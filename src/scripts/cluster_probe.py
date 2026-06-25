#!/usr/bin/env python3

"""cluster_probe.py -- Discover a site's cluster dispatch settings.

This is the companion *tool* to the ``clusterrc.py`` settings file.
``clusterrc.py`` is pure data -- the parameters the dispatch generator
reads; this script *creates a starter copy* of that file by asking the
scheduler what it knows.  It runs ``sinfo`` (the queues and the size of
their nodes) and ``sacctmgr`` (the accounts you may charge), and writes
a ``clusterrc.py`` with what it could learn already filled in and a
brief plain-language note on every setting.

What the scheduler cannot report is convention and policy: how a worker
brings up the Imago environment (``worker_init``) and which account to
charge are site-specific and human-chosen, so those are left as clearly
marked ``FILL IN`` blanks.

Two honesty rules matter here (DESIGN 6.2.11):

  * **Only the scheduler is trusted, not the login node.**  ``sinfo``
    reports the *compute* nodes' core and memory counts; this tool does
    NOT read the login node's own CPU layout, which would describe the
    wrong machine -- so no login-node facts ever reach the file.

  * **A heterogeneous cluster is not guessed at.**  When the cluster's
    nodes differ in core count or memory, the tool does not silently
    pick one value; it leaves that setting blank and lists the values
    it saw in a comment, so the user chooses deliberately.

This tool is *self-contained*: it carries its own copy of the schema
(:func:`_starter_schema`) and only ever *writes* a ``clusterrc.py`` --
it never reads one.  A test (``test_cluster_probe``) keeps this schema
identical to ``clusterrc.parameters_and_defaults()`` so they cannot
drift apart.
"""

import argparse
import os
import re
import subprocess
import sys
import textwrap
from datetime import datetime


# The cluster-settings schema this tool writes, as (key, default,
#   one-line plain-language note) triples.  The keys, order, and
#   defaults match clusterrc.parameters_and_defaults() (a test enforces
#   it); the note is shown beside each setting in the generated file so
#   a newcomer can read the file without cross-referencing anything.
_SETTINGS = [
    ("partitions", None,
     "The scheduler queues you can submit to (first is default)."),
    ("worker_init", None,
     "Shell commands that set up Imago before a job runs."),
    ("account", None,
     "The account charged for compute time, if one is required."),
    ("cores_per_node", None,
     "How many CPU cores one node has."),
    ("workers_per_node", None,
     "Calculations to run at once on a node (blank = auto)."),
    ("cores_per_worker", 1,
     "How many cores each calculation uses (one for now)."),
    ("nodes", 1,
     "How many nodes to request per job."),
    ("walltime", "01:00:00",
     "The most time one job may run, as HH:MM:SS."),
    ("default_topology", "slurm-per-job",
     "Default work spread: one job per calc, or one allocation."),
    ("max_blocks", 1,
     "Most allocations to grow to at once when work piles up."),
    ("memory_per_node", None,
     "How much memory one node has, in megabytes."),
    ("memory_per_worker", None,
     "How much memory each calculation needs, in megabytes."),
    ("launcher", "single",
     "How one calculation starts ('single' = a serial run)."),
    ("ranks_per_worker", 1,
     "Parallel pieces per calculation (one for now)."),
    ("threads_per_rank", 1,
     "Threads per parallel piece (one for now)."),
    ("binding", None,
     "How to pin a calculation to particular cores (off for now)."),
    ("omp_places", None,
     "Fine control over where threads run (off for now)."),
    ("omp_proc_bind", None,
     "Fine control over how threads are pinned (off for now)."),
    ("gpus_per_node", 0,
     "How many GPUs per node to request (zero = none)."),
    ("queue_overrides", {},
     "Settings that should differ for one particular queue."),
    ("profiles", {},
     "Named setting groups, one per cluster you use."),
    ("extra_scheduler_options", [],
     "Extra raw scheduler directives to add."),
]

# The required-core keys: shipped blank, because no scheduler query can
#   supply them -- the user must.
_REQUIRED_KEYS = ("partitions", "worker_init")

# Per-node numbers the scheduler can fill, each mapped to the fact key
#   that holds the distinct values seen when the cluster's nodes
#   disagree, plus the column in a parsed sinfo row and a label for the
#   "nodes vary" note.
_PER_NODE = {
    "cores_per_node":  ("core_options", "cores", "core counts"),
    "memory_per_node": ("mem_options", "memory_mb",
                        "memory sizes (MB)"),
    "gpus_per_node":   ("gpu_options", "gpus", "GPU counts"),
}


def _starter_schema():
    """Return this tool's copy of the settings as a ``{key: default}``
    dict.  Derived from :data:`_SETTINGS`; a test keeps it identical to
    ``clusterrc.parameters_and_defaults()`` so the two cannot drift."""
    return {key: default for key, default, _ in _SETTINGS}


# ================================================================
#  Scheduler queries and their parsers
#
#  The command runner is kept separate from the pure text parsers so
#  the parsers can be unit-tested on captured output without any
#  scheduler present.  Each query is guarded so a missing tool or an
#  unparseable line drops that one fact rather than aborting the probe
#  -- "fill what I could" behaviour.
# ================================================================

def run_query(command_words):
    """Run an external query command and return its stdout, or None.

    A discovery probe must never crash, so every failure mode -- the
    tool not being installed, a non-zero exit, a timeout -- collapses
    to None, which the caller reads as "this fact is unavailable" and
    simply omits.

    Parameters
    ----------
    command_words : list[str]
        The command and its arguments, e.g. ``["sinfo", "-h"]``.

    Returns
    -------
    str or None
        The captured standard output, or None on any failure.
    """
    try:
        completed = subprocess.run(command_words, capture_output=True,
                                   text=True, timeout=15, check=True)
    except (OSError, subprocess.SubprocessError):
        return None
    return completed.stdout


def parse_gres_gpu_count(gres_field):
    """Return the GPU count encoded in a SLURM GRES field.

    SLURM reports generic resources (GRES) in fields such as
    ``gpu:a100:4`` or ``gpu:4(S:0-1)`` or ``(null)`` when none are
    present.  This pulls out the integer count for the ``gpu`` resource
    and returns 0 when no GPU resource is named.

    Parameters
    ----------
    gres_field : str
        One node's GRES column from ``sinfo``.

    Returns
    -------
    int
        The number of GPUs per node, or 0 if none are advertised.
    """
    if not gres_field or gres_field.lower().startswith("(null)"):
        return 0
    # Match ``gpu`` optionally followed by a type, then the count, e.g.
    #   gpu:4, gpu:a100:4, gpu:a100:4(S:0-1).  The count is the last
    #   colon-separated integer before any trailing parenthesis.
    match = re.search(r"gpu(?::[^:,()]+)*:(\d+)", gres_field)
    return int(match.group(1)) if match else 0


def parse_sinfo_rows(sinfo_text):
    """Parse ``sinfo`` partition rows into a list of node-fact dicts.

    The probe runs ``sinfo -h -o "%R %c %m %G"`` so each row carries the
    partition name, the cores per node, the memory per node in
    megabytes, and the GRES (accelerator) field.  A heterogeneous
    cluster yields several rows -- one per distinct node configuration.
    SLURM appends ``+`` to a value that is a lower bound (mixed node
    sizes); the ``+`` is stripped before the integer is read.  A row
    that cannot be parsed is skipped rather than aborting the parse.

    Parameters
    ----------
    sinfo_text : str
        The captured stdout of the ``sinfo`` query.

    Returns
    -------
    list[dict]
        One dict per row with keys ``partition``, ``cores``,
        ``memory_mb``, and ``gpus``.
    """
    rows = []
    for line in sinfo_text.splitlines():
        fields = line.split()
        if len(fields) < 3:
            continue
        partition = fields[0]
        try:
            cores = int(fields[1].rstrip("+"))
            memory_mb = int(fields[2].rstrip("+"))
        except ValueError:
            continue
        gres_field = fields[3] if len(fields) > 3 else ""
        rows.append({"partition": partition, "cores": cores,
                     "memory_mb": memory_mb,
                     "gpus": parse_gres_gpu_count(gres_field)})
    return rows


def parse_sacctmgr_accounts(sacctmgr_text):
    """Parse the accounts a user is permitted to use.

    The probe runs ``sacctmgr -nP show assoc`` whose parseable (``-P``)
    output is ``|``-delimited; the first field is the account.  The
    distinct accounts are returned in first-seen order as a *hint* for
    the ``account`` field -- which one to actually charge remains the
    user's policy decision, so the probe never picks for them.

    Parameters
    ----------
    sacctmgr_text : str
        The captured stdout of the ``sacctmgr`` query.

    Returns
    -------
    list[str]
        The distinct account names, in first-seen order.
    """
    accounts = []
    for line in sacctmgr_text.splitlines():
        fields = line.split("|")
        account = fields[0].strip() if fields else ""
        if account and account not in accounts:
            accounts.append(account)
    return accounts


def probe_site():
    """Ask the scheduler what it knows about this site.

    Runs only scheduler queries -- never the login node's own hardware,
    which would describe the wrong machine.  For each per-node number
    (cores, memory, GPUs): if every node agrees, the value is filled in;
    if the nodes disagree (a heterogeneous cluster), the value is *not*
    guessed -- the distinct values seen are returned under an
    ``*_options`` key for the starter to list, and the setting is left
    for the user.

    Returns
    -------
    dict
        Some of: ``partitions``; ``cores_per_node`` /
        ``memory_per_node`` / ``gpus_per_node`` (only when the nodes
        agree); ``core_options`` / ``mem_options`` / ``gpu_options``
        (the distinct values seen, when they do not); and ``accounts``.
    """
    facts = {}

    # Queues and per-node sizes, from sinfo (the compute nodes).
    sinfo_text = run_query(["sinfo", "-h", "-o", "%R %c %m %G"])
    if sinfo_text:
        rows = parse_sinfo_rows(sinfo_text)
        if rows:
            # Ordered, de-duplicated queue list; the user reorders so
            #   the preferred default queue sits first.
            partitions = []
            for row in rows:
                if row["partition"] not in partitions:
                    partitions.append(row["partition"])
            facts["partitions"] = partitions
            # Each per-node number: fill it in when the nodes agree;
            #   otherwise record the distinct values and leave it blank.
            for setting, (opt_key, row_key, _) in _PER_NODE.items():
                values = sorted({row[row_key] for row in rows})
                if len(values) == 1:
                    facts[setting] = values[0]
                else:
                    facts[opt_key] = values

    # Accounts the user may charge, from sacctmgr -- a hint only.
    user_name = os.getenv("USER", "")
    if user_name:
        sacctmgr_text = run_query(["sacctmgr", "-nP", "show", "assoc",
                                   "user=" + user_name,
                                   "format=Account"])
        if sacctmgr_text:
            accounts = parse_sacctmgr_accounts(sacctmgr_text)
            if accounts:
                facts["accounts"] = accounts

    return facts


# ================================================================
#  Rendering and writing the starter clusterrc.py
# ================================================================

def _comment_lines(text):
    """Wrap a comment into ``        # ...`` lines, each within the
    80-character convention, so even a long "nodes vary" note stays
    readable in the generated file."""
    body_width = 78 - len("        # ")
    pieces = textwrap.wrap(text, width=body_width) or [""]
    return [f"        # {piece}" for piece in pieces]


def _value_lines(key, value, suffix):
    """Render a ``'key': value,`` line, wrapping a long list value
    across several indented lines so nothing exceeds the line-length
    convention (a scalar that is still too long is left as one line)."""
    single = f"        {key!r}: {value!r},{suffix}"
    if len(single) <= 78 or not isinstance(value, (list, tuple)):
        return [single]
    lines = [f"        {key!r}: ["]
    elements = ", ".join(repr(item) for item in value)
    for piece in textwrap.wrap(elements, width=66):
        lines.append(f"            {piece}")
    lines.append(f"        ],{suffix}")
    return lines


def render_starter_clusterrc(discovered_facts):
    """Render a starter ``clusterrc.py`` source file as text.

    Emits every setting in :data:`_SETTINGS`, each preceded by its
    plain-language note, with the discovered value where the scheduler
    could supply one.  A ``FILL IN`` marker is added to anything the
    user must still decide -- the required core, plus any per-node
    number the cluster disagreed on (a "nodes vary" note then lists the
    values seen).  The account hint from ``sacctmgr`` rides along on the
    ``account`` line.  No login-node facts are written.

    Parameters
    ----------
    discovered_facts : dict
        The output of :func:`probe_site` (possibly empty).

    Returns
    -------
    str
        The full text of a starter ``clusterrc.py``.
    """
    facts = discovered_facts or {}
    lines = [
        "#!/usr/bin/env python3",
        "",
        '"""clusterrc.py -- starter written by cluster_probe.py.',
        "",
        "Review every value below.  A FILL IN setting could not be",
        "discovered -- complete it from your site's documentation.",
        '"""',
        "",
        "",
        "def parameters_and_defaults():",
        '    """Return this site\'s cluster dispatch settings."""',
        "    return {",
    ]
    for key, default, note in _SETTINGS:
        value = facts.get(key, default)
        comments = [note]
        # A per-node number the nodes disagreed on is left blank, with
        #   the values seen listed so the user picks deliberately.
        options = None
        if key in _PER_NODE:
            options = facts.get(_PER_NODE[key][0])
        if options is not None:
            label = _PER_NODE[key][2]
            seen = ", ".join(str(value_seen) for value_seen in options)
            comments.append(f"Nodes vary -- {label} seen: {seen}. "
                            "Choose one.")
        # The discovered accounts are a hint for the account choice.
        if key == "account" and facts.get("accounts"):
            comments.append("You may use: "
                            + ", ".join(facts["accounts"]) + ".")
        blank_required = key in _REQUIRED_KEYS and not value
        for comment in comments:
            lines.extend(_comment_lines(comment))
        suffix = "  # FILL IN" if (blank_required or options) else ""
        lines.extend(_value_lines(key, value, suffix))
    lines.extend(["    }", ""])
    return "\n".join(lines) + "\n"


def write_starter_clusterrc(output_path, discovered_facts=None, *,
                            overwrite=False):
    """Probe the machine and write a starter ``clusterrc.py``.

    Parameters
    ----------
    output_path : str
        Where to write the starter file.
    discovered_facts : dict, optional
        Pre-supplied facts (used by the tests); when omitted the
        machine is probed live via :func:`probe_site`.
    overwrite : bool
        When False (the default) an existing file is left untouched and
        a message is printed, so a probe never clobbers a settings file
        the user has already edited.

    Returns
    -------
    bool
        True if a file was written, False if it already existed and was
        kept.
    """
    if os.path.exists(output_path) and not overwrite:
        print(f"cluster_probe: {output_path} already exists; "
              "not overwriting (use --force).")
        return False
    facts = probe_site() if discovered_facts is None else discovered_facts
    with open(output_path, "w") as handle:
        handle.write(render_starter_clusterrc(facts))
    discovered = ", ".join(sorted(facts)) if facts else "nothing"
    print(f"cluster_probe: wrote starter {output_path} "
          f"(discovered: {discovered}).")
    print("  Complete the FILL IN fields before dispatching.")
    return True


def build_arg_parser():
    """Build the command-line parser for the cluster-probe tool."""
    parser = argparse.ArgumentParser(
        description="Probe the scheduler and write a starter "
                    "clusterrc.py settings file.")
    parser.add_argument(
        "-o", "--output", default="clusterrc.py",
        help="path for the starter file "
             "(default: clusterrc.py in the current directory).")
    parser.add_argument(
        "--force", action="store_true",
        help="overwrite an existing output file "
             "(default: refuse and keep the existing file).")
    return parser


def main(argv=None):
    """Entry point: probe the machine and write a starter clusterrc.py.

    Parameters
    ----------
    argv : list[str], optional
        Argument list (defaults to ``sys.argv[1:]``); accepted so the
        tests can drive the CLI directly.
    """
    args = build_arg_parser().parse_args(argv)
    wrote = write_starter_clusterrc(args.output, overwrite=args.force)
    return 0 if wrote else 1


def record_command():
    """Append the issued command line to a file named "command" in
    the current directory, so the exact invocation can be recovered
    later.  This is a standing project convention: each run appends
    a dated block, so the file builds up a history of how the script
    was called."""

    with open("command", "a") as cmd:
        now = datetime.now()
        stamp = now.strftime("%b. %d, %Y: %H:%M:%S")
        cmd.write(f"Date: {stamp}\n")
        cmd.write("Cmnd:")
        for argument in sys.argv:
            cmd.write(f" {argument}")
        cmd.write("\n\n")


if __name__ == "__main__":
    record_command()
    sys.exit(main())
