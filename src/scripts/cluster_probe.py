#!/usr/bin/env python3

"""cluster_probe.py -- Discover a site's cluster dispatch settings.

This is the companion *tool* to the ``clusterrc.py`` settings file.
``clusterrc.py`` is pure data -- the parameters the dispatch generator
reads; this script *generates a starter copy* of that file by reading
what the machine can report.  It queries the scheduler (``sinfo``,
``scontrol``) and the node hardware (``lscpu``, ``numactl``) and writes
a ``clusterrc.py`` with every fact it could learn already filled in --
the performance and advanced tiers plus the partition list -- leaving
only the convention-and-policy fields (``worker_init`` and
``account``) as clearly marked blanks for the user to complete from
local documentation.

The fact-versus-policy split is deliberate (DESIGN 6.2.11): counts and
topology are machine-knowable, but how a worker brings up the imago
environment and which account to charge are site convention and
policy, which no query can report.  The tool is therefore best-effort
and scheduler-specific (SLURM today); every query is wrapped so a
missing tool or an unparseable line drops that one fact rather than
aborting, and its output is a draft the user reviews and edits, never
an authority.

This tool is *self-contained*: it carries its own copy of the schema
(:func:`_starter_schema`) and only ever *writes* a ``clusterrc.py`` --
it never reads one.  That keeps the division of labour clean
(``cluster_config`` reads the settings file, ``cluster_probe`` creates
it) and means the tool needs no ``clusterrc.py`` to exist and no
directory lookup at all.  The cost is that the key list lives in two
files; a test (``test_cluster_probe``) asserts this schema stays
identical to ``clusterrc.parameters_and_defaults()`` so they cannot
drift apart.
"""

import argparse
import os
import re
import subprocess
import sys
from datetime import datetime


def _starter_schema():
    """Return the cluster-settings schema this tool writes out.

    A self-contained copy of ``clusterrc.parameters_and_defaults()`` --
    the same keys in the same order with the same defaults -- so the
    tool can generate a complete starter file without importing (or even
    needing the existence of) a ``clusterrc.py``.  Kept identical to the
    real settings file by ``test_cluster_probe.test_schema_matches_*``;
    when a key is added to ``clusterrc.py`` it is added here too.
    """
    return {
        # Required core -- shipped as None, left for the user to fill.
        "partitions":        None,
        "worker_init":       None,
        "account":           None,
        # Performance tuning.
        "cores_per_node":    None,
        "workers_per_node":  None,
        "cores_per_worker":  1,
        "nodes":             1,
        "walltime":          "01:00:00",
        "default_topology":  "slurm-per-job",
        "max_blocks":        1,
        "memory_per_node":   None,
        "memory_per_worker": None,
        # Advanced / forward-looking.
        "launcher":          "single",
        "ranks_per_worker":  1,
        "threads_per_rank":  1,
        "binding":           None,
        "omp_places":        None,
        "omp_proc_bind":     None,
        "gpus_per_node":     0,
        "queue_overrides":   {},
        "profiles":          {},
        "extra_scheduler_options": [],
    }


# ================================================================
#  Machine queries and their parsers
#
#  The command runners are kept separate from the pure text parsers
#  so the parsers can be unit-tested on captured output without any
#  scheduler present.  Each query is guarded so a missing tool or an
#  unparseable line drops that one fact rather than aborting the
#  probe -- "fill what I could" behaviour.
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
    megabytes, and the GRES (accelerator) field.  SLURM appends ``+`` to
    a value that is a lower bound (mixed node sizes); the ``+`` is
    stripped before the integer is read.  A row that cannot be parsed is
    skipped rather than aborting the parse.

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


def parse_lscpu_topology(lscpu_text):
    """Parse ``lscpu`` output into the node's CPU topology.

    The socket / core / thread / NUMA layout reported here is exactly
    what the ``binding`` and ``omp_*`` knobs describe, so the probe
    records it to guide a power user even though the serial launcher
    does not use it yet.  Only the keys ``lscpu`` actually prints are
    returned; anything missing is simply absent.

    Parameters
    ----------
    lscpu_text : str
        The captured stdout of ``lscpu``.

    Returns
    -------
    dict
        Any of ``cpus``, ``sockets``, ``cores_per_socket``,
        ``threads_per_core``, ``numa_nodes`` that could be read.
    """
    label_to_key = {
        "CPU(s)":           "cpus",
        "Socket(s)":        "sockets",
        "Core(s) per socket": "cores_per_socket",
        "Thread(s) per core": "threads_per_core",
        "NUMA node(s)":     "numa_nodes",
    }
    topology = {}
    for line in lscpu_text.splitlines():
        if ":" not in line:
            continue
        label, _, value = line.partition(":")
        key = label_to_key.get(label.strip())
        if key is None:
            continue
        try:
            topology[key] = int(value.strip())
        except ValueError:
            continue
    return topology


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
    """Read every cluster fact the machine can report.

    Orchestrates the individual queries, each guarded by
    :func:`run_query` so a missing tool drops one fact rather than
    failing the probe.  The returned dictionary holds only the facts
    that were actually discovered; :func:`render_starter_clusterrc`
    turns it into a starter settings file.

    Returns
    -------
    dict
        Discovered facts.  Possible keys: ``partitions`` (list of queue
        names), ``cores_per_node``, ``memory_per_node``,
        ``gpus_per_node``, ``topology`` (the lscpu layout), and
        ``accounts`` (the sacctmgr hint).
    """
    facts = {}

    # Scheduler queues and per-node size, from sinfo.
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
            # Per-node facts representative of the first row; flagged in
            #   the starter file so the user reviews mixed-node queues.
            facts["cores_per_node"] = rows[0]["cores"]
            facts["memory_per_node"] = rows[0]["memory_mb"]
            # Advertise GPUs if any queue has them.
            facts["gpus_per_node"] = max(row["gpus"] for row in rows)

    # CPU topology (sockets / cores / threads / NUMA), from lscpu.
    lscpu_text = run_query(["lscpu"])
    if lscpu_text:
        topology = parse_lscpu_topology(lscpu_text)
        if topology:
            facts["topology"] = topology

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

# Discovered facts that map directly onto a settings key (the rest --
#   topology, accounts -- are hints surfaced as comments only).
_DISCOVERABLE_KEYS = ("partitions", "cores_per_node",
                      "memory_per_node", "gpus_per_node")

# The required-core keys: shipped as None, left as marked blanks until
#   the user supplies them (no machine query can).
_REQUIRED_KEYS = ("partitions", "worker_init")


def render_starter_clusterrc(discovered_facts):
    """Render a starter ``clusterrc.py`` source file as text.

    The starter offers every key in :func:`_starter_schema` (this
    tool's own copy of the settings, kept identical to the real file by
    a test), so the user sees and can tune the whole surface in one
    file.  The discoverable hardware facts overlay their defaults; the
    required-core fields the machine cannot know (``worker_init``, and
    ``partitions`` when no queue was found) are left as ``None`` with a
    ``FILL IN`` comment, so the dispatch generator's required-field
    check refuses an unfinished starter rather than dispatching it.  The
    discovered CPU topology and account hints ride along as comments,
    since they inform human choices (binding, which account) rather than
    dictate them.

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

    # Start from this tool's own schema, then overlay the hardware facts
    #   the probe can fill directly.  No clusterrc.py is read -- the test
    #   suite keeps _starter_schema() identical to the real settings.
    settings = _starter_schema()
    for key in _DISCOVERABLE_KEYS:
        if key in facts:
            settings[key] = facts[key]

    lines = [
        "#!/usr/bin/env python3",
        "",
        '"""clusterrc.py -- starter written by cluster_probe.py.',
        "",
        "Review every value below.  The discoverable fields were read",
        "off the machine; the FILL IN fields could not be -- complete",
        "them from your site's documentation before dispatching.  See",
        "the shipped clusterrc.py for the full meaning of each key.",
        '"""',
        "",
    ]

    # Surface the topology and account hints as comments.
    topology = facts.get("topology")
    if topology:
        layout = ", ".join(f"{key}={value}"
                           for key, value in topology.items())
        lines.append(f"# Discovered CPU topology: {layout}")
        lines.append("#   (informs the binding / omp_* knobs below.)")
    accounts = facts.get("accounts")
    if accounts:
        lines.append("# Accounts you may charge (sacctmgr): "
                     + ", ".join(accounts))
        lines.append("#   Set 'account' to one if the cluster "
                     "requires it.")
    if topology or accounts:
        lines.append("")

    # Emit the full parameter list -- every key the schema defines, so
    #   the user sees and can tune the whole surface in one file.
    lines.extend([
        "",
        "def parameters_and_defaults():",
        '    """Return this site\'s cluster dispatch settings."""',
        "    return {",
    ])
    for key, value in settings.items():
        blank = key in _REQUIRED_KEYS and not value
        suffix = "  # FILL IN" if blank else ""
        lines.append(f"        {key!r}: {value!r},{suffix}")
    lines.extend([
        "    }",
        "",
    ])

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
        description="Probe the scheduler and node hardware and write a "
                    "starter clusterrc.py settings file.")
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
