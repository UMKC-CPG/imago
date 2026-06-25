#!/usr/bin/env python3

"""clusterrc.py -- Per-site cluster dispatch settings for Imago.

This resource control file describes the *cluster* a flight is
dispatched to: which scheduler queues exist, what one worker must do
to bring up the imago environment, how many cores and how much memory
a node has, and the finer knobs that tune parallel placement.  Like
every other ``*rc.py`` file in Imago, it exposes a single
``parameters_and_defaults()`` function returning one dictionary, and a
user overrides the defaults by placing a customized copy in the
current working directory or in the directory pointed to by
``$IMAGO_RC``.

Why a separate file (and not part of ``imagorc``)?  These facts are
stable for a given cluster-and-account and differ from user to user,
so they belong in their own settings file rather than mixed into the
per-run input.  The generator in ``kaleidoscope`` reads this file plus
a few per-run command-line choices and assembles the Parsl ``Config``
that actually reaches the scheduler.

The settings are *tiered* so the file is approachable for a newcomer
and rewarding for a power user:

  * **Required core** -- the few facts the machine cannot report on its
    own, so they must be supplied: the queue list, the worker bring-up
    commands, and (where the cluster demands one) the account string.
    A required key carries the :data:`REQUIRED` sentinel until filled.

  * **Performance tuning** -- optional knobs that improve throughput
    (worker packing, default nodes/walltime/topology, memory guards);
    every one has a sensible built-in default.

  * **Advanced / forward-looking** -- power-user controls, several of
    them seams for a future parallel imago (the MPI/OpenMP split, CPU
    and NUMA binding, accelerator requests, named multi-cluster
    profiles, and a raw scheduler-directive escape hatch).

Discovering the settings.  Much of this file can be read straight off
the machine.  Run the discovery helper::

    python clusterrc.py --probe -o clusterrc.py

and it queries the scheduler (``sinfo``, ``scontrol``) and the node
hardware (``lscpu``, ``numactl``) and writes a *starter* settings file
with everything it could learn already filled in -- the performance
and advanced tiers plus the partition list -- leaving only the
convention-and-policy fields (``worker_init`` and ``account``) as
clearly marked blanks for you to complete from local documentation.
The helper is best-effort and SLURM-specific; its output is a draft to
review and edit, never an authority.  With no options the file simply
prints its defaults, like the other ``*rc.py`` files.
"""

import argparse
import os
import re
import subprocess
import sys


class _RequiredSentinel:
    """Marker for a settings key the user MUST supply.

    A few facts -- the queue list and the worker bring-up commands --
    cannot be guessed or read off the machine reliably, so the shipped
    defaults leave them set to the single shared :data:`REQUIRED`
    instance.  The dispatch generator (``cluster_config``) treats a
    value that is still :data:`REQUIRED` (or empty) as a configuration
    error raised up front, never a crash mid-flight.  The distinct
    ``repr`` makes an unfilled key obvious both in the defaults dump
    and in any traceback.
    """

    def __repr__(self):
        return "<REQUIRED>"


# The one shared sentinel instance.  Compared by identity, so the
#   generator must import THIS object rather than construct its own.
REQUIRED = _RequiredSentinel()


def parameters_and_defaults():
    """Return the per-site cluster dispatch settings dictionary.

    Keys are grouped into the three tiers described in the module
    docstring.  A newcomer fills only the required core; every other
    key already holds a usable default, so an omitted key is never an
    error.

    Returns
    -------
    dict
        Keys are setting names; values are the defaults (with the two
        non-discoverable required fields left as :data:`REQUIRED`).
    """
    settings = {

        # ============================================================
        #  Required core -- enough to dispatch at all.
        # ============================================================

        # The scheduler queues available to this account.  The first
        #   entry is treated as the default queue.  A starter file from
        #   --probe fills this from ``sinfo``; review the order so the
        #   queue you want as default sits first.
        "partitions":        REQUIRED,   # list[str]; [0] = default

        # The shell commands one worker runs before imago, so it can
        #   find the executable: activate the environment, load
        #   modules, set paths.  Pure site convention -- no query can
        #   report it -- so it must be supplied as a list of lines.
        "worker_init":       REQUIRED,   # list[str] of shell lines

        # The scheduler account to charge.  Required only where the
        #   cluster demands one; leave as None where it does not.
        "account":           None,

        # ============================================================
        #  Performance tuning -- optional; improves throughput.
        # ============================================================

        # Cores on one node.  Lets the generator pack several workers
        #   onto a node; when None the generator falls back to one
        #   worker per node.
        "cores_per_node":    None,

        # How many calculations run at once on a node, and how many
        #   cores each is given.  With today's serial imago a worker
        #   uses one core, so the generator packs as many workers as a
        #   node has cores.  ``workers_per_node`` None means "derive
        #   from cores_per_node // cores_per_worker".
        "workers_per_node":  None,
        "cores_per_worker":  1,

        # Defaults for the per-run choices, so a common run needs no
        #   command-line options at all.  ``default_topology`` is the
        #   dispatch shape used when --dispatch is not given.
        "nodes":             1,
        "walltime":          "01:00:00",
        "default_topology":  "slurm-per-job",

        # How many allocations the pooled shape may grow to when work
        #   backs up (the auto-scaling cap).
        "max_blocks":        1,

        # Memory guards (in megabytes) so a run neither overflows a
        #   node's memory nor silently under-requests it.  None means
        #   "let the scheduler apply its own default".
        "memory_per_node":   None,
        "memory_per_worker": None,

        # ============================================================
        #  Advanced / forward-looking -- power users; future imago.
        # ============================================================

        # How a single calculation is launched across cores or ranks.
        #   "single" is the serial launcher used today; this is the
        #   seam where an MPI / GPU launcher plugs in once imago runs
        #   in parallel.
        "launcher":          "single",

        # The hybrid-parallel split for a future parallel imago: how
        #   many MPI ranks one calculation spawns and how many OpenMP
        #   threads each rank drives.  Their product is the number of
        #   cores the calculation occupies, so the two together trade
        #   message-passing breadth against shared-memory threading.
        "ranks_per_worker":  1,
        "threads_per_rank":  1,

        # How ranks and threads are pinned to the hardware: to cores,
        #   to sockets, or to NUMA memory domains.  Pinning keeps a
        #   rank's threads on cores sharing a cache and a memory
        #   controller, which on a multi-socket node is often the
        #   difference between scaling and stalling on remote memory
        #   traffic.  None lets the scheduler place threads itself.
        "binding":           None,

        # Finer OpenMP thread placement (spread across sockets versus
        #   packed onto neighbouring cores), for tuning thread locality
        #   beyond the coarse ``binding`` choice.
        "omp_places":        None,
        "omp_proc_bind":     None,

        # Accelerators per node, for the future GPU path.  Zero means
        #   no GPU request is added to the scheduler submission.
        "gpus_per_node":     0,

        # Per-queue overrides: a mapping from a queue name to a small
        #   dictionary of settings that differ for that queue only.
        "queue_overrides":   {},

        # Named site profiles: a mapping from a profile name to a
        #   partial settings dictionary, so a user with access to
        #   several clusters selects one by name.
        "profiles":          {},

        # A raw passthrough of arbitrary scheduler directives -- the
        #   final escape hatch for anything the schema above does not
        #   name, so a power user is never blocked.
        "extra_scheduler_options": [],
    }

    return settings


# ================================================================
#  Discovery helper (--probe)
#
#  These functions read what the machine can report and write a
#  starter settings file.  Each external query is wrapped so a
#  missing tool or an unparseable line drops that one fact rather
#  than aborting the probe -- "fill what I could" behaviour.  The
#  command runners are kept separate from the pure text parsers so
#  the parsers can be unit-tested on captured output without any
#  scheduler present.
# ================================================================

def run_query(command_words):
    """Run an external query command and return its stdout, or None.

    A discovery probe must never crash the helper, so every failure
    mode -- the tool not being installed, a non-zero exit, a timeout --
    collapses to None, which the caller reads as "this fact is
    unavailable" and simply omits.

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


def render_starter_clusterrc(discovered_facts):
    """Render a starter ``clusterrc.py`` source file as text.

    Builds a complete, editable settings file: the performance and
    advanced tiers and the partition list are filled from the probe,
    while the two non-discoverable fields are left as obvious blanks --
    ``worker_init`` as an empty list and ``account`` as None, each with
    a ``FILL IN`` comment.  Until the user supplies ``worker_init`` the
    dispatch generator's required-field check refuses the file, so an
    incomplete starter never silently dispatches.  The discovered CPU
    topology and account hints are written as guiding comments rather
    than as live settings, since they inform human choices (binding,
    which account) rather than dictate them.

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

    # Header comments record what the probe found and what the user
    #   must still complete, so the starter file documents itself.
    lines = [
        "#!/usr/bin/env python3",
        "",
        '"""clusterrc.py -- starter file written by '
        '`clusterrc.py --probe`.',
        "",
        "Review every value below.  The performance and advanced",
        "tiers and the partition list were read off the machine; the",
        "FILL IN fields could not be -- complete them from your site's",
        "documentation before dispatching.",
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
        lines.append("#   Choose one for 'account' if the cluster "
                     "requires it.")
    if topology or accounts:
        lines.append("")

    # The settings function.  Only the keys worth showing in a starter
    #   are emitted; everything else keeps its built-in default.
    partitions = facts.get("partitions")
    # The trailing comma must precede any comment, so the discovered
    #   and blank cases build the whole line rather than a bare value.
    if partitions:
        partitions_line = f"        \"partitions\":        " \
                          f"{partitions!r},"
    else:
        partitions_line = "        \"partitions\":        []," \
                          "  # FILL IN: your scheduler queues"

    lines.extend([
        "",
        "def parameters_and_defaults():",
        '    """Return this site\'s cluster dispatch settings."""',
        "    return {",
        "        # --- Required core ---",
        partitions_line,
        "        \"worker_init\":       [],  # FILL IN: shell lines "
        "to find imago",
        "        \"account\":           None,  # FILL IN if the "
        "cluster requires it",
        "        # --- Performance tuning ---",
    ])
    for key in ("cores_per_node", "memory_per_node", "gpus_per_node"):
        if key in facts:
            lines.append(f"        \"{key}\":  {facts[key]!r},")
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
        print(f"clusterrc.py: {output_path} already exists; "
              "not overwriting (use --force).")
        return False
    facts = probe_site() if discovered_facts is None else discovered_facts
    with open(output_path, "w") as handle:
        handle.write(render_starter_clusterrc(facts))
    discovered = ", ".join(sorted(facts)) if facts else "nothing"
    print(f"clusterrc.py: wrote starter {output_path} "
          f"(discovered: {discovered}).")
    print("  Complete the FILL IN fields before dispatching.")
    return True


def build_arg_parser():
    """Build the command-line parser for the clusterrc utility.

    With no options the file prints its defaults, matching every other
    ``*rc.py``.  With ``--probe`` it discovers what it can and writes a
    starter settings file.
    """
    parser = argparse.ArgumentParser(
        description="Show the default cluster dispatch settings, or "
                    "probe the machine to write a starter clusterrc.py.")
    parser.add_argument(
        "--probe", action="store_true",
        help="query the scheduler and node hardware and write a "
             "starter settings file (default: just print defaults).")
    parser.add_argument(
        "-o", "--output", default="clusterrc.py",
        help="path for the starter file written by --probe "
             "(default: clusterrc.py in the current directory).")
    parser.add_argument(
        "--force", action="store_true",
        help="allow --probe to overwrite an existing output file "
             "(default: refuse and keep the existing file).")
    return parser


def main(argv=None):
    """Entry point: print defaults, or probe and write a starter file.

    Parameters
    ----------
    argv : list[str], optional
        Argument list (defaults to ``sys.argv[1:]``); accepted so the
        tests can drive the CLI directly.
    """
    args = build_arg_parser().parse_args(argv)
    if args.probe:
        wrote = write_starter_clusterrc(args.output,
                                        overwrite=args.force)
        return 0 if wrote else 1
    # Default behaviour, like the other *rc.py files: dump the
    #   defaults so a user can see every key and its built-in value.
    for key, value in parameters_and_defaults().items():
        print(f"  {key:24s} = {value!r}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
