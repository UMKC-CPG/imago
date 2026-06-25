#!/usr/bin/env python3

"""clusterrc.py -- Per-site cluster dispatch settings for Imago.

This resource control file describes the *cluster* a flight is
dispatched to: which scheduler queues exist, what one worker must do
to bring up the imago environment, how many cores and how much memory
a node has, and the finer knobs that tune parallel placement.  Like
every other ``*rc.py`` file in Imago, it is pure data -- a single
``parameters_and_defaults()`` function returning one dictionary -- and
a user overrides the defaults by placing a customized copy in the
current working directory or in the directory pointed to by
``$IMAGO_RC``.  Running the file with no arguments simply prints the
defaults, exactly like the other ``*rc.py`` files.

Why a separate file (and not part of ``imagorc``)?  These facts are
stable for a given cluster-and-account and differ from user to user,
so they belong in their own settings file rather than mixed into the
per-run input.  The generator in ``kaleidoscope.cluster_config`` reads
this file plus a few per-run command-line choices and assembles the
Parsl ``Config`` that actually reaches the scheduler.

The settings are *tiered* so the file is approachable for a newcomer
and rewarding for a power user:

  * **Required core** -- the few facts the machine cannot report on its
    own, so they must be supplied: the queue list, the worker bring-up
    commands, and (where the cluster demands one) the account string.
    A required key is left as ``None`` with a ``REQUIRED`` comment until
    filled; the dispatch generator refuses to run until it is.

  * **Performance tuning** -- optional knobs that improve throughput
    (worker packing, default nodes/walltime/topology, memory guards);
    every one has a sensible built-in default.

  * **Advanced / forward-looking** -- power-user controls, several of
    them seams for a future parallel imago (the MPI/OpenMP split, CPU
    and NUMA binding, accelerator requests, named multi-cluster
    profiles, and a raw scheduler-directive escape hatch).

Filling it in.  Much of this file can be read straight off the
machine: the companion tool ``cluster_probe.py`` queries the scheduler
(``sinfo``, ``scontrol``) and the node hardware (``lscpu``,
``numactl``) and writes a *starter* copy of this file with everything
it could learn already filled in, leaving only the
convention-and-policy fields (``worker_init`` and ``account``) as
clearly marked blanks for you to complete from local documentation.
"""


def parameters_and_defaults():
    """Return the per-site cluster dispatch settings dictionary.

    Keys are grouped into the three tiers described in the module
    docstring.  A newcomer fills only the required core; every other
    key already holds a usable default, so an omitted key is never an
    error.

    Returns
    -------
    dict
        Keys are setting names; values are the defaults.  The two
        non-discoverable required fields ship as ``None`` and must be
        filled before a cluster dispatch will run.
    """
    settings = {

        # ============================================================
        #  Required core -- enough to dispatch at all.
        # ============================================================

        # The scheduler queues available to this account; the first
        #   entry is the default queue.  REQUIRED: a list of queue
        #   names, e.g. ["general", "gpu"].  cluster_probe.py fills this
        #   from sinfo; review the order so your preferred default queue
        #   sits first.
        "partitions":        None,       # REQUIRED: list[str]

        # The shell commands one worker runs before imago, so it can
        #   find the executable: activate the environment, load modules,
        #   set paths.  Pure site convention -- no query can report it.
        #   REQUIRED: a list of shell lines, e.g.
        #   ["module load imago", "source <venv>/bin/activate"].
        "worker_init":       None,       # REQUIRED: list[str]

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


if __name__ == "__main__":
    for key, value in parameters_and_defaults().items():
        print(f"  {key:24s} = {value!r}")
