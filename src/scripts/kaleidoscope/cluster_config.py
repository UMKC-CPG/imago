"""kaleidoscope.cluster_config -- the dispatch Config generator
(DESIGN 6.2.11; PSEUDOCODE 13.7).

A flight reaches SLURM purely through the Parsl ``Config`` it carries
(DESIGN 6.2.3); *only the Config changes* between a laptop and a
cluster.  This module is where a client obtains that Config.  It reads
two of the three configuration layers -- the per-site settings file
(``clusterrc.py``, the stable facts) and the per-run choices (the CLI
flags) -- and assembles either a Parsl ``Config`` for a cluster shape
or nothing at all for the ``local`` opt-out.  The third layer,
per-unit right-sizing, is deferred (DESIGN 6.2.11, decision 3); every
unit gets the same slice for now.

The generator lives here, in the dispatcher package every flight
client already imports, so the producer and all future flights share
one copy rather than each writing its own executor builder (DESIGN
6.2.11, decision 4).

Two cluster shapes, two Config shapes:

  * **slurm-pooled** -- one (optionally auto-scaled) allocation whose
    workers stream many units.  Best for many small, similar units --
    the convergence sweeps and the database seed.
  * **slurm-per-job** -- one scheduler submission per unit, each unit
    in its own one-worker block.  Best for large or heterogeneous
    units.  This is the command-line default (DESIGN 6.2.11).

Parallel imago is a *deferred seam*.  The advanced settings that only
matter once a single calculation spans multiple cores -- a non-serial
``launcher``, the ``ranks_per_worker`` x ``threads_per_rank`` split,
and the ``binding`` / ``omp_*`` placement knobs -- have no effect
today (imago is serial).  Rather than silently ignore them, the
builders raise a clear error if any is set away from its serial
default, so the seam is named and honest until a parallel imago and a
real MPI launcher land (TODO C100 / C81).

Parsl is imported lazily, inside the SLURM builders, so importing this
module -- and the whole ``local`` path -- never requires Parsl to be
installed (the same discipline as ``dispatch.py``).
"""

import os
import sys
from types import SimpleNamespace

from .workspace import toml_line


class ConfigError(Exception):
    """A cluster dispatch configuration is unusable -- a missing
    required field, an unknown profile, or an unknown dispatch shape.
    Raised up front, before any flight is launched, so a configuration
    mistake is a clear message rather than a crash deep in a run (the
    strict-contract discipline the producer already follows, DESIGN
    6.3.1)."""
    pass


# ----------------------------------------------------------------
#  Layer 1: the per-site settings file
# ----------------------------------------------------------------

def _clusterrc_search_dirs():
    """The directories searched for ``clusterrc.py``, highest
    precedence first.

    A ``clusterrc.py`` in the current working directory -- a per-run
    override placed beside a particular campaign -- wins; otherwise the
    one installed to ``$IMAGO_RC`` is used as the convenient global
    default.  So a routine user populates the global file once and every
    run picks it up, while a campaign that needs different queues or
    walltime drops a local copy that takes precedence for that run only.
    """
    search_dirs = [os.getcwd()]
    imago_rc = os.getenv("IMAGO_RC")
    if imago_rc:
        search_dirs.append(imago_rc)
    return search_dirs


def _load_clusterrc_module():
    """Import the active ``clusterrc`` module and return it.

    Resolves the settings file by the precedence of
    :func:`_clusterrc_search_dirs` -- current directory first, then
    ``$IMAGO_RC``.  The module is pure data -- a single
    ``parameters_and_defaults()`` -- so only that dictionary is read;
    the starter generator lives separately in ``cluster_probe.py``.
    """
    # Earlier entries on sys.path win, so insert the search directories
    #   in reverse precedence -- the highest-precedence one ends up
    #   first (at index 0).
    for candidate in reversed(_clusterrc_search_dirs()):
        if candidate and candidate not in sys.path:
            sys.path.insert(0, candidate)
    import clusterrc
    return clusterrc


def _is_empty(value):
    """True for a value that is present but carries nothing usable --
    None, or an empty string / list / dict.  An empty required field is
    treated exactly like an unfilled one."""
    if value is None:
        return True
    if isinstance(value, (str, list, tuple, dict)) and len(value) == 0:
        return True
    return False


def merge_settings(base, overlay):
    """Merge one overlay onto a settings dict, per key, one level down.

    This is the single merge every overlay uses -- the named profile,
    the per-queue override, and the per-run flags alike (DESIGN
    6.2.11, decision 1).  Most settings are a single value, and
    overlaying one simply replaces it.  But a setting may itself be a
    *block* of settings -- ``orchestrator`` and ``md`` are the two,
    holding respectively the driver's cores, memory and walltime, and
    the MD job's ranks, memory, walltime and bring-up -- and there the
    overlay names only the keys it means to change, leaving the rest
    as the layer beneath gave them.

    Replacing a whole block instead would silently discard facts the
    curator never mentioned.  A user who writes, of the debug queue,
    "the driver needs only two gigabytes there" would also lose the
    driver's core count and time limit, and would get not an error but
    two plausible-looking fallbacks in their place.

    The descent stops at one level: a block of settings holds plain
    values, never further blocks, so there is nothing deeper to merge
    and a deeper rule would only obscure what an overlay can reach.
    """
    merged = dict(base)
    for key, value in overlay.items():
        beneath = base.get(key)
        if isinstance(beneath, dict) and isinstance(value, dict):
            merged[key] = {**beneath, **value}
        else:
            merged[key] = value
    return merged


def _require_core(site):
    """Refuse a settings file whose non-discoverable core is unfilled.

    A field shipped as None (the unfilled default) or left empty --
    including a probe starter the user has not completed -- is a
    configuration error raised up front, never a crash mid-flight.
    """
    for name in ("partitions", "worker_init"):
        if _is_empty(site.get(name)):
            raise ConfigError(
                f"cluster settings file is missing required field "
                f"{name!r}; fill it in clusterrc.py "
                f"(generate a starter with cluster_probe.py).")


def load_site_config(profile=None, partition=None):
    """Load the per-site settings, with every overlay already applied.

    Reading the settings file and overlaying it are **one operation**
    (DESIGN 6.2.11, decision 1), which is why this takes the queue.
    Three layers resolve here, most general first: the built-in
    defaults, then the named profile when one is selected, then the
    override for the queue this run will use.  Only the per-run
    command-line flags, resolved by the caller, rank above them.

    The single-operation shape is deliberate.  Were the queue overlay
    a separate step, every reader of this file would have to remember
    to take it, and a reader that forgot would receive settings that
    look complete and are quietly wrong -- the cluster-wide walltime
    where a queue's cap belongs, which a scheduler answers with a
    rejected or silently truncated job rather than an error naming the
    cause.  There is more than one reader (the per-unit dispatch, and
    the driver's own batch submission).  Making un-overlaid settings
    unobtainable is what keeps them honest.

    Parameters
    ----------
    profile : str, optional
        A key into the ``profiles`` mapping; when given, that profile's
        settings overlay the base dictionary.
    partition : str, optional
        The queue this run uses.  Defaults to the first entry of the
        (profile-overlaid) ``partitions`` list, matching how
        :func:`resolve_choices` defaults it.

    Returns
    -------
    dict
        The fully resolved and overlaid site settings.
    """
    clusterrc = _load_clusterrc_module()
    site = dict(clusterrc.parameters_and_defaults())

    # Overlay 1: a named profile (advanced tier) overlays the base
    #   dict, so a user with several clusters selects one by name.
    if profile is not None:
        profiles = site.get("profiles", {})
        if profile not in profiles:
            raise ConfigError(f"unknown cluster profile {profile!r}")
        site = merge_settings(site, profiles[profile])

    # Checked BEFORE the queue overlay, because picking the default
    #   queue reads `partitions`.
    _require_core(site)

    # Overlay 2: the selected queue's own settings.  The queue is a
    #   per-run choice, defaulting to the first entry of the list.
    queue = partition or site["partitions"][0]
    site = apply_queue_overrides(site, queue)

    # An override may legitimately set worker_init, so re-check the
    #   core it could have emptied.  (It may not set partitions.)
    _require_core(site)
    return site


#: Settings a per-queue override may not set, because they choose
#:   WHICH overlay applies -- an overlay rewriting them would refer
#:   to itself (DESIGN 6.2.11, decision 1).
_OVERLAY_SELECTING_KEYS = ("partitions", "profiles")


def apply_queue_overrides(site, partition):
    """Overlay the selected queue's settings onto the site dict.

    A setting may legitimately differ by queue -- a debug queue with
    a short walltime cap, a large-memory queue with a different
    per-node capacity -- so ``queue_overrides`` maps a queue name to
    the settings that differ there.  This is the third of four
    overlays (DESIGN 6.2.11, decision 1): built-in defaults, then the
    named profile, then *this*, then the per-run command-line flags.

    It lives here rather than inside :func:`load_site_config` because
    it needs to know which queue the run uses, and the queue is
    itself a per-run choice.  The caller resolves the queue from the
    profile-overlaid file, applies this, and only then lets the
    remaining choices take their defaults from the result.

    A file may carry overrides for every queue on the cluster; the
    ones this run does not use are simply not applied.

    Raises
    ------
    ConfigError
        If an override names a setting that does not exist -- at the
        top level or inside a block -- which is almost always a typo,
        and a silently ignored typo in a resource request is exactly
        what this settings file exists to prevent; or if it tries to
        set ``partitions`` or ``profiles``.
    """
    override = site.get("queue_overrides", {}).get(partition)
    if not override:
        return site

    for key in override:
        if key in _OVERLAY_SELECTING_KEYS:
            raise ConfigError(
                f"queue override for {partition!r} may not set "
                f"{key!r}: it selects which overlay applies.")
        if key not in site:
            raise ConfigError(
                f"queue override for {partition!r} names unknown "
                f"setting {key!r}; check the spelling against "
                f"clusterrc.py.")
        # The guard descends one level, exactly as far as the merge
        #   below does, so the merge cannot reach a place the guard
        #   cannot see (DESIGN 6.2.11).  A typo inside a block is the
        #   quieter fault: it leaves the real key standing at its old
        #   value beside the stray one, so the run uses the number the
        #   curator meant to change and nothing says otherwise -- an
        #   override reading `rank` for `ranks` runs the job at the
        #   site's width while its author believes they widened it.
        if isinstance(site[key], dict) and isinstance(
                override[key], dict):
            for inner_key in override[key]:
                if inner_key not in site[key]:
                    raise ConfigError(
                        f"queue override for {partition!r} names "
                        f"unknown setting {key}.{inner_key}; check "
                        f"the spelling against clusterrc.py.")

    # Per key, one level down: a queue naming only the driver's memory
    #   keeps the site's driver cores and walltime (DESIGN 6.2.11).
    return merge_settings(site, override)


# ----------------------------------------------------------------
#  Layer 2: the per-run choices
# ----------------------------------------------------------------

def resolve_choices(site, cli):
    """Resolve the four per-run choices, each defaulting from the site.

    A flag the user did not give falls back to the site default, so a
    fully configured site needs no per-run options at all (the dispatch
    shape falls back to ``default_topology``).  ``cli`` is any object
    exposing ``dispatch`` / ``partition`` / ``nodes`` / ``walltime``
    attributes -- the parsed argparse namespace in practice.

    Returns
    -------
    dict
        Keys ``dispatch``, ``partition``, ``nodes``, ``walltime``.
    """
    return {
        "dispatch":  getattr(cli, "dispatch", None)
                     or site["default_topology"],
        "partition": getattr(cli, "partition", None)
                     or site["partitions"][0],
        "nodes":     getattr(cli, "nodes", None) or site["nodes"],
        "walltime":  getattr(cli, "walltime", None)
                     or site["walltime"],
    }


#: The keys of the site's ``orchestrator`` block, each overridable
#:   per run by an ``--orchestrator-<key>`` command-line option.
ORCHESTRATOR_KEYS = ("cores", "memory", "walltime")


def resolve_orchestrator(site, cli):
    """Resolve the driver's own resource shape (DESIGN 6.2.11).

    The site's ``orchestrator`` block is a *default* shape, not a
    fixed one: a run raises or lowers any of its keys from the
    command line, which is what keeps the settings file bounded --
    a second orchestrator with different needs overrides the shape
    for its own run rather than earning a block of its own
    (ARCHITECTURE 9.4).

    The merge is **per key**.  Overriding the memory leaves the
    site's cores and walltime standing, where replacing the whole
    block would silently discard site facts the curator never meant
    to touch.  A key nobody sets stays absent, and
    :func:`build_orchestrator_sbatch` decides what absent means:
    cores and memory simply go unrequested, while walltime falls
    back once more to the run's resolved ``--walltime`` so a driver
    job always carries a time limit.

    Note the worker-sizing flags do not reach here.  ``--walltime``
    and ``--nodes`` size one calculation; a curator shortening
    ``--walltime`` to clear a short queue is speaking about the
    calculations, not about the process that submits them.

    Parameters
    ----------
    site : dict
        The loaded per-site settings.
    cli : object
        Any object exposing ``orchestrator_cores`` /
        ``orchestrator_memory`` / ``orchestrator_walltime`` -- the
        parsed argparse namespace in practice.  A client that
        exposes none of them gets the site block unchanged.

    Returns
    -------
    dict
        The merged shape, ready for :func:`build_orchestrator_sbatch`.
    """
    shape = dict(site.get("orchestrator", {}))
    for key in ORCHESTRATOR_KEYS:
        override = getattr(cli, "orchestrator_" + key, None)
        if override is not None:
            shape[key] = override
    return shape


# ----------------------------------------------------------------
#  Assembling the Config
# ----------------------------------------------------------------

def build_dispatch_config(site, choices):
    """Turn (site facts + per-run choices) into a Parsl ``Config``.

    Returns None for the ``local`` opt-out -- the driver then runs the
    flight in process, one unit at a time (PSEUDOCODE 13.5).  Each
    cluster shape delegates to its own builder.

    Raises
    ------
    ConfigError
        If the dispatch shape is not one of the three known values.
    """
    dispatch_shape = choices["dispatch"]
    if dispatch_shape == "local":
        return None
    if dispatch_shape == "slurm-pooled":
        return build_pooled_config(site, choices)
    if dispatch_shape == "slurm-per-job":
        return build_per_job_config(site, choices)
    raise ConfigError(f"unknown dispatch shape {dispatch_shape!r}")


def _require_serial_only(site):
    """Guard the deferred parallel seam (DESIGN 6.2.11).

    Today imago is serial: one calculation runs on one core through the
    single-node launcher.  The advanced knobs that describe a parallel
    calculation -- a non-serial launcher, the MPI/OpenMP rank-and-thread
    split, and the CPU/NUMA binding and placement -- cannot yet be
    realized, so setting any of them away from its serial default is a
    clear error rather than a setting that is silently dropped.  When a
    parallel imago and a real MPI launcher land (TODO C100 / C81) this
    guard is replaced by the launcher that honours those knobs.
    """
    parallel_knobs = {
        "launcher":         "single",
        "ranks_per_worker": 1,
        "threads_per_rank": 1,
        "binding":          None,
        "omp_places":       None,
        "omp_proc_bind":    None,
    }
    for knob, serial_default in parallel_knobs.items():
        if site.get(knob, serial_default) != serial_default:
            raise NotImplementedError(
                f"cluster setting {knob!r} configures parallel "
                f"execution, which is the deferred parallel-imago seam "
                f"(DESIGN 6.2.11; TODO C100). Leave it at its serial "
                f"default ({serial_default!r}) for now.")


def make_launcher(site):
    """Return the Parsl launcher for one calculation.

    Serial today: the single-node launcher starts one process per
    worker.  The MPI/OpenMP launcher that would honour the
    ``ranks_per_worker`` x ``threads_per_rank`` split and the binding /
    ``omp_*`` placement is the deferred seam (guarded by
    :func:`_require_serial_only`), so a non-serial launcher raises.
    """
    if site.get("launcher", "single") == "single":
        from parsl.launchers import SingleNodeLauncher
        return SingleNodeLauncher()
    raise NotImplementedError(
        f"launcher {site['launcher']!r}: the MPI/OpenMP launcher is "
        f"the deferred parallel-imago seam (DESIGN 6.2.11; TODO C100). "
        f"Use launcher='single' for now.")


def scheduler_options(site, workers_per_block=1):
    """Build the raw ``#SBATCH`` directives the site settings imply.

    Assembles the directives that are meaningful for today's serial
    runs -- a memory guard, a core request, and a GPU request -- then
    appends ``extra_scheduler_options`` verbatim so a power user is
    never blocked by the schema.  (CPU/NUMA binding is applied by the
    launcher in the deferred parallel path, not as a batch directive,
    so it is not emitted here; see :func:`_require_serial_only`.)

    The memory guard is derived, not copied.  ``memory_per_worker`` is
    the memory ONE calculation needs (gigabytes -- the per-job request),
    whereas SLURM's ``--mem`` is a per-NODE figure.  A node runs
    ``workers_per_block`` calculations at once (one for the per-job
    shape, the node's packed worker count for the pooled shape), so the
    node's request is ``memory_per_worker * workers_per_block``
    gigabytes.  The separate ``memory_per_node`` is deliberately NOT
    spent here: it records the node's physical capacity and is reserved
    as a ceiling for future packing / estimation checks, not a request.

    The core request is derived the same way and for the same reason:
    the block asks for its own workers' slices and nothing wider, so
    that sibling blocks (and other users) may share the node beside it
    -- DESIGN 6.2.11, "a block asks for its slice, not for the node."
    Stating it is not optional.  A block that omits the count takes
    SLURM's one-core default, which is right for a one-worker block
    only by luck and would leave a packed pool's workers contending
    for a single core.  One task per node holds the whole worker pool
    (the dispatch driver submits with one task per node), so the
    node's cores are that task's cores and ``--cpus-per-task`` carries
    the count.  This pairs with the provider's ``exclusive=False`` in
    :func:`slurm_provider`: that declines the whole node, and this
    names what we need of it.  Neither is correct alone.

    Returns
    -------
    str
        The directives joined by newlines (Parsl prepends this string
        to the submit script), or an empty string when none apply.
    """
    directives = []
    memory_per_worker = site.get("memory_per_worker")
    if memory_per_worker:
        node_memory_gb = memory_per_worker * workers_per_block
        directives.append(f"#SBATCH --mem={node_memory_gb}G")
    node_cores = site["cores_per_worker"] * workers_per_block
    directives.append(f"#SBATCH --cpus-per-task={node_cores}")
    gpus_per_node = site.get("gpus_per_node", 0) or 0
    if gpus_per_node > 0:
        directives.append(f"#SBATCH --gres=gpu:{gpus_per_node}")
    directives.extend(site.get("extra_scheduler_options", []))
    return "\n".join(directives)


def workers_per_node(site):
    """How many workers (calculations at once) to pack onto a node.

    An explicit ``workers_per_node`` wins; otherwise it is derived from
    the node's cores and the per-worker core count; otherwise it falls
    back to one worker per node (the no-``cores_per_node`` default).
    """
    if site.get("workers_per_node") is not None:
        return site["workers_per_node"]
    if site.get("cores_per_node") is not None:
        return max(1, site["cores_per_node"] // site["cores_per_worker"])
    return 1


def slurm_provider(site, choices, *, nodes_per_block, init_blocks,
                   min_blocks, max_blocks, workers_per_block=1):
    """Build the SLURM provider shared by both cluster shapes.

    The two shapes are the *same* provider wiring with different block
    geometry; the builders below differ only in how blocks map to
    units.  The worker bring-up script lets a worker find imago; the
    account, partition, and walltime come from the resolved choices;
    the memory, GPU, and CPU knobs ride along as scheduler directives.
    ``workers_per_block`` is how many calculations a node runs at once,
    so the per-node memory and core requests both scale with it (one
    for the per-job shape, the packed worker count for the pooled
    shape).

    ``exclusive`` is stated rather than left to Parsl, whose default
    claims the entire node.  That default would undo the slice
    :func:`scheduler_options` just asked for (DESIGN 6.2.11): a
    one-core block would hold every core on the node, and sibling
    blocks would each queue for a node of their own instead of sharing
    one -- which is precisely the independent scheduling the per-job
    shape promises.  Declining the node and naming the cores are two
    halves of one request; see :func:`scheduler_options` for the other.
    """
    _require_serial_only(site)
    from parsl.providers import SlurmProvider
    return SlurmProvider(
        partition=choices["partition"],
        account=site.get("account"),
        walltime=choices["walltime"],
        nodes_per_block=nodes_per_block,
        init_blocks=init_blocks,
        min_blocks=min_blocks,
        max_blocks=max_blocks,
        worker_init="\n".join(site["worker_init"]),
        launcher=make_launcher(site),
        exclusive=False,
        scheduler_options=scheduler_options(site, workers_per_block),
    )


def build_pooled_config(site, choices):
    """One (optionally auto-scaled) allocation; many units stream
    through its workers.  The block is sized by the per-run nodes and
    the site's per-node worker packing; ``max_blocks`` lets the pool
    grow when work backs up."""
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    # One node packs this many calculations; the same count sizes both
    #   the executor's worker cap and the block's per-node memory request.
    packed_workers = workers_per_node(site)
    provider = slurm_provider(
        site, choices,
        nodes_per_block=choices["nodes"],
        init_blocks=1, min_blocks=1,
        max_blocks=site.get("max_blocks", 1),
        workers_per_block=packed_workers)
    executor = HighThroughputExecutor(
        label="imago-pooled",
        provider=provider,
        cores_per_worker=site["cores_per_worker"],
        max_workers_per_node=packed_workers)
    return Config(executors=[executor])


def build_per_job_config(site, choices):
    """One scheduler submission per unit: each unit maps to its own
    one-node, one-worker block, so calculations queue and run
    independently.  ``max_blocks`` bounds how many run at once."""
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    provider = slurm_provider(
        site, choices,
        nodes_per_block=1,
        init_blocks=0, min_blocks=0,
        max_blocks=site.get("max_blocks", 1),
        workers_per_block=1)      # one calc per node -> one worker's mem
    executor = HighThroughputExecutor(
        label="imago-per-job",
        provider=provider,
        max_workers_per_node=1)   # exactly one unit per block
    return Config(executors=[executor])


# ----------------------------------------------------------------
#  Client-side wrappers (PSEUDOCODE 13.7 run_flight)
#
#  Every flight client -- the producer, the validation harness,
#  future flights -- turns its dispatch choice into a flight Config
#  the SAME way, so the two steps live here rather than being copied
#  into each client.  The client itself keeps only what is its own:
#  building the flight and, after dispatch, harvesting the results.
# ----------------------------------------------------------------

def resolve_dispatch(dispatch_shape, partition=None, nodes=None,
                     walltime=None, profile=None):
    """Turn a client's dispatch choice into a flight Parsl ``Config``.

    The ``local`` opt-out reads no settings file and returns no config
    at all, so the driver runs every unit in process; a cluster shape
    reads the per-site settings, resolves the per-run choices against
    them, and returns a real ``Config``.  A missing settings file or
    unfilled required field surfaces as a :class:`ConfigError` here,
    up front -- never a quiet local fall-back (DESIGN 6.2.11,
    decision 2).

    Returns
    -------
    tuple
        ``(parsl_config, choices)``.  Both are None for ``local``; for
        a cluster shape ``choices`` is the resolved per-run dict, which
        :func:`write_resolved_dispatch` records beside the run.
    """
    if dispatch_shape == "local":
        return None, None
    # Passing the queue is what makes the loader overlay it; the
    #   remaining per-run choices then default from the overlaid site
    #   (DESIGN 6.2.11: defaults -> profile -> queue -> flags).
    site = load_site_config(profile, partition)
    choices = resolve_choices(site, SimpleNamespace(
        dispatch=dispatch_shape, partition=partition,
        nodes=nodes, walltime=walltime))
    return build_dispatch_config(site, choices), choices


def write_resolved_dispatch(run_dir, choices, profile):
    """Write the resolved dispatch choices beside the run.

    A small human-readable record (DESIGN 6.2.11, decision 2) so a
    cluster run is reproducible from one file: the dispatch shape, the
    queue, the node count, and the time limit actually used, plus the
    profile when one was selected.  The stable site facts are not
    duplicated here -- they live in ``clusterrc.py``, and the profile
    name pins which overlay fed this run.
    """
    os.makedirs(run_dir, exist_ok=True)
    path = os.path.join(run_dir, "resolved_dispatch.toml")
    with open(path, "w") as handle:
        handle.write("# Resolved cluster dispatch for this run.\n")
        if profile is not None:
            handle.write(toml_line("profile", profile))
        for key in ("dispatch", "partition", "nodes", "walltime"):
            handle.write(toml_line(key, choices[key]))


# ----------------------------------------------------------------
#  The orchestrator's own batch job (DESIGN 6.2.11)
# ----------------------------------------------------------------

def build_orchestrator_sbatch(site, choices, command,
                              orchestrator=None):
    """Build the text of an ``sbatch`` script that runs the producer
    (the orchestrator/driver) as its own batch job (DESIGN 6.2.11).

    The driver is a single process, so the header requests one node
    with the ``orchestrator`` resource shape -- a distinct job class
    from the per-worker sizing (ARCHITECTURE 9.4).  Under a fan-out
    dispatch it is modest (a core or two); under ``--dispatch local``
    the driver runs the SCFs in process, and that run raises the
    shape with ``--orchestrator-*`` rather than editing the site
    file.  ``account`` comes from the site and ``partition`` from the
    resolved choices; the ``worker_init`` bring-up runs first so the
    batch job finds imago, then ``command`` -- the producer
    re-invoked with ``--dispatch`` and no ``--submit``, structures
    already materialized -- runs.

    ``orchestrator`` is the shape :func:`resolve_orchestrator` merged
    from the site block and this run's flags; a caller that passes
    none gets the site block alone.  An absent ``cores`` or ``memory``
    goes unrequested, letting the scheduler apply its own default; an
    absent ``walltime`` falls back to the run's resolved walltime, so
    the driver's job always carries a time limit.

    ``command`` is the already-quoted command line to execute.
    Returns the script text; the caller writes and submits it.
    """
    if orchestrator is None:
        orchestrator = site.get("orchestrator", {})
    cores = orchestrator.get("cores", 1)
    memory = orchestrator.get("memory")
    walltime = orchestrator.get("walltime") or choices["walltime"]

    # A login shell, so the worker_init below runs in a shell whose
    #   profile has been read.  Where that bring-up uses ``module``, a
    #   plain shell would work only when the submitting shell happened
    #   to be set up already, and would fail from cron, from a workflow
    #   driver, or under ``sbatch --export=NONE`` (DESIGN 6.2.11).
    lines = ["#!/bin/bash -l",
             "#SBATCH --job-name=imago-orchestrator"]
    if site.get("account"):
        lines.append(f"#SBATCH --account={site['account']}")
    lines.append(f"#SBATCH --partition={choices['partition']}")
    lines.append("#SBATCH --nodes=1")
    lines.append(f"#SBATCH --cpus-per-task={cores}")
    if memory:
        lines.append(f"#SBATCH --mem={memory}")
    lines.append(f"#SBATCH --time={walltime}")
    # The passthrough directives are already complete lines (each
    #   carries its own "#SBATCH"), exactly as scheduler_options
    #   forwards them to Parsl, so they are copied verbatim.
    lines.extend(site.get("extra_scheduler_options", []))

    lines.append("")
    lines.extend(site["worker_init"])
    lines.append("")
    lines.append(command)
    lines.append("")
    return "\n".join(lines)


def build_md_sbatch(site, choices, command, md=None):
    """Build the text of an ``sbatch`` script that runs an external
    molecular-dynamics program under MPI (DESIGN 6.2.11).

    This is the second submission-file generator, and it sits beside
    :func:`build_orchestrator_sbatch` so a reader who finds one finds
    the other and the two stay alike where they can -- both open with
    a login shell, for the reason given there.  It differs in four
    ways, each following from the job being many MPI ranks of an
    outside program rather than one driver process: it asks for
    ``--ntasks`` instead of ``--cpus-per-task``; it runs the md
    block's own bring-up instead of ``worker_init``, and refuses to
    write a file at all when that bring-up is missing; and it pins one
    thread per rank.

    ``md`` is the shape to use; a caller that passes none gets the
    site block alone.  The sizing keys fall back where the site left
    them unset -- ranks from the node, walltime from the run's
    resolved value, memory simply unrequested -- but the bring-up does
    not, so a site with no md block at all is refused rather than
    handed a job that cannot start.

    ``command`` is the already-quoted command line to execute.
    Returns the script text; the caller writes it.

    Raises
    ------
    ConfigError
        When the md block records no bring-up.
    """
    # Test for ABSENCE, not falsiness, exactly as the orchestrator
    #   generator does: a caller passing an empty shape is saying
    #   "request nothing but the fallbacks" and must not silently
    #   re-inherit the site block.
    shape = site.get("md", {}) if md is None else md
    memory = shape.get("memory")                  # None -> no --mem
    walltime = shape.get("walltime") or choices["walltime"]

    # The bring-up is required HERE rather than in the loader's
    #   required core.  Without it nothing puts the MD program on the
    #   path and the job cannot start -- but a flight that never
    #   condenses must not be refused over it, so this generator is
    #   the one that insists.  Same emptiness test the loader uses,
    #   so None and [] are both unfilled.
    bring_up = shape.get("init")
    if _is_empty(bring_up):
        raise ConfigError(
            "cluster settings file records no md bring-up "
            "('init' in the 'md' block), so the generated job "
            "would have no way to find the MD program; fill it in "
            "clusterrc.py (generate a starter with cluster_probe.py).")

    # Ranks come from the node, never from a number written into the
    #   source.  Where the site recorded no core count there is
    #   nothing to derive from, so ask for ONE rank and say so in the
    #   file itself: a one-rank MD job is visibly wrong to whoever
    #   opens it, whereas a guessed count would run, and run wrong,
    #   without ever announcing that the site was never configured.
    ranks = shape.get("ranks") or site.get("cores_per_node")
    unsized = ranks is None
    if unsized:
        ranks = 1

    lines = ["#!/bin/bash -l", "#SBATCH --job-name=lmp"]
    if site.get("account"):                       # some sites need none
        lines.append(f"#SBATCH --account={site['account']}")
    lines.append(f"#SBATCH --partition={choices['partition']}")
    lines.append("#SBATCH --nodes=1")
    lines.append(f"#SBATCH --ntasks={ranks}")
    if memory:
        lines.append(f"#SBATCH --mem={memory}")
    lines.append(f"#SBATCH --time={walltime}")
    # The passthrough directives are already complete lines, exactly
    #   as the orchestrator generator forwards them.
    lines.extend(site.get("extra_scheduler_options", []))

    # The directives are done, so a comment here cannot swallow one.
    if unsized:
        lines.append("")
        lines.append("# One rank only: this site's settings file")
        lines.append("#   records no cores_per_node, so there was")
        lines.append("#   nothing to size this job from.  Set it,")
        lines.append("#   or the md block's ranks, and rerun.")

    # The bring-up runs first so the batch job can find the MD
    #   program.  The thread pin comes AFTER it, so that a module
    #   setting a thread count of its own cannot overwrite it: the
    #   ranks were sized to fill the node, so each must hold one
    #   core, and a threaded BLAS left to itself would start a thread
    #   per core in EVERY rank -- on a forty-core node, sixteen
    #   hundred threads contending for forty cores.
    lines.append("")
    lines.extend(bring_up)
    lines.append("")
    lines.append("export OMP_NUM_THREADS=1")
    lines.append("")
    lines.append(command)
    lines.append("")
    # Deliberately NOT written, both to stay alike with the
    #   orchestrator generator: --output/--error, leaving the
    #   scheduler's own default naming; and a cd to the submit
    #   directory, which SLURM has already done.
    return "\n".join(lines)


def condense_write_submission(site):
    """Write the SLURM submission file for a condensation run.

    Called at the close of ``condense.py``'s ``create_lammps_files``,
    which has already entered the ``lammps/`` directory, so the file
    lands beside the input it submits.

    The script reads ``condenserc.py`` for the settings that are its
    own business and the site file for the cluster facts that are not
    (DESIGN 6.2.11), and it has no dispatch flags of its own -- so
    every choice falls through to the site default.  Those defaults
    are reached through the same resolver a flight uses, handed a flag
    set with nothing in it, rather than by reading ``partitions[0]``
    out of the site directly: that is the case :func:`resolve_choices`
    already serves for a fully configured site passing no flags, and
    going through it means the day ``condense.py`` grows a
    ``--partition`` of its own, the resolution needs no redesign.

    The site arrives already loaded, from the script's settings time,
    and nothing is read from disk here.  Reading it at this point
    instead would be wrong twice over: an unconfigured site would be
    refused only after the bonds, the angles and the whole LAMMPS
    input had been computed, throwing that work away; and by now the
    script has entered ``lammps/``, so a loader that searches the
    current directory first could resolve a different settings file
    than the one the run started under.

    Parameters
    ----------
    site : dict
        The fully overlaid site settings, as returned by
        :func:`load_site_config` when the script read its settings.

    Raises
    ------
    ConfigError
        When the site records no md bring-up.  The refusal for an
        unfilled site fires earlier, at settings time.
    """
    choices = resolve_choices(site, SimpleNamespace())
    # The one line naming the MD program; everything above it in the
    #   generated file came from the site.
    command = 'mpirun -np "$SLURM_NTASKS" lmp -in lammps.in'
    with open("slurm", "w") as submission_file:
        submission_file.write(build_md_sbatch(site, choices, command))
