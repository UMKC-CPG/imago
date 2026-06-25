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

def _load_clusterrc_module():
    """Import the active ``clusterrc`` module and return it.

    Follows the same override search the other Imago scripts use: a
    customized ``clusterrc.py`` in the current directory or in
    ``$IMAGO_RC`` takes precedence over the shipped default.  The module
    is pure data -- a single ``parameters_and_defaults()`` -- so only
    that dictionary is read; the starter generator lives separately in
    ``cluster_probe.py``.
    """
    for candidate in (os.getcwd(), os.getenv("IMAGO_RC")):
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


def load_site_config(profile=None):
    """Load the per-site settings, optionally overlaying a profile.

    Reads ``clusterrc.parameters_and_defaults()``; if a named profile
    is requested (the advanced multi-cluster tier) its partial settings
    are overlaid on the base.  The required core -- the queue list and
    the worker bring-up commands -- must be present and non-empty; a
    gap is a :class:`ConfigError` raised here, up front, never a crash
    mid-flight.

    Parameters
    ----------
    profile : str, optional
        A key into the ``profiles`` mapping; when given, that profile's
        settings overlay the base dictionary.

    Returns
    -------
    dict
        The fully resolved site settings.
    """
    clusterrc = _load_clusterrc_module()
    site = dict(clusterrc.parameters_and_defaults())

    # A named profile (advanced tier) overlays the base dict, so a user
    #   with several clusters selects one by name.
    if profile is not None:
        profiles = site.get("profiles", {})
        if profile not in profiles:
            raise ConfigError(f"unknown cluster profile {profile!r}")
        site = {**site, **profiles[profile]}

    # The required core must be present.  A field shipped as None (the
    #   unfilled default) or left empty -- including a probe starter the
    #   user has not completed -- is a configuration error raised here.
    for name in ("partitions", "worker_init"):
        if _is_empty(site.get(name)):
            raise ConfigError(
                f"cluster settings file is missing required field "
                f"{name!r}; fill it in clusterrc.py "
                f"(generate a starter with cluster_probe.py).")
    return site


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


def scheduler_options(site):
    """Build the raw ``#SBATCH`` directives the site settings imply.

    Assembles the directives that are meaningful for today's serial
    runs -- a memory guard and a GPU request -- then appends
    ``extra_scheduler_options`` verbatim so a power user is never
    blocked by the schema.  (CPU/NUMA binding is applied by the
    launcher in the deferred parallel path, not as a batch directive,
    so it is not emitted here; see :func:`_require_serial_only`.)

    Returns
    -------
    str
        The directives joined by newlines (Parsl prepends this string
        to the submit script), or an empty string when none apply.
    """
    directives = []
    memory_per_node = site.get("memory_per_node")
    if memory_per_node:
        directives.append(f"#SBATCH --mem={memory_per_node}")
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
                   min_blocks, max_blocks):
    """Build the SLURM provider shared by both cluster shapes.

    The two shapes are the *same* provider wiring with different block
    geometry; the builders below differ only in how blocks map to
    units.  The worker bring-up script lets a worker find imago; the
    account, partition, and walltime come from the resolved choices;
    the memory, GPU, and CPU knobs ride along as scheduler directives.
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
        scheduler_options=scheduler_options(site),
    )


def build_pooled_config(site, choices):
    """One (optionally auto-scaled) allocation; many units stream
    through its workers.  The block is sized by the per-run nodes and
    the site's per-node worker packing; ``max_blocks`` lets the pool
    grow when work backs up."""
    from parsl.config import Config
    from parsl.executors import HighThroughputExecutor
    provider = slurm_provider(
        site, choices,
        nodes_per_block=choices["nodes"],
        init_blocks=1, min_blocks=1,
        max_blocks=site.get("max_blocks", 1))
    executor = HighThroughputExecutor(
        label="imago-pooled",
        provider=provider,
        cores_per_worker=site["cores_per_worker"],
        max_workers_per_node=workers_per_node(site))
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
        max_blocks=site.get("max_blocks", 1))
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
    site = load_site_config(profile)
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
