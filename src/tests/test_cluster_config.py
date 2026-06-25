"""Tests for kaleidoscope.cluster_config -- the dispatch ``Config``
generator (DESIGN 6.2.11; PSEUDOCODE 13.7).

The pure layers -- loading and validating the site settings, resolving
the per-run choices, the ``local`` short-circuit, the worker-packing
arithmetic, the scheduler-directive string, and the deferred-parallel
guard -- are tested without Parsl.  The two cluster shapes additionally
build a real Parsl ``Config`` and inspect its block geometry; those
tests skip cleanly if Parsl is not installed.  The import works because
conftest.py inserts ``src/scripts`` on ``sys.path``.
"""

import tomllib
from types import SimpleNamespace

import pytest

from kaleidoscope import cluster_config as cc

pytestmark = pytest.mark.unit

# Parsl is only needed by the two Config-building shapes; skip those
#   tests gracefully where it is absent.
try:
    import parsl  # noqa: F401
    HAVE_PARSL = True
except Exception:
    HAVE_PARSL = False
needs_parsl = pytest.mark.skipif(not HAVE_PARSL,
                                 reason="Parsl is not installed")


def _filled_site(**overrides):
    """A fully filled site dictionary, as a user's clusterrc would
    yield, with optional per-test overrides."""
    site = {
        "partitions": ["general", "gpu"],
        "worker_init": ["module load imago", "source venv/bin/activate"],
        "account": "rulisp-lab",
        "cores_per_node": 32, "workers_per_node": None,
        "cores_per_worker": 1, "nodes": 1, "walltime": "02:00:00",
        "default_topology": "slurm-per-job", "max_blocks": 4,
        "memory_per_node": None, "memory_per_worker": None,
        "launcher": "single", "ranks_per_worker": 1,
        "threads_per_rank": 1, "binding": None, "omp_places": None,
        "omp_proc_bind": None, "gpus_per_node": 0,
        "queue_overrides": {}, "profiles": {},
        "extra_scheduler_options": [],
    }
    site.update(overrides)
    return site


# ----------------------------------------------------------------
#  Layer 1: load_site_config and its required-field guard
# ----------------------------------------------------------------

def test_load_site_config_reads_shipped_defaults_and_rejects_them():
    """The shipped clusterrc leaves the required core as None, so a
    bare load is a ConfigError -- the user must fill the file first."""
    with pytest.raises(cc.ConfigError) as excinfo:
        cc.load_site_config()
    assert "partitions" in str(excinfo.value) or \
        "worker_init" in str(excinfo.value)


def test_load_site_config_accepts_a_filled_file(monkeypatch):
    """A clusterrc whose required core is filled loads cleanly and the
    settings come straight through."""
    filled = _filled_site()
    fake_rc = SimpleNamespace(
        parameters_and_defaults=lambda: dict(filled))
    monkeypatch.setattr(cc, "_load_clusterrc_module", lambda: fake_rc)
    site = cc.load_site_config()
    assert site["partitions"] == ["general", "gpu"]
    assert site["account"] == "rulisp-lab"


def test_load_site_config_rejects_empty_required_field(monkeypatch):
    """An empty required field is rejected -- the case of a probe
    starter the user has not completed (worker_init still blank)."""
    broken = _filled_site(worker_init=[])
    fake_rc = SimpleNamespace(
        parameters_and_defaults=lambda: dict(broken))
    monkeypatch.setattr(cc, "_load_clusterrc_module", lambda: fake_rc)
    with pytest.raises(cc.ConfigError):
        cc.load_site_config()


def test_load_site_config_applies_named_profile(monkeypatch):
    """A named profile overlays the base settings; an unknown profile
    is a ConfigError."""
    base = _filled_site(profiles={"hpc2": {"account": "other-acct",
                                           "partitions": ["batch"]}})
    fake_rc = SimpleNamespace(
        parameters_and_defaults=lambda: dict(base))
    monkeypatch.setattr(cc, "_load_clusterrc_module", lambda: fake_rc)
    site = cc.load_site_config(profile="hpc2")
    assert site["account"] == "other-acct"
    assert site["partitions"] == ["batch"]
    with pytest.raises(cc.ConfigError):
        cc.load_site_config(profile="nope")


# ----------------------------------------------------------------
#  Layer 2: resolve_choices
# ----------------------------------------------------------------

def test_resolve_choices_defaults_every_flag_from_the_site():
    """When no flag is given each choice falls back to the site, so the
    dispatch shape comes from default_topology and the queue from the
    first partition."""
    cli = SimpleNamespace(dispatch=None, partition=None, nodes=None,
                          walltime=None)
    choices = cc.resolve_choices(_filled_site(), cli)
    assert choices == {"dispatch": "slurm-per-job",
                       "partition": "general", "nodes": 1,
                       "walltime": "02:00:00"}


def test_resolve_choices_lets_flags_win():
    """An explicit flag overrides the site default."""
    cli = SimpleNamespace(dispatch="slurm-pooled", partition="gpu",
                          nodes=8, walltime="12:00:00")
    choices = cc.resolve_choices(_filled_site(), cli)
    assert choices == {"dispatch": "slurm-pooled", "partition": "gpu",
                       "nodes": 8, "walltime": "12:00:00"}


# ----------------------------------------------------------------
#  build_dispatch_config dispatch on the shape
# ----------------------------------------------------------------

def test_build_dispatch_config_local_is_none():
    """The local opt-out builds no Config at all."""
    site, choices = _filled_site(), {"dispatch": "local"}
    assert cc.build_dispatch_config(site, choices) is None


def test_build_dispatch_config_unknown_shape_raises():
    """An unrecognized dispatch shape is a ConfigError."""
    with pytest.raises(cc.ConfigError):
        cc.build_dispatch_config(_filled_site(),
                                 {"dispatch": "sideways"})


# ----------------------------------------------------------------
#  workers_per_node arithmetic
# ----------------------------------------------------------------

def test_workers_per_node_override_then_derive_then_fallback():
    """An explicit count wins; else cores // cores_per_worker; else one
    worker per node."""
    assert cc.workers_per_node(_filled_site(workers_per_node=8)) == 8
    assert cc.workers_per_node(_filled_site(cores_per_node=32,
                                            cores_per_worker=4)) == 8
    assert cc.workers_per_node(_filled_site(cores_per_node=None)) == 1


# ----------------------------------------------------------------
#  scheduler_options string assembly
# ----------------------------------------------------------------

def test_scheduler_options_emits_mem_gpu_and_extra_but_not_binding():
    """Memory and GPU become #SBATCH directives and the raw passthrough
    rides along; binding does NOT appear -- it is a launch-time concern
    in the deferred parallel path, not a batch directive."""
    site = _filled_site(memory_per_node=192000, gpus_per_node=2,
                        binding="socket",
                        extra_scheduler_options=["#SBATCH --exclusive"])
    # binding alone would trip the serial guard at build time, but
    #   scheduler_options itself must simply never emit a bind line.
    text = cc.scheduler_options(site)
    assert "#SBATCH --mem=192000" in text
    assert "#SBATCH --gres=gpu:2" in text
    assert "#SBATCH --exclusive" in text
    assert "bind" not in text.lower()


def test_scheduler_options_empty_when_nothing_applies():
    """A site with no memory/GPU/extra yields an empty directive
    string."""
    assert cc.scheduler_options(_filled_site()) == ""


# ----------------------------------------------------------------
#  The deferred parallel seam
# ----------------------------------------------------------------

@pytest.mark.parametrize("knob,value", [
    ("launcher", "mpi"), ("ranks_per_worker", 2),
    ("threads_per_rank", 4), ("binding", "socket"),
    ("omp_places", "cores"), ("omp_proc_bind", "spread"),
])
def test_serial_guard_refuses_each_parallel_knob(knob, value):
    """Setting any parallel knob away from its serial default raises a
    clear NotImplementedError naming the deferred seam."""
    site = _filled_site(**{knob: value})
    with pytest.raises(NotImplementedError) as excinfo:
        cc._require_serial_only(site)
    assert knob in str(excinfo.value)


def test_serial_guard_passes_for_a_serial_site():
    """A site at its serial defaults passes the guard silently."""
    assert cc._require_serial_only(_filled_site()) is None


def test_make_launcher_single_then_deferred():
    """The single launcher builds; any other launcher raises the
    deferred-seam error."""
    if HAVE_PARSL:
        from parsl.launchers import SingleNodeLauncher
        assert isinstance(cc.make_launcher(_filled_site()),
                          SingleNodeLauncher)
    with pytest.raises(NotImplementedError):
        cc.make_launcher(_filled_site(launcher="mpi"))


# ----------------------------------------------------------------
#  The two cluster shapes build real Parsl Configs
# ----------------------------------------------------------------

@needs_parsl
def test_per_job_config_block_geometry():
    """slurm-per-job maps each unit to its own one-node, one-worker
    block (init_blocks 0, one worker per node) and threads the site and
    choices through to the provider."""
    site = _filled_site()
    choices = {"dispatch": "slurm-per-job", "partition": "general",
               "nodes": 1, "walltime": "02:00:00"}
    config = cc.build_dispatch_config(site, choices)
    executor = config.executors[0]
    assert executor.label == "imago-per-job"
    assert executor.max_workers_per_node == 1
    assert executor.provider.init_blocks == 0
    assert executor.provider.nodes_per_block == 1
    assert executor.provider.partition == "general"
    assert executor.provider.account == "rulisp-lab"
    assert "module load imago" in executor.provider.worker_init


@needs_parsl
def test_pooled_config_block_geometry():
    """slurm-pooled builds one growable allocation whose workers stream
    units: init_blocks 1, max_blocks from the site, and workers packed
    by the per-node cores."""
    site = _filled_site(max_blocks=4, cores_per_node=32,
                        cores_per_worker=1)
    choices = {"dispatch": "slurm-pooled", "partition": "general",
               "nodes": 2, "walltime": "02:00:00"}
    config = cc.build_dispatch_config(site, choices)
    executor = config.executors[0]
    assert executor.label == "imago-pooled"
    assert executor.max_workers_per_node == 32
    assert executor.provider.init_blocks == 1
    assert executor.provider.max_blocks == 4
    assert executor.provider.nodes_per_block == 2


@needs_parsl
def test_cluster_build_refuses_parallel_knob():
    """Building a cluster Config with a parallel knob set raises through
    the provider's serial guard."""
    site = _filled_site(threads_per_rank=4)
    choices = {"dispatch": "slurm-pooled", "partition": "general",
               "nodes": 1, "walltime": "02:00:00"}
    with pytest.raises(NotImplementedError):
        cc.build_dispatch_config(site, choices)


# ----------------------------------------------------------------
#  The client-side wrappers (run_flight, shared by every client)
# ----------------------------------------------------------------

def test_resolve_dispatch_local_builds_no_config():
    """The local opt-out reads no settings file and produces neither a
    Config nor a resolved-choices record."""
    assert cc.resolve_dispatch("local") == (None, None)


def test_resolve_dispatch_cluster_threads_choices_through(monkeypatch):
    """A cluster shape loads the site, resolves the per-run choices
    against it (the four inputs reaching resolve_choices verbatim), and
    returns the generator's Config alongside those choices."""
    monkeypatch.setattr(cc, "load_site_config",
                        lambda profile: {"_profile": profile})
    captured = {}

    def fake_resolve(site, cli):
        captured["cli"] = cli
        return {"dispatch": cli.dispatch, "partition": cli.partition,
                "nodes": cli.nodes, "walltime": cli.walltime}

    monkeypatch.setattr(cc, "resolve_choices", fake_resolve)
    monkeypatch.setattr(cc, "build_dispatch_config",
                        lambda site, choices: ("CONFIG", choices))

    config, choices = cc.resolve_dispatch(
        "slurm-pooled", "gpu", 4, "06:00:00", "hpc2")
    assert config[0] == "CONFIG"
    assert captured["cli"].dispatch == "slurm-pooled"
    assert choices == {"dispatch": "slurm-pooled", "partition": "gpu",
                       "nodes": 4, "walltime": "06:00:00"}


def test_write_resolved_dispatch_records_the_choices(tmp_path):
    """The reproducible record carries the resolved choices and the
    profile, omits the stable site facts, and round-trips as TOML."""
    choices = {"dispatch": "slurm-per-job", "partition": "general",
               "nodes": 2, "walltime": "02:00:00"}
    cc.write_resolved_dispatch(str(tmp_path), choices, "hpc2")
    with open(tmp_path / "resolved_dispatch.toml", "rb") as handle:
        data = tomllib.load(handle)
    assert data["dispatch"] == "slurm-per-job"
    assert data["partition"] == "general"
    assert data["nodes"] == 2
    assert data["profile"] == "hpc2"
