"""Tests for clusterrc.py -- the per-site cluster dispatch settings
file (DESIGN 6.2.11; PSEUDOCODE 13.7).

``clusterrc.py`` is pure data, like every other ``*rc.py``: a single
``parameters_and_defaults()`` returning one dictionary, with a trivial
``__main__`` that dumps it.  The discovery tool that *generates* a
starter copy lives separately in ``cluster_probe.py`` (tested in
``test_cluster_probe.py``).  The import works because conftest.py
inserts ``src/scripts`` on ``sys.path``.
"""

import pytest

import clusterrc

pytestmark = pytest.mark.unit


def test_required_core_ships_unfilled():
    """The two non-discoverable required fields ship as None, so the
    dispatch generator refuses to run until the user supplies them."""
    settings = clusterrc.parameters_and_defaults()
    assert settings["partitions"] is None
    assert settings["worker_init"] is None
    # account is optional (only some clusters require one).
    assert settings["account"] is None


def test_per_job_is_the_default_topology():
    """The default dispatch shape is slurm-per-job (DESIGN 6.2.11,
    decision 2)."""
    assert clusterrc.parameters_and_defaults()["default_topology"] == \
        "slurm-per-job"


def test_every_tier_key_is_present_with_a_usable_default():
    """Beyond the required core, every performance and advanced key is
    present and already holds a usable default."""
    settings = clusterrc.parameters_and_defaults()
    expected_defaults = {
        "cores_per_worker": 1, "nodes": 1, "walltime": "01:00:00",
        "max_blocks": 1, "launcher": "single", "ranks_per_worker": 1,
        "threads_per_rank": 1, "gpus_per_node": 0,
        "queue_overrides": {}, "profiles": {},
        "extra_scheduler_options": [],
    }
    for key, value in expected_defaults.items():
        assert settings[key] == value
    # The optional-but-unset knobs are present as None.
    for key in ("cores_per_node", "workers_per_node",
                "memory_per_node", "memory_per_worker", "binding",
                "omp_places", "omp_proc_bind"):
        assert settings[key] is None


def test_module_is_pure_data_no_cli_machinery():
    """The file is data only: it exposes parameters_and_defaults but
    none of the probe/CLI machinery, which now lives in
    cluster_probe.py."""
    assert hasattr(clusterrc, "parameters_and_defaults")
    for moved in ("probe_site", "render_starter_clusterrc",
                  "write_starter_clusterrc", "main", "REQUIRED"):
        assert not hasattr(clusterrc, moved)
