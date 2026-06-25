"""Tests for cluster_probe.py -- the discovery tool that generates a
starter ``clusterrc.py`` (DESIGN 6.2.11; PSEUDOCODE 13.7).

These are pure-computation tests: the parsers run on captured-style
scheduler output and the live queries are monkeypatched, so nothing
here touches a real scheduler.  The import works because conftest.py
inserts ``src/scripts`` on ``sys.path``.
"""

import ast

import pytest

import cluster_probe

pytestmark = pytest.mark.unit


# ----------------------------------------------------------------
#  The pure discovery parsers
# ----------------------------------------------------------------

def test_parse_gres_gpu_count_variants():
    """The GPU count is pulled from the several GRES spellings, and a
    GRES-free node reports zero."""
    assert cluster_probe.parse_gres_gpu_count("gpu:a100:4") == 4
    assert cluster_probe.parse_gres_gpu_count("gpu:4") == 4
    assert cluster_probe.parse_gres_gpu_count("gpu:a100:4(S:0-1)") == 4
    assert cluster_probe.parse_gres_gpu_count("(null)") == 0
    assert cluster_probe.parse_gres_gpu_count("") == 0


def test_parse_sinfo_rows_strips_plus_and_reads_gpus():
    """Each sinfo row yields partition, cores, memory, and GPU count;
    the SLURM ``+`` lower-bound suffix is stripped, and an unparseable
    row is skipped rather than aborting the parse."""
    sinfo_text = ("general 32 192000 (null)\n"
                  "gpu 40+ 384000+ gpu:a100:4\n"
                  "bigmem 64 1536000 (null)\n"
                  "garbage line with no numbers\n")
    rows = cluster_probe.parse_sinfo_rows(sinfo_text)
    assert [row["partition"] for row in rows] == \
        ["general", "gpu", "bigmem"]
    assert rows[1]["cores"] == 40           # "40+" -> 40
    assert rows[1]["memory_mb"] == 384000   # "384000+" -> 384000
    assert rows[1]["gpus"] == 4
    assert rows[0]["gpus"] == 0


def test_parse_sacctmgr_accounts_dedupes_in_order():
    """The account hint lists distinct accounts in first-seen order."""
    text = ("rulisp-lab|general\n"
            "rulisp-lab|gpu\n"
            "other-acct|bigmem\n")
    assert cluster_probe.parse_sacctmgr_accounts(text) == \
        ["rulisp-lab", "other-acct"]


# ----------------------------------------------------------------
#  probe_site() orchestration and graceful degradation
# ----------------------------------------------------------------

def test_probe_site_degrades_when_no_tools(monkeypatch):
    """With every query unavailable, the probe returns an empty dict
    rather than failing -- the best-effort 'fill what I could' rule."""
    monkeypatch.setattr(cluster_probe, "run_query",
                        lambda command: None)
    assert cluster_probe.probe_site() == {}


def test_probe_site_homogeneous_fills_per_node_values(monkeypatch):
    """When every node agrees on cores/memory/GPUs, the per-node values
    are filled in directly, alongside the queue list and account
    hint."""
    def fake_query(command):
        tool = command[0]
        if tool == "sinfo":
            return ("general 64 384000 (null)\n"
                    "debug 64 384000 (null)\n")
        if tool == "sacctmgr":
            return "rulisp-lab|general\n"
        return None

    monkeypatch.setattr(cluster_probe, "run_query", fake_query)
    monkeypatch.setenv("USER", "rulisp")
    facts = cluster_probe.probe_site()
    assert facts["partitions"] == ["general", "debug"]
    assert facts["cores_per_node"] == 64
    assert facts["memory_per_node"] == 384000
    assert facts["gpus_per_node"] == 0
    assert facts["accounts"] == ["rulisp-lab"]
    # No "options" keys when the nodes agree.
    assert "core_options" not in facts


def test_probe_site_heterogeneous_records_options_not_a_guess(
        monkeypatch):
    """When the nodes disagree, the per-node value is NOT guessed: the
    distinct values seen are recorded under *_options and the setting
    itself is left unset for the user."""
    def fake_query(command):
        if command[0] == "sinfo":
            return ("general 32 192000 (null)\n"
                    "bigmem 64 1536000 (null)\n")
        return None

    monkeypatch.setattr(cluster_probe, "run_query", fake_query)
    facts = cluster_probe.probe_site()
    assert facts["partitions"] == ["general", "bigmem"]
    # Cores and memory disagree -> options recorded, value left unset.
    assert "cores_per_node" not in facts
    assert facts["core_options"] == [32, 64]
    assert "memory_per_node" not in facts
    assert facts["mem_options"] == [192000, 1536000]
    # GPUs agree (both zero) -> filled, no options.
    assert facts["gpus_per_node"] == 0


def test_probe_site_reads_only_the_scheduler_never_the_login_node(
        monkeypatch):
    """The probe runs only scheduler queries (sinfo, sacctmgr) -- never
    lscpu/numactl, which would describe the login node, the wrong
    machine."""
    seen_tools = []

    def fake_query(command):
        seen_tools.append(command[0])
        return None

    monkeypatch.setattr(cluster_probe, "run_query", fake_query)
    monkeypatch.setenv("USER", "rulisp")
    cluster_probe.probe_site()
    assert "lscpu" not in seen_tools
    assert "numactl" not in seen_tools
    assert set(seen_tools) <= {"sinfo", "sacctmgr"}


# ----------------------------------------------------------------
#  Rendering and writing a starter clusterrc.py
# ----------------------------------------------------------------

def _starter_settings(text):
    """Compile a rendered starter file and return its settings dict."""
    ast.parse(text)                       # must be valid Python
    namespace = {}
    exec(compile(text, "<starter>", "exec"), namespace)
    return namespace["parameters_and_defaults"]()


def test_starter_schema_matches_clusterrc():
    """cluster_probe is self-contained -- it carries its own copy of the
    schema and never reads clusterrc.py -- so this test is the guard
    that keeps the two from drifting: same keys, same order, same
    defaults."""
    import clusterrc
    canonical = clusterrc.parameters_and_defaults()
    schema = cluster_probe._starter_schema()
    assert schema == canonical
    assert list(schema) == list(canonical)   # same insertion order too


def test_render_starter_offers_the_full_schema():
    """The starter offers every key the tool's own schema defines,
    fills the discovered facts, carries a plain-language note on each
    setting, and leaves the required blanks flagged FILL IN."""
    schema_keys = set(cluster_probe._starter_schema())

    facts = {"partitions": ["general", "gpu"], "cores_per_node": 32,
             "memory_per_node": 192000, "gpus_per_node": 4,
             "accounts": ["rulisp-lab"]}
    text = cluster_probe.render_starter_clusterrc(facts)
    settings = _starter_settings(text)
    # Full key parity with the tool's schema.
    assert set(settings) == schema_keys
    # Discovered facts filled.
    assert settings["partitions"] == ["general", "gpu"]
    assert settings["cores_per_node"] == 32
    assert settings["gpus_per_node"] == 4
    # The non-discoverable required field stays blank and is flagged.
    assert settings["worker_init"] is None
    assert "# FILL IN" in text
    # Every setting carries a plain-language note (spot-check two).
    assert "# How many CPU cores one node has." in text
    assert "scheduler queues you can submit to" in text
    # The account hint rides along on the account line.
    assert "rulisp-lab" in text


def test_render_starter_writes_no_login_node_facts():
    """No login-node CPU topology (or any 'discovered topology' note)
    is ever written into the file."""
    text = cluster_probe.render_starter_clusterrc(
        {"partitions": ["general"], "cores_per_node": 64})
    lowered = text.lower()
    # The login-node hardware markers must be absent.  ("topology" by
    #   itself is not checked -- it legitimately appears in the
    #   'default_topology' setting key.)
    assert "cpu topology" not in lowered
    assert "lscpu" not in lowered
    assert "socket" not in lowered
    assert "numa" not in lowered


def test_render_starter_lists_heterogeneous_options_not_a_guess():
    """A per-node value the cluster disagreed on is left blank, flagged
    FILL IN, with the values seen listed -- never silently guessed."""
    facts = {"partitions": ["general", "bigmem"],
             "core_options": [32, 64],
             "mem_options": [192000, 1536000]}
    text = cluster_probe.render_starter_clusterrc(facts)
    settings = _starter_settings(text)
    # The disagreed-on settings are left unset (their schema default).
    assert settings["cores_per_node"] is None
    assert settings["memory_per_node"] is None
    # ...and the file tells the user what was seen.
    assert "Nodes vary -- core counts seen: 32, 64." in text
    assert "memory sizes (MB) seen: 192000, 1536000." in text


def test_render_starter_with_no_facts_blanks_partitions(monkeypatch):
    """When the probe found nothing, the rendered file is still valid
    Python with both required fields blank and flagged."""
    text = cluster_probe.render_starter_clusterrc({})
    settings = _starter_settings(text)
    assert settings["partitions"] is None
    assert settings["worker_init"] is None
    assert text.count("# FILL IN") == 2


def test_write_starter_writes_then_refuses_to_clobber(tmp_path):
    """write_starter_clusterrc writes the file once; a second call
    without overwrite keeps the existing file, and force replaces
    it."""
    target = tmp_path / "clusterrc.py"
    facts = {"partitions": ["general"]}

    assert cluster_probe.write_starter_clusterrc(
        str(target), discovered_facts=facts) is True
    assert target.exists()
    first_text = target.read_text()

    # Second call without overwrite must not clobber.
    assert cluster_probe.write_starter_clusterrc(
        str(target), discovered_facts={"partitions": ["other"]}) \
        is False
    assert target.read_text() == first_text

    # With overwrite the file is replaced.
    assert cluster_probe.write_starter_clusterrc(
        str(target), discovered_facts={"partitions": ["other"]},
        overwrite=True) is True
    assert _starter_settings(target.read_text())["partitions"] == \
        ["other"]


# ----------------------------------------------------------------
#  The command-line entry point
# ----------------------------------------------------------------

def test_main_probes_and_writes_starter(tmp_path, monkeypatch):
    """The CLI probes the machine and writes a starter at the requested
    path, returning zero on success."""
    monkeypatch.setattr(cluster_probe, "probe_site",
                        lambda: {"partitions": ["general"]})
    target = tmp_path / "out_clusterrc.py"
    assert cluster_probe.main(["-o", str(target)]) == 0
    assert target.exists()
    assert _starter_settings(target.read_text())["partitions"] == \
        ["general"]
