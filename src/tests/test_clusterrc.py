"""Tests for clusterrc.py -- the per-site cluster dispatch settings
file and its ``--probe`` discovery helper (DESIGN 6.2.11; PSEUDOCODE
13.7).

These are pure-computation tests: the discovery parsers run on
captured-style scheduler output and the live queries are monkeypatched,
so nothing here touches a real scheduler or needs ``$IMAGO_DATA``.  The
import works because conftest.py inserts ``src/scripts`` on ``sys.path``.
"""

import ast

import pytest

import clusterrc

pytestmark = pytest.mark.unit


# ----------------------------------------------------------------
#  The defaults dictionary and the REQUIRED sentinel
# ----------------------------------------------------------------

def test_defaults_have_required_core_and_per_job_default():
    """The shipped defaults leave the two non-discoverable fields as
    the REQUIRED sentinel, default the dispatch shape to slurm-per-job
    (DESIGN 6.2.11, decision 2), and otherwise hold usable values."""
    settings = clusterrc.parameters_and_defaults()
    assert settings["partitions"] is clusterrc.REQUIRED
    assert settings["worker_init"] is clusterrc.REQUIRED
    assert settings["account"] is None
    assert settings["default_topology"] == "slurm-per-job"
    assert settings["cores_per_worker"] == 1
    # Every tier key from the design is present.
    for key in ("launcher", "ranks_per_worker", "threads_per_rank",
                "binding", "omp_places", "omp_proc_bind",
                "gpus_per_node", "queue_overrides", "profiles",
                "extra_scheduler_options", "max_blocks"):
        assert key in settings


def test_required_sentinel_repr_is_obvious():
    """An unfilled key reads clearly in the dump and in tracebacks."""
    assert repr(clusterrc.REQUIRED) == "<REQUIRED>"


# ----------------------------------------------------------------
#  The pure discovery parsers
# ----------------------------------------------------------------

def test_parse_gres_gpu_count_variants():
    """The GPU count is pulled from the several GRES spellings, and a
    GRES-free node reports zero."""
    assert clusterrc.parse_gres_gpu_count("gpu:a100:4") == 4
    assert clusterrc.parse_gres_gpu_count("gpu:4") == 4
    assert clusterrc.parse_gres_gpu_count("gpu:a100:4(S:0-1)") == 4
    assert clusterrc.parse_gres_gpu_count("(null)") == 0
    assert clusterrc.parse_gres_gpu_count("") == 0


def test_parse_sinfo_rows_strips_plus_and_reads_gpus():
    """Each sinfo row yields partition, cores, memory, and GPU count;
    the SLURM ``+`` lower-bound suffix is stripped, and an unparseable
    row is skipped rather than aborting the parse."""
    sinfo_text = ("general 32 192000 (null)\n"
                  "gpu 40+ 384000+ gpu:a100:4\n"
                  "bigmem 64 1536000 (null)\n"
                  "garbage line with no numbers\n")
    rows = clusterrc.parse_sinfo_rows(sinfo_text)
    assert [row["partition"] for row in rows] == \
        ["general", "gpu", "bigmem"]
    assert rows[1]["cores"] == 40           # "40+" -> 40
    assert rows[1]["memory_mb"] == 384000   # "384000+" -> 384000
    assert rows[1]["gpus"] == 4
    assert rows[0]["gpus"] == 0


def test_parse_lscpu_topology_reads_layout():
    """The socket / core / thread / NUMA layout is read; a missing
    label is simply absent from the result."""
    lscpu_text = ("Architecture:            x86_64\n"
                  "CPU(s):                  64\n"
                  "Socket(s):               2\n"
                  "Core(s) per socket:      16\n"
                  "Thread(s) per core:      2\n"
                  "NUMA node(s):            2\n")
    topology = clusterrc.parse_lscpu_topology(lscpu_text)
    assert topology == {"cpus": 64, "sockets": 2,
                        "cores_per_socket": 16,
                        "threads_per_core": 2, "numa_nodes": 2}


def test_parse_sacctmgr_accounts_dedupes_in_order():
    """The account hint lists distinct accounts in first-seen order."""
    text = ("rulisp-lab|general\n"
            "rulisp-lab|gpu\n"
            "other-acct|bigmem\n")
    assert clusterrc.parse_sacctmgr_accounts(text) == \
        ["rulisp-lab", "other-acct"]


# ----------------------------------------------------------------
#  probe_site() orchestration and graceful degradation
# ----------------------------------------------------------------

def test_probe_site_degrades_when_no_tools(monkeypatch):
    """With every query unavailable, the probe returns an empty dict
    rather than failing -- the best-effort 'fill what I could' rule."""
    monkeypatch.setattr(clusterrc, "run_query", lambda command: None)
    assert clusterrc.probe_site() == {}


def test_probe_site_assembles_discovered_facts(monkeypatch):
    """A working scheduler/host yields the partition list, the
    representative per-node facts, the topology, and the account
    hint."""
    def fake_query(command):
        tool = command[0]
        if tool == "sinfo":
            return ("general 32 192000 (null)\n"
                    "gpu 40 384000 gpu:a100:4\n")
        if tool == "lscpu":
            return "Socket(s):  2\nNUMA node(s):  2\n"
        if tool == "sacctmgr":
            return "rulisp-lab|general\n"
        return None

    monkeypatch.setattr(clusterrc, "run_query", fake_query)
    monkeypatch.setenv("USER", "rulisp")
    facts = clusterrc.probe_site()
    assert facts["partitions"] == ["general", "gpu"]
    assert facts["cores_per_node"] == 32      # first row, representative
    assert facts["memory_per_node"] == 192000
    assert facts["gpus_per_node"] == 4         # max across queues
    assert facts["topology"]["sockets"] == 2
    assert facts["accounts"] == ["rulisp-lab"]


# ----------------------------------------------------------------
#  Rendering and writing a starter clusterrc.py
# ----------------------------------------------------------------

def _starter_settings(text):
    """Compile a rendered starter file and return its settings dict."""
    ast.parse(text)                       # must be valid Python
    namespace = {}
    exec(compile(text, "<starter>", "exec"), namespace)
    return namespace["parameters_and_defaults"]()


def test_render_starter_fills_discovered_and_blanks_the_rest():
    """A discovered fact set renders a valid file that fills the
    partition list and per-node facts, while leaving worker_init and
    account as blanks the user must complete."""
    facts = {"partitions": ["general", "gpu"], "cores_per_node": 32,
             "memory_per_node": 192000, "gpus_per_node": 4,
             "topology": {"sockets": 2, "numa_nodes": 2},
             "accounts": ["rulisp-lab"]}
    text = clusterrc.render_starter_clusterrc(facts)
    settings = _starter_settings(text)
    assert settings["partitions"] == ["general", "gpu"]
    assert settings["cores_per_node"] == 32
    assert settings["gpus_per_node"] == 4
    # The non-discoverable fields are present but empty (so the
    #   generator's required-check refuses an unfinished file).
    assert settings["worker_init"] == []
    assert settings["account"] is None
    # The topology and account hints ride along as comments.
    assert "Discovered CPU topology" in text
    assert "rulisp-lab" in text


def test_render_starter_with_no_facts_is_valid_with_blank_partitions():
    """Even when the probe found nothing, the rendered file is valid
    Python with an empty partition list awaiting the user."""
    text = clusterrc.render_starter_clusterrc({})
    settings = _starter_settings(text)
    assert settings["partitions"] == []
    assert settings["worker_init"] == []


def test_write_starter_writes_then_refuses_to_clobber(tmp_path):
    """write_starter_clusterrc writes the file once; a second call
    without overwrite keeps the existing file, and --force replaces
    it."""
    target = tmp_path / "clusterrc.py"
    facts = {"partitions": ["general"]}

    assert clusterrc.write_starter_clusterrc(
        str(target), discovered_facts=facts) is True
    assert target.exists()
    first_text = target.read_text()

    # Second call without overwrite must not clobber.
    assert clusterrc.write_starter_clusterrc(
        str(target), discovered_facts={"partitions": ["other"]}) \
        is False
    assert target.read_text() == first_text

    # With overwrite the file is replaced.
    assert clusterrc.write_starter_clusterrc(
        str(target), discovered_facts={"partitions": ["other"]},
        overwrite=True) is True
    assert _starter_settings(target.read_text())["partitions"] == \
        ["other"]


# ----------------------------------------------------------------
#  The command-line entry point
# ----------------------------------------------------------------

def test_main_default_dumps_settings(capsys):
    """With no arguments the file prints its defaults (the *rc.py
    convention) and exits zero."""
    assert clusterrc.main([]) == 0
    captured = capsys.readouterr().out
    assert "partitions" in captured
    assert "default_topology" in captured


def test_main_probe_writes_starter(tmp_path, monkeypatch):
    """``--probe`` probes the machine and writes a starter file at the
    requested path, returning zero on success."""
    monkeypatch.setattr(clusterrc, "probe_site",
                        lambda: {"partitions": ["general"]})
    target = tmp_path / "out_clusterrc.py"
    assert clusterrc.main(["--probe", "-o", str(target)]) == 0
    assert target.exists()
    assert _starter_settings(target.read_text())["partitions"] == \
        ["general"]
