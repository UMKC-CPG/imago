"""kaleidoscope -- the high-throughput Imago campaign runner
(VISION Goal 4; ARCHITECTURE 9.4, 9.6; DESIGN 6.2; PSEUDOCODE
13).

kaleidoscope drives a *set* of Imago calculations: it dispatches
the per-structure work, tracks each one's outcome, caches
completed runs so a campaign can be resumed by simply re-running
it, and surfaces a report the client uses to decide acceptance.
Per VISION Principle 9 it is ordinary, domain-agnostic
scientific Python -- it never interprets what a run computed.
The Imago-specific behavior lives in the runner
(``kaleidoscope.runners``) below it, and in the client that
reads each run directory after the campaign (the harvest) above
it.

Typical use::

    from kaleidoscope import Campaign, CalcUnit, run_campaign

    campaign = Campaign(root="/work/my_campaign", units=[...])
    report = run_campaign(campaign)
    for entry in report.with_detail("converged"):
        ...   # client harvests entry.run_dir/result.toml

This ``__init__`` is the package's public façade; the supporting
modules (``model``, ``workspace``, ``cache``, ``runners``,
``dispatch``) carry the implementation.
"""

from .model import (
    KaleidoscopeError,
    KeyFile,
    KeyFields,
    CalcUnit,
    Campaign,
    RunOutcome,
    ReportEntry,
    CampaignReport,
)
from .runners import (
    Runner,
    ImagoRunner,
    register_runner,
    resolve_runner,
    RUNNERS,
)
from .dispatch import (
    run_campaign,
    LocalExecutor,
    ParslExecutor,
    TaskLost,
    report_entry_from_status,
)
from .workspace import (
    unit_run_dir,
    validate_campaign,
    read_status,
    write_status,
)
from .cache import (
    is_cache_hit,
    cache_key_matches,
    write_cache_key,
)

__all__ = [
    "KaleidoscopeError", "KeyFile", "KeyFields", "CalcUnit",
    "Campaign", "RunOutcome", "ReportEntry", "CampaignReport",
    "Runner", "ImagoRunner", "register_runner", "resolve_runner",
    "RUNNERS", "run_campaign", "LocalExecutor", "ParslExecutor",
    "TaskLost", "report_entry_from_status", "unit_run_dir",
    "validate_campaign", "read_status", "write_status",
    "is_cache_hit", "cache_key_matches", "write_cache_key",
]
