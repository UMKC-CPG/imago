"""kaleidoscope -- the high-throughput Imago campaign dispatcher
(VISION Goal 4; ARCHITECTURE 9.4, 9.6; DESIGN 6.2; PSEUDOCODE
13).

kaleidoscope drives a *set* of Imago calculations: it dispatches
the per-structure work, tracks each one's outcome, caches
completed runs so a campaign can be resumed by simply re-running
it, and surfaces a report the client uses to decide acceptance.
Per VISION Principle 9 it is ordinary, domain-agnostic
scientific Python -- it never interprets what a run computed.
The Imago-specific behavior lives in the wingbeat
(``kaleidoscope.wingbeats``) below it, and in the client that
reads each run directory after the campaign (the harvest) above
it.

Typical use::

    from kaleidoscope import Campaign, CalcUnit, dispatch

    campaign = Campaign(root="/work/my_campaign", units=[...])
    report = dispatch(campaign)
    for entry in report.with_detail("converged"):
        ...   # client harvests entry.wingbeat_dir/result.toml

This ``__init__`` is the package's public façade; the supporting
modules (``model``, ``workspace``, ``cache``, ``wingbeats``,
``dispatch``) carry the implementation.
"""

from .model import (
    KaleidoscopeError,
    KeyFile,
    KeyFields,
    CalcUnit,
    Campaign,
    WingbeatOutcome,
    ReportEntry,
    CampaignReport,
)
from .wingbeats import (
    Wingbeat,
    ImagoWingbeat,
    register_wingbeat,
    resolve_wingbeat,
    WINGBEATS,
)
from .dispatch import (
    dispatch,
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
    "Campaign", "WingbeatOutcome", "ReportEntry", "CampaignReport",
    "Wingbeat", "ImagoWingbeat", "register_wingbeat", "resolve_wingbeat",
    "WINGBEATS", "dispatch", "LocalExecutor", "ParslExecutor",
    "TaskLost", "report_entry_from_status", "unit_run_dir",
    "validate_campaign", "read_status", "write_status",
    "is_cache_hit", "cache_key_matches", "write_cache_key",
]
