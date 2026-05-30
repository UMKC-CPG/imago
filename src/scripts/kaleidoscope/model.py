"""kaleidoscope.model -- the flight data model
(DESIGN 6.2.1, 6.2.2, 6.2.6; PSEUDOCODE 13.1, 13.2, 13.6).

These are plain, domain-agnostic records.  Per VISION
Principle 9 the flight layer is ordinary scientific Python:
kaleidoscope dispatches, tracks, and caches calculations
without knowing what any of them computes.  The Imago-specific
meaning of a run lives one level down in the wingbeat
(``wingbeats.py``) and one level up in the client that reads each
run directory after the flight (the harvest, DESIGN 6.2.6).
"""

from dataclasses import dataclass, field
from typing import Any, Callable, Optional


class KaleidoscopeError(Exception):
    """A flight-construction or workspace fault -- for
    example an ``id`` that is not a filesystem-safe slug, two
    units colliding on one run directory, or a request for an
    unregistered wingbeat.  Distinct from ``imago.ImagoError``,
    which concerns the running of a single calculation."""
    pass


@dataclass
class KeyFile:
    """One file that contributes to a unit's cache identity
    (DESIGN 6.2.5).  ``name`` is the path, relative to the run
    directory, of the staged copy left by a prior run;
    ``source`` is the absolute path of the current input the
    cache compares against it, byte-for-byte.  Naming both
    explicitly keeps kaleidoscope from guessing how a client's
    inputs map onto staged files."""
    name: str
    source: str


@dataclass
class KeyFields:
    """The client-declared cache identity of a calculation
    (DESIGN 6.2.5).  ``scalars`` are compared verbatim,
    field-by-field (e.g. a k-point spec, a convergence
    threshold, an Imago commit).  ``files`` are compared by
    byte-comparison against their staged copies.  Only the
    client knows which inputs define identity for its
    calculations, so it supplies these; kaleidoscope never
    guesses (a too-broad key risks false hits and wrong
    science, a too-narrow key risks needless re-runs)."""
    scalars: dict = field(default_factory=dict)
    files: list = field(default_factory=list)   # list[KeyFile]


@dataclass
class CalcUnit:
    """One unit of work: a single calculation to run
    (DESIGN 6.2.1).

    - ``id``        : the stable per-structure key (a COD id or
                      a curation reference_id); must be a slug.
    - ``structure`` : path to an ``imago.skl`` for the run.
    - ``options``   : the makeinput options used to build the
                      run directory (passed to the wingbeat).
    - ``calc``      : the per-axis directory components for this
                      unit -- a tuple carrying one
                      ``<axis>-<value>`` slug per varied sweep axis
                      (DESIGN 6.2.1).  The empty tuple ``()`` means
                      "no second directory level", which is the
                      common single-calculation case; the
                      predict-then-verify helper (DESIGN 6.2.8)
                      fills it with one component per swept knob.
    - ``wingbeat``    : the wingbeat name, or None for the flight
                      default.
    - ``key_fields``: the cache identity (DESIGN 6.2.5).
    """
    id: str
    structure: str
    options: dict = field(default_factory=dict)
    calc: tuple = ()
    wingbeat: Optional[str] = None
    key_fields: KeyFields = field(default_factory=KeyFields)


@dataclass
class SweepRecord:
    """Description of the option-axis sweep a flight performs
    (DESIGN 6.2.1).  It is attached only when the flight was
    built by the predict-then-verify helper (DESIGN 6.2.8);
    a flight assembled by hand leaves ``Flight.sweep`` None.

    - ``varied_axes`` : the axis names that change from one unit
                        to the next, in the SAME order they appear
                        as directory components in each
                        ``CalcUnit.calc`` tuple.  In v1 this is the
                        single ``("kpt-density",)``.
    - ``fixed_axes``  : a mapping from each axis held constant
                        across the whole flight to its single
                        value (e.g. ``basis``, ``functional``).

    The harvest step (DESIGN 7.8) reads this to learn what was
    swept without having to reverse-engineer run-directory paths.
    """
    varied_axes: tuple = ()
    fixed_axes: dict = field(default_factory=dict)


@dataclass
class Flight:
    """A set of calculations plus global options
    (DESIGN 6.2.1).

    - ``root``           : the workspace root directory.
    - ``units``          : the list of CalcUnit to run.
    - ``default_wingbeat`` : the wingbeat for units naming none.
    - ``parsl_config``   : an optional Parsl ``Config``.  When
                           present the flight dispatches
                           through Parsl (cluster parallelism);
                           when None it runs locally and
                           serially.  Keeping Parsl behind this
                           optional field keeps it central
                           (VISION Goal 4) without making it a
                           hard import-time dependency.
    - ``on_outcome``     : an optional callback invoked with each
                           unit's ReportEntry as it reaches a
                           terminal state, for streaming
                           progress.
    - ``sweep``          : a SweepRecord describing the varied and
                           fixed axes, set only when the flight was
                           built by the predict-then-verify helper
                           (DESIGN 6.2.8); None for a hand-built
                           flight.  When present, serialize_flight
                           emits a ``[flight.sweep]`` block.
    - ``metadata``       : opaque per-key tables that
                           serialize_flight round-trips VERBATIM as
                           ``[flight.<key>]`` blocks and the
                           dispatch core never reads (Principle 9).
                           This is the seam by which the 6.2.8
                           helper stashes its PredictionRecord
                           (under key ``"prediction"``) without the
                           domain-agnostic core having to interpret
                           it.
    """
    root: str
    units: list = field(default_factory=list)
    default_wingbeat: str = "imago"
    parsl_config: Any = None
    on_outcome: Optional[Callable] = None
    sweep: Optional[SweepRecord] = None
    metadata: dict = field(default_factory=dict)


@dataclass
class WingbeatOutcome:
    """The domain-agnostic result a wingbeat returns
    (DESIGN 6.2.2).  ``ok`` says whether the unit *completed*,
    not whether it "succeeded scientifically".  ``detail`` is an
    opaque string the wingbeat chooses (e.g. "converged",
    "not_converged") that kaleidoscope records verbatim and
    never interprets -- which is what lets the flight layer
    surface convergence in ``status.toml`` while staying
    ignorant of what convergence means."""
    ok: bool
    detail: str = ""
    runtime_seconds: float = 0.0
    message: str = ""


@dataclass
class ReportEntry:
    """One unit's terminal record in the flight report
    (DESIGN 6.2.6).  Mirrors the generic ``status.toml`` fields
    and carries nothing domain-specific; the client reads each
    ``wingbeat_dir`` itself to harvest results.  ``calc`` carries
    the unit's per-axis directory tuple verbatim (``()`` when the
    unit has no second directory level)."""
    id: str
    calc: tuple
    status: str
    detail: Optional[str]
    wingbeat_dir: str
    runtime_seconds: Optional[float]
    message: Optional[str]


@dataclass
class FlightReport:
    """The result of ``dispatch``: one ReportEntry per unit,
    in unit order, plus convenience views (DESIGN 6.2.6).
    kaleidoscope never decides whether the aggregate is
    acceptable -- that judgment belongs to the client (VISION
    Principle 10)."""
    entries: list = field(default_factory=list)

    def by_status(self, status):
        """Entries whose lifecycle status equals ``status``
        (one of queued/running/done/failed/lost)."""
        return [e for e in self.entries if e.status == status]

    def with_detail(self, detail):
        """Entries whose wingbeat-supplied ``detail`` equals
        ``detail`` (e.g. "converged").  This is how a client
        such as the potential-DB producer selects the units it
        will harvest."""
        return [e for e in self.entries if e.detail == detail]

    def failures(self):
        """Entries that did not complete cleanly (status
        ``failed`` or ``lost``)."""
        return [e for e in self.entries
                if e.status in ("failed", "lost")]
