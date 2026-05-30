"""kaleidoscope.builders.predict_verify -- the predict-then-verify
k-point-density flight builder (DESIGN 6.2.8 / 7.7;
PSEUDOCODE 15.6).

This is the FIRST option-axis builder in the flight-builder split
(DESIGN 6.2; see ``kaleidoscope.builders``).  It turns "a
structure + a dict of fixed makeinput options + a loaded
historical-guidance ``Dataspace``" into a ``Flight`` of
``CalcUnit``s laid out as a *verification grid* around the
predictor's predicted operating point, plus a
``PredictionRecord`` the later harvest step recovers.

It is one *strategy* for one *axis*: it sweeps the k-point density
and does so by PREDICTING a converged value from the guidance
database and then VERIFYING around it.  Sibling builders for other
axes (an XANES target-atom sweep, a basis-size sweep) are plain
enumerations and live beside this module; structure-axis sweeps
are out of kaleidoscope's scope entirely (DESIGN 6.2).

**Why this is a builder, not part of the dispatch core.**  The
core (``model`` / ``workspace`` / ``cache`` / ``dispatch``) is
deliberately domain-agnostic: it runs and tracks calculations
without knowing what they compute (VISION Principle 9/12).  This
builder is the opposite -- it consults the guidance database,
knows the swept knob is a k-point density, and knows how to spell
the per-grid-point directory tag.  To keep the dumb core's import
graph free of the physics layer, a client imports it explicitly::

    from kaleidoscope.builders.predict_verify import predict_settings

Importing it pulls in ``guidance_db`` and ``structure_control``;
the core never does.

The verification-grid idea (DESIGN 7.7): rather than make the
user guess a converged k-point density up front, look up what
worked for chemically-similar systems, predict a density, and
run a small spread of calculations *around* that prediction --
tight when the predictor is confident, wide when it is not -- so
the truly converged point is captured even if the guess was off.
"""

import math
import os
import re
from dataclasses import asdict, dataclass

from ..model import (
    CalcUnit,
    Flight,
    KaleidoscopeError,
    KeyFields,
    KeyFile,
    SweepRecord,
)

# The physics layer.  These are top-level scripts installed flat
#   under bin/ (on a script's sys.path), so an absolute import is
#   correct -- they are siblings of the kaleidoscope package, not
#   members of it.
from guidance_db import compute_signature, predict
from structure_control import StructureControl


# ------------------------------------------------------------------
#  The prediction record persisted alongside the flight
# ------------------------------------------------------------------

@dataclass
class PredictionRecord:
    """Everything about the prediction that drove the grid, kept
    so the harvest step (DESIGN 7.8) can recover the confidence
    score and the neighbor entries without re-running the
    predictor (PSEUDOCODE 15.6).  ``serialize_flight`` writes it
    verbatim as the ``[flight.prediction]`` table.

    - ``policy``             : which grid path was taken
                               (``trust_no_verify`` /
                               ``wide_grid_no_prior`` /
                               ``verify_around_prediction``).
    - ``predicted_value``    : the predicted converged k-density.
    - ``confidence``         : the predictor's combined confidence
                               in [0, 1] that drove the grid width.
    - ``is_under_trained``   : True when the dataspace was too thin
                               for a trustworthy prediction.
    - ``neighbor_entry_ids`` : the guidance entries the prediction
                               drew on (provenance).
    - ``predicted_gap`` /
      ``predicted_spin_pol`` : the intermediate stage-1 quantities
                               (None for non-crystalline systems).
    - ``system_type``        : the declared system type.
    - ``feature_vector``     : the query Signature the prediction
                               was made for.
    """
    policy: str
    predicted_value: float | None
    confidence: float
    is_under_trained: bool
    neighbor_entry_ids: tuple = ()
    predicted_gap: float | None = None
    predicted_spin_pol: float | None = None
    system_type: str = ""
    feature_vector: object = None


# ------------------------------------------------------------------
#  Grid construction (DESIGN 7.7 / 7.9)
# ------------------------------------------------------------------

def default_wide_kpoint_density_grid():
    """The fixed 8-point bracket spanning a factor of 16, used
    when no usable predictor exists (DESIGN 7.9).  It lives here
    rather than in the dataspace because an empty crystalline
    subtree has no content to consult -- the chicken-and-egg the
    seed flight (C75) is meant to break."""
    return [25.0, 50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 400.0]


def logspace(low, high, num_points):
    """Return ``num_points`` geometrically (log-) spaced values
    from ``low`` to ``high`` inclusive.  Log spacing gives
    symmetric multiplicative coverage around a center, which is
    what the verification grid wants.  A request for a single
    point returns the geometric midpoint."""
    if num_points <= 1:
        return [math.sqrt(low * high)]
    log_step = (math.log(high) - math.log(low)) / (num_points - 1)
    return [math.exp(math.log(low) + index * log_step)
            for index in range(num_points)]


def build_verification_grid(center, confidence):
    """Lay out the k-density values to sweep around ``center``,
    with both the multiplicative span and the number of points
    scaling inversely with predictor ``confidence`` (DESIGN 7.7).

    At ``confidence = 1`` the grid is a tight 3 points spanning
    ``[center/1.2, center*1.2]`` -- a quick check that the prior
    still holds.  At ``confidence = 0`` it widens to ~7 points
    spanning ``[center/2.7, center*2.7]``, broad enough to bracket
    the true converged point even if the prediction was a poor
    guide.  The constants are starting heuristics (DESIGN 7.7);
    keeping them in this one function makes post-seed calibration
    a single-file change."""
    # Multiplicative half-span: 1.2 means center/1.2 .. center*1.2.
    width = 1.2 + 1.5 * (1.0 - confidence)
    # Point count grows as confidence falls.
    num_points = round(3 + 4 * (1.0 - confidence))
    return logspace(center / width, center * width, num_points)


# ------------------------------------------------------------------
#  Calc-tag encoding (DESIGN 6.2.4)
# ------------------------------------------------------------------

_AXIS_NAME_RE = re.compile(r"^[a-z0-9-]+$")


def encode_axis_value(value):
    """Encode one swept axis value as a slug-safe token
    (DESIGN 6.2.4 rule 3).  An integer-valued float renders as a
    plain decimal integer; otherwise the compact decimal has its
    ``.`` rewritten to ``p`` and a leading ``-`` to ``m`` so the
    result is always ``[a-z0-9-]`` and parses back unambiguously.
    In v1 the k-density is rounded to an integer before tagging,
    so this is just the decimal integer; the general encoder is
    kept for future non-integer axes."""
    if value == round(value):
        text = str(int(round(value)))
    else:
        # Compact decimal: trim trailing zeros and any bare point.
        text = f"{value:.6f}".rstrip("0").rstrip(".")
    negative = text.startswith("-")
    if negative:
        text = text[1:]
    text = text.replace(".", "p")
    return ("m" + text) if negative else text


def build_calc_tag(calc_axes):
    """Turn an ORDERED ``{axis: value}`` mapping into the unit's
    ``calc`` tuple: one ``"<axis>-<encoded-value>"`` directory
    component per varied axis, in mapping order (DESIGN
    6.2.1/6.2.4).  The order must match
    ``SweepRecord.varied_axes``.  In v1 this is the single
    component ``("kpt-density-<int>",)``.  Each axis name must
    itself be a slug so it is safe as a directory level."""
    components = []
    for axis, value in calc_axes.items():
        if not _AXIS_NAME_RE.match(axis):
            raise KaleidoscopeError(
                f"calc axis name is not a slug "
                f"([a-z0-9-]+): {axis!r}")
        components.append(f"{axis}-{encode_axis_value(value)}")
    return tuple(components)


# ------------------------------------------------------------------
#  Cache identity and structure handling
# ------------------------------------------------------------------

# The producer's cache identity for a converged-potential run
#   (DESIGN 6.2.1): the scalar settings that define "the same
#   calculation" plus the structure file, byte-compared.
_KEY_SCALAR_NAMES = ("scf_threshold", "imago_commit")


def standard_key_fields(structure, options):
    """Build the ``KeyFields`` cache identity for a unit
    (DESIGN 6.2.1/6.2.5): the structure file byte-compared, plus
    the scalar settings that define run identity.

    The scalar names are taken from ``options`` when present
    (DESIGN-6.2.1 lists ``scf_threshold`` and ``imago_commit``).
    ``scf_threshold`` is naturally in the makeinput options; the
    build-identity ``imago_commit`` is producer-injected -- it is
    carried in ``options`` when the producer (C74) supplies it and
    silently omitted otherwise, so this helper never has to learn
    how the build stamps its own commit (a C78 concern).

    The single key file ``"structure"`` is byte-compared against
    the staged copy under the run directory; its source is the
    structure path (the str argument, or a loaded
    StructureControl's ``imago_skl`` path)."""
    source = (structure if isinstance(structure, str)
              else getattr(structure, "imago_skl", "imago.skl"))
    scalars = {name: options[name]
               for name in _KEY_SCALAR_NAMES if name in options}
    return KeyFields(
        scalars=scalars,
        files=[KeyFile(name="structure", source=source)])


def _load_structure(structure):
    """Resolve the ``structure`` argument to an object
    ``compute_signature`` can read.  A string is treated as a
    path to an ``imago.skl`` and loaded into a StructureControl;
    an already-loaded StructureControl (or any duck-typed
    structure carrying ``num_atoms`` / ``atom_element_name`` /
    ``space_group_num``) is returned unchanged so a caller that
    already has it parsed need not re-read the file."""
    if isinstance(structure, str):
        loaded = StructureControl()
        loaded.read_imago_skl(structure)
        return loaded
    return structure


def _slug_from_path(path):
    """Derive a filesystem-safe unit id from a structure path: the
    base filename without its extension, lower-cased, with any run
    of non-slug characters collapsed to a single hyphen.  Used
    only when the caller passes no explicit ``id`` (PSEUDOCODE
    15.6); a non-path structure with no id is an error because no
    stable id can be derived from an in-memory object."""
    if not isinstance(path, str):
        raise KaleidoscopeError(
            "predict_settings needs an explicit id when structure "
            "is not a path (cannot derive a slug from an object)")
    stem = os.path.splitext(os.path.basename(path))[0]
    slug = re.sub(r"[^a-z0-9_-]+", "-", stem.lower()).strip("-")
    if not slug:
        raise KaleidoscopeError(
            f"cannot derive a slug id from structure path {path!r}")
    return slug


# ------------------------------------------------------------------
#  Attaching the prediction without teaching the core about it
# ------------------------------------------------------------------

def attach_prediction_record(flight, record):
    """Stash a PredictionRecord on the flight as a plain dict under
    ``metadata["prediction"]`` (PSEUDOCODE 15.6).  The
    domain-agnostic core only round-trips this opaque table
    (serialize_flight emits it as ``[flight.prediction]``); the
    harvest step reads it back.  ``asdict`` recurses through the
    nested ``feature_vector`` Signature so the whole record is
    plain scalars, arrays, and a sub-table."""
    flight.metadata["prediction"] = asdict(record)


# ------------------------------------------------------------------
#  The builder
# ------------------------------------------------------------------

def predict_settings(structure, options, dataspace, system_type,
                     verify=True, id=None, root=""):
    """Build a predict-then-verify ``Flight`` for one structure
    (DESIGN 6.2.8 / 7.7; PSEUDOCODE 15.6).

    Parameters
    ----------
    structure
        An ``imago.skl`` path, or an already-loaded
        StructureControl.
    options
        The fixed (non-swept) makeinput options, held constant
        across every grid point.  MUST carry ``basis``,
        ``functional``, and ``kpoint_integration`` -- together
        they select the predictor sub-model (DESIGN 7.6 step 2),
        so a prediction never mixes incompatible settings.
    dataspace
        The historical-guidance ``Dataspace`` loaded by
        ``guidance_db.load`` (DESIGN 7.4).
    system_type
        One of the four valid system types (DESIGN 7.2),
        declared by the caller.
    verify
        When False, *trust mode*: a length-1 grid at the
        predicted value, no widening (DESIGN 6.2.1).
    id
        The flight-level unit id (a slug).  Derived from the
        structure path when omitted.
    root
        The workspace root for the returned flight.  PSEUDOCODE
        15.6 leaves this to the caller; it is exposed here as an
        optional argument so the producer can set it directly
        rather than mutating the returned flight.

    Returns
    -------
    (Flight, PredictionRecord)
        The flight is ready to dispatch (once ``root`` is set);
        the record is also stashed in ``flight.metadata`` so the
        harvest recovers it.
    """
    # The three sub-model-selecting options are mandatory.
    for required in ("basis", "functional", "kpoint_integration"):
        if required not in options:
            raise KaleidoscopeError(
                f"predict_settings options must carry "
                f"{required!r} (it selects the predictor "
                f"sub-model)")

    # 1. Feature signature for the query structure (DESIGN 7.4).
    resolved = _load_structure(structure)
    query_sig = compute_signature(
        resolved, system_type, dataspace.group_table)

    # 2. Ask the predictor.  It always returns a PredictionResult
    #    (never None), so no None-guards are needed below.
    result = predict(
        dataspace, query_sig, options["basis"],
        options["functional"], options["kpoint_integration"])

    # 3. Choose the grid of k-densities and record which path we
    #    took (DESIGN 7.7 step 3).
    if not verify:
        # Trust mode: one point at the predicted value.
        grid_values = [result.predicted_kpoint_density]
        policy = "trust_no_verify"
    elif result.is_under_trained:
        # No usable prior: fall back to the wide-grid default.
        grid_values = default_wide_kpoint_density_grid()
        policy = "wide_grid_no_prior"
    else:
        # Predict-then-verify with confidence-aware widening.
        grid_values = build_verification_grid(
            result.predicted_kpoint_density, result.confidence)
        policy = "verify_around_prediction"

    # 4. Round + dedupe to integer k-densities so the on-disk tag
    #    parses back to exactly the swept value (DESIGN 6.2.4); a
    #    grid where rounding merged two close points collapses.
    kpd_grid = sorted(set(round(value) for value in grid_values))

    unit_id = id if id is not None else _slug_from_path(structure)
    units = []
    for kpd_int in kpd_grid:
        unit_options = dict(options)
        unit_options["kpd"] = kpd_int        # makeinput options key
        # The display name in the calc-tag tree (DESIGN 6.2.4);
        #   the helper is the translation point between the
        #   makeinput "kpd" key and the "kpt-density" axis name.
        calc_axes = {"kpt-density": kpd_int}
        units.append(CalcUnit(
            id=unit_id,
            calc=build_calc_tag(calc_axes),
            structure=structure,
            options=unit_options,
            wingbeat="imago",
            key_fields=standard_key_fields(structure, options)))

    # 5. Record the sweep shape so serialize_flight emits
    #    [flight.sweep] and the harvest recovers the varied axis
    #    without parsing run-directory paths (DESIGN 6.2.8 step 6).
    sweep = SweepRecord(
        varied_axes=("kpt-density",),
        fixed_axes={
            "basis": options["basis"],
            "functional": options["functional"],
            "kpoint_integration": options["kpoint_integration"]})

    # 6. Assemble the full prediction record (DESIGN 7.7 step 5).
    record = PredictionRecord(
        policy=policy,
        predicted_value=result.predicted_kpoint_density,
        confidence=result.confidence,
        is_under_trained=result.is_under_trained,
        neighbor_entry_ids=result.neighbor_entry_ids,
        predicted_gap=result.predicted_gap,
        predicted_spin_pol=result.predicted_spin_pol,
        system_type=system_type,
        feature_vector=query_sig)

    flight = Flight(root=root, units=units, sweep=sweep)
    attach_prediction_record(flight, record)
    return flight, record
