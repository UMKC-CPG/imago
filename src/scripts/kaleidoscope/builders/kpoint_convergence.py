"""kaleidoscope.builders.kpoint_convergence -- the k-point
convergence builder the adaptive mesh climb uses (DESIGN 7.7 /
3.12; PSEUDOCODE 4e.7).

This is an option-axis builder in the flight-builder split (DESIGN
6.2; see ``kaleidoscope.builders``).  It provides the two pieces
the historical-guidance producer needs to converge a solid's
k-point sampling: ``predict_kpoint_density`` looks up a converged
k-point density for "a structure + a loaded historical-guidance
``Dataspace`` + a sub-model" (predicting from chemically-similar
systems, laying no grid), and ``build_mesh_unit`` builds one
explicit-mesh convergence ``CalcUnit``.  It also defines the
``PredictionRecord`` the later harvest recovers and the calc-tag
encoding for the density and mesh axes.

The producer's climb (``mesh_climb``; DESIGN 3.12) drives these:
it seeds from the predicted density, then walks a sequence of
symmetry-compatible meshes -- one ``build_mesh_unit`` per rung --
until the total energy flattens.  This module supplies the
prediction and the per-mesh unit; the search itself lives in the
climb.

**Why this is a builder, not part of the dispatch core.**  The
core (``model`` / ``workspace`` / ``cache`` / ``dispatch``) is
deliberately domain-agnostic: it runs and tracks calculations
without knowing what they compute (VISION Principle 9/12).  This
module is the opposite -- it consults the guidance database, knows
the swept knob is a k-point mesh, and knows how to spell the
per-calc directory tag.  To keep the dumb core's import graph free
of the physics layer, a client imports it explicitly::

    from kaleidoscope.builders.kpoint_convergence \
        import predict_kpoint_density, build_mesh_unit

Importing it pulls in ``guidance_db`` and ``structure_control``;
the core never does.
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

    - ``policy``             : which prediction / grid path was
                               taken (``trust_no_verify`` /
                               ``wide_grid_no_prior`` /
                               ``verify_around_prediction`` /
                               ``curator_override`` for the density
                               grid; ``predict_then_climb`` for the
                               adaptive mesh climb, 7.7 / 3.12).
    - ``predicted_kpoint_density`` : the predicted converged
                               k-density (or the curator-pinned
                               value in override mode).
    - ``confidence``         : the predictor's combined confidence
                               in [0, 1] that drove the grid width.
    - ``is_under_trained``   : True when the dataspace was too thin
                               for a trustworthy prediction.
    - ``neighbor_entry_ids`` : the guidance entries the prediction
                               drew on (provenance).
    - ``predicted_gap`` /
      ``predicted_magnetization`` : the intermediate stage-1
                               quantities -- the predicted gap (eV)
                               and the predicted intensive moment
                               (Bohr magnetons per atom); None for
                               non-crystalline systems.
    - ``system_type``        : the declared system type.
    - ``feature_vector``     : the query ``Signature`` the
                               prediction was keyed on (recorded
                               for forensics; the harvest recomputes
                               its own from the .skl rather than
                               reading this back).  A guidance_db
                               ``Signature`` bundles four feature
                               fields: ``system_type`` (the
                               four-way label), the 13-long
                               ``composition_vector`` (atom-fraction
                               weight per element group, summing to
                               1), ``lattice_family`` (the crystal
                               family string -- empty when non-
                               crystalline), and its 6-long
                               ``lattice_onehot`` (all zeros when
                               non-crystalline).  Typed ``object``
                               (not ``Signature``) only to spare
                               this builder importing that type;
                               ``asdict`` serializes it to a
                               sub-table either way.
    - ``basis`` / ``functional`` / ``kpoint_integration`` : the
                               sub-model the run uses.  Carried on
                               the per-structure record (and ONLY
                               here -- never duplicated into
                               ``sweep.fixed_axes``) so a combined
                               multi-structure flight whose
                               structures differ in sub-model is
                               still harvestable: each structure's
                               harvest reads its own sub-model back
                               from its own record (DESIGN 6.2.9 /
                               7.8 step 3f).
    """
    policy: str
    predicted_kpoint_density: float | None
    confidence: float
    is_under_trained: bool
    neighbor_entry_ids: tuple = ()
    predicted_gap: float | None = None
    predicted_magnetization: float | None = None
    system_type: str = ""
    feature_vector: object = None
    basis: str = ""
    functional: str = ""
    kpoint_integration: str = ""


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


def encode_mesh_value(mesh):
    """Encode an axial-count mesh as a slug-safe calc-tag value
    (DESIGN 6.2.4 / 7.7): the three counts joined by hyphens, e.g.
    ``[4, 4, 4] -> "4-4-4"``.  Axial counts are always positive
    integers, so the token needs none of ``encode_axis_value``'s
    sign / decimal escaping and stays ``[a-z0-9-]``.  Paired with
    ``decode_mesh_value`` so the harvest reads the mesh back off the
    calc tag (DESIGN 7.7)."""
    return "-".join(str(int(count)) for count in mesh)


def decode_mesh_value(token):
    """Invert :func:`encode_mesh_value`: ``"4-4-4" -> [4, 4, 4]``.
    The hyphen-separated counts parse straight back to a list of
    three integers (DESIGN 6.2.4 / 7.7)."""
    return [int(part) for part in token.split("-")]


# ------------------------------------------------------------------
#  Cache identity and structure handling
# ------------------------------------------------------------------

# The producer's cache identity for a converged-potential run
#   (DESIGN 6.2.1): the scalar settings that define "the same
#   calculation" plus the structure file, byte-compared.
_KEY_SCALAR_NAMES = ("converg", "imago_commit")


def standard_key_fields(structure, options):
    """Build the ``KeyFields`` cache identity for a unit
    (DESIGN 6.2.1/6.2.5): the structure file byte-compared, plus
    the scalar settings that define run identity.

    The scalar names are taken from ``options`` when present
    (DESIGN 6.2.1/6.2.10 list ``converg`` and ``imago_commit``).
    ``converg`` (the SCF convergence limit) is naturally in the
    makeinput options; the build-identity ``imago_commit`` is
    producer-injected -- it is
    carried in ``options`` when the producer (C74) supplies it and
    silently omitted otherwise, so this helper never has to learn
    how the build stamps its own commit (a C78 concern).

    The single key file ``"structure.dat"`` is byte-compared
    against the staged copy under the run directory.  It is
    makeinput's RESOLVED output, not the raw skeleton, so it
    bakes in every input that changes the result (the type/species
    assignment, basis, functional, potential); any of those
    changing misses the cache on its own.  The ``source`` set here
    is provisional (the skeleton path); the producer's prepare
    step re-points it at the built ``structure.dat`` before the
    hit-test (DESIGN 6.2.5, Model A)."""
    source = (structure if isinstance(structure, str)
              else getattr(structure, "imago_skl", "imago.skl"))
    scalars = {name: options[name]
               for name in _KEY_SCALAR_NAMES if name in options}
    return KeyFields(
        scalars=scalars,
        files=[KeyFile(name="structure.dat", source=source)])


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


# ------------------------------------------------------------------
#  Predict-then-climb: the two halves the producer uses
#  (DESIGN 7.7; PSEUDOCODE 4e.7)
#
#  The adaptive climb (DESIGN 3.12) seeds from a prediction but lays
#  its OWN rungs -- explicit meshes, one per calc -- so the producer
#  reaches for two independent pieces here: predict_kpoint_density
#  makes the prediction (no grid), and build_mesh_unit builds one
#  explicit-mesh convergence unit that the producer's round adapter
#  dispatches.  The climb (mesh_climb) chooses which meshes to build.
# ------------------------------------------------------------------

def build_mesh_unit(structure, options, mesh, id):
    """Build one explicit-mesh convergence ``CalcUnit`` for the
    adaptive climb (DESIGN 7.7; PSEUDOCODE 4e.7).

    The climb searches in mesh space, so a rung is dispatched as an
    explicit mesh rather than a density.  ``scfkp`` is makeinput's
    key for an explicit axial-count mesh (a style-code-1 k-point
    file: axial counts, shift, and point operations); imago resolves
    its own symmetry shift and irreducible-wedge reduction from those
    counts, so the requested mesh is honoured exactly (DESIGN 7.7).
    ``kpt-mesh`` is the calc-tag axis the count triple renders under
    (``kpt-mesh-<a>-<b>-<c>``), mirroring the density path's
    ``kpt-density``.

    The cache identity is the SAME one the density units used
    (:func:`standard_key_fields`; DESIGN 6.2.1), so a mesh re-run in
    a later climb round is a cache hit and costs nothing.

    Parameters
    ----------
    structure
        An ``imago.skl`` path or a loaded StructureControl, copied
        verbatim onto the unit as the density builder does.
    options
        The fixed run settings in each tool's coded vocabulary; a
        copy gains the ``scfkp`` mesh so the caller's dict is left
        untouched.
    mesh
        The axial counts ``[a, b, c]`` to dispatch.
    id
        The unit id.  The producer passes the material's reference
        id, which the round adapter reads back to route the energy.
    """
    unit_options = dict(options)
    unit_options["scfkp"] = list(mesh)
    # Assembled like build_calc_tag but with the mesh encoder: that
    #   helper's encode_axis_value handles a scalar, not a count
    #   triple.  "kpt-mesh" is a static slug, safe as a directory
    #   level (DESIGN 6.2.4).
    calc = ("kpt-mesh-" + encode_mesh_value(mesh),)
    return CalcUnit(
        id=id,
        calc=calc,
        structure=structure,
        options=unit_options,
        wingbeat="imago",
        key_fields=standard_key_fields(structure, options))


def predict_kpoint_density(structure, dataspace, system_type,
                           submodel, center=None):
    """Predict the converged k-point density for one structure,
    laying no grid (DESIGN 7.7; PSEUDOCODE 4e.7).

    This is the prediction half of the producer's convergence
    search: the query signature, the predictor call, and the
    per-structure ``PredictionRecord``.  The adaptive climb (DESIGN
    3.12) seeds from the returned density and picks its dispatch
    mode and persistence from the confidence, then lays its own
    rungs (:func:`build_mesh_unit` per rung), so no grid is built
    here -- only the prediction.

    A curator-pinned ``center`` (the 5.7 ``kpoint_spec`` density
    override) BYPASSES the predictor: the density is the pinned
    value at full confidence, and the record documents the override.

    Returns
    -------
    (density, confidence, is_under_trained, record)
        ``density`` seeds the climb (DESIGN 3.12.4); ``confidence``
        and ``is_under_trained`` drive the confidence policy (3.12.6 /
        7.9); ``record`` is the ``PredictionRecord`` the harvest
        recovers (7.8).
    """
    for required in ("basis", "functional", "kpoint_integration"):
        if required not in submodel:
            raise KaleidoscopeError(
                f"predict_kpoint_density submodel must carry "
                f"{required!r} (it selects the predictor sub-model)")

    resolved = _load_structure(structure)
    query_sig = compute_signature(
        resolved, system_type, dataspace.group_table)

    if center is not None:
        # Curator override (5.7 kpoint_spec.density): the predictor
        #   is never consulted; the density is the pinned value.
        record = PredictionRecord(
            policy="curator_override",
            predicted_kpoint_density=float(center),
            confidence=1.0,
            is_under_trained=False,
            neighbor_entry_ids=(),
            predicted_gap=None,
            predicted_magnetization=None,
            system_type=system_type,
            feature_vector=query_sig,
            basis=submodel["basis"],
            functional=submodel["functional"],
            kpoint_integration=submodel["kpoint_integration"])
        return float(center), 1.0, False, record

    result = predict(
        dataspace, query_sig, submodel["basis"],
        submodel["functional"], submodel["kpoint_integration"])
    record = PredictionRecord(
        policy="predict_then_climb",
        predicted_kpoint_density=result.predicted_kpoint_density,
        confidence=result.confidence,
        is_under_trained=result.is_under_trained,
        neighbor_entry_ids=result.neighbor_entry_ids,
        predicted_gap=result.predicted_gap,
        predicted_magnetization=result.predicted_magnetization,
        system_type=system_type,
        feature_vector=query_sig,
        basis=submodel["basis"],
        functional=submodel["functional"],
        kpoint_integration=submodel["kpoint_integration"])
    return (result.predicted_kpoint_density, result.confidence,
            result.is_under_trained, record)
