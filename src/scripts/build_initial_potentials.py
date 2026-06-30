#!/usr/bin/env python3
"""build_initial_potentials.py -- the augmented potential
database *producer* (DESIGN 5.7; PSEUDOCODE 11.4).

Role in the pipeline
--------------------
This script is the *producer* half of the library / producer /
consumer split documented in DESIGN 5.4.  It builds the augmented
per-element potential database files
(``share/atomicPDB/<elem>/s_gaussian_pot.toml``) that the consumer
(``makeinput.py``, C47) later reads.  Its inputs are a curated set
of reference solids; for each it refreshes the atomSCF-derived
``"isolated"`` baseline and -- once the SCF machinery lands -- runs
(or reuses a cached) Imago SCF, harvests the converged potential at
named atom sites, and writes the results via
``initial_potential_db.save()``.

The curation manifest -- what it is and why
-------------------------------------------
The pipeline's primary input is one human-readable file: the
**curation manifest** (``share/atomicBDB/manifest.toml``).  It
declares *which solids, which atom sites, under which labels, with
which SCF settings*.  Its three jobs (DESIGN 5.7):

1. **Declare the curation set.**  The manifest *is* the curation
   strategy made explicit.  Adding a reference solid means adding a
   manifest entry; reviewing the curation set means reading the
   manifest -- nothing is hidden inside the script.
2. **Tell the pipeline what to harvest.**  For each reference
   solid, which atom sites' converged potentials enter the
   database, and under what labels.
3. **Record the SCF settings used.**  k-points, convergence
   threshold, etc., copied into the provenance fields of DESIGN 5.2
   so every database entry carries the conditions of its reference
   run.

VISION Principle 5 ("the database must be regeneratable from the
curated set, not a hand-edited artifact") rules out hardcoding the
curation set inside this script and rules out folder-of-files
conventions that lose metadata.  The manifest is the smallest piece
of structured data that closes the gap: every curation choice
captured in one version-controlled file alongside the structure
files it points at, so regeneration becomes a deterministic
function of (manifest, structure files, Imago build).

Build analogy
-------------
The manifest is the **build configuration** for the database.  The
structure files are the **source**.  The Imago build is the
**toolchain**.  ``build_initial_potentials.py`` is the **build
script**.  The augmented database is the **compiled output**.  Same
role ``pyproject.toml`` plays for a Python package, or a Makefile
plays for a binary.

Reproducibility (layered)
-------------------------
- **Emitter determinism (bit-level, strict).**  The TOML emitter
  (DESIGN 5.5, in ``initial_potential_db``) writes byte-identical
  bytes for a fixed in-memory ``ElementDatabase``.
- **Pipeline numerical output (precision-level, loose).**  Given
  the same manifest, ``pot1``/``coeff1`` files, and Imago build, the
  numerical outputs agree at the precision the SCF / fit chain can
  reach -- bit-identity is *not* promised (accumulation order,
  threading, library versions can shift the last bits).
- **Provenance metadata (free).**  Timestamps and commit SHAs
  refresh every run and carry no reproducibility guarantee.

Implementation status
----------------------
This file is being built incrementally (C48):

* **C48.1 (landed) -- the manifest reader.**
  :func:`load_manifest_v2` and the manifest dataclasses.  Enforces
  validation rules 1-8 of DESIGN 5.7; rule 9 (method must be a
  registered matcher) is gated on an optional ``known_methods``
  argument and skipped when it is ``None``, exactly as
  ``initial_potential_db.load`` handles its own rule 9 -- the
  matcher registry lives in ``makeinput.py`` (C54) and does not
  exist yet.
* **C48.2 (landed) -- the ``"isolated"`` baseline refresh.**
  :func:`refresh_isolated_entries` (step 1 of the pipeline) loads
  or creates each element's database and rebuilds its
  ``"isolated"`` entry from the current ``pot1``/``coeff1`` files
  via :func:`build_isolated_entry`; :func:`save_databases` writes
  them back.  This alone forms a working producer->consumer loop
  with no SCF: every element's file gets its rule-6 baseline,
  default-tagged when the manifest curates nothing else for it.
* **C48.3 -- SCF orchestration, the content-keyed cache, COD fetch,
  the run log, and the CLI flags** (next; needs a live Imago
  toolchain to exercise).
* **C60 -- fingerprint harvest** (Phase 2): compute and attach the
  ``[[reference_solid.entry.fingerprint]]`` records.  The Python-side
  (reduce) harvest is implemented (the light half): it computes the
  shell code in process from the run's expanded structure and stores
  an element-only ``shell_code`` (DESIGN 5.2 / 5.7).  The Fortran-side
  (bispectrum) harvest -- which needs a dispatched ``-loen`` unit --
  is deferred to C55/C58, and a Fortran-side declaration is refused
  rather than silently dropped.
"""

import argparse
import os
import re
import subprocess
import sys
import tomllib
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any

import initial_potential_db as ipdb
import guidance_db
import guidance_harvest
# The curation-manifest schema (dataclasses + readers) lives in its
#   own neutral library so the producer and the authoring tool
#   (expand_manifest) share one definition; the producer imports the
#   names it runs a manifest with (DESIGN 5.7).
from curation_manifest import (
    ReferenceEntry, ReferenceSolid, CurationManifest,
    load_manifest_v2, load_structure_sources)
from kaleidoscope import CalcUnit, Flight, SweepRecord, dispatch
from kaleidoscope.builders.kpoint_convergence import (
    build_kpoint_convergence, standard_key_fields)
from kaleidoscope.cluster_config import (
    resolve_dispatch, write_resolved_dispatch)
from kaleidoscope.workspace import read_status, toml_line
# The Phase-2 matcher registry (ARCHITECTURE 8.9) lives in the neutral
#   matchers module; the fingerprint harvest dispatches reduce
#   fingerprints through it.  StructureControl reads the run's expanded
#   structure for those shells.
from matchers import MATCHERS
from structure_control import StructureControl
# makegroups owns the two loen-side helpers the producer shares: the
#   sub_spec -> -loeninput value mapping (so loen units and the
#   makegroups grouping flow emit identical LOEN blocks) and the
#   <edge>_loen<basis>.plot descriptor finder (validated live, DESIGN
#   5.10.3).  Keeping them in one place stops the producer and
#   makegroups from drifting on the loen seam.
import makegroups


# ============================================================
#  Isolated-baseline refresh (DESIGN 5.7 step 2; PSEUDOCODE
#  11.4 build_isolated_entry / read_pot1 / read_coeff1)
# ============================================================
#
# Step 1 of the pipeline always rebuilds every element's
# "isolated" entry from the *current* legacy pot1/coeff1 files,
# so any change in atomSCF output propagates into the augmented
# database on the next run.  The isolated entry is also the
# rule-6 baseline every per-element file must carry, and -- for
# elements the manifest does not curate -- it is the file's
# sole entry and therefore its default-tagged one (rule 7).
#
# The legacy pot1/coeff1 text format is the same fixed-line
# layout the C47 consumer materializes; these parsers are the
# inverse of makeinput._write_legacy_pot_files_from_entry.


@dataclass
class _PotFileData:
    """The scalar fields parsed out of a legacy ``pot`` file.

    The per-term coefficients and alphas live in the companion
    ``coeff`` file (see :func:`_parse_coeff_file`); this holds
    only the per-entry scalars the ``pot`` file carries.
    """

    nuclear_z: float
    nuclear_alpha: float
    covalent_radius: float
    num_gaussians: int
    alpha_min: float
    alpha_max: float


def _parse_pot_file(path: str) -> _PotFileData:
    """Parse a legacy ``pot`` file's fixed eight-line layout.

    The layout (atomSCF output, and what the C47 consumer
    regenerates) is positional::

        0  NUCLEAR_CHARGE__ALPHA
        1  <Z> <nuclear_alpha>
        2  COVALENT_RADIUS
        3  <covalent_radius>
        4  NUM_ALPHAS
        5  <num_gaussians>
        6  ALPHAS
        7  <alpha_min> <alpha_max>

    ``Z`` is written by atomSCF as a float (``79.000000``) and is
    kept as the real ``nuclear_z`` the schema (DESIGN 5.2) uses --
    nominally an integer, but Imago consumes Z as a real number.
    The four tag lines are checked so a malformed file fails with
    a clear message rather than a silent misparse.
    """

    with open(path) as handle:
        lines = [ln.rstrip("\n") for ln in handle
                 if ln.strip() != ""]

    if len(lines) < 8:
        raise ValueError(
            f"{path}: malformed pot file (expected 8 content "
            f"lines, found {len(lines)})")
    for index, tag in ((0, "NUCLEAR_CHARGE__ALPHA"),
                       (2, "COVALENT_RADIUS"),
                       (4, "NUM_ALPHAS"), (6, "ALPHAS")):
        if lines[index].strip() != tag:
            raise ValueError(
                f"{path}: malformed pot file (line {index} "
                f"expected tag {tag!r}, found "
                f"{lines[index].strip()!r})")

    charge_alpha = lines[1].split()
    alpha_range = lines[7].split()
    return _PotFileData(
        nuclear_z=float(charge_alpha[0]),
        nuclear_alpha=float(charge_alpha[1]),
        covalent_radius=float(lines[3].split()[0]),
        num_gaussians=int(lines[5].split()[0]),
        alpha_min=float(alpha_range[0]),
        alpha_max=float(alpha_range[1]))


def _parse_coeff_file(path: str) -> tuple[list[float],
                                          list[float]]:
    """Parse a HEADERLESS ``coeff1`` file into (coefficients,
    alphas) -- one element's potential, the atomicPDB layout.

    The file is a bare count line followed by one line per Gaussian
    term; the term lines carry five whitespace-separated columns
    of which only column 1 (the coefficient) and column 2 (the
    alpha) are meaningful here -- columns 3-5 are the placeholder
    fields Imago ignores (see C47).  The count line is
    cross-checked against the number of term lines so a truncated
    or padded file is caught at parse time.

    This is NOT the converged ``scfV`` output layout: that file
    carries every potential type behind a ``NUM_TYPES`` header and
    per-channel tags (``TOTAL__OR__SPIN_UP`` / ``SPIN_DN``), parsed
    by ``_parse_scfv_type_block`` instead (DESIGN 5.7).
    """

    with open(path) as handle:
        lines = [ln for ln in handle if ln.strip() != ""]

    declared_count = int(lines[0].split()[0])
    coefficients: list[float] = []
    alphas: list[float] = []
    for term_line in lines[1:]:
        tokens = term_line.split()
        coefficients.append(float(tokens[0]))
        alphas.append(float(tokens[1]))

    if len(coefficients) != declared_count:
        raise ValueError(
            f"{path}: coeff file count line says {declared_count} "
            f"terms but {len(coefficients)} term lines follow")
    return coefficients, alphas


def element_path(pdb_root: str, elem: str) -> str:
    """Return the per-element database path under ``pdb_root``.

    Matches both the consumer's lookup (C47) and PSEUDOCODE
    11.4 ``element_path``: the element directory name is
    lower-cased, so ``Au`` and ``au`` resolve to the same
    ``share/atomicPDB/au/s_gaussian_pot.toml``.
    """

    return os.path.join(pdb_root, elem.lower(),
                        "s_gaussian_pot.toml")


def is_isolated_default_for(elem: str,
                            manifest: CurationManifest) -> bool:
    """True iff the manifest curates no default entry for ``elem``.

    The isolated baseline carries the per-element file's
    ``default = true`` tag exactly when the manifest declares no
    other default-tagged customization for the element (PSEUDOCODE
    11.4 ``is_isolated_default_for``).  A ``default = true``
    customization that names this element wins over the baseline; an
    element no such customization names falls back to its isolated
    entry as the default.  A customization may omit its ``element``
    (DESIGN 5.7 rule 3) -- the harvest attributes it once the run
    reveals the site's element -- so a default customization without
    an element cannot be credited to any element here and is skipped;
    the interim build honours only element-named default
    customizations.  Element symbols are compared case-insensitively
    because the manifest uses proper case (``Au``) while directory
    names are lower case (``au``).
    """

    for solid in manifest.reference_solids:
        for entry in solid.entries:
            if (entry.element is not None
                    and entry.element.lower() == elem.lower()
                    and entry.default):
                return False
    return True


def build_isolated_entry(pdb_root: str, elem: str, commit: str,
                         timestamp: str,
                         manifest: CurationManifest
                         ) -> ipdb.PotentialEntry:
    """Build the ``"isolated"`` entry from current pot1/coeff1.

    Reads ``<pdb_root>/<elem>/pot1`` and ``coeff1``, cross-checks
    that the term counts agree, and returns a fresh
    :class:`initial_potential_db.PotentialEntry` tagged
    ``"isolated"`` with atomSCF-source provenance.  The
    ``default`` flag is computed from the manifest via
    :func:`is_isolated_default_for`.  No fingerprints are
    attached -- the baseline never participates in environment-
    scheme matching.
    """

    elem_dir = os.path.join(pdb_root, elem.lower())
    pot = _parse_pot_file(os.path.join(elem_dir, "pot1"))
    coefficients, alphas = _parse_coeff_file(
        os.path.join(elem_dir, "coeff1"))

    if not (len(coefficients) == len(alphas)
            == pot.num_gaussians):
        raise ValueError(
            f"{elem}: pot1/coeff1 disagree on term count "
            f"(pot num_gaussians={pot.num_gaussians}, "
            f"coeff coefficients={len(coefficients)}, "
            f"alphas={len(alphas)})")

    symbol = elem.capitalize()
    return ipdb.PotentialEntry(
        label="isolated",
        default=is_isolated_default_for(elem, manifest),
        description=(f"Single isolated {symbol} atom "
                     f"(from atomSCF)."),
        num_gaussians=pot.num_gaussians,
        alpha_min=pot.alpha_min,
        alpha_max=pot.alpha_max,
        coefficients=coefficients,
        # A single observed potential: the running mean is just
        # this potential, its spread is zero, and the dedup
        # counts start at 1 (DESIGN 5.2.3).
        coefficient_std=[0.0] * pot.num_gaussians,
        alphas=alphas,
        provenance={
            "source": "atomSCF",
            "commit": commit,
            "generated_at": timestamp},
        fingerprints=[])


def list_element_dirs(pdb_root: str) -> list[str]:
    """Return the element directory names under ``pdb_root``.

    Only directories that actually carry a ``pot1`` file are
    returned, so non-element siblings (a ``cache`` directory, a
    stray ``manifest.toml``) are skipped.  The list is sorted for
    deterministic processing order.
    """

    names = []
    for name in os.listdir(pdb_root):
        if os.path.isfile(os.path.join(pdb_root, name, "pot1")):
            names.append(name)
    return sorted(names)


def refresh_isolated_entries(pdb_root: str,
                             manifest: CurationManifest,
                             commit: str, timestamp: str,
                             elements: list[str] | None = None
                             ) -> dict[str, ipdb.ElementDatabase]:
    """Step 1 of the pipeline: reset every element database to a
    fresh isolated baseline (the REGENERATE model, DESIGN 5.7).

    For each element (all element directories under ``pdb_root``,
    or just ``elements`` when given), load the existing
    ``s_gaussian_pot.toml`` if present or create a fresh
    :class:`initial_potential_db.ElementDatabase` from the
    element's ``pot1`` scalars, then **discard all of its prior
    entries** and seed it with a single ``"isolated"`` entry rebuilt
    from the current pot1/coeff1.  The harvest phase appends this
    run's solid entries afterward, so the produced file is a pure
    function of the current pot1/coeff1 and the manifest -- no entry
    from a removed solid lingers, and re-running cannot double-count
    a dedup entry's multiplicity.  Returns the in-memory databases
    keyed by element directory name; the caller saves them (see
    :func:`save_databases`).

    Implements PSEUDOCODE 11.4 step 1.  The existing-file load
    passes ``known_methods=None`` (rule 9 is skipped) because the
    matcher registry is not threaded in here; the loaded entries are
    dropped regardless, so the file is rebuilt rather than edited.
    """

    if elements is None:
        elements = list_element_dirs(pdb_root)

    databases: dict[str, ipdb.ElementDatabase] = {}
    for elem in elements:
        path = element_path(pdb_root, elem)
        if os.path.isfile(path):
            database = ipdb.load(path, known_methods=None)
        else:
            pot = _parse_pot_file(
                os.path.join(pdb_root, elem.lower(), "pot1"))
            database = ipdb.ElementDatabase(
                schema_version=2,
                element_symbol=elem.capitalize(),
                nuclear_z=pot.nuclear_z,
                nuclear_alpha=pot.nuclear_alpha,
                covalent_radius=pot.covalent_radius,
                potentials=[])

        # REGENERATE (DESIGN 5.7): drop *every* prior entry -- the
        # isolated baseline and all previously harvested solid
        # entries alike -- and start the file from a single fresh
        # isolated baseline.  Rebuilding from scratch each run keeps
        # the database a pure function of the current pot1/coeff1 and
        # the manifest: a solid dropped from the manifest leaves no
        # orphan, and re-running never inflates the dedup
        # multiplicity/model_count of an entry (DESIGN 5.2.3).  The
        # harvest phase then appends this run's solid entries.
        database.potentials = [build_isolated_entry(
            pdb_root, elem, commit, timestamp, manifest)]
        databases[elem] = database

    return databases


def save_databases(databases: dict[str, ipdb.ElementDatabase],
                   pdb_root: str) -> None:
    """Write each element database to its on-disk path.

    Implements PSEUDOCODE 11.4 step 3: every affected
    :class:`initial_potential_db.ElementDatabase` is written via
    the deterministic emitter (DESIGN 5.5) to its
    ``element_path``.  The element directory already exists (it
    held the ``pot1`` we read), so no directory creation is
    needed.
    """

    for elem, database in databases.items():
        ipdb.save(database, element_path(pdb_root, elem))


# ============================================================
#  Build identity, workspace, and structure materialization
#  (DESIGN 5.7; PSEUDOCODE 11.4)
# ============================================================

def _git_sha() -> str:
    """The current HEAD commit, or ``"unknown"`` when git is
    unavailable.  Injected into the run options as ``imago_commit``
    so kaleidoscope's run-reuse cache key (DESIGN 6.2.5) and every
    produced entry's provenance both record which build the
    *producer* believed it ran.  This records the producer's belief,
    which can drift from the binary actually executed; C84 hardens
    it by having Imago stamp its own build commit."""

    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True)
        return completed.stdout.strip() or "unknown"
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def _now_iso8601_utc() -> str:
    """The build timestamp in the schema's ISO-8601 UTC form."""

    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def curation_workspace_root(pdb_root: str) -> str:
    """The kaleidoscope workspace the producer dispatches its
    combined flight into (PSEUDOCODE 11.4).  It sits beside the
    databases under the shared data root so its run-reuse cache
    (DESIGN 6.2.5) persists across producer runs and reference
    solids dedupe against earlier builds."""

    data_root = os.path.dirname(pdb_root.rstrip("/"))
    return os.path.join(data_root, "curation", "workspace")


def _cod_extension(ref: ReferenceSolid) -> str:
    """The on-disk extension for a fetched COD structure.  COD
    serves CIF, so v1 always writes ``.cif``; kept as a hook so a
    future format negotiation has a single place to change."""

    return ".cif"


def _fetch_cod_structure(cod_id: int, cod_revision: str,
                         dest: str) -> None:
    """Fetch one pinned COD revision to ``dest`` (DESIGN 5.7,
    Option A).  Delegates to ``cod_fish.fetch_cif``, the canonical
    strict COD fetch (ARCHITECTURE 9.5), so the producer and the
    standalone ``cod_fish`` tool share one implementation.  The
    pinned revision is verified against the served CIF: a mismatch,
    a network outage, or a COD outage raises -- the fetch NEVER
    falls back to a different revision, because a silent fallback
    would desync the reproducible build from the pinned manifest.

    NOTE (C74 end-to-end): the live COD fetch needs network access
    and the COD per-revision API; it is exercised only on the
    cluster.  A ``structure_path`` manifest avoids it entirely for
    offline / unit-test runs."""

    import cod_fish

    try:
        payload = cod_fish.fetch_cif(cod_id, revision=cod_revision)
    except cod_fish.CodFishError as exc:
        raise RuntimeError(
            f"COD fetch failed for cod_id={cod_id} (pinned "
            f"revision {cod_revision!r}): {exc}.  The build pins "
            f"this revision and refuses to fall back.") from exc
    with open(dest, "wb") as handle:
        handle.write(payload)


def structure_cache_dir(pdb_root: str) -> str:
    """The directory the producer caches fetched and converted
    reference structures in -- ``<data_root>/atomicBDB/cache/
    structures``, beside the databases -- keyed by ``reference_id``.
    A full producer run and a ``--materialize-only`` pre-flight share
    it by default, so a pre-flight's converted skeletons are reused by
    the real run with no second fetch.  A pre-flight may redirect to
    its own location (a cache mirror) with ``--materialize-dir``."""

    data_root = os.path.dirname(pdb_root.rstrip("/"))
    return os.path.join(data_root, "atomicBDB", "cache", "structures")


def materialize_structure(ref: ReferenceSolid, manifest_dir: str,
                          pdb_root: str,
                          cache_dir: str | None = None) -> str:
    """Guarantee the reference solid's structure exists as a local
    file and return its path (Option A; DESIGN 5.7 / PSEUDOCODE
    11.4).  This is the producer's ONLY network access and is
    deliberately decoupled from any run cache: it carries no SCF
    state and makes no hit/miss decision.  Recompute avoidance
    belongs to kaleidoscope's run-reuse cache (DESIGN 6.2.5), which
    keys on this file's bytes.

    A ``structure_path`` ref is a plain disk read, resolved under
    the manifest directory (rule 4 already validated it exists); it
    is expected to be an ``imago.skl``.  A ``cod_id`` ref is fetched
    once as a CIF and converted to an ``imago.skl`` with its space
    group preserved (``cif2skl``), because the run consumes a
    skeleton and a crystal's Brillouin-zone integration samples the
    irreducible wedge using that space group (ARCHITECTURE 9.5).
    Both the fetched CIF and the converted skeleton are cached in
    ``cache_dir`` (the shared ``structure_cache_dir`` beside the
    databases when not given); the skeleton is the returned artifact.
    A CIF whose space group ``cif2skl`` cannot resolve is a hard
    error -- the curator then converts it by hand (with ``cif2skl``'s
    ``--space`` override) and supplies the result as a
    ``structure_path`` instead."""

    if ref.structure_path is not None:
        return os.path.join(manifest_dir, ref.structure_path)

    import cif2skl

    if cache_dir is None:
        cache_dir = structure_cache_dir(pdb_root)
    cif_path = os.path.join(
        cache_dir, ref.reference_id + _cod_extension(ref))
    skl_path = os.path.join(cache_dir, ref.reference_id + ".skl")
    if not os.path.exists(skl_path):
        os.makedirs(cache_dir, exist_ok=True)
        if not os.path.exists(cif_path):
            _fetch_cod_structure(ref.cod_id, ref.cod_revision, cif_path)
        try:
            cif2skl.convert(cif_path, skl_path, title=ref.reference_id)
        except cif2skl.CifConversionError as exc:
            raise RuntimeError(
                f"COD structure for {ref.reference_id!r} "
                f"(cod_id={ref.cod_id}) could not be converted to a "
                f"skeleton with its space group preserved: {exc}.  "
                f"Resolve it by hand with cif2skl (its --space "
                f"override) and supply the result as a structure_path "
                f"entry instead.") from exc
    return skl_path


# ============================================================
#  Structure pre-flight (DESIGN 5.7): fetch + convert only
# ============================================================

def materialize_only(manifest_path: str, pdb_root: str,
                     cache_dir: str | None = None
                     ) -> list[dict[str, Any]]:
    """Pre-flight every reference structure in the manifest, then
    STOP -- no SCF dispatch (DESIGN 5.7).  This is Phase 1's
    ``materialize_structure`` loop run on its own, over the relaxed
    ``load_structure_sources`` view, so a curator can validate a
    freshly pinned structure set -- and catch any space group
    ``cif2skl`` cannot resolve -- before filling in the run and
    harvest details.  Because it is the EXACT path a full run uses, a
    clean pre-flight is a guarantee the run will not trip on
    acquisition, and the converted skeletons it caches are reused by
    that run (unless ``cache_dir`` redirects them to a mirror).

    Each solid is attempted independently: a failure is captured and
    the loop continues, so one bad space group never hides the status
    of the rest.  Returns one report row per solid (also printed by
    the CLI), each carrying ``reference_id``, a human-readable
    ``source``, an ``ok`` flag, the ``skl_path`` on success, and the
    error ``message`` on failure."""

    sources = load_structure_sources(manifest_path)
    manifest_dir = os.path.dirname(os.path.abspath(manifest_path))
    if cache_dir is None:
        cache_dir = structure_cache_dir(pdb_root)

    report: list[dict[str, Any]] = []
    for ref in sources:
        if ref.structure_path is not None:
            source = f"structure_path {ref.structure_path}"
        else:
            source = (f"cod_id {ref.cod_id} "
                      f"(revision {ref.cod_revision})")
        row: dict[str, Any] = {
            "reference_id": ref.reference_id, "source": source}
        try:
            skl = materialize_structure(
                ref, manifest_dir, pdb_root, cache_dir=cache_dir)
            row.update(ok=True, skl_path=skl, message="")
        except (RuntimeError, FileNotFoundError, OSError) as exc:
            # A failed solid is reported, not raised: the curator
            #   wants every problem in one pass, not the first only.
            row.update(ok=False, skl_path=None, message=str(exc))
        report.append(row)

    return report


# ============================================================
#  Per-solid options + the (deferred) loen fingerprint units
# ============================================================

# Translation from the manifest's human-readable physics names to
#   the coded option values each tool's argparse dest expects
#   (DESIGN 6.2.10, decision 2).  The producer owns this translation
#   because it is where the physics intent is known; the wingbeat
#   downstream only routes keys by namespace and never value-codes.

# Exchange-correlation functional name -> makeinput ``xccode``
#   integer (the codes in ``share/xc_code.dat``).  v1 covers the
#   non-relativistic, spin-restricted functionals a reference-solid
#   manifest selects; extend this map as that vocabulary grows.
_FUNCTIONAL_TO_XCCODE = {
    "wigner": 100,
    "ceperley-alder": 101,
    "hedin-lundqvist": 102,
    "pbe": 200,
}

# k-point integration method name -> makeinput ``scfkpint`` integer
#   (0 = Gaussian broadening / histogram, 1 = the linear analytic
#   tetrahedron).
_KPOINT_INTEGRATION_TO_SCFKPINT = {
    "gaussian": 0,
    "linear-tetrahedral": 1,
}


def _scfkpint_for(token: str) -> int:
    """Map a k-point integration token to makeinput's ``scfkpint``
    code.  ``linear-tetrahedral`` -> 1; ``gaussian`` -> 0.  A
    smeared Gaussian may name its smearing width in the token
    (e.g. ``gaussian-0.1``); the integration *code* is still 0 --
    the width is a separate concern read by :func:`_thermsmear_for`
    and forwarded as makeinput's thermal smearing sigma, so the two
    travel as distinct options.  An unknown token raises at the
    manifest boundary rather than reaching makeinput as a mystery
    value."""
    if token in _KPOINT_INTEGRATION_TO_SCFKPINT:
        return _KPOINT_INTEGRATION_TO_SCFKPINT[token]
    if token.startswith("gaussian-"):
        return 0
    raise ValueError(
        f"unknown kpoint_integration {token!r}; the producer "
        f"understands {sorted(_KPOINT_INTEGRATION_TO_SCFKPINT)} "
        f"(or 'gaussian-<smearing>')")


def _thermsmear_for(token: str) -> float | None:
    """Extract the thermal (electronic) smearing sigma named in a
    smeared Gaussian integration token, or ``None`` when none is
    named.

    A bare ``gaussian`` (or ``linear-tetrahedral``) names no width,
    so the producer leaves the smearing unset and makeinput keeps
    its rc-sourced ``therm_smear_main`` default.  A
    ``gaussian-<sigma>`` token (e.g. ``gaussian-0.1``) pins that
    sigma for the run; the producer forwards it as makeinput's
    ``thermsmear`` option, which becomes the THERMAL_SMEARING_SIGMA
    field of the SCF input.  The width is an electron-volt
    broadening applied at the Fermi level, matching the engine's
    smearing field.  A malformed width -- a ``gaussian-`` with no
    number, or a non-numeric tail -- raises at the manifest
    boundary rather than reaching makeinput as a mystery value."""
    if not token.startswith("gaussian-"):
        return None
    width_text = token[len("gaussian-"):]
    try:
        return float(width_text)
    except ValueError:
        raise ValueError(
            f"malformed kpoint_integration {token!r}; the smearing "
            f"width after 'gaussian-' must be a number in eV "
            f"(e.g. 'gaussian-0.1')")


def make_producer_options(ref: ReferenceSolid,
                          imago_commit: str) -> dict[str, Any]:
    """The fixed (non-swept) run settings for a reference solid's
    convergence flight, in each tool's own coded vocabulary
    (DESIGN 6.2.10 / PSEUDOCODE 11.4).

    The manifest records physics in human-readable form
    (``functional = "wigner"``, ``kpoint_integration =
    "linear-tetrahedral"``, ``basis = "fb"``); this function
    translates that into the dest-keyed, coded options the tools'
    ``from_options`` APIs require, so the wingbeat (DESIGN 6.2.2 /
    13.2) can route each key purely by namespace without
    value-coding:

      - ``functional``         -> ``xccode``   (wigner = 100, ...)
      - ``kpoint_integration`` -> ``scfkpint`` (LAT = 1, ...)
      - ``basis``              -> ``scf_basis`` (an imago run-time
                                  selection; imago codes ``fb -> 2``)
      - ``scf_threshold``      -> ``converg``  (the makeinput SCF
                                  convergence limit)
      - ``gaussian-<sigma>``   -> ``thermsmear`` (the makeinput
                                  THERMAL_SMEARING_SIGMA, in eV;
                                  emitted only when the integration
                                  token names a smearing width)
      - ``kpoint_spec.shift``  -> ``kpshift``

    Also carried: ``imago_commit``, the build identity that (with
    ``converg``) forms kaleidoscope's run-reuse cache key (DESIGN
    6.2.5); the wingbeat drops it before forwarding.

    This dict carries NO physics-name keys.  The human sub-model
    (basis / functional / kpoint_integration) travels to the
    builder in its own ``submodel`` dict (DESIGN 6.2.8 / 6.2.10),
    never mixed in here, so makeinput never sees a name it would
    reject.  The swept k-density is also absent: the builder adds
    it per grid point (the ``kpd`` option), so a fixed value here
    would collide with the sweep."""

    try:
        xccode = _FUNCTIONAL_TO_XCCODE[ref.functional]
    except KeyError:
        raise ValueError(
            f"unknown functional {ref.functional!r}; the producer "
            f"understands {sorted(_FUNCTIONAL_TO_XCCODE)}")

    options: dict[str, Any] = {
        "scf_basis": ref.basis,
        "xccode": xccode,
        "scfkpint": _scfkpint_for(ref.kpoint_integration),
        "converg": ref.scf_threshold,
        "imago_commit": imago_commit,
    }
    shift = ref.kpoint_spec.get("shift")
    if shift is not None:
        options["kpshift"] = shift
    # A smeared Gaussian token (``gaussian-<sigma>``) pins the
    #   thermal smearing for the run; a bare token names no width,
    #   so the key is omitted and makeinput keeps its rc default.
    smearing_sigma = _thermsmear_for(ref.kpoint_integration)
    if smearing_sigma is not None:
        options["thermsmear"] = smearing_sigma
    return options


def _slug_safe(text: Any) -> str:
    """Lowercase ``text`` and replace every character a kaleidoscope
    slug forbids (anything outside ``[a-z0-9_-]``) with ``_``.

    Run-directory components must be slugs (DESIGN 6.2.4), but a
    sub_spec value formatted for a tag -- a float like ``9.0`` -> ``9``
    or ``0.85`` -> ``0.85`` -- can carry a dot, so it is sanitized here
    before it ever reaches the directory name."""

    return re.sub(r"[^a-z0-9_-]", "_", str(text).lower())


def _sub_spec_slug(sub_spec: dict[str, Any]) -> str:
    """Deterministic slug for a fingerprint sub_spec (PSEUDOCODE 11.4
    ``sub_spec_slug``).

    Keys are taken in alphabetical order and joined as ``key_value``
    segments, hyphen-separated; floats format as ``%.6g`` (long enough
    to tell apart the parameters humans pick, short enough for a
    directory name).  Both halves are slug-sanitized.  The same
    (method, sub_spec) therefore always produces the same slug, so the
    loen unit a declaration builds (in :func:`build_loen_units`) and the
    descriptor the harvest reads (in :func:`harvest_loen_fingerprint`)
    resolve to one and the same run directory."""

    parts = []
    for key in sorted(sub_spec):
        value = sub_spec[key]
        text = f"{value:.6g}" if isinstance(value, float) else str(value)
        parts.append(f"{_slug_safe(key)}_{_slug_safe(text)}")
    return "-".join(parts)


def _loen_calc_tag(method: str, sub_spec: dict[str, Any]) -> str:
    """The single ``calc`` directory component of a loen unit (DESIGN
    6.2.4): ``loen-<method>-<sub_spec slug>``.  Encoding the method and
    every sub_spec key in the tag means two declarations that differ in
    any parameter land in different run directories by construction, so
    one loen run is never reused for the wrong descriptor."""

    return f"loen-{method}-{_sub_spec_slug(sub_spec)}"


def build_loen_units(ref: ReferenceSolid, struct_path: str,
                     options: dict[str, Any]) -> list:
    """Structure-only ``imago -loen -scf no`` units, one per distinct
    Fortran-side fingerprint declaration (PSEUDOCODE 11.4; DESIGN 5.10
    producer half).

    A loen-side (bispectrum) descriptor is computed by the Imago Fortran
    engine, so each such declaration needs its own dispatched
    ``-loen -scf no`` run.  The bispectrum is geometry-only, so these
    runs need no converged SCF -- they share the solid's structure and
    build options but override the job to ``loen`` and the SCF basis to
    ``no``, and add the ``-loeninput`` LOEN block that
    :func:`makegroups.loen_input_values` derives from the declaration's
    sub_spec (the same mapping the makegroups grouping flow uses, so the
    descriptors are comparable, DESIGN 5.10.5).

    One run serves every site that shares a (method, sub_spec): the
    descriptor table holds one row per atom, so declarations are
    de-duplicated by their calc tag and at most one unit is built per
    distinct tag.  Each unit is tagged ``kind="fingerprint"`` so the
    convergence harvest (:func:`pick_converged_unit`, which keeps only
    ``"convergence"``) ignores it and only the fingerprint harvest reads
    its descriptor.  Python-side (reduce) declarations need no unit --
    they are computed in process during the harvest -- so they are
    skipped here.

    ``options`` is the solid's ``make_producer_options`` dict; the loen
    overrides are layered on a copy so the convergence units are
    untouched."""

    units = []
    seen_tags: set[str] = set()
    for entry in ref.entries:
        for declaration in entry.fingerprints:
            matcher = MATCHERS[declaration.method]()
            if not matcher.needs_loen_run:
                continue            # Python-side: harvested in process.
            calc_tag = _loen_calc_tag(
                declaration.method, declaration.sub_spec)
            if calc_tag in seen_tags:
                continue            # One run already covers this sub_spec.
            seen_tags.add(calc_tag)

            loen_options = dict(options)
            loen_options["job"] = "loen"
            loen_options["scf_basis"] = "no"
            loen_options["loeninput"] = makegroups.loen_input_values(
                matcher, declaration.sub_spec)
            units.append(CalcUnit(
                id=ref.reference_id,
                structure=struct_path,
                options=loen_options,
                calc=(calc_tag,),
                kind="fingerprint",
                key_fields=standard_key_fields(
                    struct_path, loen_options)))
    return units


def _read_structure_with_distances(structure_path: str,
                                   cutoff: float) -> StructureControl:
    """Read the run's expanded full-cell structure and build its
    minimum-image distance matrix out to ``cutoff`` (DESIGN 5.7).

    ``structure_path`` is the run's ``outputs["structure"]`` -- the
    ``imago.fract-mi`` makeinput wrote, every atom explicit in space
    group 1 at the run's sorted (dat) numbering.  ``read_input_file``
    dispatches on the filename to the skeleton reader; the file is
    already a full P1 cell, so no further space-group expansion
    happens.  ``set_limit_dist`` sizes the periodic search to the
    reduce cutoff before ``create_min_dist_matrix`` populates
    ``sc.min_dist`` -- the same minimum-image geometry makeinput's
    grouping pass uses, so periodic boundary conditions enter exactly
    once and exactly here."""

    structure = StructureControl()
    structure.read_input_file(structure_path)
    structure.set_limit_dist(cutoff)
    structure.create_min_dist_matrix()
    return structure


def read_skeleton_to_dat_map(result_toml: dict
                             ) -> dict[int, tuple[int, str]]:
    """Read the run's ``datSkl.map`` into ``{skeleton_atom:
    (dat_atom, element)}`` (DESIGN 5.7 / ARCH 9.7).

    The reduce harvest indexes the expanded structure by the run's
    sorted (dat) numbering, while a manifest ``atom_site`` is a
    *skeleton* index, so this returns the skeleton-to-dat mapping
    (plus each site's element symbol, used to guard that the
    structure row and the map agree).  Columns are ``DAT#
    SKELETON#  ELEMENT  SPECIES  TYPE``; the header line is
    skipped."""

    path = result_toml["outputs"]["datSkl_map"]
    mapping: dict[int, tuple[int, str]] = {}
    with open(path) as handle:
        rows = [line for line in handle if line.strip()]
    for row in rows[1:]:                       # skip the header line
        columns = row.split()
        dat_atom = int(columns[0])
        skeleton_atom = int(columns[1])
        element = columns[2]
        mapping[skeleton_atom] = (dat_atom, element)
    return mapping


def harvest_fingerprints(flight: Flight, ref: ReferenceSolid,
                         atom_site: int, overrides: list,
                         result_toml: dict,
                         characterization: list) -> list:
    """Build the ``FingerprintRecord`` list for one environment
    (PSEUDOCODE 11.4 ``harvestFingerprints``).

    Every environment harvests the database-wide ``[characterization]``
    preferred recipe (one ``sub_spec`` per method, each
    ``preferred = true``) plus any rare per-entry override the
    customization added (extra ``preferred = false`` ``sub_spec``\\ s,
    DESIGN 5.7).  Each declaration already carries its ``preferred``
    flag -- ``True`` for a characterization record, ``False`` for an
    override -- so the two lists simply concatenate; a curator is
    expected to give overrides ``sub_spec``\\ s distinct from the
    recipe.  An environment with an empty recipe and no override
    harvests nothing and never reads the structure.

    Each declaration is dispatched by its matcher's family.  Python-side
    matchers (``reduce``) compute in process from the run's *expanded*
    full-cell structure (``outputs["structure"]``), which carries the
    geometry the shells need and the run's numbering; Fortran-side
    matchers (``bispectrum``) read the descriptor of the
    ``-loen -scf no`` unit ``build_loen_units`` already dispatched for
    this solid (:func:`harvest_loen_fingerprint`), which is why
    ``flight`` and ``ref`` are needed.  ``atom_site`` (a skeleton index)
    is mapped to the run's dat numbering through ``datSkl.map`` once and
    shared by both families.  The neighbor multiset is element-only, so
    the stored fingerprint transfers across structures (DESIGN 5.2)."""

    # The preferred recipe applies to every environment; the entry's
    #   own fingerprints (if any) ride along as non-preferred overrides.
    declarations = list(characterization) + list(overrides)
    if not declarations:
        return []

    # atom_site is a skeleton index; the run's structure and descriptor
    #   are both in sorted (dat) numbering, so resolve it once here and
    #   reuse the (dat row, element) for every declaration.  The element
    #   is the cross-check both families guard against.
    skeleton_to_dat = read_skeleton_to_dat_map(result_toml)
    dat_atom, map_element = skeleton_to_dat[atom_site]

    # The expanded structure is read only when a Python-side
    #   declaration needs it -- a loen-only entry never touches it -- and
    #   then only once, sized to the largest cutoff those declarations
    #   request (each matcher trims to its own sub_spec cutoff).
    python_decls = [declaration for declaration in declarations
                    if not MATCHERS[declaration.method]().needs_loen_run]
    structure = None
    if python_decls:
        cutoff = max(declaration.sub_spec["cutoff"]
                     for declaration in python_decls)
        structure = _read_structure_with_distances(
            result_toml["outputs"]["structure"], cutoff)
        # Guard the numbering assumption: the structure row and the map
        #   must name the same element, or the expansion and the map
        #   have desynced and the fingerprint would describe the wrong
        #   atom.  (The loen branch guards the same way against the
        #   descriptor's own identity column.)
        structure_element = structure.atom_element_name[dat_atom]
        if structure_element.lower() != map_element.lower():
            raise ValueError(
                f"site {atom_site}: datSkl.map names element "
                f"{map_element!r} but the expanded structure row "
                f"{dat_atom} is {structure_element!r}; numbering desync")

    records = []
    for declaration in declarations:
        matcher = MATCHERS[declaration.method]()
        if matcher.needs_loen_run:
            payload = harvest_loen_fingerprint(
                flight, ref, dat_atom, map_element, matcher,
                declaration.sub_spec)
        else:
            vectors = matcher.compute_query(
                structure, declaration.sub_spec)
            payload = matcher.build_payload(vectors[dat_atom])
        records.append(ipdb.FingerprintRecord(
            method=declaration.method,
            sub_spec=declaration.sub_spec,
            preferred=declaration.preferred,
            payload=payload))
    return records


def harvest_loen_fingerprint(flight: Flight, ref: ReferenceSolid,
                             dat_atom: int, element: str,
                             matcher, sub_spec: dict) -> dict:
    """Read one site's bispectrum payload from its loen unit's
    descriptor (PSEUDOCODE 11.4 ``harvestLoenFingerprint``).

    No engine run happens here: ``build_loen_units`` already dispatched
    the ``-loen -scf no`` unit for this ``(solid, method, sub_spec)``,
    and kaleidoscope's run-reuse cache (DESIGN 6.2.5) owns recompute
    avoidance.  The unit's run directory is reconstructed from the same
    calc tag the build used -- ``<root>/wingbeats/<reference_id>/
    loen-<method>-<slug>`` -- and the descriptor is the
    ``<edge>_loen<basis>.plot`` file makegroups' finder locates (DESIGN
    5.10.3, validated against a live run).

    The descriptor table has one row per atom in dat order; ``dat_atom``
    is the row this site maps to (the caller resolved it through
    ``datSkl.map``).  The row's self-describing ``element`` column is
    cross-checked against the map's element so a numbering desync fails
    loudly rather than storing the wrong atom's fingerprint.  The chosen
    vector is wrapped through ``matcher.build_payload`` so the on-disk
    payload field (``values`` for bispectrum) matches what the consumer
    reads (DESIGN 5.2)."""

    calc_tag = _loen_calc_tag(matcher.name, sub_spec)
    run_dir = os.path.join(
        flight.root, "wingbeats", ref.reference_id, calc_tag)
    descriptor_path = makegroups.find_loen_descriptor(run_dir)
    rows = matcher.parse_loen_output(descriptor_path, sub_spec)

    rows_by_site = {row.site: row for row in rows}
    if dat_atom not in rows_by_site:
        raise ValueError(
            f"loen descriptor {descriptor_path!r} has no row for dat "
            f"site {dat_atom} (it holds sites "
            f"{sorted(rows_by_site)}); the descriptor and the run's "
            f"numbering have desynced")
    row = rows_by_site[dat_atom]
    if row.element.lower() != element.lower():
        raise ValueError(
            f"dat site {dat_atom}: datSkl.map names element "
            f"{element!r} but the loen descriptor row is "
            f"{row.element!r}; numbering desync")
    return matcher.build_payload(row.vector)


# ============================================================
#  Harvest helpers: pick the converged run, read the potential
#  (DESIGN 5.7 / 7.8; PSEUDOCODE 11.4)
# ============================================================

def _read_unit_result(workspace_root: str, unit) -> dict:
    """Parse one dispatched unit's ``result.toml`` from its run
    directory (the ``kaleidoscope.workspace.unit_run_dir`` layout:
    ``<root>/wingbeats/<id>/<calc...>/result.toml``)."""

    path = os.path.join(workspace_root, "wingbeats", unit.id,
                        *unit.calc, "result.toml")
    with open(path, "rb") as handle:
        return tomllib.load(handle)


def _unit_completed(workspace_root: str, unit) -> bool:
    """True iff a dispatched unit ran to completion -- its
    ``status.toml`` records status ``"done"`` (the DESIGN 6.2.4
    lifecycle: queued / running / done / failed / lost).

    A failed / lost / still-queued unit, or one with no
    ``status.toml`` at all, did not complete and has written no
    ``result.toml`` to harvest.  Checking this BEFORE reaching for
    ``result.toml`` is the DESIGN 6.2.10 robustness fix: when every
    unit fails at the makeinput/imago seam, the old harvest opened a
    missing ``result.toml`` and crashed with ``FileNotFoundError``
    instead of recording a non-converged solid."""
    wingbeat_dir = os.path.join(
        workspace_root, "wingbeats", unit.id, *unit.calc)
    status = read_status(wingbeat_dir)
    return bool(status) and status.get("status") == "done"


def pick_converged_unit(flight: Flight, reference_id: str,
                        workspace_root: str, scf_threshold: float):
    """Return ``(unit, result_toml)`` for the converged grid point
    of one solid's convergence sweep, or ``None`` when the energy is
    still moving at the top of the grid -- or when nothing ran to
    completion (PSEUDOCODE 11.4).

    Uses the same two-sided delta-below-threshold rule as the
    guidance harvest (DESIGN 7.8 step 3c), reused via
    ``guidance_harvest.pick_converged``: the smallest interior
    k-density whose total energy is within ``scf_threshold`` of both
    neighbours.  A single-point grid (a trust-mode run or a
    single-point curator override) has no interior point to judge,
    so its lone run IS the deliverable -- the producer trusts the
    pinned/predicted point and harvests its potential directly.

    Only units that ran to completion are considered: a failed or
    result-less unit is treated as non-converged and dropped before
    its ``result.toml`` is opened (DESIGN 6.2.10).  When no unit
    completed -- e.g. every unit failed at the makeinput/imago seam
    -- the solid did not converge and the function returns ``None``,
    which the caller logs and skips (5.7)."""

    units = [unit for unit in flight.units
             if unit.id == reference_id
             and unit.kind == "convergence"]
    units = sorted(
        units,
        key=lambda unit: guidance_harvest.swept_value_of(
            unit, "kpt-density"))

    # Keep only the units that completed; a failed/result-less unit
    # is not a convergence data point and has no result.toml to read.
    completed = [unit for unit in units
                 if _unit_completed(workspace_root, unit)]
    if not completed:
        return None

    results = [_read_unit_result(workspace_root, unit)
               for unit in completed]

    if len(completed) == 1:
        return completed[0], results[0]

    energies = [result["total_energy"] for result in results]
    index = guidance_harvest.pick_converged(energies, scf_threshold)
    if index is None:
        return None
    return completed[index], results[index]


def _parse_scfv_type_block(path: str, type_number: int
                           ) -> tuple[list[float], list[float]]:
    """Parse the converged ``scfV`` output and return one potential
    type's ``(coefficients, alphas)`` (DESIGN 5.7 / ARCH 9.7).

    The file holds every OLCAO potential type, not one bare block::

        NUM_TYPES   <n>
        TOTAL__OR__SPIN_UP          <- total / spin-up channel
          <count_1>                 <- type 1: a count line ...
          <count_1 term lines>      <-   ... then that many terms
          <count_2> / <terms> ...   <- types 2..n follow in order
        SPIN_DN                     <- a redundant copy (non-spin)
          ...

    Each term line carries five columns; only column 1 (the
    coefficient) and column 2 (the alpha) are meaningful -- columns
    3-5 are the placeholder fields Imago ignores (as in ``coeff1``).
    The producer runs non-spin, so the ``TOTAL__OR__SPIN_UP``
    channel IS the total potential; the ``SPIN_DN`` channel and the
    trailing +UJ section are not read (spin handling is deferred).

    ``type_number`` is 1-based (the OLCAO potential-type number from
    ``datSkl.map``).  The block is found by walking the
    count-delimited blocks in order, so the per-type term counts
    need not be known ahead of time; each block's count line is
    cross-checked against the term lines that follow."""

    with open(path) as handle:
        lines = [line for line in handle if line.strip() != ""]

    header = lines[0].split()
    if not header or header[0] != "NUM_TYPES":
        raise ValueError(
            f"{path}: expected a 'NUM_TYPES' header, found "
            f"{lines[0].strip()!r}")
    num_types = int(header[1])
    if not 1 <= type_number <= num_types:
        raise ValueError(
            f"{path}: requested type {type_number} but the file "
            f"declares NUM_TYPES = {num_types}")

    channel_tag = lines[1].strip()
    if channel_tag != "TOTAL__OR__SPIN_UP":
        raise ValueError(
            f"{path}: expected the 'TOTAL__OR__SPIN_UP' channel "
            f"tag on line 2, found {channel_tag!r}")

    # Walk the count-delimited type blocks in order and keep the
    #   one for `type_number`.  Each block is a count line followed
    #   by that many five-column term lines.
    cursor = 2
    for current_type in range(1, num_types + 1):
        declared_count = int(lines[cursor].split()[0])
        cursor += 1
        term_lines = lines[cursor:cursor + declared_count]
        cursor += declared_count
        if len(term_lines) != declared_count:
            raise ValueError(
                f"{path}: type {current_type} count line says "
                f"{declared_count} terms but only "
                f"{len(term_lines)} follow")
        if current_type == type_number:
            coefficients = [float(line.split()[0])
                            for line in term_lines]
            alphas = [float(line.split()[1]) for line in term_lines]
            return coefficients, alphas

    # Unreachable: the range guard guarantees `type_number` lies in
    #   1..num_types, so the loop always returns first.
    raise AssertionError("scfV type block not found")


def extract_potential(result_toml: dict, atom_site: int
                      ) -> tuple[list[float], list[float]]:
    """Harvest the converged Gaussian potential for one atom site
    from a dispatched run's result (DESIGN 5.7 / ARCHITECTURE 9.7;
    PSEUDOCODE 11.4).

    The converged ``scfV`` output (``result.outputs["scfV"]``, the
    ``<edge>_scfV-<basis>.dat`` Imago writes from ``fort.8``) holds
    the self-consistent potential for EVERY OLCAO potential type in
    the material -- not one bare coefficient block.  Its layout is a
    ``NUM_TYPES`` header, then a ``TOTAL__OR__SPIN_UP`` channel that
    lists each type as a count line plus that many Gaussian-term
    lines, followed by a redundant ``SPIN_DN`` channel.  The
    producer runs non-spin, so the ``TOTAL__OR__SPIN_UP`` channel is
    the total potential and the ``SPIN_DN`` copy is ignored (spin
    handling is deferred).

    This selects the type block for ``atom_site``: the site's type
    number comes from the run's ``datSkl.map`` (DESIGN 5.2.1 / ARCH
    9.7, the same map the storage label is built from), and the
    block's columns 1 and 2 are taken together as the coefficient
    and its alpha.  Those alphas equal the basis input the producer
    fed makeinput -- the consistency "converged ``scfV`` matches
    input ``scfV``" (5.7) names.

    (Contrast ``_parse_coeff_file``, which parses a single
    headerless ``coeff1`` block -- one element's potential -- and is
    NOT interchangeable with this multi-type output.)"""

    scfv_path = result_toml["outputs"]["scfV"]
    # The site's potential-type number (1-based) is read from the
    #   run's datSkl.map -- the same map the storage label is built
    #   from (DESIGN 5.2.1 / ARCH 9.7) -- and selects the type block
    #   within the multi-type scfV output.
    type_number = read_site_identity_map(result_toml)[atom_site][2]
    return _parse_scfv_type_block(scfv_path, type_number)


def read_site_identity_map(result_toml: dict
                           ) -> dict[int, tuple[str, int, int]]:
    """Read the run's ``datSkl.map`` into ``{skeleton_atom:
    (element, species, type)}`` (DESIGN 5.2.1 / 5.7; ARCHITECTURE
    9.7).

    makeinput writes this file during input preparation, recording
    -- per atom -- the sorted-dat number, the original skeleton
    number, and the site's element symbol, OLCAO species number,
    and potential-type number (C87).  The producer keys harvest by
    ``atom_site``, which is a *skeleton* numbering index, so this
    returns the map keyed by the skeleton column (``values[1]``).

    The columns are ``DAT#  SKELETON#  ELEMENT  SPECIES  TYPE``; the
    header line is skipped.  Reading the whole file once per
    converged solid lets every site's label be assembled without a
    second pass."""

    path = result_toml["outputs"]["datSkl_map"]
    identity: dict[int, tuple[str, int, int]] = {}
    with open(path) as handle:
        rows = [line for line in handle if line.strip()]
    for row in rows[1:]:                       # skip the header line
        columns = row.split()
        skeleton_atom = int(columns[1])
        element = columns[2]
        species = int(columns[3])
        type_number = int(columns[4])
        identity[skeleton_atom] = (element, species, type_number)
    return identity


def assemble_entry_label(reference_id: str, element: str,
                         species: int, type_number: int,
                         atom_site: int) -> str:
    """Assemble the DESIGN 5.2.1 entry label
    ``<reference_id>-<element><species>-t<type>-a<site>``.

    The element symbol is lowercased so it fuses with the species
    number into the OLCAO species token the CLI speaks (``si1``),
    and the whole label is lowercase.  ``reference_id`` is already
    label-safe (manifest rule 5), so no further escaping is
    needed."""

    return (f"{reference_id}-{element.lower()}{species}"
            f"-t{type_number}-a{atom_site}")


def compose_auto_description(source_description: str | None,
                             reference_id: str, element: str,
                             species: int, atom_site: int) -> str:
    """Compose an entry's description when the customization omits one
    (DESIGN 5.7): the reference solid's ``source_description``,
    qualified by the harvested site's species and 1-based
    ``atom_site``, so an auto-discovered environment still carries
    readable prose rather than an empty field.  When the solid has no
    ``source_description`` either, fall back to the bare species/site
    identity plus the ``reference_id`` that produced it."""

    symbol = element.capitalize()
    qualifier = f"{symbol} species {species}, site {atom_site}"
    if source_description:
        base = source_description.rstrip().rstrip(".")
        return f"{base} ({qualifier})."
    return f"{qualifier}, from reference solid {reference_id}."


@dataclass(frozen=True)
class Environment:
    """One distinct local environment discovered in a converged
    reference run (PSEUDOCODE 11.4 ``discover_environments``).

    The per-element database stores *distinct environments, not
    atoms* (DESIGN 5.2.3), so the harvest visits one representative
    site per environment rather than every atom.  Within one
    crystalline run the assigning partition is the potential *type*:
    all atoms sharing an ``(element, species, type)`` are
    symmetry-equivalent and carry the same converged potential, so a
    single representative speaks for the whole group.

    Fields:
      ``element``     -- the site's element symbol.
      ``atom_site``   -- skeleton index of the representative site,
                         whose potential and fingerprints are
                         harvested.
      ``species``     -- OLCAO species number (used in the label).
      ``type_number`` -- potential-type number (used in the label).
      ``label``       -- the entry label: a curator customization's
                         explicit label, else the derived
                         ``<ref_id>-<elem><species>-t<type>-a<site>``
                         (DESIGN 5.2.1).
      ``default``     -- whether this is the element's default entry
                         (a customization's flag, else False).
      ``description`` -- the entry description: a customization's
                         text, else auto-composed from the solid's
                         ``source_description`` (DESIGN 5.7).
      ``overrides``   -- the customization's per-entry fingerprint
                         declarations (non-preferred overrides), or
                         an empty list for an auto-discovered
                         environment."""

    element: str
    atom_site: int
    species: int
    type_number: int
    label: str
    default: bool
    description: str
    overrides: list


def discover_environments(result_toml: dict, ref: ReferenceSolid,
                          *, identity_fn=read_site_identity_map
                          ) -> list[Environment]:
    """Yield one :class:`Environment` per distinct environment in a
    converged run (PSEUDOCODE 11.4; DESIGN 5.7 / 5.2.3).

    The run's site-identity map (``datSkl.map``) partitions every
    atom by ``(element, species, type)``.  Atoms sharing that key are
    equivalent under the run's assigning method (symmetry for a
    crystalline reference), so the harvest keeps one representative
    per group rather than one entry per atom.  The representative is
    chosen order-independently -- the lowest skeleton index in the
    group -- so the discovered set never depends on the map's row
    order (DESIGN 5.6.5).

    A manifest customization annotates the environment that contains
    its pinned ``atom_site``, supplying the curator's label, default
    flag, description, fingerprint overrides, and the representative
    site to harvest.  A *site-less* customization cannot yet be
    matched to an environment and is skipped (an interim limitation;
    matching by label or element is later work).

    ``identity_fn`` is injected so the orchestration can be tested
    with the ``datSkl.map`` reader mocked."""

    site_identity = identity_fn(result_toml)

    # Partition the run's atoms by environment key.  Sorting the keys
    #   below makes the discovered order deterministic regardless of
    #   the map's row order.
    sites_by_env: dict[tuple[str, int, int], list[int]] = {}
    for skeleton_atom, identity in site_identity.items():
        sites_by_env.setdefault(identity, []).append(skeleton_atom)

    # Match each site-pinned customization to the environment it
    #   annotates.  Two customizations on one environment is
    #   ambiguous (which label/default wins?), so it is a hard error.
    custom_by_env: dict[tuple[str, int, int], ReferenceEntry] = {}
    for spec in ref.entries:
        if spec.atom_site is None:
            continue        # site-less: not yet matchable (interim)
        if spec.atom_site not in site_identity:
            raise ValueError(
                f"{ref.reference_id}: customization pins site "
                f"{spec.atom_site}, which the converged run does "
                f"not contain")
        key = site_identity[spec.atom_site]
        if key in custom_by_env:
            raise ValueError(
                f"{ref.reference_id}: two customizations annotate "
                f"the same environment {key}; at most one is allowed")
        custom_by_env[key] = spec

    environments: list[Environment] = []
    for key in sorted(sites_by_env):
        site_element, species, type_number = key
        spec = custom_by_env.get(key)

        # The representative site: the curator's pinned site when a
        #   customization names one, else the order-independent
        #   lowest skeleton index in the group.
        if spec is not None and spec.atom_site is not None:
            atom_site = spec.atom_site
        else:
            atom_site = min(sites_by_env[key])

        # element: an explicit customization element is cross-checked
        #   against the site; omitted, it is the site's own element
        #   the run discovered (DESIGN 5.7 rule 3).
        if spec is not None and spec.element is not None:
            if spec.element.lower() != site_element.lower():
                raise ValueError(
                    f"{ref.reference_id}: customization names element "
                    f"{spec.element!r} but site {atom_site} is "
                    f"{site_element!r} in the converged run")
            element = spec.element
        else:
            element = site_element

        # label: the curator's explicit override, else derived from
        #   the run's site identity (DESIGN 5.2.1).
        if spec is not None and spec.label is not None:
            label = spec.label
        else:
            label = assemble_entry_label(
                ref.reference_id, element, species,
                type_number, atom_site)

        # description: the explicit override, else auto-composed from
        #   the solid's source_description (DESIGN 5.7).
        if spec is not None and spec.description is not None:
            description = spec.description
        else:
            description = compose_auto_description(
                ref.source_description, ref.reference_id,
                element, species, atom_site)

        default = spec.default if spec is not None else False
        overrides = (list(spec.fingerprints)
                     if spec is not None else [])

        environments.append(Environment(
            element=element, atom_site=atom_site, species=species,
            type_number=type_number, label=label, default=default,
            description=description, overrides=overrides))

    return environments


def make_imago_provenance(commit: str, timestamp: str,
                          ref: ReferenceSolid, atom_site: int,
                          scf_iterations) -> dict[str, Any]:
    """The ``[potential.provenance]`` block for a harvested
    Imago-source entry (DESIGN 5.2 / 5.7).

    Carries the ``source = "Imago"`` discriminant and every field
    ``initial_potential_db.require_provenance`` demands of an Imago
    entry (``reference_id``, ``atom_site``, ``kpoint_spec``,
    ``scf_threshold``, ``scf_iterations``) so the 5.8 validation
    harness can re-run the originating SCF, plus ``system_type``
    recorded for forensics (5.7 rule 2)."""

    return {
        "source": "Imago",
        "commit": commit,
        "generated_at": timestamp,
        "reference_id": ref.reference_id,
        "system_type": ref.system_type,
        "atom_site": atom_site,
        "kpoint_spec": dict(ref.kpoint_spec),
        "scf_threshold": ref.scf_threshold,
        "scf_iterations": scf_iterations,
    }


# ============================================================
#  Run log (DESIGN 5.7 / 5.8; PSEUDOCODE 11.4 write_run_log)
# ============================================================

def make_run_log_entry(ref: ReferenceSolid, unit,
                       result_toml: dict) -> dict[str, Any]:
    """One converged-solid row for the run log: the reference id,
    the converged k-density (read off the chosen unit's calc tag),
    and the SCF iteration count the 5.8 harness reads."""

    return {
        "reference_id": ref.reference_id,
        "converged": True,
        "converged_kpoint_density": guidance_harvest.swept_value_of(
            unit, "kpt-density"),
        "scf_iterations": result_toml.get("scf_iterations"),
    }


def make_nonconverged_log_entry(ref: ReferenceSolid
                                ) -> dict[str, Any]:
    """One non-converged-solid row for the run log: the sweep never
    flattened, so no potential was harvested and the curator must
    widen the grid (DESIGN 7.9)."""

    return {"reference_id": ref.reference_id, "converged": False}


def write_run_log(path: str, imago_commit: str, timestamp: str,
                  per_run_log: list[dict[str, Any]]) -> None:
    """Write the producer's run log (PSEUDOCODE 11.4): a manifest
    snapshot header (the Imago commit + timestamp) followed by one
    ``[[run]]`` block per reference solid.  The 5.8 validation
    harness reads this to know which solids converged and in how
    many SCF iterations.  Emitted with the kaleidoscope
    ``toml_line`` helper so scalars/arrays format consistently."""

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as handle:
        handle.write(toml_line("schema_version", 1))
        handle.write(toml_line("imago_commit", imago_commit))
        handle.write(toml_line("generated_at", timestamp))
        for entry in per_run_log:
            handle.write("\n[[run]]\n")
            for key, value in entry.items():
                if value is not None:
                    handle.write(toml_line(key, value))


# ============================================================
#  The producer pipeline (DESIGN 5.7; PSEUDOCODE 11.4)
# ============================================================

def build_initial_potentials(manifest_path: str, pdb_root: str,
                             data_root: str, *, force: bool = False,
                             single_element: str | None = None,
                             dispatch_shape: str = "local",
                             partition: str | None = None,
                             nodes: int | None = None,
                             walltime: str | None = None,
                             profile: str | None = None,
                             save_config: bool = False,
                             dispatch_fn=dispatch,
                             extract_fn=extract_potential,
                             identity_fn=read_site_identity_map,
                             fingerprint_fn=harvest_fingerprints
                             ) -> list[dict[str, Any]]:
    """The three-phase producer (DESIGN 5.7; PSEUDOCODE 11.4):
    *build* a single combined flight, *dispatch* it through
    kaleidoscope, then *harvest* each solid's converged potential
    and contribute the same converged point back to the guidance
    dataspace.  Returns the per-run log (also written to disk).

    ``dispatch_shape`` selects where the flight runs (DESIGN 6.2.11):
    ``local`` (the default here, so tests and in-process callers stay
    serial without any cluster settings file) runs every unit in
    process; ``slurm-pooled`` and ``slurm-per-job`` build a Parsl
    ``Config`` from the per-site settings file and the per-run
    ``partition`` / ``nodes`` / ``walltime`` / ``profile`` choices.
    ``save_config`` records the resolved cluster choices beside the
    run.  ``force`` bypasses the run-reuse cache (DESIGN 6.2.5).

    ``dispatch_fn``, ``extract_fn``, ``identity_fn``, and
    ``fingerprint_fn`` are injected (defaulting to the real
    kaleidoscope dispatch, the real scfV reader, the real
    ``datSkl.map`` reader, and the real fingerprint harvest) so the
    orchestration can be unit-tested with the toolchain seam mocked:
    end-to-end dispatch, the per-site ``scfV`` read, the site-identity
    read, and the fingerprint harvest (which needs the run's expanded
    structure and loen descriptor) all need a live Imago run (C74)."""

    # Pass the matcher registry so the loader enforces rule 9: every
    #   declared fingerprint method must be a registered matcher (C54).
    manifest = load_manifest_v2(
        manifest_path, known_methods=set(MATCHERS))
    manifest_dir = os.path.dirname(manifest.manifest_path)
    guidance_root = os.path.join(data_root, "historicalGuidanceDB")
    dataspace = guidance_db.load(guidance_root)
    imago_commit = _git_sha()
    timestamp = _now_iso8601_utc()
    workspace = curation_workspace_root(pdb_root)

    # ----- Phase 1: build.  Refresh the isolated baselines, then
    # accumulate every solid's convergence units (and any loen
    # units) into ONE combined flight so the whole producer run
    # dispatches as a single flat parallel batch.
    elements = None if single_element is None else [single_element]
    databases = refresh_isolated_entries(
        pdb_root, manifest, imago_commit, timestamp, elements)

    all_units: list = []
    predictions: dict[str, Any] = {}
    for ref in manifest.reference_solids:
        struct = materialize_structure(ref, manifest_dir, pdb_root)
        options = make_producer_options(ref, imago_commit)
        # The predictor and the PredictionRecord speak the human
        # physics names, not the codes, so the sub-model travels to
        # the builder in its OWN dict -- never mixed into the
        # tool-facing options (DESIGN 6.2.8 / 6.2.10), which would
        # both duplicate the basis and make makeinput reject
        # "functional" / "kpoint_integration".
        submodel = {
            "basis": ref.basis,
            "functional": ref.functional,
            "kpoint_integration": ref.kpoint_integration,
        }
        # One builder call per solid (DESIGN 6.2.9): a pinned
        # kpoint_spec.density is the curator override (predictor
        # bypassed); otherwise full predict-then-verify.
        flight_i, record_i = build_kpoint_convergence(
            struct, options, dataspace, ref.system_type, submodel,
            id=ref.reference_id,
            center=ref.kpoint_spec.get("density"))
        all_units.extend(flight_i.units)
        all_units.extend(build_loen_units(ref, struct, options))
        # Store the plain dict (metadata must be TOML-serializable).
        predictions[ref.reference_id] = asdict(record_i)

    flight = Flight(
        root=workspace, units=all_units,
        sweep=SweepRecord(varied_axes=("kpt-density",),
                          fixed_axes={}),
        metadata={"predictions": predictions})

    # ----- Phase 2: dispatch.  Kaleidoscope runs and tracks every
    # unit through the wingbeat seam and its run-reuse cache.  The
    # dispatch shape becomes the flight's Parsl Config (None for the
    # local opt-out); the driver then auto-selects Local vs Parsl and
    # honours the force cache-bypass (DESIGN 6.2.11; PSEUDOCODE 13.5).
    flight.parsl_config, dispatch_choices = resolve_dispatch(
        dispatch_shape, partition, nodes, walltime, profile)
    if save_config and dispatch_choices is not None:
        write_resolved_dispatch(workspace, dispatch_choices, profile)
    dispatch_fn(flight, force=force)

    # ----- Phase 3: harvest.  Per solid, pick the converged grid
    # point, extract the potential at each named site, and record
    # the run.  Non-converged solids are logged and skipped.
    per_run_log: list[dict[str, Any]] = []
    for ref in manifest.reference_solids:
        converged = pick_converged_unit(
            flight, ref.reference_id, workspace, ref.scf_threshold)
        if converged is None:
            per_run_log.append(make_nonconverged_log_entry(ref))
            continue
        unit, result_toml = converged
        per_run_log.append(
            make_run_log_entry(ref, unit, result_toml))
        scf_iterations = result_toml.get("scf_iterations")

        # Discover one representative per distinct environment in the
        #   converged run (DESIGN 5.2.3 / 5.7).  Each carries its
        #   resolved label, default flag, description, and fingerprint
        #   overrides, with any site-pinned customization already
        #   layered on; the harvest only extracts the potential and
        #   inserts the entry.
        for env in discover_environments(
                result_toml, ref, identity_fn=identity_fn):
            elem_key = env.element.lower()
            if elem_key not in databases:
                continue        # filtered out by --element
            coefficients, alphas = extract_fn(
                result_toml, env.atom_site)

            new_entry = ipdb.PotentialEntry(
                label=env.label,
                default=env.default,
                description=env.description,
                num_gaussians=len(coefficients),
                alpha_min=min(alphas),
                alpha_max=max(alphas),
                coefficients=coefficients,
                # The entry keeps this representative's harvested
                # potential; a zero spread fills the column the
                # current schema still carries.  The insert-or-skip of
                # the B phase (C88) keeps the first representative on a
                # duplicate (DESIGN 5.2.3).
                coefficient_std=[0.0] * len(coefficients),
                alphas=alphas,
                provenance=make_imago_provenance(
                    imago_commit, timestamp, ref, env.atom_site,
                    scf_iterations),
                # Every environment harvests the database-wide
                #   [characterization] recipe (preferred) plus any
                #   per-entry override (non-preferred).
                fingerprints=fingerprint_fn(
                    flight, ref, env.atom_site, env.overrides,
                    result_toml, manifest.characterization))
            database = databases[elem_key]
            # INTERIM (C88): cross-run bispectrum dedup (insert-or-
            #   skip) is not built yet, so replace-by-label collapses
            #   only a same-run label collision; a new label appends.
            database.potentials = [
                entry for entry in database.potentials
                if entry.label != env.label]
            database.potentials.append(new_entry)

    # ----- Phase 3b: guidance contribution.  The same converged
    # grid points feed the historical-guidance dataspace staging,
    # so every solid the producer converges sharpens the predictor.
    guidance_harvest.harvest_flight(
        workspace, guidance_root, dataspace)

    # ----- Write outputs: every affected element file, plus the
    # run log the 5.8 validation harness reads.
    save_databases(databases, pdb_root)
    write_run_log(
        os.path.join(data_root, "curation", "run_log.toml"),
        imago_commit, timestamp, per_run_log)
    return per_run_log


# ============================================================
#  Command-line interface
# ============================================================

def _default_pdb_root() -> str:
    """``$IMAGO_DATA/atomicPDB`` (DESIGN 5.4 layout), or empty when
    $IMAGO_DATA is unset so the parser can demand ``--pdb-root``."""

    data_dir = os.environ.get("IMAGO_DATA", "")
    return os.path.join(data_dir, "atomicPDB") if data_dir else ""


def _print_materialize_report(report: list[dict[str, Any]]) -> None:
    """Print the per-solid pre-flight result -- one line per solid
    plus a final tally -- in the style of the producer's summary."""

    for row in report:
        mark = "ok  " if row["ok"] else "FAIL"
        print(f"  [{mark}] {row['reference_id']}: {row['source']}")
        if row["ok"]:
            print(f"          -> {row['skl_path']}")
        else:
            print(f"          {row['message']}")
    ready = sum(1 for row in report if row["ok"])
    print(f"materialize: {ready}/{len(report)} reference structures "
          f"fetched and converted")


def main(argv=None) -> int:
    """CLI entry point: run the producer over a curation manifest
    (DESIGN 5.7).  ``--element`` restricts the run to one element's
    database; ``--force`` bypasses kaleidoscope's run-reuse cache so
    every reference run re-executes.  ``--dispatch`` chooses where the
    flight runs -- one scheduler job per unit by default, ``local``
    to run in process without a cluster settings file -- with
    ``--partition`` / ``--nodes`` / ``--walltime`` / ``--profile``
    tuning a cluster dispatch and ``--save-config`` recording the
    resolved choices (DESIGN 6.2.11).  ``--materialize-only`` runs
    just the structure fetch-and-convert pre-flight (no SCF),
    optionally redirecting the output with ``--materialize-dir``."""

    parser = argparse.ArgumentParser(
        description="Build the augmented initial-potential database "
                    "from a curation manifest.")
    parser.add_argument(
        "--manifest", default="manifest.toml",
        help="path to the curation manifest, schema v2 (default: "
             "manifest.toml in the current working directory)")
    parser.add_argument(
        "--pdb-root", default=_default_pdb_root(),
        help="the atomicPDB root (default: $IMAGO_DATA/atomicPDB)")
    parser.add_argument(
        "--element", default=None,
        help="restrict the build to this one element's database "
             "(default: build every element named in the manifest)")
    parser.add_argument(
        "--force", action="store_true",
        help="bypass the run-reuse cache so every reference run "
             "re-executes (default: reuse any cached runs)")
    parser.add_argument(
        "--dispatch", default="slurm-per-job",
        choices=["local", "slurm-pooled", "slurm-per-job"],
        help="where to run the flight: 'local' runs every unit in "
             "process and needs no cluster settings file; "
             "'slurm-pooled' streams units through one shared "
             "allocation; 'slurm-per-job' submits one scheduler job "
             "per unit (default: slurm-per-job)")
    parser.add_argument(
        "--partition", default=None,
        help="scheduler queue for a cluster dispatch (default: the "
             "first partition in the cluster settings file)")
    parser.add_argument(
        "--nodes", type=int, default=None,
        help="nodes per allocation for a cluster dispatch (default: "
             "the cluster settings file's nodes value)")
    parser.add_argument(
        "--walltime", default=None,
        help="time limit for a cluster dispatch, e.g. 02:00:00 "
             "(default: the cluster settings file's walltime value)")
    parser.add_argument(
        "--profile", default=None,
        help="named cluster profile to select from the settings "
             "file (default: the file's base settings)")
    parser.add_argument(
        "--save-config", action="store_true",
        help="write the resolved cluster dispatch choices beside the "
             "run for a reproducible record (default: do not write "
             "the record)")
    parser.add_argument(
        "--materialize-only", action="store_true",
        help="fetch and convert every reference structure named in "
             "the manifest, then stop without running any SCF -- a "
             "pre-flight to validate a freshly pinned structure set "
             "(default: run the full build)")
    parser.add_argument(
        "--materialize-dir", default=None,
        help="directory the --materialize-only pre-flight writes "
             "fetched CIFs and converted skeletons into (default: "
             "the shared structure cache beside the databases, which "
             "a later full build reuses)")
    args = parser.parse_args(argv)

    if not args.pdb_root:
        parser.error("--pdb-root not given and $IMAGO_DATA is unset")
    if args.materialize_dir and not args.materialize_only:
        parser.error("--materialize-dir applies only with "
                     "--materialize-only")
    data_root = os.path.dirname(args.pdb_root.rstrip("/"))

    # The pre-flight short-circuits before any SCF dispatch: it only
    #   needs the structure sources, so it never touches the run or
    #   harvest fields the full build requires.
    if args.materialize_only:
        report = materialize_only(
            args.manifest, args.pdb_root,
            cache_dir=args.materialize_dir)
        _print_materialize_report(report)
        return 0 if all(row["ok"] for row in report) else 1

    per_run_log = build_initial_potentials(
        args.manifest, args.pdb_root, data_root,
        force=args.force, single_element=args.element,
        dispatch_shape=args.dispatch, partition=args.partition,
        nodes=args.nodes, walltime=args.walltime,
        profile=args.profile, save_config=args.save_config)
    converged = sum(1 for row in per_run_log if row["converged"])
    print(f"producer: {converged}/{len(per_run_log)} reference "
          f"solids converged and harvested")
    return 0


def record_command():
    """Append the issued command line to a file named "command" in
    the current directory, so the exact invocation can be recovered
    later.  This is a standing project convention: each run appends
    a dated block, so the file builds up a history of how the script
    was called."""

    with open("command", "a") as cmd:
        now = datetime.now()
        stamp = now.strftime("%b. %d, %Y: %H:%M:%S")
        cmd.write(f"Date: {stamp}\n")
        cmd.write("Cmnd:")
        for argument in sys.argv:
            cmd.write(f" {argument}")
        cmd.write("\n\n")


if __name__ == "__main__":
    record_command()
    sys.exit(main())
