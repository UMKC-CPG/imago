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
import shlex
import subprocess
import sys
import tomllib
from collections import namedtuple
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any

import initial_potential_db as ipdb
import guidance_db
import guidance_harvest
import mesh_climb
# The curation-manifest schema (dataclasses + readers) lives in its
#   own neutral library so the producer and the authoring tool
#   (expand_manifest) share one definition; the producer imports the
#   names it runs a manifest with (DESIGN 5.7).
from curation_manifest import (
    ReferenceEntry, ReferenceSolid, CurationManifest,
    load_manifest_v2, load_structure_sources, resolve_settings)
from kaleidoscope import (
    CalcUnit, Flight, SweepRecord, dispatch, send_off, collect_next,
    make_executor)
from kaleidoscope.builders.kpoint_convergence import (
    build_mesh_unit, predict_kpoint_density, standard_key_fields)
from kaleidoscope.cluster_config import (
    resolve_dispatch, write_resolved_dispatch,
    load_site_config, resolve_choices, resolve_orchestrator,
    build_orchestrator_sbatch)
from kaleidoscope.workspace import toml_line
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
        # The isolated baseline stores the single atomSCF
        # potential verbatim, like any other environment
        # (DESIGN 5.2.3).
        coefficients=coefficients,
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
    """Step 1 of the pipeline: load each element database and
    refresh only its isolated baseline (the incremental model,
    DESIGN 5.7).

    For each element (all element directories under ``pdb_root``,
    or just ``elements`` when given), load the existing
    ``s_gaussian_pot.toml`` if present -- **preserving every
    environment harvested by earlier runs** -- or create a fresh
    :class:`initial_potential_db.ElementDatabase` from the
    element's ``pot1`` scalars.  Then refresh only the ``"isolated"``
    entry from the current pot1/coeff1 (so atomSCF changes
    propagate), leaving all harvested entries untouched.  The harvest
    phase inserts-or-skips this run's solids on top, so the database
    grows incrementally and re-running an unchanged manifest moves
    nothing (DESIGN 5.2.3).  Returns the in-memory databases keyed by
    element directory name; the caller saves them (see
    :func:`save_databases`).

    Implements PSEUDOCODE 11.4 step 1.  The existing-file load passes
    ``known_methods=None`` (rule 9 is skipped) because the matcher
    registry is not threaded in here; the harvested fingerprints were
    validated when the file was first written.
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

        # INCREMENTAL (DESIGN 5.7): keep every harvested entry the
        # loaded file already holds -- there is no reset -- and
        # refresh only the "isolated" baseline from the current
        # pot1/coeff1 (so atomSCF changes propagate): drop the old
        # isolated entry by label and append the fresh one.  The
        # harvest phase then inserts-or-skips this run's solids on
        # top; re-running an unchanged manifest moves nothing because
        # skip-on-match drops every duplicate (DESIGN 5.2.3).
        fresh_isolated = build_isolated_entry(
            pdb_root, elem, commit, timestamp, manifest)
        database.potentials = [
            entry for entry in database.potentials
            if entry.label != "isolated"]
        database.potentials.append(fresh_isolated)
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
                     options: dict[str, Any],
                     characterization: list) -> list:
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

    The declarations built here must be exactly the set the harvest will
    read, or a harvested fingerprint would find no dispatched run.  The
    harvest (:func:`harvest_fingerprints`) reads the database-wide
    ``[characterization]`` recipe PLUS each entry's own overrides, so
    this build side unions the same two sources: the ``characterization``
    list (the preferred recipe, applied to every environment) and every
    entry's ``fingerprints`` (rare per-entry overrides).  A Si default
    manifest puts its bispectrum in ``[characterization]`` with no
    per-entry overrides, so omitting the recipe here would build no loen
    unit at all and the harvest would fail on a missing descriptor.

    One run serves every site that shares a (method, sub_spec): the
    descriptor table holds one row per atom, so declarations across both
    sources are de-duplicated by their calc tag and at most one unit is
    built per distinct tag.  Each unit is tagged ``kind="fingerprint"``
    so the producer dispatches it in the separate loen pre-flight batch
    (never in a climb convergence round), and only the fingerprint
    harvest reads its descriptor.  Python-side (reduce) declarations
    need no unit -- they are computed in process during the harvest --
    so they are skipped here.

    ``options`` is the solid's ``make_producer_options`` dict; the loen
    overrides are layered on a copy so the convergence units are
    untouched."""

    # Mirror the harvest's declaration stream (recipe first, then each
    #   entry's overrides) so every Fortran-side fingerprint the harvest
    #   will read has a dispatched run; the calc-tag dedup below collapses
    #   any (method, sub_spec) shared across the two sources.
    declarations = list(characterization)
    for entry in ref.entries:
        declarations.extend(entry.fingerprints)

    units = []
    seen_tags: set[str] = set()
    for declaration in declarations:
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


def _preferred_bispectrum(entry: ipdb.PotentialEntry):
    """Return ``entry``'s preferred bispectrum fingerprint record,
    or None when it carries no bispectrum descriptor.

    The dedup keys on the database-wide ``[characterization]``
    bispectrum -- the transferable descriptor every harvested entry
    shares (DESIGN 5.2.3).  Per-entry overrides are non-preferred
    and never key the dedup, and the isolated baseline carries no
    fingerprints at all, so this returns None for it."""

    for fingerprint in entry.fingerprints:
        if fingerprint.method == "bispectrum" and fingerprint.preferred:
            return fingerprint
    return None


def find_bispectrum_duplicate(database: ipdb.ElementDatabase,
                              new_entry: ipdb.PotentialEntry,
                              similarity_floor: float | None = None
                              ) -> ipdb.PotentialEntry | None:
    """Return the stored entry whose bispectrum descriptor matches
    ``new_entry``'s within the similarity floor, or None (DESIGN
    5.2.3).

    The dedup keys on the preferred bispectrum descriptor at its
    ``sub_spec`` -- the transferable one every harvested entry
    carries.  An entry with no preferred bispectrum record (the
    isolated baseline, or a recipe with no bispectrum) has no key,
    so it never matches and is always treated as novel.  The
    comparison is the same L2 match the consumer uses
    (``BispecMatcher.match_distance``, DESIGN 5.6.5), and the floor
    defaults to the matcher's ``default_similarity_floor`` -- the
    producer-side mirror of the consumer's floor (DESIGN 5.2.3,
    TODO C61).  When several stored entries match, the nearest is
    returned."""

    matcher = MATCHERS["bispectrum"]()
    new_record = _preferred_bispectrum(new_entry)
    if new_record is None:
        return None
    new_vector = matcher.extract_query_vector(new_record.payload)
    floor = (similarity_floor if similarity_floor is not None
             else matcher.default_similarity_floor)

    nearest_entry = None
    nearest_distance = None
    for entry in database.potentials:
        try:
            stored = ipdb.find_fingerprint(
                entry, "bispectrum", new_record.sub_spec)
        except KeyError:
            continue        # not comparable at this sub_spec
        distance = matcher.match_distance(
            new_vector, matcher.extract_query_vector(stored.payload))
        if distance <= floor and (nearest_distance is None
                                  or distance < nearest_distance):
            nearest_entry, nearest_distance = entry, distance
    return nearest_entry


def insert_or_skip(database: ipdb.ElementDatabase,
                   new_entry: ipdb.PotentialEntry,
                   similarity_floor: float | None = None) -> None:
    """Insert ``new_entry`` into ``database`` or skip it as a
    duplicate (DESIGN 5.2.3 insert-or-skip).

    An entry whose label already exists is replaced in place: that
    is a re-harvested solid (same derived label) or a curator
    customization whose explicit label names an entry to override.
    Otherwise the bispectrum dedup decides: a new environment whose
    descriptor matches a stored one within the similarity floor is
    SKIPPED (the first representative's potential stands), and a
    genuinely novel environment is appended.  The leaner model
    stores nothing extra on a skip -- no counts, no merge;
    reconciling duplicates into a statistical mean is the deferred
    C103 upgrade (DESIGN 5.2.3)."""

    for index, entry in enumerate(database.potentials):
        if entry.label == new_entry.label:
            database.potentials[index] = new_entry
            return
    if find_bispectrum_duplicate(
            database, new_entry, similarity_floor) is None:
        database.potentials.append(new_entry)
    # else: a stored entry already covers this environment; skip it.


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
    recorded for forensics (5.7 rule 2).

    ``type_assignment`` names the scheme that drew the run's
    species/type partition; from it each fingerprint's native vs
    witness role is derived (DESIGN 5.2.2: method M is native iff
    M == type_assignment).  The producer assigns types
    crystallographically for the references it currently handles,
    so every harvested entry is ``"symmetry"``-assigned and both
    its reduce and bispectrum fingerprints are (exact) witnesses;
    the bispectrum-assigned path for disordered references is
    deferred."""

    return {
        "source": "Imago",
        "commit": commit,
        "generated_at": timestamp,
        "reference_id": ref.reference_id,
        "system_type": ref.system_type,
        "type_assignment": "symmetry",
        "atom_site": atom_site,
        "kpoint_spec": dict(ref.kpoint_spec),
        "scf_threshold": ref.scf_threshold,
        "scf_iterations": scf_iterations,
    }


# ============================================================
#  Run log (DESIGN 5.7 / 5.8; PSEUDOCODE 11.4 write_run_log)
# ============================================================

def make_run_log_entry(ref: ReferenceSolid, harvest_inputs: dict,
                       result_toml: dict) -> dict[str, Any]:
    """One converged-solid row for the run log: the reference id,
    the converged mesh AND its k-density, and the SCF iteration
    count the 5.8 harness reads (PSEUDOCODE 11.4).

    ``harvest_inputs`` is :func:`record_converged`'s output for this
    solid, so the mesh and density are the exact converged rung's --
    the climb searches in mesh space, and a mesh does not round-trip
    from its calc tag the way a swept density did, so both are read
    from the recorded rung rather than the unit's tag (DESIGN
    3.12.4)."""

    return {
        "reference_id": ref.reference_id,
        "converged": True,
        "converged_mesh": harvest_inputs["converged_mesh"],
        "converged_kpoint_density":
            harvest_inputs["converged_kpoint_density"],
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

def apply_manifest_defaults(manifest: CurationManifest) -> None:
    """Fold the top-level ``[defaults]`` and ``[harvest]`` blocks into
    each reference solid, in place, so the rest of the producer reads
    one fully resolved value per field -- the five run settings
    (``basis``, ``functional``, ``kpoint_integration``,
    ``kpoint_spec``, ``scf_threshold``) and the harvest setting
    (``kpoint_convergence_threshold``) -- and never has to consult the
    shared blocks again (DESIGN 5.7).

    A solid that names its own value keeps it; a solid that omits a
    run setting inherits the ``[defaults]`` value, and the harvest
    setting inherits ``[harvest]`` or the built-in default.  The
    loader already proved every run setting is resolvable for every
    solid (manifest rule 2), so after this pass no field is left
    ``None`` -- the downstream ``refresh_isolated_entries`` call and
    the two per-solid loops can treat every setting as present.
    """

    manifest.reference_solids = [
        resolve_settings(solid, manifest.defaults, manifest.harvest)
        for solid in manifest.reference_solids]


def prepare_units(flight: Flight, workspace: str, units=None) -> None:
    """Build each unit's staged inputs in the DRIVER, before
    dispatch (DESIGN 6.2.5 "Prepare before the hit-test";
    PSEUDOCODE 11.4 Phase 1b).

    ``units`` selects which of the flight's units to prepare; it
    defaults to every unit (a one-shot flight) but the climb passes
    just the rungs it has newly decided (PSEUDOCODE 4e.7), so an
    accreting flight re-prepares nothing it already staged.

    makeinput's resolved output, ``structure.dat``, is the cache
    key file (:func:`standard_key_fields`), so it must exist before
    kaleidoscope's hit-test can read it.  Each unit is built into
    its own per-unit staging directory under
    ``<workspace>/prepare/<id>/<calc...>`` -- deliberately SEPARATE
    from the run directory, so a prior run's staged ``structure.dat``
    is not clobbered before the byte-compare (the "must not clobber"
    rule).  The unit's ``prepared_dir`` records the staging directory
    (the wingbeat commits it into the run directory on a cache miss),
    and the ``structure.dat`` KeyFile's ``source`` is re-pointed at
    the freshly built copy so the hit-test compares the right file.

    Running makeinput here, in the driver, is also what lets a cache
    hit be decided from local files: a hit never reaches the
    scheduler, so the surviving misses are the only units that cost
    a calculation.
    """

    # The wingbeat's option split is the single source of truth for
    #   which keys build the inputs (makeinput) versus ride as imago
    #   run-time settings; the prepare pass, which replaces the
    #   wingbeat's own build, reuses it so the two cannot drift.
    from kaleidoscope.wingbeats import _partition_options
    import makeinput

    targets = flight.units if units is None else units
    for unit in targets:
        staging = os.path.join(
            workspace, "prepare", unit.id, *unit.calc)
        makeinput_options, _ = _partition_options(unit.options)
        makeinput.build_run_dir(
            unit.structure, makeinput_options, staging)
        unit.prepared_dir = staging
        # Re-point the structure.dat KeyFile source at the freshly
        #   built copy (standard_key_fields left it provisional).
        for key_file in unit.key_fields.files:
            if key_file.name == "structure.dat":
                key_file.source = os.path.join(
                    staging, "structure.dat")


# ==================================================================
#  The adaptive k-point mesh climb (DESIGN 3.12; PSEUDOCODE 4e.3 /
#  4e.5)
#
#  The producer converges each reference solid's k-point sampling by
#  CLIMBING through symmetry-compatible meshes rather than sweeping a
#  fixed density grid: it launches every solid's opening rung at once,
#  then collects rungs as they land, and for each landing decides that
#  one solid's next mesh and launches it -- so a solid climbs the
#  instant its own rung lands and no solid waits on another, until
#  every material goes flat (converged) or reaches a ceiling
#  (non-converged).  This is the iterative-client shape ARCHITECTURE
#  9.7 reserves for the producer -- the dispatch core stays dumb
#  (Principle 12), running whatever meshes it is handed; the loop that
#  reads energies and chooses the next mesh lives here.
#
#  The mesh arithmetic (rung moves, the ceiling, first-round seeding,
#  the confidence policy) lives in ``mesh_climb``; the energy
#  flatness test is ``guidance_harvest.pick_converged_climb``, single-
#  sourced with the harvest's own convergence rule so the two cannot
#  drift.  What lives here is the control loop tying them together.
#
#  How a chosen mesh becomes an actual dispatched run -- the climb
#  ``dispatcher`` with its ``send`` and ``next_rung`` -- is INJECTED
#  into ``converge_by_climb``, exactly as ``dispatch_fn`` /
#  ``prepare_fn`` are on ``build_initial_potentials`` below, so this
#  loop is unit-testable with a synthetic dispatcher.  The real
#  dispatcher (``make_climb_dispatcher``: the explicit-mesh CalcUnit,
#  its calc tag, the one accreting flight, and the mesh-keyed energy
#  read-back) is below.
# ==================================================================

# One rung of a climb: a resolved mesh (axial counts ``[a, b, c]``)
#   and the total-cell energy the engine returned for it (raw
#   hartree, the basis ``pick_converged_climb`` normalizes to eV per
#   atom before judging flatness, DESIGN 7.8).
Rung = namedtuple("Rung", ["mesh", "energy"])


## Everything one material's climb needs, gathered per material
##   (DESIGN 3.12).  Geometry comes from the loaded cell, the policy
##   fields from ``resolve_climb_policy`` (mesh_climb, keyed on the
##   predictor's confidence), and the energy/ceiling knobs from the
##   manifest characterisation block:
##     classes            axis-class labels (mesh_climb 4c.1 / 4c.7)
##     recip_mag          the three reciprocal-axis magnitudes |b_i|
##     recip_cell_volume  reciprocal cell volume (density<->mesh map)
##     mode               PARALLEL_GRID, BRACKET_REFINE, or UNIT_STEP
##                        (mesh_climb) -- how the ladder is sampled
##     flat_needed        consecutive flat rungs the stop test needs
##     grid_width         rungs each side of the seed (grid mode)
##     start_offset       rungs below the seed a climb starts from
##     max_stride         largest geometric stride the bracket phase
##                        may take (bracket-refine only, DESIGN
##                        3.12.3)
##     cell_atom_count    atoms per cell (per-atom energy normalizer)
##     threshold          per-atom eV flatness tolerance (DESIGN 7.8)
##     max_count          per-axis ceiling backstop (DESIGN 3.12.3)
ClimbConfig = namedtuple(
    "ClimbConfig",
    ["classes", "recip_mag", "recip_cell_volume", "mode",
     "flat_needed", "grid_width", "start_offset", "max_stride",
     "cell_atom_count", "threshold", "max_count"])


# The three verdicts a climb step returns (DESIGN 3.12.3 / 3.12.5):
#   run one more mesh, stop converged at a rung, or stop at the
#   ceiling having never gone flat (non-converged).  Every search
#   shape -- the unit-step climb (``climb_action``), the bracket-
#   refine climb (``bracket_refine_next``), and the grid continuation
#   -- returns one of these through ``climb_next`` (PSEUDOCODE 4e.3).
_ACTION_RUN = "run"
_ACTION_CONVERGED = "converged"
_ACTION_CEILING = "ceiling"

# A climb-step result.  ``kind`` is one of the three verdicts above;
#   ``rung`` is the converged ``Rung`` (set for CONVERGED only);
#   ``mesh`` is the next mesh to run (set for RUN only).  CONVERGED
#   carries the rung itself, not an index, so the caller need not
#   track which ladder the index points into -- the bracket-refine
#   climb judges a filled sub-block, whose indices differ from the
#   full ladder's (PSEUDOCODE 4e.3).  The fields a verdict does not
#   use are ``None``.
ClimbAction = namedtuple("ClimbAction", ["kind", "rung", "mesh"])

# The outcome a material carries when it stops at the ceiling without
#   ever going flat (DESIGN 3.12.3).  A distinct sentinel, not a
#   ``Rung``, so the harvest tells a non-converged material apart
#   from a converged one at a glance.
NON_CONVERGED = "non_converged"

# The marker the climb dispatcher's ``next_rung`` returns in place of
#   a ``Rung`` when a requested mesh did not run to completion
#   (PSEUDOCODE 4e.5 / 4e.7; DESIGN 7.7).  It is distinct from
#   ``NON_CONVERGED`` (a whole material's verdict): this marks one
#   failed rung, which the climb loop turns into a run-failure stop.
_RUN_FAILED = "run_failed"


def _mesh_point_count(mesh):
    """The full-mesh k-point count of an axial-count mesh -- the
    product of its three counts.  This orders a climb: a rung with
    more points sits higher, so it is the sort key for the ladder."""
    return mesh[0] * mesh[1] * mesh[2]


def _sort_by_mesh(rungs):
    """Return ``rungs`` sorted ascending by mesh size, so the stop
    test reads a monotone ladder.  The full-mesh point count orders
    them; the mesh tuple breaks any tie deterministically."""
    return sorted(
        rungs,
        key=lambda rung: (_mesh_point_count(rung.mesh),
                          tuple(rung.mesh)))


def _merge_distinct(existing, incoming):
    """Return ``existing`` and ``incoming`` rungs combined, dropping
    any incoming rung whose mesh is already present, and re-sorted
    ascending.  A climb only ever appends meshes above its current
    top, but the distinct merge keeps the ladder clean even if a
    round hands back a mesh already run -- a manufactured zero-delta
    pair would otherwise fool the two-sided stop test (DESIGN 3.11)."""
    seen = {tuple(rung.mesh) for rung in existing}
    merged = list(existing)
    for rung in incoming:
        if tuple(rung.mesh) not in seen:
            merged.append(rung)
            seen.add(tuple(rung.mesh))
    return _sort_by_mesh(merged)


def climb_action(rungs, config):
    """Decide the unit-step climb's next move from its ladder
    (PSEUDOCODE 4e.3, ``climbAction``; DESIGN 3.12.3 / 3.12.5).

    The unit-step climb walks one rung at a time, testing the whole
    accumulated ladder.  Returns a ``ClimbAction``:

    - ``converged`` with the converged ``Rung``, once the per-atom
      energy has gone flat over ``config.flat_needed`` consecutive
      interior rungs -- the same two-sided test the harvest applies
      (``guidance_harvest.pick_converged_climb``);
    - ``ceiling`` when the top rung has reached the per-axis backstop
      without converging, so the search stops non-converged;
    - ``run`` with the next mesh up the climb otherwise
      (``mesh_climb.climb_one_rung``).

    This is also the continuation of a confident ``PARALLEL_GRID``
    whose opening grid did not converge: the grid becomes a unit-step
    climb from its top rung (``climb_next``, PSEUDOCODE 4e.3).

    ``rungs`` are this material's distinct meshes run so far, sorted
    ascending, each a ``Rung(mesh, energy)``."""
    index = guidance_harvest.pick_converged_climb(
        [rung.energy for rung in rungs],
        config.cell_atom_count, config.threshold, config.flat_needed)
    if index is not None:
        return ClimbAction(_ACTION_CONVERGED, rungs[index], None)
    if mesh_climb.at_ceiling(rungs[-1].mesh, config.max_count):
        return ClimbAction(_ACTION_CEILING, None, None)
    next_mesh = mesh_climb.climb_one_rung(
        rungs[-1].mesh, config.classes, config.recip_mag)
    return ClimbAction(_ACTION_RUN, None, next_mesh)


# ==================================================================
#  The bracket-refine climb: the stateful search (PSEUDOCODE 4e.3;
#  DESIGN 3.12.3)
#
#  This is the loop that reads energies and chooses the next mesh --
#  the concern ARCHITECTURE 9.7 keeps in the producer -- so it lives
#  here, beside ``climb_action``, not in the pure-geometry
#  ``mesh_climb``.  The default cold/moderate search does not walk the
#  ladder rung by rung; it BRACKETS with a geometric stride, then
#  FILLS the small bracket it lands in and re-judges the now-
#  consecutive block with the two-sided test.  Its per-material state
#  is threaded by ``converge_by_climb`` (4e.5); the grid and the
#  unit-step climb carry an empty state they never read.
# ==================================================================

# The two phases of the bracket-refine search (DESIGN 3.12.3): stride
#   to bracket the convergence, then fill and re-judge that bracket.
_PHASE_BRACKET = "bracket"
_PHASE_REFINE = "refine"

## One material's bracket-refine search state (PSEUDOCODE 4e.3).
##   Immutable -- each step returns a fresh state via ``_replace`` --
##   so a landing advances exactly one material's search and never
##   touches another's:
##     phase      _PHASE_BRACKET or _PHASE_REFINE
##     stride     the current geometric stride (bracket phase)
##     endpoints  the bracket endpoint meshes computed so far,
##                ascending; endpoints[0] is the seed, endpoints[-1]
##                the newest
##     lo, hi     the interval the refine phase fills, as meshes
##                (None until the refine phase begins)
##     from_cap   True iff this bracket runs up to the ceiling, so an
##                empty refine is a CEILING stop, not a false bracket
BracketRefineState = namedtuple(
    "BracketRefineState",
    ["phase", "stride", "endpoints", "lo", "hi", "from_cap"])


def new_bracket_refine_state(seed_mesh):
    """The opening bracket-refine state for a climb seeded at
    ``seed_mesh`` (PSEUDOCODE 4e.3, ``newBracketRefineState``).

    The search opens in the bracket phase with a stride of one and
    the seed as its only computed endpoint; the first call to
    ``bracket_refine_next`` launches the first stride from it."""
    return BracketRefineState(
        phase=_PHASE_BRACKET, stride=1, endpoints=[seed_mesh],
        lo=None, hi=None, from_cap=False)


def new_search_state(config, opening_meshes):
    """The per-material search state ``converge_by_climb`` threads
    (PSEUDOCODE 4e.4, ``newSearchState``).

    Only the bracket-refine climb is stateful; its state is seeded
    from its single opening rung.  The grid and the unit-step climb
    carry an empty state (``None``) they never read -- a grid seeds no
    state, and the unit-step climb is stateless in its ladder."""
    if config.mode == mesh_climb.BRACKET_REFINE:
        return new_bracket_refine_state(opening_meshes[0])
    return None


def _stride_up(state, stride, config):
    """Launch the next bracket endpoint ``stride`` rungs above the
    current top, recording it on the endpoint list (PSEUDOCODE 4e.3,
    ``strideUp``; DESIGN 3.12.2 / 3.12.3).

    Returns ``(run_action, next_state)``: the ``RUN`` for the new
    endpoint mesh and the state carrying it and the stride just
    taken."""
    next_mesh = mesh_climb.climb_n_rungs(
        state.endpoints[-1], stride, config.classes, config.recip_mag)
    updated = state._replace(
        endpoints=state.endpoints + [next_mesh], stride=stride)
    return ClimbAction(_ACTION_RUN, None, next_mesh), updated


def _enter_refine(rungs, state, lo_mesh, hi_mesh, from_cap, config):
    """Switch a search into the refine phase on ``[lo_mesh, hi_mesh]``
    and launch its first fill mesh (PSEUDOCODE 4e.3, ``enterRefine``;
    DESIGN 3.12.3).

    If the interval is already fully computed -- its endpoints with
    nothing between -- there is nothing to fill, so judge it at once
    by re-entering ``bracket_refine_next`` in the refine phase."""
    updated = state._replace(
        phase=_PHASE_REFINE, lo=lo_mesh, hi=hi_mesh, from_cap=from_cap)
    gap = mesh_climb.next_fill_mesh(
        rungs, lo_mesh, hi_mesh, config.classes, config.recip_mag)
    if gap is None:
        return bracket_refine_next(rungs, updated, config)
    return ClimbAction(_ACTION_RUN, None, gap), updated


def bracket_refine_next(rungs, state, config):
    """Decide the bracket-refine climb's next move (PSEUDOCODE 4e.3,
    ``bracketRefineNext``; DESIGN 3.12.3).

    A two-phase state machine.  In the BRACKET phase it strides by a
    geometric step until a stride reads flat (its endpoints within the
    per-atom threshold), then brackets the LAST non-flat interval and
    switches to refine.  In the REFINE phase it fills the bracket one
    ladder position at a time and applies the authoritative two-sided
    test to the now-consecutive block, returning the smallest mesh
    that passes.  A coincidentally flat stride (an oscillating near-
    metal energy) is caught here -- no rung in a falsely bracketed
    interval passes the two-sided test -- and the search resumes
    striding from the top of the bracket.

    Returns ``(action, next_state)``.  ``rungs`` is the material's
    computed ladder, ascending; ``state`` its bracket-refine state."""
    if state.phase == _PHASE_BRACKET:
        top = state.endpoints[-1]
        if len(state.endpoints) == 1:
            # Only the seed is computed; launch the first stride
            #   (stride 1).  There is no flatness to test yet.
            return _stride_up(state, 1, config)

        # Two or more endpoints: test whether the last stride is flat.
        prev = state.endpoints[-2]
        if guidance_harvest.stride_is_flat(
                mesh_climb.rung_at(rungs, prev),
                mesh_climb.rung_at(rungs, top),
                config.cell_atom_count, config.threshold):
            # First flat stride.  The converged rung lies at or just
            #   above the bottom of the flat stride, prev.  Fill
            #   flat_needed + 1 rungs above prev, so the persistence
            #   test has flat_needed interior candidates each with a
            #   computed neighbour on both sides (the two-sided test
            #   excludes a block's endpoints).  The + 1 matters because
            #   prev itself need not be settled: its own lower
            #   neighbour may still be moving, so the first rung the
            #   test can confirm is often prev + 1, and confirming
            #   flat_needed rungs from there needs a computed neighbour
            #   up through prev + flat_needed + 1.  Filling fewer would
            #   expose too few interior candidates, and a two-
            #   consecutive-flat search could never confirm.  If prev
            #   is the seed (the very first stride was flat), there is
            #   no lower endpoint either, so lo extends one rung below
            #   prev.
            hi_mesh = mesh_climb.climb_n_rungs(
                prev, config.flat_needed + 1, config.classes,
                config.recip_mag)
            if len(state.endpoints) >= 3:
                lo_mesh = state.endpoints[-3]
            else:
                lo_mesh = mesh_climb.descend_one_rung(
                    prev, config.classes, config.recip_mag)
            return _enter_refine(
                rungs, state, lo_mesh, hi_mesh, False, config)

        # Not flat: grow the stride geometrically and step up, unless
        #   the next endpoint would pass the ceiling -- then refine
        #   from the top endpoint up to the ceiling (DESIGN 3.12.3), so
        #   a convergence a stride jumped over just below the cap is
        #   still found.
        next_stride = min(2 * state.stride, config.max_stride)
        candidate = mesh_climb.climb_n_rungs(
            top, next_stride, config.classes, config.recip_mag)
        if mesh_climb.at_ceiling(candidate, config.max_count):
            ceiling = mesh_climb.ceiling_mesh(
                top, config.classes, config.recip_mag,
                config.max_count)
            return _enter_refine(
                rungs, state, top, ceiling, True, config)
        return _stride_up(state, next_stride, config)

    # REFINE: fill [lo, hi] lowest-first, TESTING the consecutive block
    #   after each fill so a convergence low in the bracket stops the
    #   fill before the wide rungs above it are computed (DESIGN
    #   3.12.3).  Test only the CONSECUTIVE run anchored at lo -- not
    #   the whole [lo, hi] range -- because mid-fill the range still
    #   has gaps (a sparse bracket endpoint sitting above an unfilled
    #   rung), and the two-sided test would compare non-neighbours
    #   across such a gap and could read a false convergence.  lo is
    #   always computed by the time the first fill lands
    #   (next_fill_mesh fills it first), so the anchor is safe.  The
    #   fill climbs from the bottom and pick_converged_climb returns
    #   the SMALLEST passing rung, so the converged mesh is exactly the
    #   one a full fill would find -- only the rungs above it go
    #   uncomputed.
    block = mesh_climb.consecutive_block(
        rungs, mesh_climb.rung_at(rungs, state.lo),
        config.classes, config.recip_mag)
    index = guidance_harvest.pick_converged_climb(
        [rung.energy for rung in block],
        config.cell_atom_count, config.threshold, config.flat_needed)
    if index is not None:
        return ClimbAction(_ACTION_CONVERGED, block[index], None), state

    # Not verified yet: fill the next-lowest gap if any remains.
    gap = mesh_climb.next_fill_mesh(
        rungs, state.lo, state.hi, config.classes, config.recip_mag)
    if gap is not None:
        return ClimbAction(_ACTION_RUN, None, gap), state    # keep filling

    # Interval fully filled and still nothing verified.
    if state.from_cap:
        # Still visibly steep even at the cap: a genuine non-converged
        #   ceiling stop, not a false bracket (DESIGN 3.12.3).
        return ClimbAction(_ACTION_CEILING, None, None), state
    # A coincidentally flat stride (an oscillating energy): no rung
    #   verified.  Resume striding from hi (DESIGN 3.12.3).
    return bracket_refine_next(
        rungs, new_bracket_refine_state(state.hi), config)


def climb_next(rungs, state, config):
    """Decide one material's next move, dispatching on its search mode
    (PSEUDOCODE 4e.3, ``climb_next``).

    The bracket-refine climb threads its per-material ``state``; the
    unit-step climb and a grid continuation are stateless in the
    ladder, so they run ``climb_action`` and pass the state through
    untouched.  Returns ``(action, next_state)`` so the caller can
    persist the advanced state."""
    if config.mode == mesh_climb.BRACKET_REFINE:
        return bracket_refine_next(rungs, state, config)
    return climb_action(rungs, config), state


def converge_by_climb(materials, configs, seed_densities,
                      dispatcher, on_non_converged=None):
    """Drive every material through the climb to a verdict
    (PSEUDOCODE 4e.5; DESIGN 3.12.5).

    Serial within a material -- rung N+1's mesh is not known until
    rung N's energy is judged -- but concurrent across materials, and
    NO material waits on another.  The producer launches every
    material's opening rung at once, then collects rungs as they land
    and advances only the one material each landing belongs to.  A
    material that finishes a rung early climbs on immediately instead
    of idling until some slower material's rung completes, and it
    leaves the active set the moment it converges or hits its ceiling,
    so a late, expensive material never holds back the ones already
    done.

    Parameters
    ----------
    materials
        The material identifiers to converge.
    configs
        ``{material: ClimbConfig}`` -- each material's geometry,
        resolved policy, and energy / ceiling knobs.
    seed_densities
        ``{material: float}`` -- the seed density for the opening
        rung (the guidance prediction, or the wide-grid floor for an
        under-trained bootstrap, DESIGN 7.9).  Passed as a density
        rather than a prediction object because that is all the
        seeding needs; the policy was already resolved into
        ``configs``.
    dispatcher
        The injected climb dispatcher (``make_climb_dispatcher``,
        4e.7).  It owns the in-flight set, so this loop tracks only
        its per-material ladders (Principle 12), and exposes two
        calls:

        - ``dispatcher.send(mesh_lists)`` launches one calc per
          ``(material, mesh)`` in ``{material: [mesh, ...]}`` WITHOUT
          waiting (send_off, 13.5).
        - ``dispatcher.next_rung()`` blocks until the next rung lands
          and returns ``(material, result)``, where ``result`` is a
          ``Rung`` or the ``_RUN_FAILED`` marker for a rung that did
          not complete (DESIGN 7.7).
    on_non_converged
        Optional ``on_non_converged(material)`` callback, invoked
        when a material stops non-converged -- at the ceiling, or
        because a rung failed to run -- so the caller can tag the
        prediction mismatch (DESIGN 7.8 step 3d).  Injected and
        defaulting to a no-op so this loop stays free of the
        workspace; the producer wires it when it drives the climb.

    Returns
    -------
    (outcomes, rungs)
        ``outcomes[material]`` is the converged ``Rung`` or the
        ``NON_CONVERGED`` sentinel (a ceiling stop or a run failure);
        ``rungs[material]`` is the full distinct-mesh ladder that
        material climbed, ascending -- the flatness trace the harvest
        re-judges (4e.6)."""

    rungs = {material: [] for material in materials}
    search = {}                 # per-material search state (4e.4)
    outcomes = {}
    active = set(materials)
    in_air = {}                 # rungs still in flight, per material
    opening = set(materials)    # still in the opening (grid) phase

    def retire(material, verdict):
        # Record a material's outcome and drop it from the active
        #   set; a non-converged stop tags the mismatch (7.8 3d).
        outcomes[material] = verdict
        if verdict == NON_CONVERGED and on_non_converged is not None:
            on_non_converged(material)
        active.discard(material)

    def judge(material):
        # Read the next action from a material's ladder and either
        #   retire it or launch its single next rung.  climb_next
        #   threads the per-material search state, so the bracket-
        #   refine phase persists across landings; the grid and the
        #   unit-step climb pass an empty state through untouched.
        action, search[material] = climb_next(
            rungs[material], search[material], configs[material])
        if action.kind == _ACTION_CONVERGED:
            # A converged rung is the outcome itself (not a
            #   mismatch), so it leaves active with no failure tag.
            outcomes[material] = action.rung
            active.discard(material)
        elif action.kind == _ACTION_CEILING:
            retire(material, NON_CONVERGED)              # 7.8 3d
        else:                                            # _ACTION_RUN
            dispatcher.send({material: [action.mesh]})
            in_air[material] += 1

    # Seed every material's opening rung or grid at once: one rung for
    #   a climb, a small grid for the confident mode
    #   (mesh_climb.initial_meshes, 4e.4).  The bracket-refine climb's
    #   search state is seeded from its opening rung; the grid and the
    #   unit-step climb carry an empty state (new_search_state, 4e.4).
    first = {}
    for material in materials:
        config = configs[material]
        policy = mesh_climb.ClimbPolicy(
            config.mode, config.flat_needed, config.grid_width,
            config.start_offset, config.max_stride)
        first[material] = mesh_climb.initial_meshes(
            seed_densities[material], policy, config.classes,
            config.recip_mag, config.recip_cell_volume)
        search[material] = new_search_state(config, first[material])
        in_air[material] = len(first[material])
    dispatcher.send(first)

    # Collect rungs as they land, in landing order, until nothing is
    #   in flight.  Each landing advances exactly the one material it
    #   belongs to; the others are untouched, so no material is paced
    #   by another.
    while any(count > 0 for count in in_air.values()):
        material, result = dispatcher.next_rung()
        in_air[material] -= 1
        if result is not _RUN_FAILED:
            rungs[material] = _merge_distinct(
                rungs[material], [result])

        if material in opening:
            # The confident mode's opening grid is judged as a group,
            #   so wait until the WHOLE grid has resolved, then judge
            #   on whatever landed.  A material whose entire opening
            #   failed has no rung to stand on (run failure, 7.7); a
            #   climb's opening is a single rung, judged at once.
            if in_air[material] > 0:
                continue
            opening.discard(material)
            if not rungs[material]:
                retire(material, NON_CONVERGED)          # run failure
            else:
                judge(material)
        else:
            # A continuation is exactly one rung.  A failed rung means
            #   the climb cannot advance, so stop the material rather
            #   than re-dispatch it forever (7.7); otherwise judge the
            #   extended ladder.
            if result is _RUN_FAILED:
                retire(material, NON_CONVERGED)          # run failure
            else:
                judge(material)

    return outcomes, rungs


def _unit_key(unit):
    """A hashable identity for a climb unit -- its id and calc tag --
    keying the dispatcher's origin map back to (material, mesh)."""
    return (unit.id, tuple(unit.calc))


class _ClimbDispatcher:
    """The send / collect adapter ``converge_by_climb`` drives
    (DESIGN 7.7; PSEUDOCODE 4e.7).

    ONE flight spans the whole climb -- its root is the workspace --
    and its unit list ACCRETES as rungs are decided, so ``flight.toml``
    records every rung asked for (DESIGN 7.7) and any mesh already run
    is a cache hit (DESIGN 6.2.5).  The dispatcher owns the in-flight
    set and a small map from each launched unit to its
    ``(material, mesh)``, so an energy routes straight back without
    re-decoding the calc tag.  Two calls drive it:

    - ``send(mesh_lists)`` builds one explicit-mesh unit per requested
      mesh (:func:`build_mesh_unit`), appends it to the one flight,
      prepares just the new units, and launches them WITHOUT waiting
      (``send_off``).
    - ``next_rung()`` blocks until the next rung lands
      (``collect_next``), reads its completed unit's ``(mesh, energy)``
      into a ``Rung``, and returns ``(material, Rung)`` -- or
      ``(material, _RUN_FAILED)`` for a unit that did not complete
      (DESIGN 7.7).

    The material key doubles as the unit id (materials ARE the
    reference ids the producer already uses).  ``prepare_fn`` /
    ``send_off_fn`` / ``collect_next_fn`` / ``read_fn`` are injected --
    defaulting to the real driver-side prepare, kaleidoscope's real
    send-off and collect, and the real result reader -- so a caller can
    unit-test the climb with the toolchain seam mocked (each live run
    needs a real imago, C74).  ``force`` bypasses the run-reuse cache.
    """

    def __init__(self, structures, options_by_material, workspace,
                 flight, executor, force, prepare_fn, send_off_fn,
                 collect_next_fn, read_fn):
        self._structures = structures
        self._options = options_by_material
        self._workspace = workspace
        self._flight = flight
        self._executor = executor
        self._force = force
        self._prepare_fn = prepare_fn
        self._send_off_fn = send_off_fn
        self._collect_next_fn = collect_next_fn
        self._read_fn = read_fn
        self._outstanding = []      # (unit, future) still in flight
        self._origin = {}           # _unit_key -> (material, mesh)

    def send(self, mesh_lists):
        """Launch one calc per requested mesh WITHOUT waiting.  Each
        new unit is remembered by its origin, appended to the growing
        flight, prepared (only the new units), and handed to
        ``send_off``, which re-serializes the whole flight so
        ``flight.toml`` records every rung."""
        new_units = []
        for material, meshes in mesh_lists.items():
            for mesh in meshes:
                unit = build_mesh_unit(
                    self._structures[material],
                    self._options[material], mesh, material)
                self._origin[_unit_key(unit)] = (material, list(mesh))
                new_units.append(unit)
                self._flight.units.append(unit)
        self._prepare_fn(self._flight, self._workspace, new_units)
        launched = self._send_off_fn(
            self._flight, new_units, self._executor, self._force)
        self._outstanding.extend(launched)

    def next_rung(self):
        """Block until the next rung lands and translate it to
        ``(material, Rung)`` -- or ``(material, _RUN_FAILED)`` for a
        unit that did not complete (the climb stops that material as a
        run failure, 4e.5).  A landed mesh must equal the one
        requested (DESIGN 7.7): a resolved mesh that differs means
        makeinput or imago silently changed it, so fail loudly rather
        than record the wrong rung."""
        unit, entry, remaining = self._collect_next_fn(
            self._flight, self._outstanding)
        self._outstanding = remaining
        material, mesh = self._origin[_unit_key(unit)]
        if entry.status != "done":
            return material, _RUN_FAILED
        result = self._read_fn(self._workspace, unit)
        resolved = result.get("kpoint_mesh")
        if resolved is not None and list(resolved) != list(mesh):
            raise RuntimeError(
                f"requested mesh {list(mesh)} but the run for "
                f"{material!r} resolved {list(resolved)}; an "
                f"explicit scfkp mesh must be honoured exactly "
                f"(DESIGN 7.7)")
        return material, Rung(list(mesh), result["total_energy"])


def make_climb_dispatcher(structures, options_by_material, workspace,
                          *, parsl_config=None, executor=None,
                          prepare_fn=prepare_units,
                          send_off_fn=send_off,
                          collect_next_fn=collect_next,
                          read_fn=_read_unit_result, force=False):
    """Build the climb dispatcher ``converge_by_climb`` drives
    (DESIGN 7.7; PSEUDOCODE 4e.7).

    It closes over each material's ``structure`` and coded ``options``
    (keyed by the material id, which is also the unit id so an energy
    routes straight back), the ``workspace`` the run directories live
    under, the resolved dispatch ``parsl_config`` recorded on the
    climb's one flight (``None`` for the local opt-out, a real Parsl
    ``Config`` for a cluster shape), and the ONE shared ``executor``
    every send runs beneath.  That single executor is what the loen
    pre-flight and every climb rung run under, so the whole run rides
    one warm pool (DESIGN 6.2.11) and its units land in one tree.  See
    :class:`_ClimbDispatcher` for the send / next_rung contract and the
    injected seams."""
    flight = Flight(
        root=workspace, units=[],
        parsl_config=parsl_config,
        sweep=SweepRecord(varied_axes=("kpt-mesh",), fixed_axes={}))
    return _ClimbDispatcher(
        structures, options_by_material, workspace, flight, executor,
        force, prepare_fn, send_off_fn, collect_next_fn, read_fn)


def record_converged(rung, rungs, config):
    """Build the density / mesh / grid harvest inputs for a
    converged climb material (PSEUDOCODE 4e.6; DESIGN 3.12.4).

    The guidance dataspace is keyed on a DENSITY, but the climb
    converges a MESH, so the converged rung is recorded both ways.
    The density a mesh represents is its full-mesh volume density,
    ``product(mesh) / recip_cell_volume`` -- self-consistent with the
    count selection (``mesh_climb.select_axial_counts``), so a future
    prediction of this density reproduces this mesh in this cell --
    and the exact mesh is stored beside it (DESIGN 3.12.4 / 7.2).

    ``rung`` is the converged rung; ``rungs`` its ascending
    distinct-mesh ladder.  The stored flatness trace is the
    CONSECUTIVE block of that ladder around the converged rung -- the
    rungs the two-sided test actually compared
    (``mesh_climb.consecutive_block``, 4e.6).  For the unit-step climb
    and the grid that is the whole ladder; for the bracket-refine
    climb it is the filled bracket, dropping the sparse stride
    endpoints below it (search scaffolding, DESIGN 3.12.3).  Recording
    only the consecutive block is what lets the curator's auto-promote
    rule re-judge on adjacent meshes -- a sparse endpoint left in could
    make the two-sided test read a false early convergence.

    Returns the density / mesh / grid fields only; the producer's
    climb harvest threads them into
    :func:`guidance_harvest.build_entry`, which adds the gap,
    magnetization, sub-model, and provenance (4e.6 / 15.7)."""
    volume = config.recip_cell_volume
    trace = mesh_climb.consecutive_block(
        rungs, rung, config.classes, config.recip_mag)
    return {
        "converged_kpoint_density":
            _mesh_point_count(rung.mesh) / volume,
        "converged_mesh": list(rung.mesh),
        # Ascending because `trace` is and the point count rises with
        #   each rung; raw total-cell energies (Option B), which the
        #   consumer normalizes per atom (DESIGN 7.8).
        "grid_values": [_mesh_point_count(one.mesh) / volume
                        for one in trace],
        "grid_energies": [one.energy for one in trace],
    }


def _lattice_rows(one_indexed_lattice):
    """Extract a plain 3x3 (rows a/b/c, columns x/y/z) from a
    ``StructureControl`` 1-indexed 4x4 lattice.

    ``StructureControl`` stores lattices Perl-port style: index 0 is
    an unused ``None`` sentinel, so the three lattice vectors live at
    rows 1..3 and their Cartesian components at columns 1..3.
    ``mesh_climb.axis_classes_for_cell`` wants the lattice vectors as
    ordinary 0-indexed rows, which this returns."""

    return [[one_indexed_lattice[axis][component]
             for component in range(1, 4)]
            for axis in range(1, 4)]


def build_climb_config(ref, structure, confidence, under_trained,
                       thresholds, max_count):
    """Assemble one reference solid's :data:`ClimbConfig` (PSEUDOCODE
    11.4; DESIGN 3.12 / 5.7) from three sources: the loaded cell's
    reciprocal geometry (the rung mechanics search there), the
    confidence-derived climb policy, and the run's energy / ceiling
    knobs.

    Called once per solid in the build phase, BEFORE any run -- the
    climb needs a cell's axis classes to seed and step, so they are
    recomputed here from the cell's own space-group operations
    (``cell.point_ops`` -> the shared reader, PSEUDOCODE 4b.4), never
    read back from a run.

    Parameters
    ----------
    ref
        The ``ReferenceSolid``, for its resolved per-atom k-point
        flatness threshold (the climb's stop tolerance).
    structure
        The solid's local ``.skl`` path.
    confidence, under_trained
        The predictor's confidence and bootstrap flag, which shape
        the search (parallel grid vs. serial climb; DESIGN 3.12.6).
    thresholds, max_count
        The run-wide climb policy resolved once from the manifest
        (:func:`mesh_climb.climb_policy_from_manifest`, 4e.4): the
        confidence-to-mode ``thresholds`` bundle and the per-axis
        ceiling ``max_count``.
    """

    cell = guidance_harvest.load_structure(structure)

    # Recompute the reciprocal lattice from the FINAL loaded cell
    #   (primitive for a prim reduction), so the reciprocal-axis
    #   magnitudes and cell volume match the mesh the climb builds
    #   (DESIGN 3.12.4).  make_inv_or_recip_lattice fills
    #   recip_lattice and recip_cell_volume (the 2*pi convention);
    #   it is idempotent, so calling it here is safe whatever the
    #   loader already computed.
    #
    # Unlike real_lattice (whose ROWS are the cell vectors), the
    #   reciprocal is stored as an inverse, so its reciprocal
    #   vectors b_i are its COLUMNS -- component index first, vector
    #   index second.  The magnitude of reciprocal axis i is thus
    #   the norm DOWN column i.  (Reading it row-wise transposes the
    #   vectors and silently mis-scales every non-cubic cell -- a
    #   cubic cell hides it because its row and column norms match.)
    cell.make_inv_or_recip_lattice(make_recip=True)
    recip_lattice = cell.recip_lattice
    recip_mag = [
        (recip_lattice[1][axis] ** 2 + recip_lattice[2][axis] ** 2
         + recip_lattice[3][axis] ** 2) ** 0.5
        for axis in range(1, 4)]

    # Axis classes (4c.7 / DESIGN 2.7): which reciprocal axes must
    #   share a k-point count.  The rotations come from the cell's
    #   own space group via the shared reader, so the classes the
    #   climb seeds from are derived from the operations imago will
    #   run under.  cell_mode is "full" when the loaded cell IS the
    #   conventional cell (do_full_cell), else "prim".
    cell_mode = "full" if cell.do_full_cell else "prim"
    classes = mesh_climb.axis_classes_for_cell(
        cell.point_ops(),
        _lattice_rows(cell.real_lattice),
        _lattice_rows(cell.full_cell_real_lattice),
        cell_mode)

    # Confidence -> the shape of the search (4e.4): dispatch mode,
    #   flatness persistence, grid width, climb start offset.  The
    #   knob thresholds were resolved once for the run from the
    #   manifest (climb_policy_from_manifest).
    policy = mesh_climb.resolve_climb_policy(
        confidence, under_trained, thresholds)

    return ClimbConfig(
        classes=classes,
        recip_mag=recip_mag,
        recip_cell_volume=cell.recip_cell_volume,
        mode=policy.mode,
        flat_needed=policy.flat_needed,
        grid_width=policy.grid_width,
        start_offset=policy.start_offset,
        max_stride=policy.max_stride,
        cell_atom_count=cell.num_atoms,
        threshold=ref.kpoint_convergence_threshold,
        max_count=max_count)


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
                             prepare_fn=prepare_units,
                             extract_fn=extract_potential,
                             identity_fn=read_site_identity_map,
                             fingerprint_fn=harvest_fingerprints
                             ) -> list[dict[str, Any]]:
    """The three-phase producer (DESIGN 5.7; PSEUDOCODE 11.4):
    *build* each solid's ClimbConfig from a predicted seed density
    (no grid), *converge* every solid through the adaptive mesh
    climb (predict in density, search in mesh), then *harvest* each
    converged solid's potential and contribute the same converged
    rung back to the guidance dataspace.  Returns the per-run log
    (also written to disk).

    ``dispatch_shape`` selects where the flight runs (DESIGN 6.2.11):
    ``local`` (the default here, so tests and in-process callers stay
    serial without any cluster settings file) runs every unit in
    process; ``slurm-pooled`` and ``slurm-per-job`` build a Parsl
    ``Config`` from the per-site settings file and the per-run
    ``partition`` / ``nodes`` / ``walltime`` / ``profile`` choices.
    ``save_config`` records the resolved cluster choices beside the
    run.  ``force`` bypasses the run-reuse cache (DESIGN 6.2.5).

    ``dispatch_fn``, ``prepare_fn``, ``extract_fn``, ``identity_fn``,
    and ``fingerprint_fn`` are injected (defaulting to the real
    kaleidoscope dispatch, the real driver-side prepare pass, the
    real scfV reader, the real ``datSkl.map`` reader, and the real
    fingerprint harvest) so the orchestration can be unit-tested
    with the toolchain seam mocked:
    end-to-end dispatch, the per-site ``scfV`` read, the site-identity
    read, and the fingerprint harvest (which needs the run's expanded
    structure and loen descriptor) all need a live Imago run (C74)."""

    # Pass the matcher registry so the loader enforces rule 9: every
    #   declared fingerprint method must be a registered matcher (C54).
    manifest = load_manifest_v2(
        manifest_path, known_methods=set(MATCHERS))
    apply_manifest_defaults(manifest)
    manifest_dir = os.path.dirname(manifest.manifest_path)
    guidance_root = os.path.join(data_root, "historicalGuidanceDB")
    dataspace = guidance_db.load(guidance_root)
    imago_commit = _git_sha()
    timestamp = _now_iso8601_utc()
    workspace = curation_workspace_root(pdb_root)

    # Resolve the climb policy ONCE for the run (PSEUDOCODE 4e.4 /
    #   DESIGN 3.12.6): the optional [harvest.kpoint_climb] knobs
    #   merged over the provisional defaults into the confidence
    #   `thresholds` bundle and the per-axis `max_count` ceiling every
    #   solid's ClimbConfig reads.  An empty sub-table yields the
    #   built-in policy; a mistyped knob already failed loudly at load.
    thresholds, max_count = mesh_climb.climb_policy_from_manifest(
        manifest.harvest.get("kpoint_climb", {}))

    # ----- Phase 1: build.  Refresh the isolated baselines, then per
    # solid: materialize the structure, PREDICT its seed k-point
    # density (no grid), and build its ClimbConfig.  The convergence
    # units are NOT built here -- the climb builds each round's meshes
    # as it runs (Phase 2).  Only the geometry-only fingerprint units
    # are built now, collected into one pre-flight batch (they are
    # mesh-independent, so they belong to no climb round; DESIGN 5.7).
    elements = None if single_element is None else [single_element]
    databases = refresh_isolated_entries(
        pdb_root, manifest, imago_commit, timestamp, elements)

    struct_of: dict[str, str] = {}       # reference_id -> .skl path
    options_of: dict[str, Any] = {}      # reference_id -> coded options
    configs: dict[str, ClimbConfig] = {}      # reference_id -> config
    seed_densities: dict[str, float] = {}     # reference_id -> density
    predictions: dict[str, Any] = {}     # reference_id -> record dict
    loen_units: list = []                # every solid's fingerprint units
    for ref in manifest.reference_solids:
        struct = materialize_structure(ref, manifest_dir, pdb_root)
        struct_of[ref.reference_id] = struct
        options = make_producer_options(ref, imago_commit)
        options_of[ref.reference_id] = options
        # The predictor and the PredictionRecord speak the human
        # physics names, not the codes, so the sub-model travels in
        # its OWN dict -- never mixed into the tool-facing options
        # (DESIGN 6.2.8 / 6.2.10), which would duplicate the basis and
        # make makeinput reject "functional" / "kpoint_integration".
        submodel = {
            "basis": ref.basis,
            "functional": ref.functional,
            "kpoint_integration": ref.kpoint_integration,
        }
        # PREDICTION ONLY (PSEUDOCODE 4e.7): the climb seeds from the
        #   density and picks its mode / persistence from the
        #   confidence (DESIGN 3.12.4 / 3.12.6).  A pinned
        #   kpoint_spec.density is the curator override (predictor
        #   bypassed); otherwise None and predict runs, flagging
        #   under-trained when the dataspace has no useful prior (7.9).
        density, confidence, under_trained, record = \
            predict_kpoint_density(
                struct, dataspace, ref.system_type, submodel,
                center=ref.kpoint_spec.get("density"))
        seed_densities[ref.reference_id] = density
        # Everything the climb needs for this solid, gathered once: the
        #   reciprocal geometry the rung mechanics read, the confidence-
        #   derived policy, and the energy / ceiling knobs.
        configs[ref.reference_id] = build_climb_config(
            ref, struct, confidence, under_trained,
            thresholds, max_count)
        # Store the record as a plain dict (metadata must be TOML-
        #   serializable), and stamp the resolved per-atom k-point
        #   flatness tolerance onto it: the guidance harvest reads
        #   both, and the tolerance is a manifest/resolved fact absent
        #   from any run's result.toml (DESIGN 7.8 / 5.7).
        predictions[ref.reference_id] = asdict(record)
        predictions[ref.reference_id][
            "kpoint_convergence_threshold"] = (
                ref.kpoint_convergence_threshold)
        # Geometry-only fingerprint units: one structure-only
        #   `-loen -scf no` unit per Fortran-side declaration.  The
        #   bispectrum fingerprint depends on geometry alone, so these
        #   need not wait for a converged mesh; they dispatch in the
        #   pre-flight below and their run dirs persist for the harvest.
        loen_units.extend(build_loen_units(
            ref, struct, options, manifest.characterization))

    # ----- Phase 1b: loen pre-flight (DESIGN 5.7 / 11.4).  Resolve
    # the dispatch Config ONCE for the whole run (DESIGN 6.2.11):
    # local -> None (the driver runs in process); a cluster shape -> a
    # real Parsl Config.  Both this pre-flight and every climb round
    # dispatch under it; `force` bypasses the run-reuse cache.  The
    # fingerprint units are mesh-independent, so they dispatch once
    # here as a small flat batch (makeinput runs driver-side: the
    # cache keys on structure.dat, DESIGN 6.2.5); their run dirs
    # persist for the Phase 3 fingerprint harvest.  The convergence
    # units are prepared per round inside the climb instead.
    parsl_config, dispatch_choices = resolve_dispatch(
        dispatch_shape, partition, nodes, walltime, profile)
    if save_config and dispatch_choices is not None:
        write_resolved_dispatch(workspace, dispatch_choices, profile)

    # One executor for the whole run's dispatch (DESIGN 6.2.11's
    #   pooled shape): the loen pre-flight and every climb round share
    #   one warm pool.  Build it ONCE (make_executor) and close it
    #   ONCE, in a finally so a mid-climb error still releases the
    #   SLURM allocation.  Phase 3 harvest reads run dirs only and
    #   needs no dispatch, so the pool is freed as convergence ends.
    executor = make_executor(parsl_config)
    try:
        loen_flight = Flight(
            root=workspace, units=loen_units,
            parsl_config=parsl_config,
            sweep=SweepRecord(varied_axes=(), fixed_axes={}))
        if loen_units:
            prepare_fn(loen_flight, workspace)
            dispatch_fn(loen_flight, executor=executor, force=force)

        # ----- Phase 2: converge.  Drive every solid through the
        # adaptive mesh climb (converge_by_climb): serial within a
        # solid, concurrent across, each solid climbing the instant
        # its own rung lands so no solid waits on another.  The climb
        # dispatcher closes over each solid's structure, coded options,
        # and the shared executor, and owns one flight whose unit list
        # accretes as rungs are decided; its send launches new rungs
        # and its next_rung collects whichever lands first, until each
        # solid flattens (its converged rung) or hits the max_count
        # ceiling (NON_CONVERGED).  A mesh re-run later is a cache hit
        # (DESIGN 6.2.5).  on_non_converged tags the solid's workspace
        # with a prediction mismatch (DESIGN 7.8 step 3d); it is
        # injected so the climb loop stays free of the workspace.
        dispatcher = make_climb_dispatcher(
            struct_of, options_of, workspace,
            parsl_config=parsl_config, executor=executor,
            prepare_fn=prepare_fn, force=force)
        materials = [ref.reference_id
                     for ref in manifest.reference_solids]
        outcomes, rungs = converge_by_climb(
            materials, configs, seed_densities, dispatcher,
            on_non_converged=lambda material:
                guidance_harvest.tag_prediction_mismatch(
                    workspace, material))
    finally:
        executor.close()

    # ----- Phase 3: harvest.  Per solid, locate the converged rung's
    # run, extract the potential at each named site, record the run,
    # and contribute the same converged rung back to the guidance
    # dataspace in memory.  A NON_CONVERGED solid (a ceiling stop or a
    # run failure) is logged and skipped -- no potential, no entry.
    per_run_log: list[dict[str, Any]] = []
    for ref in manifest.reference_solids:
        struct = struct_of[ref.reference_id]
        outcome = outcomes[ref.reference_id]
        if outcome is NON_CONVERGED:
            per_run_log.append(make_nonconverged_log_entry(ref))
            continue

        # record_converged turns this solid's converged rung and its
        #   ladder into the guidance density, exact mesh, and flatness
        #   trace in memory (DESIGN 3.12.4).  Both the run log and the
        #   guidance entry read it, so it is computed once, here.
        harvest_inputs = record_converged(
            outcome, rungs[ref.reference_id],
            configs[ref.reference_id])

        # The converged rung names the mesh whose run carries the
        #   converged potential (DESIGN 3.12.4 / Q4).  Rebuild the SAME
        #   mesh unit the climb dispatched -- a cache hit, so it costs
        #   nothing -- to point at that run dir, and read its
        #   result.toml for the iteration count and measured character.
        converged = build_mesh_unit(
            struct, options_of[ref.reference_id], outcome.mesh,
            ref.reference_id)
        converged_result = _read_unit_result(workspace, converged)
        per_run_log.append(
            make_run_log_entry(ref, harvest_inputs, converged_result))
        scf_iterations = converged_result.get("scf_iterations")

        # Discover one representative per distinct environment in the
        #   converged run (DESIGN 5.2.3 / 5.7).  Each carries its
        #   resolved label, default flag, description, and fingerprint
        #   overrides, with any site-pinned customization already
        #   layered on; the harvest only extracts the potential and
        #   inserts the entry.
        for env in discover_environments(
                converged_result, ref, identity_fn=identity_fn):
            elem_key = env.element.lower()
            if elem_key not in databases:
                continue        # filtered out by --element
            coefficients, alphas = extract_fn(
                converged_result, env.atom_site)

            new_entry = ipdb.PotentialEntry(
                label=env.label,
                default=env.default,
                description=env.description,
                num_gaussians=len(coefficients),
                alpha_min=min(alphas),
                alpha_max=max(alphas),
                # The entry keeps this representative's harvested
                # potential verbatim; insert_or_skip keeps the first
                # representative on a duplicate (DESIGN 5.2.3).
                coefficients=coefficients,
                alphas=alphas,
                provenance=make_imago_provenance(
                    imago_commit, timestamp, ref, env.atom_site,
                    scf_iterations),
                # Every environment harvests the database-wide
                #   [characterization] recipe (preferred) plus any
                #   per-entry override (non-preferred).  Fortran-side
                #   matchers read the pre-flight loen run.
                fingerprints=fingerprint_fn(
                    loen_flight, ref, env.atom_site, env.overrides,
                    converged_result, manifest.characterization))
            # Insert-or-skip (DESIGN 5.2.3): replace a same-label
            #   entry (a re-harvested solid or a curator override),
            #   skip a bispectrum duplicate (the stored representative
            #   stands), or append a genuinely novel environment.
            insert_or_skip(databases[elem_key], new_entry)

        # In-memory guidance contribution (DESIGN 5.7 / 11.4).  The
        #   climb already holds this solid's converged rung and its
        #   ascending ladder, so the entry is built in place -- no
        #   re-read of the workspace.  record_converged turns the rung
        #   into the entry's density and flatness ladder; those chosen
        #   facts, with the converged run's result.toml (gap /
        #   magnetization / SCF threshold / exact mesh / commit), feed
        #   the SHARED build_entry -- the same builder the standalone
        #   density harvest uses -- and save_entry stages it.  A
        #   converged climb always carries >= 3 distinct rungs (the
        #   two-sided stop test, DESIGN 3.12.3), so every converged
        #   solid contributes an entry.
        entry = guidance_harvest.build_entry(
            workspace, struct, predictions[ref.reference_id],
            dataspace, guidance_harvest.load_structure(struct),
            ref.kpoint_convergence_threshold,
            harvest_inputs["grid_values"],
            harvest_inputs["grid_energies"],
            harvest_inputs["converged_kpoint_density"],
            converged_result)
        guidance_harvest.save_entry(entry, guidance_root)

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


def submit_orchestrator_batch(argv: list[str], args,
                              data_root: str) -> str:
    """Materialize-then-submit (DESIGN 6.2.11): build the
    orchestrator's sbatch script and submit it, returning the SLURM
    job id.

    Called after the login-node materialize pre-flight has fetched
    every structure into the shared cache.  The batch job re-invokes
    this producer with the SAME arguments minus ``--submit``, so it
    runs the full build under ``--dispatch`` with the structures
    already materialized (no network on the compute node).  The
    driver's own resources come from the site's ``orchestrator``
    block (ARCHITECTURE 9.4), sized to the dispatch shape, with any
    ``--orchestrator-*`` flag overriding that block for this run.

    ``argv`` is the flag vector :func:`main` parsed, with no program
    name in it.  It is threaded in rather than read from the process,
    because a library caller may drive ``main(argv)`` with a vector
    that has nothing to do with the process's own arguments -- and
    reading the process would then submit a batch job that re-runs
    whatever launched us.  The script names itself by path for the
    same reason.
    """
    # The driver's own job is sized from the same overlaid site its
    #   units are, so the queue rides into the loader here too: a
    #   debug queue that caps walltime caps the driver's job as well
    #   (DESIGN 6.2.11).
    site = load_site_config(args.profile, args.partition)
    choices = resolve_choices(site, args)
    orchestrator = resolve_orchestrator(site, args)
    # Re-run this producer in the batch job, dropping --submit so the
    #   batch invocation runs the build itself (not another submit).
    #   The --orchestrator-* flags may ride along untouched: they are
    #   read only when submitting, and the inner run does not submit.
    inner = [item for item in argv if item != "--submit"]
    command = " ".join(
        shlex.quote(item) for item in
        [sys.executable, os.path.abspath(__file__), *inner])
    script_text = build_orchestrator_sbatch(site, choices, command,
                                            orchestrator)
    script_path = os.path.join(data_root, "orchestrator.sbatch")
    with open(script_path, "w") as handle:
        handle.write(script_text)
    completed = subprocess.run(
        ["sbatch", script_path],
        capture_output=True, text=True, check=True)
    # sbatch prints "Submitted batch job <id>".
    return completed.stdout.strip().split()[-1]


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
    optionally redirecting the output with ``--materialize-dir``.
    ``--submit`` materializes on the login node and then submits the
    producer as its own batch job (DESIGN 6.2.11)."""

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
        "--dispatch", default=None,
        choices=["local", "slurm-pooled", "slurm-per-job"],
        help="where to run the flight: 'local' runs every unit in "
             "process and needs no cluster settings file; "
             "'slurm-pooled' streams units through one shared "
             "allocation; 'slurm-per-job' submits one scheduler job "
             "per unit (default: the cluster settings file's "
             "default_topology value)")
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
    parser.add_argument(
        "--submit", action="store_true",
        help="materialize structures on the login node, then submit "
             "the producer as its OWN batch job that runs the build "
             "under --dispatch (DESIGN 6.2.11); the driver job's "
             "resources come from the clusterrc 'orchestrator' block "
             "(default: run the build in this process)")
    # The driver's own resources, overriding the site's orchestrator
    #   block key by key for this run.  Only read when submitting.
    parser.add_argument(
        "--orchestrator-cores", type=int, metavar="N",
        help="cores to request for the --submit driver job "
             "(default: the clusterrc 'orchestrator' block's cores)")
    parser.add_argument(
        "--orchestrator-memory", metavar="SIZE",
        help="memory to request for the --submit driver job, as a "
             "scheduler size such as 8G (default: the clusterrc "
             "'orchestrator' block's memory)")
    parser.add_argument(
        "--orchestrator-walltime", metavar="HH:MM:SS",
        help="time limit for the --submit driver job; it must "
             "outlast the flight it supervises (default: the "
             "clusterrc 'orchestrator' block's walltime, else the "
             "run's --walltime)")
    # The flag vector this run was given, program name excluded.  It is
    #   what --submit re-invokes inside the batch job, so it must be the
    #   vector we actually parsed, not the process's arguments (a
    #   library caller may drive main(argv) from a different process).
    flags = list(argv) if argv is not None else sys.argv[1:]
    args = parser.parse_args(flags)

    if not args.pdb_root:
        parser.error("--pdb-root not given and $IMAGO_DATA is unset")
    if args.materialize_dir and not args.materialize_only:
        parser.error("--materialize-dir applies only with "
                     "--materialize-only")
    if args.submit and args.materialize_only:
        parser.error("--submit and --materialize-only are mutually "
                     "exclusive modes")
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

    # Materialize-then-submit (DESIGN 6.2.11): the orchestrator runs
    #   as its own batch job.  Fetch every structure on the login
    #   node first (the only network step -- compute nodes may lack
    #   internet), then submit a batch job that re-runs this producer
    #   under --dispatch, reading the materialized cache with no
    #   further network.
    if args.submit:
        report = materialize_only(args.manifest, args.pdb_root)
        _print_materialize_report(report)
        if not all(row["ok"] for row in report):
            print("producer: materialize failed; not submitting the "
                  "orchestrator batch job")
            return 1
        job_id = submit_orchestrator_batch(flags, args, data_root)
        print("producer: submitted orchestrator batch job "
              + job_id)
        return 0

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
