# Imago

Imago is an all-electron electronic structure code implementing the
OLCAO (Orthogonalized Linear Combination of Atomic Orbitals) method.
It applies to a broad range of material systems: crystals, amorphous
solids, nanoparticles, molecules, interfaces, and grain boundaries.

Key characteristics:

- **All-electron** -- no pseudopotentials; core electrons are
  treated explicitly.
- **Periodic boundary conditions** -- used throughout, even for
  non-periodic systems.
- **LCAO basis** -- wavefunctions expressed as linear combinations
  of atomic orbitals.
- **Orthogonalization** -- valence orbitals are orthogonalized to
  the core, shrinking the secular (eigenvalue) equation.

Imago is the successor to the long-running OLCAO codebase. The
method retains its established name; only the code has been
rebranded.

## Building

See [BUILD.md](BUILD.md) for the production build and the opt-in
debugging builds.

## Documentation

Project documentation is in [`docs/`](docs/). The design document
chain -- vision, architecture, design, and pseudocode -- lives in
[`dev/`](dev/) and is the authoritative description of how and why
the code works the way it does.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md) for licensing, credit, and
the design chain discipline that governs changes to `src/`.

## License

Imago is licensed under the Educational Community License, Version
2.0 (ECL-2.0). See [LICENSE](LICENSE) for the full text and
[NOTICE](NOTICE) for attribution notices that must accompany
redistribution.

Copyright (c) 2026 Paul Rulis.

## Citation

If you use Imago in published work, please cite it. Machine-readable
metadata is in [CITATION.cff](CITATION.cff); GitHub renders it as a
"Cite this repository" panel, and `cffconvert` will turn it into
BibTeX or other formats.

A DOI will be minted at the first tagged release and added to
`CITATION.cff`. Until then, cite the repository directly.

---

*Computational Physics Group, University of Missouri--Kansas City.*
