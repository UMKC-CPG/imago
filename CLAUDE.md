# Imago Project

> **Session setup:** Run `/color blue` and `/rename Imago` at the
> start of each session.

## Document Hierarchy

This project follows a five-level document chain. All design
documents live in `dev/`. Read them in order when starting work:

1. `dev/VISION.md` -- goals, principles, non-negotiables
2. `dev/ARCHITECTURE.md` -- layout, modules, dependencies
3. `dev/DESIGN.md` -- algorithms, data structures, math
4. `dev/PSEUDOCODE.md` -- language-agnostic algorithm specs
5. Source code in `src/`

`dev/TODO.md` tracks tasks organized by level. Use `/focus` to
start a session and `/refine` to check consistency across the
chain.

## Level Awareness

During development conversations, the programmer may shift
between levels of the design chain without explicitly noticing.
For example, a discussion about a code fix may drift into
questioning an algorithm's design, or a design discussion may
surface a conflict with a core principle.

When you notice the conversation has moved to a different level
than where it started, say so briefly. For example: "This sounds
like it's becoming an ARCHITECTURE question -- should we capture
it there before continuing with the code?" The goal is awareness,
not interruption. Let the programmer decide whether to switch
context, propagate the change to the appropriate document, or
stay focused and defer.

Do not enforce rigid boundaries. The levels exist to organize
thinking, not to prevent it. A developer who is on a productive
train of thought should not be stopped -- but when the thought
resolves, help them recognize which documents it touches so
nothing is left inconsistent.

## Coding Style

### Line Length

Lines MUST NOT exceed 80 characters. This is a hard limit.

Short statements whose content is naturally brief (`implicit
none`, `endif`, `integer :: i`) are fine at their natural length
-- there is nothing to fill them with. The rule applies when
content is available: long expository comments, argument lists,
complex expressions, etc.

**Common failure mode -- do not do this:**
```fortran
   ! Guard each deallocation because
   !   this subroutine is called from
   !   multiple program paths.
```
Those three ~35-character lines have plenty of content to fill
longer lines. Write them as:
```fortran
   ! Guard each deallocation because this subroutine is called
   !   from multiple program paths.
```
The idea: let each line run toward 80 before wrapping. Do not
break at natural phrase boundaries when the line is only at
40-50 characters. The result will be fewer, fuller lines -- not
lines padded with filler.

### Documentation and Naming

CRITICAL: All program code must include rich, expressive
documentation so that students reading the source can easily
follow what is happening. Every function, class, and non-trivial
block should carry clear explanatory comments or docstrings that
describe purpose, inputs, outputs, and any relevant physics or
math.

Variable names must be readable and self-documenting. Avoid
cryptic one- or two-letter abbreviations. Prefer concise but
meaningful names that a reader can understand without
cross-referencing a legend. Slightly-too-long names are far
better than opaque short ones.

Good naming examples:
- `elec_mom` instead of `em` or `electron_momentum`
- `nuc_pot` instead of `np` or `nuclear_potential`
- `grid_spacing` instead of `gs` or
  `the_spacing_between_grid_points`

The goal is a middle ground: short enough to keep expressions
tidy, long enough that any student can read the code cold and
follow the logic without guessing what a variable holds.

### Structured Comment Blocks

Some comments contain *structured* content whose visual layout
is itself meaningful: equations with aligned `=` signs, ASCII
tables, multi-line derivations, sub-lists keyed by hand-aligned
labels, or commented-out code that you may re-enable later.
These blocks must be marked so the reflow tools
(`rewrap_prose.py`) leave them alone.  Without the marker the
tool will treat them as flowing prose and mash the lines into a
paragraph, destroying the layout.

The marker convention is per-language but is conceptually
identical: the comment opener gets *doubled* to signal
"structured content, do not reflow."

**Python** -- prefix the comment lines with `##` (double hash)
instead of `#`.  Any line that begins with `##` is treated as
protected and is never reflowed.  Use this for:
  - Commented-out code blocks
  - Structured equation blocks (multi-line derivations, aligned
    `=`)
  - ASCII tables and aligned sub-lists
  - Any block where the visual layout carries meaning

Example:
```python
# This prose comment may be reflowed by rewrap_prose.py.

## K_theta = 0.15 * sqrt(K_arm1 * K_arm2) * scale
##         = 0.15 * sqrt(400 * 900) * 1.0
##         = 90.0
```

**Fortran** -- prefix the comment lines with `!!` (double
exclamation, with content after) for structured prose, OR use
`!` immediately followed by content (no space) for
commented-out code.  Both forms are recognised by the reflow
engine as "do not touch."

Example:
```fortran
! This prose comment may be reflowed by rewrap_prose.py.

!! K_theta = 0.15 * sqrt(K_arm1 * K_arm2) * scale
!!         = 0.15 * sqrt(400 * 900) * 1.0
!!         = 90.0

!do i = 1, n           ! commented-out code: no space after !
!  call compute(i)
!end do
```

When writing a comment, ask: "if a reflow tool joined these
lines into one paragraph, would I lose meaningful structure?"
If yes, use the doubled form (`##` for Python, `!!` for
Fortran).  Plain `#`/`! ` (single, with space) is the right
choice only for genuine free-flowing prose paragraphs that have
no internal structure.

## What This Is

Imago is an electronic structure code that implements the OLCAO
(Orthogonalized Linear Combination of Atomic Orbitals) method.
It is applicable to a broad range of material systems: crystals,
amorphous solids, nanoparticles, molecules, interfaces, grain
boundaries, and more. Key characteristics:
- **All-electron**: no pseudopotentials; core electrons are
  treated explicitly
- **Periodic boundary conditions**: used throughout, even for
  non-periodic systems
- **LCAO basis**: wavefunctions expressed as linear combinations
  of atomic orbitals
- **Orthogonalization**: valence orbitals are orthogonalized to
  the core, reducing the size of the secular (eigenvalue)
  equation

Imago is the successor to the long-running OLCAO codebase.  The
method retains its established name; only the code has been
rebranded.

## Working Convention

Use `/focus <name>` to load context for a specific subprogram or
script before working on it. This avoids reading the whole
codebase unnecessarily.

## Repository Layout

```
src/                  Source code
  imago/              Primary electronic-structure engine
                      (HDF5-based)
  atomSCF/            Atomic SCF (generates basis/potential)
  makeKPoints/        k-point mesh generation
  gaussFit/           Gaussian fitting of charge density
  contract/           Basis set contraction
  applySpaceGroup/    Space group symmetry operations
  scripts/            Utility scripts (Python, bash)
  tests/              Test suite
  data/               Source databases (unpacked to share/ at
                      install time)
  kinds.f90           Shared: precision kinds
  constants.f90       Shared: physical/mathematical constants
  elementData.f90     Shared: periodic table data
  readData.f90        Shared: input file parsing
  writeData.f90       Shared: output file writing
  radialGrid.f90      Shared: radial grid definitions
bin/                  Installed executables and scripts
                      (generated)
build/                CMake out-of-source build directory
                      (generated)
share/                Basis set and potential databases
                      (generated; built by unpackIMAGODB)
dev/                  Design document chain
.imago/               Resource-control files for the user's
                      installed scripts
rappture/             nanoHUB / Rappture integration
middleware/           nanoHUB invoke wrapper
skl/                  Example skeleton input files
```

## Language and Build

- Fortran 90 for all programs; Python and bash for scripts
- CMake build system; out-of-source build in `build/`
- Compiler set via `$FC` env var (must match the HDF5-compiled
  Fortran compiler)
- Install prefix set via `$IMAGO_DIR` env var
- HDF5 is used for primary I/O in `imago`
- Supported compilers: gfortran, ifort

## Key Conventions

- All Fortran files use `implicit none` and the shared
  `kinds.f90` precision types
- Shared modules in `src/` are compiled once and used across
  subprograms
- Scripts are standalone tools

## Documentation Policy

This is an academic codebase used by students who frequently
need to read and understand the source code.  When refactoring
or writing any code (scripts, Fortran programs, Python
modules):
- **Preserve all existing documentation.**  Every comment
  block, usage note, option explanation, and conceptual
  description must be carried over when restructuring code.
  Do not summarize or abbreviate.
- **Use appropriate formats.**  In Python, use module-level
  docstrings, class/method docstrings, argparse help text, and
  inline comments.  In Fortran, use header comment blocks and
  inline comments.
- **Explain the "why", not just the "what".**  Students benefit
  from knowing the physics/chemistry motivation, not just the
  code mechanics.
