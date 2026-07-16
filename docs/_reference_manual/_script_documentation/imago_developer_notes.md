---
# Developer-oriented notes for the imago.py calculation options.
#
# `published: false` keeps this file in the repository but tells
# Jekyll NOT to build it into the site. It is a holding place for
# engine-internal details (Fortran source files, subroutine names,
# control flags, scratch-file unit numbers) that were pulled out of
# the user-facing imago.md page. If a developer-oriented section of
# the documentation is created later, this content can seed it.
published: false
script: imago.py
description: Engine-internal implementation notes for imago.py calculation options (not user-facing).
---

# imago.py -- Developer Notes

These notes record how each calculation option is implemented in the
Fortran engine and the `imago.py` driver. They are intended for
people developing imago, not for people running it, and are kept out
of the rendered site.

## Optical Properties (-optc / -scfoptc)

Two-stage pipeline:

1. **Optical conductivity and epsilon2.** The Fortran engine
   (`optc.F90` / `optcPrint.F90`) evaluates the interband optical
   conductivity and epsilon2 from the momentum matrix elements
   between occupied and unoccupied states, broadening each
   transition by a Gaussian (width `sigma`) weighted by k-point
   weight and summing over the Brillouin zone. Conductivity is
   written to scratch unit `fort.40`, epsilon2 to `fort.50`.
2. **Kramers-Kronig conversion.** The `imagoKKc` helper reads
   `fort.50` and produces epsilon1 (`fort.110`), the energy loss
   function (`fort.120`), refractive index (`fort.130`), extinction
   coefficient (`fort.140`), the alternate imaginary-integrand
   epsilon1 (`fort.150`), reflectivity (`fort.160`), absorption
   coefficient (`fort.170`), and a combined epsilon1/epsilon2/ELF
   file (`fort.100`).

`imago.py:_manage_optc_output` renames these scratch files into the
job directory. Decomposed files are written in Fortran `5e15.7`
format (five columns); the combined file is `4e15.7` (four columns).

**Partial optical properties (POPTC)** are triggered when the run
produces `fort.240` / `fort.209`. `imago.py` then calls
`processPOPTC` (which repeatedly invokes `imagoKKc`) to produce the
`.p.*` partial files. Relevant handlers: `_manage_poptc_spin`,
`_manage_poptc_nonspin`.

Job IDs: `scfoptc` = 104, `optc` = 204. Default bases in `JOB_DEFS`:
`-optc` scf FB / pscf EB; `-scfoptc` EB.

## Non-Linear Optical Properties (-nlop / -scfnlop)

Uses the same engine as the linear optical properties. `optc.F90`
runs in nonlinear mode (`doOPTC == 4`), reading its control
parameters (energy cutoff, max transition energy, delta, sigma) from
the `NLOP_INPUT_DATA` block via `readNlopControl` in `input.f90`. The
imaginary part chi2 is written to `fort.50`.

Because `job_id % 100 == 6`, `imago.py` invokes `imagoKKc` on
`fort.50` exactly as for `-optc`, giving chi1 (`fort.110`) and the
combined file (`fort.100`). `imagoKKc` labels its output headers
`Epsilon1` / `Epsilon2` / `ELF` regardless of whether the input was a
dielectric function or a susceptibility -- hence the "misleading"
headers in the chi files. Handler: `_manage_nlop_output`
(`imago.py:2164`); note it reuses the `optc` name tag.

Job IDs: `scfnlop` = 106, `nlop` = 206. Default bases: `-nlop` scf
FB / pscf EB; `-scfnlop` EB.

## Photo-Absorption Cross Section (-pacs / -scfpacs)

`optc.F90` in PACS mode (`doOPTC == 2`), broadened with `sigmaPACS`
in steps of `deltaPACS`. The energy onset is set from the
total-energy difference between the ground and core-excited states
(`energyMin = totalEnergyDiffPACS - mod(totalEnergyDiffPACS,
onsetEnergySlackPACS)`), not the orbital energy difference. The
spectrum is printed via `printSpectrum(0, ...)` (specType 0 =
XANES/ELNES) to `fort.50`. No Kramers-Kronig step
(`job_id % 100 == 5`). Handler: `_manage_pacs_output`
(`imago.py:2142`). Output format: `5e15.7`, header
`Energy totalXANES xXANES yXANES zXANES`.

Job IDs: `scfpacs` = 105, `pacs` = 205. Default bases: `-pacs` scf
FB / pscf EB; `-scfpacs` EB. An edge argument is required.

## Sigma(E) Curve (-sige / -scfsige)

`optc.F90` in Sigma(E) mode (`doOPTC == 3`) via the dedicated
`computeSigmaE` routine (`optc.F90:1850`). The conductivity
accumulator is scaled into `(micro-ohm-cm)^-1` and written to
`fort.50` with format `5e20.8e4` and **no header line**
(`optc.F90:2292`). A scalar total conductivity is written to the log
unit `fort.20` ("The total sigma is: ..."). No Kramers-Kronig step.
Handler: `_manage_sige_output` (`imago.py:2205`).

Job IDs: `scfsige` = 107, `sige` = 207. Default bases: FB / FB.

## Symmetric Band Structure (-sybd / -scfsybd)

After the eigenvalues are computed, `imago.py` (for `job_id % 100 ==
8`) runs `makeSYBD.py -dat fort.5 -out fort.20 -raw fort.31 -plot
fort.41` to format the raw band data into a `.raw` and a plottable
`.plot`. Spin-polarized runs make a second call with `-raw fort.32
-plot fort.42`. The k-path itself is set up earlier by `makeinput.py`
(`-sybdpath`).

Handler `_manage_sybd_output` (`imago.py:2228`): moves `fort.41` →
`[edge]_sybd-[basis].plot`; keeps `fort.31`/`fort.33` as
`[edge]_sybd-[basis].raw.31` / `.raw.33` in the working dir; and
copies `fort.33` to the phase dir as `vdim-[basis].raw` (the
`fn.vdim` tag). In the kyanite sample the `.plot` had 300 path
points × 318 band columns and `vdim-fb.raw` had 316 rows; the band
count and valence-basis dimension differ slightly.

Job IDs: `scfsybd` = 108, `sybd` = 208. Default bases: FB / FB.

## Quantum and Electronic Properties (-field / -scffield)

Handled at `job_id % 100 == 10` by `_manage_field_output`
(`imago.py:2280`), which moves scratch files only if they exist:
`fort.30`/`fort.31`/`fort.32` → `[edge]_field-[basis].prof-a/-b/-c.dat`
(1-D profiles along a, b, c); `fort.56` →
`[edge]_field-[basis].rho.dat` (charge-center file); `fort.78` →
`[edge]_field-[basis].xdmf3`.

The volumetric data are produced through the HDF5-based field path
(`hdf5PSCFField.F90`). The `.xdmf3` file is an XML index that
references a companion `[edge]_pscf+field-[basis].hdf5` file (which
holds the actual 3D arrays). HDF5 groups seen in the sample:
`meshGroup/mesh`, `psiRGroup`, `psiIGroup`, `wavGroup`, `rhoGroup`,
`potGroup`; each field has attributes `_live_up+dn`, `_live_up-dn`,
`_diff_up+dn`, `_neutral`. The kyanite sample used a 10x10x10 mesh;
`psi_i` was all zeros because it was a Gamma-point (real) run, and
the `.hdf5` file itself was not retained in the job directory.

Profile column header (sample):
`aPos psi_r_live_up+dn psi_i_live_up+dn wav_live_up+dn rho_live_up+dn pot_live_up+dn`.

Job IDs: `scffield` = 110, `field` = 210. Default bases: FB / FB.

## Local Environment Analysis (-loen)

Dedicated engine module `loen.f90` (`module O_LocalEnv`), entry
`bispec` → `computeBispectrumComponent`. Neighbors are gathered to a
radial cutoff with a switching function `f_c` (`neighSwitchFn`), and
`bsComp` / `bsCompSum` hold the per-channel and summed bispectrum
components. Written to `fort.21` (`loen.f90:879`): header
`  site#  element  species  type_sp  type_flat` + `2j_NNN` columns
(count `twoj2+1`, values `(twoj1+twoj2) - (i-1)*2`) + `total`; data
rows formatted `i7,x,a8,x,i8,x,i8,x,i10,x` then `e15.7` per component
and `e15.4` for the sum.

`imago.py` job ID is **311** (`JOB_DEFS` "loen", both bases None). It
shares property code 11 (`job_id % 100`) with mtop, so the dispatch
guards on `job_id == 311` / `job_id < 300` to avoid the mtop branch
(`_manage_mtop_output` would hunt for a `fort.180` loen never
writes). Handler `_manage_loen_output` (`imago.py:2331`) copies
`fort.21` → `[edge]_loen-[basis].plot.21` and moves it →
`[edge]_loen-[basis].plot`. No sample output exists in the repo.

## Inter-Atomic Forces (-force / -scfforce)

**Experimental / incomplete.** Engine module `forces.F90`
(`module O_Force`) builds the force integral matrices `valeValeF`,
`coreCoreF`, `coreValeF` (the position derivatives of the
Hamiltonian/overlap integrals). The only output is a diagnostic dump
(`forces.F90:834`) written to units `97+j` -- `fort.98` (spin up /
default), `fort.99` (spin down) -- as labeled blocks: `i(x,y,z)=`
(Cartesian component), `j(spin)=`, `k(kp)=`, then a `valeDim x
valeDim` matrix of `real(valeValeF(m,l,k,j,i))` in `e13.5`. There is
no finalized per-atom force-table format yet.

Handler `_manage_force_output` (`imago.py:2260`): `fort.98` →
`[edge]_force-[basis].dat` (spin: `.up.dat`), `fort.99` →
`.dn.dat`, with raw `.dat.98` / `.dat.99` copies. No sample output
exists in the repo.

Job IDs: `scfforce` = 109, `force` = 209. Default bases: FB / FB.
