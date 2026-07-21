#!/usr/bin/bash
for cell in full prim; do
  root=jobs/si_fingerprint/c109/$cell
  mkdir -p $root/atomicBDB/cache/structures

  # Clean element databases: pot1/coeff1 only, no harvested entries.
  cp -r share/atomicPDB $root/
  rm -f $root/atomicPDB/*/s_gaussian_pot.toml

  # Pre-seed the CIFs so neither run needs COD.  Only *.cif --
  # the skeletons are what we are varying.
  cp share/atomicBDB/cache/structures/*.cif \
     $root/atomicBDB/cache/structures/

  # The manifest: identical but for the cell.
  sed 's/^kpoint_integration = "gaussian"$/&\ncell = "'$cell'"/' \
      jobs/si_fingerprint/seed/manifest.toml > $root/manifest.toml
done
