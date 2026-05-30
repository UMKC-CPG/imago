"""kaleidoscope.builders -- the option-axis flight builders
(DESIGN 6.2 "flight-builder split"; 6.2.8).

A practical corollary of VISION Principle 12: the dispatch core
(``kaleidoscope.model`` / ``workspace`` / ``cache`` / ``dispatch``)
stays dumb -- it runs and tracks calculations without knowing what
they compute -- while the knowledge of HOW to lay out a particular
*option-axis sweep* lives in a domain-aware builder here.  Each
builder turns a structure plus fixed options into a ``Flight`` of
``CalcUnit``s sweeping one option axis.

Anticipated members of this family (DESIGN 6.2):

- ``predict_verify`` -- sweep k-point density as a predict-then-
  verify verification grid, consulting the historical-guidance
  database (the first builder, DESIGN 6.2.8 / 7.7).
- (future) an XANES target-atom sweep.
- (future) a basis-size sweep.

By contrast, *structure-axis* sweeps (supercell expansion,
LAMMPS-snapshot splitting, defect-site enumeration) are NOT in
kaleidoscope's scope -- they live upstream in structure control
and acquisition (DESIGN 6.2).

This ``__init__`` is intentionally lightweight: it imports none of
its builders, so importing the package costs nothing and pulling a
builder's physics-layer dependencies (e.g. ``guidance_db``) is an
explicit, opt-in act by the client::

    from kaleidoscope.builders.predict_verify import predict_settings
"""
