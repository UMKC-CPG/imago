"""test_makeinput_gamma.py -- the Gamma-point k-point sentinel.

Locks the rule from DESIGN 3.5: a shift is meaningful only when an
axis is sampled by more than one k-point, so Gamma is requested with
an explicit ``0`` sentinel rather than inferred from a ``1 1 1`` mesh.

  * ``-kp 0 0 0`` (and ``-scfkp``/``-pscfkp 0 0 0``) and ``-kpd 0``
    (and ``-scfkpd``/``-pscfkpd 0``) mark a group as Gamma: a single
    k-point at the origin with no shift, written canonically as a
    1x1x1 style-1 mesh so ``imago.py``'s ``check_gamma_kp`` routes it
    to the gamma-specialized executable (imagoG).
  * ``-kp 1 1 1`` is NOT Gamma -- it is one *shifted* point (a
    mean-value sample) that keeps the shift and runs on the general
    complex executable.
  * A zero mixed with positive counts (e.g. ``0 0 1``) is a fatal
    typo for the all-zero sentinel.

conftest.py's ``SCRIPTS_DIR`` insertion lets us import ``makeinput``
directly.  Each test builds a real ``ScriptSettings`` (its rc-sourced
defaults give the k-point state), overlays the relevant CLI option
through ``_args_from_options``, and runs ``reconcile`` -- the same
path a bare ``makeinput`` invocation takes.
"""

import pytest

import makeinput


def _reconcile(options):
    """Build a default ScriptSettings, overlay ``options`` onto an
    otherwise-bare argparse namespace, run reconcile, and return the
    settings so the test can inspect the resolved k-point state."""
    settings = makeinput.ScriptSettings()
    settings.reconcile(settings._args_from_options(options))
    return settings


def test_kp_000_is_gamma():
    """``-kp 0 0 0`` marks both groups Gamma."""
    settings = _reconcile({"kp": [0, 0, 0]})
    assert settings.kp_gamma[1] is True
    assert settings.kp_gamma[2] is True
    assert settings.kp_note[1] == "(Gamma)"
    assert settings.kp_note[2] == "(Gamma)"


def test_kp_111_is_not_gamma():
    """``-kp 1 1 1`` is a single *shifted* point, not Gamma, so it
    stays on the general executable (note ``(General)``)."""
    settings = _reconcile({"kp": [1, 1, 1]})
    assert settings.kp_gamma[1] is False
    assert settings.kp_gamma[2] is False
    assert settings.kp_note[1] == "(General)"
    assert settings.kp_note[2] == "(General)"


def test_per_group_scfkp_000_only_scf_gamma():
    """``-scfkp 0 0 0`` makes only the SCF group Gamma; the unset
    PSCF group keeps its general default."""
    settings = _reconcile({"scfkp": [0, 0, 0]})
    assert settings.kp_gamma[1] is True
    assert settings.kp_note[1] == "(Gamma)"
    assert settings.kp_gamma[2] is False
    assert settings.kp_note[2] == "(General)"


def test_mixed_zero_mesh_is_fatal():
    """A zero mixed with positive counts is rejected, not silently
    treated as Gamma."""
    with pytest.raises(SystemExit):
        _reconcile({"kp": [0, 0, 1]})


def test_kpd_0_is_gamma():
    """``-kpd 0`` is the density-mode Gamma sentinel."""
    settings = _reconcile({"kpd": 0.0})
    assert settings.kp_gamma[1] is True
    assert settings.kp_gamma[2] is True
    assert settings.kp_note[1] == "(Gamma)"
    assert settings.kp_note[2] == "(Gamma)"


def test_kpd_positive_is_density_not_gamma():
    """A positive ``-kpd`` is genuine density mode, not Gamma."""
    settings = _reconcile({"kpd": 5.0})
    assert settings.kp_gamma[1] is False
    assert settings.kp_gamma[2] is False
    assert settings.kp_note[1] == "(Density)"
    assert settings.kp_note[2] == "(Density)"
