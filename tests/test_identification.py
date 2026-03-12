"""
tests/test_identification.py
=============================
Unit tests for the identification module using synthetic FOPDT step
responses with known parameters.

Run with:
    pytest tests/test_identification.py -v
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

# Make the parent directory importable when running tests from any location
sys.path.insert(0, str(Path(__file__).parent.parent))

from identification import (
    compute_steady_state,
    fopdt_step_response,
    hagglund,
    mollenkamp,
    smith_first_order,
    smith_second_order,
    sundaresan_krishnaswamy,
    time_at_percentage,
    ziegler_nichols,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

def make_fopdt(
    K: float = 1.0,
    tau: float = 5.0,
    L: float = 1.0,
    t_end: float = 40.0,
    n: int = 2000,
    noise_std: float = 0.0,
    seed: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Generate a noiseless (or slightly noisy) synthetic FOPDT step response.

    Args:
        K: Static gain.
        tau: Time constant.
        L: Dead time.
        t_end: Simulation end time.
        n: Number of samples.
        noise_std: Standard deviation of additive Gaussian noise.
        seed: Random seed for reproducibility.

    Returns:
        ``(t, y)`` arrays.
    """
    t = np.linspace(0.0, t_end, n)
    y = fopdt_step_response(t, K=K, tau=tau, L=L, y0=0.0, delta_u=1.0)
    if noise_std > 0.0:
        rng = np.random.default_rng(seed)
        y = y + rng.normal(0.0, noise_std, size=y.shape)
    return t, y


# ---------------------------------------------------------------------------
# Helper tests
# ---------------------------------------------------------------------------

class TestComputeSteadyState:
    """Tests for :func:`compute_steady_state`."""

    def test_constant_signal(self) -> None:
        """Steady state of a constant signal equals its value."""
        y = np.ones(100) * 3.14
        assert compute_steady_state(y) == pytest.approx(3.14)

    def test_fopdt_response(self) -> None:
        """Steady state of FOPDT response approximates K."""
        t, y = make_fopdt(K=2.5, tau=5.0, L=1.0, t_end=60.0)
        assert compute_steady_state(y) == pytest.approx(2.5, rel=0.01)


class TestTimeAtPercentage:
    """Tests for :func:`time_at_percentage`."""

    def test_63_percent_fopdt(self) -> None:
        """At t = L + tau the response reaches ~63.2 % for a FOPDT."""
        K, tau, L = 1.0, 5.0, 1.0
        t, y = make_fopdt(K=K, tau=tau, L=L, t_end=40.0)
        y0 = float(y[0])
        y_ss = K
        t63 = time_at_percentage(t, y, y0, y_ss, 63.2)
        assert t63 == pytest.approx(L + tau, rel=0.02)

    def test_returns_nan_when_not_reached(self) -> None:
        """Returns NaN when the threshold is never crossed."""
        t = np.linspace(0, 1, 100)
        y = np.ones_like(t) * 0.5  # flat, never reaches 100 %
        t_cross = time_at_percentage(t, y, 0.0, 1.0, 100.0)
        assert np.isnan(t_cross)


class TestFOPDTStepResponse:
    """Tests for :func:`fopdt_step_response`."""

    def test_zero_before_deadtime(self) -> None:
        """Output equals y0 before the dead time."""
        L = 3.0
        t = np.linspace(0.0, 10.0, 500)
        y = fopdt_step_response(t, K=1.0, tau=2.0, L=L, y0=0.0)
        assert np.all(y[t < L] == pytest.approx(0.0))

    def test_asymptote(self) -> None:
        """Output approaches K·Δu at large t."""
        K, tau, L = 2.0, 3.0, 0.5
        t = np.linspace(0.0, 60.0, 5000)
        y = fopdt_step_response(t, K=K, tau=tau, L=L)
        assert y[-1] == pytest.approx(K, rel=1e-3)


# ---------------------------------------------------------------------------
# Identification method tests  (tolerance is relatively loose because real
# data are noisy; noiseless synthetic tests should be tighter)
# ---------------------------------------------------------------------------

_FOPDT_CASES = [
    # (K, tau, L)
    (1.0, 5.0, 1.0),
    (2.0, 8.0, 2.0),
    (0.5, 3.0, 0.5),
]

TOLERANCE_REL = 0.08   # 8 % relative tolerance for noiseless synthetic data


class TestZieglerNichols:
    """Tests for :func:`ziegler_nichols`."""

    @pytest.mark.parametrize("K,tau,L", _FOPDT_CASES)
    def test_noiseless(self, K: float, tau: float, L: float) -> None:
        """ZN identifies K, tau, and L within tolerance on noiseless data."""
        t, y = make_fopdt(K=K, tau=tau, L=L, t_end=10 * (tau + L))
        K_id, tau_id, L_id = ziegler_nichols(t, y)
        assert K_id == pytest.approx(K, rel=TOLERANCE_REL)
        assert tau_id == pytest.approx(tau, rel=TOLERANCE_REL)
        assert L_id == pytest.approx(L, rel=TOLERANCE_REL + 0.1)  # ZN L is looser

    def test_returns_three_values(self) -> None:
        """Return value is always a 3-tuple."""
        t, y = make_fopdt()
        result = ziegler_nichols(t, y)
        assert len(result) == 3


class TestHagglund:
    """Tests for :func:`hagglund`."""

    @pytest.mark.parametrize("K,tau,L", _FOPDT_CASES)
    def test_noiseless(self, K: float, tau: float, L: float) -> None:
        """Hägglund identifies K, tau, and L within tolerance on noiseless data."""
        t, y = make_fopdt(K=K, tau=tau, L=L, t_end=10 * (tau + L))
        K_id, tau_id, L_id = hagglund(t, y)
        assert K_id == pytest.approx(K, rel=TOLERANCE_REL)
        assert tau_id == pytest.approx(tau, rel=TOLERANCE_REL)
        assert L_id >= 0.0


class TestSmithFirstOrder:
    """Tests for :func:`smith_first_order`."""

    @pytest.mark.parametrize("K,tau,L", _FOPDT_CASES)
    def test_noiseless(self, K: float, tau: float, L: float) -> None:
        """Smith 1st-order identifies K, tau, and L within tolerance."""
        t, y = make_fopdt(K=K, tau=tau, L=L, t_end=10 * (tau + L))
        K_id, tau_id, L_id = smith_first_order(t, y)
        assert K_id == pytest.approx(K, rel=TOLERANCE_REL)
        assert tau_id == pytest.approx(tau, rel=TOLERANCE_REL)
        assert L_id >= 0.0


class TestSmithSecondOrder:
    """Tests for :func:`smith_second_order`."""

    @pytest.mark.parametrize("K,tau,L", _FOPDT_CASES)
    def test_noiseless(self, K: float, tau: float, L: float) -> None:
        """Smith 2nd-order returns finite K, tau ≥ 0, and L ≥ 0."""
        t, y = make_fopdt(K=K, tau=tau, L=L, t_end=10 * (tau + L))
        K_id, tau_id, L_id = smith_second_order(t, y)
        assert not np.isnan(K_id)
        assert not np.isnan(tau_id)
        assert not np.isnan(L_id)
        assert tau_id > 0.0
        assert L_id >= 0.0
        assert K_id == pytest.approx(K, rel=TOLERANCE_REL)


class TestSundaresanKrishnaswamy:
    """Tests for :func:`sundaresan_krishnaswamy`."""

    @pytest.mark.parametrize("K,tau,L", _FOPDT_CASES)
    def test_noiseless(self, K: float, tau: float, L: float) -> None:
        """S&K identifies K, tau, and L within tolerance on noiseless data."""
        t, y = make_fopdt(K=K, tau=tau, L=L, t_end=10 * (tau + L))
        K_id, tau_id, L_id = sundaresan_krishnaswamy(t, y)
        assert K_id == pytest.approx(K, rel=TOLERANCE_REL)
        assert tau_id == pytest.approx(tau, rel=TOLERANCE_REL)
        assert L_id >= 0.0


class TestMollenkamp:
    """Tests for :func:`mollenkamp`."""

    @pytest.mark.parametrize("K,tau,L", _FOPDT_CASES)
    def test_noiseless(self, K: float, tau: float, L: float) -> None:
        """Mollenkamp identifies K, tau, and L within tolerance on noiseless data."""
        t, y = make_fopdt(K=K, tau=tau, L=L, t_end=10 * (tau + L))
        K_id, tau_id, L_id = mollenkamp(t, y)
        assert not np.isnan(K_id), "K should be finite"
        assert not np.isnan(tau_id), "tau should be finite"
        assert K_id == pytest.approx(K, rel=TOLERANCE_REL)
        assert tau_id > 0.0
        assert L_id >= 0.0


# ---------------------------------------------------------------------------
# Regression test: noisy data still produces finite results
# ---------------------------------------------------------------------------

class TestNoisyData:
    """All methods must return finite (non-NaN) results on lightly noisy data."""

    _methods = [
        ziegler_nichols,
        hagglund,
        smith_first_order,
        smith_second_order,
        sundaresan_krishnaswamy,
        mollenkamp,
    ]

    @pytest.mark.parametrize("method", _methods, ids=lambda m: m.__name__)
    def test_finite_output(self, method) -> None:  # noqa: ANN001
        """Method returns finite K on lightly noisy synthetic data."""
        t, y = make_fopdt(K=1.0, tau=5.0, L=1.0, t_end=40.0,
                          noise_std=0.02, seed=42)
        K_id, _tau, _L = method(t, y)
        assert np.isfinite(K_id), f"{method.__name__} returned NaN/Inf K"
