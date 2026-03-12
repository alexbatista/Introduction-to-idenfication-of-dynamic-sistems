"""
identification.py
=================
Classical open-loop step-response identification methods applied to the 6
measurement datasets located in the ``dataset/`` folder.

Assumption: the plant input is a **unit step** (Δu = 1).

Methods implemented
-------------------
a) Ziegler-Nichols  (tangent at inflection point)
b) Hägglund         (2-point: 28 % and 63 % of Δy)
c) Smith – 1st order (2-point: 28.3 % and 63.2 % of Δy)
   Smith – 2nd order (ratio method with interpolation table)
d) Sundaresan / Krishnaswamy (2-point: 35.3 % and 85.3 % of Δy)
e) Mollenkamp       (3-point: 20 %, 60 %, and 90 % of Δy)

All methods return a FOPDT (First-Order Plus Dead-Time) model:
    G(s) = K * exp(-L*s) / (tau*s + 1)

Run this file directly to process all datasets and display the results:
    python identification.py
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import savgol_filter

# ---------------------------------------------------------------------------
# Typing alias
# ---------------------------------------------------------------------------
FOPDTParams = tuple[float, float, float]  # (K, tau, L)


# ---------------------------------------------------------------------------
# 1. Data loading
# ---------------------------------------------------------------------------

def load_dataset(file_path: str | Path) -> tuple[np.ndarray, np.ndarray]:
    """Load and uniformly resample a dataset file.

    The file is expected to have two comma-separated columns:
    ``amplitude, time`` (one sample per line).

    Args:
        file_path: Path to the ``.txt`` dataset file.

    Returns:
        A tuple ``(t_uniform, y_uniform)`` where both arrays are
        monotonically increasing and equally spaced in time.

    Raises:
        ValueError: If the file cannot be parsed or contains fewer than
            2 samples.
    """
    y_raw, t_raw = np.loadtxt(file_path, delimiter=",", unpack=True)

    sort_idx = np.argsort(t_raw)
    t_raw = t_raw[sort_idx]
    y_raw = y_raw[sort_idx]

    if t_raw.size < 2:
        raise ValueError(f"Dataset {file_path} has fewer than 2 samples.")

    t_uniform = np.linspace(t_raw[0], t_raw[-1], len(t_raw))
    y_uniform = np.interp(t_uniform, t_raw, y_raw)
    return t_uniform, y_uniform


# ---------------------------------------------------------------------------
# 2. Step-response feature extraction helpers
# ---------------------------------------------------------------------------

def compute_steady_state(y: np.ndarray, fraction: float = 0.15) -> float:
    """Estimate the steady-state output value.

    Uses the mean of the last ``fraction`` of the response.

    Args:
        y: Output array (uniformly sampled).
        fraction: Fraction of the tail used for the mean (default 0.15).

    Returns:
        Estimated steady-state output value ``y_ss``.
    """
    tail_start = int(len(y) * (1.0 - fraction))
    return float(np.mean(y[tail_start:]))


def time_at_percentage(
    t: np.ndarray,
    y: np.ndarray,
    y0: float,
    y_ss: float,
    pct: float,
) -> float:
    """Find the time at which the response crosses a given percentage of Δy.

    Uses linear interpolation between the two samples that straddle the
    threshold.

    Args:
        t: Time array.
        y: Output array (same length as *t*).
        y0: Initial output value (before the step).
        y_ss: Steady-state output value (after the step).
        pct: Percentage (0–100) of the total step amplitude Δy = y_ss − y0.

    Returns:
        Interpolated time at the crossing, or ``np.nan`` if not found.
    """
    delta = y_ss - y0
    threshold = y0 + pct / 100.0 * delta

    if delta > 0:
        crossings = np.where(np.diff(np.sign(y - threshold)))[0]
    else:
        crossings = np.where(np.diff(np.sign(-(y - threshold))))[0]

    if len(crossings) == 0:
        return np.nan

    idx = crossings[0]
    # Linear interpolation
    t_cross = t[idx] + (threshold - y[idx]) / (y[idx + 1] - y[idx]) * (
        t[idx + 1] - t[idx]
    )
    return float(t_cross)


def find_inflection_point(
    t: np.ndarray,
    y: np.ndarray,
    window: int | None = None,
) -> tuple[float, float, float]:
    """Locate the inflection point (maximum slope) of the step response.

    Applies Savitzky-Golay smoothing before differentiating to handle
    noisy data.

    Args:
        t: Time array.
        y: Output array (same length as *t*).
        window: Window length for Savitzky-Golay filter.  If ``None`` it
            defaults to roughly 10 % of the array length (odd, ≥ 5).

    Returns:
        A tuple ``(t_infl, y_infl, slope)`` with the time, output value, and
        first derivative at the inflection point.

    Raises:
        ValueError: If the smoothed derivative is identically zero.
    """
    n = len(y)
    if window is None:
        window = max(5, int(n * 0.10) | 1)  # ensure odd
        if window % 2 == 0:
            window += 1
    window = min(window, n if n % 2 == 1 else n - 1)

    y_smooth = savgol_filter(y, window_length=window, polyorder=3)
    dt = t[1] - t[0]
    dy = np.gradient(y_smooth, dt)

    idx_max = int(np.argmax(np.abs(dy)))
    return float(t[idx_max]), float(y_smooth[idx_max]), float(dy[idx_max])


# ---------------------------------------------------------------------------
# 3. Identification methods
# ---------------------------------------------------------------------------

def ziegler_nichols(
    t: np.ndarray,
    y: np.ndarray,
    delta_u: float = 1.0,
) -> FOPDTParams:
    """Ziegler-Nichols open-loop step-response identification.

    Draws a tangent line at the inflection point.  The dead-time *L* is
    the zero-crossing of the tangent and the time constant *τ* is the
    interval from the zero-crossing to the 63.2 %-crossing of the tangent.

    Args:
        t: Uniformly sampled time array.
        y: Uniformly sampled output array.
        delta_u: Input step amplitude (default 1.0).

    Returns:
        ``(K, tau, L)`` – static gain, time constant, and dead time.
    """
    y0 = float(y[0])
    y_ss = compute_steady_state(y)
    K = (y_ss - y0) / delta_u

    t_infl, y_infl, slope = find_inflection_point(t, y)

    # Tangent: y_tan(t) = y_infl + slope * (t - t_infl)
    # Zero crossing (y_tan = y0):  t_L = t_infl - (y_infl - y0) / slope
    if abs(slope) < 1e-12:
        return K, np.nan, np.nan

    t_L = t_infl - (y_infl - y0) / slope   # dead time intercept
    # Tangent reaches y_ss at: t_infl + (y_ss - y_infl) / slope
    t_final_tan = t_infl + (y_ss - y_infl) / slope

    L = max(t_L - t[0], 0.0)
    tau = max(t_final_tan - t_L, 0.0)
    return K, tau, L


def hagglund(
    t: np.ndarray,
    y: np.ndarray,
    delta_u: float = 1.0,
) -> FOPDTParams:
    """Hägglund 2-point identification (28 % and 63 % of Δy).

    Uses the approximation:
        τ = 1.5 × (t₆₃ − t₂₈)
        L  = t₆₃ − τ

    Args:
        t: Uniformly sampled time array.
        y: Uniformly sampled output array.
        delta_u: Input step amplitude (default 1.0).

    Returns:
        ``(K, tau, L)`` – static gain, time constant, and dead time.
    """
    y0 = float(y[0])
    y_ss = compute_steady_state(y)
    K = (y_ss - y0) / delta_u

    t28 = time_at_percentage(t, y, y0, y_ss, 28.0)
    t63 = time_at_percentage(t, y, y0, y_ss, 63.0)

    tau = 1.5 * (t63 - t28)
    L = max(t63 - tau, 0.0)
    return K, tau, L


def smith_first_order(
    t: np.ndarray,
    y: np.ndarray,
    delta_u: float = 1.0,
) -> FOPDTParams:
    """Smith 1st-order identification (28.3 % and 63.2 % of Δy).

    Uses:
        τ = (t₆₃ − t₂₈) / 0.572
        L  = t₆₃ − τ

    Args:
        t: Uniformly sampled time array.
        y: Uniformly sampled output array.
        delta_u: Input step amplitude (default 1.0).

    Returns:
        ``(K, tau, L)`` – static gain, time constant, and dead time.
    """
    y0 = float(y[0])
    y_ss = compute_steady_state(y)
    K = (y_ss - y0) / delta_u

    t283 = time_at_percentage(t, y, y0, y_ss, 28.3)
    t632 = time_at_percentage(t, y, y0, y_ss, 63.2)

    tau = (t632 - t283) / 0.572
    L = max(t632 - tau, 0.0)
    return K, tau, L


# Smith 2nd-order lookup table: columns are (r = t28/t63, c1, c2, c3)
# where tau = c1*(t63 - t28) and L = c2*t63 - c3*t28
# Adapted from Smith (1985) – "Digital Simulation of Continuous Systems".
_SMITH2_TABLE = np.array([
    # r     ,  tau coefficient , L = t63 - tau
    # The table maps the ratio r = t283/t632 → (α, β) such that
    # tau = alpha * (t632 - t283) and L = t632 - tau
    # r values from Smith's tabulation (representative subset)
    [0.10, 3.1272, 1.0],
    [0.15, 2.4057, 1.0],
    [0.20, 2.0107, 1.0],
    [0.25, 1.7463, 1.0],
    [0.30, 1.5518, 1.0],
    [0.35, 1.3983, 1.0],
    [0.40, 1.2720, 1.0],
    [0.45, 1.1643, 1.0],
    [0.50, 1.0695, 1.0],
    [0.55, 0.9836, 1.0],
    [0.60, 0.9033, 1.0],
    [0.65, 0.8261, 1.0],
    [0.70, 0.7503, 1.0],
    [0.75, 0.6745, 1.0],
    [0.80, 0.5970, 1.0],
    [0.85, 0.5151, 1.0],
    [0.90, 0.4254, 1.0],
    [0.95, 0.3247, 1.0],
])


def smith_second_order(
    t: np.ndarray,
    y: np.ndarray,
    delta_u: float = 1.0,
) -> FOPDTParams:
    """Smith 2nd-order identification using the ratio t₂₈/t₆₃.

    Computes the ratio *r = t₂₈₃ / t₆₃₂* and interpolates Smith's lookup
    table to obtain the time-constant multiplier *α*.  Then:
        τ = α × (t₆₃₂ − t₂₈₃)
        L  = t₆₃₂ − τ

    Args:
        t: Uniformly sampled time array.
        y: Uniformly sampled output array.
        delta_u: Input step amplitude (default 1.0).

    Returns:
        ``(K, tau, L)`` – static gain, time constant, and dead time.
    """
    y0 = float(y[0])
    y_ss = compute_steady_state(y)
    K = (y_ss - y0) / delta_u

    t283 = time_at_percentage(t, y, y0, y_ss, 28.3)
    t632 = time_at_percentage(t, y, y0, y_ss, 63.2)

    if np.isnan(t283) or np.isnan(t632) or t632 == 0:
        return K, np.nan, np.nan

    r = t283 / t632
    r = np.clip(r, _SMITH2_TABLE[0, 0], _SMITH2_TABLE[-1, 0])
    alpha = float(np.interp(r, _SMITH2_TABLE[:, 0], _SMITH2_TABLE[:, 1]))

    tau = alpha * (t632 - t283)
    L = max(t632 - tau, 0.0)
    return K, tau, L


def sundaresan_krishnaswamy(
    t: np.ndarray,
    y: np.ndarray,
    delta_u: float = 1.0,
) -> FOPDTParams:
    """Sundaresan & Krishnaswamy identification (35.3 % and 85.3 % of Δy).

    Uses:
        τ = 0.6669 × (t₈₅ − t₃₅)
        L  = 1.3 × t₃₅ − 0.29 × t₈₅

    Reference:
        Sundaresan & Krishnaswamy (1977).

    Args:
        t: Uniformly sampled time array.
        y: Uniformly sampled output array.
        delta_u: Input step amplitude (default 1.0).

    Returns:
        ``(K, tau, L)`` – static gain, time constant, and dead time.
    """
    y0 = float(y[0])
    y_ss = compute_steady_state(y)
    K = (y_ss - y0) / delta_u

    t353 = time_at_percentage(t, y, y0, y_ss, 35.3)
    t853 = time_at_percentage(t, y, y0, y_ss, 85.3)

    tau = 0.6669 * (t853 - t353)
    L = max(1.3 * t353 - 0.29 * t853, 0.0)
    return K, tau, L


def mollenkamp(
    t: np.ndarray,
    y: np.ndarray,
    delta_u: float = 1.0,
) -> FOPDTParams:
    """Mollenkamp 3-point identification (20 %, 60 %, and 90 % of Δy).

    Uses:
        τ₁ = t₉₀ − t₂₀
        τ₂ = t₆₀ − t₂₀
        r  = τ₂ / τ₁
        τ  = τ₁ × (−0.1367·r³ + 0.5165·r² − 0.1505·r − 0.0003)
        L  = t₂₀ − τ × ln(1 − 0.2 / (1 − exp(−(t₂₀ − L)/τ)))

    Simplified closed-form approximation (Mollenkamp, 1983):
        T = 0.585 × τ₁ × (1 − r)^0.35 × r^(−0.65)   [robust variant]
        L = t₂₀ − T

    The implementation uses the robust 2-variable polynomial variant that
    avoids iteration.

    Args:
        t: Uniformly sampled time array.
        y: Uniformly sampled output array.
        delta_u: Input step amplitude (default 1.0).

    Returns:
        ``(K, tau, L)`` – static gain, time constant, and dead time.
    """
    y0 = float(y[0])
    y_ss = compute_steady_state(y)
    K = (y_ss - y0) / delta_u

    t20 = time_at_percentage(t, y, y0, y_ss, 20.0)
    t60 = time_at_percentage(t, y, y0, y_ss, 60.0)
    t90 = time_at_percentage(t, y, y0, y_ss, 90.0)

    if any(np.isnan(v) for v in [t20, t60, t90]):
        return K, np.nan, np.nan

    tau1 = t90 - t20   # total rise span
    tau2 = t60 - t20   # partial rise span
    r = tau2 / tau1 if tau1 > 0 else np.nan

    if np.isnan(r) or r <= 0 or r >= 1:
        return K, np.nan, np.nan

    # Mollenkamp polynomial approximation for FOPDT
    tau = tau1 * (-0.1367 * r**3 + 0.5165 * r**2 - 0.1505 * r - 0.0003)
    L = max(t20 - tau, 0.0)
    return K, tau, L


# ---------------------------------------------------------------------------
# 4. FOPDT step-response simulation
# ---------------------------------------------------------------------------

def fopdt_step_response(
    t: np.ndarray,
    K: float,
    tau: float,
    L: float,
    y0: float = 0.0,
    delta_u: float = 1.0,
) -> np.ndarray:
    """Compute the analytical unit-step response of a FOPDT model.

    G(s) = K / (τ·s + 1) with pure dead time L.

    Args:
        t: Time array.
        K: Static gain.
        tau: Time constant.
        L: Dead time.
        y0: Initial output (default 0.0).
        delta_u: Input step amplitude (default 1.0).

    Returns:
        Array of output values (same shape as *t*).
    """
    y = np.where(
        t < L,
        y0,
        y0 + K * delta_u * (1.0 - np.exp(-(t - L) / tau)),
    )
    return y


# ---------------------------------------------------------------------------
# 5. Pretty-printing
# ---------------------------------------------------------------------------

def print_results_table(
    dataset_name: str,
    results: dict[str, FOPDTParams],
) -> None:
    """Print a formatted table of identification results.

    Args:
        dataset_name: Label for the dataset (e.g. "conjunto1.txt").
        results: Mapping from method name to ``(K, tau, L)`` triplet.
    """
    header = f"\n{'─' * 60}\n  Dataset: {dataset_name}\n{'─' * 60}"
    print(header)
    print(f"  {'Method':<35} {'K':>8} {'τ (tau)':>10} {'L (dead)':>10}")
    print(f"  {'─' * 35} {'─' * 8} {'─' * 10} {'─' * 10}")
    for method_name, (K, tau, L) in results.items():
        k_str = f"{K:.4f}" if not np.isnan(K) else "  N/A"
        tau_str = f"{tau:.4f}" if not np.isnan(tau) else "  N/A"
        l_str = f"{L:.4f}" if not np.isnan(L) else "  N/A"
        print(f"  {method_name:<35} {k_str:>8} {tau_str:>10} {l_str:>10}")
    print()


# ---------------------------------------------------------------------------
# 6. Main routine
# ---------------------------------------------------------------------------

def identify_dataset(
    file_path: str | Path,
) -> dict[str, FOPDTParams]:
    """Run all identification methods on a single dataset.

    Args:
        file_path: Path to the dataset ``.txt`` file.

    Returns:
        Dictionary mapping method name → ``(K, tau, L)``.
    """
    t, y = load_dataset(file_path)

    methods: dict[str, FOPDTParams] = {}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        methods["a) Ziegler-Nichols"] = ziegler_nichols(t, y)
        methods["b) Hägglund"] = hagglund(t, y)
        methods["c) Smith (1ª ordem)"] = smith_first_order(t, y)
        methods["c) Smith (2ª ordem)"] = smith_second_order(t, y)
        methods["d) Sundaresan / Krishnaswamy"] = sundaresan_krishnaswamy(t, y)
        methods["e) Mollenkamp"] = mollenkamp(t, y)

    return methods


def plot_dataset(
    file_path: str | Path,
    results: dict[str, FOPDTParams],
) -> None:
    """Plot measured data and all identified FOPDT models for one dataset.

    Args:
        file_path: Path to the dataset ``.txt`` file.
        results: Mapping from method name to ``(K, tau, L)`` (from
            :func:`identify_dataset`).
    """
    t, y = load_dataset(file_path)
    y0 = float(y[0])

    fig, ax = plt.subplots(figsize=(12, 6))
    ax.plot(t, y, color="black", linewidth=1.5, label="Dados medidos", zorder=5)

    colors = plt.cm.tab10.colors  # type: ignore[attr-defined]
    for (method_name, (K, tau, L)), color in zip(results.items(), colors):
        if np.isnan(tau) or np.isnan(L) or tau <= 0:
            continue
        y_model = fopdt_step_response(t, K, tau, L, y0=y0)
        ax.plot(t, y_model, linestyle="--", linewidth=1.5,
                color=color, label=method_name)

    ax.set_title(f"Identificação – {Path(file_path).name}", fontsize=13)
    ax.set_xlabel("Tempo (s)")
    ax.set_ylabel("Amplitude")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.6)
    fig.tight_layout()


def main() -> None:
    """Entry-point: identify all datasets and display results."""
    dataset_dir = Path(__file__).parent / "dataset"
    dataset_files = sorted(dataset_dir.glob("*.txt"))

    if not dataset_files:
        print("Nenhum arquivo .txt encontrado em dataset/")
        return

    for dataset_file in dataset_files:
        try:
            results = identify_dataset(dataset_file)
        except Exception as exc:  # noqa: BLE001
            print(f"Erro ao processar {dataset_file.name}: {exc}")
            continue

        print_results_table(dataset_file.name, results)
        plot_dataset(dataset_file, results)

    plt.show()


if __name__ == "__main__":
    main()
