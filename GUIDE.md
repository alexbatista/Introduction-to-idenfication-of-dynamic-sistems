# 📘 Beginner's Guide — `identification.py` & `generate_report.py`

> **Audience:** Complete Python beginner — every line is explained.
>
> **Goal:** After reading this guide you will understand *what* each line does, *why* it exists, and the control-engineering concepts behind it.

---

## Table of Contents

1. [Background Concepts](#1-background-concepts)
2. [`identification.py` — Line-by-Line](#2-identificationpy--line-by-line)
   - [Module Docstring](#21-module-docstring-lines-1-23)
   - [Imports](#22-imports-lines-25-32)
   - [Type Alias](#23-type-alias-line-37)
   - [Data Loading — `load_dataset`](#24-data-loading--load_dataset-lines-44-72)
   - [Step-Response Helpers](#25-step-response-feature-extraction-helpers-lines-79-171)
   - [Identification Methods](#26-identification-methods-lines-178-431)
   - [FOPDT Simulation — `fopdt_step_response`](#27-fopdt-step-response-simulation-lines-438-466)
   - [Pretty-Printing — `print_results_table`](#28-pretty-printing--print_results_table-lines-473-492)
   - [Main Routines](#29-main-routines-lines-499-583)
3. [`generate_report.py` — Line-by-Line](#3-generate_reportpy--line-by-line)
   - [Module Docstring](#31-module-docstring-lines-1-18)
   - [Imports](#32-imports-lines-20-32)
   - [Error Metrics — `compute_error_metrics`](#33-error-metrics--compute_error_metrics-lines-40-71)
   - [Method Descriptions Dictionary](#34-method-descriptions-dictionary-lines-78-116)
   - [Helper `_fmt`](#35-helper-_fmt-lines-124-136)
   - [Report Generator — `generate_report`](#36-report-generator--generate_report-lines-139-394)
   - [Main Routine](#37-main-routine-lines-397-407)
4. [Glossary](#4-glossary)

---

## 1. Background Concepts

Before diving into the code, here are the core ideas you will encounter:

### FOPDT Model

A **First-Order Plus Dead-Time** model is a simple mathematical way to describe how a system (like a heater, motor, or chemical reactor) responds to a sudden change in input (a *step*). Its transfer function is:

$$G(s) = \frac{K \cdot e^{-Ls}}{\tau s + 1}$$

| Symbol | Name | Meaning |
|--------|------|---------|
| **K** | Static gain | How much the output changes per unit of input change. |
| **τ** (tau) | Time constant | How fast the system responds (larger τ → slower). |
| **L** | Dead time | The delay before the system starts reacting. |

### Step Response

When you apply a sudden constant input (like flipping a switch), the output follows a characteristic curve that gradually rises from its initial value to a new steady-state value. Different identification methods measure specific points on this curve to estimate K, τ, and L.

### What the Code Does (Big Picture)

1. **`identification.py`** — Loads measurement data, applies six classical identification methods, computes the FOPDT parameters (K, τ, L), and plots the results.
2. **`generate_report.py`** — Re-uses `identification.py` to run all methods, computes error metrics (MSE, IAE, ISE, ITAE), and generates a Markdown report with tables and a comparative conclusion.

---

## 2. `identification.py` — Line-by-Line

### 2.1 Module Docstring (Lines 1–23)

```python
"""
identification.py
=================
Classical open-loop step-response identification methods applied to the 6
measurement datasets located in the ``dataset/`` folder.
...
"""
```

* **What:** A multi-line string enclosed in triple quotes (`"""`). It is the *docstring* for the entire module (file).
* **Why:** It describes the purpose of the file so that anyone reading the code (or using `help(identification)` in a Python shell) immediately understands what the module does.
* **Key info inside:** Lists the six methods implemented, states the FOPDT model equation, and shows how to run the file (`python identification.py`).

---

### 2.2 Imports (Lines 25–32)

```python
from __future__ import annotations
```

* **`from __future__ import annotations`** — This is a *future import*. It tells Python to treat all type hints as **strings** instead of evaluating them immediately. This lets you use modern type-hint syntax (like `str | Path`) even in older Python versions (3.9 and below). It must always be the **very first import** in a file.

```python
import warnings
```

* **`import warnings`** — Loads Python's built-in `warnings` module, which lets you suppress or filter non-critical warning messages. Used later to temporarily silence math warnings during identification.

```python
from pathlib import Path
```

* **`from pathlib import Path`** — `pathlib` is Python's built-in module for working with file system paths. `Path` is a class that lets you manipulate file paths in an OS-independent way (works on Windows, Linux, macOS). For example, `Path("dataset") / "file.txt"` produces the correct path for your operating system.

```python
import matplotlib.pyplot as plt
```

* **`import matplotlib.pyplot as plt`** — [Matplotlib](https://matplotlib.org/) is Python's most popular plotting library. `pyplot` is its MATLAB-like interface. The alias `plt` is the universal convention. You use it to create charts (e.g., `plt.plot(x, y)`).

```python
import numpy as np
```

* **`import numpy as np`** — [NumPy](https://numpy.org/) is the fundamental library for numerical computing in Python. It provides fast arrays (`np.ndarray`) and mathematical functions. The alias `np` is the universal convention.

```python
from scipy.signal import savgol_filter
```

* **`from scipy.signal import savgol_filter`** — [SciPy](https://scipy.org/) is a scientific computing library built on top of NumPy. `savgol_filter` is the *Savitzky-Golay filter*, a method for smoothing noisy data while preserving the shape of the signal (peaks, slopes). It is used here before computing derivatives to avoid amplifying noise.

---

### 2.3 Type Alias (Line 37)

```python
FOPDTParams = tuple[float, float, float]  # (K, tau, L)
```

* **`tuple[float, float, float]`** — A *type alias*. It says: "whenever you see `FOPDTParams` in the code, it means a tuple containing exactly three floating-point numbers."
* **Why:** Instead of writing `tuple[float, float, float]` in every function signature, you write `FOPDTParams`. This makes the code cleaner and communicates that the three numbers always represent **(K, τ, L)**.
* **`# (K, tau, L)`** — An inline comment clarifying the order of values.

---

### 2.4 Data Loading — `load_dataset` (Lines 44–72)

```python
def load_dataset(file_path: str | Path) -> tuple[np.ndarray, np.ndarray]:
```

* **`def`** — Keyword to define a function.
* **`load_dataset`** — The function name. Chosen to clearly describe what it does.
* **`file_path: str | Path`** — The function accepts one argument named `file_path`. The type hint `str | Path` means it can be either a regular string (`"dataset/conjunto1.txt"`) or a `Path` object.
* **`-> tuple[np.ndarray, np.ndarray]`** — The *return type hint*. It says the function returns a tuple of two NumPy arrays.

```python
    """Load and uniformly resample a dataset file.
    ...
    """
```

* **Docstring** adhering to PEP 257 (Google style). Contains `Args:`, `Returns:`, and `Raises:` sections — the standard way to document a Python function.

```python
    y_raw, t_raw = np.loadtxt(file_path, delimiter=",", unpack=True)
```

* **`np.loadtxt()`** — Reads a text file and returns its contents as NumPy arrays.
  * **`file_path`** — The file to read.
  * **`delimiter=","`** — Tells NumPy that values on each line are separated by commas (CSV format).
  * **`unpack=True`** — Transposes the result: instead of getting rows, you get **columns**. The file has two columns: `amplitude, time`. So `y_raw` receives the first column (amplitude) and `t_raw` receives the second column (time).
* **`y_raw, t_raw = ...`** — *Tuple unpacking*: the two arrays returned are assigned to two separate variables in one line.

```python
    sort_idx = np.argsort(t_raw)
    t_raw = t_raw[sort_idx]
    y_raw = y_raw[sort_idx]
```

* **`np.argsort(t_raw)`** — Returns an array of *indices* that would sort `t_raw` in ascending order. For example, if `t_raw = [3, 1, 2]`, then `argsort` returns `[1, 2, 0]` (the index of the smallest comes first).
* **`t_raw[sort_idx]`** — *Fancy indexing*: reorders the array according to the sorted indices.
* **Why:** The raw data might not be in chronological order. Sorting ensures time values go from smallest to largest.

```python
    if t_raw.size < 2:
        raise ValueError(f"Dataset {file_path} has fewer than 2 samples.")
```

* **`t_raw.size`** — The total number of elements in the array.
* **`raise ValueError(...)`** — Raises an *exception* (an error). If the file has fewer than 2 data points, the function cannot do anything useful, so it stops and tells the caller what went wrong.
* **`f"...{file_path}..."`** — An *f-string*: a string where `{expression}` is replaced by the value of that expression at runtime.

```python
    t_uniform = np.linspace(t_raw[0], t_raw[-1], len(t_raw))
```

* **`np.linspace(start, stop, num)`** — Creates an array of `num` **equally spaced** values between `start` and `stop` (inclusive on both ends).
  * `t_raw[0]` — First (smallest) time value.
  * `t_raw[-1]` — Last (largest) time value. The index `-1` means "the last element."
  * `len(t_raw)` — Number of elements (same count as the original data).
* **Why:** Real measurement data may not be evenly spaced in time. Many algorithms assume uniform spacing, so we create a new uniformly spaced time grid.

```python
    y_uniform = np.interp(t_uniform, t_raw, y_raw)
```

* **`np.interp(x, xp, fp)`** — *Linear interpolation*. For each value in `x` (our new uniform time points), it finds the corresponding `y` value by drawing straight lines between the original (`xp`, `fp`) data points.
* **Result:** `y_uniform` is the resampled output at uniform time intervals.

```python
    return t_uniform, y_uniform
```

* **`return`** — Sends the two arrays back to whoever called `load_dataset`. The caller receives them as a tuple.

---

### 2.5 Step-Response Feature Extraction Helpers (Lines 79–171)

These are utility functions used by the identification methods.

#### `compute_steady_state` (Lines 79–92)

```python
def compute_steady_state(y: np.ndarray, fraction: float = 0.15) -> float:
```

* **`fraction: float = 0.15`** — A parameter with a *default value*. If you call `compute_steady_state(y)` without specifying `fraction`, it defaults to `0.15` (15 %).

```python
    tail_start = int(len(y) * (1.0 - fraction))
```

* Calculates the index where the "tail" (last portion) of the array starts. For a 1000-element array with `fraction=0.15`, `tail_start = int(1000 * 0.85) = 850`.
* **`int(...)`** — Converts to an integer because array indices must be whole numbers.

```python
    return float(np.mean(y[tail_start:]))
```

* **`y[tail_start:]`** — *Slicing*: takes all elements from index `tail_start` to the end.
* **`np.mean(...)`** — Computes the arithmetic average.
* **Why:** The steady-state value (`y_ss`) is the value the system settles at after its transient dies out. Averaging the last 15 % of the data gives a robust estimate, especially if there is noise.

#### `time_at_percentage` (Lines 95–133)

```python
def time_at_percentage(
    t: np.ndarray,
    y: np.ndarray,
    y0: float,
    y_ss: float,
    pct: float,
) -> float:
```

* Multiple parameters:
  * `t` — Time array.
  * `y` — Output array.
  * `y0` — Initial output (before the step).
  * `y_ss` — Steady-state output.
  * `pct` — Percentage (0–100) of the total change at which we want to find the crossing time.

```python
    delta = y_ss - y0
    threshold = y0 + pct / 100.0 * delta
```

* **`delta`** — The total change in output (amplitude of the step response).
* **`threshold`** — The actual `y` value that corresponds to `pct` % of the total change. For example, 63 % of a `delta` of 10 starting from `y0=0` gives a threshold of 6.3.

```python
    if delta > 0:
        crossings = np.where(np.diff(np.sign(y - threshold)))[0]
    else:
        crossings = np.where(np.diff(np.sign(-(y - threshold))))[0]
```

* **`y - threshold`** — Subtracts the threshold from every element. Values below the threshold become negative; values above become positive.
* **`np.sign(...)`** — Returns -1, 0, or +1 for each element.
* **`np.diff(...)`** — Computes differences between consecutive elements. When the sign changes (from negative to positive or vice-versa), the difference is non-zero — that is a *crossing*.
* **`np.where(...)[0]`** — Returns the indices where the condition is true (i.e., where a crossing occurred).
* **The `if/else`:** If the response goes *up* (`delta > 0`) we look for upward crossings; if it goes *down* we negate the expression so the logic works the same way.

```python
    if len(crossings) == 0:
        return np.nan
```

* If no crossing was found, return `np.nan` (*Not a Number*) — a special float value meaning "undefined" or "missing."

```python
    idx = crossings[0]
    t_cross = t[idx] + (threshold - y[idx]) / (y[idx + 1] - y[idx]) * (
        t[idx + 1] - t[idx]
    )
    return float(t_cross)
```

* **`crossings[0]`** — Takes the *first* crossing index (the earliest one in time).
* **Linear interpolation formula:** Between the two samples that straddle the threshold (`y[idx]` below and `y[idx+1]` above), we compute the exact fractional time where the threshold is crossed using the slope between those two points. This gives a more accurate time than simply picking the nearest sample.

#### `find_inflection_point` (Lines 136–171)

```python
def find_inflection_point(
    t: np.ndarray,
    y: np.ndarray,
    window: int | None = None,
) -> tuple[float, float, float]:
```

* **`window: int | None = None`** — The Savitzky-Golay filter window size. `None` means "auto-calculate."

```python
    n = len(y)
    if window is None:
        window = max(5, int(n * 0.10) | 1)  # ensure odd
        if window % 2 == 0:
            window += 1
```

* **`int(n * 0.10) | 1`** — `|` here is the *bitwise OR* operator. ORing with `1` forces the last bit to 1, which makes any even number odd. The Savitzky-Golay filter requires an **odd** window length.
* **`max(5, ...)`** — Ensures the window is at least 5 elements wide.
* The next `if` is a safety check: if the bitwise trick didn't work (it always does, but the author was cautious), explicitly make it odd.

```python
    window = min(window, n if n % 2 == 1 else n - 1)
```

* The window cannot be larger than the data. If the data length `n` is even, use `n - 1` (the largest odd number ≤ n).

```python
    y_smooth = savgol_filter(y, window_length=window, polyorder=3)
```

* **`savgol_filter`** — Smooths the data by fitting a polynomial of degree `polyorder=3` (cubic) to successive windows of `window_length` data points. This removes high-frequency noise while preserving the overall shape.

```python
    dt = t[1] - t[0]
```

* **`dt`** — The time step between consecutive samples (uniform spacing assumed).

```python
    dy = np.gradient(y_smooth, dt)
```

* **`np.gradient(array, spacing)`** — Computes the *numerical derivative* (rate of change) of the smoothed signal. Each element of `dy` is the slope of the signal at that point.

```python
    idx_max = int(np.argmax(np.abs(dy)))
```

* **`np.abs(dy)`** — Absolute value of each element (we want the steepest slope regardless of sign).
* **`np.argmax(...)`** — Returns the **index** of the maximum value. This is the inflection point — the point of maximum rate of change.

```python
    return float(t[idx_max]), float(y_smooth[idx_max]), float(dy[idx_max])
```

* Returns:
  1. The **time** at the inflection point.
  2. The **smoothed output value** at that point.
  3. The **slope** (first derivative) at that point.

---

### 2.6 Identification Methods (Lines 178–431)

All six methods follow the same pattern: compute K, find characteristic time points, then calculate τ and L using method-specific formulas.

#### a) `ziegler_nichols` (Lines 178–214)

```python
def ziegler_nichols(
    t: np.ndarray,
    y: np.ndarray,
    delta_u: float = 1.0,
) -> FOPDTParams:
```

* **`delta_u: float = 1.0`** — The amplitude of the input step. Default is 1 (unit step).
* **`-> FOPDTParams`** — Returns `(K, tau, L)`.

```python
    y0 = float(y[0])
    y_ss = compute_steady_state(y)
    K = (y_ss - y0) / delta_u
```

* **`y0`** — Initial output (first sample).
* **`y_ss`** — Steady-state output (estimated from the tail).
* **`K`** — Static gain: total output change divided by input change.

```python
    t_infl, y_infl, slope = find_inflection_point(t, y)
```

* Finds the inflection point: the time, output value, and slope at the steepest part of the response.

```python
    if abs(slope) < 1e-12:
        return K, np.nan, np.nan
```

* **Safety check:** If the slope is essentially zero (flat line), the tangent-line method cannot work. Returns `NaN` for τ and L.
* **`1e-12`** — Scientific notation for `0.000000000001`. A very small number used as a threshold to detect "practically zero."

```python
    t_L = t_infl - (y_infl - y0) / slope   # dead time intercept
```

* **Tangent line equation:** At the inflection point, the tangent line is `y_tan(t) = y_infl + slope × (t − t_infl)`. Setting `y_tan = y0` (the initial value) and solving for `t` gives the dead-time intercept: the time when the tangent crosses the initial output level.

```python
    t_final_tan = t_infl + (y_ss - y_infl) / slope
```

* Setting `y_tan = y_ss` (steady state) gives the time when the tangent reaches the final value.

```python
    L = max(t_L - t[0], 0.0)
    tau = max(t_final_tan - t_L, 0.0)
    return K, tau, L
```

* **`L`** — Dead time = time from the start of the step until the tangent crosses `y0`. `max(..., 0.0)` ensures it cannot be negative (physically, dead time ≥ 0).
* **`tau`** — Time constant = time from the zero-crossing to the steady-state crossing of the tangent.

#### b) `hagglund` (Lines 217–245)

```python
    t28 = time_at_percentage(t, y, y0, y_ss, 28.0)
    t63 = time_at_percentage(t, y, y0, y_ss, 63.0)
```

* Finds the times when the response reaches 28 % and 63 % of its total change.

```python
    tau = 1.5 * (t63 - t28)
    L = max(t63 - tau, 0.0)
```

* These are the Hägglund formulas. The coefficient `1.5` was derived empirically for first-order plus dead-time systems.

#### c) `smith_first_order` (Lines 248–276)

```python
    t283 = time_at_percentage(t, y, y0, y_ss, 28.3)
    t632 = time_at_percentage(t, y, y0, y_ss, 63.2)
```

* Uses 28.3 % and 63.2 % — these come from the exact time-constant definition for a pure first-order system (1 − e⁻¹ ≈ 0.632).

```python
    tau = (t632 - t283) / 0.572
    L = max(t632 - tau, 0.0)
```

* **`0.572`** — A coefficient specific to Smith's first-order method. It converts the time difference between the two reference points into the actual time constant.

#### c) `smith_second_order` (Lines 279–344)

```python
_SMITH2_TABLE = np.array([
    [0.10, 3.1272, 1.0],
    [0.15, 2.4057, 1.0],
    ...
    [0.95, 0.3247, 1.0],
])
```

* **Lookup table:** A NumPy 2D array. Each row contains `[r, alpha, 1.0]`. This table comes from Smith's 1985 publication. For a given ratio `r = t₂₈₃ / t₆₃₂`, the table gives the coefficient `alpha` used to compute τ.

```python
    r = t283 / t632
    r = np.clip(r, _SMITH2_TABLE[0, 0], _SMITH2_TABLE[-1, 0])
```

* **`np.clip(value, min, max)`** — Constrains `r` to stay within the table's valid range (0.10 to 0.95). If `r` is below 0.10 it becomes 0.10; if above 0.95 it becomes 0.95.

```python
    alpha = float(np.interp(r, _SMITH2_TABLE[:, 0], _SMITH2_TABLE[:, 1]))
```

* **`_SMITH2_TABLE[:, 0]`** — Slicing: all rows, column 0 (the `r` values).
* **`_SMITH2_TABLE[:, 1]`** — All rows, column 1 (the `alpha` values).
* **`np.interp`** — Linearly interpolates between table rows to find `alpha` for our specific `r`.

```python
    tau = alpha * (t632 - t283)
    L = max(t632 - tau, 0.0)
```

#### d) `sundaresan_krishnaswamy` (Lines 347–378)

```python
    t353 = time_at_percentage(t, y, y0, y_ss, 35.3)
    t853 = time_at_percentage(t, y, y0, y_ss, 85.3)
```

* Uses 35.3 % and 85.3 % — specific to this method.

```python
    tau = 0.6669 * (t853 - t353)
    L = max(1.3 * t353 - 0.29 * t853, 0.0)
```

* The coefficients `0.6669`, `1.3`, and `0.29` are from the original 1977 paper by Sundaresan & Krishnaswamy.

#### e) `mollenkamp` (Lines 381–431)

```python
    t20 = time_at_percentage(t, y, y0, y_ss, 20.0)
    t60 = time_at_percentage(t, y, y0, y_ss, 60.0)
    t90 = time_at_percentage(t, y, y0, y_ss, 90.0)
```

* A 3-point method: 20 %, 60 %, and 90 %.

```python
    if any(np.isnan(v) for v in [t20, t60, t90]):
        return K, np.nan, np.nan
```

* **`any(...)`** — Returns `True` if at least one element is `True`.
* **Generator expression `np.isnan(v) for v in [t20, t60, t90]`** — Checks each of the three values. If any is `NaN` (crossing not found), the method cannot proceed.

```python
    tau1 = t90 - t20   # total rise span
    tau2 = t60 - t20   # partial rise span
    r = tau2 / tau1 if tau1 > 0 else np.nan
```

* **Conditional expression** (`x if condition else y`): if `tau1` is positive, compute the ratio; otherwise set it to `NaN` to avoid division by zero.

```python
    if np.isnan(r) or r <= 0 or r >= 1:
        return K, np.nan, np.nan
```

* Guards against invalid ratios. Physically, `r` should be between 0 and 1.

```python
    tau = tau1 * (-0.1367 * r**3 + 0.5165 * r**2 - 0.1505 * r - 0.0003)
    L = max(t20 - tau, 0.0)
```

* **Polynomial approximation** from Mollenkamp (1983). The expression `-0.1367·r³ + 0.5165·r² − 0.1505·r − 0.0003` is a cubic fitted to mapping tables, avoiding the need for iterative solving.
* **`r**3`** — Python's exponentiation operator (`**`). Computes r³.

---

### 2.7 FOPDT Step-Response Simulation (Lines 438–466)

```python
def fopdt_step_response(
    t: np.ndarray,
    K: float,
    tau: float,
    L: float,
    y0: float = 0.0,
    delta_u: float = 1.0,
) -> np.ndarray:
```

* Computes what the theoretical FOPDT model predicts, so we can compare it against the actual data.

```python
    y = np.where(
        t < L,
        y0,
        y0 + K * delta_u * (1.0 - np.exp(-(t - L) / tau)),
    )
```

* **`np.where(condition, value_if_true, value_if_false)`** — An element-wise "if-else" for arrays:
  * **If `t < L`** (time is before the dead time ends): the output stays at `y0` — the system hasn't reacted yet.
  * **Otherwise:** the first-order step response formula:  `y0 + K·Δu·(1 − e^(−(t−L)/τ))`.
* **`np.exp(x)`** — The exponential function eˣ.

```python
    return y
```

---

### 2.8 Pretty-Printing — `print_results_table` (Lines 473–492)

```python
def print_results_table(
    dataset_name: str,
    results: dict[str, FOPDTParams],
) -> None:
```

* **`-> None`** — The function does not return anything; it only prints to the console.
* **`dict[str, FOPDTParams]`** — A dictionary where keys are strings (method names) and values are `(K, tau, L)` tuples.

```python
    header = f"\n{'─' * 60}\n  Dataset: {dataset_name}\n{'─' * 60}"
    print(header)
```

* **`'─' * 60`** — String repetition: creates a line of 60 "─" characters. This is used as a visual separator.

```python
    print(f"  {'Method':<35} {'K':>8} {'τ (tau)':>10} {'L (dead)':>10}")
    print(f"  {'─' * 35} {'─' * 8} {'─' * 10} {'─' * 10}")
```

* **`{'Method':<35}`** — f-string formatting: left-align (`<`) the text `"Method"` in a field 35 characters wide.
* **`{'K':>8}`** — Right-align (`>`) `"K"` in 8 characters.

```python
    for method_name, (K, tau, L) in results.items():
```

* **`results.items()`** — Iterates over key-value pairs of the dictionary.
* **`(K, tau, L)`** — Destructuring: the tuple value is unpacked into three variables in one step.

```python
        k_str = f"{K:.4f}" if not np.isnan(K) else "  N/A"
```

* **`{K:.4f}`** — Formats the number with 4 decimal places (e.g., `1.2345`).
* **Conditional expression:** If `K` is `NaN`, show `"N/A"` instead.

---

### 2.9 Main Routines (Lines 499–583)

#### `identify_dataset` (Lines 499–523)

```python
def identify_dataset(
    file_path: str | Path,
) -> dict[str, FOPDTParams]:
```

* Runs all six identification methods on one dataset and returns a dictionary of results.

```python
    t, y = load_dataset(file_path)
```

* Loads and resamples the data.

```python
    methods: dict[str, FOPDTParams] = {}
```

* Creates an empty dictionary to store results.

```python
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        methods["a) Ziegler-Nichols"] = ziegler_nichols(t, y)
        methods["b) Hägglund"] = hagglund(t, y)
        ...
```

* **`with warnings.catch_warnings():`** — A *context manager*. Any warning settings changed inside this `with` block are automatically restored when the block ends.
* **`warnings.simplefilter("ignore")`** — Suppresses all warnings. Some numerical operations (like dividing by very small numbers) can trigger warnings that are not errors — we silence them here.
* Each identification method is called and its result stored in the dictionary.

```python
    return methods
```

#### `plot_dataset` (Lines 526–556)

```python
    fig, ax = plt.subplots(figsize=(12, 6))
```

* **`plt.subplots()`** — Creates a new figure and axes (the area where the plot lives).
* **`figsize=(12, 6)`** — Width 12 inches, height 6 inches.
* **`fig`** — The figure object (container).
* **`ax`** — The axes object (the actual plot area). You call `ax.plot(...)` to draw on it.

```python
    ax.plot(t, y, color="black", linewidth=1.5, label="Dados medidos", zorder=5)
```

* Plots the raw measured data as a black line.
* **`label="Dados medidos"`** — The legend label (Portuguese for "Measured data").
* **`zorder=5`** — Drawing order: higher numbers are drawn on top of lower numbers.

```python
    colors = plt.cm.tab10.colors
```

* **`plt.cm.tab10`** — A *colormap* (a set of 10 distinct colors). `.colors` gives the actual RGB values.

```python
    for (method_name, (K, tau, L)), color in zip(results.items(), colors):
```

* **`zip(...)`** — Pairs each method's result with a color from the colormap.

```python
        if np.isnan(tau) or np.isnan(L) or tau <= 0:
            continue
```

* **`continue`** — Skips the rest of this loop iteration. If the method returned invalid results, don't plot it.

```python
        y_model = fopdt_step_response(t, K, tau, L, y0=y0)
        ax.plot(t, y_model, linestyle="--", linewidth=1.5,
                color=color, label=method_name)
```

* Computes the theoretical response and plots it as a dashed (`"--"`) line.

```python
    ax.set_title(f"Identificação – {Path(file_path).name}", fontsize=13)
    ax.set_xlabel("Tempo (s)")
    ax.set_ylabel("Amplitude")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(True, linestyle="--", alpha=0.6)
    fig.tight_layout()
```

* **`set_title`/`set_xlabel`/`set_ylabel`** — Set the plot's title and axis labels.
* **`ax.legend()`** — Adds a legend box using the `label` strings from each `plot()` call.
* **`ax.grid(True, ...)`** — Enables grid lines. `alpha=0.6` makes them semi-transparent.
* **`fig.tight_layout()`** — Adjusts spacing so nothing overlaps.

#### `main` (Lines 559–582)

```python
def main() -> None:
    dataset_dir = Path(__file__).parent / "dataset"
```

* **`Path(__file__)`** — The path to the current Python file (`identification.py`).
* **`.parent`** — The directory containing the file.
* **`/ "dataset"`** — Appends `"dataset"` to the path (the `/` operator is overloaded for `Path` objects). Result: the `dataset` folder in the same directory.

```python
    dataset_files = sorted(dataset_dir.glob("*.txt"))
```

* **`.glob("*.txt")`** — Finds all files matching the pattern `*.txt` (any `.txt` file).
* **`sorted(...)`** — Sorts the list alphabetically.

```python
    if not dataset_files:
        print("Nenhum arquivo .txt encontrado em dataset/")
        return
```

* If no files were found, print a message and exit.

```python
    for dataset_file in dataset_files:
        try:
            results = identify_dataset(dataset_file)
        except Exception as exc:
            print(f"Erro ao processar {dataset_file.name}: {exc}")
            continue
```

* **`try/except`** — *Exception handling*. If `identify_dataset` raises any error, the program doesn't crash; instead it prints the error and moves on to the next file.
* **`except Exception as exc`** — Catches any exception and stores it in the variable `exc`.

```python
        print_results_table(dataset_file.name, results)
        plot_dataset(dataset_file, results)
    plt.show()
```

* For each file: print the results table, then create the plot.
* **`plt.show()`** — Displays all plots at once in interactive windows. Called once at the end so all plots appear together.

```python
if __name__ == "__main__":
    main()
```

* **`if __name__ == "__main__":`** — A Python idiom. When you run `python identification.py` directly, `__name__` is set to `"__main__"`, so `main()` runs. If another file *imports* this module, `__name__` will be the module name instead, and `main()` will **not** run automatically — this lets other files use the functions without executing the main script.

---

## 3. `generate_report.py` — Line-by-Line

### 3.1 Module Docstring (Lines 1–18)

```python
"""
generate_report.py
==================
Runs all identification methods from ``identification.py`` on every dataset in
``dataset/`` and produces a Markdown report ...
"""
```

* Same structure as `identification.py`. Documents the purpose, the error metrics, and how to run (`python generate_report.py`).

---

### 3.2 Imports (Lines 20–32)

```python
from __future__ import annotations
from pathlib import Path
import numpy as np
```

* Same as in `identification.py` (see [Section 2.2](#22-imports-lines-25-32)).

```python
from identification import (
    FOPDTParams,
    compute_steady_state,
    fopdt_step_response,
    identify_dataset,
    load_dataset,
)
```

* **`from identification import ...`** — Imports specific names from the `identification` module (the other file in this project). This avoids code duplication: instead of copying functions, we reuse them.
* The parentheses let you split the import across multiple lines for readability.

---

### 3.3 Error Metrics — `compute_error_metrics` (Lines 40–71)

```python
def compute_error_metrics(
    t: np.ndarray,
    y_measured: np.ndarray,
    y_model: np.ndarray,
) -> dict[str, float]:
```

* Returns a dictionary mapping metric name → value.

```python
    error = y_measured - y_model
    abs_error = np.abs(error)
```

* **`error`** — Element-wise difference between measured and modeled data. Positive means the model *underestimates*.
* **`abs_error`** — Absolute values (removes sign, keeping magnitude only).

```python
    _trapz = getattr(np, "trapezoid", getattr(np, "trapz", None))
```

* **`getattr(object, name, default)`** — Gets an attribute (function) from an object by name. If the attribute doesn't exist, returns `default`.
* **Why the fallback?** NumPy ≥ 2.0 renamed `np.trapz` to `np.trapezoid`. This line checks for the new name first; if it doesn't exist (older NumPy), falls back to the old name. This makes the code compatible with both old and new versions.
* **`_trapz`** — The underscore prefix is a convention meaning "this variable is for internal use."

```python
    mse = float(np.mean(error ** 2))
```

* **MSE (Mean Squared Error):** Square each error, then take the mean. Penalizes large errors heavily (because of squaring).

```python
    iae = float(_trapz(abs_error, t))
```

* **IAE (Integral of Absolute Error):** Numerically integrates |error| over time using the *trapezoidal rule* — a method that approximates the area under a curve by breaking it into trapezoids.

```python
    ise = float(_trapz(error ** 2, t))
```

* **ISE (Integral of Squared Error):** Like IAE but integrates error² — penalizes large errors even more.

```python
    itae = float(_trapz(t * abs_error, t))
```

* **ITAE (Integral of Time × Absolute Error):** Multiplies |error| by time before integrating. This means errors that occur *later* in time are penalized more. A model that settles quickly will have low ITAE.

```python
    return {"MSE": mse, "IAE": iae, "ISE": ise, "ITAE": itae}
```

---

### 3.4 Method Descriptions Dictionary (Lines 78–116)

```python
METHOD_DESCRIPTIONS: dict[str, str] = {
    "a) Ziegler-Nichols": (
        "O método de **Ziegler-Nichols** (malha aberta) traça uma reta tangente ..."
    ),
    ...
}
```

* **`dict[str, str]`** — A dictionary where both keys and values are strings.
* Each entry maps a method name to its Portuguese description (with Markdown bold markup `**...**`).
* The parentheses around the string values allow multi-line strings without triple quotes. Python automatically concatenates adjacent string literals.

---

### 3.5 Helper `_fmt` (Lines 124–136)

```python
def _fmt(value: float, decimals: int = 6) -> str:
```

* **`_fmt`** — A helper function for formatting numbers. The leading underscore (`_`) is a convention meaning "private use only" — it's not intended to be called from outside this file.

```python
    if np.isnan(value):
        return "N/A"
    return f"{value:.{decimals}f}"
```

* **`f"{value:.{decimals}f}"`** — Nested f-string formatting. The inner `{decimals}` is replaced first (e.g., `6`), resulting in `f"{value:.6f}"`, which formats the number with 6 decimal places.

---

### 3.6 Report Generator — `generate_report` (Lines 139–394)

This is the largest function. It builds the entire Markdown report as a list of strings.

```python
def generate_report() -> str:
```

* Returns the complete report as a single string.

```python
    dataset_dir = Path(__file__).parent / "dataset"
    dataset_files = sorted(dataset_dir.glob("*.txt"))
```

* Same pattern as in `identification.py` — finds all `.txt` dataset files.

```python
    lines: list[str] = []
```

* **Strategy:** Build the report by appending strings to a list. Joining a list at the end (`"".join(lines)`) is much faster than concatenating strings with `+` repeatedly in Python.

#### Report Header (Lines 150–169)

```python
    lines.append("# Relatório de Identificação de Sistemas Dinâmicos\n\n")
```

* **`lines.append(...)`** — Adds a string to the end of the list.
* **`#`** — Markdown heading level 1.
* **`\n\n`** — Two newlines create a blank line in Markdown (paragraph break).

```python
    lines.append("$$G(s) = \\frac{K \\cdot e^{-Ls}}{\\tau s + 1}$$\n\n")
```

* **`$$...$$`** — LaTeX math block in Markdown. Renders the FOPDT equation.
* **`\\frac{...}{...}`** — LaTeX fraction. The double backslash (`\\`) is needed because Python itself uses `\` as an escape character — `\\` becomes a single `\` in the output string.

#### Method Descriptions (Lines 174–178)

```python
    for method_name, description in METHOD_DESCRIPTIONS.items():
        lines.append(f"### {method_name}\n\n")
        lines.append(f"{description}\n\n")
```

* Loops over the descriptions dictionary and adds each method as a subsection.

#### Error Metrics Description (Lines 183–205)

```python
    lines.append(
        "| Critério | Fórmula | Descrição |\n"
        "|----------|---------|----------|\n"
        ...
    )
```

* Markdown table syntax:
  * `|` separates columns.
  * The second row (`|---|---|---|`) defines column alignment (dashes = default left-align).
  * **Adjacent string literals** (no comma between them) are automatically concatenated by Python.

#### Per-Dataset Results (Lines 210–269)

```python
    summary_rows: list[dict[str, str | float]] = []
```

* A list to collect summary data for the comparative section later. Each element is a dictionary with keys like `"dataset"`, `"method"`, `"MSE"`, etc.

```python
    for ds_file in dataset_files:
        ds_name = ds_file.name
```

* **`.name`** — Property of a `Path` object that returns just the filename (e.g., `"conjunto1.txt"` without the directory).

```python
        t, y = load_dataset(ds_file)
        y0 = float(y[0])
        results = identify_dataset(ds_file)
```

* Loads data and runs all identification methods.

```python
        for method_name, (K, tau, L) in results.items():
            lines.append(
                f"| {method_name} | {_fmt(K, 4)} | {_fmt(tau, 4)} | {_fmt(L, 4)} |\n"
            )
```

* Fills in each row of the parameter table using `_fmt` to format numbers.

```python
        for method_name, (K, tau, L) in results.items():
            if np.isnan(tau) or np.isnan(L) or tau <= 0:
                lines.append(
                    f"| {method_name} | N/A | N/A | N/A | N/A |\n"
                )
                continue
```

* Skips error computation for invalid models.

```python
            y_model = fopdt_step_response(t, K, tau, L, y0=y0)
            metrics = compute_error_metrics(t, y, y_model)
```

* Simulates the model and computes all four error metrics.

```python
            summary_rows.append({
                "dataset": ds_name,
                "method": method_name,
                "MSE": metrics["MSE"],
                ...
            })
```

* Stores results for the comparative summary.

#### Comparative Summary (Lines 274–310)

```python
    wins: dict[str, int] = {}
```

* Tracks how many times each method "wins" (has the lowest error).

```python
    for ds_file in dataset_files:
        ds_name = ds_file.name
        ds_rows = [r for r in summary_rows if r["dataset"] == ds_name]
```

* **List comprehension** `[r for r in summary_rows if r["dataset"] == ds_name]` — Creates a new list containing only the rows for the current dataset. It reads: "take each `r` from `summary_rows` if that row's dataset name matches."

```python
        for metric_key in ("MSE", "IAE", "ISE", "ITAE"):
            best_row = min(ds_rows, key=lambda r, mk=metric_key: r[mk])
```

* **`min(..., key=...)`** — Finds the element with the smallest value. The `key` function tells `min` *what* to compare.
* **`lambda r, mk=metric_key: r[mk]`** — An anonymous (unnamed) function. It takes a row `r` and returns the metric value `r[mk]`. The trick `mk=metric_key` captures the current value of `metric_key` into the lambda (this is a default-argument binding pattern in Python).

```python
            winner = str(best_row["method"])
            best[metric_key] = winner
            wins[winner] = wins.get(winner, 0) + 1
```

* **`wins.get(winner, 0)`** — Gets the current win count; if the method hasn't won before, returns `0`.
* **`+ 1`** — Increments the count.

#### Conclusion (Lines 317–392)

```python
    if wins:
        overall_best = max(wins, key=lambda m: wins[m])
```

* **`max(wins, key=...)`** — Finds the method name with the highest win count.
* `wins` is a dictionary; iterating over it yields the keys (method names). The `key` function maps each key to its value (win count).

```python
        overall_best_wins = wins[overall_best]
        total_decisions = sum(wins.values())
```

* **`sum(wins.values())`** — Sum of all win counts. `wins.values()` gives all the values in the dictionary.

```python
        all_methods = list(METHOD_DESCRIPTIONS.keys())
        worst = min(all_methods, key=lambda m: wins.get(m, 0))
```

* Finds the method with the *fewest* wins. Uses `wins.get(m, 0)` because a method might not appear in `wins` at all (zero wins).

```python
    metric_winners: dict[str, dict[str, int]] = {
        "MSE": {}, "IAE": {}, "ISE": {}, "ITAE": {}
    }
```

* **Nested dictionary:** For each metric, stores which methods won and how many times.

```python
    for row in best_per_dataset:
        for mk in ("MSE", "IAE", "ISE", "ITAE"):
            m = row.get(mk, "")
            if m:
                metric_winners[mk][m] = metric_winners[mk].get(m, 0) + 1
```

* Tallies wins per metric per method.

```python
    for mk, label in (
        ("MSE",  "Erro Médio Quadrático (MSE)"),
        ...
    ):
        if metric_winners[mk]:
            top = max(metric_winners[mk], key=lambda m: metric_winners[mk][m])
```

* For each metric, finds and reports which method won most often.

```python
    return "".join(lines)
```

* **`"".join(lines)`** — Concatenates all strings in the list into one single string with no separator. This is the complete Markdown report.

---

### 3.7 Main Routine (Lines 397–407)

```python
def main() -> None:
    report = generate_report()
    output_path = Path(__file__).parent / "relatorio_identificacao.md"
    output_path.write_text(report, encoding="utf-8")
    print(f"Relatório salvo em: {output_path}")
```

* **`generate_report()`** — Builds the entire report string.
* **`output_path.write_text(report, encoding="utf-8")`** — Writes the string to a file. `encoding="utf-8"` ensures special characters (like `τ`, `ä`, accented letters) are saved correctly.
* **`print(...)`** — Tells the user where the file was saved.

```python
if __name__ == "__main__":
    main()
```

* Same guard as in `identification.py` — only runs when executed directly.

---

## 4. Glossary

| Term | Meaning |
|------|---------|
| **Array (`np.ndarray`)** | A NumPy multi-dimensional container for numbers. Faster and more memory-efficient than Python lists. |
| **Context manager (`with`)** | A Python construct that automatically sets up and tears down resources (e.g., file handles, warning filters). |
| **Dead time (L)** | The delay between the input change and the first visible change in output. |
| **Docstring** | A string literal placed at the start of a module, class, or function to document it. |
| **Exception** | An error object raised (thrown) when something goes wrong. `try/except` handles it. |
| **f-string** | A Python string prefixed with `f` that embeds expressions in `{curly braces}`. |
| **FOPDT** | First-Order Plus Dead-Time — a simple dynamic model. |
| **IAE** | Integral of Absolute Error — area under |error| curve. |
| **ISE** | Integral of Squared Error — area under error² curve. |
| **ITAE** | Integral of Time × Absolute Error — penalizes late errors. |
| **Lambda** | An anonymous (unnamed) function defined inline: `lambda x: x + 1`. |
| **Linear interpolation** | Estimating a value between two known data points by assuming a straight line between them. |
| **List comprehension** | A compact syntax to create lists: `[expr for item in iterable if condition]`. |
| **MSE** | Mean Squared Error — average of error² across all samples. |
| **NaN** | *Not a Number* — a special floating-point value meaning "undefined." |
| **Savitzky-Golay filter** | A digital smoothing filter that fits local polynomials to overlapping subsets of data. |
| **Slicing** | Extracting parts of an array with `array[start:stop:step]`. |
| **Static gain (K)** | The ratio of output change to input change at steady state. |
| **Step response** | The system's output when the input suddenly changes from one constant level to another. |
| **Steady state** | The final, settled value the output reaches after all transients decay. |
| **Time constant (τ)** | The time needed for the output to reach about 63.2 % of its total change (in a pure first-order system). |
| **Trapezoidal rule** | A numerical integration method that approximates the area under a curve using trapezoids. |
| **Tuple** | An immutable ordered collection in Python: `(a, b, c)`. |
| **Type hint** | An annotation like `x: int` that documents a variable's expected type. |
