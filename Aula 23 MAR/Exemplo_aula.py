import matplotlib.pyplot as plt
import numpy as np


def simulate_arx(t: np.ndarray, u: np.ndarray, e: np.ndarray) -> np.ndarray:
    n = len(t)
    y = np.zeros(n)

    for k in range(1, n):
        y[k] = 0.85 * y[k - 1] + 0.15 * u[k - 1]  # + e[k]

    return y


def plot_arx(t: np.ndarray, u: np.ndarray, e: np.ndarray) -> None:

    y = simulate_arx(t, u, e)

    fig, ax = plt.subplots(figsize=(10, 4))

    ax.plot(t, y, label="y[t]", color="tab:green")
    ax.set_ylabel("Output y[t]")
    ax.set_xlabel("Time")
    ax.legend()
    ax.grid(True)

    fig.suptitle("ARX Model: y[t] = 0.85·y[t-1] + 0.15·u[t-1] + e[t]")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    N = 30
    t = np.arange(N)
    u = np.ones(N)
    # u[:50] = 0
    e = np.random.normal(0, 0.1, N)

    plot_arx(t, u, e)
