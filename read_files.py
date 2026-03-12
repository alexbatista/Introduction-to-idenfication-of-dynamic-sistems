import control as ct
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


def load_and_resample_dataset(file_path):
    # File format: amplitude,time
    u_data, t_data = np.loadtxt(file_path, delimiter=',', unpack=True)

    # forced_response requires exactly equally spaced timestamps.
    # Build a uniform time base and interpolate the measured input onto it.
    sort_idx = np.argsort(t_data)
    t_data = t_data[sort_idx]
    u_data = u_data[sort_idx]
    t_uniform = np.linspace(t_data[0], t_data[-1], len(t_data))
    u_uniform = np.interp(t_uniform, t_data, u_data)
    return t_data, u_data, t_uniform, u_uniform

# 2. Define the Transfer Functions
s = ct.tf('s')

# G1(s) = 20 / ((s+2)(s+10))
g1 = 20 / ((s + 2) * (s + 10))

# G2(s) = 2 / (s+2)
g2 = 2 / (s + 2)

# G3(s) = 40 / ((s+2)(s^2 + 10s + 20))
g3 = 40 / ((s + 2) * (s**2 + 10 * s + 20))

systems = [g1, g2, g3]
labels = [
    r'$g_1(s) = \frac{20}{(s+2)(s+10)}$',
    r'$g_2(s) = \frac{2}{s+2}$',
    r'$g_3(s) = \frac{40}{(s+2)(s^2+10s+20)}$'
]

# 3. Simulate every dataset file in separate figures
dataset_files = sorted(Path('dataset').glob('*.txt'))

if not dataset_files:
    print("No .txt files found in dataset/")
else:
    for dataset_file in dataset_files:
        try:
            t_data, u_data, t_uniform, u_uniform = load_and_resample_dataset(dataset_file)
        except Exception as e:
            print(f"Error loading {dataset_file}: {e}")
            continue

        plt.figure(figsize=(12, 8))
        plt.plot(t_data, u_data, label='Input (Measured)', color='black', linestyle=':', alpha=0.5)
        plt.plot(t_uniform, u_uniform, label='Input (Resampled)', color='gray', alpha=0.8)
        plt.xlabel('Time (Tempo)')
        plt.ylabel('Amplitude')
        plt.title(f'Measured Input Data ({dataset_file.name})')
        
        # for sys, label in zip(systems, labels):
        #     # Source: https://python-control.readthedocs.io/en/0.10.2/generated/control.forced_response.html
        #     time, response = ct.forced_response(sys, T=t_uniform, U=u_uniform)
        #     plt.plot(time, response, label=f'Output: {label}', linewidth=2)

        # plt.title(f'System Response to Measured Input Data ({dataset_file.name})')
        # plt.xlabel('Time (Tempo)')
        # plt.ylabel('Amplitude')
        # plt.legend(loc='lower right')
        # plt.grid(True, linestyle='--', alpha=0.7)
        # plt.tight_layout()

    plt.show()
