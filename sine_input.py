import control as ct
import matplotlib.pyplot as plt
import numpy as np

# 1. Define the Laplace variable 's'
s = ct.tf('s')

# Define the Transfer Functions
g1 = 20 / ((s + 2) * (s + 10))
g2 = 2 / (s + 2)
g3 = (2*15*200) / ((s + 2) * (s+15) * (s + 200))

# 2. Setup Time and Input Signal (Sine Wave)
t = np.linspace(0, 10, 1000)  # Time from 0 to 10 seconds
frequency = 1.0  # Frequency in Hz
u = np.sin(2 * np.pi * frequency * t)  # Sine wave input

# 3. Simulation and Plotting
plt.figure(figsize=(12, 8))

# Plot the Input for reference
plt.plot(t, u, 'k--', label='Input: Sine Wave (1Hz)', alpha=0.5)

systems = [g1, g2, g3]
labels = [
    r'$g_1(s) = \frac{20}{(s+2)(s+10)}$',
    r'$g_2(s) = \frac{2}{s+2}$',
    r'$g_3(s) = \frac{40}{(s+2)(s^2+10s+20)}$'
]

for sys, label in zip(systems, labels):
    # Forced response calculates output for the specific input 'u'
    # Source: https://python-control.readthedocs.io/en/0.10.2/generated/control.forced_response.html
    time, response = ct.forced_response(sys, T=t, U=u)
    plt.plot(time, response, label=f'Output: {label}')

plt.title('System Response to Sine Wave Input')
plt.xlabel('Time (seconds)')
plt.ylabel('Amplitude')
plt.grid(True)
plt.legend(loc='upper right')
plt.show()