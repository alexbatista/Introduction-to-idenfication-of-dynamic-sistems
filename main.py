import control as ct
import matplotlib.pyplot as plt
import numpy as np

# Define the Laplace variable 's' for easier algebraic creation
s = ct.tf('s')

# 1. g1(s) = 20 / ((s+2)(s+10))
# Expanded: 20 / (s^2 + 12s + 20)
g1 = 20 / ((s + 2) * (s + 10))

# 2. g2(s) = 2 / (s+2)
g2 = 2 / (s + 2)

# 3. g3(s) = 40 / ((s+2)(s^2 + 10s + 20))
# Expanded: 40 / (s^3 + 12s^2 + 40s + 40)
g3 = 400 / ((s + 2) * (s**2 + 100 * s + 200))

# List of systems and labels for plotting
systems = [g1, g2, g3]
labels = [
    r'$g_1(s) = \frac{20}{(s+2)(s+10)}$',
    r'$g_2(s) = \frac{2}{s+2}$',
    r'$g_3(s) = \frac{400}{(s+2)(s^2+100s+200)}$'
]

# Plotting the Step Responses
plt.figure(figsize=(10, 6))

T = np.linspace(0, 3, 1000)

for sys, label in zip(systems, labels):
    time, response = ct.step_response(sys, T)
    plt.plot(time, response, label=label)

plt.title('Step Response of Linear Systems')
plt.xlabel('Time (seconds)')
plt.ylabel('Amplitude')
plt.xlim(left=0)
plt.grid(True)
plt.legend()
plt.show()