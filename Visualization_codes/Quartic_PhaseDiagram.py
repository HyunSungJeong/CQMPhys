import matplotlib.pyplot as plt
import numpy as np

ymax = 0.4
xmax = 0.8

# Phase boundary data (in increasing order of lambda_z)


# Define the boundary as two arrays (clockwise order)
x = np.array([0.735, 0.685/2, 0.565/2, 0.525/2, 0.485/2, 0.445/2, 0.405/2, 0.375/2, 0.335/2, 0, 0, 1])
y = np.array([0.001, 0.01, 0.1, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, ymax, 0, 0])

# Plot and fill the region
plt.fill(x, y, color='green')  # Fill the region with red
plt.xlim(0, xmax)
plt.ylim(0, ymax)
plt.gca().set_aspect('equal')
plt.show()