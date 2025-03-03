#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

x_min = 6.5e-3       # 6.5 mm
x_max = 14e-3        # 14 mm
y_min_top = 9.5e-3   # Flipped: 14 - 4.5 = 9.5 mm
y_max_top = 10.25e-3 # Flipped: 14 - 3.75 = 10.25 mm
y_min_bottom = 7.75e-3 # Flipped: 14 - 6.25 = 7.75 mm
y_max_bottom = 8.5e-3  # Flipped: 14 - 5.5 = 8.5 mm

xspace = np.linspace(0, 14e-3, 140)  # 0 to 15 mm, 151 points
yspace = np.linspace(0, 14e-3, 140)    # 0 to 15 mm, 151 points

sigma_varied = np.zeros((len(yspace), len(xspace)), dtype=float)
eps_varied = np.ones((len(yspace), len(xspace)), dtype=float)

# Create masks for the waveguide regions
x_mask = (xspace >= x_min) & (xspace <= x_max)
y_top_mask = (yspace >= y_min_top) & (yspace <= y_max_top)
y_bottom_mask = (yspace >= y_min_bottom) & (yspace <= y_max_bottom)

# Create meshgrid for proper indexing
X, Y = np.meshgrid(xspace, yspace)

# Apply conductivity to the waveguide regions
mask = ((X >= x_min) & (X <= x_max)) & (
    ((Y >= y_min_top) & (Y <= y_max_top)) |
    ((Y >= y_min_bottom) & (Y <= y_max_bottom))
)
sigma_varied[mask] = 1e100

eps0 = 8.854e-12
eps_varied *= eps0

plt.figure(figsize=(6, 4))
plt.imshow(
    sigma_varied,
    extent=[xspace[0], xspace[-1], yspace[0], yspace[-1]],
    origin='upper',
    aspect='auto'
)

plt.colorbar(label="σ (S/m)")
plt.title("Waveguide Design")
plt.xlabel("x (m)")
plt.ylabel("y (m)")

plt.tight_layout()
plt.show()
