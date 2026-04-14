import numpy as np

# ─────────────────────────────────────────
# USER INPUTS
# ─────────────────────────────────────────
half_side   = 1.0        # half length of square side
R_outer     = 60.0       # outer circular boundary radius
cx, cy      = 0.0, 0.0   # center
N_points    = 350        # total boundary points (multiple of 4 recommended)
output_file = "cylinder_points.txt"
# ─────────────────────────────────────────


# -------------------------------------------------
# Generate Square Boundary (Angle Based)
# -------------------------------------------------
theta = np.linspace(0, 2*np.pi, N_points, endpoint=False)

x_inner = []
y_inner = []

for t in theta:
    dx = np.cos(t)
    dy = np.sin(t)

    # Scale ray to hit square boundary
    scale = half_side / max(abs(dx), abs(dy))

    x_inner.append(cx + scale * dx)
    y_inner.append(cy + scale * dy)

x_inner = np.array(x_inner)
y_inner = np.array(y_inner)


# -------------------------------------------------
# Generate Outer Circular Boundary
# -------------------------------------------------
x_outer = cx + R_outer * np.cos(theta)
y_outer = cy + R_outer * np.sin(theta)


# -------------------------------------------------
# Write File
# Format:
# N N
# x_inner y_inner x_outer y_outer
# -------------------------------------------------
with open(output_file, "w") as f:

    f.write(f"{N_points} {N_points}\n")

    for xi, yi, xo, yo in zip(x_inner, y_inner, x_outer, y_outer):
        f.write(f"{xi:.8f} {yi:.8f} {xo:.8f} {yo:.8f}\n")


print(f"\nDone! File '{output_file}' created.")
print(f"Total boundary points: {N_points}")
print("Each row: Square boundary point + corresponding circular boundary point")
