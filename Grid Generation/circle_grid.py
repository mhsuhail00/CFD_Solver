import numpy as np

# ─────────────────────────────────────────
#  USER INPUTS
# ─────────────────────────────────────────
R_inner     = 1.0        # cylinder (inner) radius
R_outer     = 60.0        # outer domain radius
cx, cy      = 0.0, 0.0   # center of cylinder
N_points    = 450         # number of points on each circle
output_file = "cylinder_points.txt"
# ─────────────────────────────────────────

theta = np.linspace(0, 2 * np.pi, N_points, endpoint=False)

# --- Generate inner and outer circle points at same angles ---
inner_pts = [(cx + R_inner * np.cos(t), cy + R_inner * np.sin(t)) for t in theta]
outer_pts = [(cx + R_outer * np.cos(t), cy + R_outer * np.sin(t)) for t in theta]

# --- Write to file: 1 inner + 1 outer point per row ---
# Format: x_inner, y_inner, x_outer, y_outer
with open(output_file, "w") as f:

    #f.write(f"# 2D Circular Cylinder Points\n")
    #f.write(f"# Inner Radius: {R_inner} | Outer Radius: {R_outer} | Center: ({cx},{cy})\n")
    #f.write(f"# N_points: {N_points}\n")
    f.write(f"{N_points} {N_points}\n")
    # Format: x_inner, y_inner, x_outer, y_outer

    for (xi, yi), (xo, yo) in zip(inner_pts, outer_pts):
        f.write(f"{xi:.8f} {yi:.8f} {xo:.8f} {yo:.8f}\n")

print(f"Done! {N_points} rows written to '{output_file}'")
print(f"Each row: 1 inner circle point + 1 outer domain point (same angle)")
