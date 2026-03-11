# Lipschitz Constant Estimation via Gradient Descent
A MATLAB implementation that uses gradient descent to estimate the Lipschitz 
constant of a function, with graph construction and geometric slit visualization.

## Files
- `drivingScript.m` — Main script to run the project
- `gd_compute_lip.m` — Computes global Lipschitz constant of the filling map over all graph edges
- `gd_local_lip.m` — Computes worst Lipschitz ratio over edges incident to a single node (used inside descent loop)
- `gd_dlookup.m` / `gd_get_row.m` — Distance cache utilities
- `hl_build_graph_from_points.m` — Constructs graph from point data
- `hl_draw_slits_iter.m` — Iterative slit drawing
- `hl_overlay_grid_with_slit_holes.m` — Grid and slit overlay visualization
- `hl_plot_graph_edges.m` — Graph edge plotting

## Presentation
Slides in progress — [view draft](./Lipschitz_Gradient_Descent_Project.pdf)

## Requirements
MATLAB (R2021a or later recommended)
