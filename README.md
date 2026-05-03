# Lipschitz Gradient Descent on Hakobyan–Li Slit Domains

Numerical investigation of Lipschitz 1-connectedness in a Hakobyan–Li slit
domain $X_N$ with defining sequence $r_n = 1/n$, conducted as part of an
undergraduate research project at the University of Hawaiʻi at Mānoa under
Prof. Matthew Romney.

## What this does

Discretizes the slit domain $X_N$ as a graph, constructs an initial filling
map $\bar{c}_N : Q \to X_N$ from the boundary loop, and runs gradient
descent to minimize the Lipschitz constant of the filling. Reports the
converged constant for $N = 2, 3, 4$ at fixed grid resolution.

## Quick start

Open MATLAB (R2021a or later) in the repo root and run:

```matlab
driver
```

This builds the domain graphs, computes the lane filling, runs the descent,
and saves results to `driver_results.mat`.

## Files

**Driver and pipeline**
- `driver.m` — main entry point; runs the full experiment
- `test_lane_filling.m` — sanity-check script for the initial filling

**Algorithm**
- `gd_lane.m`, `gd_top.m`, `gd_loop.m` — coordinate descent variants
- `compute_lip_lane.m` — Lipschitz constant over graph + cross-slit pairs
- `global_lip.m` — global Lipschitz constant computation
- `gd_dlookup.m`, `gd_get_row.m` — distance-cache helpers

**Domain construction**
- `hl_build_domain_graph.m`, `hl_build_graph_from_points.m` — graph builders
- `hl_lane_filling.m`, `hl_top_filling.m` — initial filling constructions
- `hl_draw_slits_iter.m`, `hl_overlay_grid_with_slit_holes.m`,
  `hl_plot_graph_edges.m` — visualization helpers

**Visualization**
- `visualize_veritcal_filling.m` — animation of the descent
- `capture_frame.m`, `draw_colored_arrows.m` — frame export and quiver plots

## Results

Select the .mp4 for the descent videos.

## Poster

A poster summarizing this work was presented at the Math 480 Spring 2026
seminar, viewable here: https://www.overleaf.com/read/cbvvpfvmpnst#9507d4