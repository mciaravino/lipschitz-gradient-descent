# Math480Project

MATLAB scripts for constructing and analyzing a Hakobyan–Li slit domain.

## Usage

Run `drivingScript.m` to execute the full pipeline. All other files are helper functions called automatically.

Set `N` and `MaxN` at the top of `drivingScript.m` to control the dyadic depth and grid resolution respectively.

## Files

| File | Description |
|------|-------------|
| `drivingScript.m` | Main script — run this |
| `hl_draw_slits_iter.m` | Generates slit geometry for dyadic stage N |
| `hl_overlay_grid_with_slit_holes.m` | Builds grid with slit nodes removed |
| `hl_build_graph_from_points.m` | Builds graph from grid points |
| `hl_plot_graph_edges.m` | Plots graph edges over the domain |
| `gd_get_row.m` | Distance cache helper |
| `gd_dlookup.m` | Distance lookup helper |

## Requirements

MATLAB with the Graph and Network Algorithms toolbox.
