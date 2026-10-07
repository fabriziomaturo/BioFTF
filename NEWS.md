# BioFTF 3.0.0

This release brings the package back to CRAN. The previous CRAN release was
1.2-0.

## Changes that affect existing code

- All functional tools are now computed on Hill-number profiles of order `q`
  (default domain `0 <= q <= 2`) instead of the Patil-Taillie beta profile on
  `-1 <= beta <= 1`. The two profiles are linked by `q = beta + 1`, but they
  are measured on different scales, so numerical results differ from those of
  version 1.2-0.
- The functions `beta()`, `beta_plot()`, `betaprime()`, `betaprime_plot()`,
  `betasecond()`, `betasecond_plot()`, `curvature()`, `curvature_plot()`,
  `radius()`, `radius_plot()`, `arc()`, `area()` and `ranking()` are deprecated.
  They issue a warning and return the corresponding Hill-number tool.
- `alltools()` returns a list and no longer draws plots. `export_plots()`
  requires the output path in `file`.
- Functions no longer print intermediate objects or attach packages.

## New features

- `hill_profile()` evaluates the profile on any grid or set of orders,
  including `Inf`; negative orders are available on request.
- `hill_derivative()`, `hill_curvature()` and `hill_radius()` use analytical
  first and second derivatives.
- `hill_cum()` and `hill_area()` compute cumulative integral functions and
  areas by the trapezoidal rule, with an optional baseline.
- `hill_total_change()` reports changes in area between consecutive rows.
- `hill_surface()`, `hill_surface_plot()` and `hill_surface_compare()`
  compute, draw and compare biodiversity surfaces over time, either as
  rotatable `plotly` figures or as static perspective plots.
- Input is validated, and row names are kept in all outputs.
- Added tests and a vignette.
