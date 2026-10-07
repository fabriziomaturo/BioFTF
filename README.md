# BioFTF

BioFTF provides functional tools for biodiversity assessment based on
Hill-number diversity profiles. The profile of each community is treated as a
function of the order `q` and summarised by its derivatives, curvature, radius
of curvature, arc length, cumulative integral function and area. When rows are
ordered in time, the profiles form a biodiversity surface that can be drawn as
a static or rotatable figure and compared with a second surface.

## Installation

```r
install.packages("BioFTF")
```

## Example

```r
library(BioFTF)

x <- matrix(c(54, 20, 17,  6,  3,
              25, 25, 25, 25,  0,
              35, 35, 30,  0,  0),
            nrow = 3, byrow = TRUE)
rownames(x) <- c("A", "B", "C")

hill_profile(x, domain = c(0, 1, 2, Inf))
hill_profile_plot(x, from = 0, to = 10, n = 501)
hill_first_derivative(x, domain = c(0, 1, 2))
hill_area(x)
hill_ranking(x)
```

`D(0)` is species richness, `D(1)` the exponential of Shannon entropy and
`D(2)` the inverse Simpson concentration. The default domain is `0 <= q <= 2`;
`from`, `to`, `n` and `domain` change it. Negative orders are rejected unless
`allow_negative = TRUE`, because for `q < 0` the Hill number exceeds species
richness and is driven by the rarest species.

## Surfaces

```r
years <- 2018:2023
z <- matrix(c(40, 30, 15, 10,  5,
              35, 30, 18, 12,  5,
              30, 28, 20, 15,  7,
              25, 25, 22, 18, 10,
              20, 24, 24, 20, 12,
              15, 22, 25, 22, 16),
            nrow = length(years), byrow = TRUE)
rownames(z) <- years

hill_total_change(z)
hill_surface_plot(z, from = 0, to = 10, n = 121, plot = "static")
hill_surface_plot(z, from = 0, to = 10, n = 121)   # rotatable, needs plotly
```

## Versions 1.x

Versions 1.x computed the same tools on the Patil-Taillie beta profile over
`-1 <= beta <= 1`, which corresponds to `0 <= q <= 2` through `q = beta + 1`.
The old function names (`beta()`, `area()`, `arc()`, ...) are deprecated and
now return the Hill-number tools; see `?"BioFTF-deprecated"`.

## References

Di Battista, T., Fortuna, F. and Maturo, F. (2017). BioFTF: An R package for
biodiversity assessment with the functional data analysis approach.
*Ecological Indicators*, 73, 726-732. <https://doi.org/10.1016/j.ecolind.2016.10.032>

Maturo, F. and Di Battista, T. (2018). A functional approach to Hill's numbers
for assessing changes in species variety of ecological communities over time.
*Ecological Indicators*, 84, 70-81. <https://doi.org/10.1016/j.ecolind.2017.08.016>

Maturo, F. (2018). Unsupervised classification of ecological communities ranked
according to their biodiversity patterns via a functional principal component
decomposition of Hill's numbers integral functions. *Ecological Indicators*,
90, 305-315. <https://doi.org/10.1016/j.ecolind.2018.03.013>
