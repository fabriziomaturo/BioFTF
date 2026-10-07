## ---------------------------------------------------------------------------
## 3. Data summaries
## ---------------------------------------------------------------------------

datirel <- function(x) {
  .relative_abundance(x)
}

summary_species <- function(x) {
  x <- .check_abundance_matrix(x)
  data.frame(
    richness = rowSums(x > 0),
    total_abundance = rowSums(x),
    min_abundance = apply(x, 1L, min),
    max_abundance = apply(x, 1L, max),
    row.names = rownames(x)
  )
}

summary_species_relative <- function(x) {
  p <- .relative_abundance(x)
  data.frame(
    richness = rowSums(p > 0),
    min_relative_abundance = apply(p, 1L, min),
    max_relative_abundance = apply(p, 1L, max),
    row.names = rownames(p)
  )
}

