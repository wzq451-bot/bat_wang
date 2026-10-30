# Recompute descriptive counts, candidate support and coefficient source tables.
# See README.md for numerical conventions and scope.

f <- read("fixes.csv")

v <- read("vertical_analysis.csv", "results")

m <- read("steps.csv")

d <- read("activity_rows.csv.gz")

raw <- read("raw_records.csv.gz")

s1 <- do.call(rbind, lapply(sort(unique(m$ind)), function(id) {
    dt <- d$raw_dt_min[d$ind == id] * 60
    data.frame(ind = paste0("HQ", id), filtered_fixes = sum(f$ind == id), vertical_fixes = sum(v$ind == 
        id), selection_steps = sum(m$ind == id), movement_records = sum(d$ind == id), raw_records = sum(raw$ind == 
        id), raw_dt_min_seconds = min(dt), raw_dt_median_seconds = median(dt), raw_dt_max_seconds = max(dt))
}))

write(s1, "Table_S1.csv")

s2 <- read("selection_coefficients.csv", "results")

write(s2[s2$model != "primary", ], "Table_S2.csv")

write(read("movement_checks.csv", "results"), "Table_S3_movement.csv")

write(read("activity_interactions.csv", "results"), "Table_S3_activity_terms.csv")

candidates <- read("choices.csv.gz")

candidates <- candidates[candidates$used == 0, ]

support <- do.call(rbind, lapply(split(candidates, candidates$step_id), function(z) {
    data.frame(step_id = z$step_id[1], tpi_sd = sd(z$tpi), elev_min = min(z$elev), elev_max = max(z$elev), 
        tpi_min = min(z$tpi), tpi_max = max(z$tpi), elev_range = diff(range(z$elev)))
}))

write(support, "candidate_support.csv")

