# Python/R comparison. Single-cluster Python zero SEs are explicitly checked and corrected to NA.
# See README.md for numerical conventions and scope.

reports = list()

check = function(actual, ref, keys, cols, tol = 1e-06) {
    a = read(actual, "results")
    b = read(ref, "reference")
    aa = a[keys]
    bb = b[keys]
    for (k in keys) if (is.numeric(aa[[k]])) {
        aa[[k]] = sprintf("%.8f", aa[[k]])
        bb[[k]] = sprintf("%.8f", bb[[k]])
    }
    ka = do.call(paste, c(aa, sep = "|"))
    kb = do.call(paste, c(bb, sep = "|"))
    stopifnot(!anyDuplicated(ka), !anyDuplicated(kb), all(kb %in% ka))
    a = a[match(kb, ka), ]
    for (col in cols) {
        x = a[[col]]
        y = b[[col]]
        corrected = 0L
        if (col == "cluster_se" && "n_clusters" %in% names(b)) {
            single = b$n_clusters == 1
            stopifnot(all(is.na(x[single])), all(y[single] == 0 | is.na(y[single])))
            corrected = sum(single & !is.na(y))
            y[single] = NA_real_
        }
        ok = identical(is.na(x), is.na(y))
        delta = max(c(0, abs(x - y)), na.rm = TRUE)
        pass = ok && delta <= tol
        reports[[length(reports) + 1]] <<- data.frame(output = actual, column = col, n = length(y), 
            max_absolute_difference = delta, tolerance = tol, documented_NA_corrections = corrected, 
            passed = pass)
    }
}

check("vertical_refit.csv", "vertical_refit.csv", "model", c("beta", "cluster_se", "cluster_lo", 
    "cluster_hi", "cluster_p", "n", "n_bats"))

check("vertical_segments.csv", "vertical_segments.csv", "seg", c("n", "beta", "ground_sd"))

check("vertical_null.csv", "vertical_null.csv", "replicate", c("mean_SD_AGL", "within_beta"))

check("movement_checks.csv", "movement_checks.csv", c("outcome", "model"), c("beta", "cluster_se", 
    "cluster_lo", "cluster_hi", "cluster_p", "n", "n_bats", "wild_p"))

check("individual_slopes.csv", "individual_slopes.csv", c("outcome", "ind"), c("n", "beta"))

check("shape_predictions.csv", "shape_predictions.csv", c("outcome", "model", "TRI"), c("fit", "lo", 
    "hi", "response_mean", "response_sd"))

check("shape_bins.csv", "shape_bins.csv", c("outcome", "bin"), c("TRI", "adjusted", "n"))

check("activity_interactions.csv", "activity_interactions.csv", c("model", "term"), c("beta", "cluster_se", 
    "cluster_lo", "cluster_hi", "cluster_p"))

a = read("selection_coefficients.csv", "results")

names(a)[match(c("cluster_lo", "cluster_hi", "cluster_p", "n", "n_bats"), names(a))] = c("cluster_ci_lo", 
    "cluster_ci_hi", "cluster_p_t", "n_steps", "n_clusters")

write(a, "selection_comparison.csv")

for (ref in c("selection_existing_checks.csv", "selection_new_checks.csv")) check("selection_comparison.csv", 
    ref, c("model", "term"), c("beta", "cluster_se", "cluster_ci_lo", "cluster_ci_hi", "cluster_p_t", 
        "n_steps", "n_clusters", "se", "p_normal", "ci_lo", "ci_hi"))

sp = read("selection_predictions.csv", "results")

names(sp)[match(c("cluster_lo", "cluster_hi"), names(sp))] = c("cluster_ci_lo", "cluster_ci_hi")

write(sp, "selection_prediction_comparison.csv")

check("selection_prediction_comparison.csv", "selection_predictions.csv", "agl_m", c("cond_slope", 
    "cluster_se", "cluster_ci_lo", "cluster_ci_hi"))

check("synchronized_window_audit.csv", "synchronized_window_audit.csv", "row", c("duration", "count", 
    "mean_TRI", "coverage"))

a = read("synchronized_window_audit.csv", "results")

b = read("synchronized_window_audit.csv", "reference")

stopifnot(identical(truth(a$valid), truth(b$valid[match(a$row, b$row)])))

check("candidate_support.csv", "candidate_support.csv", "step_id", c("tpi_sd", "elev_min", "elev_max", 
    "tpi_min", "tpi_max", "elev_range"))

check("Table_S1.csv", "Table_S1.csv", "ind", setdiff(names(read("Table_S1.csv", "reference")), "ind"))

report = do.call(rbind, reports)

write.csv(report, "validation/python_R_comparison.csv", row.names = FALSE)

print(report)

stopifnot(all(report$passed))

