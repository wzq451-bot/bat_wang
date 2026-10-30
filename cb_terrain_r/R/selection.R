# Conditional logistic models, fixed analytical choice sets, sensitivity fits and predictions.
# See README.md for numerical conventions and scope.

c = read("choices.csv.gz")

c = c[order(c$step_id, c$choice_index), ]

stopifnot(nrow(c) == 137760, all(table(c$step_id) == 41), all(tapply(c$used, c$step_id, sum) == 
    1), all(c$used[seq(1, nrow(c), 41)] == 1))

m = read("steps.csv")

m = m[order(m$step_id), ]

sc = jsonlite::fromJSON("data/selection_scaling.json")$predictors

for (k in c("tpi", "elev", "tri", "start_agl")) c[[paste0(k, "_z")]] = (c[[k]] - sc[[k]]$mean)/sc[[k]]$sd

n = nrow(m)

X = array(0, c(n, 41, 4))

for (j in 1:3) X[, , j] = matrix(c[[c("tpi_z", "elev_z", "tri_z")[j]]], n, 41, byrow = TRUE)

az = matrix(c$start_agl_z, n, 41, byrow = TRUE)

X[, , 4] = az * X[, , 1]

terms = c("tpi_z", "elev_z", "tri_z", "agl_z:tpi_z")

rows = list()

fits = list()

one = function(name, xx = X, tt = terms, mask = rep(TRUE, n)) {
    r = clogit(xx[mask, , , drop = FALSE], m$ind[mask])
    fits[[name]] <<- r
    for (j in seq_along(tt)) rows[[length(rows) + 1]] <<- cbind(model = name, term = tt[j], cbind(stats(r, 
        j), se = sqrt(r$V_model[j, j]), p_normal = 2 * pnorm(-abs(r$beta[j]/sqrt(r$V_model[j, j]))), 
        ci_lo = r$beta[j] - qnorm(0.975) * sqrt(r$V_model[j, j]), ci_hi = r$beta[j] + qnorm(0.975) * 
            sqrt(r$V_model[j, j])))
    cat("Selection", name, "\n")
}

one("M4_main_effects", X[, , 1:3, drop = FALSE], terms[1:3])

one("M5_interaction")

one("primary")

for (b in sort(unique(m$ind))) one(paste0("LOO_without_HQ", b), mask = m$ind != b)

wc = m$start_agl - ave(m$start_agl, m$ind)

xx = X

xx[, , 4] = X[, , 1] * (wc/sd(rep(wc, each = 41)))

one("within_bat_AGL", xx)

one("nonnegative_start_AGL", mask = m$start_agl >= 0)

one("valid_heights_both", mask = truth(m$start_altitude_valid) & truth(m$end_altitude_valid))

one("dt_le_2p5", mask = m$dt_min <= 2.5)

one("dt_gt_2p5", mask = m$dt_min > 2.5)

h = as.POSIXlt(m$t, tz = "UTC")

hour = h$hour + h$min/60 + h$sec/3600

hz = as.vector(scale((hour - 18)%%24))

h2 = hz^2 - mean(hz^2)

xx = array(0, c(n, 41, 10))

xx[, , 1:4] = X

for (j in 1:3) {
    xx[, , j + 4] = X[, , j] * hz
    xx[, , j + 7] = X[, , j] * h2
}

one("clock_time_terrain_interactions", xx, c(terms, paste0("clock_z:", terms[1:3]), paste0("clock_z2:", 
    terms[1:3])))

one("before_midnight", mask = hour >= 18)

one("after_midnight", mask = hour < 12)

ids = sort(unique(m$ind))

xx = array(0, c(n, 41, 22))

tt = character()

for (i in seq_along(ids)) for (j in 1:3) {
    xx[, , (i - 1) * 3 + j] = X[, , j] * (m$ind == ids[i])
    tt = c(tt, paste0("HQ", ids[i], ":", terms[j]))
}

xx[, , 22] = X[, , 1] * wc/sc$start_agl$sd

one("individual_terrain_slopes", xx, c(tt, "within_agl_z:tpi_z"))

for (b in ids) one(paste0("individual_HQ", b), mask = m$ind == b)

a = split(c[c$used == 0, ], c$step_id[c$used == 0])

sup = t(sapply(a, function(z) c(sd(z$tpi), diff(range(z$elev)))))

sup = sup[match(as.character(m$step_id), rownames(sup)), ]

xx = array(0, c(n, 41, 6))

xx[, , 1:4] = X

for (j in 1:2) xx[, , j + 4] = X[, , 1] * as.vector(scale(sup[, j]))

one("candidate_terrain_variation", xx, c(terms, "tpi_z:available_tpi_sd", "tpi_z:available_elev_range"))

write(do.call(rbind, rows), "selection_coefficients.csv")

r = fits[["M5_interaction"]]

h = seq(0, 300, 5)

L = cbind(1, 0, 0, (h - sc$start_agl$mean)/sc$start_agl$sd)

p = as.vector(L %*% r$beta)

se = sqrt(rowSums((L %*% r$V) * L))

write(data.frame(agl_m = h, cond_slope = p, cluster_se = se, cluster_lo = p - qt(0.975, 6) * se, 
    cluster_hi = p + qt(0.975, 6) * se, relative_weight = exp(-p)), "selection_predictions.csv")

write(do.call(rbind, lapply(names(fits), function(nm) {
    z = fits[[nm]]
    data.frame(model = nm, n = z$n, n_bats = z$G, loglik = z$loglik, max_mean_score = z$gradient)
})), "selection_fit_diagnostics.csv")

saveRDS(fits, "results/selection_fits.rds")

