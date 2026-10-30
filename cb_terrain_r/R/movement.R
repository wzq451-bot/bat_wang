# Linear models, exact sign-enumeration tests, fixed-knot curves and synchronized activity.
# See README.md for numerical conventions and scope.

d = prepare_movement(read("activity_rows.csv.gz"))

d$log_rate = log(d[["运动量"]]/d$raw_dt_min)

rows = list()

inds = list()

pred = list()

bins = list()

sc = jsonlite::fromJSON("data/movement_scaling.json")

kn = jsonlite::fromJSON("data/shape_knots.json")$knots_z

one = function(kind, name, b, cols = NULL, equal = FALSE, segment = FALSE, wild = FALSE) {
    y = if (kind == "straightness") 
        "straight_recomputed"
    else "log_rate"
    if (is.null(cols)) 
        cols = c("end_tri_z", if (kind == "activity") "straight_recomputed_z", "log_dt_z", "log_span_z")
    X = design(b, cols, if (segment) 
        "seg"
    else "ind")
    r = ols(X, b[[y]], b$ind, if (equal) 
        1/ave(b$ind, b$ind, FUN = length)
    else rep(1, nrow(b)))
    rr = cbind(outcome = kind, model = name, stats(r))
    rr$wild_p = NA_real_
    if (wild) {
        r0 = ols(X[, -2, drop = FALSE], b[[y]], b$ind)
        ids = sort(unique(b$ind))
        signs = as.matrix(expand.grid(rep(list(c(-1, 1)), length(ids))))
        ts = apply(signs, 1, function(s) {
            rb = ols(X, r0$fit + r0$e * s[match(b$ind, ids)], b$ind)
            abs(rb$beta[2]/sqrt(rb$V[2, 2]))
        })
        rr$wild_p = mean(ts >= abs(r$beta[2]/sqrt(r$V[2, 2])) - 1e-12)
    }
    rows[[length(rows) + 1]] <<- rr
    if (name == "primary") 
        for (id in sort(unique(b$ind))) {
            a = b[b$ind == id, ]
            xx = design(a, cols)
            bb = qr.solve(xx, a[[y]])
            inds[[length(inds) + 1]] <<- data.frame(outcome = kind, ind = id, n = nrow(a), beta = bb[2])
        }
    r
}

basis = function(x) {
    dd = function(j) (pmax(x - kn[j], 0)^3 - pmax(x - kn[4], 0)^3)/(kn[4] - kn[j])
    cbind(dd(1) - dd(3), dd(2) - dd(3))
}

for (kind in c("straightness", "activity")) {
    one(kind, "primary", d, wild = TRUE)
    one(kind, "equal_bat_weight", d, equal = TRUE)
    one(kind, "segment_intercepts", d, segment = TRUE)
    one(kind, "nonoverlapping", nonoverlap(d))
    one(kind, "raw_duration_45_75_seconds", d[d$raw_dt_min >= 0.75 & d$raw_dt_min <= 1.25, ])
    cs = c("end_tri_z", if (kind == "activity") "straight_recomputed_z", "log_dt_z", "log_span_z")
    one(kind, "clock_time", d, c(cs, "clock_sin", "clock_cos"))
    d$log_hdop = log(d$HDOP)
    one(kind, "GPS_quality_covariates", d, c(cs, "log_hdop", "卫星"))
    one(kind, "strict_window_quality", d[d$window_max_hdop <= 5 & d$window_min_sat >= 5, ])
    for (cut in c(6, 10, 20, 50)) one(kind, paste0("count_ge_", cut), d[d[["运动量"]] >= cut, 
        ])
    for (cut in c(5, 10, 20)) one(kind, paste0("rate_ge_", cut), d[d[["运动量"]]/d$raw_dt_min >= 
        cut, ])
    for (id in sort(unique(d$ind))) one(kind, paste0("without_HQ", id), d[d$ind != id, ])
    bs = basis(d$end_tri_z)
    d$spline1 = bs[, 1]
    d$spline2 = bs[, 2]
    cs = c("end_tri_z", "spline1", "spline2", if (kind == "activity") "straight_recomputed_z", "log_dt_z", 
        "log_span_z")
    y = if (kind == "straightness") 
        "straight_recomputed"
    else "log_rate"
    grid = seq(quantile(d$end_tri, 0.025), quantile(d$end_tri, 0.975), length.out = 180)
    z = (grid - sc$end_tri$mean)/sc$end_tri$sd
    bb = basis(z)
    for (name in c("primary", "equal_bat_weight", paste0("without_HQ", sort(unique(d$ind))))) {
        b = if (startsWith(name, "without_HQ")) 
            d[d$ind != as.numeric(sub("without_HQ", "", name)), ]
        else d
        X = design(b, cs)
        r = ols(X, b[[y]], b$ind, if (name == "equal_bat_weight") 
            1/ave(b$ind, b$ind, FUN = length)
        else rep(1, nrow(b)))
        L = matrix(0, length(grid), ncol(X))
        L[, 1] = 1
        L[, 2] = z
        L[, 3:4] = bb
        L[, (length(cs) + 2):ncol(L)] = 1/length(unique(b$ind))
        p = as.vector(L %*% r$beta)
        se = sqrt(pmax(0, rowSums((L %*% r$V) * L)))
        pred[[length(pred) + 1]] = data.frame(outcome = kind, model = name, TRI = grid, fit = p, 
            lo = p - qt(0.975, r$G - 1) * se, hi = p + qt(0.975, r$G - 1) * se, response_mean = mean(d[[y]]), 
            response_sd = sd(d[[y]]))
        if (name == "primary") {
            curve = r$beta[1] + sum(tail(r$beta, r$G - 1))/r$G + as.vector(X[, 2:4] %*% r$beta[2:4])
            a = data.frame(TRI = b$end_tri, adjusted = b[[y]] - r$fit + curve)
            a$bin = cut(a$TRI, unique(quantile(a$TRI, seq(0, 1, length.out = 13))), include.lowest = TRUE, 
                labels = FALSE) - 1
            ab = aggregate(cbind(TRI, adjusted) ~ bin, a, mean)
            ab$n = as.integer(table(a$bin))
            ab$outcome = kind
            bins[[length(bins) + 1]] = ab
        }
    }
}

raw = read("raw_records.csv.gz")

raw$ts = time(raw$time)

m = read("steps.csv")

ter = rbind(data.frame(ind = m$ind, ts = time(m$t), terrain = m$end_tri), data.frame(ind = m$ind, 
    ts = round(time(m$t) - m$dt_min * 60), terrain = m$start_tri))

ter = ter[!duplicated(paste(ter$ind, ter$ts)), ]

sy = list()

for (i in 1:nrow(d)) {
    r = d[i, ]
    a = raw[raw$ind == r$ind & raw$ts > time(r$window_start) & raw$ts <= time(r$window_end), ]
    tt = ter[ter$ind == r$ind & ter$ts >= time(r$window_start) & ter$ts <= time(r$window_end), ]
    dur = sum(a$raw_dt_min)
    cov = nrow(tt)/(r$window_n_steps_x + 1)
    valid = nrow(a) > 0 && all(a$raw_dt_min > 0 & a$raw_dt_min <= 2.5) && abs(dur - r$window_span_min) < 
        1e-06 && cov >= 0.8
    sy[[i]] = data.frame(row = i - 1, valid = valid, duration = dur, count = sum(a[["运动量"]]), 
        mean_TRI = mean(tt$terrain), coverage = cov)
}

sy = do.call(rbind, sy)

write(sy, "synchronized_window_audit.csv")

ok = which(sy$valid)

b = d[ok, ]

b$log_rate = log(sy$count[ok]/sy$duration[ok])

b$end_tri_z = (sy$mean_TRI[ok] - sc$end_tri$mean)/sc$end_tri$sd

one("activity", "whole_window", b)

one("activity", "whole_window_nonoverlapping", nonoverlap(b))

write(do.call(rbind, rows), "movement_checks.csv")

write(do.call(rbind, inds), "individual_slopes.csv")

write(do.call(rbind, pred), "shape_predictions.csv")

write(do.call(rbind, bins), "shape_bins.csv")

ir = list()

for (name in c("additive", "interaction")) {
    cs = c("end_tri_z", "straight_recomputed_z", "log_dt_z", "log_span_z", if (name == "interaction") "TRI_x_straightness")
    r = ols(design(d, cs), d$log_rate, d$ind)
    for (j in seq_along(cs)) ir[[length(ir) + 1]] = cbind(model = name, term = cs[j], stats(r, j + 
        1))
}

write(do.call(rbind, ir), "activity_interactions.csv")

