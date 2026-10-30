# Independently reconstruct original-record durations and centered-window geometry.
# See README.md for numerical conventions and scope.

d = read("activity_rows.csv.gz")

raw = read("raw_records.csv.gz")

raw = raw[order(raw$ind, time(raw$time)), ]

raw$duration_R = ave(time(raw$time), raw$ind, FUN = function(z) c(NA, diff(z)/60))

key = function(i, t) paste(i, t)

j = match(key(d$ind, time(d$time)), key(raw$ind, time(raw$time)))

stopifnot(!anyNA(j))

stopifnot(max(abs(raw$duration_R[j] - d$raw_dt_min)) < 1e-08, max(abs(raw[["运动量"]][j] - d[["运动量"]])) == 
    0)

f = read("fixes.csv")

f$ts = time(f$t)

rows = list()

for (i in 1:nrow(d)) {
    r = d[i, ]
    a = f[f$ind == r$ind & f$seg == r$seg & f$ts >= time(r$time) - 300 & f$ts <= time(r$time) + 
        300, ]
    a = a[order(a$ts), ]
    nn = nrow(a)
    path = sum(distance(a$lat[-1], a$lon[-1], a$lat[-nn], a$lon[-nn]))
    s = distance(a$lat[1], a$lon[1], a$lat[nn], a$lon[nn])/path
    rows[[i]] = data.frame(step_id = r$step_id, straightness = s, span = (a$ts[nn] - a$ts[1])/60, 
        path = path, n_steps = nn - 1, start = a$ts[1], end = a$ts[nn])
}

rr = do.call(rbind, rows)

stopifnot(max(abs(rr$straightness - d$straight_recomputed)) < 1e-08, max(abs(rr$span - d$window_span_min)) < 
    1e-08, max(abs(rr$path - d$window_path_m)) < 1e-05, all(rr$n_steps == d$window_n_steps_x))

write(rr, "reconstructed_geometry.csv")

write(data.frame(n = nrow(d), geometry_max_abs = max(abs(rr$straightness - d$straight_recomputed)), 
    duration_max_abs = max(abs(raw$duration_R[j] - d$raw_dt_min))), "reconstruction_checks.csv")

