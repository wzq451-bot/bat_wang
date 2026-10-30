# Vertical coupling and circular-shift randomization from the archived lag schedule.
# See README.md for numerical conventions and scope.

m = read("steps.csv")

f = read("fixes.csv")

sg = read("segments.csv")

g = rbind(data.frame(ind = m$ind, t = time(m$t), ground = m$end_elev), data.frame(ind = m$ind, t = round(time(m$t) - 
    60 * m$dt_min), ground = m$start_elev))

key = function(i, t) paste(i, t)

g = g[!duplicated(key(g$ind, g$t)), ]

f$ground = g$ground[match(key(f$ind, time(f$t)), key(g$ind, g$t))]

v = f[f$seg %in% sg$seg & truth(f$altitude_valid) & !is.na(f$ground), ]

v = v[order(v$ind, v$seg, v$t), ]

v$x = v$ground - ave(v$ground, v$seg)

v$y = v$gps_alt - ave(v$gps_alt, v$seg)

z = aggregate(cbind(ground, gps_alt) ~ ind + seg, v, mean)

out = rbind(cbind(model = "within", stats(ols(matrix(v$x, ncol = 1), v$y, v$ind), 1)), cbind(model = "between", 
    stats(ols(cbind(1, z$ground), z$gps_alt, z$ind))))

write(out, "vertical_refit.csv")

ss = do.call(rbind, lapply(split(v, v$seg), function(a) data.frame(ind = a$ind[1], seg = a$seg[1], 
    n = nrow(a), beta = sum(a$x * a$y)/sum(a$x^2), ground_sd = sd(a$ground))))

write(ss, "vertical_segments.csv")

sh = read("circular_shift_schedule.csv")

idx = split(seq_len(nrow(v)), v$seg)

score = function(a) mean(tapply(a - v$ground, v$ind, sd))

null = matrix(NA, 2000, 3)

for (rep in 1:2000) {
    a = v$gps_alt
    s = sh[sh$replicate == rep, ]
    for (j in 1:nrow(s)) {
        ii = idx[[s$seg[j]]]
        lag = s$lag[j]
        a[ii] = v$gps_alt[ii][((seq_along(ii) - 1 - lag)%%length(ii)) + 1]
    }
    null[rep, ] = c(rep, score(a), sum(v$x * a)/sum(v$x^2))
}

colnames(null) = c("replicate", "mean_SD_AGL", "within_beta")

write(as.data.frame(null), "vertical_null.csv")

write(data.frame(observed = score(v$gps_alt), median = median(null[, 2]), p = (1 + sum(null[, 2] <= 
    score(v$gps_alt)))/2001), "vertical_summary.csv")

write(v, "vertical_analysis.csv")

