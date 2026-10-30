# Numerical helpers: CR1 covariance, conditional-logistic likelihood, geometry and transformations.
# See README.md for numerical conventions and scope.

read <- function(f, where = "data") read.csv(file.path(where, f), check.names = FALSE, stringsAsFactors = FALSE)

write <- function(x, f) write.csv(x, file.path("results", f), row.names = FALSE, na = "")

time <- function(x) as.numeric(as.POSIXct(x, tz = "UTC"))

truth <- function(x) tolower(as.character(x)) == "true"

design <- function(d, cols, group = "ind") {
    ids = sort(unique(d[[group]]))
    X = cbind(Intercept = 1, as.matrix(d[, cols, drop = FALSE]))
    if (length(ids) > 1) 
        for (i in ids[-1]) X = cbind(X, as.numeric(d[[group]] == i))
    X
}

ols <- function(X, y, g, w = rep(1, length(y))) {
    n = nrow(X)
    k = ncol(X)
    G = length(unique(g))
    b = qr.solve(X * sqrt(w), y * sqrt(w))
    e = as.vector(y - X %*% b)
    bread = solve(crossprod(X, w * X))
    sc = rowsum(X * (w * e), g)
    V = if (G > 1) 
        bread %*% crossprod(sc) %*% bread * G/(G - 1) * (n - 1)/(n - k)
    else matrix(NA, k, k)
    list(beta = as.vector(b), V = V, n = n, G = G, e = e, fit = as.vector(X %*% b))
}

stats <- function(r, j = 2) {
    b = r$beta[j]
    se = sqrt(r$V[j, j])
    q = if (r$G > 1) 
        qt(0.975, r$G - 1)
    else NA
    data.frame(beta = b, cluster_se = se, cluster_lo = b - q * se, cluster_hi = b + q * se, cluster_p = if (r$G > 
        1) 
        2 * pt(-abs(b/se), r$G - 1)
    else NA, n = r$n, n_bats = r$G)
}

clogit <- function(X, g) {
    n = dim(X)[1]
    J = dim(X)[2]
    k = dim(X)[3]
    for (a in 1:k) X[, , a] = X[, , a] - rowMeans(X[, , a])
    calc = function(b, hess = FALSE) {
        eta = matrix(0, n, J)
        for (a in 1:k) eta = eta + X[, , a] * b[a]
        mx = apply(eta, 1, max)
        ee = exp(eta - mx)
        prob = ee/rowSums(ee)
        mu = matrix(0, n, k)
        for (a in 1:k) mu[, a] = rowSums(prob * X[, , a])
        score = matrix(X[, 1, ], n, k) - mu
        out = list(value = mean(log(rowSums(ee)) + mx - eta[, 1]), grad = -colMeans(score), score = score)
        if (hess) {
            H = matrix(0, k, k)
            for (a in 1:k) for (b in 1:k) H[a, b] = sum(prob * X[, , a] * X[, , b]) - sum(mu[, a] * 
                mu[, b])
            out$H = (H + t(H))/2
        }
        out
    }
    o = optim(rep(0, k), function(b) calc(b)$value, function(b) calc(b)$grad, method = "BFGS", control = list(reltol = 1e-14, 
        maxit = 1000))
    b = o$par
    for (iter in 1:8) {
        a = calc(b, TRUE)
        if (max(abs(a$grad)) < 1e-10) 
            break
        b = b - solve(a$H, a$grad * n)
    }
    a = calc(b, TRUE)
    stopifnot(max(abs(a$grad)) < 1e-07, min(eigen(a$H, symmetric = TRUE)$values) > 0)
    bread = solve(a$H)
    sc = rowsum(a$score, g)
    G = length(unique(g))
    V = if (G > 1) 
        bread %*% crossprod(sc) %*% bread * G/(G - 1) * (n - 1)/(n - k)
    else matrix(NA, k, k)
    list(beta = b, V = V, V_model = bread, n = n, G = G, gradient = max(abs(a$grad)), loglik = -n * 
        a$value)
}

nonoverlap <- function(d) {
    ix = integer()
    for (g in split(seq_len(nrow(d)), interaction(d$ind, d$seg, drop = TRUE))) {
        end = -Inf
        for (i in g[order(time(d$time[g]))]) if (time(d$window_start[i]) > end) {
            ix = c(ix, i)
            end = time(d$window_end[i])
        }
    }
    d[ix, , drop = FALSE]
}

distance <- function(la, lo, lb, lp) {
    a = la * pi/180
    b = lo * pi/180
    c = lb * pi/180
    dd = lp * pi/180
    v = sin((a - c)/2)^2 + cos(a) * cos(c) * sin((b - dd)/2)^2
    2 * 6371000 * asin(sqrt(pmin(1, pmax(0, v))))
}

prepare_movement <- function(d) {
    sc <- jsonlite::fromJSON("data/movement_scaling.json")
    rec <- list(log_dt = log(d$dt_min), log_span = log(d$window_span_min))
    for (nm in names(rec)) {
        stopifnot(max(abs(rec[[nm]] - d[[nm]])) < 1e-10)
        d[[nm]] <- rec[[nm]]
    }
    for (nm in names(sc)) {
        z <- (d[[nm]] - sc[[nm]]$mean)/sc[[nm]]$sd
        stopifnot(max(abs(z - d[[paste0(nm, "_z")]])) < 1e-09)
        d[[paste0(nm, "_z")]] <- z
    }
    tt <- as.POSIXlt(d$time, tz = "UTC")
    h <- tt$hour + tt$min/60 + tt$sec/3600
    for (nm in c("clock_sin", "clock_cos", "TRI_x_straightness")) {
        val <- switch(nm, clock_sin = sin(2 * pi * h/24), clock_cos = cos(2 * pi * h/24), TRI_x_straightness = d$end_tri_z * 
            d$straight_recomputed_z)
        stopifnot(max(abs(val - d[[nm]])) < 1e-09)
        d[[nm]] <- val
    }
    d
}

