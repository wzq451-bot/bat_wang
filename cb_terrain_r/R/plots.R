# Statistical audit PDFs, not the final cartographic or 3D figure composition.
# See README.md for numerical conventions and scope.

font <- Sys.getenv("PLOT_FONT", "sans")

dir.create("figures", showWarnings = FALSE)

open_plot <- function(name, width = 8, height = 6) {
    if (capabilities("cairo")) 
        cairo_pdf(file.path("figures", paste0(name, ".pdf")), width = width, height = height, family = font)
    else pdf(file.path("figures", paste0(name, ".pdf")), width = width, height = height, family = "Helvetica")
    par(mar = c(4.2, 4.4, 2, 1), mgp = c(2.7, 0.7, 0), las = 1, cex = 0.85, bty = "l")
}

forest <- function(z, labels, title, color = "#277C86") {
    yy <- rev(seq_len(nrow(z)))
    plot(z$beta, yy, xlim = range(c(z$cluster_lo, z$cluster_hi, 0), na.rm = TRUE), ylim = c(0.5, 
        nrow(z) + 0.5), yaxt = "n", ylab = "", xlab = "Coefficient (95% cluster CI)", pch = 16, 
        col = color, main = title)
    abline(v = 0, lty = 2, col = "grey60")
    segments(z$cluster_lo, yy, z$cluster_hi, yy, col = color)
    axis(2, at = yy, labels = labels, tick = FALSE, cex.axis = 0.7)
}

v <- read("vertical_analysis.csv", "results")

vr <- read("vertical_refit.csv", "results")

null <- read("vertical_null.csv", "results")

vs <- read("vertical_summary.csv", "results")

sg <- read("vertical_segments.csv", "results")

open_plot("vertical_audit", 10, 7)

par(mfrow = c(2, 2))

plot(v$x, v$y, pch = 16, cex = 0.2, col = adjustcolor("#277C86", 0.2), xlab = "Ground elevation deviation (m)", 
    ylab = "Flight altitude deviation (m)")

abline(a = 0, b = vr$beta[vr$model == "within"], col = "#277C86", lwd = 2)

z <- aggregate(cbind(ground, gps_alt) ~ ind + seg, v, mean)

plot(z$ground, z$gps_alt, pch = 21, xlab = "Mean ground elevation (m)", ylab = "Mean flight altitude (m)")

abline(lm(gps_alt ~ ground, z), col = "#277C86", lwd = 2)

plot(sg$ground_sd, sg$beta, pch = 16, col = "#277C86", xlab = "Segment ground elevation SD (m)", 
    ylab = "Within-segment slope")

abline(h = c(0, 1), lty = c(3, 2), col = "grey60")

hist(null$mean_SD_AGL, breaks = 35, col = "grey80", border = "white", main = "", xlim = range(c(null$mean_SD_AGL, 
    vs$observed)), xlab = "Mean within-bat AGL SD (m)")

abline(v = vs$observed, col = "#277C86", lwd = 2)

dev.off()

s <- read("selection_coefficients.csv", "results")

sp <- read("selection_predictions.csv", "results")

open_plot("selection_audit", 10, 5)

par(mfrow = c(1, 2), mar = c(4.2, 5.8, 2, 1))

a <- s[s$model == "M4_main_effects", ]

forest(a, c("TPI", "Elevation", "TRI"), "Additive terrain model")

plot(sp$agl_m, sp$cond_slope, type = "n", ylim = range(c(sp$cluster_lo, sp$cluster_hi)), xlab = "Starting height above ground (m)", 
    ylab = "Conditional TPI coefficient")

polygon(c(sp$agl_m, rev(sp$agl_m)), c(sp$cluster_lo, rev(sp$cluster_hi)), border = NA, col = adjustcolor("#277C86", 
    0.15))

lines(sp$agl_m, sp$cond_slope, col = "#277C86", lwd = 2)

abline(h = 0, lty = 2, col = "grey60")

dev.off()

open_plot("selection_sensitivity", 9, 8)

par(mar = c(4.2, 14, 2, 1))

a <- s[s$term %in% c("agl_z:tpi_z", "within_agl_z:tpi_z") & s$n_bats > 1 & s$model != "primary", 
    ]

forest(a, gsub("_", " ", a$model), "Height x TPI")

dev.off()

sh <- read("shape_predictions.csv", "results")

bins <- read("shape_bins.csv", "results")

open_plot("movement_curves", 10, 5)

par(mfrow = c(1, 2))

for (kind in c("straightness", "activity")) {
    a <- sh[sh$outcome == kind & sh$model == "primary", ]
    mu <- a$response_mean[1]
    ss <- a$response_sd[1]
    plot(a$TRI, (a$fit - mu)/ss, type = "n", ylim = range(c(a$lo, a$hi) - mu)/ss, xlab = "Terrain ruggedness index", 
        ylab = "Adjusted response (SD units)", main = kind)
    polygon(c(a$TRI, rev(a$TRI)), c(a$lo - mu, rev(a$hi - mu))/ss, border = NA, col = adjustcolor("#277C86", 
        0.15))
    for (nm in unique(sh$model)) {
        z <- sh[sh$outcome == kind & sh$model == nm, ]
        lines(z$TRI, (z$fit - mu)/ss, col = if (nm == "primary") 
            "#277C86"
        else "grey65", lty = if (nm == "equal_bat_weight") 
            2
        else 1, lwd = if (nm == "primary") 
            2
        else 0.6)
    }
    z <- bins[bins$outcome == kind, ]
    points(z$TRI, (z$adjusted - mu)/ss, pch = 21, bg = "white", col = "#277C86")
}

dev.off()

a <- read("movement_checks.csv", "results")

for (kind in c("straightness", "activity")) {
    open_plot(paste0(kind, "_sensitivity"), 9, 8)
    par(mar = c(4.2, 14, 2, 1))
    z <- a[a$outcome == kind, ]
    forest(z, gsub("_", " ", z$model), kind)
    dev.off()
}

writeLines(paste("Requested graphics font:", font), "validation/graphics_font.txt")

