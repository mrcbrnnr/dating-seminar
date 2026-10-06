# =============================================================================
# Dating Team: Shiny app for the seminar "Dating Methods"
# 14C calibration + simple Bayesian refinement (stratigraphy, TPQ/TAQ,
# dendro windows) with assumption switches, scenarios, old-wood check,
# measurement conversion (F14C), wiggle matching, SPD, OxCal export
#
# Start:   install.packages(c("shiny", "rcarbon"))
#          shiny::runApp("dating_app")
#
# The example data are FICTIONAL and meant for demonstration only.
# The model is a didactic Monte Carlo model, not a replacement for OxCal.
# =============================================================================

library(shiny)

# ---- IntCal20 laden ---------------------------------------------------------
# Order: (1) file intcal20.14c next to app.R, (2) package rcarbon,
# (3) download from intcal.org (saved next to app.R).
parse_curve <- function(df) {
  df <- as.data.frame(df)[, 1:3]
  names(df) <- c("calBP", "c14", "err")
  df[] <- lapply(df, function(x) suppressWarnings(as.numeric(x)))
  df <- df[stats::complete.cases(df), ]
  if (nrow(df) < 1000) stop("Curve incomplete")
  df
}

load_curve <- function() {
  ic <- NULL
  f <- file.path(getwd(), "intcal20.14c")

  # (1) lokale Datei
  if (file.exists(f)) {
    ic <- tryCatch(parse_curve(read.csv(f, comment.char = "#", header = FALSE)), error = function(e) NULL)
  }
  # (2) rcarbon (local only; not available in browser mode with Shinylive/webR)
  is_wasm <- grepl("emscripten", R.version$platform, fixed = TRUE)
  RC <- "rcarbon"   # name deliberately kept in a variable so that Shinylive does not try to bundle the package
  if (is.null(ic) && !is_wasm && requireNamespace(RC, quietly = TRUE)) {
    for (nm in c("intcal20", "IntCal20")) {
      x <- tryCatch(getExportedValue(RC, nm), error = function(e) NULL)
      if (!is.null(x)) {
        ic <- tryCatch(parse_curve(x), error = function(e) NULL)
        if (!is.null(ic)) break
      }
    }
  }
  # (3) Download
  if (is.null(ic) && !is_wasm) {
    ok <- tryCatch({
      download.file("https://intcal.org/curves/intcal20.14c", f, quiet = TRUE, mode = "wb")
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && file.exists(f)) {
      ic <- tryCatch(parse_curve(read.csv(f, comment.char = "#", header = FALSE)), error = function(e) NULL)
    }
  }
  if (is.null(ic)) {
    stop(paste0(
      "IntCal20 could not be loaded. Options: (a) install.packages('rcarbon') and restart the app, ",
      "(b) download https://intcal.org/curves/intcal20.14c and save it as 'intcal20.14c' in the app folder ",
      "(for the web version the file must be in the app/ folder of the repository). ",
      "Working directory: ", getwd()))
  }
  ic <- ic[ic$calBP <= 15000, ]
  ic$year <- 1950 - ic$calBP          # calendar year, negative = BC (1-year uncertainty: no year 0)
  ic[order(ic$year), ]
}
CURVE <- load_curve()

# ---- Years: input (astronomical OR "480 BC" / "101 AD") and display ----
# Internal counting as in OxCal: astronomical, year 0 = 1 BC, -479 = 480 BC
parse_year1 <- function(x) {
  if (length(x) != 1 || is.na(x)) return(NA_real_)
  if (is.numeric(x)) return(as.numeric(x))
  s <- trimws(as.character(x))
  if (!nzchar(s)) return(NA_real_)
  m <- regmatches(s, regexec("^(-?[0-9]+)\\s*(v\\.?\\s*Chr\\.?|BCE?|n\\.?\\s*Chr\\.?|CE|AD)$", s, ignore.case = TRUE))[[1]]
  if (length(m) == 3) {
    n <- abs(as.numeric(m[2]))
    era <- tolower(gsub("[ .]", "", m[3]))
    if (era %in% c("vchr", "bc", "bce")) return(1 - n)
    return(n)
  }
  m2 <- regmatches(s, regexec("^(AD|CE)\\s*([0-9]+)$", s, ignore.case = TRUE))[[1]]
  if (length(m2) == 3) return(as.numeric(m2[3]))
  suppressWarnings(as.numeric(s))
}
parse_year <- function(x) vapply(x, parse_year1, numeric(1), USE.NAMES = FALSE)

fmt_year <- function(y, mode = "astro") {
  if (identical(mode, "bcad")) ifelse(y <= 0, paste0(1 - y, " BC"), paste0(y, " AD"))
  else as.character(y)
}
xlab_years <- function(mode) if (identical(mode, "bcad")) "Calendar year" else "Calendar year (astronomical: -479 = 480 BC)"
axis_years <- function(rng, mode = "astro") {
  if (identical(mode, "bcad")) {
    # ticks on round historical numbers (e.g. 500 BC), placed internally in astronomical years
    h_rng <- ifelse(rng <= 0, rng - 1, rng)
    h <- pretty(h_rng, n = 6)
    h <- h[h != 0]
    at <- ifelse(h < 0, h + 1, h)
    keep <- at >= rng[1] & at <= rng[2]
    h <- h[keep]; at <- at[keep]
    labels <- ifelse(h < 0, paste0(-h, " BC"), paste0(h, " AD"))
  } else {
    at <- pretty(rng, n = 6)
    at <- at[at >= rng[1] & at <= rng[2]]
    labels <- as.character(at)
  }
  axis(1, at = at, labels = labels, cex.axis = 0.9)
}

# ---- Hilfsfunktionen --------------------------------------------------------
# English column names are mapped to the internal ones (the German names are accepted as well)
norm_names <- function(df) {
  map <- c(layer = "schicht", rank = "rang", type = "typ", value1 = "wert1", value2 = "wert2", offset = "abstand", error = "fehler")
  nm <- tolower(trimws(names(df)))
  idx <- match(nm, names(map))
  nm[!is.na(idx)] <- unname(map[idx[!is.na(idx)]])
  names(df) <- nm
  df
}
read_txt <- function(txt) {
  txt <- trimws(txt)
  if (!nzchar(txt)) return(NULL)
  first <- strsplit(txt, "\n")[[1]][1]
  sep <- if (grepl(";", first)) ";" else if (grepl("\t", first)) "\t" else ","
  df <- read.table(text = txt, header = TRUE, sep = sep, strip.white = TRUE,
                   stringsAsFactors = FALSE, comment.char = "",
                   na.strings = c("", "NA"), fill = TRUE)
  norm_names(df)
}

curve_at <- function(yrs) {
  list(mu = approx(CURVE$year, CURVE$c14, xout = yrs, rule = 2)$y,
       sc = approx(CURVE$year, CURVE$err, xout = yrs, rule = 2)$y)
}

# Calibration: probability per calendar year (sum = 1)
dens_c14 <- function(age, err, yrs) {
  cv <- curve_at(yrs)
  s2 <- err^2 + cv$sc^2
  d <- exp(-0.5 * (age - cv$mu)^2 / s2) / sqrt(s2)
  d / sum(d)
}

smooth5 <- function(d, k = 5) {
  s <- as.numeric(stats::filter(d, rep(1 / k, k), sides = 2))
  s[is.na(s)] <- 0
  if (sum(s) > 0) s / sum(s) else s
}

# Smoothing with edge padding (for distance distributions that start at 0)
smooth_edge <- function(d, k = 5) {
  n <- length(d)
  if (n < 6) return(d)
  d2 <- c(rev(d[1:2]), d, rev(d[(n - 1):n]))
  s <- as.numeric(stats::filter(d2, rep(1 / k, k), sides = 2))[3:(n + 2)]
  s[is.na(s)] <- 0
  if (sum(s) > 0) s / sum(s) else s
}

# 95.4% HPD ranges as text
hpd <- function(d, yrs, p = 0.954, gap = 3, mode = "astro") {
  if (sum(d) <= 0) return(NA_character_)
  d <- d / sum(d)
  o <- order(d, decreasing = TRUE)
  keep <- o[seq_len(which(cumsum(d[o]) >= p)[1])]
  y <- sort(yrs[keep])
  br <- c(0, which(diff(y) > gap), length(y))
  segs <- vapply(seq_len(length(br) - 1), function(j) {
    a <- y[br[j] + 1]; b <- y[br[j + 1]]
    if (a == b) fmt_year(a, mode) else paste0(fmt_year(a, mode), " to ", fmt_year(b, mode))
  }, character(1))
  paste(segs, collapse = "; ")
}

is_c14 <- function(typ) toupper(trimws(typ)) %in% c("14C", "C14")

# Time grid for a single 14C age (for the conversion tab)
age_grid <- function(age, err) {
  s <- sqrt(err^2 + CURVE$err^2)
  idx <- which(abs(age - CURVE$c14) < 4 * s)
  if (length(idx) == 0) stop("Age outside the calibration curve")
  seq(floor(min(CURVE$year[idx])) - 25, ceiling(max(CURVE$year[idx])) + 25)
}

# Derive the time grid from the data
make_grid <- function(smp, lay, upper) {
  pts <- upper
  for (i in seq_len(nrow(smp))) {
    if (is_c14(smp$typ[i])) {
      s <- sqrt(smp$wert2[i]^2 + CURVE$err^2)
      idx <- which(abs(smp$wert1[i] - CURVE$c14) < 4 * s)
      if (length(idx) == 0) stop(paste0("Sample ", smp$id[i], ": 14C age is outside the calibration curve."))
      pts <- c(pts, range(CURVE$year[idx]))
    } else {
      pts <- c(pts, smp$wert1[i], smp$wert2[i])
    }
  }
  pts <- c(pts, lay$tpq[!is.na(lay$tpq)], lay$taq[!is.na(lay$taq)])
  lo <- floor(min(pts)) - 25
  hi <- ceiling(max(pts)) + 25
  if (hi - lo > 5000) stop("The time span is too long (> 5000 years). Check the upper limit or the data.")
  seq(lo, hi)
}

# ---- Check and convert tables ------------------------------------------------
# Years may be entered astronomically ("-479") or as text ("480 BC").
# For dendro samples value1/value2 are years, for 14C samples age BP and error.
coerce_tables <- function(smp, lay) {
  need_s <- c("id", "schicht", "typ", "wert1", "wert2")
  need_l <- c("schicht", "rang", "tpq", "taq")
  if (is.null(smp) || !all(need_s %in% names(smp))) stop("Sample table: columns id, layer, type, value1, value2 are required.")
  if (is.null(lay) || !all(need_l %in% names(lay))) stop("Layer table: columns layer, rank, tpq, taq are required.")
  smp <- smp[!is.na(smp$id), , drop = FALSE]
  lay <- lay[!is.na(lay$schicht), , drop = FALSE]
  smp$schicht <- as.character(smp$schicht); lay$schicht <- as.character(lay$schicht)
  smp$typ <- as.character(smp$typ)
  dend <- !is_c14(smp$typ)
  w1 <- as.character(smp$wert1); w2 <- as.character(smp$wert2)
  smp$wert1 <- ifelse(dend, parse_year(w1), suppressWarnings(as.numeric(w1)))
  smp$wert2 <- ifelse(dend, parse_year(w2), suppressWarnings(as.numeric(w2)))
  lay$rang <- suppressWarnings(as.numeric(as.character(lay$rang)))
  lay$tpq <- parse_year(lay$tpq)
  lay$taq <- parse_year(lay$taq)
  if (nrow(smp) == 0) stop("Sample table is empty.")
  if (anyNA(smp$wert1) || anyNA(smp$wert2)) stop("Sample table: value1/value2 must be numbers (for dendro samples also '480 BC' is possible).")
  if (anyNA(lay$rang)) stop("Layer table: 'rank' must be a number for every layer (1 = oldest/lowest).")
  if (!all(smp$schicht %in% lay$schicht)) stop("There are samples whose layer is not in the layer table.")
  list(smp = smp, lay = lay)
}

# Switch assumptions on and off (for the before/after comparison)
apply_assumptions <- function(smp, lay, opt) {
  if (!isTRUE(opt$tpq)) lay$tpq <- NA_real_
  if (!isTRUE(opt$taq)) lay$taq <- NA_real_
  if (!isTRUE(opt$dendro)) smp <- smp[is_c14(smp$typ), , drop = FALSE]
  if (!isTRUE(opt$strat)) lay$rang <- 1
  list(smp = smp, lay = lay)
}
OPT_ALL <- list(strat = TRUE, tpq = TRUE, taq = TRUE, dendro = TRUE)
OPT_NONE <- list(strat = FALSE, tpq = FALSE, taq = FALSE, dendro = FALSE)
OFFMAX <- 500L

# ---- The model -------------------------------------------------------------
# Assumptions (didactic):
#  * Every layer has a deposition date D.
#  * Every sample (14C or dendro) dates an event S that does NOT lie after D (S <= D).
#  * Coin TPQ: D >= TPQ.  Written source/sealing TAQ: D <= TAQ.
#  * Stratigraphy: deposition of a lower (smaller rank) layer <= that of the upper one.
#  * Uniform priors; calculation by Monte Carlo drawing and rejection.
run_model <- function(smp, lay, upper, N = 40000, opt = OPT_ALL, yrs = NULL) {
  ct <- coerce_tables(smp, lay)
  smp <- ct$smp; lay <- ct$lay
  if (is.null(yrs)) yrs <- make_grid(smp, lay, upper)   # grid from ALL data so that the switches do not shift the axis
  ap <- apply_assumptions(smp, lay, opt)
  smp <- ap$smp; lay <- ap$lay
  if (nrow(smp) == 0) stop("No samples left (all dendro samples switched off?).")
  ng <- length(yrs); ns <- nrow(smp); nl <- nrow(lay)

  L <- lapply(seq_len(ns), function(i) {
    if (is_c14(smp$typ[i])) {
      dens_c14(smp$wert1[i], smp$wert2[i], yrs)
    } else {
      d <- as.numeric(yrs >= min(smp$wert1[i], smp$wert2[i]) & yrs <= max(smp$wert1[i], smp$wert2[i]))
      if (sum(d) == 0) stop(paste0("Sample ", smp$id[i], ": dendro window is not within the time grid."))
      d / sum(d)
    }
  })
  cdf <- lapply(L, cumsum)
  members <- lapply(seq_len(nl), function(k) which(smp$schicht == lay$schicht[k]))

  # Marginal distribution of the deposition D per layer: P(D=d) ~ product of the cumulative distributions
  Dm <- matrix(0L, N, nl)
  for (k in seq_len(nl)) {
    tpq <- if (is.na(lay$tpq[k])) -Inf else lay$tpq[k]
    taq <- if (is.na(lay$taq[k])) upper else lay$taq[k]
    if (is.finite(tpq) && tpq > taq) {
      stop(paste0("Conflict in layer ", lay$schicht[k], ": the TPQ (", tpq, ") lies after the TAQ (", taq,
                  "). Is the coin intrusive, is the TAQ wrong, or is the layer not a single deposit? (years counted astronomically)"))
    }
    w <- as.numeric(yrs >= tpq & yrs <= taq)
    for (i in members[[k]]) w <- w * cdf[[i]]
    if (sum(w) <= 0) {
      stop(paste0("Conflict in layer ", lay$schicht[k],
                  ": TPQ/TAQ and the samples are mutually exclusive. At least one sample is younger than the TAQ or older than allowed."))
    }
    Dm[, k] <- sample.int(ng, N, replace = TRUE, prob = w / sum(w))
  }

  # Stratigraphy: reject all draws that violate the order
  ok <- rep(TRUE, N)
  prevmax <- rep(-Inf, N)
  for (r in sort(unique(lay$rang))) {
    ks <- which(lay$rang == r)
    for (k in ks) ok <- ok & (Dm[, k] >= prevmax)
    prevmax <- pmax(prevmax, apply(Dm[, ks, drop = FALSE], 1, max))
  }
  Dacc <- Dm[ok, , drop = FALSE]
  na <- nrow(Dacc)
  if (na < 50) {
    stop("The stratigraphic order contradicts the data (almost all draws rejected). Check rank, TPQ/TAQ and samples.")
  }

  postD <- matrix(sapply(seq_len(nl), function(k) smooth5(tabulate(Dacc[, k], ng))), nrow = ng)
  postS <- matrix(0, ng, ns)
  off_med <- rep(NA_real_, ns); off_q95 <- rep(NA_real_, ns); off_dens <- vector("list", ns)   # distance to deposition
  rel_med <- rep(NA_real_, ns); rel_q95 <- rep(NA_real_, ns); rel_dens <- vector("list", ns)   # distance to the youngest sample of the layer
  for (k in seq_len(nl)) {
    mem <- members[[k]]
    if (!length(mem)) next
    Di <- Dacc[, k]
    IDX <- matrix(0L, na, length(mem))
    for (j in seq_along(mem)) {
      i <- mem[j]
      u <- runif(na) * cdf[[i]][Di]
      idx <- pmin(findInterval(u, cdf[[i]]) + 1L, Di)
      IDX[, j] <- idx
      postS[, i] <- smooth5(tabulate(idx, ng))
      off <- Di - idx                                   # years between sample event and deposition
      qs <- as.numeric(quantile(off, c(0.5, 0.95)))
      off_med[i] <- qs[1]; off_q95[i] <- qs[2]
      off_dens[[i]] <- smooth_edge(tabulate(pmin(off, OFFMAX) + 1L, OFFMAX + 1L))
    }
    if (length(mem) > 1) {
      ref <- apply(IDX, 1, max)                         # youngest sample event of the layer in each draw
      for (j in seq_along(mem)) {
        i <- mem[j]
        rel <- ref - IDX[, j]
        qs <- as.numeric(quantile(rel, c(0.5, 0.95)))
        rel_med[i] <- qs[1]; rel_q95[i] <- qs[2]
        rel_dens[[i]] <- smooth_edge(tabulate(pmin(rel, OFFMAX) + 1L, OFFMAX + 1L))
      }
    }
  }
  list(yrs = yrs, smp = smp, lay = lay, L = L, postS = postS, postD = postD,
       members = members, acc = na / N, upper = upper, opt = opt,
       off_med = off_med, off_q95 = off_q95, off_dens = off_dens,
       rel_med = rel_med, rel_q95 = rel_q95, rel_dens = rel_dens)
}

prep_tables <- function(smp, lay) coerce_tables(smp, lay)

# ---- OxCal-Code erzeugen ----------------------------------------------------
# Translation of the app model into an OxCal sequence:
#  * rank -> Sequence with Boundary between the layers (same rank: Phase group)
#  * 14C -> R_Date, dendro window -> Date with U(from,to)
#  * coin TPQ -> Date with U(tpq, upper limit) within the phase of the layer
#  * TAQ -> Date with U(lower limit, taq) directly AFTER the phase of the layer
# Note: OxCal models phases with boundaries, the app models the deposition D after
# the samples. The results are therefore similar but not identical.
oxcal_code <- function(smp, lay, upper, lo, outlier = FALSE) {
  out <- character()
  add <- function(...) out <<- c(out, paste0(...))
  q <- function(x) gsub("[\"\\\\]", "", as.character(x))
  num <- function(x) format(x, scientific = FALSE, trim = TRUE)

  layer_block <- function(k, ind) {
    sid <- which(smp$schicht == lay$schicht[k])
    nm <- q(lay$schicht[k])
    has_taq <- !is.na(lay$taq[k])
    if (has_taq) {
      add(ind, 'Sequence("', nm, ' with TAQ")')
      add(ind, "{")
      ind <- paste0(ind, " ")
    }
    add(ind, 'Phase("', nm, '")')
    add(ind, "{")
    for (i in sid) {
      if (is_c14(smp$typ[i])) {
        if (outlier) {
          add(ind, ' R_Date("', q(smp$id[i]), '",', num(smp$wert1[i]), ",", num(smp$wert2[i]), ") { Outlier(0.05); };")
        } else {
          add(ind, ' R_Date("', q(smp$id[i]), '",', num(smp$wert1[i]), ",", num(smp$wert2[i]), ");")
        }
      } else {
        add(ind, ' Date("', q(smp$id[i]), '",U(', num(min(smp$wert1[i], smp$wert2[i])), ",",
            num(max(smp$wert1[i], smp$wert2[i])), "));")
      }
    }
    if (!is.na(lay$tpq[k])) {
      add(ind, ' Date("TPQ coin ', nm, '",U(', num(lay$tpq[k]), ",", num(max(upper, lay$tpq[k] + 1)), "));")
    }
    add(ind, "};")
    if (has_taq) {
      add(ind, 'Date("TAQ ', nm, '",U(', num(lo), ",", num(lay$taq[k]), "));")
      ind <- substring(ind, 2)
      add(ind, "};")
    }
  }

  add("// OxCal code generated by the Shiny app 'Dating Team'")
  add("// Paste into https://c14.arch.ox.ac.uk/oxcal/OxCal.html (Input tab), then Run.")
  add("// Years are astronomical as in OxCal (year 0 = 1 BC, -78 = 79 BC). Curve: OxCal default (IntCal20).")
  add("Options()")
  add("{")
  add(" Resolution=1;")
  add("};")
  if (outlier) add('Outlier_Model("General",T(5),U(0,4),"t");')
  add("Plot()")
  add("{")
  add(' Sequence("Stratigraphy")')
  add(" {")
  add('  Boundary("Start");')
  ranks <- sort(unique(lay$rang))
  for (ri in seq_along(ranks)) {
    ks <- which(lay$rang == ranks[ri])
    if (length(ks) > 1) {
      add('  Phase("Rank ', num(ranks[ri]), '")')
      add("  {")
      for (k in ks) layer_block(k, "   ")
      add("  };")
    } else {
      layer_block(ks, "  ")
    }
    if (ri < length(ranks)) {
      add('  Boundary("', num(ranks[ri]), "/", num(ranks[ri + 1]), '");')
    }
  }
  add('  Boundary("End");')
  add(" };")
  add("};")
  out
}

# ---- Wiggle-Matching --------------------------------------------------------
wig_T <- function(start, d) {
  cv <- curve_at(start + d$abstand)
  sum((d$c14 - cv$mu)^2 / (d$fehler^2 + cv$sc^2))
}

# ---- Dendro crossdating (simulated reference, schematic) ---------------------
# The reference chronology is SIMULATED (the same reference for everybody). It shows the principle,
# but does not replace a real standard chronology.
with_seed <- function(seed, expr) {
  old <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) get(".Random.seed", envir = .GlobalEnv) else NULL
  on.exit({
    if (is.null(old)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", old, envir = .GlobalEnv)
  })
  set.seed(seed)
  expr
}
dendro_make_ref <- function(n = 300, seed = 2026) {
  with_seed(seed, {
    clim <- as.numeric(arima.sim(list(ar = 0.35), n = n))
    clim <- (clim - mean(clim)) / sd(clim)
    yy <- seq_len(n)
    low <- 0.25 * sin(2 * pi * yy / 170) + 0.15 * sin(2 * pi * yy / 60)
    w <- exp(0.30 * clim + low)
    w / mean(w)
  })
}
dendro_make_sample <- function(ref, n, end, noise, seed) {
  with_seed(seed, {
    seg <- ref[(end - n + 1):end]
    j <- seq_len(n)
    trend <- 0.7 + 1.3 * exp(-j / 35)            # age trend: young rings are wide
    w <- seg * trend * exp(rnorm(n, 0, noise))  # Standortrauschen
    w / mean(w)
  })
}
# Detrending: ratio to the moving average (9 years); at the edges with a shortened window,
# so that the sample keeps its full length (the last ring is the felling year!)
dendro_index <- function(x, k = 9) {
  n <- length(x); h <- (k - 1) %/% 2
  ma <- vapply(seq_len(n), function(i) mean(x[max(1, i - h):min(n, i + h)]), numeric(1))
  x / ma
}
# Statistic for a position e (= end year of the sample in the reference)
dendro_stats <- function(refidx, smpidx, e) {
  na4 <- c(r = NA_real_, t = NA_real_, glk = NA_real_, n = NA_real_)
  n <- length(smpidx)
  if (e < n || e > length(refidx)) return(na4)
  a <- refidx[(e - n + 1):e]; b <- smpidx
  ok <- !is.na(a) & !is.na(b)
  m <- sum(ok)
  if (m < 10) return(na4)
  r <- suppressWarnings(cor(a[ok], b[ok]))
  tt <- if (is.finite(r) && abs(r) < 1) r * sqrt(m - 2) / sqrt(1 - r^2) else NA_real_
  da <- diff(a); db <- diff(b)
  ok2 <- !is.na(da) & !is.na(db)
  glk <- if (any(ok2)) 100 * mean(sign(da[ok2]) == sign(db[ok2])) else NA_real_
  c(r = r, t = tt, glk = glk, n = m)
}
DENDRO_REF <- dendro_make_ref()
DENDRO_REF_IDX <- dendro_index(DENDRO_REF)
# Sapwood rings in oak, statistical ranges (after Tegel et al. 2022)
DENDRO_SAPWOOD <- list(
  "Great Britain: 10–55 sapwood rings" = c(10, 55),
  "Western Germany: 9–33 sapwood rings" = c(9, 33),
  "Northern Germany: 10–30 sapwood rings" = c(10, 30),
  "Poland: 9–23 sapwood rings" = c(9, 23))
# Felling-year window from the last measured ring
dendro_window <- function(end, mode, lo, hi, s, margin) {
  if (identical(mode, "wk")) return(c(end - margin, end + margin))
  if (identical(mode, "splint")) return(c(end + max(lo - s, 0), end + max(hi - s, 0)))
  c(end + lo, NA_real_)   # heartwood only: TPQ only
}
nz <- function(x, d) if (is.null(x) || length(x) != 1 || is.na(x)) d else x

# ---- Szenarien (alle fiktiv) -------------------------------------------------
SCEN <- list(
  "Fort Musterberg (with conflict in layer C)" = list(
    smp = "id,layer,type,value1,value2
A1,A,Dendro,101,103
A2,A,14C,1930,25
A3,A,14C,1900,25
B1,B,14C,1850,25
B2,B,14C,1820,25
B3,B,14C,1790,30
C1,C,14C,1760,25
C2,C,14C,1740,30
C3,C,14C,1900,25",
    lay = "layer,rank,tpq,taq
A,1,,
B,2,138 AD,
C,3,271 AD,260 AD",
    upper = 400,
    desc = "Timber-and-earth fort with three layers. Layer C contains a coin minted no earlier than 271, but the written source says the fort was abandoned by 260 at the latest. Task: find the conflict and test hypotheses (intrusive coin? wrong TAQ?). The 14C ages are illustrative."),
  "Fort Musterberg (conflict resolved: coin intrusive)" = list(
    smp = "id,layer,type,value1,value2
A1,A,Dendro,101,103
A2,A,14C,1930,25
A3,A,14C,1900,25
B1,B,14C,1850,25
B2,B,14C,1820,25
B3,B,14C,1790,30
C1,C,14C,1760,25
C2,C,14C,1740,30
C3,C,14C,1900,25",
    lay = "layer,rank,tpq,taq
A,1,,
B,2,138 AD,
C,3,,260 AD",
    upper = 400,
    desc = "As above, but the coin in layer C is treated as intrusive (TPQ removed). Task: switch the assumptions off one by one and observe what changes for B and C. Is C3 old wood?"),
  "Plateau example (Hallstatt period)" = list(
    smp = "id,layer,type,value1,value2
A1,A,14C,2520,30
A2,A,14C,2480,30
B1,B,14C,2450,30
B2,B,14C,2480,25
B3,B,Dendro,545 BC,540 BC
C1,C,14C,2420,30
C2,C,14C,2600,30",
    lay = "layer,rank,tpq,taq
A,1,,
B,2,,
C,3,480 BC,",
    upper = -300,
    desc = "All 14C ages lie on the plateau (approx. 800–400 BC): broad, multi-peaked distributions. Task: How strongly do stratigraphy, coin (480 BC) and the dendro window narrow the data? Try the switches one at a time."),
  "Old-wood exercise (burnt pit)" = list(
    smp = "id,layer,type,value1,value2
G1_seed,Pit,14C,1850,25
G2_seed,Pit,14C,1840,25
G3_beam,Pit,14C,2050,25",
    lay = "layer,rank,tpq,taq
Pit,1,,300 AD",
    upper = 400,
    desc = "A burnt pit with two seeds and a piece of timber. The timber is clearly older. Task: In tab “5 · Old-wood check” see how much older each sample is than the youngest sample in the pit. Which sample dates the pit, which only the wood?")
)
SCEN[["Collective grave (exercise after M. Hinz 2012)"]] <- list(
  smp = "id,layer,type,value1,value2
V1,PreUse,14C,4680,50
V2,PreUse,14C,4485,40
N1,Phase1,14C,4415,29
N2,Phase1,14C,4395,34
N3,Phase2,14C,4355,40
N4,Phase2,14C,4375,31
O1,Offering,14C,4325,40
O2,Offering,14C,4340,40
O3,Offering,14C,4335,40",
  lay = "layer,rank,tpq,taq
PreUse,1,,
Phase1,2,,
Phase2,3,,
Offering,4,,",
  upper = -2500,
  desc = "Collective grave with a pre-use phase (below the grave), two phases of use and an offering (three grains of cereal in a pot) after the use. Task: run the overall model and estimate the duration of use. Tab 3, “Combination”: the three grains come from the same event. Exercise data after M. Hinz (seminar 2012).")
SCEN_FIRST <- SCEN[[1]]

# ---- UI ---------------------------------------------------------------------
PAL <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02", "#a6761d", "#666666")

ui <- fluidPage(
  tags$head(tags$style(HTML("
    textarea { font-family: monospace; font-size: 12px; }
    .hint { color: #555; font-size: 13px; }
    .box { background: #f6f8fa; border: 1px solid #d0d7de; border-radius: 6px; padding: 8px 12px; margin-bottom: 10px; }
  "))),
  titlePanel("Dating Team: 14C, Dendro, Coins and Stratigraphy"),
  fluidRow(column(12, div(style = "text-align:right",
    radioButtons("yrfmt", "Show years as", choices = c("BC / AD" = "bcad", "astronomical (−479)" = "astro"), selected = "bcad", inline = TRUE)))),
  tabsetPanel(
    tabPanel("1 · Calibration & model",
      sidebarLayout(
        sidebarPanel(width = 4,
          h4("Scenario"),
          selectInput("scen", NULL, choices = c(names(SCEN), "Own data (leave tables unchanged)")),
          actionButton("load_scen", "Load scenario"),
          uiOutput("scen_desc"),
          h4("Samples"),
          p(class = "hint", "type = 14C: value1 = 14C age (BP), value2 = error. type = Dendro: value1/value2 = window of the felling year."),
          textAreaInput("txt_smp", NULL, SCEN_FIRST$smp, rows = 9, width = "100%"),
          fileInput("f_smp", "or upload CSV", accept = c(".csv", ".txt")),
          h4("Layers"),
          p(class = "hint", "rank: 1 = oldest/lowest layer. tpq: youngest coin. taq: written source/sealing. Years as “480 BC”, “101 AD” or astronomical (−479)."),
          textAreaInput("txt_lay", NULL, SCEN_FIRST$lay, rows = 5, width = "100%"),
          fileInput("f_lay", "or upload CSV", accept = c(".csv", ".txt")),
          div(class = "box",
            strong("Year converter"),
            fluidRow(
              column(5, numericInput("yw_n", NULL, 480, min = 0, step = 1)),
              column(7, radioButtons("yw_era", NULL, choices = c("BC" = "bc", "AD" = "ad"), inline = TRUE))),
            textOutput("yw_out")),
          numericInput("upper", "Latest conceivable deposition (if no TAQ), astronomical", SCEN_FIRST$upper, step = 10),
          div(class = "box",
            strong("Assumptions"),
            checkboxInput("a_strat", "Apply stratigraphic order", TRUE),
            checkboxInput("a_tpq", "Apply coin TPQ", TRUE),
            checkboxInput("a_taq", "Apply written source / TAQ", TRUE),
            checkboxInput("a_dendro", "Use dendro window", TRUE),
            checkboxInput("cmp", "Comparison: deposition from the 14C data only (grey)", TRUE)),
          actionButton("run", "Run model", class = "btn-primary"),
          br(), br(),
          p(class = "hint", "Dashed: calibration only. Coloured fill: sample in the model. Black: modelled deposition date. Grey: deposition from the 14C data only. Red: TPQ, blue: TAQ.")
        ),
        mainPanel(width = 8,
          textOutput("acc"),
          plotOutput("plot_model", height = "auto"),
          h4("Results (95.4% ranges)"),
          tableOutput("tab"),
          downloadButton("dl", "Table as CSV")
        )
      )
    ),
    tabPanel("2 · Wiggle matching",
      sidebarLayout(
        sidebarPanel(width = 4,
          h4("Simulate data"),
          numericInput("w_true", "True start year (simulation only), astronomical", -640, min = -1100, max = -150),
          numericInput("w_n", "Number of samples", 7, min = 3, max = 15),
          numericInput("w_step", "Years between samples (rings)", 10, min = 1, max = 50),
          numericInput("w_err", "Measurement error (14C years)", 25, min = 10, max = 60),
          actionButton("sim", "New simulation"),
          h4("Data (or paste your own)"),
          p(class = "hint", "offset = years since the first sample (ring count)."),
          textAreaInput("w_txt", NULL, "", rows = 9, width = "100%"),
          hr(),
          sliderInput("w_start", "Calendar year of the first sample (slider, astronomical)", min = -1200, max = -100, value = -800, step = 1, width = "100%"),
          checkboxInput("w_show_T", "Show fit curve (T) for all start years", FALSE),
          checkboxInput("w_show_true", "Show solution", FALSE)
        ),
        mainPanel(width = 8,
          p("Task: The samples come from a ring sequence with a known spacing. Slide the sequence along the calibration curve until the points follow the curve."),
          plotOutput("plot_wig", height = "380px"),
          textOutput("wig_text"),
          conditionalPanel("input.w_show_T", plotOutput("plot_T", height = "260px"))
        )
      )
    ),
    tabPanel("3 · SPD & combination",
      sidebarLayout(
        sidebarPanel(width = 4,
          radioButtons("sum_mode", NULL, choices = c("Sum: different events" = "sum", "Combination: same event" = "cmb")),
          conditionalPanel("input.sum_mode == 'sum'",
            checkboxInput("spd_mod", "Also show modelled distributions", FALSE),
            p(class = "hint", "The 14C samples from tab 1 are summed. A real SPD needs dozens to hundreds of dates and a permutation test (see code on the right)."),
            p(class = "hint", "Summed curves do not automatically show population or activity: sample selection, plateaus and taphonomy shape the curve.")),
          conditionalPanel("input.sum_mode == 'cmb'",
            uiOutput("cmb_layer_ui"),
            p(class = "hint", "Combination applies only if all measurements date the same event (e.g. several grains from one pot, several measurements on one skeleton). The ages are then averaged with weights and only afterwards calibrated. The χ² test checks whether the measurements agree."))
        ),
        mainPanel(width = 8,
          conditionalPanel("input.sum_mode == 'sum'",
            plotOutput("plot_spd", height = "340px"),
            h4("For larger data sets in R (rcarbon)"),
            tags$pre("library(rcarbon)
x <- calibrate(x = c14$age, errors = c14$error, calCurves = 'intcal20')
s <- spd(x, timeRange = c(2800, 2200))       # cal BP
plot(s)
# Permutation test against a null model:
# modelTest(x, errors = c14$error, bins = bins, nsim = 500,
#           timeRange = c(2800, 2200), model = 'exponential', runm = 100)")),
          conditionalPanel("input.sum_mode == 'cmb'",
            plotOutput("plot_cmb", height = "340px"),
            h4(textOutput("cmb_head")),
            textOutput("cmb_chi"),
            textOutput("cmb_hpd"))
        )
      )
    ),
    tabPanel("4 · OxCal code",
      br(),
      fluidRow(column(12,
        p("The code is generated live from the tables in tab 1 and takes the assumptions chosen there into account. Paste it into OxCal online (Input tab) and press “Run”. Then compare the results with those of the app."),
        checkboxInput("ox_out", "Add outlier model (General, 5% prior)", FALSE),
        tags$button("Copy to clipboard", class = "btn btn-default",
                    onclick = "navigator.clipboard.writeText(document.getElementById('oxcal_code').innerText)"),
        downloadButton("dl_ox", "Save as .oxcal"),
        tags$a("Open OxCal online", href = "https://c14.arch.ox.ac.uk/oxcal/OxCal.html", target = "_blank", class = "btn btn-link"),
        br(), br(),
        verbatimTextOutput("oxcal_code"),
        p(class = "hint", "Translation: rank → Sequence with boundaries; ¹⁴C → R_Date; dendro window → Date with U(from,to); coin → Date U(TPQ,upper limit) within the phase; TAQ → Date U(lower limit,TAQ) directly after the phase. OxCal models phases with boundaries, the app models deposition after the samples: differences in the results are to be expected and make a good discussion question.")
      ))
    ),
    tabPanel("5 · Old-wood check",
      sidebarLayout(
        sidebarPanel(width = 4,
          p("For each sample: how many years lie between the dated event (death, felling year) and a reference point? Large distances point to old wood or reworked material."),
          radioButtons("off_ref", "Reference point", choices = c("youngest sample of the same layer" = "rel", "deposition of the layer" = "dep"), selected = "rel"),
          sliderInput("off_max", "Show up to (years of distance)", min = 50, max = 500, value = 300, step = 25),
          sliderInput("off_thr", "Flag from median distance (years)", min = 10, max = 300, value = 60, step = 10),
          p(class = "hint", "“Youngest sample of the layer” is robust but needs at least two samples in the layer. “Deposition” depends strongly on the TAQ or the upper limit: without a fixed upper limit the distance becomes large for all samples."),
          p(class = "hint", "On plateaus the distance becomes broad even for short-lived samples because the calibration is blurred. A high median is a hint, not proof."),
          p(class = "hint", "The model from tab 1 is used (press “Run model” there).")
        ),
        mainPanel(width = 8,
          plotOutput("plot_off", height = "auto"),
          h4("Distance of the samples"),
          tableOutput("tab_off")
        )
      )
    ),
    tabPanel("6 · Measurement → age",
      sidebarLayout(
        sidebarPanel(width = 4,
          radioButtons("f_mode", "Input", choices = c("F14C → 14C age" = "f2a", "14C age → F14C" = "a2f")),
          conditionalPanel("input.f_mode == 'f2a'",
            numericInput("f14", "F14C", 0.737, min = 0.0001, step = 0.001),
            numericInput("f14s", "Measurement error of F14C (1σ)", 0.002, min = 0.00001, step = 0.0005)),
          conditionalPanel("input.f_mode == 'a2f'",
            numericInput("a_in", "14C age (BP)", 2450, step = 10),
            numericInput("a_err", "Error (1σ)", 25, min = 1, step = 5)),
          hr(),
          h4("Transfer to tab 1"),
          textInput("f_id", "Sample ID", "NEW1"),
          textInput("f_layer", "Layer", "A"),
          actionButton("f_add", "Add to sample table"),
          br(), br(),
          p(class = "hint", "The 14C age is calculated with the Libby value: age = −8033 · ln(F14C). The calibration also corrects for the difference from the true half-life.")
        ),
        mainPanel(width = 8,
          h4(textOutput("f14_text")),
          plotOutput("plot_f14", height = "320px"),
          textOutput("f14_hpd"),
          h4("For practice"),
          tableOutput("tab_f14")
        )
      )
    ),
    tabPanel("7 · Dendro crossdating",
      sidebarLayout(
        sidebarPanel(width = 4,
          p("A wood sample is slid against a (simulated) reference chronology. Where does the ring pattern fit? Then: how precise is the felling year?"),
          h4("Simulate sample"),
          numericInput("d_n", "Number of rings", 60, min = 20, max = 150, step = 5),
          sliderInput("d_noise", "Noise of the sample (site)", min = 0.05, max = 0.8, value = 0.25, step = 0.05),
          numericInput("d_end", "True end year (the solution), AD", 102, min = 20, max = 300),
          numericInput("d_seed", "Random number (same number = same sample)", 1, min = 1, step = 1),
          actionButton("d_sim", "Simulate sample"),
          hr(),
          sliderInput("d_pos", "End year of the sample in the reference (slider)", min = 20, max = 300, value = 180, step = 1, width = "100%"),
          sliderInput("d_thr", "t-value threshold (rule of thumb, lab-dependent)", min = 2, max = 6, value = 3.5, step = 0.1),
          checkboxInput("d_show_scan", "Show fit curve over all positions", FALSE),
          checkboxInput("d_show_true", "Show solution", FALSE),
          hr(),
          h4("Determine felling year"),
          radioButtons("d_mode", NULL, choices = c("Bark edge preserved" = "wk", "Sapwood, but no bark edge" = "splint", "heartwood only" = "kern")),
          selectInput("d_reg", "Region (sapwood range)", choices = names(DENDRO_SAPWOOD)),
          numericInput("d_sap", "Measured sapwood rings", 0, min = 0, max = 80),
          numericInput("d_margin", "Margin ± years (with bark edge)", 1, min = 0, max = 5),
          hr(),
          h4("Transfer to tab 1"),
          textInput("d_id", "Sample ID", "A1"),
          textInput("d_layer", "Layer", "A"),
          actionButton("d_add", "Add as dendro sample")
        ),
        mainPanel(width = 8,
          plotOutput("d_plot", height = "330px"),
          textOutput("d_stat"),
          conditionalPanel("input.d_show_scan", plotOutput("d_scan", height = "250px")),
          h4("Best positions"),
          tableOutput("d_tab"),
          h4("Felling year"),
          textOutput("d_fell_txt"),
          p(class = "hint", "The reference is simulated (schematic). t-value after Baillie-Pilcher, sign agreement (Glk) in %. Ring widths are detrended with a moving average (9 years) before comparison. Short samples quickly produce chance matches: try 25 against 100 rings. Oak sapwood ranges after Tegel et al. 2022.")
        )
      )
    ),
    tabPanel("Help & limitations",
      br(),
      h4("What the model does"),
      tags$ul(
        tags$li("Every layer has a deposition date D. Every sample dates an event S (death of the tree, death of the animal, felling year) that does not lie after D."),
        tags$li("Coins set a lower limit (TPQ) for D, written sources or sealing contexts an upper limit (TAQ)."),
        tags$li("Stratigraphy enforces D(lower layer) ≤ D(upper layer)."),
        tags$li("The calculation uses Monte Carlo (draw and reject), not MCMC. This is sufficient for the seminar, but the results are not identical to OxCal.")
      ),
      h4("Assumption switches"),
      p("With the switches in tab 1 every assumption can be switched on and off individually. The grey line shows the deposition from the ¹⁴C data only, without stratigraphy, TPQ, TAQ and dendro. This shows what each piece of additional information contributes."),
      h4("Year counting"),
      tags$ul(
        tags$li("Internally astronomical as in OxCal (year 0 = 1 BC, −78 = 79 BC). In the tables you can also enter values such as “480 BC” and “101 AD”; they are converted."),
        tags$li("The selection at the top right switches the display between BC/AD and astronomical counting. Fields such as “Latest conceivable deposition” expect astronomical years; the year converter helps.")
      ),
      h4("Dendro crossdating"),
      p("Tab 7 works with a simulated reference chronology. It shows the principle (detrending, correlation, t-value, sign agreement) but does not replace a real standard chronology."),
      h4("Not included"),
      tags$ul(
        tags$li("Outlier models (only in the OxCal export), reservoir effects, phases with boundary priors, Marine20/SHCal20."),
        tags$li("Invalid input is reported as a message. In case of conflicts (e.g. coin younger than TAQ) the model deliberately stops.")
      ),
      h4("Suggestion for comparison"),
      p("Model the same data in OxCal (Sequence, Date, Boundary, After/Before) and compare the results. Where do they differ, and why?")
    )
  )
)

# ---- Server -----------------------------------------------------------------
server <- function(input, output, session) {

  get_opt <- function() list(strat = isTRUE(input$a_strat), tpq = isTRUE(input$a_tpq),
                             taq = isTRUE(input$a_taq), dendro = isTRUE(input$a_dendro))

  # ---- Szenarien
  output$scen_desc <- renderUI({
    s <- SCEN[[input$scen]]
    if (is.null(s)) p(class = "hint", "Own data: edit the tables below or upload a CSV.")
    else p(class = "hint", s$desc)
  })
  observeEvent(input$load_scen, {
    s <- SCEN[[input$scen]]
    if (is.null(s)) return()
    updateTextAreaInput(session, "txt_smp", value = s$smp)
    updateTextAreaInput(session, "txt_lay", value = s$lay)
    updateNumericInput(session, "upper", value = s$upper)
    showNotification("Scenario loaded. Now press “Run model”.", type = "message", duration = 4)
  })

  # ---- Jahreswandler
  output$yw_out <- renderText({
    n <- input$yw_n
    if (is.null(n) || is.na(n)) return("")
    y <- if (identical(input$yw_era, "bc")) 1 - abs(n) else abs(n)
    paste0(n, if (identical(input$yw_era, "bc")) " BC" else " AD", "  →  ", y, " (astronomical)")
  })

  observeEvent(input$f_smp, {
    updateTextAreaInput(session, "txt_smp",
      value = paste(readLines(input$f_smp$datapath, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  })
  observeEvent(input$f_lay, {
    updateTextAreaInput(session, "txt_lay",
      value = paste(readLines(input$f_lay$datapath, warn = FALSE, encoding = "UTF-8"), collapse = "\n"))
  })

  # ---- Modell
  model <- eventReactive(input$run, {
    tryCatch({
      smp <- read_txt(input$txt_smp)
      lay <- read_txt(input$txt_lay)
      m <- run_model(smp, lay, upper = input$upper, opt = get_opt())
      m$base <- NULL
      if (isTRUE(input$cmp)) {
        m$base <- tryCatch(run_model(smp, lay, upper = input$upper, N = 20000, opt = OPT_NONE, yrs = m$yrs),
                           error = function(e) NULL)
      }
      m
    }, error = function(e) list(error = conditionMessage(e)))
  }, ignoreNULL = FALSE)

  nlayers <- reactive({
    lay <- tryCatch(read_txt(input$txt_lay), error = function(e) NULL)
    if (is.null(lay)) 1 else max(1, nrow(lay))
  })

  output$acc <- renderText({
    m <- model()
    if (!is.null(m$error)) return("")
    on <- c(if (m$opt$strat) "Stratigraphy", if (m$opt$tpq) "TPQ", if (m$opt$taq) "TAQ", if (m$opt$dendro) "Dendro")
    paste0("Active assumptions: ", if (length(on)) paste(on, collapse = ", ") else "none",
           "   |   Accepted draws (stratigraphy satisfied): ", round(100 * m$acc, 1), " %")
  })

  output$plot_model <- renderPlot({
    m <- model()
    validate(need(is.null(m$error), m$error))
    mode <- input$yrfmt
    nl <- nrow(m$lay)
    op <- par(mfrow = c(nl, 1), mar = c(3.4, 3, 2.2, 1), mgp = c(2, 0.7, 0))
    on.exit(par(op))
    for (k in seq_len(nl)) {
      mem <- m$members[[k]]
      has_base <- !is.null(m$base) && length(m$base$members[[k]]) > 0
      sc <- max(c(m$postD[, k], unlist(m$L[mem]), m$postS[, mem], if (has_base) m$base$postD[, k]))
      plot(NA, xlim = range(m$yrs), ylim = c(0, 1.08), yaxt = "n", xaxt = "n",
           xlab = xlab_years(mode), ylab = "rel. probability",
           main = paste0("Layer ", m$lay$schicht[k], " (rank ", m$lay$rang[k], ")"))
      axis_years(range(m$yrs), mode)
      if (has_base) lines(m$yrs, m$base$postD[, k] / sc, col = "grey55", lwd = 2.5)
      for (j in seq_along(mem)) {
        i <- mem[j]; col <- PAL[(i - 1) %% length(PAL) + 1]
        polygon(c(m$yrs, rev(m$yrs)), c(m$postS[, i] / sc, rep(0, length(m$yrs))),
                col = adjustcolor(col, 0.35), border = col)
        lines(m$yrs, m$L[[i]] / sc, col = col, lty = 2, lwd = 1.5)
      }
      lines(m$yrs, m$postD[, k] / sc, lwd = 3, col = "black")
      if (!is.na(m$lay$tpq[k])) abline(v = m$lay$tpq[k], col = "red", lwd = 2, lty = 3)
      if (!is.na(m$lay$taq[k])) abline(v = m$lay$taq[k], col = "blue", lwd = 2, lty = 3)
      leg <- c(as.character(m$smp$id[mem]), "Deposition (model)", if (has_base) "Deposition 14C only")
      cols <- c(PAL[(mem - 1) %% length(PAL) + 1], "black", if (has_base) "grey55")
      legend("topright", bty = "n", cex = 0.85, legend = leg, col = cols,
             lwd = c(rep(2, length(mem)), 3, if (has_base) 2.5))
    }
  }, height = function() max(300, 230 * nlayers()))

  results <- reactive({
    m <- model()
    validate(need(is.null(m$error), m$error))
    mode <- input$yrfmt
    s <- data.frame(
      Object = m$smp$id, Layer = m$smp$schicht,
      Type = ifelse(is_c14(m$smp$typ), "14C", "Dendro"),
      Unmodelled = vapply(seq_len(nrow(m$smp)), function(i) hpd(m$L[[i]], m$yrs, mode = mode), ""),
      Modelled = vapply(seq_len(nrow(m$smp)), function(i) hpd(m$postS[, i], m$yrs, mode = mode), ""),
      stringsAsFactors = FALSE)
    s[["14C only (deposition)"]] <- "–"
    d <- data.frame(
      Object = paste0("Deposition ", m$lay$schicht), Layer = m$lay$schicht, Type = "Model",
      Unmodelled = "–",
      Modelled = vapply(seq_len(nrow(m$lay)), function(k) hpd(m$postD[, k], m$yrs, mode = mode), ""),
      stringsAsFactors = FALSE)
    d[["14C only (deposition)"]] <- vapply(seq_len(nrow(m$lay)), function(k) {
      if (!is.null(m$base) && length(m$base$members[[k]]) > 0) hpd(m$base$postD[, k], m$yrs, mode = mode) else "–"
    }, "")
    rbind(s, d)
  })
  output$tab <- renderTable(results(), striped = TRUE, spacing = "s")
  output$dl <- downloadHandler(
    filename = function() "dating_results.csv",
    content = function(file) write.csv(results(), file, row.names = FALSE, fileEncoding = "UTF-8"))

  # ---- Altholz-Check
  off_data <- function(m) {
    if (identical(input$off_ref, "dep")) list(dens = m$off_dens, med = m$off_med, q95 = m$off_q95,
                                              xlab = "Years between sample event and deposition")
    else list(dens = m$rel_dens, med = m$rel_med, q95 = m$rel_q95,
              xlab = "Years older than the youngest sample of the layer")
  }
  output$plot_off <- renderPlot({
    m <- model()
    validate(need(is.null(m$error), m$error))
    od <- off_data(m)
    nl <- nrow(m$lay)
    xmax <- input$off_max
    op <- par(mfrow = c(nl, 1), mar = c(3.4, 3, 2.2, 1), mgp = c(2, 0.7, 0))
    on.exit(par(op))
    for (k in seq_len(nl)) {
      mem <- m$members[[k]]
      mem <- mem[!vapply(od$dens[mem], is.null, logical(1))]
      # cut off the peak at 0 (youngest sample) so that the other curves stay visible
      ymax <- max(c(0.001, unlist(lapply(mem, function(i) od$dens[[i]][6:(xmax + 1)]))))
      plot(NA, xlim = c(0, xmax), ylim = c(0, ymax * 1.15), yaxt = "n",
           xlab = od$xlab, ylab = "Probability",
           main = paste0("Layer ", m$lay$schicht[k], if (!length(mem)) " (only one sample: no comparison possible)" else ""))
      if (length(mem)) mtext("Peak at 0 cut off", side = 3, line = 0, adj = 1, cex = 0.7, col = "grey40")
      for (i in mem) {
        col <- PAL[(i - 1) %% length(PAL) + 1]
        lines(0:xmax, od$dens[[i]][1:(xmax + 1)], col = col, lwd = 2.5)
        abline(v = od$med[i], col = col, lty = 3)
      }
      if (length(mem)) legend("topright", bty = "n", cex = 0.85, legend = as.character(m$smp$id[mem]),
                              col = PAL[(mem - 1) %% length(PAL) + 1], lwd = 2.5)
    }
  }, height = function() max(300, 230 * nlayers()))

  output$tab_off <- renderTable({
    m <- model()
    validate(need(is.null(m$error), m$error))
    od <- off_data(m)
    ok <- !is.na(od$med)
    d <- data.frame(Sample = m$smp$id, Layer = m$smp$schicht,
                    `Median distance (years)` = ifelse(ok, as.character(round(od$med)), "–"),
                    `95% quantile (years)` = ifelse(ok, as.character(round(od$q95)), "–"),
                    Note = ifelse(!ok, "no comparison (only sample in the layer)",
                                     ifelse(od$med > input$off_thr, "flagged: old wood or reworked?", "not flagged")),
                    check.names = FALSE, stringsAsFactors = FALSE)
    d
  }, striped = TRUE, spacing = "s")

  # ---- Wiggle
  true_start <- reactiveVal(NA_real_)
  observeEvent(input$sim, {
    n <- max(3, input$w_n); step <- input$w_step; err <- input$w_err; ys <- input$w_true
    abst <- (0:(n - 1)) * step
    cv <- curve_at(ys + abst)
    age <- round(cv$mu + rnorm(n, 0, err))
    true_start(ys)
    updateTextAreaInput(session, "w_txt",
      value = paste(c("offset,c14,error", paste(abst, age, err, sep = ",")), collapse = "\n"))
  })

  wdat <- reactive({
    d <- tryCatch(read_txt(input$w_txt), error = function(e) NULL)
    validate(need(!is.null(d) && all(c("abstand", "c14", "fehler") %in% names(d)),
                  "Paste data with the columns offset, c14, error or start the simulation."))
    d <- d[stats::complete.cases(d[, c("abstand", "c14", "fehler")]), ]
    validate(need(nrow(d) >= 3, "At least 3 samples needed."))
    d
  })

  output$plot_wig <- renderPlot({
    d <- wdat(); st <- input$w_start; mode <- input$yrfmt
    yrs <- st + d$abstand
    lo <- st - 60; hi <- max(yrs) + 60
    cv <- CURVE[CURVE$year >= lo & CURVE$year <= hi, ]
    ylim <- range(c(cv$c14 - cv$err, cv$c14 + cv$err, d$c14 - d$fehler, d$c14 + d$fehler))
    plot(NA, xlim = c(lo, hi), ylim = ylim, xaxt = "n", xlab = xlab_years(mode),
         ylab = "14C age (BP)", main = "IntCal20 and shifted ring sequence")
    axis_years(c(lo, hi), mode)
    polygon(c(cv$year, rev(cv$year)), c(cv$c14 - cv$err, rev(cv$c14 + cv$err)),
            col = adjustcolor("grey60", 0.4), border = NA)
    lines(cv$year, cv$c14, lwd = 2, col = "grey30")
    arrows(yrs, d$c14 - d$fehler, yrs, d$c14 + d$fehler, angle = 90, code = 3, length = 0.04, col = "firebrick")
    points(yrs, d$c14, pch = 19, col = "firebrick")
    if (isTRUE(input$w_show_true) && !is.na(true_start())) abline(v = true_start(), lty = 2, col = "darkgreen")
  })

  output$wig_text <- renderText({
    d <- wdat()
    T0 <- wig_T(input$w_start, d)
    crit <- qchisq(0.95, df = nrow(d) - 1)
    paste0("T = ", round(T0, 1), "   (critical value 5%: ", round(crit, 1), ", df = ", nrow(d) - 1, ")  → ",
           if (T0 <= crit) "Fit statistically acceptable." else "Fit not acceptable.")
  })

  output$plot_T <- renderPlot({
    d <- wdat(); mode <- input$yrfmt
    starts <- seq(-1200, -100, by = 1)
    Tv <- vapply(starts, wig_T, numeric(1), d = d)
    crit <- qchisq(0.95, df = nrow(d) - 1)
    plot(starts, Tv, type = "l", lwd = 2, log = "y", xaxt = "n", xlab = xlab_years(mode), ylab = "T (log)",
         main = "Fit of the sequence as a function of start year")
    axis_years(range(starts), mode)
    abline(h = crit, col = "red", lty = 3)
    abline(v = input$w_start, col = "firebrick", lwd = 2)
    if (isTRUE(input$w_show_true) && !is.na(true_start())) abline(v = true_start(), lty = 2, col = "darkgreen")
  })

  # ---- OxCal-Export
  oxcode <- reactive({
    p <- tryCatch(prep_tables(read_txt(input$txt_smp), read_txt(input$txt_lay)),
                  error = function(e) list(error = conditionMessage(e)))
    validate(need(is.null(p$error), p$error))
    lo <- tryCatch(min(make_grid(p$smp, p$lay, input$upper)), error = function(e) NA_real_)
    if (is.na(lo)) lo <- min(c(p$lay$tpq, p$lay$taq, input$upper), na.rm = TRUE) - 3000
    ap <- apply_assumptions(p$smp, p$lay, get_opt())
    oxcal_code(ap$smp, ap$lay, input$upper, lo, outlier = isTRUE(input$ox_out))
  })
  output$oxcal_code <- renderText(paste(oxcode(), collapse = "\n"))
  output$dl_ox <- downloadHandler(
    filename = function() "model.oxcal",
    content = function(file) writeLines(oxcode(), file, useBytes = TRUE))

  # ---- Combination (same event)
  output$cmb_layer_ui <- renderUI({
    m <- model()
    if (!is.null(m$error)) return(p(class = "hint", "Please run a model in tab 1 first."))
    n14 <- vapply(seq_len(nrow(m$lay)), function(k) sum(is_c14(m$smp$typ[m$members[[k]]])), numeric(1))
    ch <- m$lay$schicht[n14 >= 2]
    if (!length(ch)) return(p(class = "hint", "No layer with at least two 14C samples."))
    selectInput("cmb_layer", "Layer (all samples = same event)", choices = ch)
  })
  cmb <- reactive({
    m <- model()
    validate(need(is.null(m$error), m$error))
    validate(need(!is.null(input$cmb_layer), "No layer with at least two 14C samples."))
    k <- which(m$lay$schicht == input$cmb_layer)
    idx <- m$members[[k]]
    idx <- idx[is_c14(m$smp$typ[idx])]
    validate(need(length(idx) >= 2, "At least two 14C samples needed in the layer."))
    ages <- m$smp$wert1[idx]; errs <- m$smp$wert2[idx]
    w <- 1 / errs^2
    mu <- sum(w * ages) / sum(w); se <- sqrt(1 / sum(w))
    chi <- sum(w * (ages - mu)^2); df <- length(idx) - 1
    list(m = m, idx = idx, mu = mu, se = se, chi = chi, df = df, crit = qchisq(0.95, df),
         dens = dens_c14(mu, se, m$yrs))
  })
  output$plot_cmb <- renderPlot({
    r <- cmb(); m <- r$m; mode <- input$yrfmt
    sc <- max(c(r$dens, unlist(m$L[r$idx])))
    plot(NA, xlim = range(m$yrs), ylim = c(0, 1.08), xaxt = "n", yaxt = "n", xlab = xlab_years(mode),
         ylab = "rel. probability", main = paste0("Combination in layer ", input$cmb_layer))
    axis_years(range(m$yrs), mode)
    for (j in seq_along(r$idx)) {
      i <- r$idx[j]; col <- PAL[(i - 1) %% length(PAL) + 1]
      lines(m$yrs, m$L[[i]] / sc, col = col, lty = 2, lwd = 1.8)
    }
    polygon(c(m$yrs, rev(m$yrs)), c(r$dens / sc, rep(0, length(m$yrs))), col = adjustcolor("black", 0.25), border = "black", lwd = 2.5)
    legend("topright", bty = "n", cex = 0.85, legend = c(as.character(m$smp$id[r$idx]), "Combination"),
           col = c(PAL[(r$idx - 1) %% length(PAL) + 1], "black"), lwd = c(rep(1.8, length(r$idx)), 2.5), lty = c(rep(2, length(r$idx)), 1))
  })
  output$cmb_head <- renderText({
    r <- cmb(); sprintf("Weighted mean: %.0f ± %.0f BP  (from %d measurements)", r$mu, r$se, length(r$idx))
  })
  output$cmb_chi <- renderText({
    r <- cmb()
    sprintf("χ² = %.1f with %d degrees of freedom (critical value 5%%: %.1f)  →  %s", r$chi, r$df, r$crit,
            if (r$chi <= r$crit) "Measurements agree: combination permissible." else "Measurements do NOT agree: not the same event, or outliers?")
  })
  output$cmb_hpd <- renderText({
    r <- cmb(); paste0("95.4% range of the combination: ", hpd(r$dens, r$m$yrs, mode = input$yrfmt))
  })

  # ---- SPD
  output$plot_spd <- renderPlot({
    m <- model()
    validate(need(is.null(m$error), m$error))
    mode <- input$yrfmt
    idx <- which(is_c14(m$smp$typ))
    validate(need(length(idx) > 0, "No 14C samples available."))
    spd_raw <- Reduce(`+`, m$L[idx])
    spd_raw <- spd_raw / max(spd_raw)
    plot(m$yrs, spd_raw, type = "l", lwd = 2, xaxt = "n", xlab = xlab_years(mode),
         ylab = "Sum (normalised)", main = paste0("Summed curve of ", length(idx), " 14C dates"), ylim = c(0, 1.05))
    axis_years(range(m$yrs), mode)
    if (isTRUE(input$spd_mod)) {
      spd_mod <- rowSums(m$postS[, idx, drop = FALSE])
      lines(m$yrs, spd_mod / max(spd_mod), lwd = 2, col = "firebrick")
      legend("topright", bty = "n", lwd = 2, col = c("black", "firebrick"), legend = c("unmodelled", "modelled"))
    }
  })

  # ---- Measurement -> age
  f14_calc <- reactive({
    if (identical(input$f_mode, "a2f")) {
      validate(need(is.finite(input$a_in) && is.finite(input$a_err) && input$a_err > 0, "Enter age and error (> 0)."))
      age <- input$a_in; err <- input$a_err
      f <- exp(-age / 8033); s <- f * err / 8033
    } else {
      validate(need(is.finite(input$f14) && input$f14 > 0, "F14C must be greater than 0."),
               need(is.finite(input$f14s) && input$f14s > 0, "The F14C error must be greater than 0."))
      f <- input$f14; s <- input$f14s
      age <- -8033 * log(f); err <- 8033 * s / f
    }
    list(f = f, s = s, age = age, err = err)
  })
  output$f14_text <- renderText({
    r <- f14_calc()
    sprintf("F14C = %.4f ± %.4f   →   14C age = %.0f ± %.0f BP", r$f, r$s, r$age, r$err)
  })
  output$plot_f14 <- renderPlot({
    r <- f14_calc()
    yrs <- tryCatch(age_grid(r$age, r$err), error = function(e) NULL)
    validate(need(!is.null(yrs), "This age lies outside the calibration curve (or is negative)."))
    mode <- input$yrfmt
    d <- dens_c14(r$age, r$err, yrs)
    plot(yrs, d / max(d), type = "n", xaxt = "n", yaxt = "n", xlab = xlab_years(mode), ylab = "rel. probability",
         main = "Calibrated date (IntCal20)", ylim = c(0, 1.05))
    axis_years(range(yrs), mode)
    polygon(c(yrs, rev(yrs)), c(d / max(d), rep(0, length(yrs))), col = adjustcolor("firebrick", 0.35), border = "firebrick")
  })
  output$f14_hpd <- renderText({
    r <- f14_calc()
    yrs <- tryCatch(age_grid(r$age, r$err), error = function(e) NULL)
    if (is.null(yrs)) return("")
    paste0("95.4% range: ", hpd(dens_c14(r$age, r$err, yrs), yrs, mode = input$yrfmt))
  })
  output$tab_f14 <- renderTable({
    f <- c(1, 0.9, 0.75, 0.5, 0.25, 0.1, 0.01)
    data.frame(F14C = f, `14C age (BP)` = round(-8033 * log(f)), check.names = FALSE)
  }, digits = 3, striped = TRUE, spacing = "s")
  observeEvent(input$f_add, {
    r <- f14_calc()
    id <- gsub("[, ]", "_", input$f_id); lay <- gsub("[, ]", "_", input$f_layer)
    line <- paste(id, lay, "14C", round(r$age), max(1, round(r$err)), sep = ",")
    cur <- sub("\\s+$", "", input$txt_smp)
    updateTextAreaInput(session, "txt_smp", value = paste0(cur, "\n", line))
    showNotification(paste0("Row added: ", line, ". The layer must appear in the layer table."), type = "message", duration = 5)
  })
  # ---- Dendro-Crossdating
  dsim <- reactiveVal(NULL)
  observeEvent(input$d_sim, {
    L <- length(DENDRO_REF)
    n <- as.integer(max(20, min(150, nz(input$d_n, 60))))
    end_true <- as.integer(min(L, max(n, round(nz(input$d_end, 102)))))
    x <- dendro_make_sample(DENDRO_REF, n, end_true, nz(input$d_noise, 0.25), nz(input$d_seed, 1))
    dsim(list(x = x, n = n, end_true = end_true))
    cand <- setdiff(n:L, (end_true - 25):(end_true + 25))
    if (length(cand)) updateSliderInput(session, "d_pos", value = cand[sample.int(length(cand), 1)])
  })

  dd <- reactive({
    s <- dsim()
    validate(need(!is.null(s), "Please press “Simulate sample”."))
    smpidx <- dendro_index(s$x)
    n <- length(s$x); L <- length(DENDRO_REF)
    ends <- n:L
    M <- t(vapply(ends, function(e) dendro_stats(DENDRO_REF_IDX, smpidx, e), numeric(4)))
    list(ends = ends, M = M, smpidx = smpidx, n = n, end_true = s$end_true)
  })
  d_pos_eff <- reactive({ d <- dd(); min(max(nz(input$d_pos, 100), d$n), length(DENDRO_REF)) })

  output$d_plot <- renderPlot({
    d <- dd(); n <- d$n; L <- length(DENDRO_REF); e <- d_pos_eff()
    z <- function(v) (v - mean(v, na.rm = TRUE)) / sd(v, na.rm = TRUE)
    zr <- z(DENDRO_REF_IDX); zs <- z(d$smpidx)
    xl <- c(max(1, e - n - 15), min(L, e + 15))
    yl <- range(c(zr[xl[1]:xl[2]], zs), na.rm = TRUE)
    plot(NA, xlim = xl, ylim = yl, xlab = "Reference year (AD)", ylab = "Ring-width index (standardised)",
         main = "Reference chronology (grey) and shifted sample (green)")
    lines(seq_len(L), zr, col = "grey55", lwd = 2)
    lines((e - n + 1):e, zs, col = "#1E8E5A", lwd = 2.5)
    abline(v = c(e - n + 1, e), lty = 3, col = "#1E8E5A")
    if (isTRUE(input$d_show_true)) abline(v = d$end_true, col = "red", lty = 2, lwd = 2)
    legend("topright", bty = "n", cex = 0.9, legend = c("Reference", "Sample", if (isTRUE(input$d_show_true)) "true end year"),
           col = c("grey55", "#1E8E5A", if (isTRUE(input$d_show_true)) "red"), lwd = c(2, 2.5, if (isTRUE(input$d_show_true)) 2),
           lty = c(1, 1, if (isTRUE(input$d_show_true)) 2))
  })

  output$d_stat <- renderText({
    d <- dd(); e <- d_pos_eff()
    st <- dendro_stats(DENDRO_REF_IDX, dendro_index(dsim()$x), e)
    if (is.na(st["t"])) return("No comparison possible at this position.")
    paste0("End year ", e, " AD:  r = ", round(st["r"], 2), ",  t = ", round(st["t"], 1),
           ",  Glk = ", round(st["glk"]), " %,  overlap ", st["n"], " rings  →  ",
           if (st["t"] >= input$d_thr) "above the threshold." else "below the threshold: not yet convincing.")
  })

  output$d_scan <- renderPlot({
    d <- dd(); tv <- d$M[, "t"]; e <- d_pos_eff()
    plot(d$ends, tv, type = "l", lwd = 2, xlab = "End year of the sample in the reference", ylab = "t-value",
         main = "Fit at all positions")
    abline(h = input$d_thr, col = "red", lty = 3)
    abline(v = e, col = "#1E8E5A", lwd = 2)
    if (any(is.finite(tv))) points(d$ends[which.max(tv)], max(tv, na.rm = TRUE), pch = 19, cex = 1.3)
    if (isTRUE(input$d_show_true)) abline(v = d$end_true, col = "red", lty = 2, lwd = 2)
  })

  output$d_tab <- renderTable({
    d <- dd()
    o <- order(d$M[, "t"], decreasing = TRUE, na.last = NA)
    o <- head(o, 5)
    out <- data.frame(`End year` = as.character(d$ends[o]), r = sprintf("%.2f", d$M[o, "r"]), `t-value` = sprintf("%.1f", d$M[o, "t"]),
                      `Glk (%)` = sprintf("%.0f", d$M[o, "glk"]), Overlap = sprintf("%.0f", d$M[o, "n"]),
                      check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[5] <- "Overlap"
    out
  }, striped = TRUE, spacing = "s")

  d_fell <- reactive({
    reg <- DENDRO_SAPWOOD[[nz(input$d_reg, names(DENDRO_SAPWOOD)[1])]]
    dendro_window(d_pos_eff(), nz(input$d_mode, "wk"), reg[1], reg[2], nz(input$d_sap, 0), nz(input$d_margin, 0))
  })
  output$d_fell_txt <- renderText({
    w <- d_fell(); mode <- input$yrfmt
    if (is.na(w[2])) paste0("Heartwood only: felling year no earlier than ", fmt_year(w[1], mode), " (TPQ only). Use as the TPQ of the layer in tab 1.")
    else if (w[1] == w[2]) paste0("Felling year: ", fmt_year(w[1], mode))
    else paste0("Felling-year window: ", fmt_year(w[1], mode), " to ", fmt_year(w[2], mode),
                "  (astronomical ", w[1], " to ", w[2], ")")
  })
  observeEvent(input$d_add, {
    w <- d_fell()
    if (is.na(w[2])) {
      showNotification("Heartwood alone yields only a TPQ, not a felling-year window. Please enter it as the TPQ in the layer table.", type = "warning", duration = 7)
      return()
    }
    id <- gsub("[, ]", "_", nz(input$d_id, "D1")); lay <- gsub("[, ]", "_", nz(input$d_layer, "A"))
    line <- paste(id, lay, "Dendro", w[1], w[2], sep = ",")
    cur <- sub("\\s+$", "", input$txt_smp)
    updateTextAreaInput(session, "txt_smp", value = paste0(cur, "\n", line))
    showNotification(paste0("Row added: ", line), type = "message", duration = 5)
  })

}

shinyApp(ui, server)
