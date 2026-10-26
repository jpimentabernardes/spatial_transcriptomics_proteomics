## =============================================================================
## R/uc_helpers.R
##
## Shared by the UC biopsy scripts (05_uc_insights.R, 05_uc_genes_eosinophils.R):
##   - Nancy score per core and low / high groups
##   - de-identified clinical table joined to the cores
##   - cell classes from the 04 labels (neutrophil, eosinophil, plasma, goblet, ...)
##   - per-core spatial features (densities, intraepithelial neutrophils, ...)
##   - statistics (core level vs Nancy, patient level vs outcomes, weighted kappa)
##   - core images on a black background
## Sourced after config/config.R and R/helpers.R.
## =============================================================================

## -----------------------------------------------------------------------
## Settings (copy into config/config.R to change them)
## -----------------------------------------------------------------------
UC_CLINICAL_CSV  <- Sys.getenv("EXP1_UC_CLINICAL", get0("UC_CLINICAL_CSV", ifnotfound = file.path(CONFIG_DIR, "clinical_uc.csv")))
UC_NANCY_COL     <- get0("UC_NANCY_COL", ifnotfound = NULL)       # NULL = Nancy_TMA, else Nancy (tma_map columns)
UC_NANCY_HIGH_MIN <- get0("UC_NANCY_HIGH_MIN", ifnotfound = 2)    # Nancy >= 2 = "high" (active), 0-1 = "low"
UC_LABEL_COL     <- get0("UC_LABEL_COL", ifnotfound = "cell_type")
UC_MIN_DENOM     <- get0("UC_MIN_DENOM", ifnotfound = 20)         # a % per core needs this many cells below the line
## clinical columns: outcomes compared per patient, Nancy columns compared with Nancy_TMA,
## eosinophil columns compared with the eosinophils counted on the TMA
UC_OUTCOME_COLS  <- get0("UC_OUTCOME_COLS", ifnotfound = c("Disease.Control", "HistoEndoscopic.Remission",
                         "Endoscopic.Remission", "Histological.Remission", "Biochemical.Remission", "Clinical.Remission"))
UC_CLIN_NANCY_COLS <- get0("UC_CLIN_NANCY_COLS", ifnotfound = c("NB.Nancy.Biopsy..TMA..single.Sigma.Score.",
                           "Nancy.histology.index", "Nancy.histology.index..calculated.", "Nancy_Pathology"))
UC_CLIN_EOS_COLS <- get0("UC_CLIN_EOS_COLS", ifnotfound = c("Presence.of.eosinophiles.in.histology.", "Eosinophiles"))
UC_CLIN_OTHER_COLS <- get0("UC_CLIN_OTHER_COLS", ifnotfound = c("Endo_mayo", "Endoscopic.Mayo.score", "Age", "Sex",
                           "Disease_Duration", "UC_Montreal", "Crp.mg.l.", "Current_Biological"))
## direct identifiers: never copied into any output
UC_IDENTIFIER_COLS <- c("Birthday", "Fallnummer.ORBIS", "Histology.journal.number", "Patient.ID..Mishra.",
                        "Sampling.Date", "Date.of.endoscopy", "Date.of.blood.withdrawal", "Samplingdate",
                        "Patient.ID.Finn.Sequencing..extended.to.full.cohort.", "BarcodeRNAseq", "Subject")

## cell class from the cell type label: first matching pattern wins (lower case,
## regular expressions). Edit in config.R if your labels are named differently.
UC_CLASS_PATTERNS <- get0("UC_CLASS_PATTERNS", ifnotfound = c(
  Eosinophil = "eosin", Neutrophil = "neutro", Mast = "mast", Plasma = "plasma",
  Inflam_mono_mac = "inflammatory_mono|inflam.*(mono|macro)", IFN_macrophage = "ifn",
  Macrophage = "macro|mono", DC = "dc$|dc[0-9_ ]|dendr",
  T_cell = "t_cell|cd4|cd8|treg|tfh|gd_t|naive|th17|_t$", NK_ILC = "^nk|ilc", B_cell = "b_cell|^b$",
  Inflamed_epi = "duox2|inflamed", Regenerative = "regenerat|^reg", Paneth = "paneth", Pyloric = "pyloric|gastric",
  Stem_TA = "stem|transit|^ta[_ ]|^ta$", Goblet = "goblet", BEST4 = "best4", EEC = "enteroendo|^eec", Tuft = "tuft",
  Colonocyte = "colonocyte|enterocyte",
  Inflam_fibroblast = "inflammatory_fib|inflam.*fib", Smooth_muscle = "muscle|myofib", Fibroblast = "fibro",
  Endothelial = "endoth|hev|lymphatic"))
EPI_LINEAGES    <- c("Epithelial", "Neural_NE")
IMMUNE_LINEAGES <- "Immune"
STROMA_LINEAGES <- c("Mesenchymal", "Endothelial", "Fibroblast", "Smooth_muscle_pericyte")

## -----------------------------------------------------------------------
## Cells, labels, classes
## -----------------------------------------------------------------------
uc_cell_table <- function(xen, label_col = UC_LABEL_COL) {
  m <- xen[[]]
  if (!label_col %in% colnames(m)) stop("Label column '", label_col, "' not in the object -- run 04 first ",
                                        "(or set UC_LABEL_COL in config.R).")
  m <- m[!is.na(m$core_id), ]
  m$label <- as.character(m[[label_col]])
  m$label[is.na(m$label) | !nzchar(m$label)] <- "unlabelled"
  ## the same name in two lineages (e.g. "Proliferating" in Immune and in Stroma) -> "Proliferating (Immune)"
  lin <- ifelse(is.na(m$lineage), "NA", as.character(m$lineage))
  multi <- tapply(lin, m$label, function(v) length(unique(v)) > 1)
  multi <- names(multi)[multi]
  if (length(multi)) {
    i <- m$label %in% multi
    m$label[i] <- paste0(m$label[i], " (", lin[i], ")")
    message("Labels used in more than one lineage, split by lineage: ", paste(multi, collapse = ", "))
  }
  m$compartment <- dplyr::case_when(m$lineage %in% EPI_LINEAGES ~ "Epithelium",
                                    m$lineage %in% IMMUNE_LINEAGES ~ "Immune",
                                    m$lineage %in% STROMA_LINEAGES ~ "Stroma", TRUE ~ "Other")
  m$class <- uc_classes(m$label)
  if (!"patient_id" %in% colnames(m)) m$patient_id <- m$core_id
  m$patient_id <- as.character(m$patient_id)
  m$patient_id[is.na(m$patient_id) | !nzchar(m$patient_id)] <- m$core_id[is.na(m$patient_id) | !nzchar(m$patient_id)]
  m
}

uc_classes <- function(label, patterns = UC_CLASS_PATTERNS) {
  lab <- tolower(label)
  out <- rep(NA_character_, length(lab))
  for (k in names(patterns)) {
    hit <- is.na(out) & grepl(patterns[[k]], lab, perl = TRUE)
    out[hit] <- k
  }
  out
}

## log-normalised expression (counts / cell total x median total), as in 04
uc_lognorm <- function(cnt, total) log1p(sweep(as.matrix(cnt), 2, pmax(total, 1), "/") * stats::median(total))

## -----------------------------------------------------------------------
## Nancy per core
## -----------------------------------------------------------------------
uc_nancy_col <- function(m) {
  if (!is.null(UC_NANCY_COL)) return(UC_NANCY_COL)
  for (c in c("Nancy_TMA", "Nancy")) if (c %in% colnames(m) && any(!is.na(suppressWarnings(as.numeric(m[[c]]))))) return(c)
  NULL
}
as_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))
nancy_group <- function(n, high_min = UC_NANCY_HIGH_MIN) {
  lv <- c(sprintf("low (Nancy < %d)", high_min), sprintf("high (Nancy >= %d)", high_min))
  factor(ifelse(is.na(n), NA, ifelse(n >= high_min, lv[2], lv[1])), levels = lv)
}
NANCY_COLOURS <- c("0" = "#2C7BB6", "1" = "#8FC1DE", "2" = "#F5C26B", "3" = "#E8743B", "4" = "#B2182B")
GROUP_COLOURS <- function(high_min = UC_NANCY_HIGH_MIN)
  stats::setNames(c("#3B6FB6", "#C0392B"), levels(nancy_group(0, high_min)))

## one row per core: slide, patient, Nancy (+ the other tma_map Nancy column), n cells
uc_core_table <- function(m) {
  nc <- uc_nancy_col(m)
  keep <- intersect(c("slide_id", "core_id", "patient_id", "organ", "diagnosis", "Nancy", "Nancy_TMA"), colnames(m))
  ct <- m[!duplicated(m$core_id), keep, drop = FALSE]
  ct$n_cells <- as.integer(table(m$core_id)[ct$core_id])
  ct$nancy <- if (is.null(nc)) NA_real_ else as_num(ct[[nc]])
  ct$nancy_group <- nancy_group(ct$nancy)
  attr(ct, "nancy_col") <- nc %||% "none"
  rownames(ct) <- ct$core_id
  ct
}

## -----------------------------------------------------------------------
## Clinical table: de-identified columns only, joined to the cores
## -----------------------------------------------------------------------
core_key <- function(sector, coord) {
  s <- tolower(gsub("[^[:alnum:]]", "", sub("(?i)^(sector|tma|panel)", "", as.character(sector), perl = TRUE)))
  co <- tolower(gsub("[^[:alnum:]]", "", as.character(coord)))
  paste0(s, "|", gsub("[^0-9]", "", co), "|", gsub("[^a-z]", "", co))
}
core_key_from_id <- function(core_id) {
  id <- as.character(core_id)
  has_sep <- grepl("_", id)
  sector <- ifelse(has_sep, sub("_.*$", "", id), sub("^([A-Za-z]+).*$", "\\1", id))
  coord <- ifelse(has_sep, sub("^[^_]*_", "", id), sub("^[A-Za-z]+", "", id))
  core_key(sector, coord)
}

read_uc_clinical <- function(path = UC_CLINICAL_CSV) {
  if (is.null(path) || !nzchar(path) || !file.exists(path)) {
    message("No clinical table (", path, ") -- clinical comparisons (Disease.Control, ...) are skipped. ",
            "Set UC_CLINICAL_CSV in config.R or EXP1_UC_CLINICAL.")
    return(NULL)
  }
  first <- readLines(path, n = 1, warn = FALSE)
  sep <- c("\t", ";", ",")[which.max(c(lengths(regmatches(first, gregexpr("\t", first))),
                                       lengths(regmatches(first, gregexpr(";", first))),
                                       lengths(regmatches(first, gregexpr(",", first)))))]
  cl <- read_config_csv(path, sep = sep, stringsAsFactors = FALSE, check.names = TRUE, fileEncoding = "latin1",
                        na.strings = c("", "NA", "na", "n.a.", "-", "#N/A"))
  ids <- intersect(UC_IDENTIFIER_COLS, colnames(cl))
  if (length(ids)) message("NOTE: the clinical table still contains direct identifiers (", paste(ids, collapse = ", "),
                           "). They are not used or written anywhere, but please keep a de-identified copy ",
                           "outside the code folder.")
  cl
}

## the clinical row of each core: by core_id, by TMA.Sector + Koordinate, or by patient_id
join_uc_clinical <- function(core_tab, cl) {
  if (is.null(cl)) return(NULL)
  want <- unique(c(UC_OUTCOME_COLS, UC_CLIN_NANCY_COLS, UC_CLIN_EOS_COLS, UC_CLIN_OTHER_COLS))
  want <- intersect(want, colnames(cl))
  how <- NULL
  if ("core_id" %in% colnames(cl)) {
    k_cl <- core_key_from_id(cl$core_id); how <- "core_id"
  } else if (all(c("TMA.Sector", "Koordinate") %in% colnames(cl))) {
    k_cl <- core_key(cl$TMA.Sector, cl$Koordinate); how <- "TMA.Sector + Koordinate"
  } else if ("patient_id" %in% colnames(cl)) {
    k_cl <- as.character(cl$patient_id); how <- "patient_id"
  } else {
    message("Clinical table: no core_id, TMA.Sector + Koordinate or patient_id column to join on -- skipped.")
    return(NULL)
  }
  k_core <- if (how == "patient_id") as.character(core_tab$patient_id) else core_key_from_id(core_tab$core_id)
  dup <- duplicated(k_cl) & !is.na(k_cl)
  if (any(dup)) message("Clinical table: ", sum(dup), " rows share a ", how, " key with an earlier row -- first kept.")
  idx <- match(k_core, k_cl[!dup])
  out <- cl[!dup, want, drop = FALSE][idx, , drop = FALSE]
  out <- cbind(core_id = core_tab$core_id, out)
  rownames(out) <- core_tab$core_id
  n_ok <- sum(!is.na(idx))
  message("Clinical table joined by ", how, ": ", n_ok, " / ", nrow(core_tab), " cores matched.")
  if (n_ok < nrow(core_tab)) {
    message("  cores without a clinical row (first 10): ", paste(utils::head(core_tab$core_id[is.na(idx)], 10), collapse = ", "),
            "\n  keys looked for: ", paste(utils::head(k_core[is.na(idx)], 5), collapse = ", "),
            "\n  keys in the clinical table: ", paste(utils::head(unique(k_cl), 5), collapse = ", "),
            "\n  (add a core_id column to the clinical table if the formats differ)")
  }
  if (n_ok == 0) return(NULL)
  out
}

## outcome as a factor: yes/no, 1/0, TRUE/FALSE -> "no" / "yes"; anything else as written
as_outcome <- function(x) {
  v <- tolower(trimws(as.character(x)))
  v[v %in% c("", "na", "nan")] <- NA
  yes <- c("1", "yes", "ja", "true", "y", "achieved", "remission")
  no  <- c("0", "no", "nein", "false", "n", "not achieved", "no remission")
  if (all(v[!is.na(v)] %in% c(yes, no))) return(factor(ifelse(v %in% yes, "yes", ifelse(v %in% no, "no", NA)),
                                                       levels = c("no", "yes")))
  factor(x)
}

## -----------------------------------------------------------------------
## Per-core spatial features
## -----------------------------------------------------------------------
## core area from the cells: occupied 20 um bins x 400 um^2 (holes in torn cores stay out)
core_area_mm2 <- function(x, y, bin = 20) {
  length(unique(paste(floor(x / bin), floor(y / bin)))) * bin^2 / 1e6
}
pct <- function(num, den, min_den = UC_MIN_DENOM) ifelse(den >= min_den, 100 * num / den, NA_real_)

## m: uc_cell_table() rows; returns one row per core
uc_core_features <- function(m, ien_k = 10, agg_radius = 20, agg_min = 3, base_um = 50) {
  classes_imm <- c("Neutrophil", "Eosinophil", "Mast", "Plasma", "B_cell", "T_cell", "NK_ILC",
                   "Macrophage", "Inflam_mono_mac", "IFN_macrophage", "DC")
  classes_epi <- c("Stem_TA", "Colonocyte", "Goblet", "BEST4", "EEC", "Tuft", "Inflamed_epi", "Regenerative",
                   "Paneth", "Pyloric")
  classes_str <- c("Fibroblast", "Inflam_fibroblast", "Smooth_muscle", "Endothelial")
  rows <- lapply(split(seq_len(nrow(m)), m$core_id), function(i) {
    d <- m[i, ]
    area <- core_area_mm2(d$x_um, d$y_um)
    n_epi <- sum(d$compartment == "Epithelium"); n_imm <- sum(d$compartment == "Immune")
    n_str <- sum(d$compartment == "Stroma")
    r <- list(core_id = d$core_id[1], n_cells = nrow(d), area_mm2 = round(area, 4),
              cells_per_mm2 = nrow(d) / area,
              Epithelium_pct = pct(n_epi, nrow(d)), Immune_pct = pct(n_imm, nrow(d)), Stroma_pct = pct(n_str, nrow(d)),
              Immune_per_mm2 = n_imm / area)
    cnt <- table(factor(d$class, levels = c(classes_imm, classes_epi, classes_str)))
    for (k in classes_imm) {
      r[[paste0(k, "_pct_immune")]] <- pct(cnt[[k]], n_imm)
      r[[paste0(k, "_per_mm2")]] <- cnt[[k]] / area
    }
    for (k in classes_epi) r[[paste0(k, "_pct_epithelium")]] <- pct(cnt[[k]], n_epi)
    for (k in classes_str) r[[paste0(k, "_pct_stroma")]] <- pct(cnt[[k]], n_str)
    r$Goblet_ratio <- pct(cnt[["Goblet"]], cnt[["Goblet"]] + cnt[["Colonocyte"]])
    ## neighbourhood features
    xy <- as.matrix(d[, c("x_um", "y_um")])
    if (nrow(d) > ien_k + 1) {
      nn <- RANN::nn2(xy, k = ien_k + 1)$nn.idx[, -1, drop = FALSE]
      epi_frac <- rowMeans(matrix(d$compartment[nn] == "Epithelium", nrow = nrow(d)))
      for (k in c("Neutrophil", "Eosinophil", "T_cell")) {
        s <- which(d$class %in% k)
        r[[paste0(k, "_intraepithelial_pct")]] <- pct(sum(epi_frac[s] >= 0.5), length(s), min_den = 5)
      }
    }
    s <- which(d$class %in% "Neutrophil")
    r$Neutrophil_aggregated_pct <- if (length(s) >= 5) {
      nb <- RANN::nn2(xy[s, , drop = FALSE], k = min(agg_min + 1, length(s)), searchtype = "radius", radius = agg_radius)$nn.idx
      100 * mean(rowSums(nb > 0) - 1 >= agg_min)
    } else NA_real_
    s <- which(d$class %in% "Plasma"); b <- which(d$class %in% "Stem_TA")
    r$Plasma_near_crypt_base_pct <- if (length(s) >= 5 && length(b) >= 5) {
      100 * mean(RANN::nn2(xy[b, , drop = FALSE], xy[s, , drop = FALSE], k = 1)$nn.dists[, 1] <= base_um)
    } else NA_real_
    as.data.frame(r, check.names = FALSE)
  })
  out <- dplyr::bind_rows(rows)
  ## classes never present in the labels: drop their all-NA / all-zero columns
  keep <- vapply(out, function(v) !is.numeric(v) || any(!is.na(v) & v != 0), logical(1))
  out[, keep, drop = FALSE]
}

## -----------------------------------------------------------------------
## Statistics
## -----------------------------------------------------------------------
## core level: each feature vs Nancy. Spearman over cores, Wilcoxon high vs low, and
## (when patients have several cores) a mixed model with the patient as random effect.
assoc_nancy <- function(df, cols, patient = "patient_id") {
  have_lme4 <- requireNamespace("lme4", quietly = TRUE)
  res <- dplyr::bind_rows(lapply(cols, function(v) {
    d <- df[!is.na(df[[v]]) & !is.na(df$nancy), ]
    if (nrow(d) < 6 || stats::sd(d[[v]]) == 0) return(NULL)
    hi <- d$nancy >= UC_NANCY_HIGH_MIN
    sp <- suppressWarnings(stats::cor.test(d[[v]], d$nancy, method = "spearman", exact = FALSE))
    wx <- if (sum(hi) >= 3 && sum(!hi) >= 3) suppressWarnings(stats::wilcox.test(d[[v]][hi], d[[v]][!hi], exact = FALSE))$p.value else NA
    p_mm <- NA
    if (have_lme4 && any(duplicated(d[[patient]])) && nrow(d) >= 8) {
      dd <- data.frame(y = d[[v]], nancy = d$nancy, pat = d[[patient]])
      p_mm <- tryCatch({
        m1 <- suppressMessages(lme4::lmer(y ~ nancy + (1 | pat), data = dd, REML = FALSE))
        m0 <- suppressMessages(lme4::lmer(y ~ 1 + (1 | pat), data = dd, REML = FALSE))
        stats::anova(m0, m1)[["Pr(>Chisq)"]][2]
      }, error = function(e) NA, warning = function(w) NA)
    }
    data.frame(feature = v, n_cores = nrow(d), n_patients = length(unique(d[[patient]])),
               spearman_rho = unname(sp$estimate), p_spearman = sp$p.value,
               median_low = stats::median(d[[v]][!hi]), median_high = stats::median(d[[v]][hi]),
               p_wilcox_high_vs_low = wx, p_mixed_model = p_mm)
  }))
  if (!nrow(res)) return(res)
  res$padj_spearman <- stats::p.adjust(res$p_spearman, "BH")
  res$padj_mixed_model <- stats::p.adjust(res$p_mixed_model, "BH")
  res[order(res$p_spearman), ]
}

## AUC = P(value of a "yes" > value of a "no") (Mann-Whitney), 0.5 = no difference
auc_mw <- function(x, y) {
  x <- x[!is.na(x)]; y <- y[!is.na(y)]
  if (!length(x) || !length(y)) return(NA_real_)
  r <- rank(c(x, y))
  (sum(r[seq_along(x)]) - length(x) * (length(x) + 1) / 2) / (length(x) * length(y))
}

## patient level: each feature vs one outcome (binary: AUC + Wilcoxon; more levels: Kruskal-Wallis)
assoc_outcome <- function(pat, cols, outcome) {
  y <- pat[[outcome]]
  res <- dplyr::bind_rows(lapply(cols, function(v) {
    ok <- !is.na(pat[[v]]) & !is.na(y)
    x <- pat[[v]][ok]; g <- droplevels(y[ok])
    if (nlevels(g) < 2 || min(table(g)) < 3) return(NULL)
    meds <- tapply(x, g, stats::median)
    out <- data.frame(outcome = outcome, feature = v, n_patients = length(x),
                      levels = paste(sprintf("%s (n=%d)", levels(g), as.integer(table(g))), collapse = " / "),
                      medians = paste(signif(meds, 3), collapse = " / "))
    if (nlevels(g) == 2) {
      out$auc <- auc_mw(x[g == levels(g)[2]], x[g == levels(g)[1]])
      out$p <- suppressWarnings(stats::wilcox.test(x ~ g, exact = FALSE))$p.value
    } else {
      out$auc <- NA; out$p <- stats::kruskal.test(x, g)$p.value
    }
    out
  }))
  if (nrow(res)) res$padj <- stats::p.adjust(res$p, "BH")
  res
}

## weighted (quadratic) Cohen's kappa between two ordinal scores
wkappa <- function(a, b) {
  ok <- !is.na(a) & !is.na(b)
  a <- a[ok]; b <- b[ok]
  if (length(a) < 3) return(NA_real_)
  lv <- sort(unique(c(a, b)))
  if (length(lv) < 2) return(NA_real_)
  O <- table(factor(a, lv), factor(b, lv)); O <- O / sum(O)
  E <- outer(rowSums(O), colSums(O))
  W <- outer(seq_along(lv), seq_along(lv), function(i, j) (i - j)^2) / (length(lv) - 1)^2
  1 - sum(W * O) / sum(W * E)
}

## -----------------------------------------------------------------------
## Core images on a black background
## -----------------------------------------------------------------------
theme_black <- function(base_size = 11) {
  theme_void(base_size = base_size) +
    theme(plot.background = element_rect(fill = "black", colour = NA),
          panel.background = element_rect(fill = "black", colour = NA),
          legend.background = element_rect(fill = "black", colour = NA),
          legend.key = element_rect(fill = "black", colour = NA),
          text = element_text(colour = "white"), legend.text = element_text(colour = "white"),
          plot.title = element_text(colour = "white", face = "bold", hjust = 0.5),
          plot.subtitle = element_text(colour = "grey80", hjust = 0.5),
          strip.text = element_text(colour = "white", face = "bold"),
          plot.margin = margin(6, 6, 6, 6))
}
## the same for charts (axes and text in light grey on black)
theme_black_chart <- function(base_size = 11) {
  theme_minimal(base_size = base_size) +
    theme(plot.background = element_rect(fill = "black", colour = NA),
          panel.background = element_rect(fill = "black", colour = NA),
          panel.grid.major = element_line(colour = "grey25", linewidth = 0.3), panel.grid.minor = element_blank(),
          text = element_text(colour = "white"), axis.text = element_text(colour = "grey85"),
          legend.background = element_rect(fill = "black", colour = NA),
          strip.text = element_text(colour = "white", face = "bold"),
          plot.title = element_text(colour = "white", face = "bold"))
}
## cell-type colours that read on black: the 04 palette, lifted towards white
lighten <- function(cols, f = 0.3) {
  rgb <- grDevices::col2rgb(cols) / 255
  stats::setNames(grDevices::rgb(t(rgb * (1 - f) + f)), names(cols))
}
## dot size (mm) so that each dot is about as wide as its cell (from cell_area), for a
## core of extent_um shown in panel_in inches
cell_dot_sizes <- function(area, extent_um, panel_in) {
  diam <- ifelse(is.na(area) | area <= 0, 10, 2 * sqrt(area / pi))
  pmax(0.15, pmin(diam, 25) * panel_in * 25.4 / extent_um / 0.75)
}
## one core: cells coloured by `colour_col` (palette pal); highlight = logical, drawn last
## with a ring; dim = logical, drawn first in dark grey
plot_core_black <- function(d, colour_col, pal, title = "", subtitle = NULL, panel_in = 5,
                            highlight = NULL, highlight_col = "#FF3EA5", dim = NULL, legend = TRUE) {
  ext <- max(diff(range(d$x_um)), diff(range(d$y_um)), 50)
  d$size <- cell_dot_sizes(d$cell_area %||% NA, ext, panel_in)
  d$col <- as.character(d[[colour_col]])
  p <- ggplot(mapping = aes(x_um, y_um))
  if (!is.null(dim) && any(dim))
    p <- p + geom_point(data = d[dim, ], aes(size = size), colour = "grey28", stroke = 0)
  front <- if (is.null(dim)) rep(TRUE, nrow(d)) else !dim
  p <- p + geom_point(data = d[front, ], aes(size = size, colour = col), stroke = 0) +
    scale_colour_manual(values = pal, name = NULL, drop = TRUE)
  if (!is.null(highlight) && any(highlight))
    p <- p + geom_point(data = d[highlight, ], aes(size = size * 1.9), shape = 21, colour = highlight_col,
                        fill = NA, stroke = 0.5)
  p <- p + scale_size_identity() + scale_y_reverse() + coord_fixed() + theme_black() +
    guides(colour = guide_legend(override.aes = list(size = 3.5), ncol = 1)) +
    labs(title = title, subtitle = subtitle) +
    ## 100 um scale bar, bottom left
    annotate("segment", x = min(d$x_um), xend = min(d$x_um) + 100, y = max(d$y_um) + 0.03 * ext,
             yend = max(d$y_um) + 0.03 * ext, colour = "white", linewidth = 0.8) +
    annotate("text", x = min(d$x_um) + 50, y = max(d$y_um) + 0.07 * ext, label = "100 um",
             colour = "white", size = 3)
  if (!legend) p <- p + theme(legend.position = "none")
  p
}
## colour by expression (counts) of one gene: zero cells dark grey, positive cells magma
plot_core_gene_black <- function(d, counts, gene, title = "", panel_in = 4) {
  ext <- max(diff(range(d$x_um)), diff(range(d$y_um)), 50)
  d$size <- cell_dot_sizes(d$cell_area %||% NA, ext, panel_in)
  d$expr <- counts
  pos <- d$expr > 0
  ggplot(mapping = aes(x_um, y_um)) +
    geom_point(data = d[!pos, ], aes(size = size), colour = "grey22", stroke = 0) +
    geom_point(data = d[pos, ][order(d$expr[pos]), ], aes(size = size * 1.15, colour = pmin(expr, 6)), stroke = 0) +
    scale_colour_viridis_c(option = "magma", begin = 0.35, end = 1, limits = c(1, 6), breaks = c(1, 2, 4, 6),
                           labels = c("1", "2", "4", "6+"), name = paste(gene, "\ncounts")) +
    scale_size_identity() + scale_y_reverse() + coord_fixed() + theme_black(10) +
    labs(title = title, subtitle = sprintf("%d / %d cells positive", sum(pos), nrow(d)))
}

## the legend of a ggplot as a stand-alone patchwork element on black (NULL if none)
legend_black <- function(p) {
  g <- ggplotGrob(p)
  i <- which(grepl("guide-box", vapply(g$grobs, function(x) x$name %||% "", character(1))))
  i <- i[vapply(g$grobs[i], function(x) !inherits(x, "zeroGrob"), logical(1))]
  if (!length(i)) return(NULL)
  patchwork::wrap_elements(full = g$grobs[[i[1]]]) +
    theme(plot.background = element_rect(fill = "black", colour = NA))
}

## how "nice" each core is for a figure: round (hull circularity), solid (tissue fills the
## hull), with both epithelium (>= 35%) and immune cells (>= 15%), and enough cells
core_shape_scores <- function(m) {
  out <- dplyr::bind_rows(lapply(split(seq_len(nrow(m)), m$core_id), function(i) {
    xy <- as.matrix(m[i, c("x_um", "y_um")])
    h <- xy[grDevices::chull(xy), , drop = FALSE]
    hx <- c(h[, 1], h[1, 1]); hy <- c(h[, 2], h[1, 2])
    A <- abs(sum(hx[-1] * hy[-length(hy)] - hx[-length(hx)] * hy[-1])) / 2
    Pm <- sum(sqrt(diff(hx)^2 + diff(hy)^2))
    data.frame(core_id = m$core_id[i[1]], n_cells = length(i), roundness = 4 * pi * A / max(Pm^2, 1),
               solidity = min(1, core_area_mm2(xy[, 1], xy[, 2]) * 1e6 / max(A, 1)),
               epithelium_pct = 100 * mean(m$compartment[i] == "Epithelium"),
               immune_pct = 100 * mean(m$compartment[i] == "Immune"))
  }))
  out$score <- with(out, roundness * solidity * pmin(epithelium_pct / 35, 1) * pmin(immune_pct / 15, 1) *
                      pmin(n_cells / 1500, 1))
  out
}
