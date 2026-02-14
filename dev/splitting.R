# ============================================================
# Unified ER splitting experiment for 4 datasets:
#   1) CORA (package cora) [gold pairs -> entity clusters]
#   2) Affiliation (csv ids + mapping entity_id)
#   3) 10K synthetic (full + duplicates pairs OR mapping)
#   4) NC voter (records + mapping OR gold pairs)
#
# Outputs:
#   res_all   : per-run results for each dataset x scheme
#   final_all : mean/sd summary by dataset x scheme
# ============================================================

suppressPackageStartupMessages({
  pkgs <- c("data.table", "stringdist", "igraph", "mclust", "aricode", "glmnet")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(data.table)
  library(stringdist)
  library(igraph)
  library(mclust)
  library(aricode)
  library(glmnet)
})

# -----------------------------
# CONFIG: EDIT ONLY THIS BLOCK
# -----------------------------
CFG <- list(
  # 1) CORA (no paths needed)
  cora = list(enabled = TRUE),

  # 2) Affiliation (you already have these)
  affiliation = list(
    enabled = TRUE,
    records_path = "D:/datasplitting/data/affiliationstrings_ids.csv",
    mapping_path = "D:/datasplitting/data/affiliationstrings_mapping.csv",
    # column names in ids file:
    ids_id_col = "id1",         # <-- change if different
    ids_text_col = "affil1",     # <-- change if different
    # column names in mapping file:
    #map_id_col = "id1",         # <-- change if different
    #map_entity_col = "id2"      # <-- change if different

    pair_id1_col = "id1",
    pair_id2_col = "id2"

  ),

  # clean
  affiliation_clean <- list(
    enabled = FALSE,
    records_path = "D:/datasplitting/data/clean_affiliations_2024_05_15.csv",
    mapping_path = "D:/datasplitting/data/affiliationstrings_mapping.csv",
    id_col = "ID",
    pair_id1_col = "id1",
    pair_id2_col = "id2",
    text_cols = c("Name1","Name2","Name3","Street1","Street2","City","State","Zipcode","Country")
  ),


  # 3) 10K synthetic (from your screenshot folder)
  # synth10k = list(
  #   enabled = FALSE,
  #   full_path = "D:/datasplitting/data/10Kfull.csv",
  #   dup_pairs_path = "D:/datasplitting/data/10Kduplicates.csv",
  #   # If full file already contains entity_id, set these:
  #   full_id_col = "Id",               # <-- change if needed
  #   # If full contains a raw string field:
  #   full_text_col = "Aggregate Value",# <-- change if needed
  #   # If duplicates file is gold pairs, set these:
  #   dup_id1_col = "id1",              # <-- change if needed
  #   dup_id2_col = "id2"               # <-- change if needed
  #   ),
  #
  #  # 4) NC voter (YOU MUST PROVIDE THESE FILES + ground truth)
  # ncvoter = list(
  #   enabled = FALSE,                  # set TRUE when you provide paths + mapping/pairs
  #   records_path = "D:/datasplitting/data/ncvoter_records.csv",
  #   mapping_or_pairs_path = "D:/datasplitting/data/ncvoter_gold.csv",
  #   # records
  #   rec_id_col = "id",                # <-- change
  #   rec_text_cols = c("first_name","last_name","address"),  # <-- change
  #   # mapping/pairs schema: choose ONE
  #   gold_type = "mapping",            # "mapping" or "pairs"
  #   map_id_col = "id",                # <-- change
  #   map_entity_col = "entity_id",     # <-- change
  #   pair_id1_col = "id1",             # <-- change
  #   pair_id2_col = "id2"              # <-- change
  #   )
)
#
# ============================================================
# Shared utilities (do not edit)
# ============================================================

norm_str <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}


# VI measures information-theoretic distance between partitions.
# aricode VI not exported
#VI_safe <- function(x, y) aricode:::VI(x, y)
# Compute VI manually
VI_manual <- function(x, y, eps = 1e-15) {
  x <- as.integer(factor(x))
  y <- as.integer(factor(y))

  n <- length(x)
  tab <- table(x, y)
  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)

  # entropies
  Hx <- -sum(px * log(px + eps))
  Hy <- -sum(py * log(py + eps))

  # mutual information
  # only sum where pxy > 0
  nz <- pxy > 0
  Ixy <- sum(pxy[nz] * log(pxy[nz] / (px[row(pxy)[nz]] * py[col(pxy)[nz]] + eps) + eps))

  Hx + Hy - 2 * Ixy
}

# FMI (Fowlkes–Mallows Index) as a pair-based substitute
FMI <- function(true_labels, pred_labels) {
  tab <- table(true_labels, pred_labels)
  # TP = sum over clusters of choose(n_ij, 2)
  tp <- sum(choose(tab, 2))
  # FP = sum_j choose(n_.j,2) - TP
  fp <- sum(choose(colSums(tab), 2)) - tp
  # FN = sum_i choose(n_i.,2) - TP
  fn <- sum(choose(rowSums(tab), 2)) - tp

  if (tp == 0) return(0)
  sqrt((tp / (tp + fp)) * (tp / (tp + fn)))
}


# Build entity_id from gold pairs using connected components
entity_from_pairs <- function(ids, pairs_dt, id1_col, id2_col) {
  ids <- sort(unique(as.integer(ids)))
  g <- make_empty_graph(n = length(ids), directed = FALSE)
  V(g)$name <- as.character(ids)

  e <- pairs_dt[get(id1_col) %in% ids & get(id2_col) %in% ids,
                .(id1 = as.character(get(id1_col)), id2 = as.character(get(id2_col)))]
  if (nrow(e) > 0) {
    g <- add_edges(g, as.vector(t(as.matrix(e))))
  }
  comp <- components(g)$membership
  data.table(id = ids, entity_id = as.integer(comp[as.character(ids)]))
}

# generate all pairs within a set of ids (no blocking)
make_all_pairs <- function(ids) {
  ids <- sort(unique(as.integer(ids)))
  n <- length(ids)
  if (n < 2) return(data.table(id1 = integer(), id2 = integer()))
  id1 <- rep(ids[1:(n - 1)], times = (n - 1):1)
  id2 <- unlist(lapply(2:n, function(i) ids[i:n]), use.names = FALSE)
  data.table(id1 = id1, id2 = id2)
}

# Splits
split_record_random <- function(d, test_frac = 0.3, seed = 1) {
  set.seed(seed)
  ids <- sample(d$id, length(d$id), replace = FALSE)
  n_test <- floor(test_frac * length(ids))
  list(test_ids = ids[1:n_test], train_ids = ids[(n_test + 1):length(ids)])
}

split_entity_disjoint <- function(d, test_frac = 0.3, seed = 1) {
  set.seed(seed)
  ents <- unique(d$entity_id)
  ents <- sample(ents, length(ents), replace = FALSE)
  n_test <- floor(test_frac * length(ents))
  test_ents <- ents[1:n_test]
  list(
    test_ids  = d[entity_id %in% test_ents, id],
    train_ids = d[!(entity_id %in% test_ents), id]
  )
}

# Features: generic text-only (multiple columns)
# Uses JW + normalized LV per column; concatenates into X
pair_features_text <- function(d, pairs, text_cols) {
  setkey(d, id)
  a <- d[pairs$id1]
  b <- d[pairs$id2]

  X_list <- list()
  for (col in text_cols) {
    xa <- a[[col]]; xb <- b[[col]]
    # JW sim
    jw <- 1 - stringdist(xa, xb, method = "jw", p = 0.1)
    # normalized LV sim
    lv_d <- stringdist(xa, xb, method = "lv")
    lv_m <- pmax(nchar(xa), nchar(xb))
    lv <- 1 - ifelse(lv_m == 0, 0, lv_d / lv_m)

    X_list[[paste0(col, "_jw")]] <- jw
    X_list[[paste0(col, "_lv")]] <- lv
  }
  X <- as.data.table(X_list)
  y <- as.integer(a$entity_id == b$entity_id)
  list(X = X, y = y)
}

# Train-only mean imputation
fit_imputer <- function(X_train) {
  mu <- vapply(X_train, function(col) mean(col, na.rm = TRUE), numeric(1))
  mu[is.nan(mu)] <- 0
  mu
}

apply_imputer <- function(X, mu) {
  X2 <- copy(X)
  for (nm in names(mu)) X2[[nm]][is.na(X2[[nm]])] <- mu[[nm]]
  X2
}

# Subsample negatives in TRAIN (critical without blocking)
subsample_train <- function(X, y, neg_ratio = 10, seed = 1) {
  set.seed(seed)
  pos <- which(y == 1)
  neg <- which(y == 0)
  if (length(pos) == 0) {
    keep <- sample(neg, min(length(neg), 5000))
    return(list(X = X[keep], y = y[keep]))
  }
  n_neg_keep <- min(length(neg), neg_ratio * length(pos))
  keep <- c(pos, sample(neg, n_neg_keep))
  keep <- sample(keep, length(keep))
  list(X = X[keep], y = y[keep])
}

# Pairwise F1 + threshold tuning (NA-safe)
pair_f1 <- function(yhat, ytrue) {
  ok <- !(is.na(yhat) | is.na(ytrue))
  yhat <- yhat[ok]; ytrue <- ytrue[ok]
  if (length(ytrue) == 0) return(NA_real_)
  tp <- sum(yhat == 1 & ytrue == 1)
  fp <- sum(yhat == 1 & ytrue == 0)
  fn <- sum(yhat == 0 & ytrue == 1)
  prec <- ifelse(tp + fp == 0, 0, tp / (tp + fp))
  rec  <- ifelse(tp + fn == 0, 0, tp / (tp + fn))
  ifelse(prec + rec == 0, 0, 2 * prec * rec / (prec + rec))
}

best_threshold_f1 <- function(proba, y_true, grid = seq(0.10, 0.90, by = 0.02)) {
  ok <- !(is.na(proba) | is.na(y_true))
  proba <- proba[ok]; y_true <- y_true[ok]
  if (length(y_true) == 0) return(list(threshold = 0.5, f1 = NA_real_))
  best_t <- 0.5; best_f <- -Inf
  for (t in grid) {
    f <- pair_f1(as.integer(proba >= t), y_true)
    if (is.na(f)) next
    if (f > best_f) { best_f <- f; best_t <- t }
  }
  if (!is.finite(best_f)) return(list(threshold = 0.5, f1 = NA_real_))
  list(threshold = best_t, f1 = best_f)
}

# Cluster from thresholded edges (connected components)
cluster_from_edges <- function(record_ids, edges_dt) {
  record_ids <- as.integer(record_ids)
  g <- make_empty_graph(n = length(record_ids), directed = FALSE)
  V(g)$name <- as.character(record_ids)
  if (nrow(edges_dt) > 0) {
    g <- add_edges(g, as.vector(t(as.matrix(edges_dt[, .(as.character(id1), as.character(id2))]))))
  }
  comp <- components(g)$membership
  as.integer(comp[as.character(record_ids)])
}

# Train classifier with glmnet (ridge) for stability
fit_predict_glmnet <- function(X_train, y_train, X_test, seed = 1) {
  set.seed(seed)
  xtr <- as.matrix(X_train)
  xte <- as.matrix(X_test)
  ytr <- as.numeric(y_train)

  # cv.glmnet provides lambda.min and lambda.1se
  cvfit <- glmnet::cv.glmnet(
    x = xtr, y = ytr,
    family = "binomial",
    alpha = 0,              # ridge
    nfolds = 5,
    standardize = TRUE
  )

  p_tr <- as.numeric(predict(cvfit, newx = xtr, type = "response", s = "lambda.1se"))
  p_te <- as.numeric(predict(cvfit, newx = xte, type = "response", s = "lambda.1se"))

  p_tr[!is.finite(p_tr)] <- NA_real_
  p_te[!is.finite(p_te)] <- NA_real_

  list(p_train = p_tr, p_test = p_te, lambda = cvfit$lambda.1se)
}


# ============================================================
# Dataset adapters (edit only if your columns differ)
# ============================================================

# CORA adapter
load_cora <- function() {
  if (!requireNamespace("cora", quietly = TRUE)) install.packages("cora")
  library(cora)
  data("cora", package = "cora")
  data("cora_gold", package = "cora")

  d <- as.data.table(cora)
  g <- as.data.table(cora_gold)

  setnames(d, "id", "id")
  # build entity_id from gold pairs
  ent <- entity_from_pairs(d$id, g, "id1", "id2")
  d <- ent[d, on = "id"]

  # choose text cols (normalize)
  d[, title_n := norm_str(title)]
  d[, authors_n := norm_str(authors)]
  d[, journal_n := norm_str(journal)]
  list(dt = d[, .(id, entity_id, title_n, authors_n, journal_n)],
       text_cols = c("title_n", "authors_n", "journal_n"))
}

# Affiliation adapter (mapping provides entity_id)
load_affiliation <- function(cfg) {
  library(data.table)

  # ----------------------------
  # 1) Load records: id1, affil1
  # ----------------------------
  d <- fread(cfg$records_path)

  stopifnot(all(c("id1", "affil1") %in% names(d)))

  setnames(d, c("id1", "affil1"), c("id", "affil_raw"))
  d[, id := as.integer(id)]

  # ----------------------------
  # 2) Load gold pairs: id1, id2
  # ----------------------------
  pairs <- fread(cfg$mapping_path)

  stopifnot(all(c("id1", "id2") %in% names(pairs)))

  pairs[, `:=`(
    id1 = as.integer(id1),
    id2 = as.integer(id2)
  )]

  # ----------------------------
  # 3) entity_id via CC (existing function)
  # ----------------------------
  ent <- entity_from_pairs(
    ids = d$id,
    pairs_dt = pairs,
    id1_col = "id1",
    id2_col = "id2"
  )

  d <- merge(d, ent, by = "id", all.x = TRUE)

  stopifnot(!anyNA(d$entity_id))

  # ----------------------------
  # 4) Normalize text
  # ----------------------------
  d[, affil_n := norm_str(affil_raw)]

  # ----------------------------
  # 5) Return standardized object
  # ----------------------------
  list(
    dt = d[, .(id, entity_id, affil_n)],
    text_cols = "affil_n"
  )
}



# Clean, multifield Affiliation adapter (mapping provides entity_id)

load_affiliation_multifield <- function(cfg) {
  library(data.table)

  # 1) records
  d <- fread(cfg$records_path)
  stopifnot(cfg$id_col %in% names(d))
  setnames(d, cfg$id_col, "id")
  d[, id := as.integer(id)]

  # 2) gold pairs (id1, id2)
  pairs <- fread(cfg$mapping_path)

  # standardize pair columns robustly
  if (!all(c("id1", "id2") %in% names(pairs))) {
    stopifnot(!is.null(cfg$pair_id1_col), !is.null(cfg$pair_id2_col))
    stopifnot(all(c(cfg$pair_id1_col, cfg$pair_id2_col) %in% names(pairs)))
    setnames(pairs, c(cfg$pair_id1_col, cfg$pair_id2_col), c("id1", "id2"))
  }
  pairs[, `:=`(id1 = as.integer(id1), id2 = as.integer(id2))]

  # 3) entity_id via connected components
  ent <- entity_from_pairs(d$id, pairs, id1_col = "id1", id2_col = "id2")

  # left-join to keep all records (explicit, robust)
  d <- merge(d, ent, by = "id", all.x = TRUE)
  stopifnot(!anyNA(d$entity_id))

  # 4) normalize multi-fields (CORA-style)
  stopifnot(all(cfg$text_cols %in% names(d)))
  for (col in cfg$text_cols) {
    newc <- paste0(col, "_n")
    setnames(d, col, newc)
    d[, (newc) := norm_str(get(newc))]
  }
  text_cols_n <- paste0(cfg$text_cols, "_n")

  list(
    dt = d[, c("id", "entity_id", text_cols_n), with = FALSE],
    text_cols = text_cols_n
  )
}




# 10K synthetic adapter:
# - If dup_pairs_path is gold pairs: build entity_id by CC
# - full_path may be pipe-delimited in one column; fix automatically
load_synth10k <- function(cfg) {
  d <- fread(cfg$full_path, showProgress = FALSE)

  # If the CSV was parsed into one column containing pipes, split it
  if (ncol(d) == 1 && grepl("\\|", names(d)[1], fixed = FALSE)) {
    # header is in the column name itself
    header <- names(d)[1]
    cols <- strsplit(header, "\\|", fixed = FALSE)[[1]]
    setnames(d, header, "V1")
    d <- d[, tstrsplit(V1, "\\|", fixed = FALSE)]
    setnames(d, names(d), cols)
  }

  # standardize id column
  setnames(d, cfg$full_id_col, "id")
  d[, id := as.integer(id)]

  # gold pairs -> entity_id
  g <- fread(cfg$dup_pairs_path, showProgress = FALSE)

  # If dup file is also pipe-delimited
  if (ncol(g) == 1 && grepl("\\|", names(g)[1], fixed = FALSE)) {
    header <- names(g)[1]
    cols <- strsplit(header, "\\|", fixed = FALSE)[[1]]
    setnames(g, header, "V1")
    g <- g[, tstrsplit(V1, "\\|", fixed = FALSE)]
    setnames(g, names(g), cols)
  }

  ent <- entity_from_pairs(d$id, g, cfg$dup_id1_col, cfg$dup_id2_col)
  d <- ent[d, on = "id"]

  # text field
  setnames(d, cfg$full_text_col, "s_raw")
  d[, s_n := norm_str(s_raw)]

  list(dt = d[, .(id, entity_id, s_n)],
       text_cols = c("s_n"))
}

# NC voter adapter (requires mapping or pairs)
load_ncvoter <- function(cfg) {
  d <- fread(cfg$records_path)
  setnames(d, cfg$rec_id_col, "id")
  d[, id := as.integer(id)]

  # normalize text cols
  for (c in cfg$rec_text_cols) {
    newc <- paste0(c, "_n")
    setnames(d, c, newc)
    d[, (newc) := norm_str(get(newc))]
  }
  text_cols <- paste0(cfg$rec_text_cols, "_n")

  g <- fread(cfg$mapping_or_pairs_path)

  if (cfg$gold_type == "mapping") {
    setnames(g, c(cfg$map_id_col, cfg$map_entity_col), c("id","entity_id"))
    d <- g[d, on = "id"]
  } else if (cfg$gold_type == "pairs") {
    ent <- entity_from_pairs(d$id, g, cfg$pair_id1_col, cfg$pair_id2_col)
    d <- ent[d, on = "id"]
  } else stop("ncvoter gold_type must be 'mapping' or 'pairs'")

  stopifnot(all(c("id","entity_id") %in% names(d)))
  list(dt = d[, c("id","entity_id", text_cols), with = FALSE],
       text_cols = text_cols)
}

# ============================================================
# Core experiment (shared)
# ============================================================

run_one <- function(d0, text_cols, scheme = c("record","entity"), test_frac = 0.3,
                    seed = 2026, neg_ratio = 10) {
  scheme <- match.arg(scheme)

  sp <- if (scheme == "record") split_record_random(d0, test_frac, seed) else split_entity_disjoint(d0, test_frac, seed)
  tr <- d0[id %in% sp$train_ids]
  te <- d0[id %in% sp$test_ids]

  pairs_tr <- make_all_pairs(tr$id)
  pairs_te <- make_all_pairs(te$id)
  if (nrow(pairs_tr) == 0 || nrow(pairs_te) == 0) return(NULL)

  ftr <- pair_features_text(tr, pairs_tr, text_cols)
  fte <- pair_features_text(te, pairs_te, text_cols)

  Xtr <- ftr$X; ytr <- ftr$y
  Xte <- fte$X; yte <- fte$y

  mu <- fit_imputer(Xtr)
  Xtr <- apply_imputer(Xtr, mu)
  Xte <- apply_imputer(Xte, mu)

  sub <- subsample_train(Xtr, ytr, neg_ratio = neg_ratio, seed = seed)
  Xtr_s <- sub$X; ytr_s <- sub$y

  #pp <- fit_predict_glmnet(Xtr_s, ytr_s, Xte)
  pp <- fit_predict_glmnet(Xtr_s, ytr_s, Xte, seed = seed)

  th <- best_threshold_f1(pp$p_train, ytr_s)$threshold

  yhat_te <- as.integer(pp$p_test >= th)
  f1_te <- pair_f1(yhat_te, yte)

  edges_te <- pairs_te[yhat_te == 1, .(id1, id2)]
  rec_ids <- te$id
  pred_cluster <- cluster_from_edges(rec_ids, edges_te)
  true_cluster <- as.integer(factor(te[match(rec_ids, id), entity_id]))

  list(
    ARI = mclust::adjustedRandIndex(true_cluster, pred_cluster),
    NMI = aricode::NMI(true_cluster, pred_cluster),
    VI  = VI_manual(true_cluster, pred_cluster),
    FMI  = FMI(true_cluster, pred_cluster),
    PairF1 = f1_te,
    threshold = th,
    n_train_records = nrow(tr),
    n_test_records  = nrow(te),
    n_train_pairs   = nrow(pairs_tr),
    n_test_pairs    = nrow(pairs_te),
    pos_rate_train  = mean(ytr),
    pos_rate_test   = mean(yte)
  )
}

run_experiment_dataset <- function(dataset_name, d0, text_cols, K = 10, test_frac = 0.30,
                                   base_seed = 2026, neg_ratio = 10) {
  out <- list()
  for (scheme in c("record","entity")) {
    for (k in seq_len(K)) {
      seed <- base_seed + k
      r <- run_one(d0, text_cols, scheme, test_frac, seed, neg_ratio)
      if (is.null(r)) next
      out[[length(out) + 1]] <- data.table(
        dataset = dataset_name,
        scheme = scheme,
        run = k,
        seed = seed,
        ARI = r$ARI, NMI = r$NMI, VI = r$VI, FMI = r$FMI,
        PairF1 = r$PairF1,
        threshold = r$threshold,
        n_train_records = r$n_train_records, n_test_records = r$n_test_records,
        n_train_pairs = r$n_train_pairs, n_test_pairs = r$n_test_pairs,
        pos_rate_train = r$pos_rate_train, pos_rate_test = r$pos_rate_test
      )
    }
  }
  rbindlist(out)
}

summarize_all <- function(res_all) {
  res_all[, .(
    n_runs = .N,
    ARI_mean = mean(ARI, na.rm = TRUE), ARI_sd = sd(ARI, na.rm = TRUE),
    NMI_mean = mean(NMI, na.rm = TRUE), NMI_sd = sd(NMI, na.rm = TRUE),
    VI_mean  = mean(VI,  na.rm = TRUE), VI_sd  = sd(VI,  na.rm = TRUE),
    FMI_mean  = mean(FMI,  na.rm = TRUE), FMI_sd  = sd(FMI,  na.rm = TRUE),
    PairF1_mean = mean(PairF1, na.rm = TRUE), PairF1_sd = sd(PairF1, na.rm = TRUE),
    thr_mean = mean(threshold, na.rm = TRUE), thr_sd = sd(threshold, na.rm = TRUE),
    median_test_pairs = median(n_test_pairs)
  ), by = .(dataset, scheme)][order(dataset, scheme)]
}

# ============================================================
# RUN ALL ENABLED DATASETS
# ============================================================

K <- 10
test_frac <- 0.30
base_seed <- 2026
neg_ratio <- 10  # reduce for speed on huge datasets

res_list <- list()

system.time(
  if (isTRUE(CFG$cora$enabled)) {
    cora_obj <- load_cora()
    res_list[["cora"]] <- run_experiment_dataset("cora", cora_obj$dt, cora_obj$text_cols,
                                                 K, test_frac, base_seed, neg_ratio)
  }
)

# K=10
# user  system elapsed
# 488.47   21.15  410.11

# K=20
# user  system elapsed
# 1714.88   38.53 1555.11


system.time(
  if (isTRUE(CFG$affiliation$enabled)) {
    aff_obj <- load_affiliation(CFG$affiliation)
    res_list[["affiliation"]] <- run_experiment_dataset("affiliation", aff_obj$dt, aff_obj$text_cols,
                                                        K, test_frac, base_seed, neg_ratio)
  }
)

# K=10
# user  system elapsed
# 25.54    1.00   15.14

# K=20
# user  system elapsed
# 47.22    2.66   27.19

system.time(
  if (isTRUE(CFG$affiliation_clean$enabled)) {
    aff_obj <- load_affiliation_multifield(CFG$affiliation_clean)
    res_list[["affiliation_clean"]] <- run_experiment_dataset("affiliation_clean", aff_obj$dt, aff_obj$text_cols,
                                                              K, test_frac, base_seed, neg_ratio)
  }
)


#### leave it to later, add blocking
# if (isTRUE(CFG$synth10k$enabled)) {
#   syn_obj <- load_synth10k(CFG$synth10k)
#   res_list[["synth10k"]] <- run_experiment_dataset("synth10k", syn_obj$dt, syn_obj$text_cols,
#                                                    K, test_frac, base_seed, neg_ratio)
# }
#
# if (isTRUE(CFG$ncvoter$enabled)) {
#   nc_obj <- load_ncvoter(CFG$ncvoter)
#   res_list[["ncvoter"]] <- run_experiment_dataset("ncvoter", nc_obj$dt, nc_obj$text_cols,
#                                                   K, test_frac, base_seed, neg_ratio)
# }


# Analysis and visualization
res_all <- rbindlist(res_list, fill = TRUE)
final_all <- summarize_all(res_all)

print(res_all)
print(final_all)


plot_er_results <- function(final_all,
                            metric = c("ARI", "NMI", "VI", "FMI", "PairF1"),
                            show_error = TRUE) {

  metric <- match.arg(metric)

  library(ggplot2)
  library(data.table)

  # ---- map metric names to columns ----
  col_map <- list(
    ARI     = c(mean = "ARI_mean",     sd = "ARI_sd",     ylab = "Adjusted Rand Index (ARI)"),
    NMI     = c(mean = "NMI_mean",     sd = "NMI_sd",     ylab = "Normalized Mutual Information (NMI)"),
    VI     = c(mean = "VI_mean",     sd = "VI_sd",     ylab = "Variation of Information (VI)"),
    FMI     = c(mean = "FMI_mean",     sd = "FMI_sd",     ylab = "Fowlkes–Mallows Index (FMI)"),
    PairF1  = c(mean = "PairF1_mean",  sd = "PairF1_sd",  ylab = "Pairwise F1")
  )

  m <- col_map[[metric]]

  dt <- as.data.table(final_all)

  p <- ggplot(
    dt,
    aes(x = scheme, y = get(m["mean"]), fill = scheme)
  ) +
    geom_col(width = 0.65, alpha = 0.85) +
    facet_wrap(~ dataset, scales = "free_y") +
    labs(
      title = paste(metric, "by split strategy and dataset"),
      subtitle = "Mean across runs (± SD)",
      x = NULL,
      y = m["ylab"]
    ) +
    theme_minimal(base_size = 13) +
    theme(
      legend.position = "none",
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    )

  if (show_error) {
    p <- p +
      geom_errorbar(
        aes(
          ymin = get(m["mean"]) - get(m["sd"]),
          ymax = get(m["mean"]) + get(m["sd"])
        ),
        width = 0.2,
        linewidth = 0.6
      )
  }

  return(p)
}

plot_er_results(final_all, metric = "ARI")
plot_er_results(final_all, metric = "NMI")
plot_er_results(final_all, metric = "VI")
plot_er_results(final_all, metric = "FMI")
plot_er_results(final_all, metric = "PairF1")


plot_thresholds <- function(final_all) {
  library(ggplot2)
  ggplot(final_all,
         aes(x = scheme, y = thr_mean, fill = scheme)) +
    geom_col(width = 0.65, alpha = 0.85) +
    geom_errorbar(
      aes(ymin = thr_mean - thr_sd, ymax = thr_mean + thr_sd),
      width = 0.2
    ) +
    facet_wrap(~ dataset, scales = "free_y") +
    labs(
      title = "Decision threshold by dataset and split strategy",
      x = NULL,
      y = "Match probability threshold"
    ) +
    theme_minimal(base_size = 13) +
    theme(legend.position = "none")
}

plot_thresholds(final_all)


