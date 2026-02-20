# ============================================================
# Unified ER splitting experiment for 4 datasets:
#   1) CORA (package cora) [gold pairs -> entity clusters]
#   2) Affiliation (csv ids + mapping entity_id)
#   3) Affiliation_clean (multi-column cleaned affiliations)
#   4) 10K synthetic (full + duplicates pairs OR mapping)
#   5) NC voter (records + mapping OR gold pairs)
#
# Outputs:
#   res_all   : per-run results for each dataset x scheme
#   final_all : mean/sd summary by dataset x scheme
# ============================================================

suppressPackageStartupMessages({
  pkgs <- c("data.table", "stringdist", "igraph", "mclust", "aricode", "glmnet", "phonics")
  for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(data.table)
  library(stringdist)
  library(igraph)
  library(mclust)
  library(aricode)
  library(glmnet)
  library(phonics)
})

# -----------------------------
# CONFIG: EDIT ONLY THIS BLOCK
# -----------------------------
CFG <- list(
  # 1) CORA (no paths needed)
  cora = list(enabled = TRUE),

  # 2) Affiliation (original single-column text)
  affiliation = list(
    enabled = TRUE,
    records_path = "D:/datasplitting/data/affiliationstrings_ids.csv",
    mapping_path = "D:/datasplitting/data/affiliationstrings_mapping.csv",
    ids_id_col = "id1",
    ids_text_col = "affil1",
    pair_id1_col = "id1",
    pair_id2_col = "id2"
  ),

  # 3) Affiliation clean (multi-column) - FIXED SYNTAX
  affiliation_clean = list(
    enabled = TRUE,
    records_path = "D:/datasplitting/data/clean_affiliations_2024_05_15.csv",
    mapping_path = "D:/datasplitting/data/affiliationstrings_mapping.csv",
    id_col = "ID",
    pair_id1_col = "id1",
    pair_id2_col = "id2",
    text_cols = c("Name1","Name2","Name3","Street1","Street2","City","State","Zipcode","Country")
  ),

  # 4) 10K synthetic dataset (D10K)
  synth10k = list(
    enabled = TRUE,  # Set to TRUE when ready to run
    records_path = "D:/datasplitting/data/10Kfull.csv",
    pairs_path = "D:/datasplitting/data/10Kduplicates.csv",
    # Column names in records file
    id_col = "Id",
    text_cols = c("Aggregate Value", "Embedded Ag.Value", "Clean Ag.Value", "Embedded Clean Ag.Value"),
    # Column names in ground truth pairs file
    pair_id1_col = "Entity1",
    pair_id2_col = "Entity2",
    # Blocking configuration
    use_blocking = TRUE,  # Set to FALSE to disable blocking
    block_on = "Aggregate Value",  # Which column to use for blocking
    block_method = "prefix",  # Options: "prefix", "soundex", "ngram"
    block_size = 3  # For prefix: num characters; for ngram: n-gram size
  )#,
  #
  # 5) NC voter
  # ncvoter = list(
  #   enabled = FALSE,
  #   records_path = "D:/datasplitting/data/ncvoter_records.csv",
  #   mapping_or_pairs_path = "D:/datasplitting/data/ncvoter_gold.csv",
  #   rec_id_col = "id",
  #   rec_text_cols = c("first_name","last_name","address"),
  #   gold_type = "mapping",
  #   map_id_col = "id",
  #   map_entity_col = "entity_id",
  #   pair_id1_col = "id1",
  #   pair_id2_col = "id2"
  # )
)

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
  nz <- pxy > 0
  Ixy <- sum(pxy[nz] * log(pxy[nz] / (px[row(pxy)[nz]] * py[col(pxy)[nz]] + eps) + eps))

  Hx + Hy - 2 * Ixy
}

# FMI (Fowlkes–Mallows Index) as a pair-based substitute
FMI <- function(true_labels, pred_labels) {
  tab <- table(true_labels, pred_labels)
  tp <- sum(choose(tab, 2))
  fp <- sum(choose(colSums(tab), 2)) - tp
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

# ============================================================
# BLOCKING FUNCTIONS
# ============================================================

# Generate blocking key based on method
generate_block_key <- function(text, method = "prefix", size = 3) {
  if (method == "prefix") {
    # Use first N characters
    return(substr(text, 1, size))
  } else if (method == "soundex") {
    # Phonetic encoding (requires phonics package)
    if (!requireNamespace("phonics", quietly = TRUE)) {
      warning("phonics package not installed. Falling back to prefix blocking.")
      return(substr(text, 1, size))
    }
    return(phonics::soundex(text))
  } else if (method == "ngram") {
    # Use first n-gram
    if (nchar(text) < size) return(text)
    return(substr(text, 1, size))
  } else {
    stop(paste("Unknown blocking method:", method))
  }
}

# Make pairs within blocks
make_blocked_pairs <- function(ids, d, block_col, method = "prefix", size = 3) {
  # Subset data to only include ids we're working with
  setkey(d, id)
  d_subset <- d[id %in% ids, .(id, text = get(block_col))]

  # Generate blocking keys
  d_subset[, block_key := generate_block_key(text, method, size)]

  # Generate pairs within each block
  pairs_list <- list()
  blocks <- unique(d_subset$block_key)

  for (key in blocks) {
    block_ids <- d_subset[block_key == key, id]
    if (length(block_ids) >= 2) {
      pairs_list[[length(pairs_list) + 1]] <- make_all_pairs(block_ids)
    }
  }

  if (length(pairs_list) == 0) {
    return(data.table(id1 = integer(), id2 = integer()))
  }

  # Combine all pairs and remove duplicates
  all_pairs <- rbindlist(pairs_list)
  all_pairs <- unique(all_pairs)

  return(all_pairs)
}

# Wrapper function that decides whether to use blocking
make_pairs_with_blocking <- function(ids, d = NULL, use_blocking = FALSE,
                                     block_col = NULL, method = "prefix", size = 3) {
  if (!use_blocking || is.null(d) || is.null(block_col)) {
    # No blocking - generate all pairs
    return(make_all_pairs(ids))
  } else {
    # Use blocking
    return(make_blocked_pairs(ids, d, block_col, method, size))
  }
}

# Calculate blocking statistics
blocking_stats <- function(d, block_col, method = "prefix", size = 3) {
  d_copy <- copy(d)
  d_copy[, block_key := generate_block_key(get(block_col), method, size)]

  stats <- list(
    n_records = nrow(d),
    n_blocks = length(unique(d_copy$block_key)),
    avg_block_size = mean(table(d_copy$block_key)),
    max_block_size = max(table(d_copy$block_key)),
    min_block_size = min(table(d_copy$block_key))
  )

  # Calculate pairs per block and sum
  block_sizes <- table(d_copy$block_key)
  pairs_per_block <- choose(block_sizes, 2)
  pairs_per_block[is.nan(pairs_per_block)] <- 0

  stats$total_pairs_blocked <- sum(pairs_per_block)
  stats$total_pairs_unblocked <- choose(nrow(d), 2)
  stats$reduction_ratio <- stats$total_pairs_blocked / stats$total_pairs_unblocked

  return(stats)
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
  if (length(yhat) == 0) return(NA_real_)
  tp <- sum(yhat == 1 & ytrue == 1)
  fp <- sum(yhat == 1 & ytrue == 0)
  fn <- sum(yhat == 0 & ytrue == 1)
  if ((tp + fp) == 0 || (tp + fn) == 0) return(0)
  prec <- tp / (tp + fp)
  rec <- tp / (tp + fn)
  if (prec + rec == 0) return(0)
  2 * prec * rec / (prec + rec)
}

best_threshold_f1 <- function(p_train, y_train) {
  ok <- !(is.na(p_train) | is.na(y_train))
  p_train <- p_train[ok]; y_train <- y_train[ok]
  if (length(p_train) < 2) return(list(threshold = 0.5, f1 = 0))

  cands <- sort(unique(p_train))
  if (length(cands) > 100) {
    cands <- quantile(p_train, probs = seq(0, 1, length.out = 100))
  }
  f1_vals <- vapply(cands, function(th) {
    yhat <- as.integer(p_train >= th)
    pair_f1(yhat, y_train)
  }, numeric(1))

  best_idx <- which.max(f1_vals)
  list(threshold = cands[best_idx], f1 = f1_vals[best_idx])
}

# Classifier (glmnet)
fit_predict_glmnet <- function(X_train, y_train, X_test, seed = 1, alpha = 1) {
  set.seed(seed)
  ok_tr <- complete.cases(X_train)
  ok_te <- complete.cases(X_test)

  Xtr <- as.matrix(X_train[ok_tr])
  ytr <- y_train[ok_tr]
  Xte <- as.matrix(X_test[ok_te])

  if (length(unique(ytr)) < 2 || nrow(Xtr) < 10) {
    return(list(
      p_train = rep(mean(ytr), nrow(X_train)),
      p_test = rep(mean(ytr), nrow(X_test))
    ))
  }

  cv_fit <- tryCatch(
    cv.glmnet(Xtr, ytr, family = "binomial", alpha = alpha, nfolds = 5),
    error = function(e) NULL
  )
  if (is.null(cv_fit)) {
    return(list(
      p_train = rep(mean(ytr), nrow(X_train)),
      p_test = rep(mean(ytr), nrow(X_test))
    ))
  }

  p_train <- numeric(nrow(X_train))
  p_train[ok_tr] <- as.vector(predict(cv_fit, Xtr, s = "lambda.min", type = "response"))
  p_train[!ok_tr] <- mean(ytr)

  p_test <- numeric(nrow(X_test))
  p_test[ok_te] <- as.vector(predict(cv_fit, Xte, s = "lambda.min", type = "response"))
  p_test[!ok_te] <- mean(ytr)

  list(p_train = p_train, p_test = p_test)
}

cluster_from_edges <- function(rec_ids, edges) {
  g <- make_empty_graph(n = length(rec_ids), directed = FALSE)
  V(g)$name <- as.character(rec_ids)
  if (nrow(edges) > 0) {
    e_chr <- edges[, .(id1 = as.character(id1), id2 = as.character(id2))]
    g <- add_edges(g, as.vector(t(as.matrix(e_chr))))
  }
  components(g)$membership
}

# ============================================================
# LOADERS
# ============================================================

load_cora <- function() {
  if (!requireNamespace("cora", quietly = TRUE)) install.packages("cora")
  library(cora)
  # Load cora data into a new environment to avoid polluting workspace
  #cora_env <- new.env()
  #data("cora", package = "cora", envir = cora_env)

  # Extract the objects from the environment
  records <- cora
  gold_pairs <- cora_gold

  if (is.null(records) || is.null(gold_pairs)) {
    stop("CORA data not loaded properly. Check package installation.")
  }

  id <- seq_len(nrow(records))
  setDT(gold_pairs)

  ent_dt <- entity_from_pairs(id, gold_pairs, "id1", "id2")
  dt <- data.table(
    id = id,
    entity_id = ent_dt$entity_id,
    title = norm_str(records$title),
    authors = norm_str(records$authors),
    journal = norm_str(records$journal)
  )
  text_cols <- c("title", "authors", "journal")
  list(dt = dt, text_cols = text_cols)
}

load_affiliation <- function(cfg) {
  rec <- fread(cfg$records_path, encoding = "UTF-8")
  setnames(rec, old = c(cfg$ids_id_col, cfg$ids_text_col), new = c("id", "text"), skip_absent = TRUE)

  pairs <- fread(cfg$mapping_path, encoding = "UTF-8")
  setnames(pairs, old = c(cfg$pair_id1_col, cfg$pair_id2_col), new = c("id1", "id2"), skip_absent = TRUE)

  ent_dt <- entity_from_pairs(rec$id, pairs, "id1", "id2")
  rec[, entity_id := ent_dt$entity_id]
  rec[, text := norm_str(text)]

  list(dt = rec[, .(id, entity_id, text)], text_cols = "text")
}

# NEW LOADER FOR MULTI-FIELD AFFILIATION
load_affiliation_multifield <- function(cfg) {
  rec <- fread(cfg$records_path, encoding = "UTF-8")

  # Rename ID column
  if (cfg$id_col %in% names(rec)) {
    setnames(rec, old = cfg$id_col, new = "id")
  } else {
    stop(paste("ID column", cfg$id_col, "not found in", cfg$records_path))
  }

  # Check all text columns exist
  missing_cols <- setdiff(cfg$text_cols, names(rec))
  if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
    cfg$text_cols <- intersect(cfg$text_cols, names(rec))
  }

  # Load pairs/mapping
  pairs <- fread(cfg$mapping_path, encoding = "UTF-8")
  setnames(pairs,
           old = c(cfg$pair_id1_col, cfg$pair_id2_col),
           new = c("id1", "id2"),
           skip_absent = TRUE)

  # Build entity_id from pairs
  ent_dt <- entity_from_pairs(rec$id, pairs, "id1", "id2")
  rec[, entity_id := ent_dt$entity_id]

  # Normalize all text columns
  for (col in cfg$text_cols) {
    set(rec, j = col, value = norm_str(rec[[col]]))
  }

  # Keep only necessary columns
  keep_cols <- c("id", "entity_id", cfg$text_cols)
  rec <- rec[, ..keep_cols]

  list(dt = rec, text_cols = cfg$text_cols)
}

# LOADER FOR 10K SYNTHETIC DATASET (D10K)
load_synth10k <- function(cfg) {
  # Load records file
  rec <- fread(cfg$records_path, encoding = "UTF-8")

  # Rename ID column
  if (cfg$id_col %in% names(rec)) {
    setnames(rec, old = cfg$id_col, new = "id")
  } else {
    stop(paste("ID column", cfg$id_col, "not found in", cfg$records_path))
  }

  # Check all text columns exist
  missing_cols <- setdiff(cfg$text_cols, names(rec))
  if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
    cfg$text_cols <- intersect(cfg$text_cols, names(rec))
  }

  if (length(cfg$text_cols) == 0) {
    stop("No valid text columns found in the records file")
  }

  # Load ground truth pairs
  pairs <- fread(cfg$pairs_path, encoding = "UTF-8")
  setnames(pairs,
           old = c(cfg$pair_id1_col, cfg$pair_id2_col),
           new = c("id1", "id2"),
           skip_absent = TRUE)

  # Build entity_id from pairs using connected components
  ent_dt <- entity_from_pairs(rec$id, pairs, "id1", "id2")
  rec[, entity_id := ent_dt$entity_id]

  # Normalize all text columns
  for (col in cfg$text_cols) {
    set(rec, j = col, value = norm_str(rec[[col]]))
  }

  # Keep only necessary columns
  keep_cols <- c("id", "entity_id", cfg$text_cols)
  rec <- rec[, ..keep_cols]

  cat(sprintf("  Text columns being used: %s\n", paste(cfg$text_cols, collapse = ", ")))

  list(dt = rec, text_cols = cfg$text_cols)
}

# ============================================================
# EXPERIMENT RUNNER
# ============================================================

run_one <- function(d0, text_cols, scheme = "record", test_frac = 0.3, seed = 1, neg_ratio = 10,
                    use_blocking = FALSE, block_col = NULL, block_method = "prefix", block_size = 3) {
  sp <- if (scheme == "record") {
    split_record_random(d0, test_frac, seed)
  } else {
    split_entity_disjoint(d0, test_frac, seed)
  }

  tr <- d0[id %in% sp$train_ids]
  te <- d0[id %in% sp$test_ids]

  # Generate pairs with or without blocking
  pairs_tr <- make_pairs_with_blocking(tr$id, tr, use_blocking, block_col, block_method, block_size)
  pairs_te <- make_pairs_with_blocking(te$id, te, use_blocking, block_col, block_method, block_size)

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
                                   base_seed = 2026, neg_ratio = 10,
                                   use_blocking = FALSE, block_col = NULL,
                                   block_method = "prefix", block_size = 3) {
  cat(sprintf("\n=== Running %s (K=%d runs) ===\n", dataset_name, K))

  # Show blocking info if enabled
  if (use_blocking) {
    cat(sprintf("  Using blocking: method=%s, column=%s, size=%d\n",
                block_method, block_col, block_size))

    # Calculate and display blocking statistics
    stats <- blocking_stats(d0, block_col, block_method, block_size)
    cat(sprintf("  Blocking stats: %d records -> %d blocks (avg size: %.1f)\n",
                stats$n_records, stats$n_blocks, stats$avg_block_size))
    cat(sprintf("  Pairs: %.2fM unblocked -> %.2fM blocked (%.1f%% reduction)\n",
                stats$total_pairs_unblocked / 1e6,
                stats$total_pairs_blocked / 1e6,
                (1 - stats$reduction_ratio) * 100))
  }

  out <- list()
  for (scheme in c("record","entity")) {
    cat(sprintf("  Scheme: %s\n", scheme))
    for (k in seq_len(K)) {
      seed <- base_seed + k
      r <- run_one(d0, text_cols, scheme, test_frac, seed, neg_ratio,
                   use_blocking, block_col, block_method, block_size)
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
      if (k %% 5 == 0) cat(sprintf("    Completed run %d/%d\n", k, K))
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
neg_ratio <- 10

res_list <- list()

cat("\n########## STARTING EXPERIMENTS ##########\n")

if (isTRUE(CFG$cora$enabled)) {
  cat("\n>>> Loading CORA dataset...\n")
  cora_obj <- load_cora()
  cat(sprintf("CORA: %d records, %d entities\n",
              nrow(cora_obj$dt),
              length(unique(cora_obj$dt$entity_id))))

  system.time({
    res_list[["cora"]] <- run_experiment_dataset("cora", cora_obj$dt, cora_obj$text_cols,
                                                 K, test_frac, base_seed, neg_ratio)
  })
}

if (isTRUE(CFG$affiliation$enabled)) {
  cat("\n>>> Loading Affiliation dataset...\n")
  aff_obj <- load_affiliation(CFG$affiliation)
  cat(sprintf("Affiliation: %d records, %d entities\n",
              nrow(aff_obj$dt),
              length(unique(aff_obj$dt$entity_id))))

  system.time({
    res_list[["affiliation"]] <- run_experiment_dataset("affiliation", aff_obj$dt, aff_obj$text_cols,
                                                        K, test_frac, base_seed, neg_ratio)
  })
}

if (isTRUE(CFG$affiliation_clean$enabled)) {
  cat("\n>>> Loading Affiliation Clean (multi-field) dataset...\n")
  aff_clean_obj <- load_affiliation_multifield(CFG$affiliation_clean)
  cat(sprintf("Affiliation Clean: %d records, %d entities, %d text fields\n",
              nrow(aff_clean_obj$dt),
              length(unique(aff_clean_obj$dt$entity_id)),
              length(aff_clean_obj$text_cols)))
  cat(sprintf("Text columns: %s\n", paste(aff_clean_obj$text_cols, collapse = ", ")))

  system.time({
    res_list[["affiliation_clean"]] <- run_experiment_dataset("affiliation_clean",
                                                              aff_clean_obj$dt,
                                                              aff_clean_obj$text_cols,
                                                              K, test_frac, base_seed, neg_ratio)
  })
}

if (isTRUE(CFG$synth10k$enabled)) {
  cat("\n>>> Loading Synthetic 10K (D10K) dataset...\n")
  synth_obj <- load_synth10k(CFG$synth10k)
  cat(sprintf("Synth10K: %d records, %d entities, %d text fields\n",
              nrow(synth_obj$dt),
              length(unique(synth_obj$dt$entity_id)),
              length(synth_obj$text_cols)))

  # Get blocking parameters from config
  use_blocking <- isTRUE(CFG$synth10k$use_blocking)
  block_col <- CFG$synth10k$block_on
  block_method <- CFG$synth10k$block_method
  block_size <- CFG$synth10k$block_size

  system.time({
    res_list[["synth10k"]] <- run_experiment_dataset("synth10k",
                                                      synth_obj$dt,
                                                      synth_obj$text_cols,
                                                      K, test_frac, base_seed, neg_ratio,
                                                      use_blocking, block_col,
                                                      block_method, block_size)
  })
}

# ============================================================
# ANALYSIS AND VISUALIZATION
# ============================================================

res_all <- rbindlist(res_list, fill = TRUE)
final_all <- summarize_all(res_all)

cat("\n########## RESULTS ##########\n\n")
print(res_all)
cat("\n")
print(final_all)

# Save results
fwrite(res_all, "er_splitting_results_detailed.csv")
fwrite(final_all, "er_splitting_results_summary.csv")
cat("\nResults saved to CSV files.\n")

# ============================================================
# PLOTTING FUNCTIONS
# ============================================================

plot_er_results <- function(final_all,
                            metric = c("ARI", "NMI", "VI", "FMI", "PairF1"),
                            show_error = TRUE) {

  metric <- match.arg(metric)

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }
  library(ggplot2)

  col_map <- list(
    ARI     = c(mean = "ARI_mean",     sd = "ARI_sd",     ylab = "Adjusted Rand Index (ARI)"),
    NMI     = c(mean = "NMI_mean",     sd = "NMI_sd",     ylab = "Normalized Mutual Information (NMI)"),
    VI      = c(mean = "VI_mean",      sd = "VI_sd",      ylab = "Variation of Information (VI)"),
    FMI     = c(mean = "FMI_mean",     sd = "FMI_sd",     ylab = "Fowlkes–Mallows Index (FMI)"),
    PairF1  = c(mean = "PairF1_mean",  sd = "PairF1_sd",  ylab = "Pairwise F1")
  )

  m <- col_map[[metric]]
  dt <- as.data.table(final_all)

  p <- ggplot(dt, aes(x = scheme, y = get(m["mean"]), fill = scheme)) +
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

plot_thresholds <- function(final_all) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' required for plotting")
  }
  library(ggplot2)

  ggplot(final_all, aes(x = scheme, y = thr_mean, fill = scheme)) +
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

# Generate plots
cat("\n########## GENERATING PLOTS ##########\n")

if (requireNamespace("ggplot2", quietly = TRUE)) {
  pdf("er_splitting_plots.pdf", width = 12, height = 8)

  print(plot_er_results(final_all, metric = "ARI"))
  print(plot_er_results(final_all, metric = "NMI"))
  print(plot_er_results(final_all, metric = "VI"))
  print(plot_er_results(final_all, metric = "FMI"))
  print(plot_er_results(final_all, metric = "PairF1"))
  print(plot_thresholds(final_all))

  dev.off()
  cat("Plots saved to er_splitting_plots.pdf\n")
} else {
  cat("ggplot2 not available - skipping plots\n")
}

cat("\n########## DONE ##########\n")
