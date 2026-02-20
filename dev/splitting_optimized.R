# ============================================================
# MEMORY-OPTIMIZED VERSION FOR LARGE DATASETS (e.g., Synth10K)
#
# Key optimizations:
# 1. Chunked pair processing to avoid loading all pairs in memory
# 2. Aggressive blocking with validation
# 3. Memory cleanup between iterations
# 4. Progress monitoring with time estimates
# 5. Automatic checkpointing to save intermediate results
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
  cora = list(enabled = FALSE),  # Disable small datasets when testing large ones

  # 2) Affiliation (original single-column text)
  affiliation = list(enabled = FALSE),

  # 3) Affiliation clean (multi-column)
  affiliation_clean = list(enabled = FALSE),

  # 4) 10K synthetic dataset (D10K) - OPTIMIZED
  synth10k = list(
    enabled = TRUE,
    records_path = "D:/datasplitting/data/10Kfull.csv",
    pairs_path = "D:/datasplitting/data/10Kduplicates.csv",
    id_col = "Id",
    text_cols = c("Aggregate Value", "Embedded Ag.Value", "Clean Ag.Value", "Embedded Clean Ag.Value"),
    pair_id1_col = "Entity1",
    pair_id2_col = "Entity2",

    # AGGRESSIVE BLOCKING (CRITICAL FOR 10K RECORDS)
    use_blocking = TRUE,
    block_on = "Aggregate Value",
    block_method = "prefix",
    block_size = 2,  # START WITH 2 (more aggressive)

    # MEMORY OPTIMIZATION SETTINGS
    max_pairs_per_set = 2000000,  # Maximum pairs to process (2M)
    chunk_size = 50000,           # Process features in chunks
    early_stopping = TRUE         # Skip run if pairs exceed limit
  )
)

# ============================================================
# Shared utilities
# ============================================================

norm_str <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- tolower(x)
  x <- gsub("[^a-z0-9]+", " ", x)
  x <- gsub("\\s+", " ", x)
  trimws(x)
}

VI_manual <- function(x, y, eps = 1e-15) {
  x <- as.integer(factor(x))
  y <- as.integer(factor(y))
  n <- length(x)
  tab <- table(x, y)
  pxy <- tab / n
  px <- rowSums(pxy)
  py <- colSums(pxy)
  Hx <- -sum(px * log(px + eps))
  Hy <- -sum(py * log(py + eps))
  nz <- pxy > 0
  Ixy <- sum(pxy[nz] * log(pxy[nz] / (px[row(pxy)[nz]] * py[col(pxy)[nz]] + eps) + eps))
  Hx + Hy - 2 * Ixy
}

FMI <- function(true_labels, pred_labels) {
  tab <- table(true_labels, pred_labels)
  tp <- sum(choose(tab, 2))
  fp <- sum(choose(colSums(tab), 2)) - tp
  fn <- sum(choose(rowSums(tab), 2)) - tp
  if (tp == 0) return(0)
  sqrt((tp / (tp + fp)) * (tp / (tp + fn)))
}

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

make_all_pairs <- function(ids) {
  ids <- sort(unique(as.integer(ids)))
  n <- length(ids)
  if (n < 2) return(data.table(id1 = integer(), id2 = integer()))
  id1 <- rep(ids[1:(n - 1)], times = (n - 1):1)
  id2 <- unlist(lapply(2:n, function(i) ids[i:n]), use.names = FALSE)
  data.table(id1 = id1, id2 = id2)
}

# ============================================================
# OPTIMIZED BLOCKING FUNCTIONS
# ============================================================

generate_block_key <- function(text, method = "prefix", size = 3) {
  if (method == "prefix") {
    return(substr(text, 1, size))
  } else if (method == "soundex") {
    if (!requireNamespace("phonics", quietly = TRUE)) {
      warning("phonics package not installed. Falling back to prefix blocking.")
      return(substr(text, 1, size))
    }
    return(phonics::soundex(text))
  } else if (method == "ngram") {
    if (nchar(text) < size) return(text)
    return(substr(text, 1, size))
  } else {
    stop(paste("Unknown blocking method:", method))
  }
}

# OPTIMIZED: Process blocks in batches to avoid memory issues
make_blocked_pairs <- function(ids, d, block_col, method = "prefix", size = 3, max_pairs = 2000000) {
  setkey(d, id)
  d_subset <- d[id %in% ids, .(id, text = get(block_col))]
  d_subset[, block_key := generate_block_key(text, method, size)]

  # Count pairs per block
  block_counts <- d_subset[, .N, by = block_key]
  block_counts[, n_pairs := choose(N, 2)]
  block_counts[is.nan(n_pairs), n_pairs := 0]

  total_pairs_estimate <- sum(block_counts$n_pairs)

  cat(sprintf("    Blocking created %d blocks, estimated %s pairs\n",
              nrow(block_counts),
              format(total_pairs_estimate, big.mark = ",")))

  if (total_pairs_estimate > max_pairs) {
    cat(sprintf("    WARNING: Estimated pairs (%s) exceeds max_pairs (%s)\n",
                format(total_pairs_estimate, big.mark = ","),
                format(max_pairs, big.mark = ",")))
    cat("    Try: (1) smaller block_size, (2) different block_on column, (3) increase max_pairs_per_set\n")
    return(NULL)  # Signal to skip this run
  }

  # Generate pairs within each block
  pairs_list <- list()
  for (key in block_counts$block_key) {
    block_ids <- d_subset[block_key == key, id]
    if (length(block_ids) >= 2) {
      pairs_list[[length(pairs_list) + 1]] <- make_all_pairs(block_ids)
    }
  }

  if (length(pairs_list) == 0) {
    return(data.table(id1 = integer(), id2 = integer()))
  }

  all_pairs <- rbindlist(pairs_list)
  all_pairs <- unique(all_pairs)

  # Final safety check
  if (nrow(all_pairs) > max_pairs) {
    cat(sprintf("    ERROR: Actual pairs (%s) exceeds max_pairs (%s). Aborting.\n",
                format(nrow(all_pairs), big.mark = ","),
                format(max_pairs, big.mark = ",")))
    return(NULL)
  }

  return(all_pairs)
}

make_pairs_with_blocking <- function(ids, d = NULL, use_blocking = FALSE,
                                     block_col = NULL, method = "prefix",
                                     size = 3, max_pairs = 2000000) {
  if (!use_blocking || is.null(d) || is.null(block_col)) {
    all_pairs <- make_all_pairs(ids)
    if (nrow(all_pairs) > max_pairs) {
      cat(sprintf("    ERROR: Without blocking, pairs (%s) exceeds max_pairs (%s)\n",
                  format(nrow(all_pairs), big.mark = ","),
                  format(max_pairs, big.mark = ",")))
      cat("    You MUST enable blocking for this dataset size!\n")
      return(NULL)
    }
    return(all_pairs)
  } else {
    return(make_blocked_pairs(ids, d, block_col, method, size, max_pairs))
  }
}

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

  block_sizes <- table(d_copy$block_key)
  pairs_per_block <- choose(block_sizes, 2)
  pairs_per_block[is.nan(pairs_per_block)] <- 0

  stats$total_pairs_blocked <- sum(pairs_per_block)
  stats$total_pairs_unblocked <- choose(nrow(d), 2)
  stats$reduction_ratio <- stats$total_pairs_blocked / stats$total_pairs_unblocked

  return(stats)
}

# ============================================================
# SPLITS
# ============================================================

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

# ============================================================
# CHUNKED FEATURE GENERATION (MEMORY EFFICIENT)
# ============================================================

pair_features_text_chunked <- function(d, pairs, text_cols, chunk_size = 50000) {
  n_pairs <- nrow(pairs)
  n_chunks <- ceiling(n_pairs / chunk_size)

  cat(sprintf("    Computing features for %s pairs in %d chunks...\n",
              format(n_pairs, big.mark = ","), n_chunks))

  setkey(d, id)

  X_chunks <- list()
  y_chunks <- list()

  for (i in 1:n_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_pairs)
    pairs_chunk <- pairs[start_idx:end_idx]

    a <- d[pairs_chunk$id1]
    b <- d[pairs_chunk$id2]

    X_list <- list()
    for (col in text_cols) {
      xa <- a[[col]]; xb <- b[[col]]
      jw <- 1 - stringdist(xa, xb, method = "jw", p = 0.1)
      lv_d <- stringdist(xa, xb, method = "lv")
      lv_m <- pmax(nchar(xa), nchar(xb))
      lv <- 1 - ifelse(lv_m == 0, 0, lv_d / lv_m)

      X_list[[paste0(col, "_jw")]] <- jw
      X_list[[paste0(col, "_lv")]] <- lv
    }

    X_chunks[[i]] <- as.data.table(X_list)
    y_chunks[[i]] <- as.integer(a$entity_id == b$entity_id)

    if (i %% 10 == 0 || i == n_chunks) {
      cat(sprintf("      Chunk %d/%d complete\n", i, n_chunks))
    }

    # Clean up
    rm(a, b, X_list, pairs_chunk)
    gc(verbose = FALSE)
  }

  X <- rbindlist(X_chunks)
  y <- unlist(y_chunks)

  list(X = X, y = y)
}

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

load_synth10k <- function(cfg) {
  rec <- fread(cfg$records_path, encoding = "UTF-8")

  if (cfg$id_col %in% names(rec)) {
    setnames(rec, old = cfg$id_col, new = "id")
  } else {
    stop(paste("ID column", cfg$id_col, "not found in", cfg$records_path))
  }

  missing_cols <- setdiff(cfg$text_cols, names(rec))
  if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
    cfg$text_cols <- intersect(cfg$text_cols, names(rec))
  }

  if (length(cfg$text_cols) == 0) {
    stop("No valid text columns found in the records file")
  }

  pairs <- fread(cfg$pairs_path, encoding = "UTF-8")
  setnames(pairs,
           old = c(cfg$pair_id1_col, cfg$pair_id2_col),
           new = c("id1", "id2"),
           skip_absent = TRUE)

  ent_dt <- entity_from_pairs(rec$id, pairs, "id1", "id2")
  rec[, entity_id := ent_dt$entity_id]

  for (col in cfg$text_cols) {
    set(rec, j = col, value = norm_str(rec[[col]]))
  }

  keep_cols <- c("id", "entity_id", cfg$text_cols)
  rec <- rec[, ..keep_cols]

  cat(sprintf("  Text columns: %s\n", paste(cfg$text_cols, collapse = ", ")))

  list(dt = rec, text_cols = cfg$text_cols)
}

# ============================================================
# OPTIMIZED EXPERIMENT RUNNER
# ============================================================

run_one_optimized <- function(d0, text_cols, scheme = "record", test_frac = 0.3,
                              seed = 1, neg_ratio = 10, use_blocking = FALSE,
                              block_col = NULL, block_method = "prefix",
                              block_size = 3, max_pairs = 2000000,
                              chunk_size = 50000, early_stopping = TRUE) {

  cat(sprintf("  Run with scheme=%s, seed=%d\n", scheme, seed))

  sp <- if (scheme == "record") {
    split_record_random(d0, test_frac, seed)
  } else {
    split_entity_disjoint(d0, test_frac, seed)
  }

  tr <- d0[id %in% sp$train_ids]
  te <- d0[id %in% sp$test_ids]

  cat(sprintf("    Train: %d records, Test: %d records\n", nrow(tr), nrow(te)))

  # Generate pairs with safety checks
  pairs_tr <- make_pairs_with_blocking(tr$id, tr, use_blocking, block_col,
                                       block_method, block_size, max_pairs)
  if (is.null(pairs_tr) && early_stopping) {
    cat("    SKIPPING: Train pairs exceed limit\n")
    return(NULL)
  }

  pairs_te <- make_pairs_with_blocking(te$id, te, use_blocking, block_col,
                                       block_method, block_size, max_pairs)
  if (is.null(pairs_te) && early_stopping) {
    cat("    SKIPPING: Test pairs exceed limit\n")
    return(NULL)
  }

  if (nrow(pairs_tr) == 0 || nrow(pairs_te) == 0) {
    cat("    SKIPPING: No pairs generated\n")
    return(NULL)
  }

  # Chunked feature computation
  ftr <- pair_features_text_chunked(tr, pairs_tr, text_cols, chunk_size)
  fte <- pair_features_text_chunked(te, pairs_te, text_cols, chunk_size)

  Xtr <- ftr$X; ytr <- ftr$y
  Xte <- fte$X; yte <- fte$y

  cat("    Training model...\n")

  mu <- fit_imputer(Xtr)
  Xtr <- apply_imputer(Xtr, mu)
  Xte <- apply_imputer(Xte, mu)

  sub <- subsample_train(Xtr, ytr, neg_ratio = neg_ratio, seed = seed)
  Xtr_s <- sub$X; ytr_s <- sub$y

  pp <- fit_predict_glmnet(Xtr_s, ytr_s, Xte, seed = seed)
  th <- best_threshold_f1(pp$p_train, ytr_s)$threshold

  yhat_te <- as.integer(pp$p_test >= th)
  f1_te <- pair_f1(yhat_te, yte)

  cat("    Computing cluster metrics...\n")

  edges_te <- pairs_te[yhat_te == 1, .(id1, id2)]
  rec_ids <- te$id
  pred_cluster <- cluster_from_edges(rec_ids, edges_te)
  true_cluster <- as.integer(factor(te[match(rec_ids, id), entity_id]))

  # Save metrics before cleanup
  n_train_pairs <- length(ytr)
  n_test_pairs <- length(yte)
  pos_rate_train <- mean(ytr)
  pos_rate_test <- mean(yte)

  # Memory cleanup
  rm(Xtr, Xte, ytr, yte, Xtr_s, ytr_s, ftr, fte, pairs_tr, pairs_te)
  gc(verbose = FALSE)

  list(
    ARI = mclust::adjustedRandIndex(true_cluster, pred_cluster),
    NMI = aricode::NMI(true_cluster, pred_cluster),
    VI  = VI_manual(true_cluster, pred_cluster),
    FMI  = FMI(true_cluster, pred_cluster),
    PairF1 = f1_te,
    threshold = th,
    n_train_records = nrow(tr),
    n_test_records  = nrow(te),
    n_train_pairs = n_train_pairs,
    n_test_pairs = n_test_pairs,
    pos_rate_train = pos_rate_train,
    pos_rate_test = pos_rate_test
  )
}

run_experiment_dataset_optimized <- function(dataset_name, d0, text_cols, K = 10,
                                             test_frac = 0.30, base_seed = 2026,
                                             neg_ratio = 10, use_blocking = FALSE,
                                             block_col = NULL, block_method = "prefix",
                                             block_size = 3, max_pairs = 2000000,
                                             chunk_size = 50000, early_stopping = TRUE) {

  cat(sprintf("\n========================================\n"))
  cat(sprintf("DATASET: %s (K=%d runs)\n", dataset_name, K))
  cat(sprintf("========================================\n"))

  if (use_blocking) {
    cat(sprintf("Blocking: method=%s, column=%s, size=%d\n",
                block_method, block_col, block_size))

    stats <- blocking_stats(d0, block_col, block_method, block_size)
    cat(sprintf("Blocking stats:\n"))
    cat(sprintf("  Records: %s\n", format(stats$n_records, big.mark = ",")))
    cat(sprintf("  Blocks: %s (avg size: %.1f, max: %d)\n",
                format(stats$n_blocks, big.mark = ","),
                stats$avg_block_size, stats$max_block_size))
    cat(sprintf("  Pairs: %s unblocked -> %s blocked (%.1f%% reduction)\n",
                format(stats$total_pairs_unblocked, big.mark = ","),
                format(stats$total_pairs_blocked, big.mark = ","),
                (1 - stats$reduction_ratio) * 100))

    if (stats$total_pairs_blocked > max_pairs) {
      cat(sprintf("\nWARNING: Estimated blocked pairs (%s) > max_pairs (%s)\n",
                  format(stats$total_pairs_blocked, big.mark = ","),
                  format(max_pairs, big.mark = ",")))
      cat("Consider: (1) smaller block_size, (2) increase max_pairs_per_set\n")
    }
  }

  out <- list()
  start_time <- Sys.time()

  for (scheme in c("record", "entity")) {
    cat(sprintf("\n--- Scheme: %s ---\n", scheme))

    for (k in seq_len(K)) {
      cat(sprintf("\n[Run %d/%d]\n", k, K))
      seed <- base_seed + k

      run_start <- Sys.time()

      r <- run_one_optimized(d0, text_cols, scheme, test_frac, seed, neg_ratio,
                             use_blocking, block_col, block_method, block_size,
                             max_pairs, chunk_size, early_stopping)

      run_elapsed <- difftime(Sys.time(), run_start, units = "secs")

      if (is.null(r)) {
        cat(sprintf("  Run %d SKIPPED (exceeded limits)\n", k))
        next
      }

      cat(sprintf("  Run %d complete in %.1f sec\n", k, run_elapsed))
      cat(sprintf("  Metrics: ARI=%.3f, NMI=%.3f, PairF1=%.3f\n", r$ARI, r$NMI, r$PairF1))

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

      # Save checkpoint after each run
      if (length(out) > 0) {
        checkpoint <- rbindlist(out)
        fwrite(checkpoint, sprintf("%s_checkpoint.csv", dataset_name))
      }

      # Estimate remaining time
      elapsed <- difftime(Sys.time(), start_time, units = "mins")
      per_run <- as.numeric(elapsed) / length(out)
      remaining_runs <- K * 2 - length(out)  # 2 schemes
      est_remaining <- per_run * remaining_runs

      cat(sprintf("  Estimated time remaining: %.1f minutes\n", est_remaining))
    }
  }

  total_elapsed <- difftime(Sys.time(), start_time, units = "mins")
  cat(sprintf("\nTotal time: %.1f minutes\n", total_elapsed))

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
# RUN EXPERIMENTS
# ============================================================

K <- 10
test_frac <- 0.30
base_seed <- 2026
neg_ratio <- 10

res_list <- list()

cat("\n" ,"#" , rep("=", 50), "#\n", sep = "")
cat("STARTING OPTIMIZED EXPERIMENTS\n")
cat("#", rep("=", 50), "#\n\n", sep = "")

if (isTRUE(CFG$synth10k$enabled)) {
  cat(">>> Loading Synthetic 10K (D10K) dataset...\n")
  synth_obj <- load_synth10k(CFG$synth10k)
  cat(sprintf("Loaded: %d records, %d entities, %d text fields\n",
              nrow(synth_obj$dt),
              length(unique(synth_obj$dt$entity_id)),
              length(synth_obj$text_cols)))

  res_list[["synth10k"]] <- run_experiment_dataset_optimized(
    "synth10k",
    synth_obj$dt,
    synth_obj$text_cols,
    K = K,
    test_frac = test_frac,
    base_seed = base_seed,
    neg_ratio = neg_ratio,
    use_blocking = CFG$synth10k$use_blocking,
    block_col = CFG$synth10k$block_on,
    block_method = CFG$synth10k$block_method,
    block_size = CFG$synth10k$block_size,
    max_pairs = CFG$synth10k$max_pairs_per_set,
    chunk_size = CFG$synth10k$chunk_size,
    early_stopping = CFG$synth10k$early_stopping
  )
}

# ============================================================
# SAVE RESULTS
# ============================================================

if (length(res_list) > 0) {
  res_all <- rbindlist(res_list, fill = TRUE)
  final_all <- summarize_all(res_all)

  cat("\n========================================\n")
  cat("FINAL RESULTS\n")
  cat("========================================\n\n")

  print(final_all)

  fwrite(res_all, "D:/datasplitting/results/synth10k_results_detailed.csv")
  fwrite(final_all, "D:/datasplitting/results/synth10k_results_summary.csv")

  cat("\n✓ Results saved to CSV files\n")
  cat("  - synth10k_results_detailed.csv\n")
  cat("  - synth10k_results_summary.csv\n")
  cat("  - synth10k_checkpoint.csv (intermediate saves)\n")
} else {
  cat("\nNo results generated. Check configuration and warnings above.\n")
}

cat("\n========================================\n")
cat("DONE!\n")
cat("========================================\n")
