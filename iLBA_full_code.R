SCA <- function(freq, B = 3) {
  if (!is.vector(freq)) freq <- as.vector(freq)
  if (!sum(freq < B)) return(freq)
  n <- length(freq)
  changer <- ifelse(stats::runif(n) < freq / B, B, 0)
  return(ifelse(freq < B, changer, freq))
}

iLBA <- function(x, B = 3) {
  nLessEqOrigin <- x[1] # K
  nEqSCA <- x[2]        # k
  SumSCA <- x[3]        # f_i
  Masked <- nEqSCA * B
  type1 <- type2 <- 0
  if (nLessEqOrigin > 1) {
    a <- (SumSCA - 1) %/% B
    Masked <- a * B + trunc(B / 2) + 1
    lbdSCA <- nEqSCA
    ubdSCA <- nEqSCA + nLessEqOrigin * (B - 1)
    lbdSystem <- a * B + 1
    ubdSystem <- (a + 1) * B
    if (a < 0) {
      Masked <- 0
    } else if (lbdSystem < lbdSCA) {
      Masked <- Masked + B
      type1 <- 1
    } else if (ubdSystem > ubdSCA) {
      Masked <- Masked - B
      type2 <- 1
    }
    if (Masked > 0 && Masked < B) Masked <- B
  }
  return(c(Masked, type1, type2))
}

savefulltb <- function(data, hkey, key = NULL, B = 3, rank = NULL, key.threshold = 100, output.path = "fulltable.rds"){
  
  dt <- data.table::as.data.table(data)
  
  # Remove rows with missing values
  na_rows <- dt[, sum(!stats::complete.cases(.SD))]
  if (na_rows > 0) {
    dt <- dt[stats::complete.cases(dt)]
    cat(paste(na_rows, "rows with missing values have been removed.\n"))
  }
  
  # Check that all hkey variables exist in data
  invalid_hkeys <- setdiff(hkey, names(dt))
  if (length(invalid_hkeys) > 0) {
    stop("The following hkey variables do not exist in data: ", paste(invalid_hkeys, collapse = ", "))
  }
  
  # If key is not provided, use all columns except hkey
  if (is.null(key)) {
    key <- setdiff(names(dt), hkey)
  } else {
    # Check that all key variables exist in data
    invalid_keys <- setdiff(key, names(dt))
    if (length(invalid_keys) > 0) {
      stop("The following key variables do not exist in data: ", paste(invalid_keys, collapse = ", "))
    }
  }
  
  # Sort hkey variables by rank if provided
  if (!is.null(rank) && length(rank) != length(hkey)) {
    stop("Length of 'rank' must be equal to length of 'hkey'")
  }
  hkey_ordered <- if (!is.null(rank)) hkey[order(rank)] else hkey
  
  # Remove key variables with too many unique values
  unique_counts <- dt[, sapply(.SD, data.table::uniqueN), .SDcols = key]
  to_remove <- names(unique_counts[unique_counts > key.threshold])
  if (length(to_remove) > 0) {
    dt[, (to_remove) := NULL]
    key <- setdiff(key, to_remove)
    cat("Removed key variables with >", key.threshold, "unique values: ",
        paste(to_remove, collapse = ", "), "\n")
  }
  
  # Create full frequency table with masked counts
  fulltb <- dt[, .(N = .N), by = c(hkey_ordered, key)][, N_SCA := SCA(N, B)]
  
  # Extract sorted unique values for hkey and key variables
  hkey_values <- lapply(dt[, ..hkey_ordered], function(x) sort(unique(x)))
  names(hkey_values) <- hkey_ordered
  key_values <- lapply(dt[, ..key], function(x) sort(unique(x)))
  names(key_values) <- key
  
  # Assemble result list
  result <- list(
    fulltb = fulltb,
    B = B,
    hkey = hkey_ordered,
    hkey_values = hkey_values,
    key = key,
    key_values = key_values
  )
  
  # Save result to RDS
  saveRDS(result, file = output.path)
  cat("Full table saved to", shQuote(output.path), "\n")
  cat("Hierarchical key variables:\n")
  for (i in seq_along(hkey_ordered)) {
    cat(i, ":", hkey_ordered[i], "\n")
  }
  cat("Key variables:", paste(key, collapse = ', '), "\n")
  cat("B:", B, "\n")
  
  invisible(result)
}

saveaggtb <- function(hkey.level, key, input.path = "fulltable.rds", output.table.path = "aggtable.csv", output.infoloss.path = "infoloss.csv") {
  
  # Check if the input file exists
  if (!file.exists(input.path)) {
    stop(paste0("The specified input file does not exist: ", shQuote(input.path)))
  }
  
  # Load saved full table
  fulltable <- readRDS(input.path)
  fulltb <- fulltable$fulltb
  B <- fulltable$B
  hkey_full <- fulltable$hkey
  key_full <- fulltable$key
  
  # Check if hkey.level is valid
  if (hkey.level > length(hkey_full)) {
    stop(paste0("Invalid hkey.level: provided level (", hkey.level,
                ") exceeds the number of hierarchical key variables in the full table (", length(hkey_full), ")."))
  }
  
  
  # Check if key is valid
  invalid_key <- setdiff(key, key_full)
  if (length(invalid_key) > 0) {
    stop("The following key variables do not exist in the full table: ", paste(invalid_key, collapse = ", "))
  }
  
  # Aggregate and summarize counts
  target <- hkey_full[seq_len(hkey.level)]
  summaryData <- fulltb %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(target, key)))) %>%
    dplyr::summarise(
      SumGrOrigin = sum(N[N > B]),
      nLessEqOrigin = sum(N <= B),
      nEqSCA = sum(N_SCA == B),
      SumOrigin = sum(N[N <= B]),
      .groups = "drop"
    )
  
  # Apply iLBA masking
  iLBAdata <- as.data.frame(summaryData)
  nCol <- ncol(iLBAdata)
  iLBACore <- iLBAdata[, (nCol - 2):nCol]
  
  Masked <- t(apply(iLBACore, 1, iLBA, B = B))
  colnames(Masked) <- c("Masked", "Type1", "Type2")
  
  FinalResult <- cbind(iLBAdata, Masked)
  FinalResult$N.Masked <- FinalResult$SumGrOrigin + FinalResult$Masked
  Loss <- FinalResult$N.Masked - FinalResult$SumGrOrigin - FinalResult$SumOrigin
  
  # Summarize information loss
  InfoLoss <- data.frame(Loss) %>%
    dplyr::group_by(Loss) %>%
    dplyr::summarise(n = dplyr::n()) %>%
    dplyr::mutate(perc = round(n / sum(n) * 100, 2))
  total_row <- data.frame(Loss = "Total", n = sum(InfoLoss$n), perc = 100)
  InfoLoss <- rbind(InfoLoss, total_row)
  
  # Retain only necessary columns and sort
  FinalResult <- FinalResult[, c(target, key, "N.Masked", "Type1", "Type2")]
  FinalResult <- data.table::as.data.table(FinalResult)
  data.table::setorderv(FinalResult, c(target, key))
  
  cat("Header of aggregated masked table via iLBA\n\n")
  print(utils::head(FinalResult), row.names = FALSE)
  cat("\nDistribution of Information Loss\n")
  print(as.data.frame(InfoLoss), row.names = FALSE)
  
  # Save results
  utils::write.csv(FinalResult, file = output.table.path, row.names = FALSE)
  cat("\nAggregated table saved to \"", output.table.path, "\"\n", sep = "")
  utils::write.csv(InfoLoss, file = output.infoloss.path, row.names = FALSE)
  cat("Information Loss saved to \"", output.infoloss.path, "\"\n", sep = "")
}

returnaggfreq <- function(hkey.level, key, hkey.value, key.value, input.path = "fulltable.rds") {
  
  # Check if input file exists
  if (!file.exists(input.path)) {
    stop("The specified input file does not exist: ", shQuote(input.path))
  }
  
  # Load saved full table
  fulltable <- readRDS(input.path)
  fulltb <- fulltable$fulltb
  B <- fulltable$B
  hkey_full <- fulltable$hkey
  hkey_values_full <- fulltable$hkey_values
  key_full <- fulltable$key
  key_values_full <- fulltable$key_values
  
  # Validate hkey.level
  if (hkey.level < 1 || hkey.level > length(hkey_full)) {
    stop(sprintf("Invalid hkey.level: %d (must be between 1 and %d)", hkey.level, length(hkey_full)))
  }
  
  # Validate hkey.value length
  if (length(hkey.value) != hkey.level) {
    stop(sprintf("`hkey.value` must be a vector of length %d, one value per hierarchical level", hkey.level))
  }
  
  # Define target hierarchical keys
  target <- hkey_full[seq_len(hkey.level)]
  
  # Validate each hierarchical key value
  for (i in seq_along(target)) {
    key_name <- target[i]
    if (!(hkey.value[i] %in% hkey_values_full[[key_name]])) {
      stop(sprintf("Invalid hkey value: '%s' not found in variable '%s'", hkey.value[i], key_name))
    }
  }
  
  # Validate key variables
  invalid_key <- setdiff(key, key_full)
  if (length(invalid_key) > 0) {
    stop("The following key variables do not exist in the full table: ", paste(invalid_key, collapse = ", "))
  }
  if (length(key.value) != length(key)) {
    stop(sprintf("`key.value` must be length %d (same as `key`)", length(key)))
  }
  for (i in seq_along(key)) {
    if (!(key.value[i] %in% key_values_full[[key[i]]])) {
      stop(sprintf("Invalid key value: '%s' not found in variable '%s'", key.value[i], key[i]))
    }
  }
  
  # Subset matching cell
  condition <- rep(TRUE, nrow(fulltb))
  # apply hierarchical key filters
  for (i in seq_along(target)) {
    condition <- condition & (fulltb[[ target[i] ]] == hkey.value[i])
  }
  # apply other key filters
  for (i in seq_along(key)) {
    condition <- condition & (fulltb[[ key[i] ]] == key.value[i])
  }
  subtb <- fulltb[condition, ]
  
  # Return 0 if no match
  if (nrow(subtb) == 0) {
    message("No matching cell found. Check if hkey.value and key.value are valid.")
    return(invisible(0))
  }
  
  # Compute masked value components
  SumGrOrigin   <- sum(subtb$N[subtb$N > B])
  nLessEqOrigin <- sum(subtb$N <= B)
  nEqSCA        <- sum(subtb$N_SCA == B)
  SumOrigin     <- sum(subtb$N[subtb$N <= B])
  
  # Run iLBA core
  result_vec <- iLBA(c(nLessEqOrigin, nEqSCA, SumOrigin), B = B)
  names(result_vec) <- c("Masked", "Type1", "Type2")
  N.Masked <- SumGrOrigin + result_vec["Masked"]
  
  # Construct result table
  out <- data.frame(
    matrix(c(hkey.value, key.value, N.Masked), nrow = 1, byrow = TRUE),
    stringsAsFactors = FALSE
  )
  colnames(out) <- c(target, key, "N.Masked")
  
  # Output
  cat("N.Masked for specified cell\n")
  cat("Hierarchical keys: ", paste(sprintf("%s = %s", target, hkey.value), collapse = ", "), "\n")
  cat("Key vars: ", paste(sprintf("%s = %s", key, key.value), collapse = ", "), "\n")
  cat("N.Masked: ", N.Masked, "\n\n")
  print(out, row.names = FALSE)
  
  return(invisible(N.Masked))
}

# Load data
load("census_full.rda")
hkey = colnames(census_full)[1:3]
key = colnames(census_full)[4:8]

# Make fulltable
savefulltb(census_full, hkey = hkey, key = key)
saveaggtb(hkey.level = 2, key = key[1:4])
