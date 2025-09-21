#' Save aggregated masked table and distribution of information loss
#'
#' Aggregates and masks counts from a previously saved full frequency table
#' via iLBA algorithm, computes the distribution of information loss,
#' and saves the results as a CSV file.
#'
#' @param hkey_level Integer indicating the level of hierarchical key at which
#' the aggregation is performed (1 indicates the top level).
#' @param key Character vector of key variable names used in the aggregated table.
#' @param input_path String path to load the RDS object
#' produced by `save_full_tb()` (default = "full_table.rds").
#' @param output_table_path String path to save the aggregated masked table
#' in CSV format (default = "agg_table.csv").
#' @param output_infoloss_path String path to save the information loss
#' distribution in CSV format (default = "infoloss.csv").
#'
#' @return
#' Saves two CSV files to the specified file paths:
#' (1) the aggregated masked table and
#' (2) the distribution of information loss.
#'

save_agg_tb <- function(hkey_level, key, input_path = "full_table.rds",
                      output_table_path = "agg_table.csv",
                      output_infoloss_path = "infoloss.csv") {
  
  # Check if the input file exists
  if (!file.exists(input_path)) {
    stop(paste0("The specified input file does not exist: ", shQuote(input_path)))
  }
  
  
  # Load saved full table
  fulltable <- readRDS(input_path)
  fulltb <- fulltable$fulltb
  mask_threshold <- fulltable$mask_threshold
  hkey_full <- fulltable$hkey
  key_full <- fulltable$key
  
  # Check if hkey_level is valid
  if (hkey_level > length(hkey_full)) {
    stop(paste0("Invalid hkey_level: provided level (", hkey_level,
                ") exceeds the number of hierarchical key variables in the full table (",
                length(hkey_full), ")."))
  }
  
  
  # Check if key is valid
  invalid_key <- setdiff(key, key_full)
  if (length(invalid_key) > 0) {
    stop("The following key variables do not exist in the full table: ",
         paste(invalid_key, collapse = ", "))
  }
  
  
  # Aggregate and summarize counts
  target <- hkey_full[seq_len(hkey_level)]
  summaryData <- fulltb %>%
    group_by(across(all_of(c(target, key)))) %>%
    summarise(
      SumGrOrigin = sum(N[N > mask_threshold]),
      nLessEqOrigin = sum(N <= mask_threshold),
      nEqSCA = sum(N_SCA == mask_threshold),
      SumOrigin = sum(N[N <= mask_threshold]),
      .groups = "drop"
    )
  
  
  # Apply iLBA masking
  iLBAdata <- as.data.frame(summaryData)
  nCol <- ncol(iLBAdata)
  iLBACore <- iLBAdata[, (nCol - 2):nCol]
  
  Masked <- t(apply(iLBACore, 1, apply_iLBA, mask_threshold = mask_threshold))
  colnames(Masked) <- c("Masked", "Type1", "Type2")
  
  FinalResult <- cbind(iLBAdata, Masked)
  FinalResult$N.Masked <- FinalResult$SumGrOrigin + FinalResult$Masked
  Loss <- FinalResult$N.Masked - FinalResult$SumGrOrigin - FinalResult$SumOrigin
  
  
  # Summarize information loss
  InfoLoss <- data.frame(Loss) %>%
    group_by(Loss) %>%
    summarise(n = n()) %>%
    mutate(perc = round(n / sum(n) * 100, 2))
  total_row <- data.frame(Loss = "Total", n = sum(InfoLoss$n), perc = 100)
  InfoLoss <- rbind(InfoLoss, total_row)
  
  
  # Retain only necessary columns and sort
  FinalResult <- FinalResult[, c(target, key, "N.Masked", "Type1", "Type2")]
  FinalResult <- as.data.table(FinalResult)
  setorderv(FinalResult, c(target, key))
  
  
  # Print summary of the results
  cat("Header of aggregated masked table via iLBA\n\n")
  print(head(FinalResult), row.names = FALSE)
  cat("\nDistribution of Information Loss\n")
  print(as.data.frame(InfoLoss), row.names = FALSE)
  
  
  # Save results
  write.csv(FinalResult, file = output_table_path, row.names = FALSE)
  cat("\nAggregated table saved to \"", output_table_path, "\"\n", sep = "")
  write.csv(InfoLoss, file = output_infoloss_path, row.names = FALSE)
  cat("Information Loss saved to \"", output_infoloss_path, "\"\n", sep = "")
}
