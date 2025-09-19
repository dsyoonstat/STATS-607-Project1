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
    group_by(across(all_of(c(target, key)))) %>%
    summarise(
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
    group_by(Loss) %>%
    summarise(n = n()) %>%
    mutate(perc = round(n / sum(n) * 100, 2))
  total_row <- data.frame(Loss = "Total", n = sum(InfoLoss$n), perc = 100)
  InfoLoss <- rbind(InfoLoss, total_row)
  
  # Retain only necessary columns and sort
  FinalResult <- FinalResult[, c(target, key, "N.Masked", "Type1", "Type2")]
  FinalResult <- as.data.table(FinalResult)
  setorderv(FinalResult, c(target, key))
  
  cat("Header of aggregated masked table via iLBA\n\n")
  print(head(FinalResult), row.names = FALSE)
  cat("\nDistribution of Information Loss\n")
  print(as.data.frame(InfoLoss), row.names = FALSE)
  
  # Save results
  write.csv(FinalResult, file = output.table.path, row.names = FALSE)
  cat("\nAggregated table saved to \"", output.table.path, "\"\n", sep = "")
  write.csv(InfoLoss, file = output.infoloss.path, row.names = FALSE)
  cat("Information Loss saved to \"", output.infoloss.path, "\"\n", sep = "")
}