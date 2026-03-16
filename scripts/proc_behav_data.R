rm(list = ls())

library(dplyr)
library(trimr) # https://cran.r-project.org/web/packages/trimr/vignettes/overview.html

overwrite <- FALSE

path.root <- file.path("/media", "data3", "Joanne_SRT_pw", "data", "behavioral")

out.path.1 <- file.path(path.root, "raw_combined.csv")
out.path.2 <- file.path(path.root, "rt_cleaned.csv")
out.path.3 <- file.path(path.root, "rt_cleaned_medians.csv")

sid.list <- c(1, 2, 4, 6, 8, 9, 10, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 24, 25)

if ( ! file.exists(out.path.3) || overwrite ) {
  
  if ( ! file.exists(out.path.2) || overwrite ) {
    
    if ( ! file.exists(out.path.1) || overwrite ) {
      i <- 0
      DF.list <- list()
      
      for ( sid in sid.list ) {
        for ( run in seq(1, 8) ) {
          i <- i + 1
          
          fn <- sprintf("%03d_run_%d_srtt_prob_new.csv", sid, run)
          fp <- file.path(path.root, "raw", fn)
          DF <- read.csv(fp) %>% 
            dplyr::rename(c(
              "key"      = "keypress", 
              "response" = "p_response_keys", 
              "correct"  = "correct.resp", 
              "rt"       = "p_response_rt"
            )) %>% 
            dplyr::select(all_of(c(
              "participant", "session", "stim_onset", 
              "rule", "rule_group", "cycle", "key", "response", "correct", "rt"
            ))) 
          DF.list[[i]] <- DF %>% 
            dplyr::mutate(across(
              "response", as.factor
            )) 
        }
      }
      
      DF.stacked <- dplyr::bind_rows(DF.list)
      write.csv(DF.stacked, row.names = FALSE, file = out.path.1)
      cat("\nSuccessfully created and saved:\n", out.path.1)
      
    } else {
      cat("\nLoading file from:\n", out.path.1)
      DF.stacked <- read.csv(out.path.1)
    }

    DF.stacked <- DF.stacked %>% 
      dplyr::mutate(across(
        .cols = c("rule", "rule_group", "cycle", "key", "response", "correct"), 
        .fns = as.factor
      )) 
    
    # DF.stacked <- DF.stacked %>% 
    #   dplyr::mutate(
    #     cond = paste(rule, key, sep = "_"), 
    #     .before = "rule"
    #   )
    
    DF.trimmed <- trimr::modifiedRecursive(
        data = DF.stacked,
        minRT = 0, 
        pptVar = "participant", 
        condVar = "rule", # "cond", 
        accVar = "correct",
        rtVar = "rt", 
        returnType = "raw"
      ) %>% 
      dplyr::arrange(
        participant, session, stim_onset
      )
      
    write.csv(DF.trimmed, row.names = FALSE, file = out.path.2)
    cat("\nSuccessfully created and saved:\n", out.path.2)
    
  } else {
    cat("\nLoading file from:\n", out.path.2)
    DF.trimmed <- read.csv(out.path.2)
  }
  
  DF.medians <- DF.trimmed %>%
    group_by(participant, session, cycle, rule) %>%
    summarise(
      rt_medians = median(rt, na.rm = TRUE), 
      .groups = "drop" # ungroup after summarizing
    )
  
  write.csv(DF.medians, row.names = FALSE, file = out.path.3)
  cat("\nSuccessfully created and saved:\n", out.path.3)
  
} else {
  cat("\n", out.path.3, "exists. Nothing is overwrited.")
}
cat("\n\n")