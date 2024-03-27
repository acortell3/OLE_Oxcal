cleanfiveunif = function(input){
  
  # Bifurcated approach for 95.4% Confidence interval of Start posterior phase distribution is split into two parts
  rng3 <- input[3]
  rng4 <- input[4]
  rng5 <- input[5]
  
  # Cleaning of the output begins
  
  rngsplit3<-unlist(strsplit(rng3, "\"    ", fixed = TRUE))
  rngsplit4<-unlist(strsplit(rng4, "\"    ", fixed = TRUE))
  rngsplit5<-unlist(strsplit(rng5, "\"    ", fixed = TRUE))
  
  keeps <- c(".")
  
  rngrange3 <- gsub(paste0(".*?($|'|", paste(paste0("\\", 
                                                    keeps), collapse = "|"), "|[^[:punct:]]).*?"), "\\1", rngsplit3[2])
  rngrange4 <- gsub(paste0(".*?($|'|", paste(paste0("\\", 
                                                    keeps), collapse = "|"), "|[^[:punct:]]).*?"), "\\1", rngsplit4[2])
  rngrange5 <- gsub(paste0(".*?($|'|", paste(paste0("\\", 
                                                    keeps), collapse = "|"), "|[^[:punct:]]).*?"), "\\1", rngsplit5[2])
  
  rngrangesplit3 <- gsub('BP','',rngrange3)
  rngrangesplit4 <- gsub('BP','',rngrange4)
  rngrangesplit5 <- gsub('BP','',rngrange5)
  
  rngrangesplit32 <- strsplit(gsub(paste0("([[:alnum:]]{",    # Apply strsplit function
                                          6,
                                          "})"),
                                   "\\1 ",
                                   rngrangesplit3),
                              " ")[[1]]
  rngrangesplit42 <- strsplit(gsub(paste0("([[:alnum:]]{",    # Apply strsplit function
                                          6,
                                          "})"),
                                   "\\1 ",
                                   rngrangesplit4),
                              
                              " ")[[1]]
  rngrangesplit52 <- strsplit(gsub(paste0("([[:alnum:]]{",    # Apply strsplit function
                                          6,
                                          "})"),
                                   "\\1 ",
                                   rngrangesplit5),
                              
                              " ")[[1]]
  rngrangesplit32_num <- as.numeric(rngrangesplit32)
  
  rngrangesplit32_num_noNA<-rngrangesplit32_num[!is.na(rngrangesplit32_num)]
  
  rngrangesplit42_num <- as.numeric(rngrangesplit42)
  
  rngrangesplit42_num_noNA<-rngrangesplit42_num[!is.na(rngrangesplit42_num)]
 
  rngrangesplit52_num <- as.numeric(rngrangesplit52)
  
  rngrangesplit52_num_noNA<-rngrangesplit52_num[!is.na(rngrangesplit52_num)]
  
  # Cleaned entires 
  
  # Filtering accurate results
  
  time.s32 <- rngrangesplit32_num_noNA[1]
  time.e32 <- rngrangesplit32_num_noNA[3]
  
  accuracy <- 0
  
  # Define value
  x1 <- 5000   
  
  # Define lower bound                 
  left1 <- time.e32  
  
  # Define upper bound                  
  right1 <- time.s32  
  
  # Apply between function
  range32 <- between(x1, left1, right1) 
  
  time.s42 <- rngrangesplit42_num_noNA[1]
  time.e42 <- rngrangesplit42_num_noNA[3]
  
  # Define value
  x2 <- 5000   
  
  # Define lower bound                 
  left2 <- time.e42  
  
  # Define upper bound                  
  right2 <- time.s42  
  
  # Apply between function
  range42 <- between(x2, left2, right2) 
  
  
  time.s52 <- rngrangesplit52_num_noNA[1]
  time.e52 <- rngrangesplit52_num_noNA[3]
  
  # Define value
  x3 <- 5000   
  
  # Define lower bound                 
  left3 <- time.e52  
  
  # Define upper bound                  
  right3 <- time.s52  
  
  # Apply between function
  range52 <- between(x2, left3, right3) 
  
  if (range32 == "TRUE" || range42 == "TRUE"|| range52 == "TRUE")
  {
    accuracy <- accuracy + 1
  }
  margin32 <-  time.s32 - time.e32
  margin42 <- time.s42 - time.e42
  margin52 <- time.s52 - time.e52
  
  margin <- margin52 + margin42 + margin32

distancefromstart <- 0
  
if (accuracy == 0){
  
  distancefromstart <- min(c(abs(5000-rngrangesplit32_num_noNA[1]),abs(5000-rngrangesplit32_num_noNA[3]),abs(5000-rngrangesplit42_num_noNA[1]),abs(5000-rngrangesplit42_num_noNA[3]),abs(5000-rngrangesplit52_num_noNA[1]),abs(5000-rngrangesplit52_num_noNA[3])))
  
}
  
  new_row <- c(margin, accuracy, distancefromstart)
  
  unif_df[nrow(unif_df) + 1, ] <<- new_row 
}
