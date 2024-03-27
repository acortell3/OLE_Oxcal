cleanthreeunif <- function (input){
  
  # print(input)
  
  input <- input[3]
  
  rngsplit3<-unlist(strsplit(input, "\"    ", fixed = TRUE))
  
  keeps <- c(".")
  
  rngrange3 <- gsub(paste0(".*?($|'|", paste(paste0("\\", 
                                                    keeps), collapse = "|"), "|[^[:punct:]]).*?"), "\\1", rngsplit3[2])
  
  rngrangesplit3 <- gsub('BP','',rngrange3)
  # print(rngrangesplit3)
  
  rngrangesplit32 <- strsplit(gsub(paste0("([[:alnum:]]{",    # Apply strsplit function
                                          6,
                                          "})"),
                                   "\\1 ",
                                   rngrangesplit3),
                              " ")[[1]]
  # print(rngrangesplit32)
  
  rngrangesplit32_num <- as.numeric(rngrangesplit32)
  
  rngrangesplit32_num_noNA<-rngrangesplit32_num[!is.na(rngrangesplit32_num)]
  
  time.s32 <- rngrangesplit32_num_noNA[1]
  time.e32 <- rngrangesplit32_num_noNA[3]
  
  # print(time.e32)
  
  accuracy <- 0
  
  # Define value
  x1 <- 5000   
  
  # Define lower bound                 
  left1 <- time.e32  
  
  # Define upper bound                  
  right1 <- time.s32  
  
  # Apply between function
  range32 <- between(x1, left1, right1) 
  
  # print(range32)
  
  if (range32 == "TRUE")
    
  {
    accuracy <- accuracy + 1
  }
    margin <-  time.s32 - time.e32

distancefromstart <- 0
    
if (accuracy == 0){
  
  distancefromstart <- min(c(abs(5000 - rngrangesplit32_num_noNA[1]),abs(5000 - rngrangesplit32_num_noNA[3])))
  
}

    new_row <- c(margin, accuracy, distancefromstart)
    
    unif_df[nrow(unif_df) + 1, ] <<- new_row
}

