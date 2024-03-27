OLE.Median = function(radioages_input, erorate) {

# radioages_input<-sort(replicate(100,rLogisticGrowth(a=5000, b=4000, k=0.01, r=0.01)), decreasing = FALSE)
# erorate <- 20
  
ero.dates <- rep(erorate,length(radioages_input)) 

# Median Approach (accounts for sampling error but not for chronological uncertainty)
input  <- calibrate(radioages_input,ero.dates)

sampledinput<-  (medCal(input) + runif(10)) |> sort() #|> round() #MedCal selects a random date from median of the prob distri. and an unif(10) is added to avert duplicates. Rounding it defeats the purpose of creating unique integers!

sampledinput <- unique(sampledinput)

res1 <- OLE(data.frame(sampledinput,rep(1,length(sampledinput))),alpha=0.05)[2:3] |> as.numeric() |> round()

accuracy <- 0

# Define value
x1 <- 5000   

# Apply between function
range32 <- between(x1, res1[1], res1[2])

if (range32 == "TRUE")
  
{
  accuracy <- accuracy + 1
}

margin <-  res1[2] - res1[1]

distancefromstart <- 0

if (accuracy == 0){
  
  distancefromstart <- min(c(abs(5000 - res1[2]),abs(5000 - res1[1])))
  
}

new_row <- c(margin, accuracy, distancefromstart)

OLE.Median.results[nrow(OLE.Median.results) + 1, ] <<- new_row
}