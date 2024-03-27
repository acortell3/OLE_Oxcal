# sim_dates = function(n,lower_earliest,upper_latest,carrying_capacity,growth_rate,error) { # create a function with the name my_function
sim_dates = function(n,lower_earliest,upper_latest,carrying_capacity,growth_rate) { # create a function with the name my_function
  calender.dates<<-replicate(n,rLogisticGrowth(a=lower_earliest, b=upper_latest, k=carrying_capacity, r=growth_rate))

  radiocarbon.ages  <<- uncalibrate(calender.dates)$rCRA
  # error.dates <<- rep(error,n) #20 years error for each date
}