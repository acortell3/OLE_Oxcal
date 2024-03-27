oxcalcleanerunif = function(input1){
  
  if (length(input1)==3) { 
    
    cleanthreeunif(input1)
    
  } else if (length(input1)==4)
  {
    
    cleanfourunif(input1)
    
  } else if (length(input1)==5)
  {
    
    cleanfiveunif(input1)
    
  } else if (length(input1)==6)
  {
    
    cleansixunif(input1)
    
  } else {
    
    # print(paste0("Uniform Date not accurate or precise. Margin: ", margin, " Accuracy:", range32))
    
    print("Date not accurate")
    
  }
  
}