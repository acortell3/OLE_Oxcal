oxcalcleaner = function(input){

if (length(input)==3) { 
  
  cleanthree(input)
  
} else if (length(input)==4)
{
  
  cleanfour(input)
 
} else if (length(input)==5)
{
  
  cleanfive(input)
  
} else if (length(input)==6)
{
  
  cleansix(input)
  
} else {
  
  # print(paste0("Trapezoid Date not accurate or precise. Margin: ", margin, " Accuracy:", range32))
  
  print("Date not accurate")
  
}

}