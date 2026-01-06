

####### Ploting results

OLE_medians_hol <- readRDS("../Results/OLE_medians_hol.rds")

c1 <- OLE_medians_hol[OLE_medians_hol$sampled_dates == 5 & OLE_medians_hol$r == 0.01 & OLE_medians_hol$Sd == 20,]

plot(x = rep(1,4), y = c(c1$start_date[1],c1$Estimate[1],c1$upperCI[1],c1$lowerCI[1]), pch = 16, col = c("blue","red","orange","orange"), xlim = c(0,51), ylim = c(0,10000))

#points(x = rep(2,4), y = c(c1$start_date[2],c1$Estimate[2],c1$upperCI[2],c1$lowerCI[2]), pch = 16, col = c("blue","red","orange","orange"))
for (i in 2:50){
	points(x = rep(i,4), y = c(c1$start_date[i],c1$Estimate[i],c1$upperCI[i],c1$lowerCI[i]), pch = 16, col = c("blue","red","orange","orange"))
}



plot(x = c1$start_date, y = rep(0,length(c1$start_date)), col = "blue", pch = 16, ylim = c(-600,600))
points(x = c1$start_date, y = c1$start_date - c1$Estimate, col = "red", pch = 16)
points(x = c1$start_date, y = c1$upperCI - c1$start_date, col = "orange", pch = 16)
