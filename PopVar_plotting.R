setwd("C:/Users/Jeff/Documents/University of Minnesota/Barley Lab/Projects/Genomic Selection/PopVar/Files")

load("PopVar_allspatial_out.RData")

Rel_eff <- read.table("../../Cross-Validation/Files/Relative_efficiencies.txt", header = TRUE)
Rel_eff <- Rel_eff[,-c(13:15)]
dim(Rel_eff)

#Vector of names from PopVar
names <- unique(unlist(strsplit(names(Popvar.out_subdf), c("_RAW","_M1","_M3","_MVA"))))

#Gather the means of the important data
Pred_val <- matrix(ncol = length(Popvar.out_subdf), nrow = 3)
colnames(Pred_val) <- unlist(strsplit(names(Popvar.out_subdf), c("_RAW","_M1","_M3","_MVA")))
rownames(Pred_val) <- colnames(Popvar.out_subdf[[1]][,c(7,9:10)])
for(i in 1:length(Popvar.out_subdf)) {
  for(r in 1:dim(Pred_val)[1]) {
    Pred_val[r,i] <- mean(unlist(Popvar.out_subdf[[i]][,(c(7,9,10)[r])]))
  }
}

pdf(file = "PopVar_plots.pdf")
#Plot loop
for(i in 1:dim(Pred_val)[1]){
  mp <- as.numeric(barplot(Pred_val[i,], space = rep(c(1,0.5,0.5,0.5),9)))[c(2:4,6:8,10:12,14:16,18:20,22:24,26:28,30:32,34:36)]
  barplot(Pred_val[i,], col = c("red", "lightblue", "yellow", "green"), axisnames = FALSE, ylim = range(pretty(range(Pred_val[i,]))), ylab = rownames(Pred_val)[1], main = paste("Predicted",rownames(Pred_val)[1], "by Spatial Adjustment", sep = " "), names.arg = NULL, space = rep(c(1,0.5,0.5,0.5),9))
  par(mar = c(6,4,4,4))
  legend("topleft", legend = c("Raw Data", "Method 1", "Method 3", "Moving Average"), fill = c("red", "lightblue", "yellow", "green"), cex = 0.75)
  axis(side = 1,las = 2, labels = names, at = seq(4,56, by = 6.5), cex.axis = 0.8)
  points(mp,Rel_eff[i,]/100, pch = 15)
  axis(side = 4, at = pretty(range(Rel_eff[1,]/100)))
  mtext("Relative Efficiency (/100)", side = 4, line = 3)
}
dev.off()
