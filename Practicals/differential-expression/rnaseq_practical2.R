# Install the library
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("maSigPro")

#install.packages("pheatmap")
#install.packages("MASS")
library("pheatmap")
library("MASS")

library("maSigPro")

# Load the accompanying data
data(data.abiotic)
data(edesign.abiotic)

design <- make.design.matrix(edesign.abiotic, degree = 2)

# Q1: Try also building the design matrix with a degree of 1 and see the difference
design.degree1 <- make.design.matrix(edesign.abiotic, degree = 1)
design.degree2 <- make.design.matrix(edesign.abiotic, degree = 2)
design.degree1[["dis"]]
design.degree2[["dis"]]
# Degree 2 adds an extra set of columns where a quadratic relationship between
# time and cold / heat / salt is added.

fit <- p.vector(data.abiotic, design, Q = 0.05, MT.adjust = "BH", min.obs = 20)

# Q2: How many genes show significance for the specified threshold Q = 0.05?
sum(fit$p.adjusted < 0.05) # 398

# Q3: And for Q = 0.005?
sum(fit$p.adjusted < 0.005) # 242
  
# Q4: What is the effect on the overall amount of significant genes when no 
# multiple testing correction is applied? Stuck? Try reading the help using 
# ?p.vector
p.adj.nsign <- round((sum(fit$p.adjusted < 0.05) / length(fit$p.adjusted))*100, digits = 4)
paste0("The percentage of adjusted p-values being significant is ", p.adj.nsign, " %")
p.nsign <- round((sum(fit$p.vector < 0.05) / length(fit$p.vector))*100, digits = 4)
paste0("The percentage of p-values being significant is ", p.nsign, " %")

tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)

# Q5: What would be the effect of increasing alfa?
# The alfa argument defines the significance level used for variable selection 
# in the stepwise regression. If we increase alfa than a test will be
# significant more easily and therefore more variables will be included in 
# the model.

sigs <- get.siggenes(tstep, rsq = 0.6, vars = "groups")

suma2Venn(sigs$summary[, c(2:4)])

# Q6: How many genes show significance between all of the timepoints but also 
# the experimental groups?
# 101

# Q7: What percentage of the genes is unique to a single experimental 
# condition?
# 1 + 9 + 1 = 11, total = 1 + 26 + 2 + 101 + 9 + 87 + 1 = 227
prcnt.unique.genes <- round((11 / 227)*100, digits = 4)
paste0("The percentage of unique genes is ", prcnt.unique.genes, " %")


pdf(width=20, height=10, file="SeeGenesOutput.pdf")
see.genes(sigs$sig.genes$ColdvsControl,
          show.fit = T,
          dis = design$dis,
          cluster.method = "hclust",
          cluster.data = 1,
          newX11 = F,
          k = 9,
          cexlab = 2)
dev.off()

pheatmap(as.matrix(na.omit(data.abiotic[as.character(sigs$summary$ColdvsControl),])))

# Q8: Do the observed clusters in the heatmap match your expectations?
# Mostly yes, except I would expect Salt_3H 1 to 3 to also cluster together
# with the other Salt and Cold samples. Perhaps the changes in Salt are not 
# yet visible at 3H but take longer to show. Perhaps that is why Salt_3H 
# clusters with controls.

# Q9: Do the results shown in the heatmap match with the previously created 
# venn diagram?
# Based on the venn diagram we can see that ColdvsControl overlaps more with
# SaltvsControl (n = 87) than HeatvsControl (n=26). We therefore expect that
# Cold and Salt samples cluster together which is also shown in the heatmap.

# Q10: Do the same for some of the other experimental groups such as 
# ‘SaltvsControl’
pheatmap(as.matrix(na.omit(data.abiotic[as.character(sigs$summary$SaltvsControl),])))
pheatmap(as.matrix(na.omit(data.abiotic[as.character(sigs$summary$HeatvsControl),])))

data(NBdata)
data(NBdesign)