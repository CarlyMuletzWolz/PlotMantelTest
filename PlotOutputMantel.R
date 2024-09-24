library(phyloseq)
data(GlobalPatterns)


## 26 samples
GlobalPatterns

sample_data(GlobalPatterns)$X.SampleID

BrayForPlotBac <- distance(GlobalPatterns, method ="bray")
BrayForMantelBac <- as.matrix(distance(GlobalPatterns, method ="jaccard", binary = T))

## Create another matrix, say each row is a different viral species and are counts 
## in each of the 26 samples (columns)
M1<-matrix(sample(0:60,156,replace=TRUE),ncol=26, nrow = 6)
M1

colnames(M1) <- sample_data(GlobalPatterns)$X.SampleID
 

M1_T = otu_table(M1, taxa_are_rows = TRUE)
physeq_vir <- phyloseq(M1_T)

physeq_vir

BrayForPlotVir <- distance(physeq_vir, method ="bray")
BrayForMantelVir <- as.matrix(distance(physeq_vir, method ="bray"))

library(vegan)
mantel(BrayForMantelBac, BrayForMantelVir, method = "pearson", permutations = 10000)

## Not significant, but let's plot!


dataP = data.frame(cbind(BrayForPlotBac, BrayForPlotVir))


library(ggplot2)
p <- ggplot(dataP, aes(BrayForPlotVir, BrayForPlotBac)) +
  geom_point()  +theme_bw()+ theme_classic() + 
  stat_smooth(method = lm, se = T, color = "black") + 
  geom_point() + theme(text = element_text(size = 16)) + 
  xlab("\nDistance between pairs in viral ASVs (Bray)") +
  ylab("Distance between pairs in bacterial ASVs (Bray)\n")


p


library(ape)
## Make plots with pc axis 1 of each
## NOTE this removes a lot, a lot of variation in the data as 
## likely PCoA only explains small amount of variation. Can plot to figure that out
## or sure some other code to figure out how much variation explained in 1st axis

p_pcoa <- pcoa(BrayForMantelBac)

str(p_pcoa$vectors)

pc1 <- p_pcoa$vectors[,1]
pc1

p_pcoa2 <- pcoa(BrayForMantelVir)

pc2 <- p_pcoa2$vectors[,1]
pc2

dfPCoA  <- cbind(pc1, pc2)

pP <-  ggplot(data = dfPCoA, aes(x = pc2, y = pc1)) +theme_bw()+ theme_classic() + 
  stat_smooth(method = lm, se = T, color = "black") + geom_point() + 
  theme(text = element_text(size = 16)) + 
  xlab("\nPCoA 1 viral ASVs") +ylab("PCoA 1 bacterial ASVs\n")


pP



