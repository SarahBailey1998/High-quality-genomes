## Analysing rata Moehau SNPs, adapted from Nat Forsdick Kuaka script by Jessie Prebble Oct 2024

## SNPRelate for PCA

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
install.packages("tinytex")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")
BiocManager::install("gdsfmt")
BiocManager::install("MASS")
install.packages("pals")
install.packages("xfun")
install.packages("gridExtra")

library(gdsfmt)
library(SNPRelate)
library(MASS)
library(gridExtra)
library(scales) # allows use of alpha to change plot opacity

sessionInfo()

setwd("T:/Lincoln/Projects A-E/DNA/Metrosideros/filtering")
getwd()

citation("gdsfmt")
citation("SNPRelate")
citation("MASS")

Venables, W. N. & Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New
York. ISBN 0-387-95457-0

##Colours##

mwlrcols <- c("#64a70b", "#898a8d", "#009cbd","#ebb700")

#import_vcfs
vcf.strict <- "rata-moehau-only_VariantCalls_20x_coverage_0site_missing_maf0.0.bcf.recode_0.8LD_VariantCalls.vcf"

showfile.gds(closeall=TRUE)
snpgdsVCF2GDS(vcf.strict, "vcf.strict.gds", method="biallelic.only")
snpgdsSummary("vcf.strict.gds")

#sanity_check
# Open the GDS file
genofile <- snpgdsOpen("vcf.strict.gds")

head(genofile)

#get population information
#hmm. looks like Nat is pulling it from a file and then selecting which column
#pop_code <- scan("T:/Lincoln/Projects A-E/DNA/Metrosideros/filtering/popmap.txt",
#                 what=character())
#table(pop_code)

# Display the first six values
#head(pop_code)

#Ah, I think I need a file with JUST the popcode to make this work, other wise could do it another way, but try that first
#pop_code <- scan("T:/Lincoln/Projects A-E/DNA/Metrosideros/filtering/pop_codes_only.txt",
#                 what=character())

#to redo with the auckland botanic garden sample re coded as Unuwhau
pop_code <- scan("T:/Lincoln/Projects A-E/DNA/Metrosideros/filtering/pop_codes_only2.txt",
                 what=character())

# Display the first six values
head(pop_code)

# Get sample id
sample.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
length(sample.id)

population <- as.factor(pop_code)
pca <- snpgdsPCA(genofile, num.thread=4,autosome.only=F)

pca$sample.id
sample.id.filename <- sub(".*/","",pca$sample.id)

# if we can assume the order of sample IDs is as the same as population codes (always double check this) 
cbind(sample.id,sample.id.filename, pop_code)

# variance proportion (%)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

# Make the data frame
tab <- data.frame(sample.id = pca$sample.id,
                  pop = as.factor(pop_code)[match(pca$sample.id, sample.id)],
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab$sample.id)
tail(tab$sample.id)
head(tab$pop)
tail(tab$pop)

# Draw

#?plot
#plot(tab$EV2, tab$EV1, col=as.integer(tab$pop), 
#     xlab="eigenvector 2", ylab="eigenvector 1")

plot(tab$EV2, tab$EV1, col=alpha(mwlrcols[tab$pop],0.6), pch=19, 
     xlab="PC 2", ylab="PC 1")
legend("topleft", legend=levels(tab$pop), 
       pch=19, col=mwlrcols[1:nlevels(tab$pop)])

tab


plot(tab$EV2, tab$EV1, col=alpha(mwlrcols[tab$pop],0.6), pch=19, 
     xlab="PC 2", ylab="PC 1")
legend("topleft", legend=levels(tab$pop), 
       pch=19, col=mwlrcols[1:nlevels(tab$pop)])
text(x=pca$eigenvect[,2], y=pca$eigenvect[,1], labels=sample.id.filename,pos=1, offset=-1)

lbls <- paste("PC", 1:5, "\n", format(pc.percent[1:5], 
                                      digits=2), "%", sep="")
pairs(pca$eigenvect[,1:5], col=alpha(mwlrcols[tab$pop],0.6), pch=19, labels=lbls)

lbls <- paste("PC", 1:2, "\n", format(pc.percent[1:2], 
                                      digits=2), "%", sep="")
pairs(pca$eigenvect[,1:2], col=alpha(mwlrcols[tab$pop],0.6), pch=19, labels=lbls)


#Parallel coordinates plot for the top principal components:
  
 # ```{r top_PCs, dev=c('pdf'), fig.path='figures-strict/', ppi=500, units="in", fig.height=5, fig.width=7}
#datpop <- factor(pop_code)[match(pca$sample.id, sample.id)]
#parcoord(pca$eigenvect[,1:5], col=alpha(y[datpop],0.6))

#To calculate the SNP correlations between eigenvectors and SNP genotypes:
  
  #```{r get_corr, dev=c('pdf'), fig.path='figures-strict/', ppi=500, units="in", fig.height=5, fig.width=7, eval=FALSE}

# Get chromosome index
##chr <- read.gdsn(index.gdsn(genofile, "snp.chromosome"))
#CORR <- snpgdsPCACorr(pca, genofile, eig.which=1:4)

#savepar <- par(mfrow=c(2,1), mai=c(0.45, 0.55, 0.1, 0.25))

#for (i in 1:2)
#{
#  plot(abs(CORR$snpcorr[i,]), ylim=c(0,1), xlab="",
#       ylab=paste("PC", i),
#       col=1:length(chr), pch="+")
#}
#```
#Can also do Weir-Cockerham Fst

#```{r fst}
# Two populations: HCB and JPT #I'm going to try Ko and Un
flag <- pop_code %in% c("Ko", "Un")

samp.sel <- sample.id[flag]
pop.sel <- pop_code[flag]
v <- snpgdsFst(genofile, sample.id=samp.sel, population=as.factor(pop.sel),
               method="W&C84", autosome.only=F)

# Weir and Cockerham weighted Fst estimate
v$Fst 
# Weir and Cockerham mean Fst estimate
v$MeanFst    
summary(v$FstSNP)
```
#looks like very low levels of differentiation

## Adegenet


```{r adegenet}
install.packages("vcfR")
install.packages("adegenet")
library(vcfR)
library(adegenet)

x <- read.vcfR("rata-moehau-only_VariantCalls_20x_coverage_0site_missing_maf0.01.bcf.recode_0.8LD_VariantCalls.vcf", verbose=F)
y <- vcfR2genind(x, ploidy=2, return.alleles=TRUE)
names(y)

#pop_code <- scan("C:/Users/ForsdickN/OneDrive - MWLR/Documents/WHDP/pop-gen/second-round-analysis/export/corr-map2.txt",
#                 what=character())
y@pop
y@pop <- as.factor(pop_code)
y@pop
```
#Now we have the inputs loaded and are begin testing the analysis. First let's max out our PCs and set 2 expected clusters. How do our individuals group?
#Jessie = should I try 3 clusters?
```{r prelimDAPC, dev=c('pdf'), fig.path='figures-strict/', ppi=500, units="in", fig.height=7, fig.width=10}
grp <- find.clusters(y, max.n.= 30, n.pca=70, n.clust=3)
#head(grp$grp, 10)
grp$size
head(grp$grp)
table(pop(y), grp$grp)
dapc1 <- dapc(y,grp$grp, n.pca = 15, n.da= 3)

table.value(table(pop(y), grp$grp), col.lab=paste("inf", 1:3),
row.lab=paste("ori", 1:3))

scatter(dapc1, scree.da=FALSE, scree.pca=TRUE, posi.pca="topleft", legend=TRUE, solid=0.5, bg="white", col=mwlrcols, clab=0)

summary(dapc1)

# unfortunately there's no argument to dictate the individual labels for assignplot().
assignplot(dapc1)

compoplot(dapc1, posi="bottomright", txt.leg=paste("Cluster", 1:2), lab="", col=mwlrcols, xlab="individuals")
```
#Jessie looking at an adgenet tutorial

install.packages("hierfstat")
library("hierfstat")
#fstat(y) # doesn't appear to exist any more
#basic.stats(y, diploid = TRUE)# looks like it will take too long

rataFST<- pairwise.fst(pop)#looks like this doesn't work anymore either. grr.
pop<-genind2genpop(y)

dist.genpop(pop)
pop_dst<-dist(pop)
ind_dist<-dist(y)
is.euclid(ind_dist)

#k-means clustering

find.clusters(y)
find.clusters(y, max.n.clust = 12)# then chose 9 PCs to retain, then 3 clusters:
$Kstat
K=1      K=2      K=3      K=4      K=5      K=6 
115.9784 116.2601 116.3151 116.1915 115.9857 115.3864 

$stat
K=3 
116.3151 

$grp
EXT049-01_S1  EXT049-02_S2  EXT049-03_S3  EXT049-04_S4  EXT049-05_S5  EXT049-06_S6  EXT049-07_S7  EXT049-08_S8 
2             2             1             1             1             3             3             1 
EXT049-09_S9 EXT049-10_S10 EXT049-11_S11 EXT049-12_S12 
1             1             3             3 
Levels: 1 2 3

$size
[1] 6 2 4

#if i chose 10 PCs to retain then 5 was the best number of clusters:
$Kstat
K=1      K=2      K=3      K=4      K=5      K=6 
116.8123 117.2552 117.5222 117.6741 117.8138 117.6987 

$stat
K=5 
117.8138 

$grp
EXT049-01_S1  EXT049-02_S2  EXT049-03_S3  EXT049-04_S4  EXT049-05_S5  EXT049-06_S6  EXT049-07_S7  EXT049-08_S8 
3             3             1             2             1             5             5             4 
EXT049-09_S9 EXT049-10_S10 EXT049-11_S11 EXT049-12_S12 
2             2             5             5 
Levels: 1 2 3 4 5

$size
[1] 2 3 2 1 4

dapc1<-dapc(y,grp$pop_code)# chose 8 PCs to retain, and 2 discriminant functions (has to be more than 1, one looks best)  
scatter(dapc1)

dapc2<-dapc(y)

scatter(dapc1, scree.da=FALSE, bg="white", pch=20, cell=0, cstar=0, col=mwlrcols, solid=.4,
       cex=3,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:4))
