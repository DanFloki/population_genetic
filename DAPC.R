# DAPC analysis:
#  - determine which alleles are more informative for clustering according to species
#  - identification of genetic clusters in our data 

#steps:
# 1) load data
# 2) enstabilish an appropriate number of PCs to use with cross-validation (it's a crucial argument to function dapc)
# 3) identification of genetic clusters (similar to structure, does not require any assumption on populations)
# 4) describe clusters inferred using dapc, running also dapc providing original groups (here species)
# 5) plot results of dapc in a STRUCTURE-like fashion
# 6) determine which alleles differentiate more between groups (useful to subset genomic data to further analysis)

# check here https://grunwaldlab.github.io/Population_Genetics_in_R/DAPC.html
# see also http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf 

setwd("elab_files_50_no_prob/")

library("adegenet")
library(vcfR)

cols <- c("#D95F02","#7570B3","#1B9E77","#E7298A")

# loading data and info ####
lonspi.filtered.VCF <- read.vcfR("lonspi_miss_25_ren_filt_no_prob.vcf")
lonspi.only.VCF<-read.vcfR("lonspi_only_50miss_noProb.VCF.gz")
load(file="gi.lonspi.only.RData") # object saved before with R command: save(object_name,file="object_name.Rdata")
load(file="gi.lonspi.filt.RData") # object saved before with R command: save(object_name,file="object_name.Rdata")

sample.data<- read.table("sampleData_lonspi_only_50miss_noProb.csv", sep=" ", header=T)
load(file="genotypes.Rdata")
lon.VCF <- read.vcfR("lon.VCF.gz")
spi.VCF <- read.vcfR("spi.VCF.gz")
sample.data.lon<- read.table("sample_data_lon.csv", sep=" ", header=T,row.names = NULL)
sample.data.spi<-read.table("sample_data_spi.csv",sep=" ", header = T,row.names = NULL)


View(sample.data.lon)
# change names to sample to make them more informative (e.g. sp47t_Hun7)
#new.names<-read.table("new_names_filt_1.txt",sep="\t", header=T)
#setdiff(new.names$ID, colnames(lonspi.filtered.VCF@gt)) #check if there is some different ID in the list, should not
#colnames(lonspi.filtered.VCF@gt)[-1]<-as.character(new.names$new_id)

# convert in genind
gi.lonspi.filt <- vcfR2genind(lonspi.filtered.VCF)
save(gi.lonspi.filt, file="gi.lonspi.filt.RData")
save(gi.lonspi.filt, file="gi.lonspi2.filt.RData")

gi.lonspi.only <-vcfR2genind(lonspi.only.VCF)
save(gi.lonspi.only, file="gi.lonspi.only.RData")
gi.lon.filt<-vcfR2genind(lon.VCF)
gi.spi.filt<-vcfR2genind(spi.VCF)
save(gi.lon.filt, file="gi.lon.filt.RData")
save(gi.spi.filt, file="gi.spi.filt.RData")
#got this warning: entirely non-type marker(s) deleted
#there were entirely non-typed locy


# convert in genlight
gl.lonspi.filt <- vcfR2genlight(lonspi.filtered.VCF)

# assign populations (according to species) and attach other info
gi.lonspi.filt$strata<-sample.data
pop(gi.lonspi.filt)<-sample.data$species

gi.lonspi.only$strata<-sample.data
pop(gi.lonspi.only)<-sample.data$population
# gl.lonspi.filt$strata<-sample.data.filtered
# pop(gl.lonspi.filt)<-sample.data.filtered$species
gi.lon.filt$strata<-sample.data.lon
pop(gi.lon.filt)<-sample.data.lon$population

gi.spi.filt$strata<-sample.data.spi
pop(gi.spi.filt)<-sample.data.spi$population

gi.spi<-gi.lonspi.filt[gi.lonspi.filt@strata$species=="spicata"]
gi.lon<-gi.lonspi.filt[gi.lonspi.filt@strata$species=="longifolia"]
gi.lonspi.only<-gi.lonspi.filt[gi.lonspi.filt@strata$species=="longifolia" | gi.lonspi.filt@strata$species=="spicata"]

# save genind object to save time
save(gi.lonspi.filt, file="gi.lonspi.filt.RData")
save(gi.lon, file="gi.lon.Rdata")
save(gi.spi, file="gi.spi.Rdata")
save(gi.lonspi.only, file="gi.lonspi.only.Rdata")

# extract genotypes with tab {adegenet} function
genotypes<-tab(gi.lonspi.filt, NA.method = "mean")
save(genotypes,file="gi.lonspi.filt.genotypes.Rdata")
gen.lon<-tab(gi.lon.filt,NA.method = "mean")
gen.spi<-tab(gi.spi.filt, NA.method = "mean")


# spicata ####
#cross-validation
load(file="gi.lon.filt.RData")
xv.spi<-xvalDapc(gen.spi, pop(gi.spi.filt),n.rep=10,ncpus = 4L) # 35 PC
save(xv.spi, file="xv.spi.Rdata")
pops<-pop(gi.spi.filt)

library(RColorBrewer)
cols <- brewer.pal(n = nPop(gi.spi.filt), #number of different colors in the palette
                   name = "Dark2")
scatter(xv.spi$DAPC,  #col = cols, 
        cex = 0.5, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,cstar=0)
grp <- find.clusters(gi.spi.filt, n.pca=35, 
                     max.n.clust=25) #evaluate up to K=40
#K=23
grp
names(grp)
head(grp$Kstat)
grp$stat
grp$size
table(pop(gi.spi.filt),grp$grp)
#table with actual groups (ori) amd inferred groups (inf)
table.value(table(gi.lonspi.only@strata$species, grp$grp), col.lab=paste("inf", 1:23),
            row.lab=paste("ori", 1:23))

#structure-like plots ####
dapc2 <- dapc(gi.spi.filt, grp$grp,n.pca = 35)

#qmatrix_K_4<-dapc2$posterior
#write.csv(qmatrix_K_4,"qmatrix_K_4.csv", sep="\t")

library(Polychrome)
cols <- palette36.colors(36)
cols<-heat.colors(23,alpha=1)
compoplot(dapc2, posi="top", col=cols,
          txt.leg=levels(grp$grp),
          xlab="individuals")

dapc2.results <- as.data.frame(dapc2$posterior)
dapc2.results$pop <- pop(gi.spi.filt)
dapc2.results$indNames <- rownames(dapc2.results)
library(reshape2)
dapc2.results <- melt(dapc2.results) # this because ggplot2 requires each observation in row and values observed in columns
colnames(dapc2.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
library(ggplot2)
p <- ggplot(dapc2.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free",space="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p


dapc1.results <- as.data.frame(dapc1$posterior)
#dapc1.results$pop <- grp$grp
dapc1.results$pop <-pop(gi.lonspi.filt)
dapc1.results$indNames <- rownames(dapc1.results)
library(reshape2)
dapc1.results <- melt(dapc1.results) # this because ggplot2 requires each observation in row and values observed in columns
colnames(dapc1.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
library(ggplot2)
cols2<-c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
p <- ggplot(dapc1.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols2) 
p <- p + facet_grid(~Original_Pop, scales = "free",space="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

# subset for spicata and longifolia only ####
gi.lonspi.filt<-gi.lonspi.filt[gi.lonspi.filt@strata$species=="longifolia" |
                                 gi.lonspi.filt@strata$species=="spicata"]
gi.lonspi.filt@strata$species
genotypes<-tab(gi.lonspi.filt, NA.method = "mean")

genotypes<-tab(gi.lonspi.only, NA.method = "mean")

# spicata and longifolia ####
# cross val
load(file="gi.lonspi.only.RData")
genotypes.lonspi<-tab(gi.lonspi.only, NA.method = "mean")
save(genotypes.lonspi,file="genotypes.lonspi.RData")
load(file="genotypes.lonspi.RData")
#gi.lonspi.only@strata<-sample.data

# to save plot on a file in cluster:
# options(bitmapType='cairo') #this option necessary to avoid X11 error
# png("xval_lonspi_50miss.png")
# dev.off()

xval <- xvalDapc(genotypes.lonspi, gi.lonspi.only@strata$species, 
                 n.pca = 1:150, training.set = 0.9, # ptop of data to be used as training set (for groups >= 10 ind.)
                 result = "overall", # groupMean = average prediction success per group; overall = overall prediction success
                 center = F, scale = FALSE,
                 n.rep = 1, xval.plot = TRUE,
                 parallel="multicore", ncpus = 2)
save(xval,file="xval_lonspi_only.Rdata")

xval <- xvalDapc(genotypes.lonspi, gi.lonspi.only@strata$species, n.pca = 1:10, training.set = 0.9, result = "overall", center = F, scale = FALSE, n.rep = 10, xval.plot = TRUE, parallel="multicore", ncpus = 2)

# load(file="xval_lonspi_only.Rdata")
xval[-1] # gives us the number of PCs with lowest error (here 8)
scatter(xval$DAPC,  col = cols, 
        cex = 0.5, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,cstar=0)

# find clusters
grp <- find.clusters(gi.lonspi.only, stat="BIC", #BIC: Bayesian Information Criterion. AIC: Aikaike's Information Criterion. WSS: within-groups sum of squares, that is, residual variance.
                     max.n.clust=50, n.iter=100) #evaluate up to K=40
# the function displays a graph of cumulated variance explained by the eigenvalues of PCA.
# here we keep all of them (we tipe 200 but the max will be anyway 156)
# after is displayed the value of BIC (Bayesian Information Criterion) using different clusters,
# we should choose a K in which we see a minimun of BIC (an elbow in the curve): 3 (somethimes is K=2, difference not much high)
# choose K=4 to have more eigenvalues to choose in next steps
grp
names(grp)
head(grp$Kstat)
grp$stat
grp$size
table(gi.lonspi.only@strata$species,grp$grp)
#table with actual groups (ori) amd inferred groups (inf)
table.value(table(gi.lonspi.only@strata$species, grp$grp), #col.lab=paste("inf", 1:2),
            row.lab=paste("ori", 1:6))

gi.lonspi.only@strata$species<-as.character(gi.lonspi.only@strata$species)
gi.lonspi.only@strata$species<-as.factor(gi.lonspi.only@strata$species)
levels(gi.lonspi.only@strata$species)
# describing clusters using DAPC
dapc1 <- dapc(gi.lonspi.only,n.pca=4, grp$grp) #groups provided are the ones identified previously with find.clusters() method
# the more PCs choose to retain (max 156) the more the scatterplot will produce concentrated 
# groups. With 156 PCs all the groups will appear as single dots

#memberhip probabilities (see pag 22 of Vignette DAPC)
dapc1$posterior
summary(dapc1) #slots assign.per.pop indicates the proportion of successful reassignment (based on the discriminant function)
#of individuals to their original clusters. Large values indicate clear-cut clusters, while low values suggest admixed group.
#see also assignplot
assignplot(dapc1,cex.lab=0.6) # heat colors represent membership probabilities (red=1, white=0),
#blue cross indicate original cluster

dapc1.results <- as.data.frame(dapc1$posterior)
dapc1.results$pop <- grp$grp
#dapc1.results$pop < sample.data$species
dapc1.results$indNames <- rownames(dapc1.results)
library(reshape2)
dapc1.results <- melt(dapc1.results) # this because ggplot2 requires each observation in row and values observed in columns
colnames(dapc1.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
library(ggplot2)
cols2 <- c("#D95F02","#A6761D","#666666","#7570B3","#66A61E","#E6AB02","#1B9E77","#E7298A")
p <- ggplot(dapc1.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols2) 
p <- p + facet_grid(~Original_Pop, scales = "free",space="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p
scatter(dapc1, 
        col = cols2, 
        cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomright", scree.pca = TRUE, posi.da = "topright",
        posi.pca = "bottomleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,cstar=0,lwd=2,lty=2)

dapc2 <- dapc(gi.lonspi.only, gi.lonspi.only@strata$species, n.pca = 10)
assignplot(dapc2,cex.lab =0.6) # heat colors represent membership probabilities (red=1, white=0),
#blue cross indicate original cluster
cols2<-c("#F8766D","#7CAE00","#00BFC4","#C77CFF")

qmatrix_K_2<-dapc2$posterior
write.csv(qmatrix_K_4,"qmatrix_K_4.csv", sep="\t")

# library(RColorBrewer)
# cols <- brewer.pal(n = nPop(gi.lonspi.filt), #number of different colors in the palette
#                    name = "Dark2")

compoplot(dapc2, posi="top", col=cols2,
          txt.leg=levels(gi.lonspi.only@strata$species),
          xlab="individuals",show.lab=T)

dapc2.results <- as.data.frame(dapc2$posterior)
dapc2.results$pop <- gi.lonspi.only@strata$species
dapc2.results$indNames <- rownames(dapc2.results)
library(reshape2)
dapc2.results <- melt(dapc2.results) # this because ggplot2 requires each observation in row and values observed in columns
colnames(dapc2.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
library(ggplot2)
p <- ggplot(dapc2.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols2) 
p <- p + facet_grid(~Original_Pop, scales = "free",space="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p


# spicata only ####
load(file="gi.spi.RData")
gi.spi<-gi.lonspi.only[gi.lonspi.only@strata$species=="spicata"]
gi.spi@strata$population<-as.character(gi.spi@strata$population)
gi.spi@strata$population<-as.factor(gi.spi@strata$population)
genotypes.spi<-tab(gi.spi, NA.method = "mean")
# cross val
xval <- xvalDapc(genotypes.spi, pop(gi.spi), 
                 n.pca = 1:40, training.set = 0.9, # ptop of data to be used as training set (for groups >= 10 ind.)
                 result = "overall", # groupMean = average prediction success per group; overall = overall prediction success
                 center = F, scale = FALSE,
                 n.rep = 10, xval.plot = TRUE,
                 parallel="multicore",ncpus = 2)
save(xval,file="xval_spi.Rdata")
xval[-1] # gives us the number of PCs with lowest error (here 37)
scatter(xval$DAPC,  col = cols, 
        cex = 0.5, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,cstar=0)

# find clusters
grp <- find.clusters(gi.spi, stat="BIC", #BIC: Bayesian Information Criterion. AIC: Aikaike's Information Criterion. WSS: within-groups sum of squares, that is, residual variance.
                     max.n.clust=45, n.iter=1000) #evaluate up to K=40
# the function displays a graph of cumulated variance explained by the eigenvalues of PCA.
# here we keep all of them (we tipe 200 but the max will be anyway 156)
# after is displayed the value of BIC (Bayesian Information Criterion) using different clusters,
# we should choose a K in which we see a minimun of BIC (an elbow in the curve): 3 (somethimes is K=2, difference not much high)
# choose K=4 to have more eigenvalues to choose in next steps
grp
names(grp)
head(grp$Kstat)
grp$stat
grp$size
table(gi.spi@strata$country,grp$grp)
#table with actual groups (ori) amd inferred groups (inf)
table.value(table(gi.spi@strata$continent, grp$grp), #col.lab=paste("inf", 1:2),
            row.lab=paste("ori", 1:6))

# describing clusters using DAPC
dapc1 <- dapc(gi.spi,n.pca=100, grp$grp) #groups provided are the ones identified previously with find.clusters() method
# the more PCs choose to retain (max 156) the more the scatterplot will produce concentrated 
# groups. With 156 PCs all the groups will appear as single dots

#memberhip probabilities (see pag 22 of Vignette DAPC)
dapc1$posterior
summary(dapc1) #slots assign.per.pop indicates the proportion of successful reassignment (based on the discriminant function)
#of individuals to thei original clusters. Large values indicate clear-cut clusters, while low values suggest admixed group.
#see also assignplot
assignplot(dapc1,cex.lab=0.6) # heat colors represent membership probabilities (red=1, white=0),
#blue cross indicate original cluster
dapc1.results <- as.data.frame(dapc1$posterior)
dapc1.results$pop <- grp$grp
#dapc1.results$pop < sample.data$species
dapc1.results$indNames <- rownames(dapc1.results)
library(reshape2)
dapc1.results <- melt(dapc1.results) # this because ggplot2 requires each observation in row and values observed in columns
colnames(dapc1.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
library(ggplot2)
cols2<-c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
p <- ggplot(dapc1.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols2) 
p <- p + facet_grid(~Original_Pop, scales = "free",space="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

# map
q_matrix_K_n<-as.qmatrix(q_matrix_K_n)
plot(q_matrix_K_n, coord, method = "map.max", interpol = FieldsKrigModel(10),  
     main = "Ancestry coefficients",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(300,300), cex = .4
     #col.palette = my.palette)
)


dapc2 <- dapc(gi.spi, gi.spi@strata$area, n.pca = 20)#60 PCs from cross-valid n.rep=100
assignplot(dapc2,cex.lab =0.6) # heat colors represent membership probabilities (red=1, white=0),
#blue cross indicate original cluster
cols2<-c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
scatter(dapc1, 
        #col = cols2, 
        cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE, posi.da = "topright",
        posi.pca = "bottomright", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,cstar=0,lwd=2,lty=2)


qmatrix_K_2<-dapc2$posterior
write.csv(qmatrix_K_4,"qmatrix_K_4.csv", sep="\t")

# library(RColorBrewer)
# cols <- brewer.pal(n = nPop(gi.lonspi.filt), #number of different colors in the palette
#                    name = "Dark2")

compoplot(dapc2, posi="top", col=cols,
          txt.leg=levels(sample.data$species),
          xlab="individuals",show.lab=T)

dapc2 <- dapc(gi.spi, gi.spi@strata$continent, n.pca = 10)#60 PCs from cross-valid n.rep=100
dapc2.results <- as.data.frame(dapc2$posterior)
dapc2.results$pop <- gi.spi@strata$continent
dapc2.results$indNames <- rownames(dapc2.results)
library(reshape2)
dapc2.results <- melt(dapc2.results) # this because ggplot2 requires each observation in row and values observed in columns
colnames(dapc2.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
library(ggplot2)
p <- ggplot(dapc2.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols2) 
p <- p + facet_grid(~Original_Pop, scales = "free",space="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p




# ALL DATA ############
# cross-validation ####
# enstablish appropriate number of principal components (PCs) 
# we can run a cross-validation test (leave out a percentage of data, run DAPC, and check 
# if data that was left out is correctly placed in the population)


# we can see we have a peak around 60 PC: we can narrow the search specifying the number of PC to
# try with n.pca, doing several replicates each
# if on the next get error, probably is because you need to assign again pop(gi.lonspi.filt)<-sample.data.filtered$species
# pramx <- xvalDapc(genotypes, pop(gi.lonspi.filt), # sample.data.filtered$species
#                   n.pca = 55:65, n.rep = 80, #with those values (n.pca = 55:65, n.rep = 20) RAM around limit (16Gb)
#                   parallel = "multicore", ncpus = 4L)

xval <- xvalDapc(genotypes, gi.lonspi.only@strata$species, 
                 n.pca = 1:200, training.set = 0.9,
                 result = "groupMean", center = F, scale = FALSE,
                 n.rep = 30, xval.plot = TRUE,
                 ncpus = 4L)

xval[-1] # gives us the number of PCs with lowest error (here 3)

save(xval,file="xval.Rdata")

scatter(xval$DAPC,  col = cols, 
        cex = 0.5, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,cstar=0)

# assessing which alleles differentiate more between groups ####
dapc10PC<-dapc(gi.lonspi.filt, var.contrib = T, #retain variable contributing to the analysis in output
              scale=F, n.pca = 13, # 60 PCs from cross-validation
              n.da = nPop(gi.lonspi.filt)-1)

scatter(dapc10PC, 
        col = cols, 
        cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "topleft", scree.pca = F,
        posi.pca = "top", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,cstar=0,lwd=2,lty=2)

set.seed(3)
contrib <- loadingplot(dapc10PC$var.contr,axis = 1, thres = 0.0004, lab.jitter = 1)
#contrib
contrib$var.names



dapc60PC$var.contr
orderedContr <- dapc60PC$var.contr[order(-dapc60PC$var.contr[,1]),] #descending order based on LD1
selectedAlleles <- orderedContr[1:300,] #select first 300 alleles (most importants for DAPC separation)
selectedAlleles
#extract names of loci
selectedLoci<-unique(substr(row.names(selectedAlleles),1,regexpr(".",row.names(selectedAlleles),fixed=T)-1)) 
write.table(as.matrix(selectedLoci),file="selecLoci.txt",sep="\t", quote=F,row.names = F,col.names = F)

dim(lonspi.filtered.VCF@fix)[1] #number of SNPs

#subset VCF file for certain SNPs
lonspi.subsetted.VCF <- lonspi.filtered.VCF[colnames(lonspi.filtered.VCF@fix) = selectedLoci,]


lonspi.subsetted.VCF@fix <- lonspi.filtered.VCF@fix[selectedLoci,]
vcf@gt <- vcf@gt[keep_snps, ]



loci<-seploc(gi.lonspi.filt) #seploc {adegenet} creates a list of individual loci
snp1<-tab(loci[["Contig59_13929338_42419"]]) #tab {adegenet} returns a matrix of genotypes
# following command find the average allele frequencies per population
freqsnp1<-apply(snp1, 2, function(e) tapply(e, pop(gi.lonspi.filt), mean, na.rm=T))


# to subset biggest data:
# data <- data[with(data,order(-Score)),]
# data <- data[1:10,]


# identification of structural SNPs using snzip #### 
# see https://rdrr.io/cran/adegenet/man/snpzip.html

dapc1 <- dapc(genotypes, pop(gi.lonspi.filt), n.da = nPop(gi.lonspi.filt)-1, n.pca = 80, var.contrib = T, scale=F)
features <- snpzip(genotypes, #SNPs matrix used as input for DAPC
                   dapc1, #factor indicating group membership of individual or dapc obj.
                   xval.plot=T,
                   loading.plot = TRUE, method = "centroid")

# generate separate loading plots for each DA
DA <- c(1:dapc1$n.da)
par(ask=F)

for(i in DA){
  title <- paste("Loading Plot for DA", i, sep=" ")
  maximus <- features$FS[[i]][[2]]
  cutoff <- abs(dapc1$var.contr[maximus,i][(which.min(dapc1$var.contr[maximus,i]))])-0.000001
  loadingplot(dapc1$var.contr[, i], threshold=cutoff, main=title)
}

scatter(dapc1, #col = other(gi.lonspi.filt)$comparePal, 
        cex = 0.5, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = TRUE,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
dapc1$var.contr
features
#names of selected alleles on DA1
features$FS[[1]]$`Names of selected alleles`
features$FS[[1]]$`Number of selected vs. unselected alleles`



########################### identification of clusters ####
# see http://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf

?find.clusters

par(mfrow=c(1,1))

indNames(gi.lonspi.filt)<-sample.data.filtered$new_id

grp <- find.clusters(gi.lonspi.only, stat="BIC", #BIC: Bayesian Information Criterion. AIC: Aikaike's Information Criterion. WSS: within-groups sum of squares, that is, residual variance.
                     max.n.clust=100, n.iter=1000) #evaluate up to K=40
# the function displays a graph of cumulated variance explained by the eigenvalues of PCA.
# here we keep all of them (we tipe 200 but the max will be anyway 156)
# after is displayed the value of BIC (Bayesian Information Criterion) using different clusters,
# we should choose a K in which we see a minimun of BIC (an elbow in the curve): 3 (somethimes is K=2, difference not much high)
# choose K=4 to have more eigenvalues to choose in next steps
grp
names(grp)
head(grp$Kstat)
grp$stat
grp$size
table(pop(gi.lonspi.filt),grp$grp)
#table with actual groups (ori) amd inferred groups (inf)
table.value(table(sample.data$species, grp$grp), col.lab=paste("inf", 1:6),
            row.lab=paste("ori", 1:6))
# it seems that to describe the data 2 clusters are usefull (K=2)
nrow(gi.lonspi.filt)

# describing clusters using DAPC #####
dapc1 <- dapc(gi.lonspi.only,n.pca=3, grp$grp) #groups provided are the ones identified previously with find.clusters() method
# the more PCs choose to retain (max 156) the more the scatterplot will produce concentrated 
# groups. With 156 PCs all the groups will appear as single dots

#memberhip probabilities (see pag 22 of Vignette DAPC)
dapc1$posterior
summary(dapc1) #slots assign.per.pop indicates the proportion of successful reassignment (based on the discriminant function)
               #of individuals to thei original clusters. Large values indicate clear-cut clusters, while low values suggest admixed group.
                #see also assignplot
assignplot(dapc1,cex.lab=0.6) # heat colors represent membership probabilities (red=1, white=0),
                  #blue cross indicate original cluster
dapc2 <- dapc(gi.lonspi.only, sample.data$species, n.pca = 10, n.da=2)#60 PCs from cross-valid n.rep=100
assignplot(dapc2,cex.lab =0.6) # heat colors represent membership probabilities (red=1, white=0),
#blue cross indicate original cluster
cols2<-c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
scatter(dapc1, 
        #col = cols2, 
        cex = 2, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomright", scree.pca = TRUE, posi.da = "topleft",
        posi.pca = "bottomleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1,cstar=0,lwd=2,lty=2)

#structure-like plots ####
qmatrix_K_2<-dapc2$posterior
write.csv(qmatrix_K_4,"qmatrix_K_4.csv", sep="\t")

library(RColorBrewer)
cols <- brewer.pal(n = nPop(gi.lonspi.filt), #number of different colors in the palette
                   name = "Dark2")

compoplot(dapc2, posi="top", col=cols,
          txt.leg=levels(sample.data$species),
          xlab="individuals",show.lab=T)

dapc2.results <- as.data.frame(dapc2$posterior)
dapc2.results$pop <- sample.data$species
dapc2.results$indNames <- rownames(dapc2.results)
library(reshape2)
dapc2.results <- melt(dapc2.results) # this because ggplot2 requires each observation in row and values observed in columns
colnames(dapc2.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
library(ggplot2)
p <- ggplot(dapc2.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols) 
p <- p + facet_grid(~Original_Pop, scales = "free",space="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p


dapc1.results <- as.data.frame(dapc1$posterior)
dapc1.results$pop <- grp$grp
#dapc1.results$pop < sample.data$species
dapc1.results$indNames <- rownames(dapc1.results)
library(reshape2)
dapc1.results <- melt(dapc1.results) # this because ggplot2 requires each observation in row and values observed in columns
colnames(dapc1.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")
library(ggplot2)
cols2<-c("#F8766D","#7CAE00","#00BFC4","#C77CFF")
p <- ggplot(dapc1.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity') 
p <- p + scale_fill_manual(values = cols2) 
p <- p + facet_grid(~Original_Pop, scales = "free",space="free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p



#look at subset of ind.
compoplot(dapc2, subset=75:150, posi="top",
          txt.leg=levels(pop(gi.lonspi.filt)), lab="",
          xlab="individuals", col=cols)

#which are the most admixed individuals (ind. with no more than 90% prob of membership in a cluster)?
admixedInd <- which(apply(dapc2$posterior,1, function(e) all(e<0.999)))
admixedInd
compoplot(dapc2, subset=admixedInd, posi="bottomright",
          txt.leg=paste("Cluster", 1:2),
          ncol=2, col=funky(4),show.lab=T)

# it seems we have overfitting (too many PCs retained?). 
# a-score ####
# alpha-score: difference between the proportion of successful reassignment of
# the analysis (observed discrimination) and values obtained using random groups (random discrimination)
# BETTER USE CROSS-VALIDATION, see http://lists.r-forge.r-project.org/pipermail/adegenet-forum/2014-September/000963.html
dapc3<-dapc(gi.lonspi.filt, pop(gi.lonspi.filt),n.da = 100, n.pca = 80)
temp <- optim.a.score(dapc3) #n.pca to actually evaluate all solutions, otherwise method evaluate few numbers of PCs in the range provided
#optimal number of PCs: 19
dapc4<-dapc(gi.lonspi.filt, pop(gi.lonspi.filt),n.da = 100, n.pca = 20)
compoplot(dapc4, posi="bottomright",
          txt.leg=paste("Cluster", 1:nPop(gi.lonspi.filt)), lab="",
          ncol=2, xlab="individuals", col=funky(nPop(gi.lonspi.filt)))
admixedInd <- which(apply(dapc4$posterior,1, function(e) all(e<0.90)))
admixedInd
compoplot(dapc4, subset=admixedInd, posi="bottomright",
          txt.leg=paste("Cluster", 1:4),
          ncol=2, col=funky(4),show.lab=T)

dapc1
dapc1$var.contr #contains contributions of alleles to each discriminant function
scatter(dapc1)

scatter(dapc1, #col = other(gi.lonspi.filt)$comparePal, 
        cex = 0.5, legend = TRUE,
        clabel = FALSE,  scree.pca = TRUE,
        posi.pca = "topleft", posi.da = "bottomleft",posi.leg = "topright",
        #cleg = 0.75,
        bg="white",
        xax = 1, yax = 2, inset.solid = 1)


# plotting DAPC ####
?scatter.dapc
scatter(dapc1, scree.da=FALSE, bg="white", 
        pch=20, 
        cell=1, 
        cstar=0, 
        solid=.4,
        cex=0.6,clab=0, leg=TRUE, txt.leg=paste("Cluster",1:4)
        ,inset.solid=1
        )

myInset <- function(){
  temp <- dapc1$pca.eig
  temp <- 100* cumsum(temp)/sum(temp)
  plot(temp, col=rep(c("black","lightgrey"),
                     c(dapc1$n.pca,1000)), ylim=c(0,100),
       xlab="PCA axis", ylab="Cumulated variance (%)",
       cex=1, pch=20, type="h", lwd=2)
}
add.scatter(myInset(), posi="top",
            inset=c(-0.03,-0.01), ratio=.28,
            bg=transp("white"))

# simulate a dataset and select SNPs ####
# see http://lists.r-forge.r-project.org/pipermail/adegenet-forum/2014-October/000970.html

# make a simulated dataset with 5 "groups"
simpop <- glSim(200, 1000, 40, k=5, sort.pop=TRUE)
snps <- as.matrix(simpop)
phen <- simpop@other$ancestral.pops
dapc1 <- dapc(snps, phen, n.pca=50, n.da=4)
scatter(dapc1)
foo <- snpzip(snps, phen, xval.plot=TRUE, loading.plot=TRUE,
              method="centroid")
# isolate the DAPC component of the snpzip results, calling it "dapc1"
dapc1 <- foo$DAPC

# specify that you want to run the following lines for all DA (ie. from
# DA=1 to DA=(k-1), where K is the number of groups in your dataset)
DA <- c(1:dapc1$n.da)
par(ask=TRUE)

# generate separate loading plots for each DA
for(i in DA){
  title <- paste("Loading Plot for DA", i, sep=" ")
  maximus <- features$FS[[i]][[2]]
  cutoff <- abs(dapc1$var.contr[maximus,i][(which.min(dapc1$var.contr[maximus,i]))])-0.000001
  loadingplot(dapc1$var.contr[, i], threshold=cutoff, main=title)
}