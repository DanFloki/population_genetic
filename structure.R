# # With R version 3.5 or greater, install Bioconductor packages using BiocManager; see https://bioconductor.org/install
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install()
# 
# BiocManager::install(c("GenomicFeatures", "AnnotationDbi", "LEA"))
library(LEA)
library(fields)

setwd("/run/user/1002/gvfs/smb-share:server=smb.uni-oldenburg.de,share=hpc_home/thesis/GBs/3 elaboration Stacks outputs/structure/")
setwd("50miss/")
source("functions/Conversion.R") #from http://membres-timc.imag.fr/Olivier.Francois/Conversion.R
source("functions/POPSutilities.R")#from http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R
source("functions/POPSutilities.R") # for maps function, perform spatial interpolation of ancestry coefficients
# to get map:
# install.packages("remotes")
# remotes::install_github("prio-data/priogrid")
#library(priogrid)


sample.data<-read.table("sample_data_lonspixm_50miss_unl.csv",sep=" ", header=T)
sample.data<-read.table("sample_data_lonspisch.csv",sep=" ", header=T)
nrow(sample.data[sample.data$species=="longifolia",])
nrow(sample.data[sample.data$species=="spicata",])
sample.data<-read.table("sample_data_spi.csv",sep=" ", header=T)
edit(sample.data)

# plot a map with populations #####
library(plyr) # for ddply
populations<-ddply(sample.data, .(population, latitude, longitude), 
                   summarize,
                   N_spi=sum(species=="spicata"),
                   N_lon=sum(species=="longifolia"),
                   N_sch=sum(species=="schmidtiana"),
                   mixed=if(N_spi>0 & N_lon>0){"lon_spi"}else{if(N_spi>0){"spi"}else{if(N_lon>0){"lon"}else{"sch"}}}
                   )
populations$mixed<-as.factor(populations$mixed)
library(mapplots)
library(maps)
coord=populations[,c(3,2)]
plot(coord,ylab="Latitude",xlab="Longitude", col="white")
map(add = T, col="grey90",
    #wrap=F,
    lforce="l",
    #myborder=c(.001,.004),
    #resolution = 1, #0=full resolution of database,
    fill=T, #xlim = c(-5,155), ylim=c(40,70),
    ylab="Latitude",xlab="Longitude")
cols <- c("#7570B3","#E7298A","#D95F02","#1B9E77")
points(coord,col=cols[populations$mixed],pch=20, cex=0.8)
text(coord,labels = populations$population, pos = 4,cex=.8)
legend("top",fill=cols,c("V. longifolia","V. longifolia and V. spicata","V. schmidtiana","V. spicata"), #c(levels(populations$mixed)),
       bty="n", cex=.8)


library(rworldmap)
mymap<-getMap(resolution="low")
par(mfrow=c(1,1))
plot(mymap,xlim = c(-0,140), ylim=c(55,55) )
points(coord,col=cols[populations$mixed],pch=4, cex=0.8)
text(coord,labels = populations$population, pos = 2,cex=.8)

axis(1)
axis(3)
axis(2)
axis(4)

library(ggmap)
map <- get_map(location = 'Eurasia', zoom = 4)

#check also https://rpubs.com/tlama/sNMFsquirrel
lonspi.filtered.VCF <- read.vcfR("/home/daniele/Dropbox/thesis/GBs/2 elaboration on enlarged set/elaboration files/lonspi_filtered_1.vcf")
colnames(lonspi.filtered.VCF@gt)[-1]
in.file<-"/home/daniele/Dropbox/thesis/GBs/2 elaboration on enlarged set/elaboration files/lonspi_filtered_1.vcf"
vcf2geno(in.file,"lonspi_filtered_1.geno")

#see file tutoRstructure.pdf
load(file="first.snmf.Rdata")
first.snmf<-snmf("lonspi_filtered_1.geno", K=1:5, 
                 ploidy=4, entropy=T,#cross-entropy to choose K later
                 repetitions=10,CPU=4,
                 #alpha=100, 
                 project="new")
save(first.snmf,file="first.snmf.Rdata")

second.snmf<-snmf("lonspi_filtered_1.geno", K=1:7, 
                 ploidy=4, entropy=T,#cross-entropy to choose K later
                 repetitions=100,CPU=4, #100 rep should take about 3,1 h
                 #alpha=100, 
                 project="new")
save(second.snmf,file="second.snmf.Rdata")


show(first.snmf)
summary(first.snmf)
# plot cross-entropy criterion for all runs in the snmf project
plot(first.snmf, col="blue4",pch=19, cex=1.2) #here we choose the number of anc. pop. that reaches plateau e.g. K=2
# select the best run for K=x
best<-which.min(cross.entropy(first.snmf, K=3))
cols <- c("#D95F02","#7570B3","#1B9E77","#E7298A")
bp <- barchart(first.snmf, K=3, run=best, border=NA, space=0, col=cols, #xlab="sample ID", 
               ylab="Ancestry proportions", main="Ancestry matrix") #barplot like structure
axis(1, at=1:length(bp$order), labels=sample.data[bp$order,]$new_id, las=2, cex.axis=.65) # put labels on x axis

#qmatrix<-Q(first.snmf, K=3) #admixture coefficients from a snmf run
#barplot(t(qmatrix), col = cols, border = NA, space = 0,
#        xlab = "Individuals", ylab = "Admixture coefficients")

# map barcharts ####

coord = read.table("coord.txt")
plot(coord,ylab="Latitude",xlab="Longitude",xlim=c(-3,150),ylim=c(42,68))
#with(sample.data, text(coord, labels=sample.data$new_id,pos=4))
#with(sample.data,points(coord, col="blue"))
map1=map(add = T, col="grey90",fill=T)
qmatrix = Q(first.snmf, K=3, run=best)
save(qmatrix,file="qmatrix.Rdata")
load("qmatrix.Rdata")
for(i in 1:nrow(qmatrix)){
  add.pie(z=qmatrix[i,], x=coord[i,1], y=coord[i,2],
          labels=sample.data[i,]$new_id,
          col=cols,radius=0.2) 
}
plot(coord,ylab="Latitude",xlab="Longitude",xlim=c(-3,50),ylim=c(44,54))
map(add = T, col="grey90",fill=T)
qmatrix = Q(first.snmf, K=3, run=best)
for(i in 1:nrow(qmatrix)){
  add.pie(z=qmatrix[i,], x=coord[i,1], y=coord[i,2],
          labels=sample.data[i,]$new_id,
          col=cols,radius=0.1)
}
#average for populations ####
n=0
for(i in 1:nrow(coord)){
  if(!(coord[i,1]==coord[i+1,1])){
    n=n+1
  }
}
n

#map geographycal colors for populations ####

# see https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html and email from Olivier Francois
library(tess3r)
q_matrix<-tess3r::as.qmatrix(qmatrix) 
# load(file="genotypes.Rdata") # extract genotypes with tab {adegenet} function: genotypes<-tab(gi.lonspi.filt, NA.method = "mean")
# tess3.obj <- tess3(X = genotypes, #genotype matrix
#                    coord = coord, K = 1:8, 
#                     method = "projected.ls", ploidy = 4, openMP.core.num = 4)



# bragraph STRUCTURE-like ####

# loading qmatrix obtained from STRUCTURE
q_matrix_K_n<-read.table("qMatr_lonspixm_50miss_K2_clumped.txt",sep="\t")
sample.data<-read.table("sample_data_lonspixm_50miss_unl.csv",sep=" ", header=T)
nrow(q_matrix_K_n)
nrow(sample.data)
row.names(q_matrix_K_n)<-sample.data$new_id


cols <- c("#D95F02","#A6761D","#666666","#7570B3","#66A61E","#E6AB02","#1B9E77","#E7298A")
cols <- c("#1B9E77","#E7298A")
x<-barplot(t(q_matrix_K_n), col = cols, border = F, space = 0,
                xlab = "Individuals", ylab = "Admixture coefficients", names.arg = row.names(q_matrix_K_n), axisnames = F #cex.names = NULL#names.arg=sample.data$new_id, angle=90
        )
text(cex=0.5, x=x+.9, y=-0.01, row.names(q_matrix_K_n), xpd=TRUE, srt=90, pos=2)
# order by probability K=2
q_matrix_ord<-q_matrix_K_n[order(q_matrix_K_n$V1),]
# order by probability K=3
q_matrix_ord<-q_matrix_ord[order(-q_matrix_ord$V4),]

x<-barplot(t(q_matrix_ord), col = cols, border = F, space = 0, 
           xlab = "", ylab = "Admixture coefficients", main="clumped K=2 STRUCTURE results 20 rep.\nunlinked SNPs - 25.73% missing data",
           names.arg = row.names(q_matrix_ord), axisnames = F #cex.names = NULL#names.arg=sample.data$new_id, angle=90
)
text(cex=0.8, x=x+.9, y=-0.01, row.names(q_matrix_ord), xpd=TRUE, srt=90, pos=2)
#select pure individuals ####
q_matr_pure_spi<-q_matrix_ord[q_matrix_ord[,1]>0.99,]
q_matr_pure_lon<-q_matrix_ord[q_matrix_ord[,2]>0.99,]
pure_spi_ids<-row.names(q_matr_pure_spi)
pure_lon_ids<-row.names(q_matr_pure_lon)


# map ####


coord = sample.data[,c(14,13)]
plot(coord,
     ylab="Latitude", xlab="Longitude", 
     xlim=c(min(coord$longitude)-2,max(coord$longitude)+2),
     ylim=c(min(coord$latitude)-2,max(coord$latitude)+2)
     #xlim=c(-3,150),ylim=c(42,68)
     )


q_matrix_K_n<-tess3r::as.qmatrix(q_matrix_K_n)
plot(q_matrix_K_n, coord, method = "map.max", 
     #interpol = FieldsKrigModel(10), 
     main = "Genetic clusters",
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(600,300), cex = .4
     
     #col.palette = cols
     )
text(coord,labels = row.names(q_matrix_K_n), pos = 3,cex=.8)

?plot

# map with average values for populations ####
q_matrix_K_n<-read.table("qMatr_spi_50miss_K5.txt",sep="\t")
sample.data<-read.table("spi_50miss_unl_noProb_sampleData.csv",sep=" ", header=T)
nrow(q_matrix_K_n)
nrow(sample.data)
colnames(q_matrix_K_n)<-c("G1","G2","G3","G4","G5")
samples_info<-cbind(sample.data,q_matrix_K_n)

library(plyr) # for ddply
pop_K_vals<-ddply(samples_info, .(population, latitude, longitude), 
                   summarize,
                   N_spi=sum(species=="spicata"),
                   N_lon=sum(species=="longifolia"),
                   N_sch=sum(species=="schmidtiana"),
                   mixed=if(N_spi>0 & N_lon>0){"lon_spi"}else{if(N_spi>0){"spi"}else{if(N_lon>0){"lon"}else{"sch"}}},
                   avg_G1=mean(G1),
                   avg_G2=mean(G2),
                   avg_G3=mean(G3),
                   avg_G4=mean(G4),
                   avg_G5=mean(G5)
)


pop_qMatrix=pop_K_vals[,c("avg_G1","avg_G2","avg_G3","avg_G4","avg_G5")]
row.names(pop_qMatrix) <- pop_K_vals$population
q_matrix_K_n<-tess3r::as.qmatrix(pop_qMatrix)
coord<-pop_K_vals[,c("longitude","latitude")]
plot(q_matrix_K_n, coord, method = "map.max", 
     #interpol = FieldsKrigModel(10), 
     main = expression("Genetic clusters identified in "*italic("V. spicata")*" populations"),
     xlab = "Longitude", ylab = "Latitude", 
     resolution = c(600,300), cex = .4,
     xlim=c(-4,47)
     #ylim=c(43,60)
     #col.palette = cols
)
text(coord,labels = row.names(q_matrix_K_n), pos = 4,cex=.8)

# to adjust size and proportions of maps work on xlim and ylim:
plot(coord,
     ylab="Latitude", xlab="Longitude", 
     # xlim=c(min(coord$longitude)-2,max(coord$longitude)+2),
     # ylim=c(min(coord$latitude)-2,max(coord$latitude)+2)
     xlim=c(-3,45)#,ylim=c(42,68)
)

# for nice map plots see https://www.r-spatial.org/r/2018/10/25/ggplot2-sf.html
install.packages(c("cowplot", "googleway", "ggplot2", "ggrepel", 
                   "ggspatial", "libwgeom", "sf"))
library(ggplot2)
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")



ggplot(data = world)

# save as map_plot_STRUCTURE_K_2

#asc.raster="http://membres-timc.imag.fr/Olivier.Francois/RasterMaps/Europe.asc"
map.raster="lowResEurope.asc"
map.raster="etopo1.asc"
grid=createGridFromAsciiRaster(map.raster)

# next command generate a memory allocation memory (probably because raster image is too big),
# try run it on the cluster (after loading qmatrix, coord, grid and libraries):
# png("map_plot.png",)
# maps(matrix = qmatrix, coord,grid, constraints, method = "max", main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude")
# dev.off()

# max option: only the cluster with the maximal local contribution to ancestry is represented at each geographic point of the map
maps(matrix = qmatrix, coord,grid, constraints, method = "max", main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude")

# see http://membres-timc.imag.fr/Olivier.Francois/TESS_Plot.html

show.key = function(cluster=1,colorGradientsList=lColorGradients){
  ncolors=length(colorGradientsList[[cluster]])
  barplot(matrix(rep(1/10,10)),col=colorGradientsList[[cluster]][(ncolors-9):ncolors],main=paste("Cluster",cluster))}

layout(matrix(c(rep(1,8),2,3,4,5), 4, 4, byrow = FALSE), widths=c(4,1), respect = T)
par(ps = 22)
maps(matrix = qmatrix, coord, grid, constraints, method = "max", 
     main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude")
par(ps = 16)
for(k in 1:4){show.key(k)}

# # lea.croc <- load.snmfProject("/Users/lc/Desktop/antalyensis_gbs_project/R/ip_crocus.snmfProject")
# qmatrix = Q(lea.croc, K=3, run=best)
# maps(matrix = qmatrix, coord, grid, constraints, method = "max",
#      main = "Ancestry coefficients", xlab = "Longitude", ylab = "Latitude", cex = .5)
# map(add = T, interior = F)