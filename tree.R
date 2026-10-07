setwd("/run/user/1002/gvfs/smb-share:server=smb.uni-oldenburg.de,share=hpc_home/thesis/GBs/3 elaboration Stacks outputs/elab_files_25_unlinked/")


library(ape)
library(phangorn)


load(file="gl.lonspi.filt.Rdata")
sample.data <- read.table("sample_data_filt.csv", sep=" ", header=T,row.names = NULL)

lonspi.VCF <- read.vcfR("lonspi_miss_25_unlinked_ren_filt.vcf")
lonspi.DNAbin<-vcfR2DNAbin(lonspi.VCF)
lonspi.phyDat<-phyDat(lonspi.DNAbin, type="DNA")
lonspi.phyDat<-phyDat(lonspi.VCF)

lonspi.phyDat <- genlight2phyDat(gl.lonspi.filt, ambiguity = NA)
# 
# lonspi.phyDat <- read.phyDat(file, format = "phylip", type = "DNA", ...)


dm<-dist.ml(lonspi.phyDat)
treeUPGMA <-upgma(dm)
treeNJ<-NJ(dm)

plot(treeUPGMA, main="UPGMA")
plot(treeNJ, "unrooted", main="NJ")

# parsimony ####
parsimony(treeUPGMA, lonspi.phyDat)
parsimony(treeNJ, lonspi.phyDat)

# rearrangements to find trees with a lower  parsimony score (implements nearest-neighbor interchanges - NNI, and subtree pruning and regrafting - SPR)
treePars <- optim.parsimony(treeUPGMA, lonspi.phyDat)
# parsimony ratchet (K. Nixon. The parsimony ratchet, a new method for rapid rarsimony analysis. Cladistics, 15:407–414, 1999)
treeRatchet <- pratchet(lonspi.phyDat, trace = 0)
parsimony(c(treePars, treeRatchet), lonspi.phyDat)
plot(treeRatchet, main = "Ratchet")

# Maximum Likelihood ####
# introduced by (oseph Felsenstein. Evolutionary trees from dna sequences: a maxumum likelihood approach. Journal of Molecular Evolution, 17:368–376, 1981.)

# compute likelihood fir a tree:
fit = pml(treeNJ, data=lonspi.phyDat)
fit
# the object fit estimated the likelihood for the tree, but brach length are not optimized for Jukes-Cantor model, so:
fitJC <- optim.pml(fit, TRUE)# optimize branch length for Jukes-Cantor model
logLik(fitJC)
# update.pml allows to change parameters
fitGTR <- update(fit, 
                 k=4, # rate categories
                 inv=0.2)
# fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
#                     rearrangement = "NNI", control = pml.control(trace = 0))


# for larg trees NNI rearrangement often get stuck in local maxima. 
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", # makes random NNI permuation to the tree, which than gets optimised to escape local optima
                    control = pml.control(trace = 0))
fitGTR
#compare nested models for JC and GTR+gamma+I
anova(fitJC, fitGTR)
SH.test(fitGTR, fitJC)

#alternatively, function modelTest
mt = modelTest(lonspi.phyDat, model="all", 
               multicore = T, mc.cores = 4)

bs = bootstrap.pml(fitJC, bs=100, optNni=TRUE,
                     control = pml.control(trace = 0))
