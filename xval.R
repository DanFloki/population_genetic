library("adegenet")

setwd("/user/posu1084/thesis/GBs/3 elaboration Stacks outputs/elab_files_50_no_prob")

load(file="gi.lonspi.only.Rdata")
genotypes.lonspi<-tab(gi.lonspi.only, NA.method = "mean")

options(bitmapType='cairo')
png("xval_lonspi_50miss.png")
xval <- xvalDapc(genotypes.lonspi, gi.lonspi.only@strata$species, 
                 n.pca = 1:150, training.set = 0.9, # ptop of data to be used as training set (for groups >= 10 ind.)
                 result = "overall", # groupMean = average prediction success per group; overall = overall prediction success
                 center = F, scale = FALSE,
                 n.rep = 100, xval.plot = TRUE,
                 parallel="multicore", ncpus = 16)
dev.off()
save(xval,file="xval_lonspi_only.Rdata")