#devtools::install_github("bwringe/parallelnewhybrid")
library(dartR)

load(file="gl.lonspischxm_50miss_unl.Rdata")
load(file="gl.lonspischxm.Rdata")

gl.lonspixm<-gl.lonspi.filt[gl.lonspi.filt@strata$species != "schmidtiana"]

gl.lonspixm<-gl.lonspischxm[gl.lonspischxm@strata$species != "schmidtiana"]
gl.lonspixm@strata$species<-as.character(gl.lonspixm@strata$species)
gl.lonspixm@strata$species<-as.factor(gl.lonspixm@strata$species)
gl.lonspixm@strata$species


# population = species

non_pure_ind<-c("lo15t_Ukr2","lo37t_Slo2","xm1t_Ukr20","lo84t_Ukr18","sp73t_Ukr19",
                "lo65x_Ger4","lo69x_Ger6","lo39x_Fra9","lo53x_Ger1","lo48t_Ukr5","lo98t_Rus20","lo24t_Hun5","lo85t_Ukr19","lo58t_Slo4","lo28t_Hun8","lo20t_Mon1","lo97t_Rus19","lo82t_Ukr17","lo74t_Rus8","lo63d_Rus4","lo8d_Hun1","lo18x_Fra5","lo83t_Ukr17",
                "sp1x_Rus1","sp77d_Pol2","sp6d_Rus2","sp17t_Rus3","fi9t_Ukr4","sp61d_Ukr8","sp76t_Pol1","sp16t_Rus3","sp20t_Fra4","sp12t_Hun3","sp11t_Hun3","sp50t_Aus6","sp13t_Hun3")

# obtained with structure resutls: non_adm_inds<-c(row.names(q_matrix_ord)[1:40],row.names(q_matrix_ord)[116:156])
non_adm_inds<-c("lo54t_Slo4",  "lo43t_Slo3",  "lo21t_Hun5",  "lo10d_Hun2",  "lo75d_Ukr13", "lo33t_Slo2" ,
                "lo44t_Slo3",  "lo56t_Slo4",  "lo22t_Hun5",  "lo27t_Hun8" , "lo7d_Hun1"  , "lo57t_Slo4" ,
                "lo12d_Hun2",  "lo23t_Hun5",  "lo46t_Slo3",  "lo29t_Hun8",  "lo25t_Hun5",  "lo30t_Hun8",
                "lo14d_Hun2",  "lo88x_cul1",  "lo89t_Ukr22", "lo34t_Slo2" , "lo64x_Ger3" , "lo36t_Slo2" ,
                "lo45t_Slo3",  "lo38t_Slo2"  ,"lo78x_Ukr15", "lo66x_Ger4" , "lo2x_Fra1" ,  "lo55t_Slo4" ,
                "lo26t_Hun8" , "lo95d_Fin2",  "lo86t_Ukr21" ,"lo49x_Fra11", "lo5d_Hun1" ,  "lo94d_Fin1" ,
                "lo3x_Fra1" ,  "lo11d_Hun2",  "lo67t_Ukr12" ,"lo32x_Fra8"  ,"sp60d_Ukr7",  "sp41t_Hun11",
                "sp84t_Pol3" , "sp18d_Fra4",  "sp5d_Swi1"   ,"sp23t_Hun4" , "sp43t_Hun11", "sp44t_Hun11",
                "sp45t_Hun11" ,"sp27t_Hun6",  "sp31t_Hun6"  ,"sp34t_Hun6" , "sp54t_Slo1",  "sp56t_Slo1" ,
                "sp65x_Ger8" , "sp47t_Aus6",  "sp48t_Aus6"  ,"sp51t_Aus6" , "sp3t_Ser21" , "sp9t_Rom2"  ,
                "sp79x_cul2" , "sp19d_Fra4" , "sp28t_Hun6" , "sp53t_Slo1" , "sp4d_Fra2"  , "sp10d_Aus2" ,
                "sp24t_Hun4" , "sp46t_Hun11" ,"fi7t_Hun10"  ,"fi2t_Hun9"  , "fi4t_Hun9" ,  "fi8t_Hun10" ,
                "sp36t_Hun7" , "sp26t_Hun6" , "sp29t_Hun6"  ,"sp32t_Hun6",  "sp35t_Hun6"  ,"sp49t_Aus6",
                "fi5t_Hun9"  , "sp57t_Slo1",  "sp39t_Aus5")

# selecting pure individuals using STRUCTURE results ####
# loading qmatrix obtained from STRUCTURE
q_matrix_K_n<-read.table("qMatr_lonspixm_50miss_K2_clumped.txt",sep="\t")
sample.data<-read.table("sample_data_lonspixm_50miss_unl.csv",sep=" ", header=T)
nrow(q_matrix_K_n)
nrow(sample.data)
row.names(q_matrix_K_n)<-sample.data$new_id

# order by probability K=2
q_matrix_ord<-q_matrix_K_n[order(q_matrix_K_n$V1),]

q_matr_pure_spi<-q_matrix_ord[q_matrix_ord[,1]>0.8,] # select as pure individuals with less than 1% of admixture coeff.
q_matr_pure_lon<-q_matrix_ord[q_matrix_ord[,2]>0.8,]
pure_spi_ids<-row.names(q_matr_pure_spi)
pure_lon_ids<-row.names(q_matr_pure_lon)
pure_inds<-c(pure_spi_ids,pure_lon_ids)

ind_species<-gl.lonspixm@strata$species
ind_species<-as.character(ind_species)
for(i in 1:length(ind_species)){
  if(!gl.lonspixm[i]@strata$new_id %in% pure_inds){
    ind_species[i]<-"admixed"
  }
}
write.table(ind_species,file="ind_species.csv")
gl.lonspixm2<-gl.lonspixm
pop(gl.lonspixm2)<-as.factor(ind_species)

# sample_cat<-data.frame(nrow = 157, ncol=2)
# sample_cat[,1]<-t(gl.lonspixm@strata$new_id)
# sample_cat[,2]<-pop(gl.lonspixm2)

gl.nhybrids(gl.lonspixm2, outfile = "newHybrids_200AvgPICSNPs4_lonspixm_50miss.txt", outpath = getwd(),
            p0 = "longifolia",#list of populations to be regarded as parental population 0 [default NULL]
            p1 = "spicata",
            method = "AvgPIC" #specifies the method (random or AvgPIC) to select 200 loci for NewHybrids [default random]
)




library(parallelnewhybrid)
#if (!require("pacman")) install.packages("pacman")
pacman::p_load(parallel, plyr, stringr, tidyr)
setwd("~/softwares/parallel_newhybrids")
data_path<-paste0(getwd(),"/unl_2/") # change last number to test different random SNPs subsets
# data_path<-paste0(getwd(),"/Ver_unlinked_from_PGDSpider/")
# data_path<-paste0(getwd(),"/Ver_unlinked/")

parallelnh_LINUX(data_path, where.NH="R/",
                 burnin=10000, sweeps=100000)
