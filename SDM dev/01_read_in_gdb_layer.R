#' @title MBM SDM workflow dev with dummy data
#' @author Eric R. Sokol
#' 
#' @references 
#' Jerosch, K., Scharf, F. K., Deregibus, D., Campana, G. L., Zacher, K., Pehlke, H., et al. (2019). Ensemble Modeling of Antarctic Macroalgal Habitats Exposed to Glacial Melt in a Polar Fjord. Front. Ecol. Evol. 7. doi:10.3389/fevo.2019.00207.
#' 
#' Wilfried Thuiller, Damien Georges, Maya Gueguen, Robin Engler and Frank Breiner (2021). biomod2: Ensemble Platform for Species Distribution Modeling. R package version 3.5.1. https://CRAN.R-project.org/package=biomod2
 

# How to read a gdb layer here: 
# https://stackoverflow.com/questions/55910398/how-do-i-get-rgdal-to-open-a-geodatabase-gdb-file

library(rgdal)

list.files("~/MBM_working/TestData.gdb")

my_gdb_path <- path.expand("~/MBM_working/TestData.gdb")
ogrListLayers(my_gdb_path)
my_layer <- readOGR(my_gdb_path,"Taylor_Valley_Test")


plot(my_layer)

## creating the sub directory for the data and output figures

dir.create("data")
dir.create("output_figure")


## loading libraries required for the SDM
library(raster)
library(dplyr)
library(rgdal)
library(sf)
library(ggplot2)
library(readxl)
library(sp)
library(biomod2)
library(ggplot2)




# loading environmental variables

pre<-raster::stack("C:/Users/Khum/Documents/mbm_updated_data/predictors_update/mat_env.grd")


# loading presence and absence data frame created in arcmap

orgmat<-read.csv("C:/Users/Khum/Documents/mbm_output/sdm_rarified/orgmat_thin_combine.csv")

## variables making ready to apply in the SDM model

# Select the name of the studied species
myRespName <- 'OrangeMat'

# Get corresponding presence/absence data
myResp <- as.numeric(orgmat[, 'presence'])

# Get corresponding XY coordinates
myRespXY <- orgmat[, c('ycentroid', 'xcentroid')]




## This is starting of formating the model,running the mdoels


# Format Data with true absences
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = pre,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName)

# plot(myBiomodData) to show presence and absence data in map

# Create default modeling options
myBiomodOptions <- BIOMOD_ModelingOptions()


# Model single models with default algorithsm
myBiomodModelOut <- BIOMOD_Modeling(bm.format = myBiomodData,
                                    modeling.id = 'AllModels',
                                    models = c('RF', 'GLM','GBM','GAM','CTA','ANN','SRE','FDA','MARS','MAXENT.Phillips', 'MAXENT.Phillips.2'),
                                    bm.options = myBiomodOptions,
                                    nb.rep = 2, # how man times to run, this is "run times"
                                    data.split.perc = 80,
                                    metric.eval = c('TSS','ROC'),
                                    var.import = 3, # permutation times to calcualte the variables important, this is "rand time"
                                    do.full.models = FALSE,
                                    seed.val = 42)

myBiomodModelOut


myBiomodEM <- BIOMOD_EnsembleModeling(bm.mod = myBiomodModelOut,
                                      models.chosen = 'all',
                                      em.by = 'all',
                                      metric.select = c('TSS'),
                                      metric.select.thresh = c(0.2),
                                      metric.eval = c('TSS', 'ROC'),
                                      var.import = 3,
                                      prob.mean = TRUE,
                                      prob.median = FALSE,
                                      prob.cv = FALSE,
                                      prob.ci = FALSE,
                                      prob.ci.alpha = 0.05,
                                      committee.averaging = TRUE,
                                      prob.mean.weight = FALSE,
                                      prob.mean.weight.decay = 'proportional',
                                      seed.val = 42)


# Project single models
myBiomodProj <- BIOMOD_Projection(bm.mod = myBiomodModelOut,
                                  proj.name = 'Current',
                                  new.env = pre,
                                  models.chosen = 'all',
                                  metric.binary = 'all',
                                  metric.filter = 'all',
                                  build.clamping.mask = TRUE,
                                  output.format = ".img")
myBiomodProj
plot(myBiomodProj)


BiomodEMProj <- BIOMOD_EnsembleForecasting(bm.em = myBiomodEM,
                                           proj.name = 'CurrentEM',
                                           new.env = pre,
                                           models.chosen = 'all',
                                           metric.binary = 'all',
                                           metric.filter = 'all',
                                           output.format = ".img")
BiomodEMProj
plot(BiomodEMProj)


## getting the outputs and plotting it

en<-raster("C:/Users/Khum/OneDrive - UCB-O365/Mats_sdm/code/OrangeMat/proj_CurrentEM/individual_projections/OrangeMat_EMmeanByTSS_mergedAlgo_mergedRun_mergedData.img")

plot(en/1000)


#### getting response curves of each predictor variables

## response curve
rescur<-bm_PlotResponseCurves(bm.out = myBiomodModelOut,
                              models.chosen = get_built_models(myBiomodModelOut)[c(1:4)],
                              fixed.var = 'median',
                              colors=NA)


# ggsave("C:/Users/Khum/OneDrive - UCB-O365/Mats_sdm/data_output/response_curve.jpg", width=10, height=8, dpi=500)


## getting importna varialbles responsible for the distribution of mats

get_value<-get_variables_importance(myBiomodModelOut, as.data.frame = TRUE)

write.csv(get_value,"C:/Users/Khum/OneDrive - UCB-O365/Mats_sdm/data_output/variable_imp_20run.csv")

var_imp<-read.csv("C:/Users/Khum/OneDrive - UCB-O365/Mats_sdm/data_output/variable_imp_fryxell.csv")

var_imp_fig<-get_value%>%
  group_by(Expl.var)%>%
  summarise(av=mean(Var.imp))%>%
  ggplot()+
  geom_bar(aes(x=reorder(Expl.var,-av),y=av*100), stat="identity")+
  labs(x="Predictor variables",y="% contribution")


# ggsave("C:/Users/Khum/OneDrive - UCB-O365/Mats_sdm/figure_output/var_imp_fryxell.png",
#        width=6,height=4,dpi=500,units="in")


## getting figures for model evaluations
bm_PlotEvalMean(bm.out=myBiomodModelOut)[2]

# ggsave("C:/Users/Khum/OneDrive - UCB-O365/Mats_sdm/figure_output/Model_fig_fryxell.png",
#        width=6,height=4,dpi=500,units="in")



