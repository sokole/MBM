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

