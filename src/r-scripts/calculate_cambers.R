# Functions required to create cambers in R

#### Loading functions and template foil ####

source("./src/r-scripts/foil_functions.R")

dir.create("./results", showWarnings = F)

# Load template foil
foil <- read.table("./data/airfoils/AS6091.dat", skip = 0)
colnames(foil) <- c("x", "y")

#### GPC parameter set ####

dir.create("./data/airfoils/gPC_camber_files", showWarnings = FALSE)
dir.create("./results/gPC_noLogRe_matfiles/", showWarnings = FALSE)
dir.create("./results/gPC_LogRe_matfiles/", showWarnings = FALSE)

# Load parameters and extract cambers
gpc_params <- read.csv("./data/parameters/gPC_Params_noLogRe.csv",
                    header = FALSE,
                    col.names = c("Re", "aoa", "camber"))
cambers <- factor(gpc_params$camber)

# Produce and save each camber file:
for(i in levels(cambers)){
  camber_new <- as.numeric(i)
  foil_midline_sm <- find_midline(foil$x, foil$y, smoothed = T, plot = F)
  new_midline <- adjust_midline(foil_midline_sm, camber_new, plot = F)
  new_foil <- create_new_foil(new_midline, plot = F)
  write.table(new_foil, file = paste0("./data/airfoils/gPC_camber_files/AS_camber_", i, ".dat"), 
              row.names = F, col.names = F)
}

#### NN parameter set ####

dir.create("./data/airfoils/NN_camber_files", showWarnings = FALSE)
dir.create("./results/NN_noLogRe_matfiles/", showWarnings = FALSE)
dir.create("./results/NN_LogRe_matfiles/", showWarnings = FALSE)

# Load parameters and extract cambers
nn_params <- read.csv("./data/parameters/NN_Params_noLogRe.csv",
                       header = FALSE,
                       col.names = c("Re", "aoa", "camber"))
cambers <- nn_params$camber

# Produce and save each camber file:
for(i in cambers){
  camber_new <- i
  foil_midline_sm <- find_midline(foil$x, foil$y, smoothed = T, plot = F)
  new_midline <- adjust_midline(foil_midline_sm, camber_new, plot = F)
  new_foil <- create_new_foil(new_midline, plot = F)
  write.table(new_foil, file = paste0("./data/airfoils/NN_camber_files/AS_camber_", i, ".dat"), 
              row.names = F, col.names = F)
}

