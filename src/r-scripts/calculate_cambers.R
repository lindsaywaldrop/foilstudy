# Functions required to create cambers in R

dir.create("./data/camber_files", showWarnings = FALSE)

source("./src/r-scripts/foil_functions.R")

# Load list of cambers
cambers <- read.csv("./data/camber_list.csv", col.names = c("number", "camber"))

# Load template foil
foil <- read.table("./data/AS6091.dat", skip = 0)
colnames(foil) <- c("x", "y")
npts <- nrow(foil)

for(i in 1:nrow(cambers)){
  camber_new <- cambers$camber[i]
  foil_midline_sm <- find_midline(foil$x, foil$y, smoothed = T)
  new_midline <- adjust_midline(foil_midline_sm, camber_new, plot = F)
  new_foil <- create_new_foil(new_midline, plot = F)
  write.table(new_foil, file = paste0("./data/camber_files/AS_camber_", i, ".dat"), 
              row.names = F, col.names = F)
}

