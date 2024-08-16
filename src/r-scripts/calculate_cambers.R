# Functions required to create cambers in R

dir.create("./data/camber_files", showWarnings = FALSE)

ASwing <- read.table("./data/NACA3320.dat")
colnames(ASwing) <- c("x", "y")
#ASwing <- read.csv("./data/AS6091.csv", header = TRUE)

camber.original <- (max(ASwing$y) - min(ASwing$y)) / (max(ASwing$x) - min(ASwing$x))

#### Maybe unneeded ####
#sink("./data/camber_files/NACA3320_base.dat")
#cat("NACA3320\ camber=", camber.original, "\n")
#sink()

#write.table(ASwing, file = "./data/camber_files/NACA3320_base.dat",
#            append = TRUE, col.names = FALSE, row.names = FALSE, sep = "  ")

camber_list <- read.csv("./data/camber_list.csv", row.names = 1)
levels.camber <- camber_list$x

#####


for (i in 1:length(levels.camber)){
  camber.new <- levels.camber[i]
  message("camber ",i,"= ",camber.new)
  ASnew.y <- ((ASwing$y - min(ASwing$y)) * camber.new )/ (max(ASwing$y) - min(ASwing$y))
  plot(ASwing$x,ASwing$y,col="blue")
  points(ASwing$x,ASnew.y,col="red")
  camber.check<-(max(ASnew.y)-min(ASnew.y))/(max(ASwing$x)-min(ASwing$x))
  message(all.equal(camber.new,camber.check))
  ASnew<-data.frame(x=ASwing$x,y=ASnew.y)
  
  sink(paste("./data/camber_files/NACA3320_cam",i,".dat",sep=""))
  cat("NACA3320\ camber=",camber.check,"\n")
  sink()
  #write.table(c("wing"),file="birdwing.dat",append=FALSE,col.names=FALSE,row.names=FALSE,sep="  ")
  write.table(ASnew, file = paste("./data/camber_files/NACA3320_cam", i, ".dat", sep = ""), 
              append = TRUE, col.names = FALSE, row.names = FALSE, sep = "  ")
}

#ASwing_cam1 <- read.table("./data/camber_files/NACA3320_base.dat", skip = 1, sep = "\t")
#colnames(ASwing_cam1) <- c("x", "y")
#plot(ASwing_cam1)
