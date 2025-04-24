# Grid functions 
# Custom functions required for generating points using the grid method

mesh_array <- function(x, y, z){
  nptsx <- length(x)
  nptsy <- length(y)
  nptsz <- length(z)
  answer <- list(x = array(x, dim = c(nptsx, nptsy, nptsz)), 
                 y = aperm(array(y, dim = c(nptsx, nptsy, nptsz)), c(2, 1, 3)), 
                 z = aperm(array(z, dim = c(nptsx, nptsy, nptsz)), c(3, 2, 1)))
  return(answer)
}
