# User specified model for use with gpc expansion code

# polynomial model (ex1 from Sudret 2008)
# Y = 1/2^N * PROD_i ( 3*x_i^2 + 1 )
#
#       NOTE: 1. x_i in [-1,1]
#             2. Validated against Sudret for 3 parameters
#
#       Sobol Indices (for 3 parameters):
#          -> 1st Order: all 0.274725274725275
#          -> 2nd Order: all 0.0549450549450549
#          -> 3rd Order: 0.010989010989011
#
#-------------------------------------------------------------------------
#   N = length(VEC);
#   prod = 1;
#   for i=1:N
#       prod = prod * ( 3*( VEC(i) )^2 + 1 );
#   end
#   val =  1/2^N * prod;
#  
#  
# -------------------------------------------------------------------------
# MODEL: Ishigami function (ex2 from Sudret 2008)
#
#       Y = sin(x1) + a*sin^2(x2) + b*x3^4*sin(x1)
#
#       NOTE: 1. a,b model parameters (not varied)
#             2. x1,x2,x3 in [-pi,pi]
#             2. Validated against Sudret for 3 parameters
#
#       ACTUAL SOBOL INDICES (for 3 parameters):
#          -> 1ST-ORDER INDICES:  [0.3138, 0.4424, 0]
#          -> 2ND-ORDER INDICES:  [0 0.2436 0]
#          -> s123-INDEX:         [0]
#          -> TOTAL-ORDER INDICES:[0.5574, 0.4424, 0.2436] 
#
# --------------------------------------------------------------------------

user_specified_model <- function(vec){
  
  # Define set parameters of Ishigami Function
  a <- 7
  b <- 0.1
  
  # Transform VEC values from [-1,1] -> [-pi,pi] = [A,B]
  A <- -pi
  B <- pi
  vec_trans <- transform_values_to_desired_interval(vec, A, B)
  
  # Evaluate Ishigami function
  val <- sin(vec_trans[1]) + a*(sin(vec_trans[2]))^2 + (b*(vec_trans[3])^4) * sin(vec_trans[1])
  
  return(val)
}

transform_values_to_desired_interval <- function(vec, A, B){
  # Setup linear system
  MAT <- matrix(c(-1, 1, 1, 1), ncol = 2, nrow = 2, byrow = TRUE)
  RHS <- matrix(c(A, B), nrow = 2, ncol = 1)
  coeffs <- solve(MAT) %*% RHS
  
  m <- coeffs[1]
  b <- coeffs[2]
  
  vec_trans <- m*vec + b
  return(vec_trans)
}