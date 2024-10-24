#gPC functions

Legendre_roots <- function(degree, test_flag = 0){
  # get roots of (n+1)^(st) Legendre polynomial
  
  #get coefficients of the nth degree Legendre polynomial
  coeffs <- Legendre_poly_coeffs(degree)
  # find roots of polynomial with those coefficients
  # Method 2
  leg_roots <- polyroots(coeffs, n = 1e-16)
  # Method 1
  #leg_roots <- roots(coeffs)
  #leg_roots <- sort(leg_roots)
  leg_roots <- leg_roots[,1]
  if(test_flag == 1){
    poly_roots <- leg_roots
    poly_roots_matlab <- read.csv("./test/gpc_matlab/poly_roots.csv", header = F)
    root_diffs <- test_matlab_match(poly_roots, poly_roots_matlab, 1e-14)
    poly_roots <- poly_roots[root_diffs]
    print(all.equal(poly_roots_matlab[,1], poly_roots, tolerance = 1e-10))
    return(poly_roots)
  } else{
    return(leg_roots)
  }
}

Legendre_poly_coeffs <- function(n){
  # get Legendre Polynomial Coefficients: 
  # L_n(x) = c0 + c1*x + c2*x^2 +...+ c_n*x^n
  # input: n <- degree of Legendre poly
  
  if(n == 0){
    coeffs <- 1   # L_0(x) = 1
  }else if(n == 1){
    coeffs = c(1, 0)   # L_1(x) = x
  }else{
    # create Legendre polynomials based off of recurrence relation
    p_nm1 <- 1
    p_n <- c(1, 0)
    for (i in 1:(n-1)){
      p_np1 <- ((2*i + 1)*c(p_n, 0) - i*c(0,0, p_nm1))/(i + 1)
      p_nm1 <- p_n
      p_n <- p_np1
    }
    coeffs <- p_np1
  }
  return(coeffs)
}

compute_all_collo_pt_combos <- function(n, poly_roots, test_flag = 0){
  # combine all plausible collocation points, i.e., combine
  # roots of Legendre polynomial --> uses same machinery
  # as creating the alphaMAT indices of the polynomial.
  
  # Start with column vector of POLY ROOTS
  root_mat <- as.matrix(poly_roots)
  num_roots <- length(poly_roots)
  
  # Loop over how many qualities we're varying 
  # (minus 1, since starting with  a root vector already)
  for(k in 1:(n-1)){
    # empty vector for new leftmost vector in root_mat
    left_add <- NULL
    for(j in 1:num_roots){
      # create a vector of all 'j' values of length of previous root matrix
      left_aux <- poly_roots[j] * rep(1, nrow(root_mat))
      # concatenate new left-most vector for root_mat
      left_add <- c(left_add, left_aux)
    }
    
    # create empty array to stack root_mat in a certain # of times
    root_stack <- NULL
    # STACK the rootMAT we already have!
    for (j in 1:num_roots){
      # stack the pre-existing root_mat from before
      root_stack <- rbind(root_stack, root_mat)
      #This is where I am stuck, look at the structure of rootMAT and how it grows in the matlab code
    }
    # REDEFINE root_mat to include new left_vector and stacked root_mat!
    root_mat <- cbind(left_add, root_stack)
  }
  if(test_flag == 1){
    param_combo_matlab <- as.matrix(read.csv("./test/gpc_matlab/param_combo.csv", header = F))
    print(all.equal(root_mat, param_combo_matlab, tolerance = 1e-10, check.attributes = F))
  }
  return(root_mat)
}

sample_parameter_combos <- function(n_subset, param_combo, test_flag){
  # subsample all parameter combinations by taking the 
  # N_subset closest number to origin (parameters are all 
  # scaled to [-1,1] here)
  distance_vec <- NULL
  for(i in 1:nrow(param_combo)){
    # compute overall distance for origin, with points in particular row
    # of parameter combination
    sum_dat <- 0
    for(j in 1:ncol(param_combo)){
      sum_dat <- sum_dat + ((param_combo[i, j])^2)
    }
    #distance from origin
    distance_vec[i] <- sum_dat
  }
  
  if(test_flag == 1){
    matlab_distanceVec <- t(as.matrix(read.csv("./test/gpc_matlab/distanceVec.csv", header = F)))
    print(all.equal(distance_vec, matlab_distanceVec[,1], tolerance = 1e-10, 
                    check.attributes = F))
  }
  
  # Re-order the data in distance_vec, get orig inds
  inds <- order(round(distance_vec, 10))
  
  if(test_flag == 1){
    matlab_inds <- t(as.matrix(read.csv("./test/gpc_matlab/inds.csv", header = F)))
    print(all.equal(inds, matlab_inds[,1], tolerance = 1e-10, check.attributes = F))
  }
  
  param_combo <- param_combo[inds,]
  
  # take first N_subset
  param_combo_subset <- param_combo[1:n_subset,]
  
  if(test_flag == 1){
    param_combo_subset_matlab <- as.matrix(read.csv("./test/gpc_matlab/param_combo_subset.csv", header = F))
    diff_subset <- test_matlab_match(param_combo_subset, param_combo_subset_matlab, 1e-10)
    print(all.equal(param_combo_subset, param_combo_subset_matlab, tolerance = 1e-10, 
                    check.attributes = F))  
  }
  
  
  return(param_combo_subset)
}

create_polynomial_ordering <- function(n, p, test_flag){
  # get convention for how to setup multivariable Legendre
  #           polynomial
  
  q <- n # number of varies parameters (qualities)
  d <- p # highest degree basis polynomial
  
  # Start with column vector of poly orders {0,1,2, ... d}
  alpha_mat <- matrix(seq(0, d, by = 1), ncol = 1)
  
  # Loop over how many qualities we're varying 
  # (minus 1, since starting with  a single alpha vector already)
  for(i in 1:(q-1)){
    # empty vector for 'new' leftmost vector in alpha_mat
    left_add <- NULL
    # CREATE left-most vector --> loop over each plausible degree
    for(j in 0:d){
      # create a vector of all 'j' values of length of previous alpha matrix
      left_aux <- j*matrix(rep(1, nrow(alpha_mat), ncol = 1))
      # concatenate new left-most vector for alphaMAT
      left_add <- rbind(left_add, left_aux)
    }
   
    # create empty array to stack alphaMAT in a certain # of times
    alpha_stack <- NULL
    # STACK the alphaMAT we already have!
    for(j in 0:d){
      #stack the pre-existing alpha_mat from before
      alpha_stack <- rbind(alpha_stack, alpha_mat)
    }
    #REDEFINE alpha_mat to include new leftVector and stacked alpha_mat!
    alpha_mat <- cbind(left_add, alpha_stack)
  }
  # CUT OUT COMBINATIONS WHERE SUM INDICES > d
  alpha_new <- NULL
  for(i in 1:nrow(alpha_mat)){
    if(sum(alpha_mat[i,]) <= p){
      alpha_new <- rbind(alpha_new, alpha_mat[i,])
    }
  }
  
  if(test_flag == 1){
    alpha_mat_matlab <- as.matrix(read.csv("./test/gpc_matlab/alpha_new.csv", header = F))
    print(all.equal(alpha_new, alpha_mat_matlab, tolerance = 1e-10, 
                    check.attributes = F))
  }
  
  return(alpha_new)
}

create_info_matrix <- function(n, p, cap_p, param_combo_subset, alpha_mat, test_flag){
  # create INFORMATION MATRIX
  info_mat <- matrix(NA, ncol = cap_p, nrow = nrow(param_combo_subset))
  # Loop over all test points
  for(j in 1:nrow(param_combo_subset)){
    param_vec <- param_combo_subset[j, ]
    # Loop over MULTIVARIABLE LEGENDRE POLY
    for(i in 0:(cap_p - 1)){
      # get Multivariable Legendre indices (i+1 bc i starts at 0)
      alpha_vec <- alpha_mat[i + 1,]
      # compute Multivariable Legendre Poly. at specific point, 
      # paramVec (i+1 bc i starts at 0)
      info_mat[j, i + 1] <- multi_dim_Legendre_poly(alpha_vec, param_vec)
    }
  }
  
  if(test_flag == 1){
    info_mat_matlab <- as.matrix(read.csv("./test/gpc_matlab/info_mat.csv", header = F))
    print(all.equal(info_mat, info_mat_matlab, tolerance = 1e-10, check.attributes = F))
  }
  return(info_mat)
}

multi_dim_Legendre_poly <- function(alpha_vec, param_vec){
  # compute multidimensional Legendre Polynomial and evaluate for 
  # parameter combination in the vector param_vec
  val <- 1
  for(k in 1:length(alpha_vec)){
    b <- alpha_vec[k] # Get degree with specific single dimensional Leg. Poly
    x <- param_vec[k] # Get value to evaluate single dimensional Leg. Poly at
    # Compute product of all Legendre Poly's evaluated at specific x and of
    # particular degree
    val <- val * Legendre_poly(b, x)
  }
  return(val)
}

Legendre_poly <- function(b, x){
  # compute nth Legendre polynomial evaluated at x
  # --> These polynomials are to be evaluated on [-1,1]
  # --> Hardcoded Legendre Polynomials are faster to evaluate than
  #     either Rodriques Formula or MATLAB's Legendre poly function
  # from https://en.wikipedia.org/wiki/Legendre_polynomials
  if(b == 0){
    val <- 1
  } else if (b == 1){
    val <- x
  } else if (b == 2){
    val <- 0.5*((3*x^2) - 1)
  } else if (b == 3){
    val <- 0.5*((5*x^3) - (3*x))
  } else if (b == 4){
    val <- 0.125*((35*x^4) - (30*x^2) + 3)
  } else if (b == 5){
    val <- 0.125*((63*x^5) - (70*x^3) + (15*x))
  } else if (b == 6){
    val <- (0.0625)*((231*x^6) - (315*x^4) + (105*x^2) - 5)
  } else if (b == 7){
    val <- (0.0625)*((429*x^7) - (693*x^5) + (315*x^3) - (35*x))
  } else if (b == 8){
    val <- (1/128)*((6435*x^8) - (12012*x^6) + (6930*x^4) - (1260*x^2) + 35)
  } else if (b == 9){
    val <- (1/128)*((12155*x^9) - (25740*x^7) + (18018*x^5) - (4620*x^3) + (315*x))
  } else if (b == 10){
    val <- (1/256)*((46189*x^10) - (109395*x^8) + (90090*x^6) - (30030*x^4) + 
                      (3465*x^2) - 63)
  } else if (b == 11){
    val <- (1/256)*(x)*(-693 + (15015*x^2) - (90090*x^4) + (218790*x^6) - (230945*x^8) + 
                          (88179*x^10))
  } else if (b == 12){
    val <- (1/1024)*(231 - (18018*x^2) + (225225*x^4) - (1021020*x^6) + (2078505*x^8) - 
                       (1939938*x^10) + (676039*x^12))
  } else if (b == 13){
    val <- (1/1024)*(x)*(3003 - (90090*x^2) + (765765*x^4) - (2771340*x^6) + (4849845*x^8) - 
                           (4056234*x^10) + (1300075*x^12))
  } else if (b == 14){
    val <- (1/2048)*(-429 + (45045*x^2) - (765765*x^4) + (4849845*x^6) - (14549535*x^8) + 
                       (22309287*x^10) - (16900975*x^12) + (5014575*x^14))
  } else {
    stop("Please see the MATLAB version to solve this polynomial")
  }
  return(val)
}

compute_validation_and_testing_error <- function(n, s_coeffs, alpha_mat, training_params,
                                                 sobol_seed = 0){
  require(spacefillr)
  # Compute Validation and Testing Error across full 3D parameter space
  # -> Samples full space using Sobol' sequence
  #  INPUTS:  n: # of uncertain parameters
  #           s_coeffs: gPC expansion coefficients
  #           alpha_mat: Legendre polynomial orderings
  #           param_combo_subset: sampled parameter combinations for training
  
  npts <- 1000
  
  # Package 1: spacefillr
  testing_params <- spacefillr::generate_sobol_set(npts, n, seed = sobol_seed)
  
  # Package 2: SobolSequence
  #testing_param <- SobolSequence::sobolSequence.points(3,25,1000)
  
  testing_params <- 2 * testing_params - 1
  
  # Training data recovery
  train_list <- compute_gpc_and_true_function_values(training_params, s_coeffs, alpha_mat)
  
  # Testing data recovery
  test_list <- compute_gpc_and_true_function_values(testing_params, s_coeffs, alpha_mat)

  return(list("training" = train_list, "testing" = test_list))
}

compute_gpc_and_true_function_values <- function(params, s_coeffs, alpha_mat){
  # compute the gPC predicted value and real function value at all parameter 
  # combinations in the training_params set
  
  # allocate space
  gpc_dat <- rep(NA, nrow(params))
  real_dat <- rep(NA, nrow(params))
  
  for(i in 1:nrow(params)){
    
    gpc_dat[i] <- eval_multidim_Legendre_poly(s_coeffs, alpha_mat, params[i,])
    
    real_dat[i] <- user_specified_model(params[i,])
  }
  return(list("gpc_dat" = gpc_dat, "real_dat" = real_dat))
}

eval_multidim_Legendre_poly <- function(s_coeffs, alpha_mat, vec){
  # Evaluate Multi-dimensional Legendre Polynomial at specific x, y, z
  sum_dat <- 0
  for(j in 1:nrow(alpha_mat)){
    prod_dat <- 1
    for(k in 1:ncol(alpha_mat)){
      prod_dat <- prod_dat * Legendre_poly(alpha_mat[j, k], vec[k])
    }
    sum_dat <- sum_dat + s_coeffs[j] * prod_dat
  }
  return(sum_dat)
}

test_gpc_expansion <- function(s_coeffs, alpha_mat){
  x_vec <- seq(-0.95, 0.95, length.out = 50)
  y_vec <- seq(-0.95, 0.95, length.out = 50)
  xy_mesh <- meshgrid(x_vec, y_vec)
  X <- xy_mesh$X
  Y <- xy_mesh$Y
  
  z_test <- matrix(NA, nrow = length(y_vec), ncol = length(x_vec))
  z_real <- matrix(NA, nrow = length(y_vec), ncol = length(x_vec))
  
  for(i in 1:length(x_vec)){
    for(j in 1:length(y_vec)){
      
      # Get specific (x,y,z) value to test in expansion 
      x <- x_vec[i]
      y <- y_vec[j]
      z <- 0.9     # choose a specific 2D subspace
      
      
      # Store in vector to pass into expansion
      vec <- c(x, y, z)
      
      # Evaluate multidimensional Legendre polynomial
      z_test[j, i] <- eval_multidim_Legendre_poly(s_coeffs, alpha_mat, vec)
      
      # Test problem
      z_real[j, i] <- user_specified_model(vec)
    }
  }
  test_df <- data.frame("x" = as.vector(X), "y" = as.vector(Y), 
                        "z" = as.vector(z_test))
  real_df <- data.frame("x" = as.vector(X), "y" = as.vector(Y), 
                        "z" = as.vector(z_real))
  z_range <- range(c(range(z_real), range(z_test)))
  p_test <- ggplot(test_df, aes(x, y, fill = z)) + 
    geom_tile() + 
    scale_fill_viridis_c(limits = z_range) +
    ggtitle("Test function") +
    theme_bw()
  p_real <- ggplot(real_df, aes(x, y, fill = z)) + 
    geom_tile() + 
    scale_fill_viridis_c(limits = z_range) +
    ggtitle("Real function") +
    theme_bw()
  return(list(p_test, p_real))
}

compute_sobol_indices <- function(s_coeffs, alpha_mat){
  # Compute Sobol Indices from aPC Expansion Coefficients!
  # INPUTS: (1) Coeffs - aPC expansion coefficients vector
  #         (2) PolynomialDegree - alpha indices for expansion
  #               - each row different expansion coefficient combo for basis
  #                 poly's, e.g., row=[1,0,5] --> L_1(x1)*L_0(x2)*L_5(x3)
  
  p <- nrow(alpha_mat)  # total number of terms in gPC expansion
  n <- ncol(alpha_mat)  # number of uncertain parameters
  
  # VAR:  SUM_{j=2}^N Coeffs_j^2 * E[ PSI_j(x1,x2,..,xN)^2] 
  #               (j index in MATLAB notation)
  
  e_psi_sqr_vec <- compute_expectation_psi_squared(alpha_mat)
  
  sigma_sqr <- sum((s_coeffs[-1])^2 * e_psi_sqr_vec[-1])
  
  # initialization
  s_first <- rep(NA, n)
  
  # compute first-order index
  for (i in 1:n){
    # Loop over every row in s_coeffs / alpha index matrix for a particular uncertain
    #   parameter, x_i
    ct <- 0
    k_vec <- NULL
    s_aux <- 0
    for(j in 1:p){
      # determine if all other indices are zero (or not) in alpha index matrix:
      use_this <- 1
      for(k in 1:n) {
        if(k != i){
          if(alpha_mat[j,k] != 0) use_this <- 0
        } 
      }
      
      # if only non-zero index corresponds to parameter of interest
      if(alpha_mat[j,i] != 0 && use_this == 1){
        s_aux <- s_aux + s_coeffs[j]^2 * e_psi_sqr_vec[j]
        k_vec <- c(k_vec, j)
        ct <- ct + 1
      }
    }
    s_aux <- s_aux / sigma_sqr
    
    s_first[i] <- s_aux
  }
  
  # Compute TOTAL-Order Index
  
  s_total <- rep(NA, n)
  
  for (i in 1:n){  # Loop over every uncertain parameter
    # Loop over every row in Coeffs / alpha index matrix for a particular
    #  uncertain parameter, x_i
    ct <- 0
    k_vec <- NULL
    s_aux <- 0
    for(j in 1:p){
      # if parameter has non-zero index in alpha_mat
      if(alpha_mat[j,i] != 0){
        s_aux <- s_aux +((s_coeffs[j]^2) * e_psi_sqr_vec[j])
        k_vec <- c(k_vec, j)
        ct <- ct + 1
      }
    }
    s_aux <- s_aux / sigma_sqr
    s_total[i] <- s_aux
  }
  
  # compute s123... Index
  
  s123 <- rep(NA, n)
  
  for(i in 1:n){
    ct <- 0
    k_vec <- NULL
    s_aux <- 0
    for(j in 1:p){
      use_this <- 1
      for(k in 1:n) if(alpha_mat[j,k] == 0) use_this <- 0
      if(use_this == 1){
        s_aux <- s_aux + ((s_coeffs[j]^2) * e_psi_sqr_vec[j])
        k_vec <- c(k_vec, j)
        ct <- ct + 1
      }
    }
    s_aux <- s_aux / sigma_sqr
    s123[i] <- s_aux
  }
  
  # Compute Second-Order Index
  inds_vec <- 1:n
  s2nd <- matrix(NA, nrow = n, ncol = n)
  for(i1 in 1:n){
    # first index
    ind1 <- i1
    start_i2 <- i1 + 1
    if(start_i2 > n) break
    for(i2 in start_i2:n){
      # second index
      ind2 <- i2
      
      # reinitiate-Auxillary indices vector
      inds_aux <- rep(0, length(inds_vec))
      
      # store indices being compared and find indices that must be zero
      inds_aux[ind1] <- ind1
      inds_aux[ind2] <- ind2
      inds_test <- which((inds_vec - inds_aux) != 0)
      
      # Loop over every row in Coeffs / alpha index matrix for a particular
      #       uncertain parameter, x_i
      ct <- 0
      k_vec <- NULL
      s_aux <- 0
      for(j in 1:p){
        use_this <- 1
        for(k in 1:length(inds_test)) if(any(alpha_mat[j,inds_test] != 0)) use_this <- 0
        
        # if only non-zero indices correspond to parameters of interest
        if(alpha_mat[j,i1] != 0 & alpha_mat[j, i2] != 0 & use_this == 1){
          s_aux <- s_aux + ((s_coeffs[j]^2) * e_psi_sqr_vec[j])
          k_vec <- c(k_vec, j)
          ct <- ct + 1
        }
      }
      s_aux <- s_aux / sigma_sqr
      s2nd[i1, i2] <- s_aux
    }
  }
  return(list("s_first" = s_first, "s_total" = s_total, "s2nd" = s2nd, 
              "s123" = s123))
}

compute_expectation_psi_squared <- function(alpha_mat){
  # compute the Expectation of PSI^2 and store each into vector
  #                          ( E[ PSI_j^2(x) ] )
  
  # Allocate space
  e_psi_sqr_vec <- rep(NA, nrow(alpha_mat))
  for(i in 1:nrow(alpha_mat)){
    #  Loop over each index in alphaMAT (loops over # of uncertain parameters)
    prod_dat <- 1
    for(j in 1:ncol(alpha_mat)){
      prod_dat <- (prod_dat * 1) / (2*alpha_mat[i,j] + 1)
    }
    e_psi_sqr_vec[i] <- prod_dat
  }
  return(e_psi_sqr_vec)
}

test_matlab_match <- function(r_values, matlab_values, tolerance, 
                              print_diff = FALSE){
  r_values <- as.matrix(r_values)
  matlab_values <- as.matrix(matlab_values)
  order_vec <- rep(NA, nrow(r_values))
  for(i in 1:nrow(r_values)){
    R_param <- r_values[i,]
    for(j in 1:nrow(matlab_values)){
      diff <- abs(R_param - matlab_values[j,])
      if(is.numeric(print_diff) & i == print_diff) print(diff)
      if(all(diff < tolerance)) order_vec[i] <- j
    }
  }
  return(order_vec)
}
