# NN generating points 
# Custom functions required for generating points using the neural network method

get_training_and_test_data <- function(num_inputs, n_train, n_test, eps_error, sobol_seed, 
                                       testing = FALSE){
  ## Get Training and Testing Data for Neural Network
  
  if(testing){
    training_data <- as.matrix(read.csv("./test/nn_matlab/Matlab_training_preeval.csv", header = F))
    test_data <- as.matrix(read.csv("./test/nn_matlab/Matlab_testing_preeval.csv", header = F))
  } else{
    # TRAINING DATA: Setting Up TRAINING SAMPLES
    ## SOBOL SEQUENCE:
    #   get Sobol' sequence in [0,1]^N
    # Package: spacefillr
    training_data <- spacefillr::generate_sobol_set(n_train, num_inputs, seed = sobol_seed)
    #   get Sobol' sequence in [-1,1]^N
    # training_data <- 2*spacefillr::generate_sobol_set(n_train, num_inputs, seed = sobol_seed) - 1
    
    # Find TEST Data Combinations
    set.seed(sobol_seed)
    test_data <- matrix(runif(n_test*num_inputs), nrow = n_test, ncol = num_inputs)
  }
  
  # Evaluate specified function (to get output values for training and testing data)
  flag_scale_output <- F
  train_output <- evaluate_function(training_data, eps_error, flag_scale_output, 0, 0)
  test_output <- evaluate_function(test_data, eps_error, flag_scale_output, 0, 0)
  
  if(testing){
    train_output_matlab <- as.matrix(read.csv("./test/nn_matlab/Matlab_training_output.csv", header = F))
    test_output_matlab <- as.matrix(read.csv("./test/nn_matlab/Matlab_testing_output.csv", header = F))
    all.equal(train_output, train_output_matlab[,1])
    all.equal(test_output, test_output_matlab[,1])
  }
  
  # Scale output data to [0, 1]
  flag_get_minmax <- T
  if(flag_get_minmax){
    min_z <- min(c(train_output, test_output))
    max_z <- max(c(train_output, test_output))
    m <- 1/(max_z - min_z)
    b <- -m*min_z
    train_output <- m * train_output + b
    test_output <- m * test_output + b
    
    if(testing){
      train_output_matlab <- as.matrix(read.csv("./test/nn_matlab/Matlab_training_output_scaled.csv", 
                                                header = F))
      test_output_matlab <- as.matrix(read.csv("./test/nn_matlab/Matlab_testing_output_scaled.csv", 
                                               header = F))
      all.equal(train_output, train_output_matlab[,1])
      all.equal(test_output, test_output_matlab[,1])
    }
  }
  return(list("training_data" = training_data, 
              "train_output" = train_output, 
              "test_data" = test_data, 
              "test_output" = test_output, 
              "min_z" = min_z, "max_z" = max_z))
}

evaluate_function <- function(dat, eps_error, flag_scale_output, min_z, max_z){
  # Model function
  # eps: small error amount to take into account in model evaluation
  
  # Assumes 1 output
  f_mat <- rep(0, nrow(dat))
  
  # Loop over all input parameter combinations in dat
  for (i in 1:nrow(dat)){
    x1 <- dat[i, 1]
    x2 <- dat[i, 2]
    x3 <- dat[i, 3]
    
    f1 <- (1 - x2^2) * cos(2 * (2*pi*x1))
    f2 <- 0.5 * sin(2 * (2*pi*x2))
    f3 <- (2*x3) + 0.2
    
    f_val = f1 + f2 + f3 + 2;
    
    f_mat[i] <- f_val*(1 + eps_error*(2*runif(1)- 1))
  }
  
  # Scale output data to [0, 1]
  if(flag_scale_output){
    m <- 1/(max_z - min_z)
    b <- -m * min_z
    f_mat <- m * f_mat + b
  }
  return(f_mat)
}

train_artificial_neural_network <- function(training_data, testing_data, hyper_param_list, 
                                            testing = FALSE){
  ## TRAINS an artificial neural network with 2 hidden layers
  flag_use_stored_weights <- F
  if(flag_use_stored_weights){
    print("Loading previously stored weights... trying to make them more accurate...")
  }else{
    print("NOT starting with previously stored weights!")
  }
  
  # Redefine Data for Consistency
  x0 <- t(training_data[,1:hyper_param_list$num_inputs])
  z0 <- t(training_data[,(hyper_param_list$num_inputs+1):ncol(training_data)])
  
  num_input_neurons <- hyper_param_list$num_inputs
  train_data_size <- ncol(x0)
  x0_test <- t(testing_data[,1:hyper_param_list$num_inputs])
  z0_test <- t(testing_data[,(hyper_param_list$num_inputs+1):ncol(training_data)])
  
  # Rule of thumb for # OF HIDDEN LAYER NEURONS (Hidden Layer Size)
  n_input <- num_input_neurons  # number of input neurons
  n_output <- length(z0[,1])    # number of output neurons
  num_train <- train_data_size  # amount of training data
  num_hid_layer_neurons <- hyper_param_list$num_hidden_layer_neurons
  
  # Get Number of DATA OUTPUTs
  num_z <- prod(dim(z0))            # number of data elements in z (training set)
  num_z_test <- prod(dim(z0_test))  # number of data elements in z (testing set)
  
  if(flag_use_stored_weights){
    
    print("Sorry, not yet implemented in R")
    
  }else{
    
    if(testing){
      w1 <- as.matrix(read.csv("./test/nn_matlab/Matlab_w1.csv", header = F))
      w2 <- as.matrix(read.csv("./test/nn_matlab/Matlab_w2.csv", header = F))
      wend <- as.matrix(read.csv("./test/nn_matlab/Matlab_wend.csv", header = F))
      b1 <- as.matrix(read.csv("./test/nn_matlab/Matlab_b1.csv", header = F))
      b2 <- as.matrix(read.csv("./test/nn_matlab/Matlab_b2.csv", header = F))
      bend <- as.matrix(read.csv("./test/nn_matlab/Matlab_bend.csv", header = F))
    }else{
      # Initialize weights/biases randomly
      coeff <- 1/sqrt(num_train)
      w1 <- coeff*(2 * matrix(runif(num_hid_layer_neurons*num_input_neurons), 
                              nrow = num_hid_layer_neurons, ncol = num_input_neurons) - 1)
      w2 <- coeff*(2 * matrix(runif(num_hid_layer_neurons*num_hid_layer_neurons), 
                              nrow = num_hid_layer_neurons, ncol = num_hid_layer_neurons) - 1)
      wend <- coeff*(2 * matrix(runif(n_output*num_hid_layer_neurons), 
                                nrow = n_output, ncol = num_hid_layer_neurons) - 1)
      b1 <- coeff*(2 * matrix(runif(num_hid_layer_neurons), nrow = num_hid_layer_neurons) - 1)
      b2 <- coeff*(2 * matrix(runif(num_hid_layer_neurons), nrow = num_hid_layer_neurons) - 1)
      bend <- coeff*(2 * matrix(runif(n_output), nrow = n_output) - 1)
    }
    
    
  }
  
  #### Initialize gradients of weight/bias matrices ####
  w1_p <- matrix(rep(0, num_hid_layer_neurons*num_input_neurons), 
                 nrow = num_hid_layer_neurons, ncol = num_input_neurons)
  dj_dw1_p <- matrix(rep(0, num_hid_layer_neurons*num_input_neurons), 
                     nrow = num_hid_layer_neurons, ncol = num_input_neurons)
  #
  w2_p <- matrix(rep(0, num_hid_layer_neurons*num_hid_layer_neurons), 
                 nrow = num_hid_layer_neurons, ncol = num_hid_layer_neurons)
  dj_dw2_p <- matrix(rep(0, num_hid_layer_neurons*num_hid_layer_neurons), 
                     nrow = num_hid_layer_neurons, ncol = num_hid_layer_neurons)
  #
  wend_p <- matrix(rep(0, n_output*num_hid_layer_neurons), 
                   nrow = n_output, ncol = num_hid_layer_neurons)
  dj_dwend_p <- matrix(rep(0, n_output*num_hid_layer_neurons), 
                       nrow = n_output, ncol = num_hid_layer_neurons)
  #
  b1_p <- matrix(rep(0, num_hid_layer_neurons), nrow = num_hid_layer_neurons)
  dj_db2_p <- matrix(rep(0, num_hid_layer_neurons), nrow = num_hid_layer_neurons)
  #
  b2_p <- matrix(rep(0, num_hid_layer_neurons), nrow = num_hid_layer_neurons)
  dj_db2_p <- matrix(rep(0, num_hid_layer_neurons), nrow = num_hid_layer_neurons)
  #
  bend_p <- matrix(rep(0, num_hid_layer_neurons), nrow = num_hid_layer_neurons)
  dj_dbend_p <- matrix(rep(0, num_hid_layer_neurons), nrow = num_hid_layer_neurons)
  
  ####  Initialize Learning Rates and Momentum ####
  
  # Weight Learning Rates
  lambda_1 <- hyper_param_list$learning_rate    # initial learning rate
  lambda_2 <- hyper_param_list$learning_rate    # initial learning rate
  lambda_end <- hyper_param_list$learning_rate  # initial learning rate
  # Bias learning rates
  lambda_b1 <- hyper_param_list$learning_rate    # initial learning rate
  lambda_b2 <- hyper_param_list$learning_rate    # initial learning rate
  lambda_bend <- hyper_param_list$learning_rate  # initial learning rate
  # Minimum lambda
  min_lam <- hyper_param_list$learning_rate      # minimum learning rate
  #
  alpha <- hyper_param_list$momentum             # momentum (not used)
  #
  flag_adaptive_step_size <- hyper_param_list$adaptive_step_size # adaptive step size w/ Barzilai_Borwein
  lam_ct <- 0
  lam_vec <- c(0.01, 0.005, 0.0025, 0.001, 0.0005, 0.0002, 0.0001, 0.00005, 0.00002)
  
  # BATCH SIZE FOR MINIBATCH GRADIENT DESCENT (vein of Stochastic grad. descent)
  batch_size_save <- hyper_param_list$batch_size
  
  #### Initialize regularization variables ####
  lam_regularize <- hyper_param_list$lam_regularize
  reg_flag <- hyper_param_list$regularization_flag
  if(reg_flag == 0){
    regularize_flag <- 'none'
  }else if (reg_flag == 1){
    regularize_flag <- 'L1'
  }else if (reg_flag == 2){
    regularize_flag <- 'L2'
  }else{
    regularize_flag <- 'none'
  }
  
  #### START TRAINING ####
  if(testing){
    num_epochs <- 100
  }else{
    num_epochs <- hyper_param_list$max_epochs   # max # of EPOCHS to TRAIN model
  }
  cost_vec <- rep(0, num_epochs)              # initialize vector storage for cost function
  cost_vec_test <- rep(0, num_epochs)         # initialize vector storage for cost function
  cost_vec_train <- rep(0, num_epochs)        # initialize vector storage for cost function
  # 
  for(epoch_iter in 1:num_epochs){
   
    # Pseudo-adaptive step size
    if(mod(epoch_iter, 251) == 0 & epoch_iter >= 200 & flag_adaptive_step_size == 1){
      lam_ct <- lam_ct + 1
      
      # scale learning rate
      min_lam <- 0.90*min_lam
      
      # use stored learning rate values
      #min_lam <- lam_vec[lam_ct]
      
      # Give computer time to rest:
      # Sys.sleep(1)
    }
    lambda_1 <- min_lam
    lambda_2 <- min_lam
    lambda_end <- min_lam
    lambda_b1 <- min_lam
    lambda_b2 <- min_lam
    lambda_bend <- min_lam
    
    # RANDOMLY SHUFFLE INDICES for MINI-BATCH RANDOM SAMPLING and RESET EPOCH PARAMETERS
    if(testing){
      inds_random <- as.matrix(read.csv(paste0("./test/nn_matlab/Matlab_indsRandom_epoch_", 
                                               epoch_iter, ".csv"), header = F))
    }else{
      inds_random <- randperm(length(1:num_train))    # Randomly shuffle training data indices for SGD
    }
    cost_sum <- 0                                   # Reset single epoch cost to 0
    batch_size <- batch_size_save                   # Reset to original batch size
    
    # Iteration number inside single epoch
    for(iter in 1:floor(num_train/batch_size_save)){
      #iter <- 5
      # RANDOMLY SHUFFLE INDICES for MINI-BATCH RANDOM SAMPLING
      if(iter != floor(num_train/batch_size)){
        inds <- inds_random[(1 + (iter - 1)*batch_size):(batch_size*iter)]
      }else{
        inds <- inds_random[(1 + (iter - 1)*batch_size):length(inds_random)]
      }
      
      # Forward propagation 
      forward_results <- forward_propagation(x0[,inds], w1, w2, wend, b1, b2, bend)
      z_hat <- forward_results$z_hat
      x2 <- forward_results$x2
      x1 <- forward_results$x1
      
      # Compute Cost Function (w/ or w/o regularization)
      # Orig. Vectorized Mean-Squared Error (MSE) Cost ("loss") Function --------
      J <- 0.5 * (z0[,inds] - z_hat)^2
      
      # Regularize the Cost Function
      if(regularize_flag == "L2"){
        cost <- (sum(J) + 0.5 * lam_regularize * (sum(w1^2)  + sum(wend^2) + sum(b1^2)))
      }else if(regularize_flag == "L1"){
        cost <- (sum(J) + 0.5 * lam_regularize * (sum(abs(w1)) + sum(abs(wend)) + sum(abs(b1))))
      }else {
        cost <- sum(J)
      }
      cost_sum <- cost_sum + cost       # cumulative cost across epoch
      
      # Compute Delta Matrices (Fall 2022 w/ Eileen Yizzi)
      #    uses f_out = x; so f'_{out}=1
      delta_end <- -(z0[,inds] - z_hat)
      act_prime_val2 <-  w2 %*% x1 + b2
      delta_2 <- (t(wend) %*% delta_end)* activation_function_prime(act_prime_val2)
      act_prime_val1 <- w1 %*% x0[,inds] + b1
      delta_1 <- (t(w2) %*% delta_2) * activation_function_prime(act_prime_val1)
      
      # Compute regularization gradients
      if(regularize_flag == "L2"){
        regularize_grad_w1 <- lam_regularize*w1
        regularize_grad_w2 <- lam_regularize*w2
        regularize_grad_wend <- lam_regularize*wend
        regularize_grad_b1 <- lam_regularize*b1
        regularize_grad_b2 <- lam_regularize*b2
        regularize_grad_bend <- lam_regularize*bend
      }else if (regularize_flag == "L1"){
        regularize_grad_w1 <- lam_regularize*sign(w1)
        regularize_grad_w2 <- lam_regularize*sign(w2)
        regularize_grad_wend <- lam_regularize*sign(wend)
        regularize_grad_b1 <- lam_regularize*sign(b1)
        regularize_grad_b2 <- lam_regularize*sign(b2)
        regularize_grad_bend <- lam_regularize*sign(bend)
      }else{
        regularize_grad_w1 <- 0
        regularize_grad_w2 <- 0
        regularize_grad_wend <- 0
        regularize_grad_b1 <- 0
        regularize_grad_b2 <- 0
        regularize_grad_bend <- 0
      }
      
      # COMPUTE GRADIENTS!  [>>> w/ Eileen Yizzi (Fall 2022) <<<]
      #  Note: 1/num_Z added from COST function coeff. 
      #    (not included in partial deriv. derivation)
      dj_dwend <- 1/batch_size * ((delta_end %*% t(x2)) + regularize_grad_wend)
      dj_dw2 <- 1/batch_size * ((delta_2 %*% t(x1)) + regularize_grad_w2)
      dj_dw1 <- 1/batch_size * ((delta_1 %*% t(x0[,inds])) + regularize_grad_w1)
      dj_dbend <- 1/batch_size * ((delta_end %*% rep(1, batch_size)) + regularize_grad_bend)
      dj_db2 <- 1/batch_size * ((delta_2 %*% rep(1, batch_size)) + regularize_grad_b2)
      dj_db1 <- 1/batch_size * ((delta_1 * rep(1, batch_size)) + regularize_grad_b1)
      
      # Perform Gradient Descent (with momentum <<-- not incorporated)
      # Weights
      w1n <- w1 - lambda_1 * (dj_dw1) # - alpha*lambda_1*dj_dw1_p
      w2n <- w2 - lambda_2 * (dj_dw2) # - alpha*lambda_2*dj_dw2_p
      wendn <- wend - lambda_end * (dj_dwend) # - alpha*lambda_end*dj_dwend_p
      
      # Biases
      b1n <- b1 - lambda_b1 * (dj_db1) # - alpha*lambda_b1 * dj_db1_p
      b2n <- b2 - lambda_b2 * (dj_db2) # - alpha*lambda_b2 * dj_db2_p
      bendn <- bend - lambda_bend * (dj_dbend) # - alpha*lambda_b2 * dj_db2_p
      
      # Use adaptive step-sizes (if flag = 1) 
      #### Not included in R code 
      
      # Increment to next
      w1 <- w1n
      w2 <- w2n
      wend <- wendn
      #
      b1 <- b1n
      b2 <- b2n
      bend <- bendn
    }
    
    # Save COST for SINGLE Epoch (TRAINING)
    cost_vec[epoch_iter] <- cost_sum
    
    # Save COST for SINGLE Epoch (TRAINING DATA)
    forward_train_results = forward_propagation(x0,w1,w2,wend,b1,b2,bend);
    z_hat_train <- forward_train_results$z_hat
    
    j_train <- sum( 0.5*(z0 - z_hat_train )^2 )
    
    cost_vec_train[epoch_iter] <- j_train/num_z
    
    # Save COST for SINGLE Epoch (TESTING DATA)
    forward_test_results = forward_propagation(x0_test,w1,w2,wend,b1,b2,bend);
    z_hat_test <- forward_test_results$z_hat
    
    j_test <- sum( 0.5*(z0_test - z_hat_test)^2 )
    cost_vec_test[epoch_iter] <- j_test/num_z_test
    
    # Compute Differences in Errors / Cost Functions and have it print the info to screen
    if(mod(epoch_iter, hyper_param_list$print_interval) == 0){
      print("#--------------------------------------#")
      print(paste0("   *** Epoch: ", epoch_iter, " ***"))
      print(paste0("Cost (BackP) = ", cost_vec[epoch_iter]))
      print(paste0("Cost (Train) = ", j_train))
      print(paste0("Cost (Test) = ", j_test))
    }
    
  }
  
  print("Model Training done!!!")
  
  if(testing){
    print("Testing R training data against MATLAB after 100 epochs.")
    test_results <- rep(NA, 7)
    Matlab_w1 <- as.matrix(read.csv("./test/nn_matlab/Matlab_w1_epoch_100.csv", header=F))
    test_results[1] <- all.equal(Matlab_w1, w1)
    Matlab_w2 <- as.matrix(read.csv("./test/nn_matlab/Matlab_w2_epoch_100.csv", header=F))
    test_results[2] <- all.equal(Matlab_w2, w2)
    Matlab_wend <- as.matrix(read.csv("./test/nn_matlab/Matlab_wend_epoch_100.csv", header=F))
    test_results[3] <- all.equal(Matlab_wend, wend)
    
    Matlab_b1 <- as.matrix(read.csv("./test/nn_matlab/Matlab_b1_epoch_100.csv", header=F))
    test_results[4] <- all.equal(Matlab_b1, b1)
    Matlab_b2 <- as.matrix(read.csv("./test/nn_matlab/Matlab_b2_epoch_100.csv", header=F))
    test_results[5] <- all.equal(Matlab_b2, b2)
    Matlab_bend <- as.matrix(read.csv("./test/nn_matlab/Matlab_bend_epoch_100.csv", header=F))
    test_results[6] <- all.equal(Matlab_bend, bend)
    
    Matlab_costvec <- read.csv("./test/nn_matlab/Matlab_costvec.csv", header = F)
    test_results[7] <- all.equal(as.vector(t(Matlab_costvec[1,])), cost_vec)
    if(all(test_results)){
      print("All tests report TRUE.")
    }
  }
  
  completed_training <- list(
    "w1_save" = w1,
    "w2_save" = w2,
    "wend_save" = wend,
    "b1_save" = b1,
    "b2_save" = b2,
    "bend_save" = bend,
    "cost_vec" = cost_vec,
    "cost_vec_test" = cost_vec_test,
    "cost_vec_train" = cost_vec_train,
    "min_cost" = cost)
  
  return(completed_training)
}

forward_propagation <- function(x0_0, w1, w2, wend, b1, b2, bend){
  # Compute z1 = w1*x0 (1st hidden layer)
  x1 <- w1 %*% x0_0 + as.vector(b1)
  
  # Apply activation function to x1 (redefining x1 in the process)
  x1 <- activation_function(x1)
  
  # Compute z2 = w2*x1
  x2 <- w2%*%x1 + as.vector(b2)
  
  # Apply activation function to x2 (redefining x2 in the process)
  x2 <- activation_function(x2)
  
  # Compute zHat = WEnd*x2 + bEnd (here, activation is thought to be sigma(x)=x)
  z_hat <- wend%*%x2 + as.vector(bend)
  
  return(list("z_hat" = z_hat, "x2" = x2, "x1" = x1))
}

activation_function <- function(x1, function_type = "ReLU"){
  if(function_type == "Sigmoid"){
    A <- 1/(1 + exp(-x1))
  }else if(function_type == "ReLU"){
    A <- matlab_max(x1, 0)
  }else {
    stop("Unknown activation function type.")
  }
  return(A)
}

activation_function_prime <- function(z, function_type = "ReLU"){
  if(function_type == "Sigmoid"){
    A <- exp(-z)/((1 + exp(-z))^2)
  }else if(function_type == "ReLU"){
    A <- matrix(rep(0, prod(dim(z))), nrow = dim(z)[1], ncol = dim(z)[2])
    A[which(z > 0)] <- 1
  }else{
    stop("Unknown activation function prime type")
  }
  return(A)
}

matlab_max <- function(x, a){
  output <- matrix(NA, nrow = dim(x)[1], ncol = dim(x)[2])
  for(i in 1:ncol(x)){
    for(j in 1:nrow(x)){
      if(x[j,i] > a){
        output[j, i] <- x[j, i]
      }else{
        output[j, i] <- a
      }
    }
  }
  return(output)
}

save_trained_model_values <- function(trained_model){
  dir.create("./src/r-scripts/nn-data/", showWarnings = F)
  for(i in names(trained_model)){
    write.table(trained_model[[i]], 
                file = paste0("./src/r-scripts/nn-data/",i,".csv"),
                row.names = F, col.names = F, sep = ",")
  }
  return(NULL)
}

print_error_information <- function(zdat, dat, report_option, num_inputs){
  # compute errors
  err_abs <- abs(zdat - dat[, (num_inputs + 1)])
  err_rel <- (err_abs / dat[, (num_inputs + 1)]) * 100
  inds <- is.infinite(err_rel)
  err_rel[inds] <- NA
  
  # print error
  print(paste0(" *** ", report_option," ERROR INFORMATION ***"))
  print(" ")
  print("Absolute Error: ")
  print(paste0("Max ABS error: ", max(err_abs, na.rm = T)))
  print(paste0("Min ABS error: ", min(err_abs, na.rm = T)))
  print(paste0("Avg ABS error: ", mean(err_abs, na.rm = T)))
  print(" ")
  #
  print("Relative Error: ")
  print(paste0("Max ABS error: ", max(err_rel, na.rm = T)))
  print(paste0("Min ABS error: ", min(err_rel, na.rm = T)))
  print(paste0("Avg ABS error: ", mean(err_rel, na.rm = T)))
  print(" ")
}

create_2D <- function(x_mesh, y_mesh, z_const, trained_model, 
                      true_dat = FALSE, train_test_output = NULL, 
                      eps_error = NULL){
  f_prediction <- matrix(NA, nrow = dim(x_mesh)[1], ncol = dim(x_mesh)[2])
  f_real <- matrix(NA, nrow = dim(x_mesh)[1], ncol = dim(x_mesh)[2])
  for(n1 in 1:nrow(f_prediction)){
    for(n2 in 1:ncol(f_prediction)){
      input_vec <- matrix(c(x_mesh[n1, n2], y_mesh[n1, n2], z_const), nrow = 3)
      # Get predicted output:
      forward_out <- forward_propagation(input_vec, trained_model$w1_save, 
                                         trained_model$w2_save, trained_model$wend_save,
                                         trained_model$b1_save, trained_model$b2_save,
                                         trained_model$bend_save)
      f_prediction[n1, n2] <- forward_out$z_hat
      
      if(true_dat){
        # Get True data value 
        f_real[n1, n2] <- evaluate_function(t(input_vec), eps_error, 1, train_test_output$min_z,
                                            train_test_output$max_z)
      }
    }
  }
  if(true_dat){
    return(list("x" = x_mesh, "y" = y_mesh, "z" = z_const,
                "f_prediction" = f_prediction, "f_real" = f_real))
  }else{
    return(list("x" = x_mesh, "y" = y_mesh, "z" = z_const,
                "f_prediction" = f_prediction))
  }
  
}

line_fxn <- function(dat){
  m <- 1/(max(dat) - min(dat))
  b <- -m*min(dat)
  return(m * dat + b)
}

scale_params_data <- function(input_params, log = F){
  input_scaled <- matrix(NA, nrow = nrow(input_params), 
                         ncol = ncol(input_params))
  if(log) input_params[,1] <- log10(input_params[,1])
  input_scaled[,1] <- line_fxn(input_params[,1])
  input_scaled[,2] <- line_fxn(input_params[,2])
  input_scaled[,3] <- line_fxn(input_params[,3])
  return(input_scaled)
}

transform_output_data <- function(dat, option = "z-score", lambda = 1){
  if(option == "box-cox"){
    dat_output <- ((data[,4]^lambda) - 1) / (lambda)
    return(list("trans_output" = dat_output, "lambda" = lambda))
  } else if(option == "box-cox-y"){
    dat_output <- ((dat^lambda) - 1) / (lambda * (mean(dat)^(lambda - 1)))
    return(list("trans_output" = dat_output, "lambda" = lambda, 
                "mean_dat" = mean(dat)))
  } else if(option == "z-score"){
    dat_output <- (dat - mean(dat)) / sd(dat)
    return(list("trans_output" = dat_output, "mean_dat" = mean(dat),
                "sd_dat" = sd(dat)))
  } else {
    stop("Unknown transformation option, please pick another")
  }
}

untransform_output_data <- function(dat, option = "z-score", 
                                    lambda = 1, dat_mean = 1, dat_sd = 1){
  if(option == "z-score"){
    output_dat <- dat*dat_sd + dat_mean
  } else if (option == "box-cox"){
    output_dat <- ((dat*lambda)  + 1)^(1/lambda) 
  } else if (option == "box-cox-y"){
    output_dat <- ((dat*lambda*(dat_mean^(lambda-1)))  + 1)^(1/lambda)
  } else {
    stop("Unknown transformation option, please pick another")
  }
  return(output_dat)
}

scale_output_data <- function(train_output, test_output){
  min_z <- min(c(train_output, test_output))
  max_z <- max(c(train_output, test_output))
  m <- 1/(max_z - min_z)
  b <- -m*min_z
  train_output <- m * train_output + b
  test_output <- m * test_output + b
  return(list("train_output" = train_output, 
              "test_output" = test_output,
              "min_z" = min_z, 
              "max_z" = max_z))
}

unscale_output_data <- function(output_data, min_z, max_z){
  m <- 1/(max_z - min_z)
  b <- -m*min_z
  unscaled_output <- (output_data - b)/m
  return(unscaled_output)
}
