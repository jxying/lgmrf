#' @title Learn sparse graph Lapalcian in Laplacian constrained Gaussian Markov Random Field (LGMRF)
#'
#' @param S sample covariance matrix, where p is the number of nodes of the graph
#' @param lmbd, eps and q hyperparameters to control the sparsity regularization by 1/( |Omega_ij| + eps)^q
#' @param backtrack whether to update the learning rate using backtrack line search
#' @param recursive whether to update the regularization weight using the estimate 
#'        obtained in the last outer loop
#' @param maxiter maximum number of iterations in the inner loop
#' @param reltol_out relative tolerance on the Frobenius norm of the estimated
#'        Laplacian matrix as a stopping criteria in the outer loop
#' @param reltol_inner relative tolerance on the Frobenius norm of the estimated
#'        Laplacian matrix as a stopping criteria in the inner loop
#' @param verbose whether or not to show a progress bar displaying the
#'        iterations
#' @return A list containing possibly the following elements:
#' \item{\code{laplacian}}{the estimated graph Lapalcian}
#' \item{\code{adjacency}}{the estimated Adjacency Matrix}
#' \item{\code{convergence}}{boolean flag to indicate whether or not the optimization converged}
#' \item{\code{elapsed_time}}{elapsed time recorded at every iteration}


learn_laplacian_alpe <- function(S, lmbd = 0.1, eps = 1e-4, q = 1, backtrack = TRUE, recursive = TRUE,
                    maxiter = 1000, reltol_out = 1e-7, reltol_inner = 1e-7, verbose = FALSE) {
  # number of nodes
  p <- nrow(S)

  # initial estimate  
   Sinv <- MASS::ginv(S)    
   w0 <- spectralGraphTopology:::Linv(Sinv)
   w0[w0 < 0] <- 0
   w <- w0 + 1e-4
    
   J <- matrix(1, p, p) / p
   H <- lmbd * (diag(p) - p * J) 
    
  if (recursive) {
      maxiter_out = 5
  }else{
      maxiter_out = 2
  }
  # set initial step size
  eta <- 0.01
  if (verbose)
    pb <- progress::progress_bar$new(format = "<:bar> :current/:total ",
                                     total = maxiter_out, clear = FALSE, width = 80)
  time_seq_1 <- c(0)
  start_time_1 <- proc.time()[3]
  w_prev <- matrix(0, p*(p-1)/2, 1)

  # Outer iteration 
  for (i in 1:maxiter_out) {
      Lw_pre <- L(w_prev)
      
      if (i > 1){
                    
     K <- S + H / ((-Lw_pre + eps)^q) 
      }else{
          K <- S
      }
      
    # Inner iteration
    for (j in 1:maxiter) {
        Lw <- L(w)
    tryCatch(
        
        {gradient <- spectralGraphTopology:::Lstar(K - spectralGraphTopology:::inv_sympd(Lw + J))},
        
        error = function(err) {
        results <- list(laplacian = L(w_prev), maxiter_out = i,
                     convergence = FALSE, elapsed_time = time_seq_1)
        return(results)
              }
              )    
    if (backtrack) {
      fun_list <- mle_pgd.obj(Lw = Lw, J = J, K = K)
      fun <- fun_list$obj_func
      while(1) {

        wi <- w - eta * gradient
        wi[wi < 0] <- 0
        Lwi <- L(wi)  
        # compute the objective function at the updated value of w
        fun_t_list <- mle_pgd.obj(Lw = Lwi, J = J, K = K)
        fun_t <- fun_t_list$obj_func
        is_NPD <- fun_t_list$is_NPD
        # check whether the previous value of the objective function is
        # smaller than the current one
        conditions <- (fun < fun_t - sum(gradient * (wi - w)) - (.5/eta)*norm(wi - w, '2')^2) | is_NPD
        if (conditions) {
          eta <- .5 * eta
            if (eta < 1e-15)
                {
                break
            }
                
        } else {
          eta <- 2 * eta
          break
        }
      }
    } else {
        wi <- w - eta * gradient
        wi[wi < 0] <- 0
        Lwi <- L(wi)
    }
    has_converged_inner <- (norm(Lwi - Lw, 'F') / norm(Lw, 'F') < reltol_inner) && (j > 10)
    if (has_converged_inner)
      break
    w <- wi   
    }

      
    if (verbose)
        pb$tick() 
         
    Lwi <- L(wi)  
    has_converged <- (norm(Lw_pre - Lwi, 'F') / norm(Lw_pre, 'F') < reltol_out)
    time_seq_1 <- c(time_seq_1, proc.time()[3] - start_time_1)
    if (has_converged)
         break
    w_prev <- wi
  
  }
 
  results <- list(laplacian = Lwi, adjacency = A(wi), convergence = has_converged, elapsed_time = time_seq_1)
  return(results)
}


   mle_pgd.obj <- function(Lw, J, K) {
   cholStatus <- try(chol_factor <- chol(Lw + J), silent = TRUE)
   cholError <- ifelse(class(cholStatus) == "try-error", TRUE, FALSE) 
   if (cholError[1]){      
       return_list <- list(obj_func = 1e16, is_NPD = TRUE)
       return(return_list)
   }else{
        return_list <- list(obj_func = sum(Lw*K) - 2*sum(log(diag(chol_factor))), is_NPD = FALSE)
        return(return_list)
   }

        
}
    
    

