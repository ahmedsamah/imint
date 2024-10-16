
# produce summing matrix of new basis time series
transform.sMat <- function(sMat, basis_set){
  m <- dim(sMat)[2]
  if (length(basis_set) != m){
    stop(simpleError(sprintf('length of basis set should be %d', m)))
  }
  S1 <- sMat[basis_set,]
  S2 <- sMat[-basis_set,]
  transitionMat <- solve(S1, diag(rep(1, m)))
  rbind(S2 %*% transitionMat, diag(rep(1, m)))
}

# forecast.reconcile <- function(base_forecasts, 
#                                sMat,
#                                weighting_matrix,
#                                immu_set=NULL,
#                                nonnegative=FALSE){
#   m <- dim(sMat)[2]
#   n <- dim(sMat)[1]
#   k <- length(immu_set)
#   if (length(immu_set) == 0){
#     weighting_matrix = solve(weighting_matrix)
#     reconciled_y = sMat %*% solve(t(sMat) %*% weighting_matrix %*% sMat) %*% t(sMat) %*% weighting_matrix %*% t(base_forecasts)
#     return(t(reconciled_y))
#   }
# 
#   # construct new basis time series
#   if (k > m) {
#     stop(simpleError(sprintf('length of basis set can not be bigger than %d', m)))
#   }
#   ## select mutable series
#   immutable_basis <- sort(immu_set)
#   candidate_basis <- setdiff((n-m+1):n, immu_set)
#   if (all(immutable_basis >= n-m+1 )){
#     mutable_basis <- candidate_basis
#   } else {
#     mutable_basis <- c()  
#     determined <- c()
#     i <- max(which(immutable_basis < n-m+1))
#     while (length(mutable_basis) != m-k) {
#       corresponding_leaves <- which(sMat[immutable_basis[i], ] != 0) + n - m
#       free_leaves <- setdiff(corresponding_leaves, c(immutable_basis, mutable_basis, determined))
#       if (length(free_leaves) == 0) stop(simpleError('the immu_set can not be used to describe the hierarchy'))
#       if (length(free_leaves) == 1) {
#         candidate_basis <- candidate_basis[candidate_basis != free_leaves[1]]
#       } else{
#         determined <- c(determined, free_leaves[1])
#         mutable_basis <- c(mutable_basis, free_leaves[2:length(free_leaves)])
#         candidate_basis <- candidate_basis[!(candidate_basis %in% free_leaves)]
#       }
#       i <- i - 1
#       if (i == 0) {
#         mutable_basis <- c(mutable_basis, candidate_basis)
#       }
#     }
#   }
#   new_basis <- c(sort(mutable_basis), immutable_basis)
#   sMat <- transform.sMat(sMat, new_basis)
#   S1 <- sMat[1:(n-k),,drop=FALSE][,1:(m-k),drop=FALSE]
#   S2 <- sMat[1:(n-m),,drop=FALSE][,(m-k+1):m,drop=FALSE]
#   determined <- setdiff(1:n, new_basis)
#   mutable_series <- c(determined, mutable_basis)
#   #mutable_weight <- solve(weighting_matrix[mutable_series,,drop=FALSE][,mutable_series,drop=FALSE])
#   weighting_matrix = solve(weighting_matrix)
#   mutable_weight <- weighting_matrix[mutable_series,,drop=FALSE][,mutable_series,drop=FALSE]
#   mutable_base <- cbind(base_forecasts[,determined,drop=FALSE] - t(S2 %*% t(base_forecasts[,immutable_basis,drop=FALSE])),
#                         base_forecasts[,sort(mutable_basis),drop=FALSE])
#   reconciled_mutable <- solve(t(S1) %*% mutable_weight %*% S1) %*% t(S1) %*% mutable_weight %*% t(mutable_base)
#   reconciled_y <- t(sMat %*% rbind(reconciled_mutable, t(base_forecasts[,immutable_basis,drop=FALSE])))
#   mutable_weight <- mutable_weight / max(diag(mutable_weight))
#   if (nonnegative){
#     for (i in 1:dim(mutable_base)[1]){
#       Dmat <- t(S1) %*% mutable_weight %*% S1
#       dvec <- as.vector(t(mutable_base[i,]) %*% mutable_weight %*% S1)
#       Amat <- diag(rep(1, dim(S1)[2]))
#       bvec <- rep(0, dim(S1)[2])
#       sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec)$solution)
#       if (is(sol, "try-error")){
#         warning(paste0("unsolvable at row ", rownames(basef)[i], " use unconstrained solution!"))
#       }else{
#         reconciled_y[i,] <- as.vector(sMat %*% c(sol, base_forecasts[i,immutable_basis,drop=FALSE]))
#       }
#     }
#   }
#   new_index <- c(determined, new_basis)
#   reconciled_y[,order(new_index)]
# }


forecast.reconcile <- function(base_forecasts, 
                               sMat,
                               weighting_matrix,
                               immu_set=NULL,
                               nonnegative=FALSE){
  m <- dim(sMat)[2]
  n <- dim(sMat)[1]
  k <- length(immu_set)
  if (length(immu_set) == 0){
    weighting_matrix = solve(weighting_matrix)
    solution <- matrix(0, dim(base_forecasts)[1], m)
    if (nonnegative){
      for (i in 1:dim(base_forecasts)[1]){
        Dmat <- t(sMat) %*% weighting_matrix %*% sMat
        dvec <- t(sMat) %*% weighting_matrix %*% base_forecasts[i,]
        Amat <- diag(rep(1, m))
        solution[i,] <- quadprog::solve.QP(Dmat, dvec, Amat)$solution
      }
      reconciled_y <- sMat %*% t(solution)
    }else{
      reconciled_y <- sMat %*% solve(t(sMat) %*% weighting_matrix %*% sMat) %*% t(sMat) %*% weighting_matrix %*% t(base_forecasts)
    }
    return(t(reconciled_y))
  }
  
  # construct new basis time series
  if (k > m) {
    stop(simpleError(sprintf('length of basis set can not be bigger than %d', m)))
  }
  ## select mutable series
  immutable_basis <- sort(immu_set)
  candidate_basis <- setdiff((n-m+1):n, immu_set)
  determined <- 1:(n-m)
  mutable_basis <- candidate_basis
  if (any(immutable_basis < n-m+1)){
    i <- max(which(immutable_basis < n-m+1))
    while (i > 0) {
      corresponding_leaves <- which(sMat[immutable_basis[i], ] != 0) + n - m
      free_leaves <- setdiff(corresponding_leaves, c(immutable_basis, determined))
      if (length(free_leaves) == 0){
        if (all(corresponding_leaves %in% immutable_basis)){
          warning(paste0("all children of ", immutable_basis[i], "th series are immutable, it is removed from the condition."))
          k <- k - 1
          immutable_basis <- immutable_basis[immutable_basis != immutable_basis[i]]
          i <- i - 1
          next
        } else{
          stop(simpleError("can not describe the hierarchy."))
        }
      }
      determined <- determined[determined != immutable_basis[i]]
      determined <- c(determined, free_leaves[1])
      mutable_basis <- mutable_basis[mutable_basis != free_leaves[1]]
      i <- i - 1
    }
  }
  new_basis <- c(sort(mutable_basis), immutable_basis)
  sMat <- transform.sMat(sMat, new_basis)
  S1 <- sMat[1:(n-k),,drop=FALSE][,1:(m-k),drop=FALSE]
  S2 <- sMat[1:(n-m),,drop=FALSE][,(m-k+1):m,drop=FALSE]
  determined <- setdiff(1:n, new_basis)
  mutable_series <- c(determined, mutable_basis)
  # mutable_weight <- solve(weighting_matrix)[mutable_series,,drop=FALSE][,mutable_series,drop=FALSE]
  weighting_matrix = solve(weighting_matrix)
  mutable_weight <- weighting_matrix[mutable_series,,drop=FALSE][,mutable_series,drop=FALSE]
  mutable_base <- cbind(base_forecasts[,determined,drop=FALSE] - t(S2 %*% t(base_forecasts[,immutable_basis,drop=FALSE])),
                        base_forecasts[,sort(mutable_basis),drop=FALSE])
  reconciled_mutable <- solve(t(S1) %*% mutable_weight %*% S1) %*% t(S1) %*% mutable_weight %*% t(mutable_base)
  reconciled_y <- t(sMat %*% rbind(reconciled_mutable, t(base_forecasts[,immutable_basis,drop=FALSE])))
  mutable_weight <- mutable_weight / max(diag(mutable_weight))
  if (nonnegative){
    for (i in 1:dim(mutable_base)[1]){
      Dmat <- t(S1) %*% mutable_weight %*% S1
      dvec <- as.vector(t(mutable_base[i,]) %*% mutable_weight %*% S1)
      Amat <- diag(rep(1, dim(S1)[2]))
      bvec <- rep(0, dim(S1)[2])
      sol <- try(quadprog::solve.QP(Dmat, dvec, Amat, bvec)$solution)
      if (is(sol, "try-error")){
        warning(paste0("unsolvable at row ", rownames(basef)[i], " use unconstrained solution!"))
      }else{
        reconciled_y[i,] <- as.vector(sMat %*% c(sol, base_forecasts[i,immutable_basis,drop=FALSE]))
      }
    }
  }
  new_index <- c(determined, new_basis)
  reconciled_y[,order(new_index)]
}

