
library(Matrix)
library(hts)
library(parallel)
library(bigstatsr)


############ immute function ##################

iMinT <- function(fcasts, weights, immute, smat, solver = c("quadprog", "osqp"), 
                    nonnegative = FALSE, parallel = FALSE, settings = NULL)
          { 
            solver <- match.arg(solver)
            nts <- nrow(smat)
            nbts <- ncol(smat)
            nagg <- nts - nbts
            seqagg <- 1L:nagg
            utmat <- cbind2(Matrix::sparseMatrix(i = seqagg, j = seqagg, x = 1), 
                            -1*smat[1L:nagg, ])
            jmat <- Matrix::sparseMatrix(i = 1L:nbts, j = (nagg + 1L):nts, 
                                 x = rep(1L, nbts), dims = c(nbts, nts))
            if (is.null(immute)) {
              immute <- 1
              message("Top level is kept immutable as immutable levels were not specified.")
            }
            if (length(immute) == 1) {
              smat_l <- Matrix(smat[immute,], ncol = dim(smat)[2], sparse = TRUE)
            }
            else {
              smat_l <- smat[immute,]  #rows of summing matrix corresponding to immutable series
            }
            if (length(immute) > nts) {
              stop("Number of immutable series should not be more than total 
                   number of time series.",
                   call. = FALSE)
            }
            if (length(immute) > nbts) {
              stop("Number of immutable series should not be more than 
                   total number of bottom-level series.",
                   call. = FALSE)
            }
            if (length(immute) != base::qr(smat_l)$rank) {
              stop("All specified immutable series cannot be immutable at the same time. 
                   Some nodes need to be free.",
                   call. = FALSE)
            }
            
            reconciledf <- matrix(0, dim(fcasts)[1], dim(fcasts)[2])
            
            # Matrix and vectors in quadratic function
            Rmat <- Matrix::chol(weights) # upper triangular part of weights matrix (R'R=weights)
            matrix1 <- kronecker(Rmat %*% t(utmat), smat) # matrix infront of vector x in the quadratic functions
            dmat <- t(matrix1) %*% matrix1
            dvec <- t(as.vector(smat %*% jmat %*% t(Rmat))) %*% matrix1
            
            reconcile.forecast <- function(yhat) {
              yhat_l <- yhat[immute]
              if (!nonnegative) {
                amat <- kronecker(t(yhat) %*% t(utmat), smat_l)
                bvec <- yhat_l - (smat_l %*% jmat %*% yhat) 
              } 
              else {
                amat <- rbind(kronecker(t(yhat) %*% t(utmat), smat_l),
                              kronecker(t(yhat) %*% t(utmat), 
                                        Matrix::sparseMatrix(i = 1:dim(smat_l)[2], 
                                                     j = 1:dim(smat_l)[2], 
                                                     x = 1, 
                                                     dims = c(dim(smat_l)[2], dim(smat_l)[2])))
                )
                bvec <- rbind(yhat_l - (smat_l %*% jmat %*% yhat),
                              -jmat %*% yhat)
              }
              if (solver == "quadprog") {
                qp_result <- quadprog::solve.QP(Dmat = dmat, dvec = -dvec, Amat = t(amat),
                                                bvec = bvec, meq = length(immute))
                xmat <- Matrix(qp_result$solution, nrow = dim(smat)[2])
              }
              else if (solver == "osqp") {
                if (is.null(settings)){
                  settings = osqp::osqpSettings(verbose = FALSE,
                                          eps_abs = 1e-8,
                                          eps_rel = 1e-8,
                                          polish_refine_iter = 100,
                                          polish = TRUE)
                }
                if (nonnegative) {
                  uvec <- rbind(yhat_l - (smat_l %*% jmat %*% yhat),
                                Matrix(rep(Inf, dim(smat)[2])))
                }
                else {
                  uvec <- bvec
                }
                qp_result <- osqp::solve_osqp(P = dmat, q = t(dvec), A = amat, 
                                              l=bvec, u=uvec, pars = settings)
                if (qp_result$info$status != "solved") {
                  cat("The problem has not converged to an optimal solution.\n")
                } 
                xmat <- Matrix(qp_result$x, nrow=dim(smat)[2])
              }
              gmat <- jmat + xmat %*% utmat
              as.vector(t(smat %*% gmat %*% yhat))
            }
            
            if (parallel) {
              cl <- parallel::makeCluster(parallel::detectCores() - 1)
              parallel::clusterExport(cl, varlist = c("fcasts", "weights", "immute", "smat", "nonnegative", 
                                            "settings", "smat_l", "utmat", "jmat", "dmat", "dvec", 
                                            "reconcile.forecast"), envir = environment())
              parallel::clusterEvalQ(cl, {
                library(Matrix)
                library(quadprog)
              })

              reconciledf <- t(parallel::parSapply(cl, 1:nrow(fcasts), function(i) reconcile.forecast(fcasts[i, ])))

              parallel::stopCluster(cl)
            }
            else {
              reconciledf <- t(sapply(1:nrow(fcasts), function(i) reconcile.forecast(fcasts[i, ])))
            }
            
            return(reconciledf)
          }


############ combinef.iMinT function ##################

combinef.iMinT <- function(fcasts, nodes = NULL, groups = NULL, variances = NULL, 
                            immute = NULL, nonnegative = FALSE, 
                            solver = c("quadprog", "osqp"), 
                            settings = NULL,
                            keep = c("all", "gts", "bottom"),
                            parallel = FALSE
                            )
              {
                if (is.null(nodes) && is.null(groups)) {
                  stop("Please specify the hierarchical or the grouping structure.", 
                       call. = FALSE)
                }
                if (!xor(is.null(nodes), is.null(groups))) {
                  stop("Please specify either nodes or groups argument, not both.", 
                       call. = FALSE)
                }
                solver <- match.arg(solver)
                keep <- match.arg(keep)
                fcasts <- stats::as.ts(fcasts)
                tspx <- stats::tsp(fcasts)
                cnames <- colnames(fcasts)
                if (nonnegative) {
                  if (any(fcasts < 0)) {
                    fcasts[fcasts < 0] <- 0
                    warning("Negative base forecasts are truncated to zero.")
                  }
                }
                if (is.null(groups)) {
                  totalts <- sum(hts:::Mnodes(nodes))
                  if (!is.matrix(fcasts)) {
                    fcasts <- t(fcasts)
                  }
                  h <- nrow(fcasts)
                  if (ncol(fcasts) != totalts) {
                    stop("Argument fcasts requires all the forecasts.",
                         call. = FALSE)
                  }
                  gmat <- hts:::GmatrixH(nodes)
                }
                else if (is.null(nodes)) {
                  rownames(groups) <- NULL
                  gmat <- hts:::GmatrixG(groups)
                  totalts <- sum(hts:::Mlevel(gmat))
                  if (!is.matrix(fcasts)) {
                    fcasts <- t(fcasts)
                  }
                  if (ncol(fcasts) != totalts) {
                    stop("Argument fcasts requires all the forecasts.", 
                         call. = FALSE)
                  }
                }
                smat <- hts:::SmatrixM(gmat)
                seqts <- 1:totalts
                if (is.null(variances)) {
                  weights <- Matrix::sparseMatrix(i=seqts, j=seqts, x=1)
                }
                else {
                  weights <- Matrix::sparseMatrix(i=seqts, j=seqts, x=variances)
                }
                reconciledf <- iMinT(fcasts=fcasts, weights=weights, 
                                       immute=immute, smat=smat, solver=solver,
                                       nonnegative=nonnegative, parallel=parallel,
                                       settings=settings)
                if (keep == "all") {
                  out <- reconciledf
                }
                else {
                  bottom <- totalts - (ncol(smat):1L) + 1L
                  bf <- reconciledf[, bottom]
                  colnames(bf) <- cnames[bottom]
                  if (keep == "gts") {
                    bf <- ts(bf, start = tspx[1L], frequency = tspx[3L])
                    if (is.null(groups)) {
                      out <- suppressMessages(hts(bf, nodes = nodes))
                    }
                    else if (is.null(nodes)) {
                      out <- suppressMessages(gts(bf, groups = groups))
                    }
                  }
                  else {
                    out <- bf
                  }
                }
              return(out)
            }



############ MinT.immute function ##################

MinT.iMinT <- function(fcasts, nodes = NULL, groups = NULL, residual, 
                        covariance = c("shr", "sam"),
                        immute = 1, nonnegative = FALSE, 
                        solver = c("quadprog", "osqp"), 
                        settings = NULL,
                        keep = c("all", "gts", "bottom"),
                        parallel = FALSE
                        )
          {
            if (is.null(nodes) && is.null(groups)) {
              stop("Please specify the hierarchical or the grouping structure.", 
                   call. = FALSE)
            }
            if (!xor(is.null(nodes), is.null(groups))) {
              stop("Please specify either nodes or groups argument, not both.", 
                   call. = FALSE)
            }
            solver <- match.arg(solver)
            keep <- match.arg(keep)
            covar <- match.arg(covariance)
            res <- residual
            fcasts <- stats::as.ts(fcasts)
            tspx <- stats::tsp(fcasts)
            cnames <- colnames(fcasts)
            if (nonnegative) {
              if (any(fcasts < 0)) {
                fcasts[fcasts < 0] <- 0
                warning("Negative base forecasts are truncated to zero.")
              }
            }
            if (missing(residual)) {
              stop("MinT needs insample residuals.", call. = FALSE)
            }
            if (covar == "sam") {
              n <- nrow(res)
              w.1 <- crossprod(res)/n
              if (hts:::is.posdef(w.1) == FALSE) {
                stop("MinT needs covariance matrix to be positive definite.", 
                     call. = FALSE)
              }
            }
            else {
              tar <- hts:::lowerD(res)
              shrink <- hts:::shrink.estim(res, tar)
              w.1 <- shrink[[1]]
              lambda <- shrink[[2]]
              if (hts:::is.posdef(w.1) == FALSE) {
                stop("MinT needs covariance matrix to be positive definite.", 
                     call. = FALSE)
              }
            }
            if (!is.null(w.1)) {
              weights <- methods::as(w.1, "sparseMatrix") # do we put 1/w.1 here as well? or doesn't matter?
            }
            if (is.null(groups)) {
              totalts <- sum(hts:::Mnodes(nodes))
              if (!is.matrix(fcasts)) {
                fcasts <- t(fcasts)
              }
              h <- nrow(fcasts)
              if (ncol(fcasts) != totalts) {
                stop("Argument fcasts requires all the forecasts.",
                     call. = FALSE)
              }
              gmat <- hts:::GmatrixH(nodes)
            }
            else if (is.null(nodes)) {
              rownames(groups) <- NULL
              gmat <- hts:::GmatrixG(groups)
              totalts <- sum(hts:::Mlevel(gmat))
              if (!is.matrix(fcasts)) {
                fcasts <- t(fcasts)
              }
              if (ncol(fcasts) != totalts) {
                stop("Argument fcasts requires all the forecasts.", 
                     call. = FALSE)
              }
            }
            smat <- hts:::SmatrixM(gmat)
            reconciledf <- iMinT(fcasts=fcasts, weights=weights, 
                                   immute=immute, smat=smat, solver=solver,
                                   nonnegative=nonnegative, parallel=parallel,
                                   settings=settings)
            if (keep == "all") {
              out <- reconciledf
            }
            else {
              bottom <- totalts - (ncol(smat):1L) + 1L
              bf <- reconciledf[, bottom]
              colnames(bf) <- cnames[bottom]
              if (keep == "gts") {
                bf <- ts(bf, start = tspx[1L], frequency = tspx[3L])
                if (is.null(groups)) {
                  out <- suppressMessages(hts(bf, nodes = nodes))
                }
                else if (is.null(nodes)) {
                  out <- suppressMessages(gts(bf, groups = groups))
                }
              }
              else {
                out <- bf
              }
            } 
          return(out)
        }

