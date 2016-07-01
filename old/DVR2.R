# DVR2 uses ecdf to create pseudo-copula data in nonparametric case, backwards compatible

# Performs D-vine Regression on given data
# Author: Daniel Kraus

# Input:
# y: response data (nx1-vector).
# X: predictor data (nxd-matrix).
# margins.known: if FALSE then margins will be estimated (1xd-vector of logicals).
# margins: how to estimate unknown margins; or classes of known margins (1xd-vector of strings).
# margin.pars: list of known margin parameters (1xd-list), unknown parameters are left empty.
# selectioncrit: penalization of the conditional loglikelihood of the model. Note that the pair-copulas are estimated using the BIC-criterion.
# use_pobs: locical, whether to use the function pobs to produce copula-data.
# correctionterm: logical, whether the cll term should be (AIC- or BIC-)corrected.
# bw.sel: how to estimate the bandwidth for nonparametric estimation of marginals.
# progress: logical, whether progress should be printed during calculation

DVR = function(y, X, margins.known=FALSE, margins="nonpar", margin.pars=list(),
               selectioncrit="BIC", familyset=NA, indeptestlevel=0.05, use_pobs=FALSE, correctionterm=TRUE,
               bw.sel="density", progress=FALSE){
  X<-as.matrix(X)
  n = length(y)
  d = ncol(X)
  p = d+1

  if(length(margins.known)==1) margins.known = rep(margins.known,p)
  if(length(margins)==1) margins = rep(margins,p)
  xdata = cbind(y,X)
  udata = matrix(0,n,p)
  for (j in (1:p)[margins.known==FALSE]){
    margin.pars[[j]] = estim.margins(xdata[,j],margins[j],bw.sel)
  }
  pobsind = (margins=="nonpar"|margins=="kcde")&use_pobs==T
  for(j in (1:p)[!pobsind]){
    udata[,j] = CDF(xdata[,j],margins[j],margin.pars[[j]])
  }
  udata[,pobsind] = pobs(xdata[,pobsind,drop=F])
  v = udata[,1]
  U = as.matrix(udata[,2:p])

  # remaining.variables = 1:d
  # my.index = NULL
  # global.max.ll = -Inf

  # initialize 1-dimensional DVine
  # currentDvine = list(RVM=list(Matrix=matrix(1,1,1),
  #                              family=matrix(0,1,1),
  #                              par=matrix(0,1,1),
  #                              par2=matrix(0,1,1)),
  #                     V=list(direct=array(v, dim = c(1, 1, n)),
  #                            indirect=array(NA, dim = c(1, 1, n)))
  # )
  #
  # for (i in 1:d){
  #   # check which variable update increases the loglikelihood of the conditional density f_V|U_I the most
  #   cll = vector(length=length(remaining.variables))
  #   DvineUpdate = list()
  #   for (j in 1:length(remaining.variables)){
  #     # update current D-Vine by adding j-th remaining variable
  #     DvineUpdate[[j]] = updateDvine(newcolumn=U[,remaining.variables[j]], familyset = familyset,
  #                                    selectioncrit = "BIC", indeptest = TRUE, level = indeptestlevel,
  #                                    currentDvine=currentDvine)
  #     # number of vine's parameters for cll calculation
  #     npar = sum(DvineUpdate[[j]]$RVM$par!=0)+sum(DvineUpdate[[j]]$RVM$par2!=0)
  #     cll[j] = sum(RVineLogLik(cbind(v,U[,my.index[1:(i-1)]],U[,remaining.variables[j]]), DvineUpdate[[j]]$RVM)$V$value[,1])+
  #       (correctionterm==T)*(-(selectioncrit=="AIC")*npar-(selectioncrit=="BIC")*npar*log(n)/2)
  #   }
  #   # pick the one with the highest cll. If none of the remaining variables increases the overall cll break the loop
  #   maxInd <- which.max(cll)
  #   if (cll[maxInd]<=global.max.ll)break
  #   my.index = c(my.index,remaining.variables[maxInd])
  #   global.max.ll = max(global.max.ll,cll[maxInd])
  #   currentDvine = DvineUpdate[[maxInd]]
  #   remaining.variables = setdiff(remaining.variables,my.index[i])
  #   if (progress) print(my.index)
  # }
  # # return the best Dvine
  # reorder = my.index
  # reorder[order(reorder)] = 1:length(my.index)
  # currentDvine$RVM$Matrix = DVineMatGen(elements=c(1,reorder+1))
  # used = sort(my.index)

  vine.fit <- RVineCopSelect(cbind(v, U), familyset = familyset, Matrix = DVineMatGen(p))
  DVR.object = list(DVM=vine.fit,margins=margins,margin.pars=margin.pars,
                    my.index=1:d,used=1:d)
  return(DVR.object)
}

# updateDvine - estimates missing pair-copulas to update the current D-vine by adding one leaf at the end.
updateDvine <- function (newcolumn, familyset = NA, selectioncrit = "BIC",
                         indeptest = TRUE, level = 0.05, trunclevel = NA, rotations = TRUE,
                         currentDvine){
  d <- n <- ncol(currentDvine$RVM$Matrix)+1
  N <- length(newcolumn)
  if (is.na(trunclevel))
    trunclevel <- n
  types <- familyset
  if (trunclevel == 0)
    types <- 0
  Matrix = DVineMatGen(d)
  M <- Matrix
  Mold <- M
  o <- diag(M)
  M <- reorderRVineMatrix(M)
  MaxMat <- createMaxMat(M)
  #CondDistr <- neededCondDistr(M)
  Types <- cbind(rbind(0,currentDvine$RVM$family),0)
  Params <- cbind(rbind(0,currentDvine$RVM$par),0)
  Params2 <- cbind(rbind(0,currentDvine$RVM$par2),0)
  V <- list()
  V$direct <- array(NA, dim = c(n, n, N))
  V$indirect <- array(NA, dim = c(n, n, N))
  for (l in 1:N){
    V$direct[ , ,l] <- cbind(rbind(NA,currentDvine$V$direct[ , ,l]),NA)
    V$indirect[ , ,l] <- cbind(rbind(NA,currentDvine$V$indirect[ , ,l]),NA)
  }
  V$direct[n, n, ] <- newcolumn
  doEst <- function(i) {
    m <- MaxMat[k, i]
    zr1 <- V$direct[k, i, ]
    zr2 <- if (m == M[k, i]) {
      V$direct[k, (d - m + 1), ]
    }else {
      V$indirect[k, (d - m + 1), ]
    }
    cfit <- BiCopSelect(zr2, zr1, familyset, selectioncrit,
                        indeptest, level, weights = NA, rotations)
    direct <- indirect <- NULL
    myHfunc = BiCopHfunc(zr1, zr2, cfit, check.pars = FALSE)
    #    if (CondDistr$direct[k - 1, i])
    direct <- myHfunc$hfunc2
    #    if (CondDistr$indirect[k - 1, i])
    indirect <- myHfunc$hfunc1
    list(direct = direct, indirect = indirect, cfit = cfit)
  }
  for (k in d:2) {
    i <- k-1
    res <- doEst(i)
    Types[k, i] <- res$cfit$family
    Params[k, i] <- res$cfit$par
    Params2[k, i] <- res$cfit$par2
    if (!is.null(res$direct))
      V$direct[k - 1, i, ] <- res$direct
    if (!is.null(res$indirect))
      V$indirect[k - 1, i, ] <- res$indirect
  }
  RVM <- RVineMatrix(Mold, family = Types, par = Params, par2 = Params2,
                     names = paste("V", 1:n, sep = ""))
  return(list(RVM=RVM, V=V))
}

### help functions

createMaxMat <- function(Matrix){

  if(dim(Matrix)[1]!=dim(Matrix)[2])
    stop("Structure matrix has to be quadratic.")

  MaxMat <- reorderRVineMatrix(Matrix)
  n  <- nrow(MaxMat)
  for(j in 1:(n-1)){
    for(i in (n-1):j){
      MaxMat[i,j] <- max(MaxMat[i:(i+1),j])
    }
  }

  tMaxMat <- MaxMat
  tMaxMat[is.na(tMaxMat)] <- 0
  oldSort <- diag(Matrix)
  oldSort <- oldSort[n:1]
  for(i in 1:n){
    MaxMat[tMaxMat == i] <- oldSort[i]
  }

  return(MaxMat)
}

reorderRVineMatrix <- function(Matrix){
  oldOrder <- diag(Matrix)
  O <- apply(t(1:nrow(Matrix)),2,"==", Matrix)
  for(i in 1:nrow(Matrix)){
    Matrix[O[,oldOrder[i]]] <- nrow(Matrix)-i+1
  }
  return(Matrix)
}
