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

DVRnp = function(y, X, margins.known=FALSE, margins="nonpar", margin.pars=list(),
               selectioncrit="BIC", familyset=NA, indeptestlevel=0.05, use_pobs=FALSE, correctionterm=TRUE,
               bw.sel="density", progress=FALSE, ...){
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

  vine.fit <- kdevinecop(cbind(v, U), matrix = DVineMatGen(p), ...)


  # return the best Dvine
  DVR.object = list(DVM=vine.fit,margins=margins,margin.pars=margin.pars,
                    my.index=1:d, used=1:d)
  return(DVR.object)
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
