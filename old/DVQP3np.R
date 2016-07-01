DVQPnp = function(X,DVR.object,alpha=0.95,uscale=F){
  X = as.matrix(X)
  N = nrow(X)
  d = ncol(X)
  if (uscale==T){
    U = X
  }else{
    U = matrix(0,N,d)
    for(j in DVR.object$used){
      U[,j] = CDF(X[,j],DVR.object$margins[j+1],DVR.object$margin.pars[[j+1]])
    }
  }
  U = U[,DVR.object$used, drop=F]
  
  RVM = DVR.object$DVM
  my.order = diag(RVM$matrix)[-1]-1
  U = U[,my.order, drop=F]
  RVM$matrix = DVineMatGen(ncol(RVM$matrix))
  
  helpfun = function(u){
    temp = cond_CDFsnp(u,RVM)
    out = vector(length=length(alpha))
    for (a in 1:length(alpha)){
      pred_quant = rkdevinecop(1, RVM, U = matrix(c(alpha[a],temp), nrow = 1))
      if (any(pred_quant<1e-9)|any(pred_quant>1-1e-9))warning("Capping occured!")
      out[a] = qFUN(pred_quant[1], DVR.object$margins[1],DVR.object$margin.pars[[1]])
    }
    return(out)
  }
  quantile = apply(U,1,helpfun)
  if(is.matrix(quantile)){
    quantile = t(quantile)
    colnames(quantile) = alpha
  }
  return(quantile)
}
