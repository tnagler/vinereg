DVQP = function(X,DVR.object,alpha=0.95,uscale=F){
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
    my.order = diag(RVM$Matrix)[-1]-1
    U = U[,my.order, drop=F]
    RVM$Matrix = DVineMatGen(ncol(RVM$Matrix))

    helpfun = function(u){
        temp = cond_CDFs(u,RVM)
        out = vector(length=length(alpha))
        for (a in 1:length(alpha)){
            pred_quant = RVineSim(1,RVM,c(alpha[a],temp))

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
