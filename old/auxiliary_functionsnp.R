
pkernel <- function(x,margin.pars){
    centers = margin.pars$x
    bw = margin.pars$bw
    n=length(centers)
    X=matrix(rep(x,n),ncol=n)
    Z=sweep(X,2,centers,FUN="-")/bw
    y=apply(pnorm(Z),1,mean)
    return(y)
}

dkernel <- function(x,margin.pars){
    centers = margin.pars$x
    bw = margin.pars$bw
    n=length(centers)
    X=matrix(rep(x,n),ncol=n)
    Z=sweep(X,2,centers,FUN="-")/bw
    y=apply(dnorm(Z),1,mean)/bw
    return(y)
}

qkernel <-
    function(q,margin.pars,eps=0.00001){
        centers = margin.pars$x
        bw = margin.pars$bw
        z=rep(0,length(q))
        # sort quantiles
        my.order<-order(q)
        q<-sort(q)
        lower<-min(centers)-4*bw
        upper<-max(centers)+4*bw
        for (i in 1:length(z)){
            # fix quantiles near 0 and 1
            if (q[i]<eps){
                z[my.order[i]]=min(centers)-4*bw
            }
            else if  ((1-q[i])<eps){
                z[my.order[i]]=max(centers) +4*bw
            }
            else{
                f=function(x) pkernel(x,margin.pars)-q[i]
                if (i>1){
                    lower<-z[my.order[i-1]]
                }
                z[my.order[i]]=uniroot(f,lower=lower,upper=upper)$root
            }
        }
        #z<-z[my.order]
        return(z)
    }


qkernel2 <-
    function(q,margin.pars,eps=0.00001){
        centers = margin.pars$x
        bw = margin.pars$h
        z=rep(0,length(q))
        # sort quantiles
        my.order<-order(q)
        q<-sort(q)
        lower<-min(centers)-4*bw
        upper<-max(centers)+4*bw
        for (i in 1:length(z)){
            # fix quantiles near 0 and 1
            if (q[i]<eps){
                z[my.order[i]]=min(centers)-4*bw
            }
            else if  ((1-q[i])<eps){
                z[my.order[i]]=max(centers) +4*bw
            }
            else{
                f=function(x) predict(margin.pars,x=x)-q[i]
                if (i>1){
                    lower<-z[my.order[i-1]]
                }
                z[my.order[i]]=uniroot(f,lower=lower,upper=upper)$root
            }
        }
        #z<-z[my.order]
        return(z)
    }

# x: vector
# margin: string indicating margin class
# margin.pars: list containing margin parameters
CDF = function(x,margin,margin.pars){
    if (margin=="unif")return(x)
    if (margin=="gauss")return(pnorm(x,margin.pars$mean,margin.pars$sd))
    if (margin=="exp")return(pexp(x,rate=margin.pars$rate))
    if (margin=="nonpar")return(pkernel(x,margin.pars))
    if (margin=="kcde")return(predict(margin.pars,x=x))
    if (margin=="ecdf")return(margin.pars$fac*ecdf(margin.pars$x)(x))
    if (margin=="t")return(pt(x,margin.pars$df))
    if (margin=="skewt")return(pst(x,margin.pars$xi,margin.pars$omega,margin.pars$alpha,margin.pars$nu))
    if (margin=="skewnormal")return(psn(x,margin.pars$xi,margin.pars$omega,margin.pars$alpha))
}

# x: vector
# margin: string indicating margin class
# margin.pars: list containing margin parameters
PDF = function(x,margin,margin.pars){
    if (margin=="unif")return(1)
    if (margin=="gauss")return(dnorm(x,margin.pars$mean,margin.pars$sd))
    if (margin=="exp")return(dexp(x,rate=margin.pars$rate))
    if (margin=="nonpar")return(dkernel(x,margin.pars))
    if (margin=="t")return(dt(x,margin.pars$df))
}

# q: vector
# margin: string indicating margin class
# margin.pars: list containing margin parameters
qFUN = function(q,margin,margin.pars){
    if (margin=="unif")return(q)
    if (margin=="gauss")return(qnorm(q,margin.pars$mean,margin.pars$sd))
    if (margin=="exp")return(qexp(q,rate=margin.pars$rate))
    if (margin=="nonpar")return(qkernel(q,margin.pars))
    if (margin=="kcde")return(qkernel2(q,margin.pars))
    if (margin=="ecdf")return(quantile(margin.pars$x,min(q/margin.pars$fac,1)))
    if (margin=="t")return(qt(q,margin.pars$df))
    if (margin=="cauchy")return(qcauchy(q,margin.pars$location,margin.pars$scale))
    if (margin=="skewt"){
        if (is.matrix(q)){
            mymat = matrix(NA,dim(q)[1],dim(q)[2])
            for (k in 1:dim(q)[1]){
                mymat[k,] = qst(q[k,],margin.pars$xi,margin.pars$omega,margin.pars$alpha,margin.pars$nu)
            }
            return(mymat)
        }
        return(qst(q,margin.pars$xi,margin.pars$omega,margin.pars$alpha,margin.pars$nu))
    }
    if (margin=="skewnormal"){
        if (is.matrix(q)){
            mymat = matrix(NA,dim(q)[1],dim(q)[2])
            for (k in 1:dim(q)[1]){
                mymat[k,] = qsn(q[k,],margin.pars$xi,margin.pars$omega,margin.pars$alpha)
            }
            return(mymat)
        }
        return(qsn(q,margin.pars$xi,margin.pars$omega,margin.pars$alpha))
    }
}


bounds = function(margin){
    if (margin=="unif")return(0:1)
    if (margin=="exp")return(c(0,Inf))
    return(c(-Inf,Inf))
}

# x: vector
# margin: string indicating margin class
# returns list of marginal parameters
estim.margins = function(x,margin,bw.sel="density"){
    if (margin=="gauss"){
        my.mean<-mean(x)
        my.sd<-sd(x)
        return(list(mean=my.mean,sd=my.sd))
    }
    #  if (margin=="nonpar")return(list(x=x, bw=density(x)$bw))
    #  if (margin=="nonpar")return(list(x=x, bw=kCDF(x)$bw))
    if (margin=="nonpar"){
        if (bw.sel=="density"){return(list(x=x, bw=density(x)$bw))}else{
            return(list(x=x, bw=kcde(x)$h))
        }
    }
    if (margin=="kcde")return(kcde(x))
    if (margin=="ecdf")return(list(x=x,fac=length(x)/(length(x)+1)))
}

# Generates a D-Vine Matrix of dimension d. An order can be specified with the elements argument.
# Author: Daniel Kraus

DVineMatGen = function(d=length(elements),elements=1:d){
    if (!setequal(1:d,elements))
        stop("'elements' has to consist of the numbers from 1 to d.")
    Mat = diag(elements)
    for (i in 2:d){
        Mat[d+2-i,1:(d+1-i)]=elements[i:d]
    }
    return(Mat)
}

cond_quant_normal = function(alpha,x_c,mean,sigma){
    d = length(mean)
    sigma_12 = sigma[1,2:d]
    sigma_11 = sigma[1,1]
    sigma_22 = sigma[2:d,2:d]
    mean_c = mean[1]+sigma_12%*%solve(sigma_22)%*%as.vector(x_c-mean[2:d])
    sigma_c = sigma_11-sigma_12%*%solve(sigma_22)%*%sigma_12
    return(qnorm(alpha,mean=mean_c, sd=sqrt(sigma_c)))
}

cond_dens_normal = function(x,x_c,mean,sigma){
    d = length(mean)
    sigma_12 = sigma[1,2:d]
    sigma_11 = sigma[1,1]
    sigma_22 = sigma[2:d,2:d]
    mean_c = mean[1]+sigma_12%*%solve(sigma_22)%*%as.vector(x_c-mean[2:d])
    sigma_c = sigma_11-sigma_12%*%solve(sigma_22)%*%sigma_12
    return(dnorm(x,mean=mean_c, sd=sqrt(sigma_c)))
}


cond_quant_t = function(alpha,x_c,mean,sigma,df){
    d = length(mean)
    df_c = df+d-1
    sigma_12 = sigma[1,2:d]
    sigma_11 = sigma[1,1]
    sigma_22 = sigma[2:d,2:d]
    mean_c = mean[1]+sigma_12%*%solve(sigma_22)%*%(x_c-mean[2:d])
    sigma_c = (df+(x_c-mean[2:d])%*%solve(sigma_22)%*%(x_c-mean[2:d]))/(df_c)*(sigma_11-sigma_12%*%solve(sigma_22)%*%sigma_12)
    return(qt(alpha,df=df_c)*sqrt(sigma_c)+mean_c)
}

cond_CDFsnp = function(u,DVM){

    d = ncol(DVM$matrix)
    if(!all(DVM$matrix==DVineMatGen(d)))warning("Input has to be D-vine")
    if(length(u)!=d-1)stop("Dimensions of u and DVM are not compatible")
    if(d==2)return(u)
    V = matrix(NA,d,d)
    V[d,-1] = u
    V2 = V
    for (j in (d-1):2){
        for (k in (d-1):j){
            # temp = BiCopHfunc(V2[k+1,j],V[k+1,j+1],Fam[k+1,j],Par[k+1,j],Par2[k+1,j])
            cfit <-  DVM[[d - k]][[j]]$c
            cfit$flip <- TRUE # ifelse(is.null(cfit$flip), TRUE, NULL)
            V2[k,j] = hkdecop(cbind(V2[k+1,j], V[k+1,j+1]), cfit, cond.var = 2)
            V[k,j] = hkdecop(cbind(V2[k+1,j], V[k+1,j+1]), cfit, cond.var = 1)
        }
    }
    return(diag(V2)[-1])
}

cond_dens <- function(object, y, X, ...){
    X = as.matrix(X)
    uev <- cbind(as.matrix(y),matrix(X,length(y),dim(object$par)[1]-1,byrow=dim(X)[2]==1))
    if(dim(uev)[2]==1)  uev <- matrix(uev, 1, dim(uev)[1])
    if(any(uev>1) || any(uev<0)) stop("Data has be in the interval [0,1].")
    n = dim(uev)[2]
    N = dim(uev)[1]
    if(dim(uev)[2] != dim(object$Matrix)[2]) stop("Dimensions of 'data'
                                                and 'object' do not match.")
    if(!is(object,"RVineMatrix")) stop("'object' has to be an
                                     RVineMatrix object")

    o = diag(object$Matrix)
    oldobject = object
    if(any(o != length(o):1)){
        object = normalizeRVineMatrix(object)
        uev = matrix(uev[, o[length(o):1]], N, n)
    }

    Family=object$family
    Params=object$par
    Params2=object$par2

    val = array(1,dim=c(n,n,N))
    V = list()
    V$direct = array(NA,dim=c(n,n,N))
    V$indirect = array(NA,dim=c(n,n,N))

    V$direct[n,,] = t(uev[,n:1])

    for(i in (n-1):1)
    {
        for(k in n:(i+1))
        {

            m = object$MaxMat[k,i]
            zr1 = V$direct[k,i,]

            if(m == object$Matrix[k,i]){
                zr2 = V$direct[k,(n-m+1),]
            }else{
                zr2 = V$indirect[k,(n-m+1),]
            }

            val[k-1,i,] <- BiCopPDF(zr2, zr1, Family[k,i], Params[k,i],
                                    Params2[k,i])

            if(object$CondDistr$direct[k-1,i])
            {
                V$direct[k-1,i,] = .C("Hfunc1",
                                      as.integer(object$family[k,i]),
                                      as.integer(length(zr1)),
                                      as.double(zr1),
                                      as.double(zr2),
                                      as.double(Params[k,i]),
                                      as.double(Params2[k,i]),
                                      as.double(rep(0,length(zr1))),
                                      PACKAGE='VineCopula')[[7]]
                V$direct[k-1,i,] <- pmin(pmax(1e-16, V$direct[k-1,i,]),
                                         1-1e-16)
            }
            if(object$CondDistr$indirect[k-1,i])
            {
                V$indirect[k-1,i,] = .C("Hfunc2",
                                        as.integer(object$family[k,i]),
                                        as.integer(length(zr2)),
                                        as.double(zr2),
                                        as.double(zr1),
                                        as.double(Params[k,i]),
                                        as.double(Params2[k,i]),
                                        as.double(rep(0,length(zr1))),
                                        PACKAGE='VineCopula')[[7]]
                V$indirect[k-1,i,] <- pmin(pmax(1e-16,
                                                V$indirect[k-1,i,]), 1-1e-16)
            }

        }
    }
    apply(val, 3, function(X){prod(X[,1])})
}

normalizeRVineMatrix = function(RVM){

    oldOrder = diag(RVM$Matrix)
    Matrix = reorderRVineMatrix(RVM$Matrix)

    return(RVineMatrix(Matrix, RVM$family, RVM$par, RVM$par2, names =
                           rev(RVM$names[oldOrder])))
}

reorderRVineMatrix = function(Matrix){
    oldOrder = diag(Matrix)

    O = apply(t(1:nrow(Matrix)),2,"==", Matrix)

    for(i in 1:nrow(Matrix)){
        Matrix[O[,oldOrder[i]]] = nrow(Matrix)-i+1
    }

    return(Matrix)
}


permutations = function (n, r, v = 1:n, set = TRUE, repeats.allowed = FALSE)
{
    if (mode(n) != "numeric" || length(n) != 1 || n < 1 || (n%%1) !=
        0)
        stop("bad value of n")
    if (mode(r) != "numeric" || length(r) != 1 || r < 1 || (r%%1) !=
        0)
        stop("bad value of r")
    if (!is.atomic(v) || length(v) < n)
        stop("v is either non-atomic or too short")
    if ((r > n) & repeats.allowed == FALSE)
        stop("r > n and repeats.allowed=FALSE")
    if (set) {
        v <- unique(sort(v))
        if (length(v) < n)
            stop("too few different elements")
    }
    v0 <- vector(mode(v), 0)
    if (repeats.allowed)
        sub <- function(n, r, v) {
            if (r == 1)
                matrix(v, n, 1)
            else if (n == 1)
                matrix(v, 1, r)
            else {
                inner <- Recall(n, r - 1, v)
                cbind(rep(v, rep(nrow(inner), n)), matrix(t(inner),
                                                          ncol = ncol(inner), nrow = nrow(inner) * n,
                                                          byrow = TRUE))
            }
        }
    else sub <- function(n, r, v) {
        if (r == 1)
            matrix(v, n, 1)
        else if (n == 1)
            matrix(v, 1, r)
        else {
            X <- NULL
            for (i in 1:n) X <- rbind(X, cbind(v[i], Recall(n -
                                                                1, r - 1, v[-i])))
            X
        }
    }
    sub(n, r, v[1:n])
}
