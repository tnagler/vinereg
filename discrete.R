library(vinereg)
library(copula)
library(kdevine)
library(rafalib)
library(beepr)
library(np)
library(sn)
library(logitnorm)

n <- 100
d <- 3
cop <- gumbelCopula(dim = d, param = 1.3)
s <- rCopula(n, cop)
qfun <- function(u) qbinom(u, 4, 0.3)
x <- qfun(s)

trueqan <- function(x, q) 1 * x[, 1] + 0.5 * (x[, 2] - 2)^2 - 0.5 * (x[, 3] - 1.5) + qnorm(q)
y <- trueqan(x, runif(n))
# pairs(cbind(y, x))
z <- x + (matrix(qlogitnorm(runif(n * d), sigma = 0.5), n) - 0.5) * 1
# pairs(cbind(y, z))

resnp <- vinereg(y, z, familyset = 1, par_1d = list(mult = 2), correction = "BIC")

x <- as.data.frame(x)
predict(resnp, x)
resnp2 <- npqreg(tydat = y,
                 txdat = do.call(data.frame, lapply(1:ncol(x), function(i) ordered(x[, i], 0:25))),
                 tol = 1e-2, itmax = 10)


news <- rCopula(n, cop)
newx <- qfun(news)
alph <- 0.01
tq <- trueqan(newx, alph)

npq <- predict(resnp, newx, alpha = alph)
npq2 <- predict(resnp2,
                tau = alph,
                exdat = do.call(data.frame, lapply(1:ncol(x), function(i) ordered(newx[, i], 0:25))),
                tol = 1e-2, itmax = 10)


ran <- range(c(tq, npq, npq2))
mypar(1, 2)
plot(tq + rnorm(length(tq)) / 25, npq + rnorm(length(tq)) / 25,
     main = formatC(mean(abs(tq - npq)), 2), xlim = ran, ylim = ran)
abline(0, 1)
plot(tq + rnorm(length(tq)) / 25, npq2 + rnorm(length(tq)) / 25,
     main = formatC(mean(abs(tq - npq2)), 2), xlim = ran, ylim = ran)
abline(0, 1)

beep(5)
