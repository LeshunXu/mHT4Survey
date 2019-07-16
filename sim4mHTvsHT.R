##########################################################################
### Simulation
### True Variance: It is the variance calculated by 100 times simulation
###                It is related to the size of domains. (e.g. for a
###                20 by 20 squre, its var is different to that of a 
###                30 by 30 one.)
###
### True Mean: It is also depends on the size of domains, and
###            it is also simulated 100 times for the mean.
###
### Author: Leshun Xu
### Last modified date: 2019-07-15
##########################################################################


# setwd("working directory")
library(MASS)

nedge <- 50        # the number of the colums and rows of the Z field
mpts <- nedge^2    # the number of individuals in the Z field

p2p.dis <- matrix(0, mpts, mpts)
for(i in 1:mpts){
  for(j in i:mpts){
    p1 <- (i-1)%/%nedge
    p2 <- (i-1)%%nedge
    q1 <- (j-1)%/%nedge
    q2 <- (j-1)%%nedge
    dis <- sqrt((p1-q1)^2+(p2-q2)^2)
    p2p.dis[i,j] <- dis
    p2p.dis[j,i] <- dis
  }
}

sigma <- exp(-0.05*p2p.dis) # transfer distance into covariances

set.seed(2019)
seqn <- mvrnorm(n = 1, rep(0, mpts), sigma)
pix <- matrix(seqn, nedge, nedge, byrow = TRUE) # generate the Z field

# plot the Z field
pixtt <- pix[,nrow(pix):1]
image(pixtt, col=heat.colors(100), axes=FALSE)


g <- function(z){
  exp(z)/(1+exp(z))
}

f <- function(z){
  z                    # case1: f(z)=z
  # z^3-3*z^2+2*z+0.01   # case2: f(z)=z^3-3*z^2+2*z+0.01
  # 2*sin(z)             # case3: f(z)=2*sin(z)
  # exp(z)-5             # case4: f(z)=exp(z)
}

# This function is for the inclusion probabilities of HT estimators.
selec.prob <- function(ind1, ind2, dn, m1, m2, k){
  p1 <- m1/dn; p2 <- m2/dn
  pi1 <- 1-(1-p1)^k
  pi2 <- 1-(1-p2)^k
  pi12 <- 2*(k-1)*pi1*pi2/(2*k-pi1-pi2)
  c(pi1, pi2, pi12)
}

field.z <- pix
field.f <- f(field.z)
field.g <- g(field.z)

##########################################################################
### Simulation for the true variance and true mean (start here)
SumSnSeq <- data.frame(matrix(ncol = nedge, nrow = 0))
colnames(SumSnSeq) <- paste("Sn", 1:nedge, sep = "")

t <- 0
for (i in 1:100) {
  t <- t+1
  set.seed(i)
  seqn <- mvrnorm(n = 1, rep(0, mpts), sigma)
  pix.t <- matrix(seqn, nedge, nedge, byrow = TRUE)
  field.z.t <- pix.t
  field.f.t <- f(field.z.t)
  vecTemp <- NULL
  for (j in 1:nedge) {
    vecTemp <- c(vecTemp, sum(field.f.t[1:j,1:j]))
  }
  SumSnSeq[t,] <- vecTemp
}
### Simulation for the true variance and true mean (end here)
##########################################################################

set.seed(2)
temp.r <- rbinom(mpts, 1, field.g)
field.r <- matrix(temp.r, nedge, nedge, byrow = FALSE)

field.x <- field.f*field.r

size <- seq(40, nedge, 1)

recorded <- data.frame(matrix(ncol = 11, nrow = 0))
colnames(recorded) <- c("SampleSize", "DnSize", "DnVar",
                        "NEW.var","NEW.err1",
                        "HT.var", "HTvar.err3",
                        "Sum.MHT", "Sum.MHT.err",
                        "Sum.HT", "Sum.HT.err"
)

t <- 0

for (s in size) {
  t <- t+1
  r01 <- matrix(0, s, nedge-s)
  r02 <- matrix(0, nedge-s, nedge)
  Rn <- field.r[1:s, 1:s]
  Rn2field <- cbind(Rn, r01)
  Rn2field <- rbind(Rn2field, r02)
  index <- which(Rn2field!=0)
  sample.cov <- sigma[index, index]
  
  Dn.1 <- matrix(1, s, s)
  temp.field <- cbind(Dn.1, r01)
  temp.field <- rbind(temp.field, r02)
  index2 <- which(temp.field!=0)
  Dn.f <- field.f[index2]
  
  sumDn.true <- sum(Dn.f)
  varDn.true <- var(SumSnSeq[,s])
  
  
  f.seq <- field.f[index]   # Observed as a sample
  g.seq <- field.g[index]
  
  n.seq <- length(index)    # it is the sample size
  Dn.size <- length(index2) # It is the finite population size
  
  varDn.est.HT <- 0         # This is H-T estimator of the variance
  sumDn.est.MHT <- 0        # This is MHT estimator of the pop total
  sumDn.est.HT <- 0         # This is H-T estimator of the pop total
  for (i in 1:n.seq) {
    for (j in 1:n.seq) {
      if (i==j) {
        selec.s <- selec.prob(index[i], index[j], s^2, 1, 1, n.seq)
        selec.i <- selec.s[1]
        varDn.est.HT <- varDn.est.HT + f.seq[i]^2*(1-selec.i)/selec.i^2
      }else{
        selec.s <- selec.prob(index[i], index[j], s^2, 1, 1, n.seq)
        selec.i <- selec.s[1]
        selec.j <- selec.s[2]
        selec.ij <- selec.s[3]
        varDn.est.HT <- varDn.est.HT + f.seq[i]*f.seq[j]*(selec.ij - selec.i * selec.j)/selec.ij/selec.i/selec.j
      }
    }
    sumDn.est.MHT <- sumDn.est.MHT + f.seq[i]/g.seq[i]
    
    selec.s <- selec.prob(index[i], index[i], s^2, 1, 1, n.seq)
    selec.i <- selec.s[1]
    sumDn.est.HT <- sumDn.est.HT + f.seq[i]/selec.i
  }
  
  field.xtmp <- field.x
  
  ##########################################################################
  ### The mHT estimator of the variance (start here)
  
  kn <- floor(s^0.7)   # the bolck size is kn*kn
  mn <- floor(50/kn)   # the number of blocks is mn*mn
  vbar <- sumDn.est.MHT/mn/mn
  figi <- field.xtmp/field.g
  sumtmp <- 0

  ### Note: In the following, blue blocks are neighbours of a red block.
  ## Case 1: Number of red blocks=1; Number of blue blocks=3;
  c1rblock <- figi[1:kn, 1:kn]
  term1 <- sum(c1rblock)-vbar
  sumtmp <- sumtmp + term1^2
  c1bblock1 <- figi[(kn+1):(2*kn), 1:kn]
  sumtmp <- sumtmp +term1 * (sum(c1bblock1) - vbar)
  c1bblock2 <- figi[1:kn, (kn+1):(2*kn)]
  sumtmp <- sumtmp +term1 * (sum(c1bblock2) - vbar)
  c1bblock3 <- figi[(kn+1):(2*kn), (kn+1):(2*kn)]
  sumtmp <- sumtmp +term1 * (sum(c1bblock3) - vbar)
  
  ## Case 2: Number of red blocks=mn-2; Number of blue blocks=5;
  for (m in 1:(mn-2)) {
    c2rblock <- figi[(m*kn+1):((m+1)*kn), 1:kn]
    term1 <- sum(c2rblock)-vbar
    sumtmp <- sumtmp + term1^2
    c2bblock1 <- figi[((m-1)*kn+1):(m*kn), 1:kn]               # on its left
    sumtmp <- sumtmp +term1 * (sum(c2bblock1) - vbar)
    c2bblock2 <- figi[((m+1)*kn+1):((m+2)*kn), 1:kn]           # on its right
    sumtmp <- sumtmp +term1 * (sum(c2bblock2) - vbar)
    c2bblock3 <- figi[((m-1)*kn+1):(m*kn), (kn+1):(kn+kn)]     # left-bottom
    sumtmp <- sumtmp +term1 * (sum(c2bblock3) - vbar)
    c2bblock4 <- figi[(m*kn+1):((m+1)*kn), (kn+1):(kn+kn)]     # under it
    sumtmp <- sumtmp +term1 * (sum(c2bblock4) - vbar)
    c2bblock5 <- figi[((m+1)*kn+1):((m+2)*kn), (kn+1):(kn+kn)] # right-bottom
    sumtmp <- sumtmp +term1 * (sum(c2bblock5) - vbar)
  }
  
  ## Case 3:  Number of red blocks=1; Number of blue blocks=3;
  c3rblock <- figi[((mn-1)*kn+1):(mn*kn), 1:kn]
  term1 <- sum(c3rblock)-vbar
  sumtmp <- sumtmp + term1^2
  c3bblock1 <- figi[((mn-2)*kn+1):((mn-1)*kn), 1:kn]
  sumtmp <- sumtmp +term1 * (sum(c3bblock1) - vbar)
  c3bblock2 <- figi[((mn-2)*kn+1):((mn-1)*kn), (kn+1):(2*kn)]
  sumtmp <- sumtmp +term1 * (sum(c3bblock2) - vbar)
  c3bblock3 <- figi[((mn-1)*kn+1):(mn*kn), (kn+1):(2*kn)]
  sumtmp <- sumtmp +term1 * (sum(c3bblock3) - vbar)
  
  ## Case 4: Number of red blocks=mn-2; Number of blue blocks=5;
  for (m in 1:(mn-2)) {
    c4rblock <- figi[1:kn, (m*kn+1):((m+1)*kn)]
    term1 <- sum(c4rblock)-vbar
    sumtmp <- sumtmp + term1^2
    c4bblock1 <- figi[1:kn, ((m-1)*kn+1):(m*kn)]               # above it
    sumtmp <- sumtmp +term1 * (sum(c4bblock1) - vbar)
    c4bblock2 <- figi[1:kn, ((m+1)*kn+1):((m+2)*kn)]           # under it
    sumtmp <- sumtmp +term1 * (sum(c4bblock2) - vbar)
    c4bblock3 <- figi[(kn+1):(kn+kn), ((m-1)*kn+1):(m*kn)]     # right-top
    sumtmp <- sumtmp +term1 * (sum(c4bblock3) - vbar)
    c4bblock4 <- figi[(kn+1):(kn+kn), (m*kn+1):((m+1)*kn)]     # right
    sumtmp <- sumtmp +term1 * (sum(c4bblock4) - vbar)
    c4bblock5 <- figi[(kn+1):(kn+kn), ((m+1)*kn+1):((m+2)*kn)] # right-bottom
    sumtmp <- sumtmp +term1 * (sum(c4bblock5) - vbar)
  }
  
  ## Case 5: Number of red blocks=(mn-2)^2; Number of blue blocks=8;
  for (m1 in 1:(mn-2)) {
    for (m2 in 1:(mn-2)) {
      c5rblock <- figi[(m1*kn+1):((m1+1)*kn), (m2*kn+1):((m2+1)*kn)]
      term1 <- sum(c5rblock)-vbar
      sumtmp <- sumtmp + term1^2
      c5bblock1 <- figi[((m1-1)*kn+1):(m1*kn), ((m2-1)*kn+1):(m2*kn)]
      sumtmp <- sumtmp +term1 * (sum(c5bblock1) - vbar)
      c5bblock2 <- figi[(m1*kn+1):((m1+1)*kn), ((m2-1)*kn+1):(m2*kn)]
      sumtmp <- sumtmp +term1 * (sum(c5bblock2) - vbar)
      c5bblock3 <- figi[((m1+1)*kn+1):((m1+2)*kn), ((m2-1)*kn+1):(m2*kn)]
      sumtmp <- sumtmp +term1 * (sum(c5bblock3) - vbar)
      
      c5bblock4 <- figi[((m1-1)*kn+1):(m1*kn), (m2*kn+1):((m2+1)*kn)]
      sumtmp <- sumtmp +term1 * (sum(c5bblock4) - vbar)
      c5bblock5 <- figi[((m1+1)*kn+1):((m1+2)*kn), (m2*kn+1):((m2+1)*kn)]
      sumtmp <- sumtmp +term1 * (sum(c5bblock5) - vbar)
      
      c5bblock6 <- figi[((m1-1)*kn+1):(m1*kn), ((m2+1)*kn+1):((m2+2)*kn)]
      sumtmp <- sumtmp +term1 * (sum(c5bblock6) - vbar)
      c5bblock7 <- figi[(m1*kn+1):((m1+1)*kn), ((m2+1)*kn+1):((m2+2)*kn)]
      sumtmp <- sumtmp +term1 * (sum(c5bblock7) - vbar)
      c5bblock8 <- figi[((m1+1)*kn+1):((m1+2)*kn), ((m2+1)*kn+1):((m2+2)*kn)]
      sumtmp <- sumtmp +term1 * (sum(c5bblock8) - vbar)
    }
  }
  
  ## Case 6: Number of red blocks=mn-2; Number of blue blocks=5;
  for (m in 1:(mn-2)) {
    c6rblock <- figi[((mn-1)*kn+1):(mn*kn), (m*kn+1):((m+1)*kn)]
    term1 <- sum(c6rblock)-vbar
    sumtmp <- sumtmp + term1^2
    c6bblock1 <- figi[((mn-1)*kn+1):(mn*kn), ((m-1)*kn+1):(m*kn)]
    sumtmp <- sumtmp +term1 * (sum(c6bblock1) - vbar)
    c6bblock2 <- figi[((mn-1)*kn+1):(mn*kn), ((m+1)*kn+1):((m+2)*kn)]
    sumtmp <- sumtmp +term1 * (sum(c6bblock2) - vbar)
    c6bblock3 <- figi[((mn-2)*kn+1):((mn-1)*kn), ((m-1)*kn+1):(m*kn)]
    sumtmp <- sumtmp +term1 * (sum(c6bblock3) - vbar)
    c6bblock4 <- figi[((mn-2)*kn+1):((mn-1)*kn), (m*kn+1):((m+1)*kn)]
    sumtmp <- sumtmp +term1 * (sum(c6bblock4) - vbar)
    c6bblock5 <- figi[((mn-2)*kn+1):((mn-1)*kn), ((m+1)*kn+1):((m+2)*kn)]
    sumtmp <- sumtmp +term1 * (sum(c6bblock5) - vbar)
  }
  
  ## Case 7: Number of red blocks=1; Number of blue blocks=3;
  c7rblock <- figi[1:kn, ((mn-1)*kn+1):(mn*kn)]
  term1 <- sum(c7rblock)-vbar
  sumtmp <- sumtmp + term1^2
  c7bblock1 <- figi[(kn+1):(2*kn), ((mn-1)*kn+1):(mn*kn)]
  sumtmp <- sumtmp +term1 * (sum(c7bblock1) - vbar)
  c7bblock2 <- figi[1:kn,  ((mn-2)*kn+1):((mn-1)*kn)]
  sumtmp <- sumtmp +term1 * (sum(c7bblock2) - vbar)
  c7bblock3 <- figi[(kn+1):(2*kn), ((mn-2)*kn+1):((mn-1)*kn)]
  sumtmp <- sumtmp +term1 * (sum(c7bblock3) - vbar)
  
  ## Case 8: Number of red blocks=mn-2; Number of blue blocks=5;
  for (m in 1:(mn-2)) {
    c8rblock <- figi[(m*kn+1):((m+1)*kn), ((mn-1)*kn+1):(mn*kn)]
    term1 <- sum(c8rblock)-vbar
    sumtmp <- sumtmp + term1^2
    c8bblock1 <- figi[((m-1)*kn+1):(m*kn), ((mn-1)*kn+1):(mn*kn)]
    sumtmp <- sumtmp +term1 * (sum(c8bblock1) - vbar)
    c8bblock2 <- figi[((m+1)*kn+1):((m+2)*kn), ((mn-1)*kn+1):(mn*kn)]
    sumtmp <- sumtmp +term1 * (sum(c8bblock2) - vbar)
    c8bblock3 <- figi[((m-1)*kn+1):(m*kn), ((mn-2)*kn+1):((mn-1)*kn)]
    sumtmp <- sumtmp +term1 * (sum(c8bblock3) - vbar)
    c8bblock4 <- figi[(m*kn+1):((m+1)*kn), ((mn-2)*kn+1):((mn-1)*kn)]
    sumtmp <- sumtmp +term1 * (sum(c8bblock4) - vbar)
    c8bblock5 <- figi[((m+1)*kn+1):((m+2)*kn), ((mn-2)*kn+1):((mn-1)*kn)]
    sumtmp <- sumtmp +term1 * (sum(c8bblock5) - vbar)
  }
  
  ## Case 9:  Number of red blocks=1; Number of blue blocks=3;
  c9rblock <- figi[((mn-1)*kn+1):(mn*kn), ((mn-1)*kn+1):(mn*kn)]
  term1 <- sum(c9rblock)-vbar
  sumtmp <- sumtmp + term1^2
  c9bblock1 <- figi[((mn-2)*kn+1):((mn-1)*kn), ((mn-1)*kn+1):(mn*kn)]
  sumtmp <- sumtmp +term1 * (sum(c9bblock1) - vbar)
  c9bblock2 <- figi[((mn-2)*kn+1):((mn-1)*kn), ((mn-2)*kn+1):((mn-1)*kn)]
  sumtmp <- sumtmp +term1 * (sum(c9bblock2) - vbar)
  c9bblock3 <- figi[((mn-1)*kn+1):(mn*kn), ((mn-2)*kn+1):((mn-1)*kn)]
  sumtmp <- sumtmp +term1 * (sum(c9bblock3) - vbar)
  
  NEW.MHT <- sumtmp - sum((1-g.seq)/g.seq/g.seq*f.seq^2)
  
  ### The mHT estimator of the variance (end here)
  ##########################################################################
  
  recorded[t,] <- 
    c(n.seq, s^2, varDn.true,
      NEW.MHT, abs((NEW.MHT-varDn.true)/varDn.true),
      varDn.est.HT, abs((varDn.est.HT-varDn.true)/varDn.true),
      sumDn.est.MHT, abs((sumDn.est.MHT-sumDn.true)/sumDn.true),
      sumDn.est.HT, abs((sumDn.est.HT-sumDn.true)/sumDn.true)
    )
  print(recorded[t,c(1,2,5,7,9,11)])
}

round(recorded[,c(1,2,5,7,9,11)],4)


yup <- max(abs(recorded$NEW.err1), abs(recorded$HTvar.err3))
ylo <- min(abs(recorded$NEW.err1), abs(recorded$HTvar.err3))

plot(abs(recorded$NEW.err1), type = "l", xlab = "Index of Domains", ylab = "ARE Errors",
     col = "red", lty=1, lwd=3, cex.lab=1.5, cex.axis=1.5,
     ylim = c(ylo,yup),xaxt = "n")
axis(1, at=1:11, labels=c(1:11))
points(abs(recorded$NEW.err1), col = "red", lty=1, lwd=2, pch=19)
lines(abs(recorded$HTvar.err3), lty=2, lwd=2, col="blue")
points(abs(recorded$HTvar.err3), col = "blue", lty=1, lwd=2, pch=19)
legend(7.6, .985, legend=c("mHT Estimator", "HT Estimator"),
       col=c("red", "blue"), lty=1:2, lwd=2, cex=1.4)

plot(abs(recorded$Sum.MHT.err), type = "l", xlab = "Index of Domains", ylab = "ARE Errors",
     ylim = c(0, max(abs(recorded$Sum.HT.err))),
     col = "red", lty=1, lwd=3, cex.lab=1.5, cex.axis=1.5,xaxt = "n")
axis(1, at=1:11, labels=c(1:11))
points(abs(recorded$Sum.MHT.err), col = "red", lty=1, lwd=2, pch=19)
lines(abs(recorded$Sum.HT.err), lty=2, lwd=2, col="blue")
points(abs(recorded$Sum.HT.err), col = "blue", lty=1, lwd=2, pch=19)
legend(7.5, .55, legend=c("mHT Estimator", "HT Estimator"),
       col=c("red", "blue"), lty=1:2, lwd=2, cex=1.5)

