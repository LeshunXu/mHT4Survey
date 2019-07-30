##########################################################################
### This file is for the estimation of NZ pop total by using small blocks.
###
### Author: Leshun Xu
### Last modified date: 2019-07-29
##########################################################################

# setwd("working directory")

emp <- read.csv("Emp.csv", header = TRUE, row.names=NULL, sep = ",")
pop <- read.csv("Pop.csv", header = FALSE, row.names=1, sep = ",")
pop <- t(pop)
pop <- as.data.frame(pop)
pop <- pop[which(pop$Year>=2006),]
emp <- emp[which(emp$Year>=2006),]
row.names(pop) <- c()
row.names(emp) <- c()

reunit <- 100000
# Setting up new unit for emp to avoid overflow in exp(Z)/(1+exp(Z))
Z <- emp[,2:18]/reunit
X <- pop[,3:19]/reunit
truepop <- pop[13,2]/reunit
region <- names(pop[,-c(1,2)])

# Calculating the g_i's based on the history data
# Being asumed in the form: g(Z) = exp(Z)/(1+exp(Z))
giAll <- exp(Z)/(1+exp(Z))

##########################################################################
### Setting up function X=f(Z) for each region by using 2006-2017 data
### (start here)
fi <- NULL          # fi is the estimated pop for each region in 2018
for (r in 1:17) {
  p <- X[1:12,r]
  z <- Z[1:12,r]
  pop2z <- data.frame(matrix(ncol = 2, nrow = 12))
  names(pop2z) <- c("Pop", "Z")
  
  pop2z$Z <- z
  pop2z$Pop <- p
  f2z <- lm(Pop ~ ., data = pop2z)
  pdata <- pop2z[1,]
  pdata$Z <- as.numeric(Z[13,r])
  fi <- c(fi, as.numeric(predict(f2z, newdata = pdata)))
}
### Setting up function X=f(Z) for each region by using 2006-2017 data
### (end here)
##########################################################################

# Setting up function g(Z) of each region for 2018 
gi <- as.numeric(giAll[13,])

# Set up E(gi) for each region by using 2006-2017 data
Egi <- NULL
for (r in 1:17) {
  Egi <- c(Egi, mean(giAll[1:12,r]))
}

##########################################################################
### Pop total estimation (start here)
set.seed(2019)
ri <- rbinom(17, 1, Egi)  # Generate ri field for 2018
region[which(ri==0)]      # regions that can not be observed

estpop <- sum(ri*fi/gi)
print(paste("True Population Total:", truepop*reunit))
print(paste("Estimated Population Total:", round(estpop,5)*reunit))
print(paste("Error of The Pop Total:", round((estpop - truepop)*reunit,0)))
print(paste("Absolute Relative Error:", abs(round((estpop - truepop)/truepop,5))))

### Pop total estimation (end here)
##########################################################################


##########################################################################
### Pop total's var estimation (start here)

# Setting up 10*10 blocks for 17 regions

xi <- matrix(0, 10,10)
xgi <- matrix(0, 10,10)

## Region 01:
xi[1:2, 1:2] <- fi[1]/4; xgi[1:2, 1:2] <- gi[1]

## Region 02:
xi[3:6, 1:2] <- fi[2]/34; xgi[3:6, 1:2] <- gi[2]
xi[1:6, 3:6] <- fi[2]/34; xgi[1:6, 3:6] <- gi[2]
xi[7, 2:3] <- fi[2]/34; xgi[7, 2:3] <- gi[2]

## Region 03:
xi[1:7, 7] <- fi[3]/10; xgi[1:7, 7] <- gi[3]
xi[7 ,4:6] <- fi[3]/10; xgi[7, 4:6] <- gi[3]

## Region 04:
xi[1, 8:10] <- fi[4]/6; xgi[1, 8:10] <- gi[4]
xi[2, 8:9] <- fi[4]/6; xgi[2, 8:9] <- gi[4]
xi[3, 8] <- fi[4]/6; xgi[3, 8] <- gi[4]

## Region 05:
xi[3, 9] <- fi[5]; xgi[3, 9] <- gi[5]

## Region 06:
xi[4, 8:9] <- fi[6]/3; xgi[4, 8:9] <- gi[6]
xi[5, 9] <- fi[6]/3; xgi[5, 9] <- gi[6]

## Region 07:
xi[8, 7:8] <- fi[7]/2; xgi[8, 7:8] <- gi[7]

## Region 08:
xi[5:7, 8] <- fi[8]/5; xgi[5:7, 8] <- gi[8]
xi[6:7, 9] <- fi[8]/5; xgi[6:7, 9] <- gi[8]

## Region 09:
xi[2:9, 10] <- fi[9]/11; xgi[2:9, 10] <- gi[9]
xi[8:9, 9] <- fi[9]/11; xgi[8:9, 9] <- gi[9]
xi[9, 8] <- fi[9]/11; xgi[9, 8] <- gi[9]

## Region 10:
xi[8, 1] <- fi[10]; xgi[8, 1] <- gi[10]

## Region 11:
xi[7, 1] <- fi[11]; xgi[7, 1] <- gi[11]

## Region 12:
xi[8, 2] <- fi[12]; xgi[8, 2] <- gi[12]

## Region 13:
xi[9, 3] <- fi[13]; xgi[9, 3] <- gi[13]

## Region 14:
xi[8, 3:6] <- fi[14]/12; xgi[8, 3:6] <- gi[14]
xi[9, 4:7] <- fi[14]/12; xgi[9, 4:7] <- gi[14]
xi[10, 7:10] <- fi[14]/12; xgi[10, 7:10] <- gi[14]

## Region 15:
xi[10, 2:6] <- fi[15]/5; xgi[10, 2:6] <- gi[15]

## Region 16:
xi[9, 1:2] <- fi[16]/2; xgi[9, 1:2] <- gi[16]

## Region 17:
xi[10, 1] <- fi[17]; xgi[10, 1] <- gi[17]

cishu <- c(1:93, 95:99, 101:102)
va <- 0; vv <- 0; ji <- 0
for (ci in cishu) {
  set.seed(ci)
  temp.r <- rbinom(100, 1, xgi)
  field.r <- matrix(temp.r, 10, 10, byrow = FALSE)

  
  estTotl <- estpop
  field.xtmp <- field.r*xi
  
  kn <- 2                 # the bolck size is kn*kn
  mn <- 5                 # the number of blocks is mn*mn
  vbar <- estTotl/mn/mn
  figi <- field.xtmp/xgi
  sumtmp <- 0
  
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
    c2bblock1 <- figi[((m-1)*kn+1):(m*kn), 1:kn] ## on its left
    sumtmp <- sumtmp +term1 * (sum(c2bblock1) - vbar)
    c2bblock2 <- figi[((m+1)*kn+1):((m+2)*kn), 1:kn] ## on its right
    sumtmp <- sumtmp +term1 * (sum(c2bblock2) - vbar)
    c2bblock3 <- figi[((m-1)*kn+1):(m*kn), (kn+1):(kn+kn)] ## its left-bottom
    sumtmp <- sumtmp +term1 * (sum(c2bblock3) - vbar)
    c2bblock4 <- figi[(m*kn+1):((m+1)*kn), (kn+1):(kn+kn)] ## under it
    sumtmp <- sumtmp +term1 * (sum(c2bblock4) - vbar)
    c2bblock5 <- figi[((m+1)*kn+1):((m+2)*kn), (kn+1):(kn+kn)] ## its right-bottom
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
    c4bblock1 <- figi[1:kn, ((m-1)*kn+1):(m*kn)] ## above it
    sumtmp <- sumtmp +term1 * (sum(c4bblock1) - vbar)
    c4bblock2 <- figi[1:kn, ((m+1)*kn+1):((m+2)*kn)] ## under it
    sumtmp <- sumtmp +term1 * (sum(c4bblock2) - vbar)
    c4bblock3 <- figi[(kn+1):(kn+kn), ((m-1)*kn+1):(m*kn)] ## its right-top
    sumtmp <- sumtmp +term1 * (sum(c4bblock3) - vbar)
    c4bblock4 <- figi[(kn+1):(kn+kn), (m*kn+1):((m+1)*kn)] ## its right
    sumtmp <- sumtmp +term1 * (sum(c4bblock4) - vbar)
    c4bblock5 <- figi[(kn+1):(kn+kn), ((m+1)*kn+1):((m+2)*kn)] ## its right-bottom
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
  
  NEW.MHT <- sumtmp
  
  ### Pop total's var estimation (end here)
  ##########################################################################
  
  sd <- sqrt(NEW.MHT)
  
  print(paste("The Standard Deviation is:", round(sd*reunit,0)))
  print(paste("The Interval is:", round(estTotl-1.96*sd,5)*reunit,
              round(estTotl+1.96*sd,5)*reunit))
  print(paste("The variation of this interval:", round(sd/estTotl,5)))
  if (sd/estTotl!='NaN'){
    ji <- ji+1
    va[ji] <- sd
    vv[ji] <- sd/estTotl
  }

}
print(paste("The mean Standard Deviation is:", mean(va)))
print(paste("The mean variation of this interval:", mean(vv)))

