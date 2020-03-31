# Laura Jimenez

# Using the Bayesian Statistical Model for the Estimation of the Fundamental Niche

# Source directory ---------
setwd("My working directory")
source("Nf_Model_functions.R")

# ------------- M A I N --------------------------------------------------------

# Read the environmental variables and occurrence data -----------
# Geo-referenced data (longitud,latitude) + Environmental observations (Comp1,Comp2)

# read background data
envall <- read.csv("Background.csv",header=T)
head(envall)

# read ccurrence data
data <- read.csv("SpD_50_20.csv",header=T)
head(data)
rotule <- "Species name"
# select species color
spcol <- "#FF00AAFF"

# Fixing the species and environmental variables to work with ---------
# Indicate which columns from the occurrence and background files contain the environmental variables
DefineSp( env = envall, data.sp = data, Comp1=c(3,3), Comp2=c(4,4))
# Now we have fixed: Comp1, Comp2, env.d, env.sp, n, m, N, Et

# Number of environmental variables in the study
dd <- 2

# Define a valid interval for mu, depending on the rage of the environmental variables
mu.lim <<- c(min(Et[,1])-0.5,max(Et[,1])+0.5,min(Et[,2])-0.5,max(Et[,2])+0.5)

# A priori tolerance limits ---------
# x-axis=(a1,b1) and y-axis=(a2,b2)
a1 <- 0 #xleft
b1 <- 1.2 #xright
a2 <- -1 #ybottom
b2 <- 0.8 #ytop

# First Plots ----------------------
### Plot the geographical and environmental spaces with the occurrence of the species on top
x11()
par(mfrow=c(1,2))
## Geographical Space:
# Plot the geographical locations of reported presences of the species, on top of the global environmental variables
plot(envall[,1],envall[,2],pch=".",col=1,xlab="Longitude",ylab="Latitude",main="Geographical space")
points(data[,1],data[,2],pch=19,col=spcol)

## Environmental Space:
# Plot environmental variables of species data and the location of reported presences of the species on top
plot(env.d, pch=".", col=1,xlab="Bio1WHStnd",ylab="Bio12WHStnd",main="Environmental space")
points(env.sp, pch=19, col=spcol)
rect(xleft=a1,xright=b1,ybottom=a2,ytop=b2,border="gold",lwd=2)
legend("topleft",legend=c("Species presences:",rotule,"Tolerance ranges"),pch=c(19,NA,NA),
       bty="n",lty=c(0,0,1),col=c(spcol,"white","gold"),lwd=2)

### Statistical Model ----------------

# Define the prior parameters ----------
# Get estimations for the sample mean and variance of every environmental variable
s1.prior <- ((b1-a1)/6)^2
s2.prior <- ((b2-a2)/6)^2
mu1.prior <- a1 + (b1-a1)/2
mu2.prior <- a2 + (b2-a2)/2

# For the centroid, this is a bivariate normal distribution with:
mu0 <- c(mu1.prior,mu2.prior)     #mean vector mu0
Sigma0 <- matrix(c(s1.prior,0,0,s2.prior),nrow=2,ncol=2) #variance covariance matrix Sigma0
CholSigma0 <- chol(Sigma0)
A0 <- chol2inv(CholSigma0)    # precision matrix
A0mu0 <- A0 %*% mu0
mu0
Sigma0
# Plot the contour line in the environmental space:
el<-ellipse::ellipse(x=Sigma0,centre = mu0,level=0.95)
lines(el,col="gold",lwd=2)

# For the precision matrix, it is a Wishart distribution with
alpha <- 2 #shape parameter, default, do not move for now
W <- alpha*A0  #scale matrix W.
CholW <- chol(W)
Winv <- chol2inv(CholW)
W
# The prior expected value for the precision matrix is E(A) = W


# Energy = - log ( posterior ) ----------
wi <<- rep(1,n)
# This is called right after Supp: mu, A and detA are already defined, th is ignored
Energy <- function(th)
{ # This is called right after Supp: mu, A and detA are already defined, th is ignored
  ax1 <- (mu - mu0)
  ax2 <- apply( env.sp, 1, function(xi) { ax<-as.matrix(xi - mu); t(ax) %*% A %*% ax })
  # first two terms are generic in the posterior due to the normal model
  S <- 0.5*sum(wi*ax2) + n*log(suma.Et)
  # these terms correspond to the priors:
  S <- S + 0.5*( t(ax1) %*% A0 %*% ax1 + sum(diag(A %*% Winv)) - (alpha-dd-1)*log(detA) )
  
  S # + 10**6 numerical artifact
}

# Run the MCMC to produce simulations from the posterior ---------
library(Rtwalk)
Run <- function(Tr=20000)
{
  info <- Runtwalk( Tr=Tr, dim=5, Obj=Energy, Supp=Supp, x0=Initth(), xp0=Initth())
  x11()
  PlotIterations( info )
  
  info
}

ptm <- proc.time()
info <- Run(Tr=7000)
proc.time() - ptm

# Plotting results:
mu
chol2inv(chol(A))
# save al the values of mu and A

x11()
PlotIterations(info,col=spcol,main=paste(rotule,"7K",sep="_"))
# a priori ellipse
el<-ellipse::ellipse(x=Sigma0,centre = mu0,level=0.95)
lines(el,col="gold",lwd=2)

fname <- paste0(paste(rotule,"7K","output",sep="_"),".csv")
save.all(info,2000,200,paste0(".\\My folder\\",fname))

### END
