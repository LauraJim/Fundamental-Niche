# Laura Jimenez
# Statistical Model for Estimating the Fundamental Niche of a Species
# from presence and physiological data

# Simple case: d=2, it means that we only consider two environmental variables in the analysis

# Functions -----------------------------------------------------------------------------

## Define the species and environmental variables to work with,
# columns Comp1[1] and Comp2[1] in the environmental csv file containing the climate data
# columns Comp1[2] and Comp2[2] in the species csv file containing the climate data
DefineSp <- function(env, data.sp, Comp1, Comp2, duplicate=F)
{
  Comp1 <<- Comp1
  Comp2 <<- Comp2
  # Define the global environmental variables
  env.d <<- as.matrix(env[ , c(Comp1[1],Comp2[1])])
  # Define the species environmental data (occurrences)
  env.sp <<- as.matrix(data.sp[ , c(Comp1[2],Comp2[2])])
  if(duplicate == T) env.sp  <<- subset(env.sp,duplicated(env.sp[,Comp1])==F)  # remove duplicated data
  n <<- dim(env.sp)[1]
  m <<- dim(env.d)[1]
  N <<- n+m
  Et <<- rbind(env.d,env.sp)
}

## In the main model mu is the mean vector and A is the precision matrix
# We define the support (ie. = TRUE if vaild values for mu and Sigma)
# This function is essential for using the t-walk algorithm
# Variable mu.lim must be defined before using this function
Supp <- function(th)
{ # If we have two environmental variables, th has 2 (mu) + 2 (A diag) + 1 (A off diag) = 5 parameters
  mu <<- th[1:2]
  rt <- (mu.lim[1] < mu[1]) & (mu[1] < mu.lim[2])
  rt <- rt & ((mu.lim[3] < mu[2]) & (mu[2] < mu.lim[4]))
  mu <<- matrix( mu, nrow=2, ncol=1)
  if (rt) {     # Test if A positive definite
    A <<- matrix( c( th[3], th[5], th[5], th[4]), nrow=2, ncol=2)
    ev <- eigen(A)$values
    detA <<- prod(ev)
    suma.Et <<- sum( apply( Et, 1, function(yi) { ax<-as.matrix(yi - mu); exp(-0.5 * (t(ax) %*% A %*% ax)) }))
    # TRUE if all eigenvalues are greater than zero (A is positive definite) and suma.Et is positive
    all(ev > 0) & (suma.Et > 0)
  }
  else
    FALSE  # Return value: T or F
}  # Now we have fixed: mu, A and detA

## Produce initial values, simulate them from prior ... here is not done right!!!
# This function is essential for using the t-walk algorithm
Initth <- function()
{
  mu <- mu0 + CholSigma0 %*% rnorm(2)
  
  X1 <- CholW %*% rnorm(2)
  X2 <- CholW %*% rnorm(2)
  
  A <- X1 %*% t(X1) + X2 %*% t(X2)  ### This works for alpha = 2
  
  c( mu[1], mu[2], A[1,1], A[2,2], A[1,2])
  
  #dg <- rgamma(1, 0.1, rate=1)
  #A <- solve(matrix( c(	rgamma(1, alpha, rate=W[1,1]), dg, dg, rgamma(1, alpha, rate=W[2,2])), nrow=2, ncol=2))
  #c( rnorm(1, mean=mu0[1], sd=Sigma0[1,1]), rnorm(1, mean=mu0[2], sd=Sigma0[2,2]), A[1,1], A[2,2], A[1,2])
}

## Plot environmental variables Comp1 and Comp2 of species data env.sp
# If th not NULL, then also plot a Bivariate Normal with parameters th (See function Supp above)
# It's posible to add other valid parameters for the function contour or change the color (See function PlotBN above)
PlotXYEnvVars <- function( th=NULL, col="red", lev=0.95, ...)
{
  plot( env.d, pch=".", col="black", ...)
  points( env.sp, pch=19, col=col)
  if (!(is.null(th)))
    if (Supp(th))
    {
      el<-ellipse::ellipse(x=chol2inv(chol(A)),centre = mu,level=lev)
      lines(el,col="blue",lwd=2)
    }
}

## After running the MCMC this will plot a selection (indices) of the MCMC iterations
# from gives the start iteration and thin the separation between iterations
PlotIterations <- function( info, from=2000, thin=200, lev=0.95,...)
{
  post <- exp(-info$Us - max(-info$Us) + 500)/exp(500) ## This is a normalized posterior from 0 to 1
  post1 <- sort(post)
  
  ix <- which( post1 == max(post))[1] ## looking for the MAP
  th <- info$output[ix,]
  PlotXYEnvVars(th,xlab="Comp1",ylab="Comp2",lev=lev,...)
  
  indices <- seq(from, info$Tr, thin)
  #colo <- colorRamp(col.ran)
  
  for (i in indices)
  {
    if(Supp(info$output[i,])) ##
    {
      el<-ellipse::ellipse(x=chol2inv(chol(A)),centre = mu,level=lev)
      lines(el,col=grey(1-post[i]),lwd=2)
    }
  }
  
  if(Supp(th)) ## after this step we kept the MAP in the variables mu and A
  {
    el<-ellipse::ellipse(x=chol2inv(chol(A)),centre = mu,level=lev)
    lines(el,col="blue",lwd=3)
  }
}

## Function used to save all the values of mu and Sigma simulated from the posterior ----------
save.all <- function(info,from,thin,filename){
  # select only the iterations chosen to plot the ellipses
  indices <- seq(from, info$Tr, thin)
  th <- info$output[indices,]
  save <- matrix(0,nrow=nrow(th),ncol=6)
  for(j in 1:nrow(th))
    # save values in the right order
    save[j,] <- c(th[j,1:3],th[j,5],th[j,5],th[j,4])
  write.csv(save,file=filename,row.names = F)
}

## END