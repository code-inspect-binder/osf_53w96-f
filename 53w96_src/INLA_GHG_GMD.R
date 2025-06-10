library(INLA)
library(ncdf4)
setwd("/PATH/TO/FOLDER/CONTAINING/DATA") #Update this as needed
set.seed(2)
#
# Example script using simulated data to estimate
# UK methane emissions (without boundary emissions)
# using the R-INLA package
#

readncdf <- function(ncname, variables){
  #' Read data from netcdf
  #' Inputs:
  #' ncname (string): name of netcdf file (without .nc)
  #' variables (list): list of variable names
  #' Output
  #' list of variable names and values
  ncfname <- paste(ncname, ".nc", sep = "")
  ncin <- nc_open(ncfname)
  output <- lapply(variables,  function(var){ncvar_get(ncin, var)})
  nc_close(ncin)
  names(output) <- variables
  return(output)
}
precision.ar1 = function(N, rho){
  #' Generate AR1 precision matrix x_t = rho*x_t-1 + sigma
  #' Input:
  #' N (int): Size if matrix
  #' rho (double): correlation 
  #' Output:
  #' Matrix containing AR1 precision matrix
  
  Q = diag(1+rho^2, N)
  for (i in 1:(N-1)) {
    Q[i, i+1] = -rho
    Q[i+1, i] = -rho
  }
  Q[1,1] = 1
  Q[N,N] = 1
  return(as(Q, 'dgCMatrix') )
}
llarea <-function(lon, lat)
  #' Calculate area of gridsquare from lat/lon
  #' Input:
  #' lon (vector): Longitudes
  #' lat (vector): Latitudes
  #' Output:
  #' area of grid squares
{
  re <- 6367500.0 #radius of Earth
  pi <- 3.141592653589793 #Pi
  dlon <- abs(mean(lon[-1] - lon[-length(lon)]))*pi/180.
  dlat <- abs(mean(lat[-1] - lat[-length(lat)]))*pi/180.
  theta <- pi*(90.-lat)/180.
  
  area <- matrix(nrow = length(lon) , ncol = length(lat))
  for(latI in 1:length(lat)){
    if (theta[latI] == 0 | theta[latI] == pi) {
      area[,latI] = (re^2)*abs(cos(dlat/2.) - 1)*dlon
    } else {
      lat1 = theta[latI] - dlat/2.
      lat2 = theta[latI] + dlat/2.
      area[,latI] = (re^2)*(cos(lat1) - cos(lat2))*dlon
    }
  }
  (area)
}

#############################################
#####BEGIN INLA/SPDE FOR GHG emissions ######
#############################################
print(paste0("Start at ", Sys.time()))


#############################################
############## READ IN DATA #################
#############################################

# #Get locations for UK and its lat and lons
EURvars <- readncdf("country_EUROPE", list("lon", "lat", "country"))
lon <- EURvars$lon
lat <- EURvars$lat
UKinds <- which((EURvars$country == 101) | (EURvars$country == 49), arr.ind = TRUE)

data <- readRDS("INLA_GHG_GMD_data.rds") #Read in footprint and prior mean data
grbf <- readRDS("INLA_GHG_GMD_grbf.rds") #Read in spatial domain data
nt <- length(data) #Number of time steps

print(paste0("Read data at ", Sys.time()))
#############################################
############## SET UP SPDE ##################
#############################################

#Get prior emission locations for mesh
llgrid <- expand.grid(lon, lat)
sdllgrid<- llgrid[grbf==0,]

# Get sub domain for mesh
subgridllix <- which(grbf == 0, arr.ind = TRUE)
subgridllon <- lon[min(subgridllix[,1])]
subgridulon <- lon[max(subgridllix[,1])]
subgridllat <- lat[min(subgridllix[,2])]
subgridulat <- lat[max(subgridllix[,2])]

#Make a mesh using UK as nodes
locs = cbind(lon[UKinds[seq(1,length(UKinds[,1]), 2),1]], lat[UKinds[seq(1,length(UKinds[,1]), 2),2]])
inner <- cbind(c(subgridllon,subgridulon,subgridulon,subgridllon), 
               c(subgridllat,subgridllat,subgridulat,subgridulat))
outer <- cbind(c(subgridllon-10,subgridulon+10,subgridulon+10,subgridllon-10), 
               c(subgridllat-10,subgridllat-10,subgridulat+10,subgridulat+10))
mesh <- inla.mesh.2d(loc=locs, boundary = list(inner, outer) , max.edge = c(1.75, 8), cutoff=c(0.25,2))
A = inla.spde.make.A(mesh, loc = as.matrix(sdllgrid))
nx <- nt*dim(A)[2]

#Operator on x, rescale to be in y unit
Hlist <- list()
for(i in seq(nt)){
  Hlist[[paste0("H",i)]] <- with(data, t(apply(data[[i]]$sdfp, 1, function(x) x*data[[i]]$prior)) %*% A)
}
H<-bdiag(unlist(Hlist))
strl <- function(x){
  x[x < max(x)/1000] <- 0
  return(x)
}
H <- as(t(apply(H, 1, strl)), 'dgCMatrix')  #Make sparser  

#Make pseudo-data
sigma0 <- 0.5
range0 <- 3.25
kappa0 <- sqrt(8)/range0
tau0 <- 1/(sqrt(4 * pi) * kappa0 * sigma0)
meshtrue <- inla.mesh.2d(loc=sdllgrid, boundary = list(inner), max.edge = 1)
At = inla.spde.make.A(meshtrue, loc = as.matrix(sdllgrid))
model.sim <- inla.spde2.matern(meshtrue,   
                               B.tau = matrix(log(tau0)),
                               B.kappa = matrix(log(kappa0)))
Qs0=inla.spde.precision(model.sim)
## Time covariance 
rho0 <- 0.8
Qt0 <- precision.ar1(nt, rho0)
Q0 <- kronecker(Qt0,Qs0)
pPrior <- inla.qsample(n=1,Q0,seed=2)
nxp <- length(pPrior)
Htlist <- list()
for(i in seq(nt)){
  Htlist[[paste0("H",i)]] <- with(data, t(apply(data[[i]]$sdfp, 1, function(x) x*data[[i]]$prior)) %*% At)
}
Ht<-bdiag(unlist(Htlist))
y <- as.vector(Ht %*% as.vector(pPrior))  #Prior measurement
y <- y + rnorm(length(y), 0, abs(rowSums(H))*0.15)

## prior Matern parameters
sigma0 <- 0.1
range0 <- 5
spde <- inla.spde2.pcmatern(
  mesh=mesh, alpha=2, ### mesh and smoothness parameter
  prior.range=c(range0, 0.5), ### P(practic.range<range0)=0.5
  prior.sigma=c(sigma0, 0.01))

#Generate a list 'mesh.index'
mesh.index = inla.spde.make.index(name = "field", n.spde = spde$n.spde, n.repl = nt)
for(i in seq(nt)){
  mesh.index$field.group[((i-1)*spde$n.spde+1):(i*spde$n.spde)] = i
}


#Collect predictor information for the observed data                
st.est = inla.stack(data = list(y = y),   #Observations
                    A = list(H),   #This should be my H matrix
                    effects = list(c(mesh.index)), #Whether they're fixed or random (include everything)
                    tag = "est")
#Linear predictor for UK
LLarea <- llarea(lon,lat) 
sdcs <- EURvars$country[grbf==0]
areas <- LLarea[grbf==0]
areas[-which(sdcs == 101, arr.ind = TRUE)] <- 0
npr <- (data$data1$prior*areas*1e-9*(365.*24.*3600.)*16.04*1e-12)
Anpr <- bdiag(npr %*%A,npr %*%A,npr %*%A,npr %*%A) 
st.pred = inla.stack(data=list(UKtot=NA), A=list(Anpr), effects=list(c(mesh.index)), tag="pred")

stack = inla.stack(st.est, st.pred)


##Formula for inference:
formula <- y ~ -1 + f(field, model=spde, 
                      control.group = list(model = "ar1",
                                           hyper=list(rho=list(prior='betacorrelation',param=c(6.5,0.1), initial=-4))),
                      group=field.group)

rm(data)
print(paste0("SPDE set up at ", Sys.time()))

# Inference using INLA
print(paste0("INLA starting at ", Sys.time()))
sigypdf <- c(-5,1)  #Hyperprior for precision likelihood
resi = inla(formula,
            data = inla.stack.data(stack, spde = spde),
            family = "normal",
            control.family = list(hyper=list(theta=list(prior='normal', param=sigypdf))),
            control.predictor=list(A=inla.stack.A(stack), compute=TRUE))



print(paste0("Inference completed at ", Sys.time()))


#Get UK totals
mfx <- resi$summary.linear.predictor$mode[(length(y)+1):(length(y)+nt)]
lower <- resi$summary.linear.predictor[['0.025quant']][(length(y)+1):(length(y)+nt)]
upper <-  resi$summary.linear.predictor[['0.975quant']][(length(y)+1):(length(y)+nt)]
sd <- resi$summary.linear.predictor$sd[(length(y)+1):(length(y)+nt)]

true <- rep(NA,nt)
estimated <- rep(NA,nt)
for(tt in seq(nt)){
  true[tt] <- as.double( npr %*%At %*% pPrior[(1+(tt-1)*(length(pPrior)/nt)):(tt*length(pPrior)/nt)] )
}
par(pty="s")
plot(true, mfx, ylim=c(min(lower),max(upper)),xlim=c(min(lower),max(upper) ), 
     pch=4, xlab='True deviation (Tg/yr)', ylab='Estimated deviation (Tg/yr)', asp=1)
lines(c(-10,10),c(-10,10), col="red")
arrows(true, lower, true, upper, length=0.05, angle=90, code=3)

