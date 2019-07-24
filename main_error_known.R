
# This program is designed to reproduce table 1 in the manuscript
# i.e., to analyze the four different types of testing outcomes
# which are individual test, Master pool test, Dorfman test and Array test
# Here we assume that the testing error is known
#############################################################################################




library(mvtnorm)
library(msm)
library(ltsa)
library(geoR)
library(Matrix)
library(coda)

source("BNR_GP.txt")
source("BNR_GP_FUN.txt")
source("BNR_PP.txt")
source("BNR_PP_FUN.txt")
source("Testing functions.txt")
Rcpp::sourceCpp("SampLatent.cpp")
Rcpp::sourceCpp("ErrorUpdate.cpp")


##########################################################################################
# Specify the number of Gibbs iterates and prior parameters
###################
# Simulation settings
N<-5000                  # Number of Individuals
n.keep <- 1000           # number of posterior sample
n.burn <- 2000           # burn in 
mod    <- 2              # Choose Data Generating Model (1 for model 1; 2 for model 2)
n.step <- 50             # step to plot
n.thin <- 5              # thinning parameter
Se.t   <-c(0.95,0.98)    # True Sensitivity
Sp.t   <-c(0.98,0.99)    # True Specificity
cj     <-5               # Pool size, note array testing uses cjxcj arrays
kap    <- 2.0            # smooth parameter in matern function
trphi_tune <- 0.1        # tune parameter for decay parameter phi
Q      <- 2              # number of unknown functions
L_1    <- 100            # number of knots for Gaussian predictive process associated with function 1 (GPP_1)
L_2    <- 100            # number of knots for Gaussian predictive process associated with function 2 (GPP_2)
lbd_1  <- -3;rbd_1 <- 3  # left and right endpoints for GPP_1
lbd_2  <- -3;rbd_2 <- 3  # left and right endpoints for GPP_2

##########################################################################################
# Define the four prevalence functions
# X range [-3,3]

aepd <- function(x,k=1.2,sig=2.5,a=2){
  0.7*exp((-(x>0)*k^a-(x<0)/k^a)*abs(x/sig)^a)
}

norm_mix <- function(x,mu1=-1.5,mu2=1.5,sig1=.6,sig2=.8,wt=.6){
  wt*exp(-(x-mu1)^2/(2*sig1^2))+(1-wt)*exp(-(x-mu2)^2/(2*sig2^2))
}

damp_sin<- function(x){0.2*(sin(pi*(x+0.3)/2.5)+2)*(exp(-(x+(x-0.3)^2*(x>0.3))/6))}

f4 <- function(x,a=1,b=1.5,c=4){
  c*exp(a+b*x)/(6+6*exp(a+b*x))
}

##########################################################################################
# Start the simulation
lrow   <- 7248          # Length of row; will use it to fill with NA in case error happened
nIter  <- 500           # Number of simulated dataset


sim    <- 1

################################################################################################
### <---------- Start simulation ---------- > ###
while(sim <=nIter){
	
	#########################################################################
	# Generate 2 continuous covariates that are nonlinearly related to logit(p)
	# e.g., age(rounded) and biomarker level
	x_1 <- round(runif(N,lbd_1,rbd_1),2)
	x_2 <- round(runif(N,lbd_2,rbd_2),2)
	
	# Generate 2 covariates in linear predictor
	x_3 <- rnorm(N) 
	x_4 <- rbinom(N,1,0.5)
	V   <- cbind(1,x_3,x_4)
	VtV <- crossprod(V,V)
	#########################################################################
	# To form the pools, use random pooling
	# True value of intercept, regression coefficients for x3 and x4
	rbeta  <- c(-1.8,0.5,0.5)
	id_eta <- rep(NA,N)
	#########################################################################
	# Generate true statuses of each individual under specified models (mod=1 or 2)
	# Combine function 1 and 2
	######################################################################### model 1
	if(mod == 1){ # combine of function aepd and mix normal
	  int_1 <- integrate(aepd,lbd_1,rbd_1)$value/(rbd_1-lbd_1)
	  int_2 <- integrate(norm_mix,lbd_2,rbd_2)$value/(rbd_2-lbd_2)
	  fun_1 <- aepd(x_1) - int_1
	  fun_2 <- norm_mix(x_2) - int_2
	  p_x = pnorm(fun_1+fun_2+V%*%rbeta)
	  Y_true<-rbinom(N,1,p_x)
	  
	  dat_1 = data.frame(table(x_1))
	  mcx_1 = as.numeric(as.character(factor(dat_1[,1], levels = unique(dat_1[,1]))))
	  d_1 = as.numeric(dat_1[,2])
	  K_1 = length(d_1)
	  D_1 = Diagonal(x=d_1)
	  for(i in 1:K_1){
		id_eta[which(x_1==mcx_1[i])] = i
	  }
	  M_mat_1 = Matrix(0,N,K_1)
	  M_mat_1[cbind(1:N,id_eta)] = 1
	  ###################################################### configure for function 2
	  dat_2 = data.frame(table(x_2))
	  mcx_2 = as.numeric(as.character(factor(dat_2[,1], levels = unique(dat_2[,1]))))
	  d_2 = as.numeric(dat_2[,2])
	  K_2 = length(d_2)
	  D_2 = Diagonal(x=d_2)
	  for(i in 1:K_2){
		id_eta[which(x_2==mcx_2[i])] = i
	  }
	  M_mat_2 = Matrix(0,N,K_2)
	  M_mat_2[cbind(1:N,id_eta)] = 1
	}

	########################################################################## model 2
	if(mod == 2){ # combine of function damp_sin and f4
	  int_1 <- integrate(damp_sin,lbd_1,rbd_1)$value/(rbd_1-lbd_1)
	  int_2 <- integrate(f4,lbd_2,rbd_2)$value/(rbd_2-lbd_2)
	  fun_1 <- damp_sin(x_1) - int_1
	  fun_2 <- f4(x_2) - int_2
	  p_x = pnorm(fun_1+fun_2+V%*%rbeta)
	  Y_true<-rbinom(N,1,p_x)
	  
	  dat_1 = data.frame(table(x_1))
	  mcx_1 = as.numeric(as.character(factor(dat_1[,1], levels = unique(dat_1[,1]))))
	  d_1 = as.numeric(dat_1[,2])
	  K_1 = length(d_1)
	  D_1 = Diagonal(x=d_1)
	  for(i in 1:K_1){
		id_eta[which(x_1==mcx_1[i])] = i
	  }
	  M_mat_1 = Matrix(0,N,K_1)
	  M_mat_1[cbind(1:N,id_eta)] = 1
	  ###################################################### configure for function 2
	  dat_2 = data.frame(table(x_2))
	  mcx_2 = as.numeric(as.character(factor(dat_2[,1], levels = unique(dat_2[,1]))))
	  d_2 = as.numeric(dat_2[,2])
	  K_2 = length(d_2)
	  D_2 = Diagonal(x=d_2)
	  for(i in 1:K_2){
		id_eta[which(x_2==mcx_2[i])] = i
	  }
	  M_mat_2 = Matrix(0,N,K_2)
	  M_mat_2[cbind(1:N,id_eta)] = 1
	}

	
	############################################################# Configure setting for function 1
	mcw_1 = seq(lbd_1,rbd_1,length.out=L_1)   # Knots for Predictive Process 1
	x_dist_1 = as.matrix(dist(mcx_1,mcx_1, method = "manhattan", diag = FALSE, upper = FALSE))
	knots_dist_1 = as.matrix(dist(mcw_1,mcw_1, method = "manhattan", diag = FALSE, upper = FALSE))
	cross_dist_1 = matrix(NA,K_1,L_1)
	for(j in 1:L_1){
	  cross_dist_1[,j] <- abs(mcx_1-mcw_1[j])
	}
	id0_mcw_1 <- which.min(abs(mcw_1))
	id0_mcx_1 <- which(abs(mcx_1)<0.001)
	if(length(id0_mcx_1)==0){id0_mcx_1<-0}

	############################################################# Configure setting for function 2
	mcw_2 = seq(lbd_2,rbd_2,length.out=L_2)   # Knots for Predictive Process 2
	x_dist_2 = as.matrix(dist(mcx_2,mcx_2, method = "manhattan", diag = FALSE, upper = FALSE))
	knots_dist_2 = as.matrix(dist(mcw_2,mcw_2, method = "manhattan", diag = FALSE, upper = FALSE))
	cross_dist_2 = matrix(NA,K_2,L_2)
	for(j in 1:L_2){
	  cross_dist_2[,j] <- abs(mcx_2-mcw_2[j])
	}
	id0_mcw_2 <- which.min(abs(mcw_2))
	id0_mcx_2 <- which(abs(mcx_2)<0.001)
	if(length(id0_mcx_2)==0){id0_mcx_2<-0}

	#############################################################
	### determine the range of decay parameter phi for function 1
	x_length_1 = rbd_1-lbd_1
	phi_range = seq(0.01,100,by=0.01)
	lw = rep(NA,length(phi_range))
	up = lw
	for(i in 1:length(phi_range)){
	  lw[i] = matern(x_length_1*2/30,phi=phi_range[i],kappa = kap)
	  up[i] = matern(x_length_1*2/3,phi=phi_range[i],kappa = kap)
	}
	phi_min_1 = phi_range[which.min(abs(lw-0.05))]
	phi_max_1 = phi_range[which.min(abs(up-0.05))]
	#############################################################
	### determine the range of decay parameter phi for function 2
	x_length_2 = rbd_2-lbd_2
	phi_range = seq(0.01,100,by=0.01)
	lw = rep(NA,length(phi_range))
	up = lw
	for(i in 1:length(phi_range)){
	  lw[i] = matern(x_length_2*2/30,phi=phi_range[i],kappa = kap)
	  up[i] = matern(x_length_2*2/3,phi=phi_range[i],kappa = kap)
	}
	phi_min_2 = phi_range[which.min(abs(lw-0.05))]
	phi_max_2 = phi_range[which.min(abs(up-0.05))]
	################################################################################## MCMC configuration
	phi_ini <- matrix(NA,Q,2)
	phi_ini[1,] <- c(phi_min_1,phi_max_1)
	phi_ini[2,] <- c(phi_min_2,phi_max_2)
	sample.parameter <- c(n.burn,n.keep,n.thin,n.step)
	
	#=================================================================================
	# Start data simulation 
	#
	#=================================================================================
	
	################################################################################## 
	# Individual Test
	#
	##################################################################################

	###### Simulates Individual testing
	test.res <- Pool.test(Y_true,Se.t[2],Sp.t[2],cj=1)
	Z.it     <- test.res$Z
	Y.it     <- test.res$Y
	Y.it_obs <- Z.it[,1]
	
	#################################################################################
	# Fit Individual Test data using Gaussian Predictive Process
	#
	# ===============================================================================
	t1 <- proc.time()
	res2<-Bayes.PP(Z=Z.it,Y=Y.it,Y.it_obs,D_1,mcx_1,mcw_1,
				   D_2,mcx_2,mcw_2,phi_ini,
				   trphi_tune,sample.parameter,kap,na=1, 
				   Se=Se.t[2],Sp=Sp.t[2],
				   est.error=FALSE,visual=FALSE)
	t2 <- proc.time()
	
	########################################## save results
	if(class(res2)!='try-error'){
		Theta <- rbind(res2$beta,res2$tau_1,res2$tau_2,res2$phi_1,res2$phi_2)
		Theta.mean <- apply(Theta,1,mean)
		Theta.med <- apply(Theta,1,median)
		Theta.sd <- apply(Theta,1,sd)
		Theta.bd <- HPDinterval(as.mcmc(t(Theta)))

		Eta_1.new <- res2$X_new_1
		Eta_1 <- res2$eta_new_1
		Eta_1.mean <- apply(Eta_1,1,mean)
		Eta_1.med <- apply(Eta_1,1,median)
		Eta_1.sd <- apply(Eta_1,1,sd)
		Eta_1.bd <- HPDinterval(as.mcmc(t(Eta_1)))

		Eta_2.new <- res2$X_new_2
		Eta_2 <- res2$eta_new_2
		Eta_2.mean <- apply(Eta_2,1,mean)
		Eta_2.med <- apply(Eta_2,1,median)
		Eta_2.sd <- apply(Eta_2,1,sd)
		Eta_2.bd <- HPDinterval(as.mcmc(t(Eta_2)))

		IT_PP <- c(as.numeric(t2-t1)[3],Theta.mean,Theta.med,Theta.sd,Theta.bd,
				Eta_1.new,Eta_1.mean,Eta_1.med,Eta_1.sd,Eta_1.bd,
				Eta_2.new,Eta_2.mean,Eta_2.med,Eta_2.sd,Eta_2.bd)}else{IT_PP<-rep(NA,lrow)}
				
	#################################################################################
	# Fit Individual Test data using Regular Gaussian Process
	#
	# ===============================================================================

	t1 <- proc.time()
	res1<-try(Bayes.GP(Z=Z.it,Y=Y.it,Y.it_obs,D_1,mcx_1,D_2,mcx_2,
				   phi_ini,trphi_tune,sample.parameter,kap,na=1,
				   Se=Se.t[2],Sp=Sp.t[2],
				   est.error=FALSE,visual=FALSE),silent=TRUE)
	t2 <- proc.time()
	
	########################################## save results
	if(class(res1)!='try-error'){
		Theta <- rbind(res1$beta,res1$tau_1,res1$tau_2,res1$phi_1,res1$phi_2)
		Theta.mean <- apply(Theta,1,mean)
		Theta.med <- apply(Theta,1,median)
		Theta.sd <- apply(Theta,1,sd)
		Theta.bd <- HPDinterval(as.mcmc(t(Theta)))

		Eta_1.new <- res1$X_new_1
		Eta_1 <- res1$eta_new_1
		Eta_1.mean <- apply(Eta_1,1,mean)
		Eta_1.med <- apply(Eta_1,1,median)
		Eta_1.sd <- apply(Eta_1,1,sd)
		Eta_1.bd <- HPDinterval(as.mcmc(t(Eta_1)))

		Eta_2.new <- res1$X_new_2
		Eta_2 <- res1$eta_new_2
		Eta_2.mean <- apply(Eta_2,1,mean)
		Eta_2.med <- apply(Eta_2,1,median)
		Eta_2.sd <- apply(Eta_2,1,sd)
		Eta_2.bd <- HPDinterval(as.mcmc(t(Eta_2)))
		IT_GP <- c(as.numeric(t2-t1)[3],Theta.mean,Theta.med,Theta.sd,Theta.bd,
			Eta_1.new,Eta_1.mean,Eta_1.med,Eta_1.sd,Eta_1.bd,
			Eta_2.new,Eta_2.mean,Eta_2.med,Eta_2.sd,Eta_2.bd)}else{IT_GP<-rep(NA,lrow)}

	################################################################################## 
	# Master Pool Test
	# Data generated from individual test data
	##################################################################################
	
	# Simulates Master Pool testing
	test.res<-Pool.test(Y_true,Se.t[1],Sp.t[1],cj=cj)
	Z.pool<-test.res$Z
	Y.pool<-test.res$Y
	Y.master_obs <- rep(NA,N)
	for(i in 1:nrow(Z.pool)){
	  Y.master_obs[seq((i-1)*cj+1,i*cj)] <- Z.pool[i,1]
	}
	#################################################################################
	# Fit Master Pool Test data using Gaussian Predictive Process
	#
	# ===============================================================================
	t1 <- proc.time()
	res2<-try(Bayes.PP(Z=Z.pool,Y=Y.pool,Y.master_obs,D_1,mcx_1,mcw_1,
				   D_2,mcx_2,mcw_2,phi_ini,
				   trphi_tune,sample.parameter,kap,na=1, 
				   Se=Se.t[1],Sp=Sp.t[1],
				   est.error=FALSE,visual=FALSE),silent=TRUE)
	t2 <- proc.time()
	
	########################################## Save results
	if(class(res2)!='try-error'){
		Theta <- rbind(res2$beta,res2$tau_1,res2$tau_2,res2$phi_1,res2$phi_2)
		Theta.mean <- apply(Theta,1,mean)
		Theta.med <- apply(Theta,1,median)
		Theta.sd <- apply(Theta,1,sd)
		Theta.bd <- HPDinterval(as.mcmc(t(Theta)))

		Eta_1.new <- res2$X_new_1
		Eta_1 <- res2$eta_new_1
		Eta_1.mean <- apply(Eta_1,1,mean)
		Eta_1.med <- apply(Eta_1,1,median)
		Eta_1.sd <- apply(Eta_1,1,sd)
		Eta_1.bd <- HPDinterval(as.mcmc(t(Eta_1)))

		Eta_2.new <- res2$X_new_2
		Eta_2 <- res2$eta_new_2
		Eta_2.mean <- apply(Eta_2,1,mean)
		Eta_2.med <- apply(Eta_2,1,median)
		Eta_2.sd <- apply(Eta_2,1,sd)
		Eta_2.bd <- HPDinterval(as.mcmc(t(Eta_2)))
		MPT_PP <- c(as.numeric(t2-t1)[3],Theta.mean,Theta.med,Theta.sd,Theta.bd,
			Eta_1.new,Eta_1.mean,Eta_1.med,Eta_1.sd,Eta_1.bd,
			Eta_2.new,Eta_2.mean,Eta_2.med,Eta_2.sd,Eta_2.bd)}else{MPT_PP<-rep(NA,lrow)}

	#################################################################################
	# Fit Master Pool Test data using Regular Gaussian Process
	#
	# ===============================================================================

	t1 <- proc.time()
	res1<-try(Bayes.GP(Z=Z.pool,Y=Y.pool,Y.master_obs,D_1,mcx_1,D_2,mcx_2,
				   phi_ini,trphi_tune,sample.parameter,kap,na=1,
				   Se=Se.t[1],Sp=Sp.t[1],
				   est.error=FALSE,visual=FALSE),silent=TRUE)
	t2 <- proc.time()

	########################################## Save results
	if(class(res1)!='try-error'){
		Theta <- rbind(res1$beta,res1$tau_1,res1$tau_2,res1$phi_1,res1$phi_2)
		Theta.mean <- apply(Theta,1,mean)
		Theta.med <- apply(Theta,1,median)
		Theta.sd <- apply(Theta,1,sd)
		Theta.bd <- HPDinterval(as.mcmc(t(Theta)))

		Eta_1.new <- res1$X_new_1
		Eta_1 <- res1$eta_new_1
		Eta_1.mean <- apply(Eta_1,1,mean)
		Eta_1.med <- apply(Eta_1,1,median)
		Eta_1.sd <- apply(Eta_1,1,sd)
		Eta_1.bd <- HPDinterval(as.mcmc(t(Eta_1)))

		Eta_2.new <- res1$X_new_2
		Eta_2 <- res1$eta_new_2
		Eta_2.mean <- apply(Eta_2,1,mean)
		Eta_2.med <- apply(Eta_2,1,median)
		Eta_2.sd <- apply(Eta_2,1,sd)
		Eta_2.bd <- HPDinterval(as.mcmc(t(Eta_2)))
		MPT_GP<-c(as.numeric(t2-t1)[3],Theta.mean,Theta.med,Theta.sd,Theta.bd,
			Eta_1.new,Eta_1.mean,Eta_1.med,Eta_1.sd,Eta_1.bd,
			Eta_2.new,Eta_2.mean,Eta_2.med,Eta_2.sd,Eta_2.bd)}else{MPT_GP<-rep(NA,lrow)}
			
	################################################################################## 
	# Dorfman Test
	# Data generated from individual test data
	##################################################################################
	

	# Simulates Dorfman testing
	test.res<-Dorfman.decode.diff.error(Y_true,Se.t,Sp.t,cj=cj)
	Z.dorf<-test.res$Z
	Y.dorf<-test.res$Y
	Y.dorf_obs <- rep(NA,N)
	Y.dorf_obs[Y.dorf[,2]==1] <- 0
	id_2 <- (Y.dorf[,2]==2)
	Y.dorf_obs[id_2] <- Z.dorf[Y.dorf[id_2,4]]
	
	#################################################################################
	# Fit Dorfman Test data using Gaussian Predictive Process
	#
	# ===============================================================================

	t1 <- proc.time()
	res2<-try(Bayes.PP(Z=Z.dorf,Y=Y.dorf,Y.dorf_obs,D_1,mcx_1,mcw_1,
				   D_2,mcx_2,mcw_2,phi_ini,
				   trphi_tune,sample.parameter,kap,na=length(Se.t), 
				   Se=Se.t,Sp=Sp.t,
				   est.error=FALSE,visual=FALSE),silent=TRUE)
	t2 <- proc.time()
	########################################## save results
	if(class(res2)!='try-error'){
		Theta <- rbind(res2$beta,res2$tau_1,res2$tau_2,res2$phi_1,res2$phi_2)
		Theta.mean <- apply(Theta,1,mean)
		Theta.med <- apply(Theta,1,median)
		Theta.sd <- apply(Theta,1,sd)
		Theta.bd <- HPDinterval(as.mcmc(t(Theta)))

		Eta_1.new <- res2$X_new_1
		Eta_1 <- res2$eta_new_1
		Eta_1.mean <- apply(Eta_1,1,mean)
		Eta_1.med <- apply(Eta_1,1,median)
		Eta_1.sd <- apply(Eta_1,1,sd)
		Eta_1.bd <- HPDinterval(as.mcmc(t(Eta_1)))

		Eta_2.new <- res2$X_new_2
		Eta_2 <- res2$eta_new_2
		Eta_2.mean <- apply(Eta_2,1,mean)
		Eta_2.med <- apply(Eta_2,1,median)
		Eta_2.sd <- apply(Eta_2,1,sd)
		Eta_2.bd <- HPDinterval(as.mcmc(t(Eta_2)))
		DT_PP<-c(as.numeric(t2-t1)[3],Theta.mean,Theta.med,Theta.sd,Theta.bd,
			Eta_1.new,Eta_1.mean,Eta_1.med,Eta_1.sd,Eta_1.bd,
			Eta_2.new,Eta_2.mean,Eta_2.med,Eta_2.sd,Eta_2.bd)}else{DT_PP<-rep(NA,lrow)}

	#################################################################################
	# Fit Dorfman Test data using Regular Gaussian Process
	#
	# ===============================================================================

	t1 <- proc.time()
	res1<-try(Bayes.GP(Z=Z.dorf,Y=Y.dorf,Y.dorf_obs,D_1,mcx_1,D_2,mcx_2,
				   phi_ini,trphi_tune,sample.parameter,kap,na=length(Se.t),
				   Se=Se.t,Sp=Sp.t,
				   est.error=FALSE,visual=FALSE),silent=TRUE)
	t2 <- proc.time()

	if(class(res1)!='try-error'){
		Theta <- rbind(res1$beta,res1$tau_1,res1$tau_2,res1$phi_1,res1$phi_2)
		Theta.mean <- apply(Theta,1,mean)
		Theta.med <- apply(Theta,1,median)
		Theta.sd <- apply(Theta,1,sd)
		Theta.bd <- HPDinterval(as.mcmc(t(Theta)))

		Eta_1.new <- res1$X_new_1
		Eta_1 <- res1$eta_new_1
		Eta_1.mean <- apply(Eta_1,1,mean)
		Eta_1.med <- apply(Eta_1,1,median)
		Eta_1.sd <- apply(Eta_1,1,sd)
		Eta_1.bd <- HPDinterval(as.mcmc(t(Eta_1)))

		Eta_2.new <- res1$X_new_2
		Eta_2 <- res1$eta_new_2
		Eta_2.mean <- apply(Eta_2,1,mean)
		Eta_2.med <- apply(Eta_2,1,median)
		Eta_2.sd <- apply(Eta_2,1,sd)
		Eta_2.bd <- HPDinterval(as.mcmc(t(Eta_2)))
		DT_GP<-c(as.numeric(t2-t1)[3],Theta.mean,Theta.med,Theta.sd,Theta.bd,
			Eta_1.new,Eta_1.mean,Eta_1.med,Eta_1.sd,Eta_1.bd,
			Eta_2.new,Eta_2.mean,Eta_2.med,Eta_2.sd,Eta_2.bd)}else{DT_GP<-rep(NA,lrow)}
		

	################################################################################## 
	# Array Test
	# Data generated from individual test data
	##################################################################################

	# Simulates array testing  
	test.res<-Array.decode.diff.error(Y_true, Se.t, Sp.t, cj=cj)
	Z.array<-test.res$Z
	Y.array<-test.res$Y
	Y.array_obs <- rep(NA,N)
	Y.array_obs[Y.array[,2]==2] <- 0
	id_2 <- Y.array[,2]==3 
	Y.array_obs[id_2] <- Z.array[Y.array[id_2,5]]
	
	#################################################################################
	# Fit Array Test data using Gaussian Predictive Process
	#
	# ===============================================================================

	t1 <- proc.time()
	res2<-try(Bayes.PP(Z=Z.array,Y=Y.array,Y.array_obs,D_1,mcx_1,mcw_1,
				   D_2,mcx_2,mcw_2,phi_ini,
				   trphi_tune,sample.parameter,kap,na=length(Se.t), 
				   Se=Se.t,Sp=Sp.t,
				   est.error=FALSE,visual=FALSE),silent=TRUE)
	t2 <- proc.time()
	########################################## save results
	if(class(res2)!='try-error'){
		Theta <- rbind(res2$beta,res2$tau_1,res2$tau_2,res2$phi_1,res2$phi_2)
		Theta.mean <- apply(Theta,1,mean)
		Theta.med <- apply(Theta,1,median)
		Theta.sd <- apply(Theta,1,sd)
		Theta.bd <- HPDinterval(as.mcmc(t(Theta)))

		Eta_1.new <- res2$X_new_1
		Eta_1 <- res2$eta_new_1
		Eta_1.mean <- apply(Eta_1,1,mean)
		Eta_1.med <- apply(Eta_1,1,median)
		Eta_1.sd <- apply(Eta_1,1,sd)
		Eta_1.bd <- HPDinterval(as.mcmc(t(Eta_1)))

		Eta_2.new <- res2$X_new_2
		Eta_2 <- res2$eta_new_2
		Eta_2.mean <- apply(Eta_2,1,mean)
		Eta_2.med <- apply(Eta_2,1,median)
		Eta_2.sd <- apply(Eta_2,1,sd)
		Eta_2.bd <- HPDinterval(as.mcmc(t(Eta_2)))
		AT_PP<-c(as.numeric(t2-t1)[3],Theta.mean,Theta.med,Theta.sd,Theta.bd,
			Eta_1.new,Eta_1.mean,Eta_1.med,Eta_1.sd,Eta_1.bd,
			Eta_2.new,Eta_2.mean,Eta_2.med,Eta_2.sd,Eta_2.bd)}else{AT_PP<-rep(NA,lrow)}

	#################################################################################
	# Fit Array Test data using Regular Gaussian Process
	#
	# ===============================================================================

	t1 <- proc.time()
	res1<-try(Bayes.GP(Z=Z.array,Y=Y.array,Y.array_obs,D_1,mcx_1,D_2,mcx_2,
				   phi_ini,trphi_tune,sample.parameter,kap,na=length(Se.t),
				   Se=Se.t,Sp=Sp.t,
				   est.error=FALSE,visual=FALSE),silent=TRUE)
	t2 <- proc.time()
	########################################## save results
	if(class(res1)!='try-error'){
		Theta <- rbind(res1$beta,res1$tau_1,res1$tau_2,res1$phi_1,res1$phi_2)
		Theta.mean <- apply(Theta,1,mean)
		Theta.med <- apply(Theta,1,median)
		Theta.sd <- apply(Theta,1,sd)
		Theta.bd <- HPDinterval(as.mcmc(t(Theta)))

		Eta_1.new <- res1$X_new_1
		Eta_1 <- res1$eta_new_1
		Eta_1.mean <- apply(Eta_1,1,mean)
		Eta_1.med <- apply(Eta_1,1,median)
		Eta_1.sd <- apply(Eta_1,1,sd)
		Eta_1.bd <- HPDinterval(as.mcmc(t(Eta_1)))

		Eta_2.new <- res1$X_new_2
		Eta_2 <- res1$eta_new_2
		Eta_2.mean <- apply(Eta_2,1,mean)
		Eta_2.med <- apply(Eta_2,1,median)
		Eta_2.sd <- apply(Eta_2,1,sd)
		Eta_2.bd <- HPDinterval(as.mcmc(t(Eta_2)))
		AT_GP<-c(as.numeric(t2-t1)[3],Theta.mean,Theta.med,Theta.sd,Theta.bd,
			Eta_1.new,Eta_1.mean,Eta_1.med,Eta_1.sd,Eta_1.bd,
			Eta_2.new,Eta_2.mean,Eta_2.med,Eta_2.sd,Eta_2.bd)}else{AT_GP<-rep(NA,lrow)}

	
	
	#write(c(mean(Y_true),nrow(Z.dorf),nrow(Z.array)),"number.txt",append=TRUE)
	
	
	###############################################################################
	# 
	# End of simulation
	# Save simulation results from one dataset
	###############################################################################
	
	RCD <- c(mean(Y_true),nrow(Z.dorf),nrow(Z.array))
	mat <- cbind(mat,c(RCD,IT_PP,IT_GP,MPT_PP,MPT_GP,DT_PP,DT_GP,AT_PP,AT_GP))
	
	# increase iteration count
	sim <- sim + 1
}

fname <- paste0(N,'_noerror_mod',mod,'.RData')

save(mat,file=fname)











