

#############################################################################################
#############################################################################################
#############################################################################################
# Read in needed functions Note: Set R directory to the file that conatins the necessary dlls


source("BNR_RealData_PP.txt")
source("Testing functions.txt")
Rcpp::sourceCpp('SampLatent.cpp')
Rcpp::sourceCpp('ErrorUpdate.cpp')


library(mvtnorm);library(msm);library(ltsa);library(geoR);library(Matrix);library(coda)

##################################################################
# Divide the data into individual and GT data sets

data.new <- read.csv("Simulated_dataset.csv")

dat <- format_realdata(data.new)

X   = dat$X
Y.C = dat$Y
Z.C = dat$Z


# dim(X)=dim(Y.C)=dim(Y.G)
# Z = A matrix of testing responses whose jth row is of the form 
#    (col1=Z_j, col2=cj, col3=assay used, col4:col(4+cj-1)=indices 
#    of the individuals assigned to the jth pool) 
# X = Covariate matrix in which the ith row contains the covariate 
#    information pertaining to the ith individual
# Y = matrix whose ith row is of the form (col1=0, col2=number 
#    of pools ith individual involved in, col3:col(3+col2-1)=pools 
#    ith individual was assigned to)

##############################################

age    <- round(X[,1],1)
age.st <- (age-mean(age))/sd(age)

# colnames(X) <- c("Age","Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact",
# "Symptom","Specimen.Type")

name.in <- c("Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact","Symptom")
X_mat<-cbind(age.st,1,X[,name.in])


##########################################################################################
# Specify the seed
set.seed(100)


###################
# Simulation settings
N<-nrow(Y.C)                # Number of Individuals
kap = 2.0                   # smooth parameter in matern function
n.thin = 10                 # thinning parameter
n.keep = 1000               # number of posterior sample
n.step = 50                 # step to plot
n.burn = 2000               # burn in 
trphi_tune = 0.1            # tune parameter for decay parameter phi
L = 100                     # number of knots for Gaussian Predictive Process
x = X_mat[,1]
lbd=min(x);rbd=max(x)       # two endpoints

##########################################################################################

id_eta = rep(NA,N)
#######################################


V = X_mat[,-1];VtV = crossprod(V,V)
####################################################
# Prepare data

##################### to get the M_mat matrix
##################### to get the M_mat matrix
x <- age
dat = data.frame(table(x))
mcx = as.numeric(as.character(factor(dat[,1], levels = unique(dat[,1]))))
d = as.numeric(dat[,2])
K = length(d)
D = Diagonal(x=d)
for(i in 1:K){
	id_eta[which(x==mcx[i])] = i
}
M_mat = Matrix(0,N,K)
M_mat[cbind(1:N,id_eta)] = 1

############################################## model 1
x = X_mat[,1]

lbd=min(x);rbd=max(x)          # two endpoints
dat = data.frame(table(x))
mcx = as.numeric(as.character(factor(dat[,1], levels = unique(dat[,1]))))

#############################################################
# Define the knots, calculate the three distance matrix
#lbd <- as.numeric(sign(min(x))*(abs(trunc(min(x)*10))+1)/10)
#rbd <- as.numeric(sign(max(x))*(abs(trunc(max(x)*10))+1)/10)  # endpoints for GPP_1
mcw = seq(lbd,rbd,length.out=100)   # Knots for Predictive Process

x_dist = as.matrix(dist(mcx,mcx, method = "manhattan", diag = FALSE, upper = FALSE))
knots_dist = as.matrix(dist(mcw,mcw, method = "manhattan", diag = FALSE, upper = FALSE))

cross_dist = matrix(NA,K,L)
for(j in 1:L){
	cross_dist[,j] <- abs(mcx-mcw[j])
}
id0_mcw <- which.min(abs(mcw))
#############################################################
# determine the lower and upper bound of decay parameter phi
x_length = rbd-lbd
phi_range = seq(0.01,100,by=0.01)
lw = rep(NA,length(phi_range))
up = lw
for(i in 1:length(phi_range)){
	lw[i] = matern(x_length*0.07,phi=phi_range[i],kappa = kap)
	up[i] = matern(x_length*0.7,phi=phi_range[i],kappa = kap)
}
phi_min = phi_range[which.min(abs(lw-0.05))]
phi_max = phi_range[which.min(abs(up-0.05))]
phi_0 = (phi_min+phi_max)/2


################################################################################## Dorfman Test
################################################################################## Dorfman Test
################################################################################## Dorfman Test
# Formatting Dorfman testing data
Y.C_obs <- rep(NA,N)
Y.C_obs[Y.C[,2]==1] <- 0
id_2 <- Y.C[,2]==2 
Y.C_obs[id_2] <- Z.C[Y.C[id_2,4]]

sum(Y.C_obs)

######################################################################################
######################################################################################
# Data from Dorfman decoding with testing error rates unknown

C.ase<-c(1,1)
C.bse<-c(1,1)
C.asp<-c(1,1)
C.bsp<-c(1,1)



t3 <- proc.time()
Res_PP_C_uninfo<-Bayes.NPReg.GT.PP(Z=Z.C, Y=Y.C, Y_obs=Y.C_obs,D,mcx,mcw,n.burn,n.keep,n.thin,n.step,trphi_tune,phi_min,phi_max,phi_0=phi_0,
																	 kap, na=2, ase=C.ase, bse=C.bse,asp=C.asp, bsp=C.bsp, est.error=TRUE)
t4 <- proc.time()

t4-t3

save.image(file="Age_st_mean_constraint_GPP_uninform.RData")

# ====================================================================================
######################################################################################

library(coda)
Sep <- rbind(Res_PP_C_uninfo$Se,Res_PP_C_uninfo$Sp)
para <- rbind(Res_PP_C_uninfo$beta,Sep)
round(cbind(apply(Res_PP_C_uninfo$beta,1,mean),HPDinterval(as.mcmc(t(Res_PP_C_uninfo$beta)))),3)[-1,]

round(cbind(apply(Sep,1,mean),HPDinterval(as.mcmc(t(Sep)))),3)


tr_beta <- c(-1.769,-0.175,0.149,0.179,0.769,0.155)
Se.C <- c(0.942,0.947) # Female Swab, Female Urine
Sp.C <- c(0.976,0.989) # Female Swab, Female Urine

# 1: Swab Pool/Individual
# 2: Urine Individual
