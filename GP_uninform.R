

##################################################################
# Divide the data into individual and GT data sets
# setwd()

data.new <- read.csv("Simulated_dataset.csv")

pool.id  <- as.numeric(data.new$Pool.ID)
data.gts <- data.new[pool.id!="NaN",]
data.ind <- data.new[pool.id=="NaN",]


##################################################################
# Creating the design matrices
X.ind<-cbind(
             as.numeric(as.character(data.ind$Age)),  
             as.numeric(data.ind$Race=="W"), 
             as.numeric(data.ind$Risk.New.Partner=="Y"),
             as.numeric(data.ind$Risk.Multiple.Partner=="Y"),
             as.numeric(data.ind$Risk.Contact=="Y"),
             as.numeric(data.ind$Symptom=="Y"),
             as.numeric(data.ind$Specimen.Type=="Swab")
)



X.pool<-cbind(
             as.numeric(as.character(data.gts$Age)),  
             as.numeric(data.gts$Race=="W"), 
             as.numeric(data.gts$Risk.New.Partner=="Y"),
             as.numeric(data.gts$Risk.Multiple.Partner=="Y"),
             as.numeric(data.gts$Risk.Contact=="Y"),
             as.numeric(data.gts$Symptom=="Y"),
             as.numeric(data.gts$Specimen.Type=="Swab")
)


###################################################################
# Building the Z and Y matrice

Z.ind.C<-matrix(-99,nrow=dim(data.ind)[1],ncol=8)
Y.ind.C<-matrix(-99,nrow=dim(data.ind)[1],ncol=4)

Z.pool.C<-matrix(-99,nrow=dim(data.gts)[1],ncol=8)
Y.pool.C<-matrix(-99,nrow=dim(data.gts)[1],ncol=4)

Z.ind.G<-matrix(-99,nrow=dim(data.ind)[1],ncol=8)
Y.ind.G<-matrix(-99,nrow=dim(data.ind)[1],ncol=4)

Z.pool.G<-matrix(-99,nrow=dim(data.gts)[1],ncol=8)
Y.pool.G<-matrix(-99,nrow=dim(data.gts)[1],ncol=4)


##################################
# For GT testing

id.ind<-1:(dim(data.gts)[1])

pool.id<-unique(data.gts$Pool.ID)
n<-length(pool.id)

track.CT<-1
track.GC<-1

#############################
#############################
### For loop starts
for(i in 1:n){

temp<-data.gts[data.gts$Pool.ID==pool.id[i],]
temp.id<-id.ind[data.gts$Pool.ID==pool.id[i]]

CT.res <- as.numeric(temp$CT.Result=="P")
GC.res <- as.numeric(temp$GC.Result=="P")

cj<-length(CT.res)

CT.retest<- sum(CT.res)>0
GC.retest<- sum(GC.res)>0

Z.pool.C[track.CT,1]<-CT.retest
Z.pool.C[track.CT,2]<-cj
Z.pool.C[track.CT,3]<-3
Z.pool.C[track.CT,4:(cj+3)]<-temp.id

Y.pool.C[temp.id,1]<-0
Y.pool.C[temp.id,2]<-1
Y.pool.C[temp.id,3]<-track.CT

Z.pool.G[track.GC,1]<-GC.retest
Z.pool.G[track.GC,2]<-cj
Z.pool.G[track.GC,3]<-3
Z.pool.G[track.GC,4:(cj+3)]<-temp.id

Y.pool.G[temp.id,1]<-0
Y.pool.G[temp.id,2]<-1
Y.pool.G[temp.id,3]<-track.GC

if(CT.retest==0){
track.CT<-track.CT+1
}

if(CT.retest>0){
tid<-(track.CT+1):(track.CT+cj)
Z.pool.C[tid,1]<-CT.res
Z.pool.C[tid,2]<-1
Z.pool.C[tid,3]<-1
Z.pool.C[tid,4]<-temp.id

Y.pool.C[temp.id,1]<-CT.res
Y.pool.C[temp.id,2]<-2
Y.pool.C[temp.id,4]<-tid
track.CT<-max(tid)+1
}

if(GC.retest==0){
track.GC<-track.GC+1
}

if(GC.retest>0){
tid<-(track.GC+1):(track.GC+cj)
Z.pool.G[tid,1]<-GC.res
Z.pool.G[tid,2]<-1
Z.pool.G[tid,3]<-1
Z.pool.G[tid,4]<-temp.id

Y.pool.G[temp.id,1]<-GC.res
Y.pool.G[temp.id,2]<-2
Y.pool.G[temp.id,4]<-tid
track.GC<-max(tid)+1
}
}

### For loop ends
#############################
#############################


Z.pool.C<-Z.pool.C[1:(track.CT-1),]
Z.pool.G<-Z.pool.G[1:(track.GC-1),]

ZnC<-dim(Z.pool.C)[1]
ZnG<-dim(Z.pool.G)[1]

YnC<-dim(Y.pool.C)[1]
YnG<-dim(Y.pool.G)[1]

#################################
# For individual level testing 

Z.ind.C[,1]<-as.numeric(data.ind$CT.Result=="P") # Test outcome for Chlamydia
Z.ind.C[,2]<-1                                   # Pool size    
Z.ind.C[,3]<- as.numeric(data.ind$Specimen.Type=="Urine")+1
Z.ind.C[,4]<-(1:(dim(data.ind)[1]))+YnC

Z.ind.G[,1]<-as.numeric(data.ind$GC.Result=="P") # Test outcome for Gonorrhea
Z.ind.G[,2]<-1                                   # Pool size    
Z.ind.G[,3]<- as.numeric(data.ind$Specimen.Type=="Urine")+1
Z.ind.G[,4]<-1:(dim(data.ind)[1])+YnG

Y.ind.C[,1]<-as.numeric(data.ind$CT.Result=="P") # Initial true status
Y.ind.C[,2]<-1
Y.ind.C[,3]<-1:(dim(data.ind)[1])+ZnC

Y.ind.G[,1]<-as.numeric(data.ind$GC.Result=="P") # Initial true status
Y.ind.G[,2]<-1
Y.ind.G[,3]<-1:(dim(data.ind)[1])+ZnG




###################################################
# Putting everything together 

X<-rbind(X.pool,X.ind)
colnames(X) <- c("Age","Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact",
"Symptom","Specimen.Type")

Z.C<-rbind(Z.pool.C,Z.ind.C)
Z.G<-rbind(Z.pool.G,Z.ind.G)

Y.C<-rbind(Y.pool.C,Y.ind.C)
Y.G<-rbind(Y.pool.G,Y.ind.G)

# dim(X)=dim(Y.C)=dim(Y.G)


##############################################

age    <- round(X[,1],1)
age.st <- (age-mean(age))/sd(age)

# colnames(X) <- c("Age","Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact",
# "Symptom","Specimen.Type")

name.in <- c("Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact","Symptom")
X_mat<-cbind(age.st,1,X[,name.in])

#############################################################################################
#############################################################################################
#############################################################################################
# Read in needed functions Note: Set R directory to the file that conatins the necessary dlls

source("BNR_RealData_GP.txt")
source("Testing functions.txt")
Rcpp::sourceCpp('SampLatent.cpp')
Rcpp::sourceCpp('ErrorUpdate.cpp')

library(mvtnorm);library(msm);library(ltsa);library(geoR);library(Matrix);library(coda)
##########################################################################################
# Specify the seed
set.seed(100)

###################
# Simulation settings
N<-nrow(Y.C)                # Number of Individuals
kap = 2.0              # smooth parameter in matern function
n.thin = 10              # thinning parameter
n.keep = 1000         # number of posterior sample
n.step = 50             # step to plot
n.burn = 2000         # burn in 
trphi_tune = 0.1      # tune parameter for decay parameter phi
L = 101                 # number of knots for Gaussian Predictive Process
x = X_mat[,1]
lbd=min(x);rbd=max(x)          # two endpoints

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

x_dist = as.matrix(dist(mcx,mcx, method = "manhattan", diag = FALSE, upper = FALSE))

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

######################################################################################
######################################################################################
# Data from Dorfman decoding with testing error rates unknown
phi_ini <- c(phi_min,phi_max)
sample.parameter <- c(n.burn,n.keep,n.thin,n.step)

######################################################################################
C.ase<-c(1,1,1)
C.bse<-c(1,1,1)
C.asp<-c(1,1,1)
C.bsp<-c(1,1,1)

set.seed(100)
t1 <- proc.time()
res1<-Bayes.GP(Z=Z.C, Y=Y.C, Y_obs=Y.C_obs,D,mcx,phi_ini,trphi_tune,sample.parameter,kap,na=3,
	             ase=C.ase, bse=C.bse, asp=C.asp, bsp=C.bsp,est.error=TRUE,visual=FALSE)
t2 <- proc.time()

save.image(file="Age_st_mean_constraint_GP_uninform.RData")


# ====================================================================================
######################################################################################

library(coda)
Sep <- rbind(res1$Se,res1$Sp)
para <- rbind(res1$beta,Sep)
round(cbind(apply(res1$beta,1,mean),HPDinterval(as.mcmc(t(res1$beta)))),3)[-1,]

round(cbind(apply(Sep,1,mean),HPDinterval(as.mcmc(t(Sep)))),3)

# 1: Swab Individual
# 2: Urine Individual
# 3: Swab Pool 


