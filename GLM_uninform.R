
# This program is designed to analyze the real group testing data, obtained by Dorfman Testing algorithm


setwd("C:/Users/liuyanyxy/Desktop/crap/GAM/Simulate")
##################################################################
# Divide the data into individual and GT data sets

data.new <- read.csv("Simulated_dataset.csv")

pool.id  <- as.numeric(data.new$Pool.ID)
data.gts <- data.new[pool.id!="NaN",]
data.ind <- data.new[pool.id=="NaN",]

dim(data.gts)
dim(data.ind)

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

dim(X.ind)
X.ind[1:10,]
###################################################################
# Building the Z and Y matrice

Z.ind.C<-matrix(-99,nrow=dim(data.ind)[1],ncol=8)
Y.ind.C<-matrix(-99,nrow=dim(data.ind)[1],ncol=4)

Z.pool.C<-matrix(-99,nrow=dim(data.gts)[1],ncol=8)
Y.pool.C<-matrix(-99,nrow=dim(data.gts)[1],ncol=4)


##################################
# For GT testing

id.ind<-1:(dim(data.gts)[1])

pool.id<-unique(data.gts$Pool.ID)
n<-length(pool.id)

track.CT<-1


#############################
#############################
### For loop starts

for(i in 1:n){

temp<-data.gts[data.gts$Pool.ID==pool.id[i],]
temp.id<-id.ind[data.gts$Pool.ID==pool.id[i]]

CT.res <- as.numeric(temp$CT.Result=="P")

cj<-length(CT.res)

CT.retest<- sum(CT.res)>0

Z.pool.C[track.CT,1]<-CT.retest
Z.pool.C[track.CT,2]<-cj
Z.pool.C[track.CT,3]<-1 # Swab Assay
Z.pool.C[track.CT,4:(cj+3)]<-temp.id

Y.pool.C[temp.id,1]<-0
Y.pool.C[temp.id,2]<-1
Y.pool.C[temp.id,3]<-track.CT


if(CT.retest==0){
track.CT<-track.CT+1
}

if(CT.retest>0){
tid<-(track.CT+1):(track.CT+cj)
Z.pool.C[tid,1]<-CT.res
Z.pool.C[tid,2]<-1
Z.pool.C[tid,3]<-1 # Swab Assay
Z.pool.C[tid,4]<-temp.id

Y.pool.C[temp.id,1]<-CT.res
Y.pool.C[temp.id,2]<-2
Y.pool.C[temp.id,4]<-tid
track.CT<-max(tid)+1
}


}

### For loop ends
#############################
#############################


Z.pool.C<-Z.pool.C[1:(track.CT-1),]

ZnC<-dim(Z.pool.C)[1]
YnC<-dim(Y.pool.C)[1]


#################################
# For individual level testing 

Z.ind.C[,1]<-as.numeric(data.ind$CT.Result=="P") # Test outcome for Chlamydia
Z.ind.C[,2]<-1                                   # Pool size    
Z.ind.C[,3]<- as.numeric(data.ind$Specimen.Type=="Urine")+1 #Urine Assay
Z.ind.C[,4]<-(1:(dim(data.ind)[1]))+YnC


Y.ind.C[,1]<-as.numeric(data.ind$CT.Result=="P") # Initial true status
Y.ind.C[,2]<-1
Y.ind.C[,3]<-1:(dim(data.ind)[1])+ZnC


###################################################
# Putting everything together 

X<-rbind(X.pool,X.ind)
colnames(X) <- c("Age","Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact",
"Symptom","Specimen.Type")

Z.C<-rbind(Z.pool.C,Z.ind.C)
Y.C<-rbind(Y.pool.C,Y.ind.C)

sum(Z.C[,1])

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

age<-X[,1]
age.st<- (age-mean(age))/sd(age)

# colnames(X) <- c("Age","Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact",
# "Symptom","Specimen.Type")

name.in <- c("Race_W","Risk.New.Partner","Risk.Multiple.Partner","Risk.Contact","Symptom")

X_mat<-cbind(1,age.st,X[,name.in])
#############################################################################################
#############################################################################################
#############################################################################################


#source("Testing functions.txt")
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
V <- X_mat; VtV = crossprod(V,V)
est.error <- TRUE
na <- 2               # number of assay accuracy sets

C.ase<-c(1,1)
C.bse<-c(1,1)
C.asp<-c(1,1)
C.bsp<-c(1,1)
ase=C.ase; bse=C.bse
asp=C.asp; bsp=C.bsp
##########################################################################################

################################################################################## Dorfman Test
################################################################################## Dorfman Test
################################################################################## Dorfman Test
# Simulates Dorfman testing
Y.C_obs <- rep(NA,N)
Y.C_obs[Y.C[,2]==1] <- 0
id_2 <- Y.C[,2]==2 
Y.C_obs[id_2] <- Z.C[Y.C[id_2,4]]

Z<-Z.C; Y<-Y.C; Y_obs<-Y.C_obs
######################################################################################
######################################################################################
# Data from Dorfman decoding with testing error rates known


Y_save <- matrix(NA,N,n.keep)
J = nrow(Z)
N = nrow(Y)
count_inverse = 0;count_inverse_new = 0
if(est.error==TRUE){
	Se.mat<-matrix(-99,nrow=na,ncol=n.keep)
	Sp.mat<-matrix(-99,nrow=na,ncol=n.keep)
}

############################### speficy the initial values of unknown parameters and random terms
a=2;b=1;tau_iter = rgamma(1,a,b)
U_iter = rep(0,N)
beta_iter =rep(0.0,ncol(V));beta_iter[1]<-qnorm(mean(Y_obs));Vb_iter = V%*%beta_iter
Sig_inv_beta=diag(rep(1/1000,length(beta_iter)))
logit <- function(x){log(x/(1-x))}
############################### create the matrix to store posterior samples
beta = matrix(NA,length(beta_iter),n.keep)
############################### initialize Y_obs
Y[,1] = Y_obs
y = Y[,1]
id0 = (y==0)
id1 = (y==1)
##########################
#GIBBS SAMPLER STARTS HERE
iter = 1
ikeep = 1
GI = n.burn + n.thin*n.keep
accept.phi = 0

while(iter<=GI){
	#################################################################################
	#################################################################################
	#################################################################################
	
	########################################################################### update Se and Sp
	if(est.error==TRUE){
		###################
		#Sampling Se and Sp
		PS<-matrix(0,nrow=na,ncol=4)
		
		res<-errorupdate(N,J,Y,Z,PS,na)
		Se.up<-matrix(res[,1:2],na,2)
		Sp.up<-matrix(res[,3:4],na,2)
		
		Se<-rbeta(na,ase+Se.up[,1],bse+Se.up[,2])
		Sp<-rbeta(na,asp+Sp.up[,1],bsp+Sp.up[,2])
	}
	
	########################################################################### update U
	U_iter[id0] = rtnorm(sum(id0), mean=(Vb_iter)[id0], sd=1, lower=-Inf, upper=0)
	U_iter[id1] = rtnorm(sum(id1), mean=(Vb_iter)[id1], sd=1, lower=0, upper=Inf)
	
	########################################################################### update beta
	beta_sig = solve(VtV+Sig_inv_beta);beta_mu=beta_sig%*%crossprod(V,U_iter)
	beta_iter = as.numeric(rmvnorm(1,beta_mu,beta_sig))
	Vb_iter = V%*%beta_iter
	
	########################################################################### update the latent variable Y
	# Sample true status Y
	p = pnorm(as.numeric(Vb_iter))
	uf<-runif(N)
	Y[,1] = SampLatent(N,p,Y=Y,Z=Z,U=uf,se=Se, sp=Sp,na=1)
	y = Y[,1]
	id0 = (y==0)
	id1 = (y==1)

	if(iter>n.burn&iter%%n.thin==0){
		

		beta[,ikeep] <- beta_iter
		Se.mat[,ikeep]<-Se
		Sp.mat[,ikeep]<-Sp
		Y_save[,ikeep] <- y
		if(ikeep%%(n.keep/10)==0){
				print(paste0("Collected ",(ikeep%/%(n.keep/10))*10,"% of required sample"))
		}
		ikeep <- ikeep + 1
		
	}
	#print(iter)
	iter = iter + 1
	
} # End Gibbs sampling 


save.image(file="Chris_method_uninform.RData")

# ====================================================================================
######################################################################################

par(mfrow=c(3,2),mar=c(2,2,2,1))
for(i in 1:3){
plot(Se.mat[i,],type="l")
plot(Sp.mat[i,],type="l")
}

par(mfrow=c(2,4),mar=c(2,2,2,1))
for(i in 1:nrow(beta)){
plot(beta[i,],type="l")
}

library(coda)
Sep <- rbind(Se.mat,Sp.mat)
para <- rbind(beta,Sep)
round(cbind(apply(beta,1,mean),HPDinterval(as.mcmc(t(beta)))),3)[-1,]

round(cbind(apply(Sep,1,mean),HPDinterval(as.mcmc(t(Sep)))),3)

# 1: Swab Individual
# 2: Urine Individual
# 3: Swab Pool 

#################################################################################
# inference of age on original scale
round(cbind(apply(para,1,mean),HPDinterval(as.mcmc(t(para))))[2,]/sd(age),3)


xp <- sort(age.st)*sd(age)+mean(age)
yp <- sort(age.st)*mean(beta[2,])

xp <- seq(min(age),max(age),0.01)
y0 <- xp*mean(beta[2,])/sd(age)
yp <- y0-mean(y0)
plot(xp,yp,type='l',lwd=2)
lines(xp,y0)



write.csv(cbind(xp,yp),'Age_parametric_plot.csv',row.names=F)



