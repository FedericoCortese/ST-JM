library(MASS)
library(cluster)
library(StatMatch)
library(caret)

simulate_observations <- function(mu=1,rho=0,phi=.8,n_states=3,P,Pcat,M,s,seed
                                  #,pNAs,typeNA=0
) {
  
  # This function simulates data from a multivariate normal distribution given the latent states sequence
  
  # Arguments:
  # mu: Mean values for data simulation (first state has mean = mu, last state has mean = -mu, and all intermediate states are equally spaced between them)
  # rho: Correlation between variables
  # n_states: Number of states
  # P: Number of features
  # Pcat: Number of categorical features
  # M: Number of spatial points
  # s: Latent states sequence
  # seed: Random seed
  # pNAs: Percentage of missing values
  # tpNA: Type of missing values (0 = random missing pattern, 1 = block (continuous) missing pattern)
  
  
  if(is.null(Pcat)){
    Pcat=floor(P/2)
  }
  
  MU=mu
  mu=seq(-mu,mu,length.out=n_states)
  
  Sigma <- matrix(rho,ncol=P,nrow=P)
  diag(Sigma)=1
  
  Sim = matrix(0, M, P * n_states)
  SimData = matrix(0, M, P)
  
  set.seed(seed)
  for(k in 1:n_states){
    u = MASS::mvrnorm(M,rep(mu[k],P),Sigma)
    Sim[, (P * k - P + 1):(k * P)] = u
  }
  
  for (i in 1:M) {
    k = s[i]
    SimData[i, ] = Sim[i, (P * k - P + 1):(P * k)]
  }
  
  if(Pcat!=0){
    SimData[,1:Pcat]=apply(SimData[,1:Pcat],2,get_cat,mc=s,mu=MU,phi=phi)
    SimData=as.data.frame(SimData)
    SimData[,1:Pcat]=SimData[,1:Pcat]%>%mutate_all(as.factor)
  }
  
  # if(pNAs>0){
  #   SimData.NA=apply(SimData,2,punct,pNAs=pNAs,type=typeNA)
  #   SimData.NA=as.data.frame(SimData.NA)
  #   SimData.NA[,1:Pcat]=SimData.NA[,1:Pcat]%>%mutate_all(as.factor)
  #   SimData.NA[,-(1:Pcat)]=SimData.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  # }
  
  #else{
  SimData.NA=SimData
  #}
  
  return(
    #list(
    SimData=SimData
    #,SimData.NA=SimData.NA
    #)
  )
  
}

generate_spatio_temporal_data <- function(M, TT, theta, beta, K = 3,
                                          mu=1,rho=0,phi=.8,
                                          P,Pcat,seed,
                                          pGap=.2,
                                          pNAs=0) {
  
  
  # Function to generate spatio-temporal data with spatial (beta) and temporal (theta) persistence
  
  # Arguments:
  # M: Number of spatial points
  # TT: Number of time points
  # theta: Spatial correlation parameter (the lower, the higher the spatial correlation)
  # beta: Temporal correlation parameter (the higher, the higher the temporal correlation)
  # K: Number of states (only K=3 available at the moment)
  # mu: Mean values for data simulation (first state has mean = mu, last state has mean = -mu, and all intermediate states are equally spaced between them)
  # rho: Correlation between variables
  # phi: Conditional probability for the categorical outcome k in state k
  # P: Number of features
  # Pcat: Number of categorical features
  # seed: Random seed
  # pGap: Percentage of time points to be removed 
  # pNAs: Percentage of missing values (only random missing pattern is available at the moment)
  
  # Value:
  # A list with the following elements:
  # S: A matrix TTxM with the simulated states
  # Y: A data frame with the complete simulated data in long format
  # Y.NA: A data frame with the simulated data with missing values in long format
  # spatial_points: A data frame with the spatial points
  # spatial_cov: The spatial covariance matrix
  # dist_matrix: The distance matrix
  
  
  
  # Increment TT by one as the first time step will be removed
  TT=TT+1
  
  
  # Generate spatial points
  spatial_points <- generate_spatial_points(M)
  
  # Create spatial covariance matrix based on distance and theta
  dist_matrix <- as.matrix(dist(spatial_points)) # Eventually substitute with Gower distance
  
  spatial_cov <- exp(-theta * dist_matrix)
  
  # Initialize data matrix
  data <- array(0, dim = c(TT, M))
  
  S=matrix(0,TT,M)
  
  # Initial time step data (from spatial process)
  data[1, ] <- mvrnorm(1, mu = rep(0, M), Sigma = spatial_cov)
  
  cluster_levels <- quantile(data[1,], probs = seq(0, 1, length.out = K + 1))
  S[1,] <- cut(data[1,], breaks = cluster_levels, labels = FALSE,
               include.lowest =T)
  
  temp=simulate_observations(mu=mu,rho=rho,phi=phi,
                             n_states=K,P=P,Pcat=Pcat,M=M,
                             s=S[1,],seed=seed+seed*1000
                             #,pNAs=pNAs,typeNA=0
  )
  
  Y=temp#$SimData
  #Y.NA=temp$SimData.NA
  
  Y=data.frame(Y)
  Y$m=1:M
  Y$t=rep(0,M)
  
  S[1,]=order_states_condMean(Y[Y$t==0,dim(Y)[2]-2],S[1,])
  
  # Y.NA=data.frame(Y.NA)
  # Y.NA$m=1:M
  # Y.NA$t=rep(0,M)
  
  # Generate data for each subsequent time step
  for (t in 2:TT) {
    eta_t <- mvrnorm(1, mu = rep(0, M), Sigma = spatial_cov)  # Spatial noise
    data[t, ] <- beta * data[t-1, ] + eta_t  # Temporal correlation
    
    #data[t, ] <- beta * data[t-1, ] + (1-beta)*eta_t  # Temporal correlation
    
    cluster_levels <- quantile(data[t,], probs = seq(0, 1, length.out = K + 1))
    S[t,] <- cut(data[t,], breaks = cluster_levels, labels = FALSE,
                 include.lowest =T)
    
    simDat=simulate_observations(mu=mu,rho=rho,phi=phi,
                                 n_states=K,P=P,Pcat=Pcat,M=M,
                                 s=S[t,],seed=seed+seed*1000+t-1
                                 #,pNAs=pNAs,typeNA=0
    )
    
    temp=data.frame(simDat#$SimData
    )
    temp$m=1:M
    temp$t=rep(t-1,M)
    Y=rbind(Y,temp)
    
    # temp=data.frame(simDat$SimData.NA)
    # temp$m=1:M
    # temp$t=rep(t-1,M)
    #Y.NA=rbind(Y.NA,temp)
    
    S[t,]=order_states_condMean(Y[Y$t==(t-1),dim(Y)[2]-2],S[t,])
    
  }
  
  data=data[-1,]
  S=S[-1,]
  Y=Y[-which(Y$t==0),]
  #Y.NA=Y.NA[-which(Y.NA$t==0),]
  
  Y <- Y %>% relocate(t,m)
  #Y.NA <- Y.NA %>% relocate(t,m)
  
  Y.NA=apply(Y[,-(1:2)],2,punct,pNAs=pNAs,type=0)
  Y.NA=as.data.frame(Y.NA)
  Y.NA[,1:Pcat]=Y.NA[,1:Pcat]%>%mutate_all(as.factor)
  Y.NA[,-(1:Pcat)]=Y.NA[,-(1:Pcat)]%>%mutate_all(as.numeric)
  Y.NA=data.frame(t=Y$t,m=Y$m,Y.NA)
  
  if(pGap>0){
    set.seed(seed)
    gaps=sort(sample(1:TT,round(TT*pGap),replace=F))
    Y=Y[-which(Y$t %in% gaps),]
    Y.NA=Y.NA[-which(Y.NA$t %in% gaps),]
  }
  
  return(list(S = S, 
              Y=Y,
              Y.NA=Y.NA,
              spatial_points = spatial_points,
              spatial_cov = spatial_cov,
              dist_matrix = dist_matrix)
  )
}

generate_spatial_points <- function(n, max_distance = 10) {
  x <- runif(n, 0, max_distance)
  y <- runif(n, 0, max_distance)
  return(data.frame(x = x, y = y))
}

order_states_condMean=function(y,s){
  
  # This function organizes states by assigning 1 to the state with the smallest conditional mean for vector y
  # and sequentially numbering each new state as 2, 3, etc., incrementing by 1 for each newly observed state.
  
  #Slong=c(t(S))
  # condMeans=sort(tapply(y,Slong,mean,na.rm=T))
  condMeans=sort(tapply(y,s,mean,na.rm=T))
  
  states_temp=match(s,names(condMeans))
  
  #states_temp=matrix(states_temp,nrow=nrow(S),byrow = T)
  
  return(states_temp)
}

STjumpDist=function(Y,n_states,
                    D,
                    jump_penalty=1e-5,
                    spatial_penalty=1e-5,
                    initial_states=NULL,
                    max_iter=10, n_init=10, tol=NULL, verbose=FALSE,timeflag=T){
  
  # This function implements the spatio-temporal jump algorithm for clustering spatio-temporal data
  # based on a distance matrix
  
  # Arguments:
  # Y is a dataframe in long format with mandatory columns t and m which are T times and M spatial indexes
  # n_states is the number of states
  # D is the distance matrix of dimension MxM
  # jump_penalty is the penalty for jumping between states
  # spatial_penalty is the penalty for jumping between spatially close points
  # initial_states is a matrix with the initial state of each point
  # max_iter is the maximum number of iterations
  # n_init is the number of initializations
  # tol is the tolerance for stopping the algorithm
  # verbose is a boolean for printing the loss at each iteration
  # timeflag is a boolean, if TRUE the times are not equally sampled
  
  # Value:
  # best_s is the best state sequence
  # Y is the imputed data
  # loss is the loss function at the optimum
  
  #
  Y=Y[order(Y$t,Y$m),]
  
  P=ncol(Y)-2
  # Time differences
  Y.orig=Y
  
  if(timeflag){
    time=sort(unique(Y.orig$t))
    dtime=diff(time)
    dtime=dtime/as.numeric(min(dtime))
    dtime=as.numeric(dtime)
    Y=subset(Y,select=-c(t,m))
    Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
  }
  
  else{
    Y=subset(Y,select=-c(t,m))
    Y=Y%>%mutate_if(is.numeric,function(x)as.numeric(scale(x)))
  }
  
  Gamma <- jump_penalty * (1 - diag(n_states))
  best_loss <- NULL
  best_s <- NULL
  
  library(dplyr)
  Y <- data.frame(t=Y.orig$t,m=Y.orig$m,Y)
  YY=subset(Y,select=-c(t,m))
  
  TT=length(unique(Y$t))
  M=length(unique(Y$m))
  
  cat.indx=which(sapply(YY, is.factor))
  cont.indx=which(sapply(YY, is.numeric))
  
  Ycont=YY[,cont.indx]
  
  Ycat=YY[,cat.indx]
  
  n_cat=length(cat.indx)
  n_cont=length(cont.indx)
  
  ###
  # Missing data imputation TBD
  # Track missings with 0 1 matrix
  Mcont=ifelse(is.na(Ycont),T,F)
  Mcat=ifelse(is.na(Ycat),T,F)
  # Initialize mu 
  mu <- colMeans(Ycont,na.rm = T)
  # Initialize modes
  mo <- apply(Ycat,2,Mode)
  # Impute missing values with mean of observed states
  for(i in 1:n_cont){
    Ycont[,i]=ifelse(Mcont[,i],mu[i],Ycont[,i])
  }
  for(i in 1:n_cat){
    x=Ycat[,i]
    Ycat[which(is.na(Ycat[,i])),i]=mo[i]
  }
  YY[,-cat.indx]=Ycont
  YY[,cat.indx]=Ycat
  
  Y[,-(1:2)]=YY
  
  ###
  
  # State initialization through kmeans++
  S=matrix(0,nrow=TT,ncol=M)
  for(m in 1:M){
    S[,m]=initialize_states(Y[which(Y$m==m),-(1:2)],n_states)
  }
  
  for (init in 1:n_init) {
    mu <- matrix(0, nrow=n_states, ncol=length(cont.indx))
    mo <- matrix(0, nrow=n_states, ncol=length(cat.indx))
    loss_old <- 1e10
    
    for (it in 1:max_iter) {
      for (i in unique(as.vector(S))) {
        mu[i,] <- colMeans(Ycont[as.vector(t(S))==i,])
        mo[i,]=apply(Ycat[as.vector(t(S))==i,],2,Mode)
      }
      
      mu=data.frame(mu)
      mo=data.frame(mo,stringsAsFactors=TRUE)
      for(i in 1:n_cat){
        mo[,i]=factor(mo[,i],levels=levels(Ycat[,i]))
      }
      mumo=data.frame(matrix(0,nrow=n_states,ncol=P))
      mumo[,cat.indx]=mo
      mumo[,cont.indx]=mu
      colnames(mumo)=colnames(YY)
      
      
      # Fit state sequence
      S_old <- S
      
      loss_v=NULL
      
      # Re-fill-in missings 
      for(i in 1:nrow(Mcont)){
        Ycont[i,]=unlist(ifelse(Mcont[i,],mu[as.vector(t(S))[i],],Ycont[i,]))
      }
      for(i in 1:nrow(Mcat)){
        Ycat[i,]=unlist(ifelse(Mcat[i,],mo[as.vector(t(S))[i],],Ycat[i,]))
      }
      
      YY[,-cat.indx]=Ycont
      YY[,cat.indx]=Ycat
      
      Y[,-(1:2)]=YY
      
      ###
      for(m in 1:M){
        loss_by_state=gower.dist(Y[which(Y$m==m),-(1:2)],mumo)
        
        for(k in 1:n_states){
          #temp <- t(t((S[,-m]==k))/D[m,-m])
          temp <- t(t((S[,-m]==k))*exp(-D[m,-m]))
          loss_by_state[,k]=loss_by_state[,k]-spatial_penalty*rowSums(temp)
        }
        
        V <- loss_by_state
        for (t in (TT-1):1) {
          #V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
          if(timeflag){
            V[t-1,] <- loss_by_state[t-1,] + apply(V[t,]/dtime[t] + Gamma, 2, min)
          }
          else{
            V[t-1,] <- loss_by_state[t-1,] + apply(V[t,] + Gamma, 2, min)
          }
        }
        
        S[1,m] <- which.min(V[1,])
        for (t in 2:TT) {
          S[t,m] <- which.min(V[t,] + Gamma[S[t-1,m],])
        }
        loss_v=c(loss_v,min(V[1,]))
      }
      if (length(unique(S)) == 1) {
        break
      }
      loss <- mean(loss_v)
      if (verbose) {
        cat(sprintf('Iteration %d: %.6e\n', it, loss))
      }
      if (!is.null(tol)) {
        epsilon <- loss_old - loss
        if (epsilon < tol) {
          break
        }
      } else if (all(S == S_old)) {
        break
      }
      loss_old <- loss
      
      ###
      
    }
    if (is.null(best_s) || (loss_old < best_loss)) {
      best_loss <- loss_old
      best_s <- S
    }
    
    for(m in 1:M){
      S[,m]=initialize_states(Y[which(Y$m==m),-(1:2)],n_states)
    }
  }
  return(list(best_s=best_s,
              Y=Y,
              K=n_states,
              lambda=jump_penalty,
              gamma=spatial_penalty,
              loss=best_loss))
  
}

get_cat=function(y,mc,mu,phi){
  # Function to simulate categorical data
  
  # Arguments:
  # y: continuous variable 
  # mc: Markov chain states
  # mu: numeric mean value
  # phi: conditional probability for the categorical outcome k in state k
  
  library(dplyr)
  
  mu=c(-mu,0,mu)
  #K=length(unique(mc))
  #mu=seq(-mu,mu,length.out=K)
  #phi1=(1-phi)/(K-1)
  phi1=(1-phi)/2
  
  TT=length(y)
  for(i in 1:TT){
    k=mc[i]
    switch(k,
           "1"={
             threshold=c(qnorm(phi1,mu[1]),qnorm(phi+phi1,mu[1]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=1
             }
             else if(y[i]<threshold[1]){
               y[i]=2
             }
             else{
               y[i]=3
             }
           },
           "2"={
             threshold=c(qnorm(phi1,mu[2]),qnorm(phi+phi1,mu[2]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=2
             }
             else if(y[i]<threshold[1]){
               y[i]=3
             }
             else{
               y[i]=1
             }
           },
           "3"={
             threshold=c(qnorm(phi1,mu[3]),qnorm(phi+phi1,mu[3]))
             if(y[i]>threshold[1]&y[i]<threshold[2]){
               y[i]=3
             }
             else if(y[i]<threshold[1]){
               y[i]=1
             }
             else{
               y[i]=2
             }
           }
    )
  }
  return(y)
  
}

punct=function(x,pNAs,typeNA){
  
  # x is a vector (column of the dataset)
  # pNAs is the percentage of missing values
  # typeNA is the type of missing values (0 for random, 1 for continuous, all other values will turn into no missing imputation)
  
  TT=length(x)
  pTT=round(TT*pNAs)
  if(typeNA==0){
    NAindx=sample(1:TT,pTT,replace = F)
    x[NAindx]=NA
  }
  else if(typeNA==1){
    NAindx=sample(1:(TT-pTT),1,replace = F)
    NAindx=seq(NAindx,NAindx+pTT)
    x[NAindx]=NA
  }
  
  return(x)
  
}

Mode <- function(x,na.rm=T) {
  if(na.rm){
    x <- x[!is.na(x)]
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
  
}

initialize_states <- function(Y, K) {
  n <- nrow(Y)
  
  ### Repeat the following few times?
  centr_indx=sample(1:n, 1)
  centroids <- Y[centr_indx, , drop = FALSE]  # Seleziona il primo centroide a caso
  
  closest_dist <- as.matrix(daisy(Y, metric = "gower"))
  closest_dist <- closest_dist[centr_indx,]
  
  for (i in 2:K) {
    prob <- closest_dist / sum(closest_dist)
    next_centr_indx <- sample(1:n, 1, prob = prob)
    next_centroid <- Y[next_centr_indx, , drop = FALSE]
    centroids <- rbind(centroids, next_centroid)
  }
  ###
  
  # init_stats=rep(0,n)
  # For cycle
  # for(i in 1:n){
  #   init_stats[i]=which.min(gower.dist(Y[i,],centroids))
  # }
  
  # Using sapply and vapply
  # init_stats2 <- sapply(1:n, function(i) which.min(gower.dist(Y[i,], centroids)))
  # init_stats3 <- vapply(1:n, function(i) which.min(gower.dist(Y[i,], centroids)), integer(1))
  
  # Faster solution 
  dist_matrix <- gower.dist(Y, centroids)
  init_stats <- apply(dist_matrix, 1, which.min)
  
  return(init_stats)
}

BAC=function(A,B,levs=3){
  
  # Compute balanced accuracy between A and B
  A=as.character(A)
  B=as.character(B)
  
  A=factor(A,levels=1:levs)
  B=factor(B,levels=1:levs)
  
  return(confusionMatrix(A,B)$overall[1])
}
