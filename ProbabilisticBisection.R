# x is a vector
# sigma (annualized volatility), r (annualized interest rate) are GBM parameters
# dt (time-step) and K (Put strike) are model parameters
bdr1 <- c(37.95088, 37.39437, 36.90874, 36.51018, 36.18202, 35.89936, 35.63883, 35.38308, 35.11630)

am.put.oracle <- function(x,K=40,boundaries=bdr1,sigma=0.2,r=0.06,dt=0.04)
{
  len <- length(x)
  M <- length(boundaries)
  
  stopPayoff <- pmax(0, K-x)
  
  curX <-  x*exp( rnorm(len)*sigma*sqrt(dt) + (r- sigma^2/2)*dt)
  
  payoff <- stopPayoff
  
  contNdx <- 1:len
  i <- 1
  
  while (i <=(M) & length(contNdx) > 0) {
    payoff[contNdx]  <- exp(-i*dt*r)*pmax(0, K- curX[contNdx])
    
    # continue if above the boundary or payoff is zero
    contNdx <- contNdx[which( curX[contNdx] > boundaries[i] | payoff[contNdx] == 0) ]
    
    curX[contNdx] <- curX[contNdx]*exp( rnorm(length(contNdx))*sigma*sqrt(dt) + (r- sigma^2/2)*dt)
    i <- i+1
  }
  payoff[contNdx]  <- exp(-i*dt*r)*pmax(0, K- curX[contNdx])
  
  return( payoff-stopPayoff)
}

EmpiricalProb <- function(noisy.function=g,n.train,n.boost){
  
  # n.train equidistant points measured n.boost times
  x <- rep(seq(0,40,length=n.train),n.boost) 
  y.prob <- (sign2(noisy.function(x)) + 1)/2
  
  # logistic regression with cubic term: success is going right
  # i.e. observing a negative sign
  p.fit <- glm(y.prob~x + I(x^2) + I(x^3) + I(x^4),family = "binomial")
  p.fit  <- step(p.fit,direction = "back") # model selection
  
  return(p.fit)
}


g <- function(x){am.put.oracle(x,K=40,boundaries=bdr1,sigma=0.2,r=0.06,dt=0.04)}



p.correct <- function(p.fit,sampling.point){
  p.x <- predict(p.fit,newdata=data.frame(x=sampling.point),type = "response")  
  p <- ifelse(p.x<1/2,1-p.x,p.x)
  
  return(p)
}

oracle <- function(noisy.function,sampling.point,increasing=TRUE){
  if(increasing==TRUE){
    y <- sign(-noisy.function(sampling.point))
  } else{
    y <- sign(noisy.function(sampling.point))
  }
  return(y)
}

PBA.am <- function(dd,N,n.boost,...){
  # Args:
  # dd: a posterior density object
  # N: Maximum number of iterations
  # x.root: true root
  # n.boost: Number of calls to the oracle at each sample
  
  post.dens <- list() # keep track of updated posterior densities
  est.roots <- vector() # keep track of estimated root (median of posterior densitty)
  sampling.point <- vector() # Sampling location obtained by minimizing cost function
  est.quant <- matrix(NA,N,6) # See how the quantiles collapse toswars the true value
  probs <- vector() # Vector of probability of moving right
  oracle.ans <- matrix(NA,N,n.boost) # Vector of probability of moving right
  sampling.quantile <- 0
  v <- list() # list of values
  
  for(i in 1:N){ # number of macro iterations
    
    if(i!=1){ # After the first iteration we use the updated pdf of X^{*}
      dd <- post.dens[[i-1]]
    } 
    
    ###
    # Sampling at argmax of mutual information
    ####
    f <- function(x){
      i.gain=thesis::In(dd,p.correct(x),x)
      return(i.gain)
    }
    
    # maximizer of mutual information
    x.sample <- optimize(f,interval=c(0,40),maximum=TRUE)$max # optimize information gain
    
    
    # Probability of correct response
    p <- p.correct(p.fit,sampling.point=x.sample)
    
    ###
    # sampling in batch
    ####
    # Generate oracle calls in batch
    y.oracle <- sign2(g(rep(x.sample,n.boost)))
    
    # Update posterior in batch
    dd <- thesis::posterior.batch(dd,x.sample,-y.oracle,p)
    
    #dd <- posterior(dd,x.sample,y.oracle,p)
    
    ###
    # Keep track of observed quantities
    ###
    
    # sampling locations
    sampling.point[i] <- x.sample
    # p correct
    probs[i] <- p
    # oracle responses
    oracle.ans[i,] <- y.oracle 
    # updated posterior
    post.dens[[i]] <- dd
    # estimated root (posterior median)
    est.roots[i] <- thesis::est.quantile(dd,.5)
    # Keep track of quantiles
    est.quant[i,] <- thesis::est.quantile(dd,c(.01,.05,.25,.75,.95,.99))
    colnames(est.quant) <- paste0("q",c(.01,.05,.25,.75,.95,.99))
    # keep track of the quantile corresponding to the sampling location
    sampling.quantile[i] <- thesis::Fn.x(dd,sampling.point[i])
    
  }
  
  return(list(post.dens=post.dens,est.roots=est.roots,sampl.points = sampling.point,
              est.quant = est.quant,p=probs,oracle.ans=oracle.ans,
              sampling.quantile=sampling.quantile))
}


SamplingScheme <- function(dd,scheme=1,...){
  
  otherParms <- list(...)
  
  if(scheme==1){ # Uniformly over [0,1]
    sampling.point <-  runif(1,min=min(dd[,1]),max=max(dd[,1])) 
  } else if(scheme==2){ # Uniformly over [0,1] in quantiles
    sampling.point <-  thesis::est.quantile(dd,runif(1,0,1)) 
  } else if(scheme==3){ # Sampling at the median of f_{n-1}
    sampling.point <- est.quantile(dd,1/2)
  } else if(scheme==4){ # Sampling at argmax of mutual information
    
    # information gain function
    f <- function(x){
      i.gain=thesis::In(dd,p.correct(otherParms$p.fit,x),x)
      return(i.gain)
    }
    print(range(dd[,1]))
    sampling.point <- optimize(f,interval=range(dd[,1]),maximum=TRUE)$max
    
  } else if(scheme==5){ 
    # Sampling outside IQR of f_{n-1} and look at the quantile with largest info gain
    
    x1 <- thesis::est.quantile(dd,0.80) # 80th quantile of f_{n}
    x2 <- thesis::est.quantile(dd,0.20) # 20th quantile of f_{n}
    
    info1 <- thesis::In(dd,p.correct(otherParms$p.fit,x1),x1)
    info2 <- thesis::In(dd,p.correct(otherParms$p.fit,x2),x2)
    
    sampling.point <- c(x1,x2)[which.max(c(info1,info2))]
  } else if(scheme==6){ # Sampling systematically both sides of the root
    
    x1 <- thesis::est.quantile(dd,0.80) # 80th quantile of f_{n}
    x2 <- thesis::est.quantile(dd,0.20) # 20th quantile of f_{n}
    
    if(otherParms$i%%2==0){
      sampling.point <- x1
    }else{
      sampling.point <- x2
    } 
  }
  return(unname(sampling.point))
}


PBA.wallClock2 <- function(dd,N,scheme=1,n.boost,noisy.function,...){
  # Args:
  # dd: a posterior density object
  # N: Maximum number of iterations
  # x.root: true root
  # n.boost: Number of calls to the oracle at each sample
  
  # wall-clock time
  wallClock = N*n.boost
  
  post.dens <- list() # keep track of updated posterior densities
  est.roots <- vector() # keep track of estimated root (median of posterior densitty)
  sampling.point <- vector() # Sampling location obtained by minimizing cost function
  est.quant <- matrix(NA,wallClock,6) # See how the quantiles collapse toswars the true value
  probs <- vector() # Vector of probability of moving right
  oracle.ans <- matrix(NA,N,n.boost) # Vector of probability of moving right
  sampling.quantile <- 0
  v <- list() # list of values
  
  macro.time <- (0:(N-1))*n.boost
  
  for(i in 1:N){ # N= total number of macro iterations
  #  print(paste("Macro iteration:",i))
    
    # Sampling location chosen according to scheme    
    sampling.point[i] <-  SamplingScheme(dd,scheme,p.fit=p.fit,i=i)
    
    # Probability of correct response
    p <- p.correct(p.fit,sampling.point[i])
    
    # call oracle as many times as n.boost for a fixed sampling location and update posterior density
    for(j in 1:n.boost){
      
      y.oracle <- oracle(noisy.function,sampling.point[i])
      oracle.ans[i,j] <- y.oracle # keep track of oracle responses
      
      # Update posterior in batch
      dd <- thesis::posterior(dd,sampling.point[i],y.oracle,p)
      
      index=macro.time[i]+j
      
      # Keep track of quantiles
      est.quant[index,] <- thesis::est.quantile(dd,c(.01,.05,.25,.75,.95,.99))
      colnames(est.quant) <- paste0("q",c(.01,.05,.25,.75,.95,.99))
      
      # estimated root (posterior median)
      est.roots[index] <- thesis::est.quantile(dd,.5)
      
      # keep track of the quantile corresponding to the sampling location
      sampling.quantile[index] <- thesis::Fn.x(dd,sampling.point[i])
      
      # keep track of posterior densities
      post.dens[[index]] <-  dd
      
    }
    
    # Keep track of probability of being correct
    probs[i] <- p
    
  }
  
  return(list(post.dens=post.dens,est.roots=est.roots,sampl.points = sampling.point,
              est.quant = est.quant,p=probs,oracle.ans=oracle.ans,
              sampling.quantile=sampling.quantile))
}


pbaData <- function(which.scheme=1:6,
         OracleCalls=c(1,5,10,20,50),
         wallClock = 100,
         noisy.function){  
  require(reshape2)
  
  # List of names of sampling schemes
  sampling.scheme = c("Uniform-x","Uniform-q","Median",
                      "Argmax InfoGain","q-max InfoGain","q-alternating")
  
  # Total number of macro iterations
  N = wallClock/OracleCalls
  
  # Generate data needed for creating the plot
  
  data.dens=NULL # To aggregate all schemes
  for(k in which.scheme){ # over sampling scheme
    #(paste("SamplingScheme",k))
    dataTemp=NULL # updated in the for loop
    for(i in 1:length(OracleCalls)){ # over number of oracle calls
     # print(paste("OracleCalls",OracleCalls[i]))
      
      pba=PBA.wallClock2(dd,N[i],scheme=k,n.boost = OracleCalls[i],noisy.function)
      
      # Run algorithm
      dataTemp=rbind(dataTemp,
                     cbind(melt(pba$est.quant[,1:3],value.name="ymin"),
                           ymax=melt(pba$est.quant[,6:4])[,3],
                           Oracle=OracleCalls[i],
                           SamplingPoint=rep(pba$sampl.points,each=OracleCalls[i]),
                           EstRoot=pba$est.roots,
                           SamplingScheme=sampling.scheme[k]))
      
    }
    data.dens=rbind(data.dens,dataTemp) # for quantiles plot
  }
  return(data.dens)
}



posterior <- function(dd,sampling.point,y.oracle,p){
  for(i in 1:length(sampling.point)){
    # CASE I: Sampling point has NOT been already included:
    x  <- sampling.point[i]
    y  <- y.oracle[i] # feb12-2015: added vectorized version of oracle answers
    
    # Get indexes for multiplication (i.e. update the posterior density): 
    if(!(x %in% dd[,1])){ 
      # bind sampling point with previous density
      dd <- rbind(dd,x) 
      # arrange domain of the posterior so it includes the new measurement
      dd <- dd[order(dd[,1]),]   
      # Index at which x is:
      (sampling.index <- which(dd[,1] == x))  
      # assign pdf of the next point mass
      # change done feb-12-2015
      (dd[sampling.index,2] <- dd[sampling.index-1,2]) 
    } 
    
    # indexes before sampling point: strictly less than sampling point
    (index.bef <- dd[,1]<x) 
    
    # indexes after the sampling point
    (index.aft <- dd[,1]>=x) 
    
    # Proportionality constant
    (const <- gammax(dd,x,p))
    
    if(y == -1 ){ # move to right
      dd[index.aft,2] <- const^{-1}*p*dd[index.aft,2] # x >= X_{n}
      dd[index.bef,2] <- const^{-1}*(1-p)*dd[index.bef,2] # x < X_{n}
      
    } else if (y == +1) { # move to left
      dd[index.aft,2] <- (1-const)^{-1}*(1-p)*dd[index.aft,2] # x >= X_{n}
      dd[index.bef,2] <- (1-const)^{-1}*p*dd[index.bef,2] # x < X_{n}
      
    } 
  }
  
  return(dd)
}


QuantileDataPlots <- function(which.scheme=1:6,
         OracleCalls=c(1,5,10,20,50),
         wallClock = 100,
         noisy.function){  
  require(reshape2)
  
  # List of names of sampling schemes
  sampling.scheme = c("Uniform-x","Uniform-q","Median",
                      "Argmax InfoGain","q-max InfoGain","q-alternating")
  
  # Total number of macro iterations
  N = wallClock/OracleCalls
  
  # Generate data needed for creating the plot
  
  data.dens=NULL # To aggregate all schemes
  for(k in which.scheme){ # over sampling scheme
    #print(paste("SamplingScheme",k))
    dataTemp=NULL # updated in the for loop
    for(i in 1:length(OracleCalls)){ # over number of oracle calls
      #print(paste("OracleCalls",OracleCalls[i]))
      
      pba=PBA.wallClock2(dd,N[i],scheme=k,n.boost = OracleCalls[i],noisy.function)
      
      # Run algorithm
      dataTemp=rbind(dataTemp,
                     cbind(melt(pba$est.quant[,1:3],value.name="ymin"),
                           ymax=melt(pba$est.quant[,6:4])[,3],
                           Oracle=OracleCalls[i],
                           SamplingPoint=rep(pba$sampl.points,each=OracleCalls[i]),
                           EstRoot=pba$est.roots,
                           SamplingScheme=sampling.scheme[k]))
      
    }
    data.dens=rbind(data.dens,dataTemp) # for quantiles plot
  }
  return(data.dens)
}





QuantilePlots <- function(data.dens){
  
  require(ggplot2)
  # plot realization
  p=ggplot(data=data.dens)
  # plot density shapes
  p= p + geom_ribbon(aes(x=Var1,group=Var2,ymin=ymin, ymax=ymax),alpha=.1)
  
  
  # plot estimated root
  p = p + geom_line(aes(x=Var1,y=EstRoot,colour="2")) + facet_grid(Oracle~SamplingScheme)
  
  # plot sampling locations
  p = p + geom_point(aes(x=Var1,y=SamplingPoint,colour="1"),
                     size=.5)
  
  # change layout for better visualization
  p = p + theme_bw()
  
  p=p + theme(legend.position = "bottom") +
    scale_colour_discrete(name="", 
                          labels = c("Sampling location","Estimated root"), 
                          guide = guide_legend(override.aes = list(alpha = 1))) 
  
  p=p+labs(title = "Sampling schemes posterior density comparison",
           x="Wall-clock time",
           y="Sampling location")
  
  return(p)
}


SchemeEvaluationData <- function(MonteCarlo = 2,
         which.scheme=1:6,
         OracleCalls=c(1,5,10,20,50),
         wallClock = 100,
         noisy.function,
         dd,...){
  
  require(plyr)
  require(dplyr)
  require(reshape2)
  # Determine which sampling scheme to plot
  sampling.scheme = c("Uniform-x","Uniform-q",
                      "Median",
                      "Argmax MutualInfo",
                      "q-argmax MutualInfo",
                      "q-argmax alternating")
  
  
  if(sum(wallClock%%OracleCalls)!=0){
    stop("Number of oracle calls should be a  multiple of total number of wall-clock iterations")
  } else{
    # Total number of macro-iterations
    Total.MacroIter=wallClock/OracleCalls
  }
  
  # matrices for summary statistics
  entropy.model <- matrix(NA,wallClock,length(OracleCalls))
  iqr.fn <- matrix(NA,wallClock,length(OracleCalls))
  p.eff<- matrix(NA,wallClock,length(OracleCalls))
  estRoot <- matrix(NA,wallClock,length(OracleCalls))
  
  # identifiers for rows and columns of performance matrix
  l <- paste0("Oracle.",OracleCalls)
  colnames(entropy.model) <- colnames(estRoot) <- l
  colnames(iqr.fn) <- colnames(p.eff) <- l
  #### 
  
  ####
  # list for keeping track about diagnostics in each method
  ####
  
  model.diag <- list()
  
  # iterate over sampling scheme
  for(scheme in which.scheme){
    
    # Reinitialize lists at each sampling scheme
    diagnostics=list()
    diagnostics$estRoot.overall=list()
    diagnostics$entropy.overall=list()
    diagnostics$iqr.overall=list()
    diagnostics$p.effective=list()
    
    #print(paste0("Sampling Scheme:",scheme))
    
    for(k in 1:MonteCarlo){ # number of replicates per sampling scheme
     # print(paste0("MonteCarlo iteration",k))
      
      # number of calls to the oracle per sampling location
      for(j in 1:length(OracleCalls)){ 
        
        # Number of oracle calls per sampling location
        n.boost <- OracleCalls[j]
        
        pba <- PBA.wallClock2(dd,N=Total.MacroIter[j],
                              scheme,n.boost,noisy.function)
        
        # Compute summary statistics
        idx <- seq(1,wallClock,by=n.boost)
        # Empirical entropy
        entropy.model[,j] <- rep(sapply(pba$post[idx],thesis::emp.entropy),
                                 each=n.boost) 
        #  IQR  
        if(length(idx)==1){
          iqr.fn[,j] <- rep(pba$est.quant[idx,4] - pba$est.quant[idx,3],each=n.boost) 
        } else{
          iqr.fn[,j] <- rep(apply(pba$est.quant[idx,],1,function(x)x[4]-x[3]),each=n.boost) 
        }
        
        # Estimated root
        estRoot[,j]  <- rep(pba$est.roots[idx],each=n.boost) 
        
        # effective p
        p.eff[,j] <- rep(pba$p,each=n.boost) 
        
        
      } # finish j loop 
      # Keep track of resulting values in a list
      diagnostics$entropy.overall[[paste0("MC.",k)]] <- entropy.model
      diagnostics$iqr.overall[[paste0("MC.",k)]] <- iqr.fn
      diagnostics$estRoot.overall[[paste0("MC.",k)]] <- estRoot
      diagnostics$p.effective[[paste0("MC.",k)]] <- p.eff
      
    } # end for loop in k
    
    model.diag[[sampling.scheme[scheme]]]  <- diagnostics
  } # end loop for sampling scheme
  
  
  # Prepare list for a data set
  #add index of wall-clock iteration
  ul0 = lapply(model.diag,function(x)lapply(x,function(x){lapply(x,function(x)cbind(wallClock=1:wallClock,x))}))
  
  # put MC simulations together per metric and per sampling scheme
  ul1 = lapply(ul0,function(x)lapply(x,function(x){ldply(x,.id="num.iter")}))
  
  # put metric together per sampling scheme
  ul2 = lapply(ul1,function(x)ldply(x,.id="PerformanceMetric"))
  
  ul3 = lapply(ul2,
               function(x){
                 x %>% 
                   dplyr::group_by(PerformanceMetric,wallClock) %>% 
                   dplyr::summarise_each(funs(mean)) %>%
                   dplyr::select(-num.iter)
               })
  
  ul4=ldply(ul3,.id="SamplingScheme")
  
  ul5=melt(ul4,id.vars=names(ul4)[1:3],variable.name="OracleCalls")
  
  return(ul5)
}


SchemeEvaluationPlot <- function(diag.data, ...){
  ellipsis <- list(...)
  require(ggplot2)
  
  # change labels of facets
  levels(diag.data$PerformanceMetric) <- c("Est. root","Entropy", "IQR","P(Correct)")
  
  p <- ggplot(data = diag.data, aes(x=as.integer(wallClock),y=value)) +
    geom_line(aes(colour=OracleCalls)) +
    facet_grid(PerformanceMetric~SamplingScheme,scales = "free_y") 
  
  p <- p + theme(legend.position = "bottom", # position of legend
                 plot.title = element_text(size = 16), # size of overall title
                 axis.title=element_text(size=14), # format of x,y lables
                 strip.text.x = element_text(size = 10),
                 strip.text.y = element_text(size = 12)) 
  
  p <- p+labs(title = paste0("Sampling scheme performance comparison. MC=",ellipsis$MC),
              x="Wall-clock time (T)",
              y="Performance Metric",
              colour="Batch-Size (K)") 
  
  temp <- gregexpr("[0-9]+", levels(diag.data$OracleCalls))
  SizeBatch <- unlist(regmatches(levels(diag.data$OracleCalls), temp))
  
  p <- p + scale_colour_discrete(labels=SizeBatch) 
  
  return(p)
  
}


sign2 <- function(x){
  ifelse(x>=0,1,-1)
}














