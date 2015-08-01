
library(thesis)
library(ggplot2)
library(tidyr)
library(dplyr)
library(plyr)
library(gridExtra)
library(grid)


# Functions
setwd("/Users/checorh/Documents/PHD Thesis/ProbabilisticBisection/")
source("ProbabilisticBisection.R") 

# Load noisy function
#source("/Users/checorh/Documents/PHD Thesis/ProbabilisticBisection/AmericanOptions/TimingValueFunction.R") # for american put oracle

# Plot function
a = 25
b = 40

sampling.locations <- seq(a,b,by=.5)
y <- sapply(sampling.locations,function(x){am.put.oracle(x,K=40,boundaries=bdr1,sigma=0.2,r=0.06,dt=0.04)})                          

p1 <- qplot(x=sampling.locations,y=y,geom="point",
            ylab="z(x)=g(x) + noise(x)",
            xlab="Sampling Locations (x)",
            main="Timing value function: g(x)=T(t,x)") + 
  geom_hline(yintercept=0,lty=2) 
p1
ggsave(file="plot1.png",p1)

  

set.seed(1) # Set seed

# Redefine noisy function
g <- function(x){am.put.oracle(x,K=40,boundaries=bdr1,sigma=0.2,r=0.06,dt=0.04)}

# 100 equidistant points measured 20 times
set.seed(1)
n.times=50
n.samples=1000
x.train <- seq(a,b,length=n.samples)
x=rep(x.train,n.times) # train data
g.mean <- tapply(g(x),x,mean)
      


y.prob=(sign2(g.mean) + 1)/2

# logistic regression with cubic term: success is 
# observing a negative sign
p.fit <- glm(y.prob ~ x.train + I(x.train^2),
             family = "binomial")
p.fit  <- step(p.fit,direction = "back") # model selection


# plot model
test.x <- seq(a,b,length=n.samples)
p.x <- predict(p.fit,newdata=data.frame(x=test.x),
               type = "response",
               se.fit = TRUE)

p2 <- qplot(x=test.x,y=p.x$fit,ylim=c(0,1),
            geom="line",
            ylab="s(x)=P(Y(x) = -1)",
            xlab="Sampling location (x)",
            main="Logistic regression estimation of s(x)") +
  geom_vline(xintercept = 35,lty="dotted") + 
  geom_hline(yintercept = 1/2,lty="dotted")
p2

# add empirical proportion 

#####
# Probablity of correct response
#####
p <- sapply(seq(a,b,length=n.samples),
            function(x)p.correct(p.fit,sampling.point=x))
p3 <- qplot(x=seq(a,b,length=100),
            y=p,geom="line",ylab="p(x)",
            xlab="Sampling location (x)",
            main="Probability of correct response")


pp=arrangeGrob(p2,p3,ncol=2)
ggsave(file="plot2.png",pp,width=12,height=6)


### run model

prior=data.frame(x=c(a,b),f=c(1/(b-a),1/(b-a))) # prior
#####

TotalwallClock=2000
#Rprof("fun.out")
diag.data2000 <- SchemeEvaluationData(MonteCarlo = 1,
                                  which.scheme=4,
                                  OracleCalls=c(1,5,25,50),
                                  wallClock=TotalwallClock,
                                  noisy.function=g,
                                  dd=prior)
#Rprof(NULL)
#summaryRprof("fun.out")
p <- SchemeEvaluationPlot(diag.data2000,MC=1)
ggsave(file="plot2000.png",p,width=10,height=10)


TotalwallClock=5000
#Rprof("fun.out")
diag.data5000 <- SchemeEvaluationData(MonteCarlo = 1,
                                      which.scheme=4,
                                      OracleCalls=c(1,5,25,50),
                                      wallClock=TotalwallClock,
                                      noisy.function=g,
                                      dd=prior)
#Rprof(NULL)
#summaryRprof("fun.out")
p <- SchemeEvaluationPlot(diag.data5000,MC=1)
ggsave(file="plot5000.png",p,width=10,height=10)


TotalwallClock=10000
#Rprof("fun.out")
diag.data10000 <- SchemeEvaluationData(MonteCarlo = 1,
                                      which.scheme=4,
                                      OracleCalls=c(1,5,25,50),
                                      wallClock=TotalwallClock,
                                      noisy.function=g,
                                      dd=prior)
#Rprof(NULL)
#summaryRprof("fun.out")
p <- SchemeEvaluationPlot(diag.data10000,MC=1)
ggsave(file="plot10000.png",p,width=10,height=10)


  

library(xtable)
library(taRifx)

diag.table2000 <-  diag.data2000 %>% 
  filter(wallClock==2000) %>%
  spread(key=OracleCalls,value="value") %>%
  dplyr::select(-wallClock)

diag.table5000 <-  diag.data5000 %>% 
  filter(wallClock==5000) %>%
  spread(key=OracleCalls,value="value") %>%
  dplyr::select(-wallClock)

diag.table10000 <-  diag.data10000 %>% 
  filter(wallClock==10000) %>%
  spread(key=OracleCalls,value="value") %>%
  dplyr::select(-wallClock)



save(diag.table2000,file="diagtable2000.RData")
save(diag.table5000,file="diagtable5000.RData")
save(diag.table10000,file="diagtable10000.RData")

print(latex.table.by(diag.table,2), 
      include.rownames = FALSE, 
      include.colnames = TRUE, 
      sanitize.text.function = force)

data.dens <- pbaData(which.scheme=1:6,
                       OracleCalls=c(1,5,10,20,50),
                       wallClock = 1000,
                       noisy.function=g)

p <- QuantilePlots(data.dens)
ggsave(file="plot4.png",p,width=12,height=10)

  
########
# posterior

TotalWallClock = 10000
n.boost = c(1,5,25,50)
N = TotalWallClock/n.boost
density.plots=list()

for(j in 1:length(n.boost)){
  print(N[j])
  pba <- PBA.wallClock2(dd=prior,N=N[j],scheme=4,
                        n.boost=n.boost[j],noisy.function=g)  
  
  print(pba$post[[TotalWallClock]][,1])
  
  density.plots[[j]]= ggplot(data = pba$post[[TotalWallClock]],
                             aes(x=x,y=f)) + geom_step() + 
    geom_rug(aes(x=x),sides="b",alpha = 0.3,col="blue") + 
    geom_vline(xintercept = pba$est.roots[TotalWallClock],
               col="red")
  
}


pp=arrangeGrob(density.plots[[1]],
               density.plots[[2]],
               density.plots[[3]],
               density.plots[[4]],
               ncol=4)

pp
ggsave(file="dens10000.png",pp,width=12,height=6)




PBA.wallClock5000 <- PBA.wallClock2(dd=prior,N=5000,scheme=4,n.boost=1,
                                    noisy.function=g)
PBA.wallClock2000 <- PBA.wallClock2(dd=prior,N=10000,scheme=4,n.boost=1,
                                    noisy.function=g)


qplot(x=pba$post[[TotalWallClock]][,1],
      y=pba$post[[TotalWallClock]][,2],
      geom="step",
      xlab="Sampling location (x)",
      ylab="Posterior density") + 
  geom_rug(aes(x=pba$sampl.points),sides="b",alpha = 0.05)




names(PBA.wallClock2000)


qplot(PBA.wallClock2000$post[[5000]])

qplot(PBA.wallClock2000$post[[10000]])
  
  
  
  
  
  