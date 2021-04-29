library(dplyr)
library(tidyverse)
library(ggplot2)

#####################################################
#####################################################
# !!!        n 1:100 took 5 hours to run        !!! #
#####################################################
#####################################################

# Set seed for reproducibility
set.seed(123)

#
## START: Define custom functions for each of the 6 required CI methods
#

# Wald interval
waldCI <- function(y, n, alpha = 0.05){
  data.frame(lower=y/n - qnorm(1-alpha/2)*sqrt(((y/n)*(1-y/n))/n),
             upper=y/n + qnorm(1-alpha/2)*sqrt(((y/n)*(1-y/n))/n)
  )
}

# Adjusted Wald interval
adjWaldCI <- function(y, n, alpha = 0.05){
  pAdj = (y+2)/(n+4)
  data.frame(lower=pAdj - qnorm(1-alpha/2)*sqrt(((pAdj)*(1-pAdj))/(n+4)),
             upper=pAdj + qnorm(1-alpha/2)*sqrt(((pAdj)*(1-pAdj))/(n+4))
  )
}

# Clopper-Pearson (exact) interval
exactCI <- function(y, n, alpha = 0.05){
  # If y=0 or y=n set CI to (0,0) or (1,1), else calculate interval from F dist.
  data.frame(lower=ifelse(y==0,0,
                          ifelse(y==n,1,
                                 1/(1+(n-y+1)/(y*qf(alpha/2,2*y,2*(n-y+1)))))),
             upper=ifelse(y==0,0,
                          ifelse(y==n,1,
                                 1/(1+(n-y)/((y+1)*qf(1-alpha/2,2*(y+1),2*(n-y))))))
  )
}

# Wilson score interval
scoreCI <- function(y, n, alpha = 0.05){
  data.frame(lower=(((y/n)+((qnorm(alpha/2)^2)/(2*n)))+(qnorm(alpha/2)*sqrt((((y/n)*(1-y/n))/n)+(qnorm(alpha/2)^2)/(4*(n^2)))))/(1+((qnorm(alpha/2)^2)/n)),
             upper=(((y/n)+((qnorm(alpha/2)^2)/(2*n)))-(qnorm(alpha/2)*sqrt((((y/n)*(1-y/n))/n)+(qnorm(alpha/2)^2)/(4*(n^2)))))/(1+((qnorm(alpha/2)^2)/n))
  )
}

# Raw percentile interval using parametric bootstrap
bootstrapRawCI <- function(y, N, n, b, alpha = 0.05){
  dataBoot <- data.frame()
  for (i in 1:length(y)){
    y_i=y[i] # Define the y being evaluated from the y table of N resamples
    pHat <- y_i/n
    yBoot <- rbinom(n=b,size=n,prob=pHat) # Generate b bootstrap samples
    pBoot <- yBoot/n # Calculate bootstrap sample proportions
    dataBoot <- bind_rows(dataBoot,
                          data.frame(# If y=0 or y=n set CI to (0,0) or (1,1), else find quantiles
                                     lower=ifelse(y_i==0,0,
                                                  ifelse(y_i==n,1,
                                                         quantile(pBoot,alpha/2))),
                                     upper=ifelse(y_i==0,0,
                                                  ifelse(y_i==n,1,
                                                         quantile(pBoot,1-alpha/2)))
                          )
    )
  } # END Loop
  return(dataBoot)
}

# Bootstrap t interval using a parametric bootstrap
bootstrapTCI <- function(y, N, n, b, alpha = 0.05){
  dataBoot <- data.frame()
  for (i in 1:length(y)){
    y_i=y[i] # Define the y being evaluated from the y table of N resamples
    pHat <- y_i/n
    yBoot <- rbinom(n=b,size=n,prob=pHat) # Generate b bootstrap samples
    pBoot <- yBoot/n # Calculate bootstrap sample proportions
    throwOut <- ifelse(pBoot==0 | pBoot==1,FALSE,TRUE) # TRUE=keep, FALSE=discard
    pBoot <- pBoot[throwOut] # Discard any bootstrap sample proportions of 0 or 1
    tBoot <- (pBoot-pHat)/sqrt((pBoot*(1-pBoot))/n) # Calculate t-values
    dataBoot <- bind_rows(dataBoot,
                          data.frame(# If y=0 or y=n set CI to (0,0) or (1,1), else find quantiles
                                     lower=ifelse(y_i==0,0,
                                                  ifelse(y_i==n,1,
                                                         pHat-quantile(tBoot,1-alpha/2)*sqrt((pHat*(1-pHat))/n))),
                                     upper=ifelse(y_i==0,0,
                                                  ifelse(y_i==n,1,
                                                         pHat-quantile(tBoot,alpha/2)*sqrt((pHat*(1-pHat))/n)))
                          )
    )
  } # END Loop
  return(dataBoot)
}

#
## END: Define custom functions
#

# Create a master dataSet to hold results for all combinations of n and p
dataSet <- data.frame()

# Define N, number of samples for rbinom() in the loops for data simulation
N <- 1500

# Define B, number of bootstrap resamples
B <-200

loopBegin <- Sys.time()

# Loops used for iterating through all combinations of n and p, with steps:
# 1. Data simulation from the binomial distribution 
# 2. Apply the 6 custom functions to the simulated data
# 3. Build the master dataSet table used for plotting
for (n in 1:100) {
  for (p in seq(from=.01,to=.99,by=.01)){
    print(paste0("Calculating CIs for n=",n," and p=",p,
                 " @ ",Sys.time()," --- ",round(nrow(dataSet)/9900*100,1),
                 "% complete.")) # Time and progress marker to console
    # Simulation with Binomial Distribution
    y <- rbinom(n=N,size=n,prob=p)
    pHat <- y/n
    # Calculation - Wald
    wald <- as.data.frame(waldCI(y,n)) # Call custom function
    wald <- wald %>% mutate(contained=ifelse(lower<=p & upper>=p,1,0))
    waldLength <- mean(wald$upper-wald$lower)
    waldContain <- sum(wald$contained)/N
    # Calculation - Adjusted Wald
    adjWald <- as.data.frame(adjWaldCI(y,n)) # Call custom function
    adjWald <- adjWald %>% mutate(contained=ifelse(lower<=p & upper>=p,1,0))
    adjWaldLength <- mean(adjWald$upper-adjWald$lower)
    adjWaldContain <- sum(adjWald$contained)/N
    # Calculation - Clopper-Pearson (exact)
    exact <- as.data.frame(exactCI(y,n)) # Call custom function
    exact <- exact %>% mutate(contained=ifelse(lower<=p & upper>=p,1,0))
    exactLength <- mean(exact$upper-exact$lower)
    exactContain <- sum(exact$contained)/N
    # Calculation - Wilson score interval
    score <- as.data.frame(scoreCI(y,n)) # Call custom function
    score <- score %>% mutate(contained=ifelse(lower<=p & upper>=p,1,0))
    scoreLength <- mean(score$upper-score$lower)
    scoreContain <- sum(score$contained)/N
    # Calculation - Raw Percentile Interval using parametric boostrap
    raw <- bootstrapRawCI(y,N,n,b=B) # Call custom function
    raw <- raw %>% mutate(contained=ifelse(lower<=p & upper>=p,1,0))
    rawLength <- mean(raw$upper-raw$lower)
    rawContain <- sum(raw$contained)/N
    # Calculation - Bootstrap T Interval using parametric boostrap
    t <- bootstrapTCI(y,N,n,b=B) # Call custom function
    t <- t %>% mutate(contained=ifelse(lower<=p & upper>=p,1,0))
    tLength <- mean(t$upper-t$lower)
    tContain <- sum(t$contained)/N
    # Build Master Data Set
    dataSet <- bind_rows(dataSet,
                         data.frame(n=n,
                                    p=p,
                                    WaldLength=waldLength,
                                    WaldContain=waldContain,
                                    AdjWaldLength=adjWaldLength,
                                    AdjWaldContain=adjWaldContain,
                                    ExactLength=exactLength,
                                    ExactContain=exactContain,
                                    ScoreLength=scoreLength,
                                    ScoreContain=scoreContain,
                                    RawLength=rawLength,
                                    RawContain=rawContain,
                                    BoottLength=tLength,
                                    BoottContain=tContain
                         )
    ) # END dataSet build
  } # END p loop
} # END n loop

loopEnd <- Sys.time()

#
## PLOTTING
#

# Split Master Data Set into two sets, one for Length and one for Contain
dataSetLength <- select(dataSet,c(n,p,ends_with("Length")))
dataSetContain <- select(dataSet,c(n,p,ends_with("Contain")))

# Gather the two data sets (convert from wide to narrow) for easy plotting
dataPlotLength <- dataSetLength %>%
  gather(variable,value,3:length(dataSetLength)) %>%
  mutate(Method=substring(variable,1,nchar(variable)-6))
dataPlotContain <- dataSetContain %>%
  gather(variable,value,3:length(dataSetContain)) %>%
  mutate(Method=substring(variable,1,nchar(variable)-7))

# Create plots for Length
lineLength <- ggplot(data=dataPlotLength,aes(x=p,y=value,color=Method))
lineLength +
  geom_line(size=1) +
  facet_grid(cols = vars(n),labeller = label_both) +
  labs(x="p",y="Average Interval Length") +
  theme_bw()

# Create plots for Contain
# Use data=subset(dataPlotContain,Method=="") to plot specific methods
# Use data=subset(dataPlotContain,n==) to plot specific n
lineContain <- ggplot(data=subset(dataPlotContain,n==100),aes(x=p,y=value,color=Method))
lineContain +
  geom_line(size=1) +
  ylim(0,1) +
  geom_hline(yintercept = .95, color="black") +
  facet_wrap(vars(Method),labeller = label_both) +
  #facet_grid(rows = vars(Method),cols = vars(n),labeller = label_both) +
  labs(x="p",y="Proportion Contained") +
  theme_bw()

# Save generated datasets -->> Called within file gifBuildScript
### ***CHECK FILENAME TO MAKE SURE NOT TO OVERWRITE ANYTHING*** ###
# save(dataSet,dataPlotContain,dataPlotLength,file="binom95CIsimulation.RData")
