# Subclassification paper analysis
# Ellie Colson
# 2/13/15

rm(list=ls())
library("MatchIt")
library("survey")
library("parallel")
setwd("C:/Users/kecolson/Google Drive/subclassification")

# Generate data function -------------------------------------------------------------

generateData <- function(n, tx.mech = 1, tx.hetero = F, ATT = F) {
  # Create population
  nn <- 1000000
  W1 <- runif(nn, min = 0.02, max = 0.7)  
  W2 <- rnorm(nn, mean = (0.2 + 0.125*W1), sd = 1) 
  W3 <- rnorm(nn, mean = -2, sd = 0.7) 
  W4 <- rbinom(nn, size = 1, p = 0.4) 
  
  if (ATT == F) {
    if (tx.mech == 1) A <- rbinom(nn, 1, plogis(-.5 + 1*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
    if (tx.mech == 2) A <- rbinom(nn, 1, plogis(-.3 + 1*W1 -1.5*I(W1^2) + 1.5*I(W1*W2) - 1.5*W2)) 
    Y.0 <- rnorm(nn, mean=(-.5 + 3*poly(W1,2) - 2*W2), sd=1)
    if (tx.hetero == F) Y.1 <- rnorm(nn, mean=(-.5 + 3*poly(W1,2) - 2*W2 + 2), sd=1)
    if (tx.hetero == T) Y.1 <- rnorm(nn, mean=(-.5 + 2*poly(W1,2) - W2 + 1.5 + 2*poly(W1,2) - W2), sd=1) 
  }
  
  if (ATT == T) {
    if (tx.mech == 1) A <- rbinom(nn, 1, plogis(-.5 + 1*W1 + .1*I(W1^2) + .5*I(W1*W2) - .5*W2))
    if (tx.mech == 2) A <- rbinom(nn, 1, plogis(-.3 + 1*W1 -1.5*I(W1^2) + 1.5*I(W1*W2) - 1.5*W2)) 
    Y.0 <-                     rnorm(nn, mean=(-.5 + 3*poly(W1,2) - 2*W2), sd=1)
    #if (tx.hetero == F) Y.1 <- rnorm(nn, mean=(-.5 + 3*poly(W1,2) - 2*W2 + 2), sd=1) # we shouldn't need to use this one-- if there is no tx hetero, there is no difference between ATE and ATT, so we're not interested in this DGM-- right?
    if (tx.hetero == T) Y.1 <- rnorm(nn, mean=(-.5 + 2*poly(W1,2) -   W2 + 1.5 + 3*poly(W1,2) - 3*W2), sd=1) 
  }
  
  Y <-  ifelse(A==1,Y.1,Y.0)
  
  # Sample from population with probability conditional on W1, W3, W4
  in_sample <- rbinom(nn, size = 1, abs((n/nn) + (n/nn*.25)*(W1 - W3 + W4 - 2.75)))
  
  survey_weights <- 1/abs((n/nn) + (n/nn*.25)*(W1 - W3 + W4 - 2.8))
  
  # NCS-A weights range between about 82 and about 40,000 and are heavily skewed. Replicate using gamma distribution and randomly assign.
  # survey_weights <- rgamma(n, 0.5, 2)*7000+82
  
  sample <- data.frame(W1,W2,A,Y,Y.1,Y.0, survey_weights)
  sample <- sample[in_sample == 1, ]
  sample
}

# DGM 1 ------------------------------------------------------------------------------
# Good overlap 
# No positivity issues
# Without treatment effect heterogeneity

# Examine treatment mechanism, propensity scores, and survey weights 
set.seed(123)
data <- generateData(n = 100000, tx.mech = 1, tx.hetero = F, ATT = F)
g <- glm(A ~ poly(W1,2) + W2 + W1:W2, family = "binomial", data = data)
pscore <- predict(g, type="response")
p1 <- hist(pscore[data$A==1], breaks=30)
p2 <- hist(pscore[data$A==0], breaks=30)
plot(p2, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot(p1, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 
summary(pscore)
hist(data$survey_weights, breaks=30)
summary(data$survey_weights)
table(data$A)
ATE1 <- mean(data$Y.1 - data$Y.0) # True ATE
ATE1
mean(data$Y[data$A==1]) - mean(data$Y[data$A==0]) # Unadjusted 

# Conduct subclassification with 10 subclasses and see how treated and control units are distributed across the classes
match <- suppressWarnings(matchit(A ~ poly(W1,2) + W2 + W1:W2, data = data, method="subclass", subclass = 10, sub.by = "all"))
m <- match.data(match)
table(m$A)
round(as.matrix(table(m$A, m$subclass))/as.matrix(rbind(table(m$subclass), table(m$subclass))),2) # Proportion of each subclass that is control v treated
round(as.matrix(table(m$A, m$subclass))/ rbind(rep(sum(m$A==0),10), rep(sum(m$A),10)),2) # Proportion of total controls (treated) in each subclass


# DGM 2 -------------------------------------------------------------------------------
# Medium overlap
# With positivity issues
# Without treatment effect heterogeneity

# Examine treatment mechanism, propensity scores, and survey weights 
set.seed(123)
data <- generateData(n = 100000, tx.mech = 2, tx.hetero = F, ATT = F)
g <- glm(A ~ poly(W1,2) + W2 + W1:W2, family = "binomial", data = data)
pscore <- predict(g, type="response")
p1 <- hist(pscore[data$A==1], breaks=30)
p2 <- hist(pscore[data$A==0], breaks=30)
plot(p2, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot(p1, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 
summary(pscore)
hist(data$survey_weights, breaks=30)
summary(data$survey_weights)
table(data$A)
ATE2 <- mean(data$Y.1 - data$Y.0) # True ATE
ATE2
mean(data$Y[data$A==1]) - mean(data$Y[data$A==0]) # Unadjusted ATE

# Conduct subclassification with 10 subclasses and see how treated and control units are distributed across the classes
match <- suppressWarnings(matchit(A ~ poly(W1,2) + W2 + W1:W2, data = data, method="subclass", subclass = 10, sub.by = "all"))
m <- match.data(match)
table(m$A)
round(as.matrix(table(m$A, m$subclass))/as.matrix(rbind(table(m$subclass), table(m$subclass))),2) # Proportion of each subclass that is control v treated
round(as.matrix(table(m$A, m$subclass))/ rbind(rep(sum(m$A==0),10), rep(sum(m$A),10)),2) # Proportion of total controls (treated) in each subclass



# DGM 3 ------------------------------------------------------------------------------
# Good overlap 
# No positivity issues
# With treatment effect heterogeneity

# Examine treatment mechanism, propensity scores, and survey weights 
set.seed(123)
data <- generateData(n = 100000, tx.mech = 1, tx.hetero = T, ATT = F)
g <- glm(A ~ poly(W1,2) + W2 + W1:W2, family = "binomial", data = data)
pscore <- predict(g, type="response")
p1 <- hist(pscore[data$A==1], breaks=30)
p2 <- hist(pscore[data$A==0], breaks=30)
plot(p2, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot(p1, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 
summary(pscore)
hist(data$survey_weights, breaks=30)
summary(data$survey_weights)
table(data$A)
ATE3 <- mean(data$Y.1 - data$Y.0) # True ATE
ATE3
mean(data$Y[data$A==1]) - mean(data$Y[data$A==0]) # Unadjusted ATE
(mean(data$Y[data$A==1]) - mean(data$Y[data$A==0]) - ATE2)/ATE2

# Conduct subclassification with 10 subclasses and see how treated and control units are distributed across the classes
match <- suppressWarnings(matchit(A ~ poly(W1,2) + W2 + W1:W2, data = data, method="subclass", subclass = 10, sub.by = "all"))
m <- match.data(match)
table(m$A)
round(as.matrix(table(m$A, m$subclass))/as.matrix(rbind(table(m$subclass), table(m$subclass))),2) # Proportion of each subclass that is control v treated
table(m$A, m$subclass)
round(as.matrix(table(m$A, m$subclass))/ rbind(rep(sum(m$A==0),10), rep(sum(m$A),10)),2) # Proportion of total controls (treated) in each subclass


# DGM 4 ------------------------------------------------------------------------------
# Medium overlap 
# With positivity issues
# With treatment effect heterogeneity

# Examine treatment mechanism, propensity scores, and survey weights 
set.seed(123)
data <- generateData(n = 100000, tx.mech = 2, tx.hetero = T, ATT = F)
g <- glm(A ~ poly(W1,2) + W2 + W1:W2, family = "binomial", data = data)
pscore <- predict(g, type="response")
p1 <- hist(pscore[data$A==1], breaks=30)
p2 <- hist(pscore[data$A==0], breaks=30)
plot(p2, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot(p1, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 
summary(pscore)
hist(data$survey_weights, breaks=30)
summary(data$survey_weights)
table(data$A)
ATE4 <- mean(data$Y.1 - data$Y.0) # True ATE
ATE4
mean(data$Y[data$A==1]) - mean(data$Y[data$A==0]) # Unadjusted ATE
(mean(data$Y[data$A==1]) - mean(data$Y[data$A==0]) - ATE2)/ATE2

# Conduct subclassification with 10 subclasses and see how treated and control units are distributed across the classes
match <- suppressWarnings(matchit(A ~ poly(W1,2) + W2 + W1:W2, data = data, method="subclass", subclass = 10, sub.by = "all"))
m <- match.data(match)
table(m$A)
round(as.matrix(table(m$A, m$subclass))/as.matrix(rbind(table(m$subclass), table(m$subclass))),2) # Proportion of each subclass that is control v treated
table(m$A, m$subclass)
round(as.matrix(table(m$A, m$subclass))/ rbind(rep(sum(m$A==0),10), rep(sum(m$A),10)),2) # Proportion of total controls (treated) in each subclass


# DGM 5 ------------------------------------------------------------------------------
# Good overlap 
# No positivity issues
# With treatment effect heterogeneity
# For ATT

# Examine treatment mechanism, propensity scores, and survey weights 
set.seed(123)
data <- generateData(n = 100000, tx.mech = 1, tx.hetero = T, ATT = T)
g <- glm(A ~ poly(W1,2) + W2 + W1:W2, family = "binomial", data = data)
pscore <- predict(g, type="response")
p1 <- hist(pscore[data$A==1], breaks=30)
p2 <- hist(pscore[data$A==0], breaks=30)
plot(p2, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot(p1, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 
summary(pscore)
hist(data$survey_weights, breaks=30)
summary(data$survey_weights)
table(data$A)
mean(data$Y[data$A==1]) - mean(data$Y[data$A==0]) # Unadjusted
ATT5 <- mean(data$Y.1[data$A==1]) - mean(data$Y.0[data$A==1])
ATT5

# Conduct subclassification with 10 subclasses and see how treated and control units are distributed across the classes
match <- suppressWarnings(matchit(A ~ poly(W1,2) + W2 + W1:W2, data = data, method="subclass", subclass = 10, sub.by = "all"))
m <- match.data(match)
table(m$A)
round(as.matrix(table(m$A, m$subclass))/as.matrix(rbind(table(m$subclass), table(m$subclass))),2) # Proportion of each subclass that is control v treated
table(m$A, m$subclass)

# z.out3 <- zelig(Y ~ A, data = match.data(match, "treat"), weights="survey_weights", model = "ls", by = "subclass")
# z.out3 <- zelig(Y ~ A, data = m, weights = "survey_weights", model = "ls", by = "subclass")
# z.out3 <- zelig(Y ~ A, data = m, weights = "weights", model = "ls", by = "subclass")
# z.out3 <- zelig(Y ~ A, data = m, model = "ls", by = "subclass")
# z.out3 <- zelig(Y ~ A, data = match.data(match, "treat"), model = "ls", by = "subclass")
# z.out3
# 
# zout <- zelig(Y ~ A, data = match.data(match), model = "ls", by = "subclass")
# xout <- setx(zout, data = match.data(match, "treat"))
# sout <- sim(zout, x = xout)
# summary(sout)
# 
# mod <- glm(Y ~ A, data = m[m$subclass == 2,])
# summary(mod)
# pred <- predict(mod, se.fit = T, type = "response")
# table(pred$se.fit)
# pred <- predict(mod, newdata = m[m$A==1 & m$subclass==2,], se.fit = T, type = "response")
# table(pred$se.fit)
# 
# design <- svydesign(~ 0, 
#                         data = m[,c("W1","W2","A","Y")], 
#                         weights = m$survey_weights)
# mod <- svyglm(Y ~ A, design = subset(design, m$subclass == 1))
# summary(mod)
# pred <- predict(mod, se.fit = TRUE, type = "response")
# table(SE(pred))
# 
# pred <- predict(mod, newdata = m[m$A==1,], se.fit = T)
# table(sqrt(vcov(pred)))
# 
# pred <- predict(mod, newdata = m[m$A==0,], se.fit = T)
# table(sqrt(vcov(pred)))

# DGM 6 ------------------------------------------------------------------------------
# Medium overlap 
# With positivity issues
# With treatment effect heterogeneity
# For ATT

# Examine treatment mechanism, propensity scores, and survey weights 
set.seed(123)
data <- generateData(n = 100000, tx.mech = 2, tx.hetero = T, ATT = T)
g <- glm(A ~ poly(W1,2) + W2 + W1:W2, family = "binomial", data = data)
pscore <- predict(g, type="response")
p1 <- hist(pscore[data$A==1], breaks=30)
p2 <- hist(pscore[data$A==0], breaks=30)
plot(p2, col=rgb(0,0,1,1/4), xlim=c(0,1)) 
plot(p1, col=rgb(1,0,0,1/4), xlim=c(0,1), add=T) 
summary(pscore)
hist(data$survey_weights, breaks=30)
summary(data$survey_weights)
table(data$A)
mean(data$Y[data$A==1]) - mean(data$Y[data$A==0]) # Unadjusted ATE
ATT6 <- mean(data$Y.1[data$A==1]) - mean(data$Y.0[data$A==1])
ATT6

# Conduct subclassification with 10 subclasses and see how treated and control units are distributed across the classes
match <- suppressWarnings(matchit(A ~ poly(W1,2) + W2 + W1:W2, data = data, method="subclass", subclass = 10, sub.by = "all"))
m <- match.data(match)
table(m$A)
round(as.matrix(table(m$A, m$subclass))/as.matrix(rbind(table(m$subclass), table(m$subclass))),2) # Proportion of each subclass that is control v treated
table(m$A, m$subclass)


## Analysis function -------------------------------------------------------

analyze <- function(iteration, tx.mech = 1, tx.hetero = F, ATT = F, cor.spec.pscore = T, n.subclass = 10, svy.wt = F, n = 3000) {
  library("MatchIt")
  set.seed(iteration)
  print(iteration)
  sample <- generateData(n = n, tx.mech = tx.mech, tx.hetero = tx.hetero, ATT = ATT)
  sub.by <- ifelse(ATT, "treat", "all")
  if (cor.spec.pscore == T) {
    match <- suppressWarnings(matchit(A ~ I(poly(W1,2)) + W2 + W1:W2, 
                                      data = sample, 
                                      sub.by = sub.by,
                                      method = "subclass", 
                                      subclass = n.subclass))
  } else {
    match <- suppressWarnings(matchit(A ~ W1 + W2, 
                                      data = sample, 
                                      sub.by = sub.by,
                                      method = "subclass", 
                                      subclass = n.subclass))
  }
  Obs <- match.data(match)
  
  ests <- var <- rep(NA, n.subclass)
   
  if (svy.wt == T) {
    design <- svydesign(~ 0, 
                        data = Obs[,c("W1","W2","A","Y")], 
                        weights = Obs$survey_weights)
    for (i in 1:n.subclass) {
      try(model <- svyglm(Y ~ A, design = subset(design, Obs$subclass == i)))
      if (sum(Obs$A[Obs$subclass == i])>0 & sum(Obs$A[Obs$subclass==i]==0)>0) {
        ests[i] <-  summary(model)$coefficients[2,1] 
        if (ATT == F) var[i]  <- (summary(model)$coefficients[2,2])^2
      }
    }
  } else {
    for (i in 1:n.subclass) {
      try(model <- glm(Y ~ A, data = Obs, subset = (Obs$subclass == i)))
      if (sum(Obs$A[Obs$subclass == i])>0 & sum(Obs$A[Obs$subclass==i]==0)>0) {
        ests[i] <-  summary(model)$coefficients[2,1] 
        if (ATT == F) var[i]  <- (summary(model)$coefficients[2,2])^2
      }
    }
  }
  
  # Get variance estimates for ATT
  if (ATT == T) {
    # Estimate the propensity score for each observation
    mod <- glm(A ~ I(poly(W1,2)) + W2 + W1:W2, data = Obs, family = "binomial")
    e <- predict(mod, type = "response")
    # Estimate the predicted outcome if tx = 1 and if tx = 0 using the true parametric model
    mod <- glm(Y ~ I(poly(W1,2)) + W2 + A + A:I(poly(W1,2)) + A:W2, data = Obs)
    tx <- ct <- Obs
    tx$A <- 1; ct$A <- 0
    mu1 <- predict(mod, newdata = tx, type = "response")
    mu0 <- predict(mod, newdata = ct, type = "response")
    # Estimate EIC for each observation within each subclass
    # e(x): the propensity score. This is specific to the observation.
    # w: treatment. 0=no, 1=yes. This is also specific to the observation.
    # y: outcome. This is also specific to the observation.
    # \mu_1 is the predicted outcome if tx=1. \mu_0 is the predicted outcome if tx=0. These values are specific to the observation. We would use the true parametric outcome model (correct variables and correct form) to get these predicted values.
    # \tau_P would be the stratum-specific ATT
    # EIC = ((w*y/e) - ((1-w)*y/(1-e))) - tau - (mu1/e + mu0/(1-e))*(w - e)
    w <- Obs$A; y <- Obs$Y; tau <- ests[Obs$subclass]
    EIC <- ((w*y/e) - ((1-w)*y/(1-e))) - tau - (mu1/e + mu0/(1-e))*(w - e)
    # Estimate variance for treated observations within each subclass:  get the sample variance of those in the stratum-specific treatment group by calculating the variance of those in the treatment group in stratum a and dividing by the number in the treatment group in stratum a.
    var <- sapply(1:n.subclass, function(x) var(EIC[Obs$subclass==x & Obs$A==1])/nrow(Obs[Obs$subclass==x & Obs$A==1,]) )
    
  }
  
  if (ATT == F) combo1 <- weighted.mean(ests[!is.na(ests)], w = as.numeric(table(Obs$subclass)[!is.na(ests)]))
  if (ATT == T) combo1 <- weighted.mean(ests[!is.na(ests)], w = as.numeric(table(Obs$subclass[Obs$A==1])[!is.na(ests)]))
  combo2 <- weighted.mean(ests[!is.na(ests)], w = 1/var[!is.na(ests)])
  c(combo1, combo2)
}

MSE.funct<- function(X, Xtrue) {
  ave <- mean(X, na.rm = T)    
  bias <- mean(X- Xtrue, na.rm = T)
  var <- var(X, na.rm = T)
  MSE <- mean( (X-Xtrue)^2, na.rm = T)
  c(ave, bias, var, MSE)
}


# Check results for each, run once
set.seed(123); analyze(iteration = 1, tx.mech = 1, tx.hetero = T, ATT = T, cor.spec.pscore = T, n.subclass = 30,  svy.wt = F, n = 5000)
ATT5

# Run simulations---------------------------------------------------------------

results <- expand.grid(n = 2000, combo = c(1,2), n.subclass = c(5,10), svy.wt = c(F,T), tx.mech = c(1,2), tx.hetero = T, ATT = T, cor.spec.pscore = c(T,F), mean = NA, bias = NA, Var = NA, MSE = NA)

results[,13:1012] <- NA

index <- 1:1000

cl <- makeCluster(4)
ex <- Filter(function(x) is.function(get(x, .GlobalEnv)), ls(.GlobalEnv)) # send Cluster user defined fns 
clusterExport(cl, ex)


for (j in seq(1,nrow(results),2)) {
  print(j)
  r <- parLapply(cl, index, analyze, tx.mech = results$tx.mech[j], 
                 tx.hetero = results$tx.hetero[j],
                 ATT = results$ATT[j],
                 cor.spec.pscore = results$cor.spec.pscore[j], 
                 n.subclass = results$n.subclass[j], 
                 svy.wt = results$svy.wt[j], n = results$n[j])
  r <- do.call(rbind, r)
  if (results$ATT[j] == F & results$tx.hetero[j] == F) Xtrue <- 2
  if (results$ATT[j] == F & results$tx.hetero[j] == T) Xtrue <- ATE3
  if (results$ATT[j] == T & results$tx.hetero[j] == T & results$tx.mech[j] == 1) Xtrue <- ATT5
  if (results$ATT[j] == T & results$tx.hetero[j] == T & results$tx.mech[j] == 2) Xtrue <- ATT6
  results[j:(j+1),9:12] <- t(apply(r, 2, MSE.funct, Xtrue = Xtrue))
  results[j,13:1012] <- r[,1]
  results[(j+1),13:1012]   <- r[,2]
}

stopCluster(cl)

write.csv(results, "results/2015_03_10_ATT_5_10.csv", row.names = F)


