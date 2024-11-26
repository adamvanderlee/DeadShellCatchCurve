################################################################################
#        1         2         3         4         5          6        7         8                   
#2345678901234567890123456789012345678901234567890123456789012345678901234567890
################################################################################
# CATCH-CURVE FROM DEAD SHELLS SIMULATIONS
# conduct simulation to test if/when you can get accurate estimates of mortality
# when conducting a catch-curve analysis on dead Unionid shells
#-------------------------------------------------------------------------------

# Clear workspace
rm(list = ls())  ### clear the workspace

# set working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# load libraries
library(pacman)
p_load(FSA)
p_load(data.table)
p_load(parallel)
p_load(ggplot2)
p_load(fitdistrplus)
p_load(simMetric)

#-------------------------------------------------------------------------------
# ggplot theme

theme_me <- theme_bw() +
  theme(axis.title = element_text(size = 11, family = "sans", face = "bold"),
        axis.text.x = element_text(size = 10, family = "sans", colour = "black"),
        axis.text.y = element_text(size = 10, family = "sans", hjust = 0.6,
                                   angle = 0, colour = "black"),
        legend.title = element_text(size = 10, family = "sans"),
        legend.text = element_text(size = 8, family = "sans"),
        strip.text = element_text(size = 11, family = "sans", face = "bold"))

#-------------------------------------------------------------------------------

# Function to generate shell remaining in system for recently deceased mussels
shell_remain_f <- function(Tmax, # Longevity
                           rate, # maximum rate shell can remain
                           type  # shape of relationship between rate and age
                                 # constant, increasing, hockey-stick, Logistic
){
  
  if(type == "constant") {          # constant remain rate
    x <- rep(rate, Tmax)
  } else if(type == "increasing"){  # linear increase
    x <- 1:(Tmax)*(rate/Tmax)
  } else if(type == "hockey_stick"){# linear increase to 50% Tmax, then constant
    x <- 1:(Tmax)*(rate/(0.5*Tmax))
    x <- ifelse(x>rate, rate, x)
  } else if(type == "logistic") {   # shape increase, inflection point at 50% Tmax
    x <- rate/(1 + exp(-0.5*(1:Tmax - Tmax/2)))
  }
  
  x
  
}

# Plot Shell remain functions
plot_data <- do.call(rbind, 
  lapply(c("constant", "increasing", "hockey_stick", "logistic"),
    function(type) {
       out <- data.table(Age = 1:30,
                         R_rate = shell_remain_f(Tmax = 30,
                                                 rate = 0.2,
                                                 type = type),
                         type = type)
                      
}))
plot_data$type <- factor(plot_data$type,
                         levels = c("constant", "increasing", "hockey_stick", 
                                    "logistic", "Live Individuals"),
                         labels = c("Constant",  "Increasing", "Hockey Stick",
                                    "Logistic", "Live Individuals")) 

png("Figures/Figure1.png", width = 3.5, height = 2.75, unit = "in", res = 300)
ggplot(plot_data)+
  geom_line(aes(Age, R_rate, colour = type), linewidth = 1.5)+
  labs(x="Age", y = "Shell Persistence Rate", colour = NULL) +
  theme_me + theme(legend.position = "inside",
                   legend.position.inside = c(0.99,0.01), 
                   legend.justification = c(1,0))
dev.off()

# Clean up
rm(plot_data)

#-------------------------------------------------------------------------------
# Run Simulation
# Estimate Z using catch curve analysis from empty shells of dead unionids
# Draw Z from uniform dist @ 0.05 increments - S ~ 40$ to 95%
# estimate Tmax from Z - add error
# draw remain rate (max) and decay rate from uniform dist.
# select type of remain shape - constant, increasing, hockey-stick, logistic
# generate stable stage distribution and abundance vector
# run simulation to generate shell assemblage available form sampling
# sample the shell assemblage for sample size N - 25, 50, 100, 200
# estimate Z with catch curve - Chapman-Robson and weighted regression
# repeat catch curve with sample from live population. 

# Number of replicates
reps <- 100 # per Z, type, & N combo

# parameter combinations - loop over with outside loop
params = tidyr::crossing(
  Z_true = seq(0.05, 0.95, 0.05), 
  type = c("constant", "increasing","hockey_stick", "logistic"),
  n_sample = c(25, 50, 100, 200))

# Simulations
no_cores <- detectCores() - 1                              # number of cores
cl <- makeCluster(no_cores)                                # create clusters
clusterExport(cl, list("shell_remain_f", "params", "reps"))# send data to clusters
clusterEvalQ(cl, {library(data.table); library(FSA)})      # load library in clusters
clusterSetRNGStream(cl, iseed = 123456)
CC_sim <- do.call(rbind, parLapply(cl, 1:nrow(params), function(p){
  
  # parameters 
  Z_true <- params$Z_true[p]     # true mortality rate
  type <- params$type[p]         # remain rate shape
  n_sample <- params$n_sample[p] # sample size for catch curve
  
  # replicates over parameters combination
  do.call(rbind, lapply(1:reps, function(r){
  #
    # LIFE-HISTORY/POPULATION DATA
    #
    K <- 2000 # population size - small b/c SAR
    
    # sample mortality rate - over approx 40%-95% suvival
    S_rate <- exp(-Z_true)
  
    # generate Tmax - predicted from Haag 2012 +/- 15% noise for variation
    Tmax <- round(rnorm(1, 
                        (Z_true/4.171)^(1/-1.070), 
                        (Z_true/4.171)^(1/-1.070)*0.15))
    Tmax <- ifelse(Tmax < 5, 5, Tmax)
    
    R_rate <- runif(1, 0.1, 0.5)  # shell remain maximal rate
    D_rate <- runif(1, 0.05, 0.3) # shell decay rate - based on Strayer and Malcom (2007) Figure 1c
  
    # Stable Stage structure
    SS <- S_rate^(1:Tmax)/sum(S_rate^(1:Tmax))
    
    # Population vector
    # include demographic stochasticity so theres a chance of older age-classes at
    # low population size
    N <- sapply(SS, function(ss) sum(rbinom(K, 1, ss)))
    
    #
    # SHELL SIMULATION
    #
    # Generate age-distrubtion of shells remaining in the system
    # run for 500 years to stabilize structure
    Shells <- rep(0, Tmax) # inititalize vector
    for(i in 1:500) {
      
      # Shells decay 
      Shells <- round(Shells * exp(-D_rate))
      
      # Number of dead shells 
      # demographic sotchasitic by age determines number that die
      Dead <- sapply(N, function(n) sum(rbinom(n, 1, 1-S_rate)))
      
      # Increment that population vector
      # Remove dead, move to next age-class; add new recruits - based on SS prop.
      N <- N - Dead                    # Surviving individuals
      N <- c(sum(rbinom(K, 1, SS[1])), # new recruits
             N[1:(Tmax - 1)])          # transition to next age (2:Tmax)
      
      # Shells that remain in the system
      shell_remain_rate <- shell_remain_f(Tmax, R_rate, type) 
      Shells <- Shells + round(Dead * shell_remain_rate)
      
    }
    
    #
    # SAMPLE SHELLS
    #
    # Take a sample of shell for population
  
    # Ensure no. of shells is > sample size
    if(sum(Shells) < n_sample) {
      
      catch_curve_CH_est <- 
        data.table(
          S_est_CH = NA,
          S_err_CH = NA,
          Z_est_CH = NA,
          Z_err_CH = NA
        )
      
      catch_curve_CH_est_N <- 
        data.table(
          S_est_CH_N = NA,
          S_err_CH_N = NA,
          Z_est_CH_N = NA,
          Z_err_CH_N = NA
        )
      
      catch_curve_CC_est <- 
        data.table(
          S_est_R = NA,
          Z_est_R = NA
        )
      
      catch_curve_CC_est_N <- 
        data.table(
          S_est_R_N = NA,
          Z_est_R_N = NA
        )
      
    } else {
    
      # sample shells
      # create vector of shell ages - length No. shells
      # table sample of n_sample from vector of shell ages
      shell_sample <- sample(rep(1:Tmax, Shells), n_sample)
      
      # Count number at each age
      sample_tab <- table(shell_sample)
      
      # organize into data frame
      shell_data <- data.table(
        Age = as.numeric(names(sample_tab)), # age
        n = as.vector(sample_tab))    # no. of shells at age
      
      # sample live individuals
      # create vector of shell ages - length No. shells
      # table sample of n_sample from vector of shell ages
      N_sample <- sample(rep(1:Tmax, N), n_sample)
      
      # Count number at each age
      sample_tab <- table(N_sample)
      
      # organize into data frame
      N_data <- data.table(
        Age = as.numeric(names(sample_tab)), # age
        n = as.vector(sample_tab))    # no. of shells at age
      
      #
      # CATCH CURVE ANALYSIS - Shells
      #
      # organize data
      
      # ID min age to include
      Age_min <- shell_data$Age[which.max(shell_data$n)]
      Age_min <- ifelse(length(Age_min) == 0, 1, Age_min)
      
      # Chapman-Robson Catch Curve
      if(Tmax - (Age_min + 1) < 3) {
        catch_curve_CH_est <- 
          data.table(
            S_est_CH = NA,
            S_err_CH = NA,
            Z_est_CH = NA,
            Z_err_CH = NA
          )
      } else {
        catch_curve_CH_est <- tryCatch({
          CC <- chapmanRobson(n~Age, data = shell_data, 
                        ages2use = shell_data[Age>=(Age_min+1)]$Age)
          data.table(
            S_est_CH = CC$est[1,1]/100, # survival rate
            S_err_CH = CC$est[1,2]/100, # survival rate SE
            Z_est_CH = CC$est[2,1],     # Instantaneous Mortality
            Z_err_CH = CC$est[2,2]      # Instantaneous Mortality SE
          )
        }, error = function(err) {      # If error fill in with NAs
          data.table(
            S_est_CH = NA,
            S_err_CH = NA,
            Z_est_CH = NA,
            Z_err_CH = NA
          )
        })
      }
      
      # Regression base catch curve
      catch_curve_CC_est <- tryCatch({
        CC <- catchCurve(n~Age, data = shell_data, 
                      ages2use = shell_data[Age>=Age_min]$Age,
                      weighted = TRUE)
        data.table(
          S_est_R = exp(CC$lm$coefficients[2]), # Survival rate
          Z_est_R = -CC$lm$coefficients[2]      # Instantaneous Mortality
        )
        
      }, error = function(err) { # If error fill in with NAs
        data.table(
          S_est_R = NA,
          Z_est_R = NA
        )
      })
      
      #
      # CATCH CURVE ANALYSIS - N
      #
      # organize data
      
      # ID min age to include
      Age_min <- N_data$Age[which.max(N_data$n)]
      Age_min <- ifelse(length(Age_min) == 0, 1, Age_min)
      
      # Chapman-Robson Catch Curve
      if(Tmax - (Age_min + 1) < 3) {
        catch_curve_CH_est_N <- 
          data.table(
            S_est_CH_N = NA,
            S_err_CH_N = NA,
            Z_est_CH_N = NA,
            Z_err_CH_N = NA
          )
      } else {
        catch_curve_CH_est_N <- tryCatch({
          CC <- chapmanRobson(n~Age, data = N_data, 
                              ages2use = N_data[Age>=(Age_min+1)]$Age)
          data.table(
            S_est_CH_N = CC$est[1,1]/100, # survival rate
            S_err_CH_N = CC$est[1,2]/100, # survival rate SE
            Z_est_CH_N = CC$est[2,1],     # Instantaneous Mortality
            Z_err_CH_N = CC$est[2,2]      # Instantaneous Mortality SE
          )
        }, error = function(err) {      # If error fill in with NAs
          data.table(
            S_est_CH_N = NA,
            S_err_CH_N = NA,
            Z_est_CH_N = NA,
            Z_err_CH_N = NA
          )
        })
      }
      
      # Regression base catch curve
      catch_curve_CC_est_N  <- tryCatch({
        CC <- catchCurve(n~Age, data = N_data, 
                         ages2use = N_data[Age>=Age_min]$Age,
                         weighted = TRUE)
        data.table(
          S_est_R_N = exp(CC$lm$coefficients[2]), # Survival rate
          Z_est_R_N = -CC$lm$coefficients[2]      # Instantaneous Mortality
        )
        
      }, error = function(err) { # If error fill in with NAs
        data.table(
          S_est_R_N = NA,
          Z_est_R_N = NA
        )
      })
    } 
    #
    # OUTPUT
    #
    data.table(
      K = K,                # Population Size
      Shells = sum(Shells), # No. of Shells in environment
      N = n_sample,         # Sample Size
      Tmax = Tmax,          # Longevity
      Z_true,               # True instantaneous Mortality
      S_rate = S_rate,      # Survival rate - actual
      R_rate = R_rate,      # Remain rate - maximum
      D_rate = D_rate,      # Decay Rate
      type = type,          # Remain curve shape
      catch_curve_CH_est,   # Chapman-robson catch curve - shells
      catch_curve_CC_est,   # Regrssion catch curve - shells
      catch_curve_CH_est_N, # Chapman-Robsen catch curve - live individuals
      catch_curve_CC_est_N  # Regression catch curve - live individuals
    )
  }))
}))
stopCluster(cl)

#-------------------------------------------------------------------------------

# Exclude samples size < 25
CC_sim <- CC_sim[N>=25, ]
CC_sim <- na.omit(CC_sim)
CC_sim <- CC_sim[Z_true<=1, ]

# Re-organize data frame to long
CC_sim_long <- rbind(
  with(CC_sim, data.table(
    Z_true = c(Z_true, Z_true),
    Z_est = c(Z_est_CH, Z_est_R),
    S_est = c(S_est_CH, S_est_R),
    K = c(K, K),
    Shells = c(Shells, Shells),
    N = c(N, N),
    Tmax = c(Tmax, Tmax),
    R_rate = c(R_rate, R_rate),
    D_rate = c(D_rate, D_rate),
    type = c(type, type),
    method = c(rep("Chapman-Robson", nrow(CC_sim)),
               rep("Regression", nrow(CC_sim)))
  )),
  with(CC_sim[sample(1:nrow(CC_sim), nrow(params)*reps/4),], data.table(
    Z_true = c(Z_true, Z_true),
    Z_est = c(Z_est_CH_N, Z_est_R_N),
    S_est = c(S_est_CH_N, S_est_R_N),
    K = c(K, K),
    Shells = c(Shells, Shells),
    N = c(N, N),
    Tmax = c(Tmax, Tmax),
    R_rate = c(R_rate, R_rate),
    D_rate = c(D_rate, D_rate),
    type = "Live Individuals",
    method = c(rep("Chapman-Robson", nrow(params)*reps/4),
               rep("Regression", nrow(params)*reps/4))
  ))
)
CC_sim_long$Tmax_cat <- round(CC_sim_long$Tmax/5)*5
CC_sim_long$type <- factor(CC_sim_long$type,
                           levels = c("constant",  "increasing", "hockey_stick",
                                      "logistic", "Live Individuals"),
                           labels = c("Constant", "Increasing", "Hockey Stick", 
                                      "Logistic", "Live Individuals")) 
# sort table
setorder(CC_sim_long, method, type, N, Z_true)

#-------------------------------------------------------------------------------
# Sample Size

CC_sim[, .("propN" = .N/7600), by = .(type)]
CC_sim[, .("propN" = .N/7600), by = .(N)]
N_data <- CC_sim[, .("propN" = .N/1600), by = .(Z_true)]
setorder(N_data,  Z_true); N_data


N_data <- CC_sim[, .("n" = .N, "propN" = .N/100), by = .(type,N,Z_true)]
setorder(N_data, type, N, Z_true)


#-------------------------------------------------------------------------------
# Calculate accuracy metrics to assess performance of CC estimates
# Follows Morris et al (2019) Using simulation studies to evaluate statistical
#   methods. Stats. in Medicine. 38:2074-2102. 
# Bias - deviation of estimate from truth
# Empirical SE - measure of precision
# Root mean square error - contains elements of bias ans precision, in units
#                          of estimate

# Overall Accuracy
accuracy_data_Z <- rbind(
  CC_sim_long[, list("estimate" = simMetric::bias(Z_true, Z_est, 
                                                  get = "bias"),
                     "MCSE" = simMetric::bias(Z_true, Z_est, 
                                              get = "bias_mcse"),
                     "n" = .N,
                     "Estimator" = "Bias"), 
              by = list(Z_true, N, type, method)],
  CC_sim_long[, list("estimate" = simMetric::empSE(Z_est, 
                                    get = "empSE"),
                     "MCSE" = sqrt(simMetric::empSE(Z_est, 
                                                    get = "empSE_mcse")),
                     "n" = .N,
                     "Estimator" = "Emp. SE"), 
              by = list(Z_true, N, type, method)],
  CC_sim_long[, list("estimate" = sqrt(simMetric::mse(Z_true, Z_est,
                                                      get = "mse")),
                     "MCSE" = sqrt(simMetric::mse(Z_true, Z_est, 
                                                  get = "mse_mcse")),
                     "n" = .N,
                     "Estimator" = "RMSE"), 
              by = list(Z_true, N, type, method)]
)

#-------------------------------------------------------------------------------

# Expected vs estimated survival
png("Figures/Figure2.png", width = 8.5, height = 4.5, unit = "in", res = 300)
ggplot(CC_sim_long,aes(Z_true, Z_est))+
  geom_point(aes(colour = as.factor(N)), alpha = 0.5, size = 0.75) +
  geom_abline(intercept = 0, slope = 1, size = 1.0) +
  geom_smooth(method = loess, se = FALSE) +
  labs(x = "True Z", y = "Estimated Z", colour = "Sample Size")+
  facet_grid(method~type) +
  theme_me + theme(legend.position = "top")+
  coord_cartesian(ylim=c(0,2))
dev.off()

# Accuracy metrics across shell loss shapes and CC estimate methods by sample
# size
png("Figures/Figure3.png", width = 8.5, height = 6, unit = "in", res = 300)
ggplot(accuracy_data_Z, aes(Z_true, estimate))+
  geom_point(aes(colour = as.factor(N), shape = method)) +
  geom_hline(yintercept=0)+
  geom_smooth(aes(linetype = method, colour = as.factor(N)),
              method = loess, se = FALSE) +
  labs(x = "True Z", y = "Metric Value", colour = "Sample Size",
       linetype = "Method", shape = "Method")+
  facet_grid(Estimator~type, scales = "free_y") +
  theme_me + theme(legend.position = "top")
dev.off()

# Sample size of simulation across: Z, method and type. - not for manuscript
ggplot(N_data, aes(Z_true, propN))+
  geom_line(aes(colour = as.factor(N))) +
  labs(x = "Expected Z", y = "Performace", colour = "Sample Size",
       linetype = "Estimatation Method", shape = "Estimatation Method")+
  facet_grid(.~type, scales = "free_y") +
  theme_me + theme(legend.position = "top")

# Mean RMSE across Z and N
accuracy_data_Z[Estimator == "RMSE", 
                mean(estimate), 
                by = .(type,method, Estimator)]

#-------------------------------------------------------------------------------
################################################################################