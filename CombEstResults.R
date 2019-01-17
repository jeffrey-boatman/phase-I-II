
################################################################################

### Purpose: Take the survival and toxcity simulated data and obtain new 
### estimates via cohort combining methods




######################## POOLED ESTIMATOR FUNCTION #####################

# PURPOSE: Takes the escalation data at the optimal dose and expansion data then 
# blindly combine them to get efficient estimation for efficacy and toxicity. 
# Although efficicent in estimation (i.e. using all the data), it may be biased.

# INPUT: 
## d.esl: Escalation death data from dose escalation trial
## d.exp: Simulated expansion death data
## t.esl: Escalation toxicity data from dose escalation trial
## t.exp: Simulated expansion toxicity data

# OUTPUT (list):
## pool: a vector of death and toxicity indicators using both escalation and 
# expansion data


Pool <- function(d.esl, d.exp, t.esl, t.exp){
	d.pool <- c(d.esl, d.exp)
	t.pool <- c(t.esl, t.exp)
	s.conf <- binconf(sum(1-d.pool),length(1-d.pool), method = "wilson")[2:3]
	t.conf <- binconf(sum(t.pool),length(t.pool), method = "wilson")[2:3]
	pooled <- c(mean(1-d.pool), s.conf, mean(t.pool), t.conf)
	return(pooled)
}


########## Rao-Blackwell Estimator Function ##############


# PURPOSE: To take the simulated data (Dose ESL and EXP) for efficacy and 
# toxicity and create a permutation distribution. Then rerun it through the 
# dose finding model and finding the average of the proportion for both efficacy 
# and toxicity for permutations where the original optimal dose equals the newly 
# found optimal dose (i.e. from the permutated data).

# INPUT: 
## d.esl: Escalation death data from dose escalation trial
## d.exp: Simulated expansion death data
## t.esl: Escalation toxicity data from dose escalation trial
## t.exp: Simulated expansion toxicity data
## dose: Dose assignments from original dose finding model
## optimal_dose: Optimal dose from the dose escalation trial
## tox_time_obs: Simulated toxicity observed time from original dose 
##  finding model
## surv_time_obs: Simulated survival observed time from original dose 
##  finding model

# OUTPUT: 
## Conditional UMVUE E[x|(permuted) expansion cohort*] for survival and toxicity

RaoBlackwell <- function(d.esl, d.exp, t.esl, t.exp, dose, optimal_dose, 
  tox_time_obs, surv_time_obs){

	# Length of death & toxicity indicators for the optimal_dose
	m <- length(d.esl[dose == optimal_dose])
	
	# Combined escalation and expansion death & toxicity vectors
	d.all <- c(d.esl[dose == optimal_dose], d.exp)
	t.all <- c(t.esl[dose == optimal_dose], t.exp)
	
	# Of all possible permutations of death indicator, a random sample (of size q) 
	# of the combined escalation & expansion
	# is taken of size (m + expansion size)
	q <- 100
	d.perm <- t(replicate(n = q, d.all[sample(1:length(d.all), replace=FALSE)]))
	t.perm <- t(replicate(n = q, t.all[sample(1:length(t.all), replace=FALSE)]))
	
	## Non-parametric Bootsrap: We will resample each permutation with replacement 
	## from the observed data from 
	## both the escalation trial and follow-up cohort. The quantiles for the 
	## distribution will be used for CI
	# d.npb <- lapply(1:q, function(i) replicate(n = q, d.perm[i,][sample(1:length(d.perm[i,]), replace = TRUE)]))
	# t.npb <- lapply(1:q, function(i) replicate(n = q, t.perm[i,][sample(1:length(t.perm[i,]), replace = TRUE)]))
	# cov <- matrix(NA, ncol = 4, nrow = q)
	# for (k in 1:q){
		
		# # Obtain each resampled matrix of death and toxicity indicators
		# d.resamp <- t(d.npb[[k]])
		# t.resamp <- t(t.npb[[k]])
		
		# # Replace the resamples for the optimal dose with the original doses (for toxicity and survival)
		# d.list <- t(sapply(1:q, function(i) replace(x = d.esl, dose == optimal_dose, d.resamp[i, 1:m])))
		# t.list <- t(sapply(1:q, function(i) replace(x = t.esl, dose == optimal_dose, t.resamp[i, 1:m])))
		
		# # Rerun samples through dose model to find optimal dose
		# dose.mod <- sapply(1:q, function(i) DOSE.MODEL(d.list[i,], t.list[i,], dose, tox_time_obs, surv_time_obs))
		
		# s.new <- sapply(1:q, function(i) ifelse(dose.mod[i] == optimal_dose, 1 - mean(d.list[i, m:length(d.list[i,])]), NA))
		# t.new <- sapply(1:q, function(i) ifelse(dose.mod[i] == optimal_dose, mean(t.list[i, m:length(t.list[i,])]), NA))
		
		# s.cov <- quantile(s.new, c(0.025, 0.975), na.rm = T)
		# t.cov <- quantile(t.new, c(0.025, 0.975), na.rm = T)
		# cov[k,] <- c(s.cov, t.cov)
		
	# }
	

	
	# Replace the new permuted death indicators for optimal dose with the orignal doses (for toxicity and survival)
	d.mat <- t(sapply(1:q, function(i) replace(x = d.esl, dose == optimal_dose, d.perm[i,1:m])))
	t.mat <- t(sapply(1:q, function(i) replace(x = t.esl, dose == optimal_dose, t.perm[i,1:m])))
	
 	# Using the dose.model function, input the death vector and toxicity vector for each q random permutation
 	dose.update <- sapply(1:q, function(i) DOSE.MODEL(d.mat[i,],t.mat[i,], dose, tox_time_obs, surv_time_obs))
		
	# Conditional UMVUE: Expectation of the sample proportion using the optimal dose expansion data conditional on
	# the permuted death outcomes for the optimal dose (escalation and expansion) and the refit model choosing 		
	# the optimal dose.  
	
	s.rb <- sapply(1:q, function(i) ifelse(dose.update[i]==optimal_dose, 1 - mean(d.perm[i, m:length(d.perm[i,])]), NA))
	
	t.rb <- sapply(1:q, function(i) ifelse(dose.update[i]==optimal_dose, mean(t.perm[i, m:length(t.perm[i,])]), NA))
	
	RB.est <- c(mean(s.rb , na.rm=TRUE), mean(t.rb , na.rm=TRUE))
	return(RB.est)
}



####### DOSE FINDING ALGORITHM FOR PERMUTATION FUNCTION #################

# PURPOSE: To be used in PERM function to find new optimal dose, d* from 
# permuted data

# INPUT:
## death: New permuted escalation death vector
## tox: New permuted escalation toxicity vector
## dose: Dose assignments from original dose finding model
## tox_time_obs: Simulated toxicity observed time from original 
##   dose finding model
## surv_time_obs: Simulated survival observed time from original dose 
##   finding model

# OUTPUT: 
## optimal_dose.perm: a vector of new Optimal dose based on the permuted 
## surivial and toxicity data

DOSE.MODEL <- function(death, tox, dose, tox_time_obs, surv_time_obs) {
		# Set total subjects to 21 subjects
		total_subjects <- max_n
		
		# Indicator for each possible scenario of death and tox (e.g. 1 is for death = 1 and tox = 1)
		scenario <- 1 * (death == 1 & tox == 1) + 2 * (death == 1 & tox == 0) + 3 * (death == 0 & tox == 1) + 4 * (death == 0 & tox == 0)	
		###update model ###
		
		jags_works <- 0
		attempt <- 1
		
		while(jags_works == 0) {

			jags_mod <- try(jags.model(	'survEffTox_revision_unif2.bug', 
								data = list('scenario' = scenario[1:total_subjects], 'dose' = dose[1:total_subjects], 'N' = total_subjects,
								'y_s' = surv_time_obs[1:total_subjects], 'h_s' = surv_T, 'y_t' = tox_time_obs[1:total_subjects], 'h_t' = tox_T,  
								'beta0t_m' = beta0t_m, 'beta1t_m' = beta1t_m, 'beta0s_m' = beta0s_m, 'beta1s_m' = beta1s_m, 'beta2s_m' = beta2s_m, 
								'beta0t_p' = beta0t_p, 'beta1t_p' = beta1t_p, 'beta0s_p' = beta0s_p, 'beta1s_p' = beta1s_p, 
								'beta2s_p' = beta2s_p, 'ones' = rep(1, total_subjects)),
								init = list('beta0t' = -2, 'beta1t' = 1, 'beta0s' = -1, 'beta1s' = 1, 'beta2s' = 0),
								n.chains = 1, quiet = TRUE))

			jags_works <- 1*(length(jags_mod) > 1 | attempt == 10)

			attempt <- attempt + 1

		}
		
		update(jags_mod, 4000, quiet = TRUE)

		coda_samples_temp <- coda.samples(jags_mod, c('pt_est', 'ps_est'), 1000, quiet = TRUE)
		##coda_samples_temp <- coda.samples(jags_mod, c('pt_est', 'ps_est', 'beta0t', 'beta1t', 'beta0s', 'beta1s', 'beta2s'), 1000)

		coda_samples <- coda_samples_temp[[1]]

		dist_samples <- matrix(data = 0, ncol = dim(coda_samples)[2]/2, nrow = dim(coda_samples)[1])
		
		for(k in 1:dim(dist_samples)[2]) {

			dist_samples[,k] <- calc_dist(ps = coda_samples[,k], pt = coda_samples[,(k + 4)])

		}
		
		###prob_ps_g_psmin is the posterior probability that the probability of survival is greater than ps_min###
		###prob_pt_l_ptmax is the posterior probability that the probability of toxicity is less than pt_max###

		prob_acceptable <- colMeans(dist_samples > 0)

		###determine acceptable doses###

		acceptable <- 1*(prob_acceptable > gatekeeper)
		accept_dose <- dose_level[acceptable  == 1]

		###Calculate desirability index for acceptable doses###

		ps_cur <- colMeans(coda_samples[,1:4])
		pt_cur <- colMeans(coda_samples[,5:8])

		dist_cur <- calc_dist(ps_cur[acceptable  == 1], pt_cur[acceptable  == 1])

		###identify optimal dose ###

		optimal_dose.perm <- accept_dose[which.max(dist_cur)]
		return(optimal_dose.perm)
}





###################### BETA COMMENSURATE PRIOR ESTIMATOR FUNCTION ###########################

# PURPOSE: The data from the escalation phase are modeled using the joint logistic regression model for efficacy and toxicity that is linear in dose for toxcity and quadratic in dose for efficacy. The prior can be found by centering the prior mean about the probability of [efficacy or toxicity] from the escalation cohort at the optimal dose-level and specifying a commensurate parameter that controls borrowing between the two cohorts.

# INPUT:
## d.esl: Escalation death data from dose escalation trial
## d.exp: Simulated expansion death data
## t.esl: Escalation toxicity data from dose escalation trial
## t.exp: Simulated expansion toxicity data
## optimal_dose: Optimal dose from the dose escalation trial
## dose: Dose assignments from original dose finding model
## tox_time_obs: Simulated toxicity observed time from original dose finding model
## surv_time_obs: Simulated survival observed time from original dose finding model

# OUTPUT: 
## Beta.est: a vector of the beta commensurate prior estimates of efficacy and toxicity

BetaComP <- function(d.esl, d.exp, t.esl, t.exp, optimal_dose, dose, surv_time_obs, tox_time_obs){
	# 
	if(mean(d.exp)==0){
		Beta.est <- rep(NA, 6)
	} else {
		# Beta commensurate prior
		
		scenario <- 1 * (d.esl == 1 & t.esl == 1) + 2 * (d.esl == 1 & t.esl == 0) + 3 * (d.esl == 0 & t.esl == 1) + 4 * (d.esl == 0 & t.esl == 0)
		
		jags_works <- 0
		attempt <- 1
		
		total_subjects <- max_n
		N2 <- dec_n
		while(jags_works == 0) {
		
			comensurate_data <- list('op_dose' =  optimal_dose, 'v.e' = 0.5, 'l.slab.e' = 1, 'u.slab.e' = 4, 'spike.e' = 400, 'v.t' = 0.5, 'l.slab.t' = 1, 'u.slab.t' = 2.25, 'spike.t' = 100, 'N2' = N2, 'surv.exp'=1-d.exp, 'tox.exp'=t.exp, 'scenario' = scenario[1:total_subjects], 'dose' = dose[1:total_subjects], 'N' = total_subjects, 'y_s' = surv_time_obs[1:total_subjects], 'h_s' = surv_T, 'y_t' = tox_time_obs[1:total_subjects], 'h_t' = tox_T, 'beta0t_m' = beta0t_m, 'beta1t_m' = beta1t_m, 'beta0s_m' = beta0s_m, 'beta1s_m' = beta1s_m, 'beta2s_m' = beta2s_m, 'beta0t_p' = beta0t_p, 'beta1t_p' = beta1t_p, 'beta0s_p' = beta0s_p, 'beta1s_p' = beta1s_p,'beta2s_p' = beta2s_p, 'ones' = rep(1, total_subjects))
		
			jags_mod.2 <- try(jags.model('commensurate_model_beta2.bug', data = comensurate_data,init = list('beta0t' = -2, 'beta1t' = 1, 'beta0s' = -1, 'beta1s' = 1, 'beta2s' = 0),
									n.chains = 1, quiet = TRUE))
		
			jags_works <- 1*(length(jags_mod.2) > 1 | attempt == 10)
		
			attempt <- attempt + 1
	
			}
		
		update(jags_mod.2, 4000, quiet = TRUE)
		
		coda_samples_temp.2 <- coda.samples(jags_mod.2, c('exp_eff', 'exp_tox'),1000, quiet = TRUE)
		coda_samples.2 <- coda_samples_temp.2[[1]]
		
		# Beta Commensurate prior data
		ps_est <- mean(coda_samples.2[,1])
		pt_est <- mean(coda_samples.2[,2])
		ps_conf <- quantile(coda_samples.2[,1], c(0.025, 0.975))
		pt_conf <- quantile(coda_samples.2[,2], c(0.025, 0.975))
		Beta.est <- c(ps_est, ps_conf, pt_est, pt_conf)
	}
		
		return(Beta.est)

}





############## GRAND ESTIMATOR FUNCTION #############################

# PURPOSE: To obtain the necessary estimators (expansion only, pooled, 
# Rao-Blackwell, beta commensurate prior, and MEM model estimates). 

# INPUT:
## d.esl: Escalation death data from dose escalation trial
## d.exp: Simulated expansion death data
## t.esl: Escalation toxicity data from dose escalation trial
## t.exp: Simulated expansion toxicity data
## optimal_dose: Optimal dose from the dose escalation trial
## dose: Dose assignments from original dose finding model
## tox_time_obs: Simulated toxicity observed time from original 
##  dose finding model
## surv_time_obs: Simulated survival observed time from original 
##  dose finding model

# OUTPUT: 
## optimal_dose: Optimal dose from the dose escalation trial
## dec: Dose expansion survival and toxicity estimates
## pooled: Dose escalation and dose expansion combined survival and 
##   toxicity estimates
## rb.est: Rao Blackwell survival and toxicity estimates
## bc.est: Beta Commensurate prior survival and toxicity estimates
## esl.est: Dose escalation survival and toxicity estimates
## mb.est: Model-based survival and toxicity estimates

ESTIMATORS <- function(d.esl, d.exp, t.esl, t.exp, optimal_dose, dose, 
  surv_time_obs, tox_time_obs, mb){
	
	if (optimal_dose != 0){
		# EXPANSION ESTIMATES ONLY
		s.dec <- mean(1-d.exp)
		t.dec <- mean(t.exp)
		s.conf <- binconf(sum(1-d.exp), length(1-d.exp), method = "wilson")[2:3]
		t.conf <- binconf(sum(t.exp), length(t.exp), method = "wilson")[2:3]
		dec <- c(s.dec, s.conf, t.dec, t.conf)
		
		# POOLED ESTIMATES
		pooled <- Pool(d.esl, d.exp, t.esl, t.exp)
		
		# RAO-BLACKWELL ESTIMATES
		rb.est <- RaoBlackwell(d.esl, d.exp, t.esl, t.exp, dose, optimal_dose, tox_time_obs, surv_time_obs)
		
		# BETA COMMENSURATE PRIOR ESTIMATES
		bc.est <- BetaComP(d.esl, d.exp, t.esl, t.exp, optimal_dose, dose, surv_time_obs, tox_time_obs)
		
		# ESCALATION ESTIMATES ONLY
		s.escl <- mean(1 - d.esl)
		t.escl <- mean(t.esl)
		s.escl_conf <- binconf(sum(1 - d.esl), length(1 - d.esl), method = "wilson")[2:3]
		t.escl_conf <- binconf(sum(t.esl), length(t.esl), method = "wilson")[2:3]
		esl.est <- c(s.escl, s.escl_conf, t.escl, t.escl_conf)
		
		# MODEL-BASED ESTIMATES (EffTox estimates and respective CI are from SurvEffTox_sim)
		mb.est <- mb
		
		# MEM MODEL ESTIMATES
		#mem.mat <- sapply(seq(0.01, 1, by = 0.01), function(i) MEM.data(1-d.esl,1-d.exp, t.esl, t.exp, dose, optimal_dose,i))
	} else {
		dec <- rep(NA, 6)
		pooled <- rep(NA, 6)
		rb.est <- rep(NA, 2)
		bc.est <- rep(NA, 6)
		esl.est <- rep(NA, 6)
		mb.est <- rep(NA, 6)
		#mem.mat <- matrix(NA, ncol = 100, nrow = 2)
	}
	return(c(optimal_dose, dec, pooled, rb.est, bc.est, esl.est, mb.est))
	
}





########## Estimator Results Function #####

# PURPOSE: To take each simulated scenario and obtain the 
# necessary estimators (expansion only, pooled, Rao-Blackwell, 
# beta commensurate prior, etc.) with different types of 
# inter-trial heterogeneity.

# INPUT:
## sim.results: Simulated results from SurvEffTox
## n_sim: Number of simulations

# OUTPUT: 
## output1: no inter-trial heterogeneity
## output2: lower inter-trial heterogeneity
## output3: upper inter-trial heterogeneity
## output4: random inter-trial heterogeneity

EstiResults <- function(sim_results, n_sim = 1000){
	
	optimal_dose <- sapply(1:n_sim, function(i) sim_results1[[i]]$dat[1])
	d.esl <- sapply(1:n_sim, function(i) sim_results[[i]]$death)
	t.esl <- sapply(1:n_sim, function(i) sim_results[[i]]$tox)
	d.exp1 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$d1) 
	t.exp1 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$t1)
	d.exp2 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$d2) 
	t.exp2 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$t2)
	d.exp3 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$d3) 
	t.exp3 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$t3)
	d.exp4 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$d4) 
	t.exp4 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$t4)
	dose <- sapply(1:n_sim, function(i) sim_results[[i]]$dose)
	surv_time_obs <- sapply(1:n_sim, function(i) sim_results[[i]]$surv_time_obs)
	tox_time_obs <- sapply(1:n_sim, function(i) sim_results[[i]]$tox_time_obs)
	mb <- sapply(1:n_sim, function(i) sim_results[[i]]$mb)

	
	
	output1 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp1[,i], t.esl[,i], t.exp1[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], tox_time_obs[,i], mb[,i]))
	
	output2 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp2[,i], t.esl[,i], t.exp2[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], tox_time_obs[,i], mb[,i]))
	
	# output3 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp3[,i], t.esl[,i], t.exp3[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], tox_time_obs[,i], mb[,i]))
	
	# output4 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp4[,i], t.esl[,i], t.exp4[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], tox_time_obs[,i], mb[,i]))
	
	return(list(output1 = output1, output2 = output2))
}








### Pacakges needed to be loaded
library(Hmisc) # For binconf (i.e. wilson binomial confidence intervals)
library(rjags)
library(parallel)  ####this needs to be parallel###
library(boot)

###set working directory###
#setwd("/Users/mohammadi02/Desktop/EffToxSIM")
#setwd("H:/Simulation Results")

###number of simulations###

n_sim <- 1000

###indicator of whether or not outcomes are treated as time-to-event variables###
###all_data = 0 implies time-to-event, all_data = 1 implies binary outcomes###

all_data <- 1

###maximum sample size/cohorts###

max_n <- 21
cohort_size <- 3
wait_time <- 2
dec_n <- 20


###assign the dose levels###

dose_level <- 1:4

###time in weeks when survival or toxicity is measured###

surv_T <- 24
tox_T <- 4
a <- 1 # Shape 1 parameter for the beta distribution to simulate survival and toxicty time


## Choose case ###
true.prob <- function(case){
	if (case == 1){
			b0t <- -3; b1t <- 1
			b0s <- -0.5; b1s <- 0.7
			
		} else if (case == 2){
			b0t <- -0.5; b1t <- 0.6
			b0s <- 1.2; b1s <- 0.3
		}else if (case == 3){
			b0t <- -3.5; b1t <- 0.9
			b0s <- -2; b1s <- 0.9
		} else if (case == 4) {
			b0t <- -1; b1t <- 0.5
			b0s <- -0.05; b1s <- 1.1
		}
		
	pt <- round(inv.logit(b0t + b1t*0:3), 2)
	ps <- round(inv.logit(b0s + b1s*0:3), 2)
	
  out <- list(pt = pt, ps = ps)
  out
}

###true probabilities of efficacy and toxicity###

###case 1###

#pt <- c(.05, .12, .27, .50) ###beta0 = -3, beta1 = 1###
#ps <- c(.38, .55, .71, .83) ###beta0 = -0.5, beta1 = 0.7###

###case 2###

##pt <- c(.38, .52, .67, .79) ###beta0 = -.5, beta1 = .6###
##ps <- c(.77, .82, .86, .89) ###beta0 = 1.2, beta1 = .3###

###case 3######

##pt <- c(.02, .07, .15, .31) ###beta0 = -3.5, beta1 = .9###
##ps <- c(.12, .25, .45, .67) ###beta0 = -2, beta1 = .9###

### case 4###

##pt <- c(0.27, 0.38,0.50, 0.62) ####beta0 = -1, beta1 = 0.5
##ps <- c(0.49, 0.74, 0.90, 0.96) #### beta0 = -0.05, beta = 1.1




###minimum efficacy and maximum toxicity###

ps_min <- .55
pt_max <- 0.50

ps_star <- .70
pt_star <- .40

###gatekeeper parameters###

gatekeeper <- .05

# ~~~ hyperparameters ~~~ ----

###mean hyperparameters###

beta0t_m <- -3
beta1t_m <- 1

beta0s_m <- -1
beta1s_m <- 1
beta2s_m <- 0

##psi_m <- .5

###sd hyperparameters###

beta0t_s <- 3
beta1t_s <- 2

beta0s_s <- 3
beta1s_s <- 2
beta2s_s <- 0.5


###translating hyperparameters for SD to precision###

beta0t_p <- 1/(beta0t_s^2)
beta1t_p <- 1/(beta1t_s^2)

beta0s_p <- 1/(beta0s_s^2)
beta1s_p <- 1/(beta1s_s^2)
beta2s_p <- 1/(beta2s_s^2)



###function to calculate desirability measure###

###calculate p using contour propossed by Thall and Cook###

find_p <- function(p) {
  (((1 - ps_star) / (1 - ps_min)) ^ p + (pt_star / pt_max) ^ p - 1) ^ 2 
}

p <- optimize(find_p, c(-100,100))$minimum


calc_dist <- function(ps, pt) {

  dist <- 1 - (((1 - ps) / (1 - ps_min)) ^ p + (pt / pt_max) ^ p ) ^ (1 / p)

  dist

}







# Calculate the new estimators (i.e. pooled, Rao Blackwell, Beta Commensurate Prior, etc.) from the 
# EffTox Simulated data from the four cases 

# ~~~ for debugging ~~~ ----
# debug(EstiResults)
# debug(Pool)
# debug(DOSE.MODEL)
# debug(BetaComP)
# debug(ESTIMATORS)

load('RData/sim_results1.RData')
output1 <- EstiResults(sim_results1)
save(output1, file = "output1.RData")

load('RData/sim_results2.RData')
output2 <- EstiResults(sim_results2)
save(output2, file = "output2.RData")

load('RData/sim_results3.RData')
output3 <- EstiResults(sim_results3)
save(output3, file = "output3.RData")

load('RData/sim_results4.RData')
output4 <- EstiResults(sim_results4)
save(output4, file = "output4.RData")





