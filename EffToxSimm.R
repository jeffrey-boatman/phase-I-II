

### Purpose: Simulate phase I-II trial and follow-up cohort to obtain survival and toxicity outcomes based on dose-finding model proposed by Koopmeiners et al. 


# Load libraries
library(rjags)
library(parallel)  ####this needs to be parallel###
library(boot)

###set working directory###
#setwd("/Users/mohammadi02/Desktop/EffToxSIM")

###number of simulations###

n_sim <- 1000

###indicator of whether or not outcomes are treated as time-to-event variables###
###all_data = 0 implies time-to-event, all_data = 1 implies binary outcomes###

all_data <- 0

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




###true probabilities of efficacy and toxicity###

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
	
return(list(pt = pt, ps = ps))
}


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

###hyperparameters###

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

return( ( ((1 - ps_star)/(1 - ps_min))^p + (pt_star/pt_max)^p - 1)^2 )

}

p <- optimize(find_p, c(-100,100))$minimum


calc_dist <- function(ps, pt) {

dist <- 1 - ( ((1 - ps)/(1 - ps_min))^p + (pt/pt_max)^p )^(1/p)

return(dist)

}


############## DOSE EXPANSION COHORT #############

# PURPOSE: To simulate binary deaths and toxicity events based on the optimal dose found in the dose-finding algorithm and the true efficacy and toxicity values. Specifically, we simulate a follow-up cohort assuming no inter-trial heterogeneity, lower inter-trial heterogeneity, upper inter-trial heterogeneity, and inter-trial heterogeniety. (We add N(0,0.25) to the logit of the probability of toxcity and efficacy)

## INPUT:
# optimal_dose: The optimal dose obtained by the dose finding model
# case: case number to get respective efficacy and toxicity coefficents from the logistic model

## OUTPUT (list):
# d1: simulated DEC death outcomes for 20 subjects assuming no inter-trial heterogeneity
# t1: simulated DEC tox outcomes for 20 subjects assuming no inter-trial heterogeneity
# d2: simulated DEC death outcomes for 20 subjects assuming lower inter-trial heterogeneity
# t2: simulated DEC tox outcomes for 20 subjects assuming lower inter-trial heterogeneity
# ps.lwr: True efficacy value at the optimal dose assuming lower inter-trial heterogeneity
# pt.lwr: True efficacy value at the optimal dose assuming lower inter-trial heterogeneity
# d3: simulated DEC death outcomes for 20 subjects assuming upper inter-trial heterogeneity
# t3: simulated DEC tox outcomes for 20 subjects assuming upperinter-trial heterogeneity
# ps.upr: True efficacy value at the optimal dose assuming upper inter-trial heterogeneity
# pt.upr: True efficacy value at the optimal dose assuming upper inter-trial heterogeneity
# d4: simulated DEC death outcomes for 20 subjects assuming random inter-trial heterogeneity
# t4: simulated DEC tox outcomes for 20 subjects assuming random inter-trial heterogeneity
# ps.lwr: True efficacy value at the optimal dose assuming random inter-trial heterogeneity
# pt.lwr: True efficacy value at the optimal dose assuming random inter-trial heterogeneity

DEC <- function(optimal_dose, case){
	
	# It optimal dose is zero (i.e. trial stopped), then return an empty list
	if(optimal_dose != 0){
		# IF else statment to call beta coefficents for the logistic model
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
		
		# The true values efficacy and toxicity (1. No inter-trial heterogeneity)
		pt <- inv.logit(b0t + b1t*0:3)[optimal_dose]
		ps <- inv.logit(b0s + b1s*0:3)[optimal_dose]
		
		# Jitter the true values efficacy and toxicity (2. lower inter-trial heterogeneity)
		pt.lwr <- inv.logit(b0t + b1t*0:3 - abs(rnorm(n = 1, sd = 0.25)))[optimal_dose]
		ps.lwr <- inv.logit(b0s + b1s*0:3 - abs(rnorm(n = 1, sd = 0.25)))[optimal_dose]
		
		# Jitter the true values efficacy and toxicity (3. upper inter-trial heterogeneity)
		pt.upr <- inv.logit(b0t + b1t*0:3 + abs(rnorm(n = 1, sd = 0.25)))[optimal_dose]
		ps.upr <- inv.logit(b0s + b1s*0:3 + abs(rnorm(n = 1, sd = 0.25)))[optimal_dose]
		
		# Jitter the true values efficacy and toxicity (4. random inter-trial heterogeneity)
		pt.ith <- inv.logit(b0t + b1t*0:3 + rnorm(n = 1, sd = 0.25))[optimal_dose]
		ps.ith <- inv.logit(b0s + b1s*0:3 + rnorm(n = 1, sd = 0.25))[optimal_dose]
		
		# Sampling 20 subjects toxicity based on true values of the optimal_dose
		t1 <- rbinom(n = dec_n, size = 1, prob = pt)
		t2 <- rbinom(n = dec_n, size = 1, prob = pt.lwr)
		t3 <- rbinom(n = dec_n, size = 1, prob = pt.upr)
		t4 <- rbinom(n = dec_n, size = 1, prob = pt.ith)
		
		# Sampling 20 subjects death status based on true values of the optimal dose
		d1 <- 1 - rbinom(n = dec_n, size = 1, prob = ps)
		d2 <- 1 - rbinom(n = dec_n, size = 1, prob = ps.lwr)
		d3 <- 1 - rbinom(n = dec_n, size = 1, prob = ps.upr)
		d4 <- 1 - rbinom(n = dec_n, size = 1, prob = ps.ith)
		} else {
			d1 = t1 = d2 = t2  = d3 = t3 = d4 = t4 = rep(NA, dec_n)
			ps.lwr = pt.lwr = ps.upr = pt.upr = ps.ith = pt.ith = NA
		}

	
	return(list(d1 = d1, t1 = t1, d2 = d2, t2 = t2, ps.lwr = ps.lwr, pt.lwr = pt.lwr, d3 = d3, t3 = t3, ps.upr = ps.upr, pt.upr = pt.upr, d4 = d4, t4 = t4, ps.ith = ps.ith, pt.ith = pt.ith))
}


##########################################################


############## Dose Escalation Clinical Trial #############

SurvEffTox_sim <- function(xxx) {


	###starting values for week, dose and max dose###

	cur_week <- 1
	cur_dose <- 1
	max_dose <- 1

	###vectors for storing data###

	dose <- rep(NA, max_n)  					###dose###
	surv <- rep(NA, max_n)  					###survival at current time###
	tox <- rep(NA, max_n)   					###toxicity at current time###
	surv_time <- rep(NA, max_n)					###survival time###
	tox_time <- rep(NA, max_n)					###toxicity time###
	start_week <- cumsum(rexp(n = max_n, rate = 1/wait_time))			###week subject starts on trial - currently 1, 3, 5, etc. - this will change###
	start_week <- start_week - start_week[1]
	surv_mult <- rep(NA, max_n)					###weight for survival outcome###
	tox_mult <- rep(1, max_n)					###weight for toxicity outcome###
	outcome <- matrix(ncol = 4, nrow = max_n, data = NA)

	X <- matrix(data = 0, nrow = max_n, ncol = 4)		###multinomial outcome###

	stop_trial <- 0
	total_subjects <- cohort_size

	while(stop_trial == 0){						

		###assign current dose to current cohort### 

		dose[(total_subjects - cohort_size + 1):total_subjects] <- cur_dose

		###simulate survival and toxicity time from uniform distribution - I'd like to try other distributions###

		surv_time[(total_subjects - cohort_size + 1): total_subjects] <- surv_T * pmax(rbinom(n = cohort_size, size = 1, prob = ps[cur_dose]), rbeta(n = cohort_size, shape1 = a, shape2 = 1))
		tox_time[(total_subjects - cohort_size + 1): total_subjects] <- tox_T * pmax((1 - rbinom(n = cohort_size, size = 1, prob = pt[cur_dose])), rbeta(n = cohort_size, shape1 = a, shape2 = 1))

		###determined how many weeks each outcome has been observed - if all_data 

		if(total_subjects < max_n & all_data == 0) {

			###determine week when next cohort will enter - this will change with stochastic enrollment time###

			cur_week <- start_week[(total_subjects + 1)]

			surv_u <- pmin(cur_week - start_week, surv_T)
			tox_u <- pmin(cur_week - start_week, tox_T)

		}else{

			surv_u <- rep(surv_T, max_n)
			tox_u <- rep(tox_T, max_n)

		}

		death <- 1*(surv_time < surv_u)
		tox <- 1*(tox_time < tox_u)

		surv_time_obs <- pmin(surv_time, surv_u)
		tox_time_obs <- pmin(tox_time, tox_u)

		##outcome[,1] <- 1 * (death == 1 & tox == 1)
		##outcome[,2] <- 1 * (death == 1 & tox == 0)
		##outcome[,3] <- 1 * (death == 0 & tox == 1)
		##outcome[,4] <- 1 * (death == 0 & tox == 0)

		scenario <- 1 * (death == 1 & tox == 1) + 2 * (death == 1 & tox == 0) + 3 * (death == 0 & tox == 1) + 4 * (death == 0 & tox == 0)

		###update model###

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

		###stop trial if all doses are unacceptable###
		###otherwise, calculate desirability index for acceptable doses###

		if(sum(acceptable) == 0) {

			stop_trial <- 1
			cur_dose <- 0

		} else {

			ps_cur <- colMeans(coda_samples[,1:4])
			pt_cur <- colMeans(coda_samples[,5:8])

			dist_cur <- calc_dist(ps_cur[acceptable  == 1], pt_cur[acceptable  == 1])

			###identify best dose and update current dose###

			best_dose <- accept_dose[which.max(dist_cur)]

			cur_dose <- min(best_dose, max_dose + 1)
			max_dose <- max(max_dose, cur_dose)

			stop_trial <- total_subjects >= max_n

		}

		###update max dose###

		total_subjects <- total_subjects + cohort_size

	}

	## If trial stopped (doses unacceptable), then set optimal dose to zero. 
	## Otherwise, set the best dose as the new optimal dose
	if(sum(acceptable) == 0){

		optimal_dose <- 0

	} else {

		optimal_dose <- best_dose

	}
		
	## indicator of whether or not outcomes are treated as time-to-event variables
	## all_data = 0 implies time-to-event, all_data = 1 implies binary outcomes
	if(all_data == 0) {

		study_duration <- max(pmax(pmin(surv_time[1:(total_subjects - cohort_size)], surv_T), pmin(tox_time[1:(total_subjects - cohort_size)], tox_T)) + start_week[1:(total_subjects - cohort_size)])

	}

	if(all_data == 1) {

		for(i in 1:max(((total_subjects/cohort_size) - 2),1) ){

		earliest_start <- max(pmax(pmin(surv_time[((i - 1) * cohort_size + 1):(i * cohort_size)], surv_T), pmin(tox_time[((i - 1) * cohort_size + 1):(i * cohort_size)], tox_T)) + start_week[((i - 1) * cohort_size + 1):(i * cohort_size)])

		start_week[(i * cohort_size + 1):(i * cohort_size + cohort_size)] <- pmax(start_week[(i * cohort_size + 1):(i * cohort_size + cohort_size)], earliest_start)

		}

	study_duration <- max(pmax(pmin(surv_time[1:(total_subjects - cohort_size)], surv_T), pmin(tox_time[1:(total_subjects - cohort_size)], tox_T)) + start_week[1:(total_subjects - cohort_size)])


	}

	## Simulate the follow-up cohort based on the optimal dose found by the dose-finding model
	## using the function DEC directly in the returned list. Four situations will be outputed (See DEC function)
	
	## Create list of data from dose finding model
	dat <- c(optimal_dose, sum(dose == 1, na.rm = T), sum(dose == 2, na.rm = T), sum(dose == 3, na.rm = T), sum(dose == 4, na.rm = T), study_duration)
	
	# model based estimates
	s.mb <- ps_cur[optimal_dose]
	t.mb <- pt_cur[optimal_dose]
	
	s.mb_conf <- quantile(coda_samples[,optimal_dose], prob = c(0.025, 0.975))
	t.mb_conf <- quantile(coda_samples[,(optimal_dose +4)], prob = c(0.025, 0.975))
	if(optimal_dose == 0){s.mb = t.mb = NA;  s.mb_conf = t.mb_conf = rep(NA, 2)}
	mb <- c(s.mb, s.mb_conf, t.mb, t.mb_conf)
	
return(list(dat = dat, death = death, tox = tox, dec = DEC(optimal_dose, case), dose = dose, surv_time_obs = surv_time_obs, tox_time_obs = tox_time_obs, mb = mb))

}


## START SIMULATION 

# SET SEEED FOR REPLICATION
RNGkind("L'Ecuyer-CMRG")
set.seed(7777442)

# Set number of cores for parallel
# **Set to one core but can be changed**
no_cores <- detectCores()-7

# start.time <- Sys.time()
# end.time <- Sys.time()
# end.time-start.time

# Simulation for Case 1
ps <- true.prob(case = 1)$ps
pt <- true.prob(case = 1)$pt
sim_results1 <- mclapply(1:n_sim, SurvEffTox_sim, mc.cores = no_cores, mc.set.seed = FALSE)

# Simulation for Case 2
ps <- true.prob(case = 2)$ps
pt <- true.prob(case = 2)$pt
sim_results2 <- mclapply(1:n_sim, SurvEffTox_sim, mc.cores = no_cores, mc.set.seed = FALSE)

# Simulation for Case 3
ps <- true.prob(case = 3)$ps
pt <- true.prob(case = 3)$pt
sim_results3 <- mclapply(1:n_sim, SurvEffTox_sim, mc.cores = no_cores, mc.set.seed = FALSE)

# Simulation for Case 4
ps <- true.prob(case = 4)$ps
pt <- true.prob(case = 4)$pt
sim_results4 <- mclapply(1:n_sim, SurvEffTox_sim, mc.cores = no_cores, mc.set.seed = FALSE)


# Save RData from all cases to Home Drive
save(sim_results1, file = "/Users/mohammadi02/Desktop/EffToxSIM/sim_results1.RData")
save(sim_results2, file = "/Users/mohammadi02/Desktop/EffToxSIM/sim_results2.RData")
save(sim_results3, file = "/Users/mohammadi02/Desktop/EffToxSIM/sim_results3.RData")
save(sim_results4, file = "/Users/mohammadi02/Desktop/EffToxSIM/sim_results4.RData")





