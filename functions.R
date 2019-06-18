


############## Simulate Dose Escalation Clinical Trial #############

SurvEffTox_sim <- function(ps, pt, keep_dose0 = FALSE) {
  all_data     <- 1
	optimal_dose <- 0
	
	while(optimal_dose == 0) {
	  
	 ###starting values for week, dose and max dose###

	  cur_week <- 1
  	cur_dose <- 1
	  max_dose <- 1

  	###vectors for storing data###

  	dose      <- rep(NA, max_n)  ###dose###
   	surv      <- rep(NA, max_n)  ###survival at current time###
	  tox       <- rep(NA, max_n)  ###toxicity at current time###
	  surv_time <- rep(NA, max_n)  ###survival time###
	  tox_time  <- rep(NA, max_n)  ###toxicity time###
	
	  ###week subject starts on trial - currently 1, 3, 5, etc. - this will change###
	  start_week <- cumsum(rexp(n = max_n, rate = 1/wait_time))
	  start_week <- start_week - start_week[1]

	  surv_mult <- rep(NA, max_n)  ###weight for survival outcome###
  	tox_mult  <- rep(1, max_n)   ###weight for toxicity outcome###
	  outcome   <- matrix(ncol = 4, nrow = max_n, data = NA)

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

		  } else {

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

		  like_scenario <- 1 * (death == 1 & tox == 1) + 2 * (death == 1 & tox == 0) + 3 * (death == 0 & tox == 1) + 4 * (death == 0 & tox == 0)

		  ###update model###

		  jags_works <- 0
		  attempt <- 1

		  while(jags_works == 0) {

			  jags_mod <- try(jags.model('bugs/survEffTox_revision_unif2.bug', 
			    data = list(
              #'scenario' = like_scenario[1:total_subjects], 
              'dose'     = dose[1:total_subjects], 
              'N'        = total_subjects,
			        tox = tox[1:total_subjects],
			        eff = 1 - death[1:total_subjects],
              #'y_s'      = surv_time_obs[1:total_subjects], 
              #'h_s'      = surv_T, 
              #'y_t'      = tox_time_obs[1:total_subjects], 
              #'h_t'      = tox_T,  
              'beta0t_m' = beta0t_m, 
              # 'beta1t_m' = beta1t_m, 
              'beta0s_m' = beta0s_m, 
              # 'beta1s_m' = beta1s_m, 
              'beta2s_m' = beta2s_m, 
              'beta0t_p' = beta0t_p, 
              #'beta1t_p' = beta1t_p, 
              'beta0s_p' = beta0s_p, 
              #'beta1s_p' = beta1s_p, 
              'beta2s_p' = beta2s_p, 
              #'ones'     = rep(1, total_subjects),
			        'beta1t_r'      = beta1_shapeT,
			        'beta1t_lambda' = beta1_rateT,
			        'beta1s_r'      = beta1_shapeE,
			        'beta1s_lambda' = beta1_rateE),
            init = list(
              'beta0t' = -2, 
              'beta1t' = 1, 
              'beta0s' = -1, 
              'beta1s' = 1, 
              'beta2s' = 0),
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
	  if(keep_dose0) break
	}

	## Create list of data from dose finding model
	dat <- c(optimal_dose, sum(dose == 1, na.rm = T), sum(dose == 2, na.rm = T), sum(dose == 3, na.rm = T), sum(dose == 4, na.rm = T), study_duration)
	
	# model based estimates
	if(optimal_dose == 0){
	  s.mb = t.mb = NA;  s.mb_conf = t.mb_conf = rep(NA, 2)
	} else {
	  s.mb <- ps_cur[optimal_dose]
	  t.mb <- pt_cur[optimal_dose]
	
	  s.mb_conf <- quantile(coda_samples[,optimal_dose], prob = c(0.025, 0.975))
	  t.mb_conf <- quantile(coda_samples[,(optimal_dose +4)], prob = c(0.025, 0.975))
	}

	mb <- c(s.mb, s.mb_conf, t.mb, t.mb_conf)
	
	if(optimal_dose == 0) {
	  p_tox_optimal <- NA
	  p_eff_optimal <- NA
	} else {
	  p_tox_optimal <- pt[optimal_dose]
	  p_eff_optimal <- ps[optimal_dose]
	}
	
  out <- list(
    # seed          = seed,
    # scenario      = scenario,
    # dat           = dat,
    p_tox_optimal  = p_tox_optimal,
    p_eff_optimal  = p_eff_optimal,
    death          = death, 
    tox            = tox, 
    # dec           = DEC(optimal_dose, scenario), 
    dose           = dose, 
    surv_time_obs  = surv_time_obs, 
    tox_time_obs   = tox_time_obs, 
    mb             = mb,
    optimal_dose   = optimal_dose,
    study_duration = study_duration)
  out

}

############## DOSE EXPANSION COHORT #############

# PURPOSE: To simulate binary deaths and toxicity events based on the optimal dose found in 
# the dose-finding algorithm and the true efficacy and toxicity values. Specifically, 
# we simulate a follow-up cohort assuming no inter-trial heterogeneity, lower inter-trial 
# heterogeneity, upper inter-trial heterogeneity, and inter-trial heterogeniety. (We add N(0,0.25) 
# to the logit of the probability of toxcity and efficacy)

## INPUT:
# optimal_dose: The optimal dose obtained by the dose finding model
# scenario: scenario number to get respective efficacy and toxicity coefficents from the logistic model

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

# DEC <- function(optimal_dose, scenario){
# 	
# 	# It optimal dose is zero (i.e. trial stopped), then return an empty list
# 	if(optimal_dose != 0){
# 		# IF else statment to call beta coefficents for the logistic model
# 		if (scenario == 1){
# 			b0t <- -3; b1t <- 1
# 			b0s <- -0.5; b1s <- 0.7
# 			
# 		} else if (scenario == 2){
# 			b0t <- -0.5; b1t <- 0.6
# 			b0s <- 1.2; b1s <- 0.3
# 		}else if (scenario == 3){
# 			b0t <- -3.5; b1t <- 0.9
# 			b0s <- -2; b1s <- 0.9
# 		} else if (scenario == 4) {
# 			b0t <- -1; b1t <- 0.5
# 			b0s <- -0.05; b1s <- 1.1
# 		}
# 		
# 		# The true values efficacy and toxicity (1. No inter-trial heterogeneity)
# 		pt <- inv.logit(b0t + b1t*0:3)[optimal_dose]
# 		ps <- inv.logit(b0s + b1s*0:3)[optimal_dose]
# 		
# 		# Jitter the true values efficacy and toxicity (2. lower inter-trial heterogeneity)
# 		pt.lwr <- inv.logit(b0t + b1t*0:3 - abs(rnorm(n = 1, sd = 0.25)))[optimal_dose]
# 		ps.lwr <- inv.logit(b0s + b1s*0:3 - abs(rnorm(n = 1, sd = 0.25)))[optimal_dose]
# 		
# 		# Jitter the true values efficacy and toxicity (3. upper inter-trial heterogeneity)
# 		pt.upr <- inv.logit(b0t + b1t*0:3 + abs(rnorm(n = 1, sd = 0.25)))[optimal_dose]
# 		ps.upr <- inv.logit(b0s + b1s*0:3 + abs(rnorm(n = 1, sd = 0.25)))[optimal_dose]
# 		
# 		# Jitter the true values efficacy and toxicity (4. random inter-trial heterogeneity)
# 		pt.ith <- inv.logit(b0t + b1t*0:3 + rnorm(n = 1, sd = 0.25))[optimal_dose]
# 		ps.ith <- inv.logit(b0s + b1s*0:3 + rnorm(n = 1, sd = 0.25))[optimal_dose]
# 		
# 		# Sampling 20 subjects toxicity based on true values of the optimal_dose
# 		t1 <- rbinom(n = dec_n, size = 1, prob = pt)
# 		t2 <- rbinom(n = dec_n, size = 1, prob = pt.lwr)
# 		t3 <- rbinom(n = dec_n, size = 1, prob = pt.upr)
# 		t4 <- rbinom(n = dec_n, size = 1, prob = pt.ith)
# 		
# 		tox_time_obs1 <- tox_T * (t1 * rbeta(dec_n, 1, 1) + (1 - t1))
# 		tox_time_obs2 <- tox_T * (t2 * rbeta(dec_n, 1, 1) + (1 - t2))
# 		tox_time_obs3 <- tox_T * (t3 * rbeta(dec_n, 1, 1) + (1 - t3))
# 		tox_time_obs4 <- tox_T * (t4 * rbeta(dec_n, 1, 1) + (1 - t4))
# 		
# 		# Sampling 20 subjects death status based on true values of the optimal dose
# 		d1 <- 1 - rbinom(n = dec_n, size = 1, prob = ps)
# 		d2 <- 1 - rbinom(n = dec_n, size = 1, prob = ps.lwr)
# 		d3 <- 1 - rbinom(n = dec_n, size = 1, prob = ps.upr)
# 		d4 <- 1 - rbinom(n = dec_n, size = 1, prob = ps.ith)
# 		
# 		surv_time_obs1 <- surv_T * (d1 * rbeta(dec_n, 1, 1) + (1 - d1))
# 		surv_time_obs2 <- surv_T * (d2 * rbeta(dec_n, 1, 1) + (1 - d2))
# 		surv_time_obs3 <- surv_T * (d3 * rbeta(dec_n, 1, 1) + (1 - d3))
# 		surv_time_obs4 <- surv_T * (d4 * rbeta(dec_n, 1, 1) + (1 - d4))
# 		
# 		} else {
# 		  pt <- NA
# 		  ps <- NA
# 			d1 = t1 = d2 = t2  = d3 = t3 = d4 = t4 = rep(NA, dec_n)
# 			tox_time_obs1 = tox_time_obs2 = tox_time_obs3 = tox_time_obs4 = rep(NA, dec_n)
# 			surv_time_obs1 = surv_time_obs2 = surv_time_obs3 = surv_time_obs4 = rep(NA, dec_n)
# 			ps.lwr = pt.lwr = ps.upr = pt.upr = ps.ith = pt.ith = NA
# 		}
# 
# 	
# 	return(list(
# 	  pt     = pt,
# 	  ps     = ps,
# 	  d1     = d1, 
# 	  t1     = t1, 
# 	  d2     = d2, 
# 	  t2     = t2, 
# 	  ps.lwr = ps.lwr, 
# 	  pt.lwr = pt.lwr, 
# 	  d3     = d3, 
# 	  t3     = t3, 
# 	  ps.upr = ps.upr, 
# 	  pt.upr = pt.upr, 
# 	  d4     = d4, 
# 	  t4     = t4, 
# 	  ps.ith = ps.ith, 
# 	  pt.ith = pt.ith,
# 	  tox_time_obs1 = tox_time_obs1,
# 	  tox_time_obs2 = tox_time_obs2,
# 	  tox_time_obs3 = tox_time_obs3,
# 	  tox_time_obs4 = tox_time_obs4,
# 	  surv_time_obs1 = surv_time_obs1,
# 	  surv_time_obs2 = surv_time_obs2,
# 	  surv_time_obs3 = surv_time_obs3,
# 	  surv_time_obs4 = surv_time_obs4))
# }

DEC <- function(optimal_dose, ps, pt, 
  data_gen_type = c("no_effect", "hiTox_loEff", "ITH")){
  
  data_gen_type <- match.arg(data_gen_type)
  
  prob_eff <- ps
  prob_tox <- pt

  if(data_gen_type == "hiTox_loEff") {
    prob_eff <- plogis(qlogis(ps) - abs(rnorm(4, sd = 1 / 4)))
    prob_tox <- plogis(qlogis(pt) + abs(rnorm(4, sd = 1 / 4)))
  } else if(data_gen_type == "ITH") {
    prob_eff <- plogis(qlogis(ps) + rnorm(4, sd = 1 / 4))
    prob_tox <- plogis(qlogis(pt) + rnorm(4, sd = 1 / 4))
  }
  
  prob_eff <- prob_eff[optimal_dose]
  prob_tox <- prob_tox[optimal_dose]
  
  n_dec <- 20
  
  death <- rbinom(n_dec, 1, 1 - prob_eff)
  tox   <- rbinom(n_dec, 1, prob_tox)
  
  tox_time_obs  <- tox_T  * (tox    * rbeta(dec_n, 1, 1) + (1 - tox))
	surv_time_obs <- surv_T * (death  * rbeta(dec_n, 1, 1) + (1 - death))

  out <- list(
    prob_eff      = prob_eff,
    prob_tox      = prob_tox,
    death         = death,
    tox           = tox,
    tox_time_obs  = tox_time_obs,
    surv_time_obs = surv_time_obs,
    data_gen_type = data_gen_type
  )
  
  out
	
}


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


Pool <- function(d.esl, d.exp, t.esl, t.exp, dose, optimal_dose){
	d.pool <- c(d.esl[dose == optimal_dose], d.exp)
	t.pool <- c(t.esl[dose == optimal_dose], t.exp)
	s.conf <- binconf(sum(1-d.pool),length(1-d.pool), method = "wilson")[2:3]
	t.conf <- binconf(sum(t.pool),length(t.pool), method = "wilson")[2:3]
	pooled <- c(mean(1-d.pool), s.conf, mean(t.pool), t.conf)
	names(pooled) <- c(
	  paste0("pool_eff_", c("est", "lwr", "upr")),
	  paste0("pool_tox_", c("est", "lwr", "upr")))
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

# for debugging UMVUE estimator
# set.seed(1)
# d.esl <- rbinom(21, 1, 0.2)
# d.exp <- rbinom(20, 1, 0.2)
# t.esl <- rbinom(21, 1, 0.35)
# t.exp <- rbinom(20, 1, 0.35)
# optimal_dose <- 1
# surv_time_obs <- d.esl * rbeta(21, 1, 1) * 24 + (1 - d.esl) * 24
# tox_time_obs  <- t.esl * rbeta(21, 1, 1) *  4 + (1 - t.esl) *  4
# surv_time_obs.exp <- d.exp * rbeta(20, 1, 1) * 24 + (1 - d.exp) * 24
# tox_time_obs.exp  <- t.exp * rbeta(20, 1, 1) *  4 + (1 - t.exp) *  4
# dose <- rep(1, 21)
# dose[c(4, 5, 6)] <- 2
# debug(RaoBlackwell)
# RaoBlackwell(d.esl, d.exp, t.esl, t.exp, dose, optimal_dose,
#   tox_time_obs, surv_time_obs, tox_time_obs.exp, surv_time_obs.exp)

RaoBlackwell <- function(d.esl, d.exp, t.esl, t.exp, dose, optimal_dose, 
  tox_time_obs, surv_time_obs, tox_time_obs.exp, surv_time_obs.exp, include.ci){
  
  # set.seed(seed)
  
  n_esl <- length(d.esl)
  n_exp <- length(d.exp)
  
  esl_seq <- seq_len(n_esl)
  
  all_death         <- c(d.esl, d.exp)
  all_tox           <- c(t.esl, t.exp)
  all_dose          <- c(dose, rep(optimal_dose, n_exp))
  all_tox_time_obs  <- c(tox_time_obs, tox_time_obs.exp)
  all_surv_time_obs <- c(surv_time_obs, surv_time_obs.exp)
  
  which_optimal_dose <- which(all_dose == optimal_dose)
  
  n_perm        <- ifelse(include.ci, 1e3, 10)
  dose_matches  <- 0
  num_tries     <- 0
  max_num_tries <- ifelse(include.ci, 1e4, 100)
  s.rb <- t.rb  <- numeric(n_perm)
  
  while(dose_matches < n_perm & num_tries <= max_num_tries) {
    death_perm         <- all_death
    tox_perm           <- all_tox
    dose_perm          <- all_dose
    tox_time_obs_perm  <- all_tox_time_obs
    surv_time_obs_perm <- all_surv_time_obs

    perm_index <- sample(which_optimal_dose)
  
    death_perm[which_optimal_dose]         <- all_death[perm_index]
    tox_perm[which_optimal_dose]           <- all_tox[perm_index]
    dose_perm[which_optimal_dose]          <- all_dose[perm_index]
    tox_time_obs_perm[which_optimal_dose]  <- all_tox_time_obs[perm_index]
    surv_time_obs_perm[which_optimal_dose] <- all_surv_time_obs[perm_index]
  
    # perm_optimal_dose <- DOSE.MODEL(
    #   death         = death_perm[esl_seq],
    #   tox           = tox_perm[esl_seq],
    #   dose          = dose_perm[esl_seq],
    #   tox_time_obs  = tox_time_obs_perm[esl_seq],
    #   surv_time_obs = surv_time_obs_perm[esl_seq]
    # )
    # ~~~~ ----
		jags_works <- 0
		attempt <- 1

	  jags_mod <- try(jags.model('bugs/survEffTox_revision_unif2.bug', 
      data = list(
            'dose'     = dose_perm[1:n_esl], 
            'N'        = n_esl,
		        'tox'      = tox_perm[1:n_esl],
		        'eff'      = 1 - death_perm[1:n_esl],
            'beta0t_m' = beta0t_m, 
            'beta0s_m' = beta0s_m, 
            'beta2s_m' = beta2s_m,
            'beta0t_p' = beta0t_p, 
            'beta0s_p' = beta0s_p, 
            'beta2s_p' = beta2s_p, 
		        'beta1t_r'      = beta1_shapeT,
		        'beta1t_lambda' = beta1_rateT,
		        'beta1s_r'      = beta1_shapeE,
		        'beta1s_lambda' = beta1_rateE),
          init = list(
            'beta0t' = -2, 
            'beta1t' = 1, 
            'beta0s' = -1, 
            'beta1s' = 1, 
            'beta2s' = 0),
          n.chains = 1, quiet = TRUE))
	  update(jags_mod, 4000, quiet = TRUE)

	  coda_samples_temp <- coda.samples(jags_mod, c('pt_est', 'ps_est'), 1000, quiet = TRUE)

	  coda_samples <- coda_samples_temp[[1]]

	  dist_samples <- matrix(data = 0, ncol = dim(coda_samples)[2]/2, nrow = dim(coda_samples)[1])

	  for(k in 1:dim(dist_samples)[2]) {

		  dist_samples[,k] <- calc_dist(ps = coda_samples[,k], pt = coda_samples[,(k + 4)])

	  }

	  prob_acceptable <- colMeans(dist_samples > 0)

	  ###determine acceptable doses###

	  acceptable <- 1*(prob_acceptable > gatekeeper)
	  accept_dose <- dose_level[acceptable  == 1]

	  if(sum(acceptable) == 0) {
  
	    perm_optimal_dose <- 0

	  } else {

		  ps_cur <- colMeans(coda_samples[,1:4])
		  pt_cur <- colMeans(coda_samples[,5:8])

		  dist_cur <- calc_dist(ps_cur[acceptable  == 1], pt_cur[acceptable  == 1])
		  ###identify best dose and update current dose###

		  perm_optimal_dose <- accept_dose[which.max(dist_cur)]

	  }
  
    if(perm_optimal_dose == optimal_dose) {
      dose_matches <- dose_matches + 1
      s.rb[dose_matches] <- mean(rev(1 - death_perm)[1:n_exp])
      t.rb[dose_matches] <- mean(rev(tox_perm)[1:n_exp])
    }
    num_tries <- num_tries + 1
    
  }

  rb_eff <- s.rb[seq_len(dose_matches)]
  rb_tox <- t.rb[seq_len(dose_matches)]
  
  mean_and_ci <- function(x) {
    m <- mean(x, na.rm = TRUE)
    s <- sd(x, na.rm = TRUE)
    out <- c(m, m + c(-1, 1) * qnorm(0.975) * s)
    names(out) <- c("mean", "lwr", "upr")
    out
  }
  
  if(include.ci) {
  	# RB.est <- c(
  	#   mean(rb_eff, na.rm=TRUE), quantile(rb_eff, c(0.025, 0.975), na.rm = TRUE),
  	#   mean(rb_tox, na.rm=TRUE), quantile(rb_tox, c(0.025, 0.975), na.rm = TRUE))
  	RB.est <- c(
  	  mean_and_ci(rb_eff),
  	  mean_and_ci(rb_tox))
  	names(RB.est) <- c(paste0("rb_eff_", c("est", "lwr", "upr")), 
  	  paste0("rb_tox_", c("est", "lwr", "upr")))
  } else{
  	RB.est <- c(
  	  mean(rb_eff, na.rm=TRUE), 
  	  mean(rb_tox, na.rm=TRUE))
  	names(RB.est) <- c("rb_eff_est", "rb_tox_est")
  }
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

#' DOSE.MODEL <- function(death, tox, dose, tox_time_obs, surv_time_obs) {
#' 		# Set total subjects to 21 subjects
#' 		total_subjects <- max_n
#' 		
#' 		# Indicator for each possible scenario of death and tox (e.g. 1 is for death = 1 and tox = 1)
#' 		like_scenario <- 1 * (death == 1 & tox == 1) + 2 * (death == 1 & tox == 0) + 3 * (death == 0 & tox == 1) + 4 * (death == 0 & tox == 0)	
#' 		###update model ###
#' 		
#' 		jags_works <- 0
#' 		attempt <- 1
#' 		
#' 		while(jags_works == 0) {
#' 
#' 			jags_mod <- try(jags.model('bugs/survEffTox_revision_unif2.bug', 
#'         data = list(
#'           'scenario' = like_scenario[1:total_subjects], 
#'           'dose'     = dose[1:total_subjects], 
#'           'N'        = total_subjects,
#'           'y_s'      = surv_time_obs[1:total_subjects], 
#'           'h_s'      = surv_T, 
#'           'y_t'      = tox_time_obs[1:total_subjects], 
#'           'h_t'      = tox_T,  
#'           'beta0t_m' = beta0t_m, 
#'           # 'beta1t_m' = beta1t_m, 
#'           'beta0s_m' = beta0s_m, 
#'           #'beta1s_m' = beta1s_m, 
#'           'beta2s_m' = beta2s_m, 
#'           'beta0t_p' = beta0t_p, 
#'           #'beta1t_p' = beta1t_p, 
#'           'beta0s_p' = beta0s_p, 
#'           #'beta1s_p' = beta1s_p, 
#'           'beta2s_p' = beta2s_p,
#'           'ones' = rep(1, total_subjects),
#'           'beta1t_r'      = beta1_shapeT,
#'           'beta1t_lambda' = beta1_rateT,
#'           'beta1s_r'      = beta1_shapeE,
#'           'beta1s_lambda' = beta1_rateE),
#'         init = list(
#'           'beta0t' = -2, 
#'           'beta1t' =  1, 
#'           'beta0s' = -1, 
#'           'beta1s' =  1, 
#'           'beta2s' =  0),
#'         n.chains = 1, quiet = TRUE))
#' 
#' 			jags_works <- 1*(length(jags_mod) > 1 | attempt == 10)
#' 
#' 			attempt <- attempt + 1
#' 
#' 		}
#' 		
#' 		update(jags_mod, 4000, quiet = TRUE)
#' 
#' 		coda_samples_temp <- coda.samples(jags_mod, c('pt_est', 'ps_est'), 1000, quiet = TRUE)
#' 		##coda_samples_temp <- coda.samples(jags_mod, c('pt_est', 'ps_est', 'beta0t', 'beta1t', 'beta0s', 'beta1s', 'beta2s'), 1000)
#' 
#' 		coda_samples <- coda_samples_temp[[1]]
#' 
#' 		dist_samples <- matrix(data = 0, ncol = dim(coda_samples)[2]/2, nrow = dim(coda_samples)[1])
#' 		
#' 		for(k in 1:dim(dist_samples)[2]) {
#' 
#' 			dist_samples[,k] <- calc_dist(ps = coda_samples[,k], pt = coda_samples[,(k + 4)])
#' 
#' 		}
#' 		
#' 		###prob_ps_g_psmin is the posterior probability that the probability of survival is greater than ps_min###
#' 		###prob_pt_l_ptmax is the posterior probability that the probability of toxicity is less than pt_max###
#' 
#' 		prob_acceptable <- colMeans(dist_samples > 0)
#' 
#' 		###determine acceptable doses###
#' 
#' 		acceptable <- 1*(prob_acceptable > gatekeeper)
#' 		if(sum(acceptable) == 0) {
#' 		  optimal_dose.perm <- 0
#' 		} else {
#'   		accept_dose <- dose_level[acceptable  == 1]
#'  
#' 	  	###Calculate desirability index for acceptable doses###
#' 
#' 		  ps_cur <- colMeans(coda_samples[,1:4])
#' 		  pt_cur <- colMeans(coda_samples[,5:8])
#' 
#' 		  dist_cur <- calc_dist(ps_cur[acceptable  == 1], pt_cur[acceptable  == 1])
#' 
#' 		  ###identify optimal dose ###
#' 
#' 		  optimal_dose.perm <- accept_dose[which.max(dist_cur)]
#' 		  
#' 		}
#' 		
#' 		optimal_dose.perm
#' }





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
	# Beta commensurate prior
		
  require(rjags)
  
	like_scenario <- 1 * (d.esl == 1 & t.esl == 1) + 2 * (d.esl == 1 & t.esl == 0) + 3 * (d.esl == 0 & t.esl == 1) + 4 * (d.esl == 0 & t.esl == 0)
	
	jags_works <- 0
	attempt <- 1
	
	# total_subjects <- max_n
	exp_n <- length(d.esl)
	dec_n <- length(d.exp)
	#N2 <- dec_n
	while(jags_works == 0 & attempt <= 10) {
		
		comensurate_data <- list(
      'optimal_dose' = optimal_dose, 
      'nu_eff'       = 0.5, 
      'u_eff'        = 4, 
      'spike_eff'    = 400, 
      'nu_tox'       = 0.5, 
      'u_tox'        = 2.25, 
      'spike_tox'    = 100, 
      # 'tox'          = t.exp, 
      'dose'         = dose, 
      # 'n_esc'        = total_subjects, 
		  'n_esc'        = exp_n, 
		  'n_exp'        = dec_n,
      'y_s_esc'      = 1 - d.esl, 
      'y_t_esc'      = t.esl, 
		  'y_s_exp'      = sum(1 - d.exp),
		  'y_t_exp'      = sum(t.exp),
      'beta0_muT'    = beta0t_m, 
      # 'beta1_muT'    = beta1t_m, 
      'beta0_muE'    = beta0s_m, 
      #'beta1_muE'    = beta1s_m, 
      'beta2_muE'    = beta2s_m, 
      'beta0_tauT'   = beta0t_p, 
      #'beta1_tauT'   = beta1t_p, 
      'beta0_tauE'   = beta0s_p, 
      #'beta1_tauE'   = beta1s_p,
      'beta2_tauE'   = beta2s_p,
		  'beta1_shapeT' = beta1_shapeT,
		  'beta1_rateT'  = beta1_rateT,
		  'beta1_shapeE' = beta1_shapeE,
		  'beta1_rateE'  = beta1_rateE)
		
		# jags_mod.2 <- try(jags.model('bugs/commensurate_model_beta2.bug', 
		#   data = comensurate_data,
		#   init = list(
		#     'beta_0T' = -2, 
		#     'beta_1T' =  1, 
		#     'beta_0E' = -1, 
		#     'beta_1E' =  1, 
		#     'beta_2E' =  0),
		# 	n.chains = 1, quiet = TRUE))
		# 
		# jags_works <- 1*(length(jags_mod.2) > 1 | attempt == 10)
		# 
		# attempt <- attempt + 1

		 jags_try <- try({
		   jags_mod.2 <- jags.model('bugs/commensurate_model_beta2.bug',
		     data = comensurate_data,
		     init = list(
		       'beta_0T' = -2,
		       'beta_1T' =  1,
		       'beta_0E' = -1,
		       'beta_1E' =  1,
		       'beta_2E' =  0),
			   n.chains = 1, quiet = TRUE)
		     update(jags_mod.2, 4000, quiet = TRUE)
		   }, silent = TRUE)

		if(class(jags_try) != "try-error") jags_works <- 1
		# jags_works <- 1*(length(jags_mod.2) > 1 | attempt == 10)

		attempt <- attempt + 1
			
	}
		
	# update(jags_mod.2, 4000, quiet = TRUE)
	
	coda_samples_temp.2 <- coda.samples(jags_mod.2, c('eff', 'tox'), 1000, 
	  quiet = TRUE)
	coda_samples.2 <- coda_samples_temp.2[[1]]
		
	# Beta Commensurate prior data
	ps_est <- mean(coda_samples.2[,1])
	pt_est <- mean(coda_samples.2[,2])
	ps_conf <- quantile(coda_samples.2[,1], c(0.025, 0.975))
	pt_conf <- quantile(coda_samples.2[,2], c(0.025, 0.975))
	Beta.est <- c(ps_est, ps_conf, pt_est, pt_conf)
		
	names(Beta.est) <- c(
	  paste0("cp_eff_", c("est", "lwr", "upr")),
	  paste0("cp_tox_", c("est", "lwr", "upr")))

	return(Beta.est)

}





### GRAND ESTIMATOR FUNCTION ###
 
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

# ESTIMATORS <- function(seed = NULL, d.esl, d.exp, t.esl, t.exp, optimal_dose, dose, 
#   surv_time_obs, tox_time_obs, surv_time_obs.exp, tox_time_obs.exp, mb){
# 	
#   if(!missing(seed)) set.seed(seed)
#   
# 	if (optimal_dose != 0){
# 		# EXPANSION ESTIMATES ONLY
# 		s.dec <- mean(1-d.exp)
# 		t.dec <- mean(t.exp)
# 		s.conf <- binconf(sum(1-d.exp), length(1-d.exp), method = "wilson")[2:3]
# 		t.conf <- binconf(sum(t.exp), length(t.exp), method = "wilson")[2:3]
# 		dec <- c(s.dec, s.conf, t.dec, t.conf)
# 		names(dec) <- c(
# 		  paste0("dec_eff_", c("est", "lwr", "upr")),
# 		  paste0("dec_tox_", c("est", "lwr", "upr"))
# 		)
# 		
# 		# POOLED ESTIMATES
# 		pooled <- Pool(d.esl, d.exp, t.esl, t.exp)
# 		
# 		# RAO-BLACKWELL ESTIMATES
# 		rb.est <- RaoBlackwell(d.esl, d.exp, t.esl, t.exp, dose, optimal_dose, 
# 		  tox_time_obs, surv_time_obs, tox_time_obs.exp, surv_time_obs.exp)
# 		
# 		# BETA COMMENSURATE PRIOR ESTIMATES
# 		bc.est <- BetaComP(d.esl, d.exp, t.esl, t.exp, optimal_dose, dose, surv_time_obs, tox_time_obs)
# 		
# 		# ESCALATION ESTIMATES ONLY
# 		s.escl <- mean(1 - d.esl)
# 		t.escl <- mean(t.esl)
# 		s.escl_conf <- binconf(sum(1 - d.esl), length(1 - d.esl), method = "wilson")[2:3]
# 		t.escl_conf <- binconf(sum(t.esl), length(t.esl), method = "wilson")[2:3]
# 		esl.est <- c(s.escl, s.escl_conf, t.escl, t.escl_conf)
# 		names(esl.est) <- c(
# 		  paste0("esl_eff_", c("est", "lwr", "upr")),
# 		  paste0("esl_tox_", c("est", "lwr", "upr"))
# 		)
# 		
# 		# MODEL-BASED ESTIMATES (EffTox estimates and respective CI are from SurvEffTox_sim)
# 		mb.est <- mb
# 		names(mb.est) <- c(
# 		  paste0("mb_eff_", c("est", "lwr", "upr")),
# 		  paste0("mb_tox_", c("est", "lwr", "upr"))
# 		)
# 		
# 		# MEM MODEL ESTIMATES
# 		# mem.mat <- sapply(seq(0.01, 1, by = 0.01), function(i) MEM.data(1-d.esl,1-d.exp, t.esl, t.exp, dose, optimal_dose,i))
# 		mem <- calc_MEM(
# 		  optimal_dose = optimal_dose,
# 		  eff_esc      = 1 - d.esl,
# 		  tox_esc      = t.esl,
# 		  dose_esc     = dose,
# 		  eff_exp      = 1 -d.exp,
# 		  tox_exp      = t.exp,
# 		  dose_exp     = rep(optimal_dose, length(t.exp))
# 		)
# 	} else {
# 		dec <- rep(NA, 6)
# 		pooled <- rep(NA, 6)
# 		rb.est <- rep(NA, 2)
# 		bc.est <- rep(NA, 6)
# 		esl.est <- rep(NA, 6)
# 		mb.est <- rep(NA, 6)
# 		#mem.mat <- matrix(NA, ncol = 100, nrow = 2)
# 		mem <- rep(NA, 6)
# 	}
# 	# return(c(optimal_dose, dec, pooled, rb.est, bc.est, esl.est, mb.est, mem))
# 	out <- c(optimal_dose = optimal_dose, dec, pooled, rb.est, bc.est, 
# 	  esl.est, mb.est, mem)
# 	out
# }

# ~~~~ estimate function ~~~~ ----
estimate <- function(esl_data, dec_data, include.ci) {
  require(Hmisc)
  
  # set.seed(seed)
  
  d.esl        <- esl_data$death
  t.esl        <- esl_data$tox
  dose         <- esl_data$dose
  sto.esl      <- esl_data$surv_time_obs
  tto.esl      <- esl_data$tox_time_obs
  mb           <- esl_data$mb
  
  optimal_dose <- esl_data$optimal_dose

  d.exp         <- dec_data$death
  t.exp         <- dec_data$tox
  sto.exp       <- dec_data$surv_time_obs
  tto.exp       <- dec_data$tox_time_obs
  data_gen_type <- dec_data$data_gen_type

  # these are the targets of inference
  prob_eff <- dec_data$prob_eff
  prob_tox <- dec_data$prob_tox
  
	# EXPANSION ESTIMATES ONLY
	s.dec <- mean(1 - d.exp)
	t.dec <- mean(t.exp)
	s.conf <- binconf(sum(1 - d.exp), length(1 - d.exp), method = "wilson")[2:3]
	t.conf <- binconf(sum(t.exp), length(t.exp), method = "wilson")[2:3]
	dec <- c(s.dec, s.conf, t.dec, t.conf)
	names(dec) <- c(
	  paste0("dec_eff_", c("est", "lwr", "upr")),
	  paste0("dec_tox_", c("est", "lwr", "upr"))
	)

	# POOLED ESTIMATES
	pooled <- Pool(d.esl, d.exp, t.esl, t.exp, dose, optimal_dose)

	# BETA COMMENSURATE PRIOR ESTIMATES
	bc.est <- BetaComP(d.esl, d.exp, t.esl, t.exp, optimal_dose, 
	  dose, sto.esl, tto.esl)
	
	# MEM MODEL ESTIMATES
	mem <- calc_MEM(
	  optimal_dose = optimal_dose,
	  eff_esc      = 1 - d.esl,
	  tox_esc      = t.esl,
	  dose_esc     = dose,
	  eff_exp      = 1 -d.exp,
	  tox_exp      = t.exp,
	  dose_exp     = rep(optimal_dose, length(t.exp))
	)  
	

	# RAO-BLACKWELL ESTIMATES
	rb.est <- RaoBlackwell(d.esl, d.exp, t.esl, t.exp, dose, optimal_dose,
		tto.esl, sto.esl, tto.exp, sto.exp, include.ci)


	# ESCALATION ESTIMATES ONLY
	s.escl <- mean(1 - d.esl)
	t.escl <- mean(t.esl)
	s.escl_conf <- binconf(sum(1 - d.esl), length(1 - d.esl), method = "wilson")[2:3]
	t.escl_conf <- binconf(sum(t.esl), length(t.esl), method = "wilson")[2:3]
	esl.est <- c(s.escl, s.escl_conf, t.escl, t.escl_conf)
	names(esl.est) <- c(
	  paste0("esl_eff_", c("est", "lwr", "upr")),
	  paste0("esl_tox_", c("est", "lwr", "upr"))
	)

	# MODEL-BASED ESTIMATES (EffTox estimates and respective CI are from SurvEffTox_sim)
	mb.est <- mb
	names(mb.est) <- c(
	  paste0("mb_eff_", c("est", "lwr", "upr")),
	  paste0("mb_tox_", c("est", "lwr", "upr"))
	)

	# return(c(optimal_dose, dec, pooled, rb.est, bc.est, esl.est, mb.est, mem))
	# out <- c(optimal_dose = optimal_dose, dec, pooled, rb.est, bc.est,
	#   esl.est, mb.est, mem, data_gen_type)
	out <- data.frame(optimal_dose, 
	  prob_eff,
	  prob_tox,
	  t(dec), 
	  t(pooled), 
	  t(rb.est), 
	  t(bc.est),
	  t(esl.est), 
	  t(mb.est), 
	  t(mem), 
	  data_gen_type)
	
	out

}

estimate_mem <- function(esl_data, dec_data, prior_weight) {
  
  # set.seed(seed)
  
  d.esl        <- esl_data$death
  t.esl        <- esl_data$tox
  dose         <- esl_data$dose
  sto.esl      <- esl_data$surv_time_obs
  tto.esl      <- esl_data$tox_time_obs
  mb           <- esl_data$mb
  
  optimal_dose <- esl_data$optimal_dose

  d.exp         <- dec_data$death
  t.exp         <- dec_data$tox
  sto.exp       <- dec_data$surv_time_obs
  tto.exp       <- dec_data$tox_time_obs
  data_gen_type <- dec_data$data_gen_type

  # these are the targets of inference
  prob_eff <- dec_data$prob_eff
  prob_tox <- dec_data$prob_tox
  
	# MEM MODEL ESTIMATES
	mem <- calc_MEM(
	  optimal_dose = optimal_dose,
	  eff_esc      = 1 - d.esl,
	  tox_esc      = t.esl,
	  dose_esc     = dose,
	  eff_exp      = 1 -d.exp,
	  tox_exp      = t.exp,
	  dose_exp     = rep(optimal_dose, length(t.exp)),
	  prior_weight = prior_weight
	)  
  
	out <- data.frame(optimal_dose, 
	  prob_eff,
	  prob_tox,
	  t(mem), 
	  data_gen_type)
	
	out

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

# EstiResults <- function(sim_results){
#   
#   n_sim <- length(sim_results)
#   
#   seed <- sapply(sim_results, '[[', "seed")
#   scenario <- sapply(sim_results, '[[', "scenario")
#   
#   dec <- lapply(sim_results, "[[", "dec")
#   
#   pt     <- sapply(dec, "[[", "pt")
#   pt.lwr <- sapply(dec, "[[", "pt.lwr")
#   pt.upr <- sapply(dec, "[[", "pt.upr")
#   pt.ith <- sapply(dec, "[[", "pt.ith")
#   
#   ps     <- sapply(dec, "[[", "ps")
#   ps.lwr <- sapply(dec, "[[", "ps.lwr")
#   ps.upr <- sapply(dec, "[[", "ps.upr")
#   ps.ith <- sapply(dec, "[[", "ps.ith")
# 	
# 	# optimal_dose <- sapply(1:n_sim, function(i) sim_results1[[i]]$dat[1])
# 	# d.esl <- sapply(1:n_sim, function(i) sim_results[[i]]$death)
# 	# t.esl <- sapply(1:n_sim, function(i) sim_results[[i]]$tox)
# 	# d.exp1 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$d1) 
# 	# t.exp1 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$t1)
# 	# d.exp2 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$d2) 
# 	# t.exp2 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$t2)
# 	# d.exp3 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$d3) 
# 	# t.exp3 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$t3)
# 	# d.exp4 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$d4) 
# 	# t.exp4 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$t4)
# 	# dose <- sapply(1:n_sim, function(i) sim_results[[i]]$dose)
# 	# surv_time_obs <- sapply(1:n_sim, function(i) sim_results[[i]]$surv_time_obs)
# 	# 
# 	# tox_time_obs1 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$tox_time_obs1)
# 	# tox_time_obs2 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$tox_time_obs2)
# 	# tox_time_obs3 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$tox_time_obs3)
# 	# tox_time_obs4 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$tox_time_obs4)
# 	# 
# 	# surv_time_obs1 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$surv_time_obs1)
# 	# surv_time_obs2 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$surv_time_obs2)
# 	# surv_time_obs3 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$surv_time_obs3)
# 	# surv_time_obs4 <- sapply(1:n_sim, function(i) sim_results[[i]]$dec$surv_time_obs4)
# 	# 
# 	# tox_time_obs <- sapply(1:n_sim, function(i) sim_results[[i]]$tox_time_obs)
# 	# mb <- sapply(1:n_sim, function(i) sim_results[[i]]$mb)
# 
# 	optimal_dose <- lapply(1:n_sim, function(i) sim_results1[[i]]$dat[1])
# 	d.esl <- lapply(1:n_sim, function(i) sim_results[[i]]$death)
# 	t.esl <- lapply(1:n_sim, function(i) sim_results[[i]]$tox)
# 	d.exp1 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$d1) 
# 	t.exp1 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$t1)
# 	d.exp2 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$d2) 
# 	t.exp2 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$t2)
# 	d.exp3 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$d3) 
# 	t.exp3 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$t3)
# 	d.exp4 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$d4) 
# 	t.exp4 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$t4)
# 	dose <- lapply(1:n_sim, function(i) sim_results[[i]]$dose)
# 	surv_time_obs <- lapply(1:n_sim, function(i) sim_results[[i]]$surv_time_obs)
# 
# 	tox_time_obs1 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$tox_time_obs1)
# 	tox_time_obs2 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$tox_time_obs2)
# 	tox_time_obs3 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$tox_time_obs3)
# 	tox_time_obs4 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$tox_time_obs4)
# 	
# 	surv_time_obs1 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$surv_time_obs1)
# 	surv_time_obs2 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$surv_time_obs2)
# 	surv_time_obs3 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$surv_time_obs3)
# 	surv_time_obs4 <- lapply(1:n_sim, function(i) sim_results[[i]]$dec$surv_time_obs4)
# 	
# 	tox_time_obs <- lapply(1:n_sim, function(i) sim_results[[i]]$tox_time_obs)
# 	mb <- lapply(1:n_sim, function(i) sim_results[[i]]$mb)
# 	
# 	# for debugging:
# 	foo <- mapply(ESTIMATORS,
# 	  d.esl             = d.esl[[8]],
# 	  d.exp             = d.exp1[[8]],
# 	  t.esl             = t.esl[[8]],
# 	  t.exp             = t.exp1[[8]],
# 	  optimal_dose      = optimal_dose[[8]],
# 	  dose              = dose[[8]],
# 	  surv_time_obs     = surv_time_obs[[8]],
# 	  tox_time_obs      = tox_time_obs[[8]],
# 	  surv_time_obs.exp = surv_time_obs1[[8]],
# 	  tox_time_obs.exp  = tox_time_obs1[[8]],
# 	  mb                = mb[[8]])
# 
# 	
# 	# no intertrial heterogeneity
# 	# output1 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp1[,i], 
# 	#   t.esl[,i], t.exp1[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], 
# 	#   tox_time_obs[,i], surv_time_obs1[,i], tox_time_obs1[,i], mb[,i]))
# 	# output1 <- mapply(ESTIMATORS,
# 	#   d.esl             = d.esl,
# 	#   d.exp             = d.exp1,
# 	#   t.esl             = t.esl,
# 	#   t.exp             = t.exp1,
# 	#   optimal_dose      = optimal_dose,
# 	#   dose              = dose,
# 	#   surv_time_obs     = surv_time_obs,
# 	#   tox_time_obs      = tox_time_obs,
# 	#   surv_time_obs.exp = surv_time_obs1,
# 	#   tox_time_obs.exp  = tox_time_obs1,
# 	#   mb                = mb)
# 	output1 <- mcmapply(failwith(NA, quiet(ESTIMATORS)),
# 	  seed              = seed,
# 	  d.esl             = d.esl, 
# 	  d.exp             = d.exp1, 
# 	  t.esl             = t.esl, 
# 	  t.exp             = t.exp1, 
# 	  optimal_dose      = optimal_dose, 
# 	  dose              = dose, 
# 	  surv_time_obs     = surv_time_obs, 
# 	  tox_time_obs      = tox_time_obs, 
# 	  surv_time_obs.exp = surv_time_obs1, 
# 	  tox_time_obs.exp  = tox_time_obs1, 
# 	  mb                = mb,
# 	  mc.cores          = n_cores)
# 
# 	# lower efficacy, higher toxicity
# 	# output2 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp2[,i], 
# 	#   t.esl[,i], t.exp3[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], 
# 	#   tox_time_obs[,i], surv_time_obs2[,i], tox_time_obs3[,i], mb[,i]))
# 	# output2 <- mapply(quiet(ESTIMATORS), 
# 	#   d.esl             = d.esl, 
# 	#   d.exp             = d.exp2, 
# 	#   t.esl             = t.esl, 
# 	#   t.exp             = t.exp3, 
# 	#   optimal_dose      = optimal_dose, 
# 	#   dose              = dose, 
# 	#   surv_time_obs     = surv_time_obs, 
# 	#   tox_time_obs      = tox_time_obs, 
# 	#   surv_time_obs.exp = surv_time_obs2, 
# 	#   tox_time_obs.exp  = tox_time_obs3, 
# 	#   mb                = mb)
# 	output2 <- mcmapply(failwith(NA, quiet(ESTIMATORS)),
# 	  seed              = seed,
# 	  d.esl             = d.esl, 
# 	  d.exp             = d.exp2, 
# 	  t.esl             = t.esl, 
# 	  t.exp             = t.exp3, 
# 	  optimal_dose      = optimal_dose, 
# 	  dose              = dose, 
# 	  surv_time_obs     = surv_time_obs, 
# 	  tox_time_obs      = tox_time_obs, 
# 	  surv_time_obs.exp = surv_time_obs2, 
# 	  tox_time_obs.exp  = tox_time_obs3, 
# 	  mb                = mb,
# 	  mc.cores          = n_cores)
# 	
# 	# intertrial heterogeneity
# 	# output3 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp4[,i], 
# 	#   t.esl[,i], t.exp4[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], 
# 	#   tox_time_obs[,i], surv_time_obs4[,i], tox_time_obs4[,i], mb[,i]))
# 	# output3 <- mapply(quiet(ESTIMATORS),
# 	#   d.esl             = d.esl, 
# 	#   d.exp             = d.exp4, 
# 	#   t.esl             = t.esl, 
# 	#   t.exp             = t.exp4, 
# 	#   optimal_dose      = optimal_dose, 
# 	#   dose              = dose, 
# 	#   surv_time_obs     = surv_time_obs, 
# 	#   tox_time_obs      = tox_time_obs, 
# 	#   surv_time_obs.exp = surv_time_obs4, 
# 	#   tox_time_obs.exp  = tox_time_obs4, 
# 	#   mb                = mb)
# 	output3 <- mcmapply(failwith(NA, quiet(ESTIMATORS)),
# 	  seed              = seed,
# 	  d.esl             = d.esl, 
# 	  d.exp             = d.exp4, 
# 	  t.esl             = t.esl, 
# 	  t.exp             = t.exp4, 
# 	  optimal_dose      = optimal_dose, 
# 	  dose              = dose, 
# 	  surv_time_obs     = surv_time_obs, 
# 	  tox_time_obs      = tox_time_obs, 
# 	  surv_time_obs.exp = surv_time_obs4, 
# 	  tox_time_obs.exp  = tox_time_obs4, 
# 	  mb                = mb,
# 	  mc.cores          = n_cores)
# 	
# 	o1 <- data.frame(do.call(rbind, output1))
# 	o2 <- data.frame(do.call(rbind, output2))
# 	o3 <- data.frame(do.call(rbind, output3))
# 	
# 	out1 <- data.frame(
# 	  seed = seed,
# 	  scenario = scenario,
# 	  data = "no_effect",
# 	  pTox = pt,
# 	  pEff = ps,
# 	  pTox_offset = 0,
# 	  pEff_offset = 0,
# 	  o1
# 	)
# 	
# 	out2 <- data.frame(
# 	  seed = seed,
# 	  scenario = scenario,
# 	  data = "hiTox_loEff",
# 	  pTox = pt.upr,
# 	  pEff = ps.lwr,
# 	  pTox_offset = qlogis(pt.upr) - qlogis(pt), # QA checking
# 	  pEff_offset = qlogis(ps.lwr) - qlogis(ps), # QA checking
# 	  o2
# 	)
# 	
# 	out3 <- data.frame(
# 	  seed = seed,
# 	  scenario = scenario,
# 	  data = "ITH",
# 	  pTox = pt.ith,
# 	  pEff = ps.ith,
# 	  pTox_offset = qlogis(pt.ith) - qlogis(pt),
# 	  pEff_offset = qlogis(ps.ith) - qlogis(ps),
# 	  o3
# 	)
# 	
# 	# output3 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp3[,i], t.esl[,i], t.exp3[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], tox_time_obs[,i], mb[,i]))
# 	
# 	# output4 <- sapply(1:n_sim, function(i) ESTIMATORS(d.esl[,i], d.exp4[,i], t.esl[,i], t.exp4[,i], optimal_dose[i], dose[,i], surv_time_obs[,i], tox_time_obs[,i], mb[,i]))
# 	
# 	# out <- list(output1 = output1, output2 = output2)
# 	
# 	out <- rbind(out1, out2, out3)
# 	
# 	out
# }


## Choose scenario ###
true.prob <- function(scenario){
	if (scenario == 1){
			b0t <- -3; b1t <- 1
			b0s <- -0.5; b1s <- 0.7
			
		} else if (scenario == 2){
			b0t <- -0.5; b1t <- 0.6
			b0s <- 1.2; b1s <- 0.3
		}else if (scenario == 3){
			b0t <- -3.5; b1t <- 0.9
			b0s <- -2; b1s <- 0.9
		} else if (scenario == 4) {
			b0t <- -1; b1t <- 0.5
			b0s <- -0.05; b1s <- 1.1
		}
		
	# pt <- round(inv.logit(b0t + b1t*0:3), 2)
	# ps <- round(inv.logit(b0s + b1s*0:3), 2)

  pt <- plogis(b0t + b1t*0:3)
	ps <- plogis(b0s + b1s*0:3)
	
  out <- list(pt = pt, ps = ps)
  out
}

###true probabilities of efficacy and toxicity###

###scenario 1###

#pt <- c(.05, .12, .27, .50) ###beta0 = -3, beta1 = 1###
#ps <- c(.38, .55, .71, .83) ###beta0 = -0.5, beta1 = 0.7###

###scenario 2###

##pt <- c(.38, .52, .67, .79) ###beta0 = -.5, beta1 = .6###
##ps <- c(.77, .82, .86, .89) ###beta0 = 1.2, beta1 = .3###

###scenario 3######

##pt <- c(.02, .07, .15, .31) ###beta0 = -3.5, beta1 = .9###
##ps <- c(.12, .25, .45, .67) ###beta0 = -2, beta1 = .9###

### scenario 4###

##pt <- c(0.27, 0.38,0.50, 0.62) ####beta0 = -1, beta1 = 0.5
##ps <- c(0.49, 0.74, 0.90, 0.96) #### beta0 = -0.05, beta = 1.1



calc_dist <- function(ps, pt) {

  dist <- 1 - (((1 - ps) / (1 - ps_min)) ^ p + (pt / pt_max) ^ p ) ^ (1 / p)

  dist

}


# mean parameters  
# beta0_muT <- -3
# beta1_muT <-  1
# beta0_muE <- -1
# beta1_muE <-  1
# beta2_muE <-  0
  
# sd parameters
# beta0_sdT <- 3
# beta1_sdT <- 2
# beta0_sdE <- 3
# beta1_sdE <- 2
# beta2_sdE <- 0.5

# precision parameters
# beta0_tauT <- 1 / (beta0_sdT ^ 2)
# beta1_tauT <- 1 / (beta1_sdT ^ 2)
# beta0_tauE <- 1 / (beta0_sdE ^ 2)
# beta1_tauE <- 1 / (beta1_sdE ^ 2)
# beta2_tauE <- 1 / (beta2_sdE ^ 2)
#   
# a_dE <- 1
# b_dE <- 1
# a_dT <- 1
# b_dT <- 1

expMean <- function(x) {
  x <- x[!is.na(x)]
  cons <- median(x)
  exp(-cons) * mean(exp(cons + x))
}


log_likelihood <- function(beta_0E, beta_1E, beta_2E, beta_0T, beta_1T, p_eff, p_tox, dose_esc, dose_exp, eff_esc, eff_exp, tox_esc, tox_exp,
  share_eff, share_tox) {
  
  # esc denotes escalation cohort
  # exp denotes expansion cohort
  
  force(share_eff)
  force(share_tox)
  
  # likelihood for efficacy
  lp_eff_esc <- c(cbind(1, dose_esc - 1, (dose_esc - 1) ^ 2) %*% c(beta_0E, beta_1E, beta_2E))
  lp_eff_exp <- c(cbind(1, dose_exp - 1, (dose_exp - 1) ^ 2) %*% c(beta_0E, beta_1E, beta_2E))
  
  prob_eff_esc <- plogis(lp_eff_esc)
  prob_eff_exp <- share_eff * plogis(lp_eff_exp)  + (1 - share_eff) * p_eff
  
  loglike_eff_esc <- sum(eff_esc * log(prob_eff_esc) + (1 - eff_esc) * log(1 - prob_eff_esc))
  loglike_eff_exp <- sum(eff_exp * log(prob_eff_exp) + (1 - eff_exp) * log(1 - prob_eff_exp))

  # likelihood for toxicty
  lp_tox_esc <- c(cbind(1, dose_esc - 1) %*% c(beta_0T, beta_1T))
  lp_tox_exp <- c(cbind(1, dose_exp - 1) %*% c(beta_0T, beta_1T))
  
  prob_tox_esc <- plogis(lp_tox_esc)
  prob_tox_exp <- share_tox * plogis(lp_tox_exp)  + (1 - share_tox) * p_tox
  
  loglike_tox_esc <- sum(tox_esc * log(prob_tox_esc) + (1 - tox_esc) * log(1 - prob_tox_esc))
  loglike_tox_exp <- sum(tox_exp * log(prob_tox_exp) + (1 - tox_exp) * log(1 - prob_tox_exp))
  
  out  <- loglike_eff_esc + loglike_eff_exp + loglike_tox_esc + loglike_tox_exp
  out
}

calc_MEM <- function(optimal_dose, eff_esc, tox_esc, dose_esc, eff_exp, tox_exp, dose_exp, prior = c("equal", "eb", "eb_constrain"), constraint=1, prior_weight = 0.5, ...){

  prior <- match.arg(prior)
  if(prior != "equal") stop("prior != 'equal' is not yet implemented")
  

  # MEM definitions ----
  # mem_1 no borrowing
	# mem_2 pool efficacy, not toxicity
	# mem_3 pool toxicity, not efficacy
	# mem_4 pool all
  
  # E is for efficacy, T is for Toxicty


	# marginal likelihood ----
	# compute marginal likelihood via Monte Carlo integration
	# using random draws from the prior

	# set.seed(1234)
  
  rand_beta_0E <- beta0s_m + 1 / (sqrt(beta0s_p)) * rnorm(1e4) 
	#rand_beta_1E <- beta1s_m + 1 / (sqrt(beta1s_p)) * rnorm(1e4)
  rand_beta_1E <- rgamma(1e4, beta1_shapeE, beta1_rateE)
	rand_beta_2E <- beta2s_m + 1 / (sqrt(beta2s_p)) * rnorm(1e4) 
	
	rand_beta_0T <- beta0t_m + 1 / (sqrt(beta0t_p)) * rnorm(1e4) 
	#rand_beta_1T <- beta1t_m + 1 / (sqrt(beta1t_p)) * rnorm(1e4) 
	rand_beta_1T <- rgamma(1e4, beta1_shapeT, beta1_rateT)
	
  # rand_p_eff <- rbeta(1e4, a_dE, b_dE)
  # rand_p_tox <- rbeta(1e4, a_dT, b_dT)
  rand_p_eff <- rbeta(1e4, 1, 1)
  rand_p_tox <- rbeta(1e4, 1, 1)
  
  # previously used gamma priors
  # rand_beta_1E <- rgamma(1e4, beta1_shapeE, beta1_rateE)
  # rand_beta_1T <- rgamma(1e4, beta1_shapeT, beta1_rateT)

	mem_1_loglike_dist <- mapply(log_likelihood,
	  beta_0E   = rand_beta_0E,
	  beta_1E   = rand_beta_1E,
	  beta_2E   = rand_beta_2E,
	  beta_0T   = rand_beta_0T,
	  beta_1T   = rand_beta_1T,
	  p_eff     = rand_p_eff,
	  p_tox     = rand_p_tox,
	  dose_esc  = list(dose_esc),
	  dose_exp  = list(dose_exp),
	  eff_esc   = list(eff_esc),
	  eff_exp   = list(eff_exp),
	  tox_esc   = list(tox_esc),
	  tox_exp   = list(tox_exp),
	  share_eff = FALSE,
	  share_tox = FALSE)
	mem_2_loglike_dist <- mapply(log_likelihood,
	  beta_0E   = rand_beta_0E,
	  beta_1E   = rand_beta_1E,
	  beta_2E   = rand_beta_2E,
	  beta_0T   = rand_beta_0T,
	  beta_1T   = rand_beta_1T,
	  p_eff     = rand_p_eff,
	  p_tox     = rand_p_tox,
	  dose_esc  = list(dose_esc),
	  dose_exp  = list(dose_exp),
	  eff_esc   = list(eff_esc),
	  eff_exp   = list(eff_exp),
	  tox_esc   = list(tox_esc),
	  tox_exp   = list(tox_exp),
	  share_eff = TRUE,
	  share_tox = FALSE)
	mem_3_loglike_dist <- mapply(log_likelihood,
	  beta_0E   = rand_beta_0E,
	  beta_1E   = rand_beta_1E,
	  beta_2E   = rand_beta_2E,
	  beta_0T   = rand_beta_0T,
	  beta_1T   = rand_beta_1T,
	  p_eff     = rand_p_eff,
	  p_tox     = rand_p_tox,
	  dose_esc  = list(dose_esc),
	  dose_exp  = list(dose_exp),
	  eff_esc   = list(eff_esc),
	  eff_exp   = list(eff_exp),
	  tox_esc   = list(tox_esc),
	  tox_exp   = list(tox_exp),
	  share_eff = FALSE,
	  share_tox = TRUE)
	mem_4_loglike_dist <- mapply(log_likelihood,
	  beta_0E   = rand_beta_0E,
	  beta_1E   = rand_beta_1E,
	  beta_2E   = rand_beta_2E,
	  beta_0T   = rand_beta_0T,
	  beta_1T   = rand_beta_1T,
	  p_eff     = rand_p_eff,
	  p_tox     = rand_p_tox,
	  dose_esc  = list(dose_esc),
	  dose_exp  = list(dose_exp),
	  eff_esc   = list(eff_esc),
	  eff_exp   = list(eff_exp),
	  tox_esc   = list(tox_esc),
	  tox_exp   = list(tox_exp),
	  share_eff = TRUE,
	  share_tox = TRUE)
	
	# summary(mem_1_loglike_dist)
	# summary(mem_2_loglike_dist)
	# summary(mem_3_loglike_dist)
	# summary(mem_4_loglike_dist)

	mem_1_like <- expMean(mem_1_loglike_dist)
	mem_2_like <- expMean(mem_2_loglike_dist)
	mem_3_like <- expMean(mem_3_loglike_dist)
	mem_4_like <- expMean(mem_4_loglike_dist)
	
	marginal_like <- c(
	  mem_1_like,
	  mem_2_like,
	  mem_3_like,
	  mem_4_like)
	
	pweights <- c(
	  (1 - prior_weight) ^ 2,
	  (1 - prior_weight) * prior_weight,
	  (1 - prior_weight) * prior_weight,
	  prior_weight ^ 2)

		
	omega <- (pweights * marginal_like) / sum(pweights * marginal_like)

	# posterior computations ----	
	eff_post_samples <- eff_post_fun(optimal_dose, eff_esc, eff_exp, dose_esc, dose_exp)
	tox_post_samples <- tox_post_fun(optimal_dose, tox_esc, tox_exp, dose_esc, dose_exp)	

	eff_posteriors <- array(dim = c(1e4, 4),
	  dimnames = list(sample = 1:1e4, mem = paste0("mem_", 1:4)))
	tox_posteriors <- array(dim = c(1e4, 4),
	  dimnames = list(sample = 1:1e4, mem = paste0("mem_", 1:4)))

	# mem_1 no borrowing
	# mem_2 pool efficacy, not toxicity
	# mem_3 pool toxicity, not efficacy
	# mem_4 pool all

	eff_posteriors[, c(2, 4)] <- eff_post_samples[, 'post_sh']
  eff_posteriors[, c(1, 3)] <- eff_post_samples[, 'post_ns']
  
  tox_posteriors[, c(1, 2)] <- tox_post_samples[, 'post_ns']
  tox_posteriors[, c(3, 4)] <- tox_post_samples[, 'post_sh']
  
  
  # head(eff_posteriors)
  # head(tox_posteriors)
 
  # apply(eff_posteriors, 2, sd)
 	# apply(tox_posteriors, 2, sd)
	
	which.mem <- cbind(1:1e4, sample(1:4, 1e4, replace = TRUE, prob = omega))
	
	eff_post <- eff_posteriors[which.mem]
  tox_post <- tox_posteriors[which.mem]	
	
  out <- c(
    mem_tox_est = mean(tox_post),
    mem_tox_lwr = unname(quantile(tox_post, 0.025)),
    mem_tox_upr = unname(quantile(tox_post, 0.975)),
    mem_eff_est = mean(eff_post),
    mem_eff_lwr = unname(quantile(eff_post, 0.025)),
    mem_eff_upr = unname(quantile(eff_post, 0.975)))
	
  out
}

eff_post_fun <- function(optimal_dose, eff_esc, eff_exp, dose_esc, dose_exp) {
  require(rjags)
  
  # jags model: posterior with shared parameters
  dlist <- list(optimal_dose = optimal_dose,
    n         = length(eff_esc) + length(eff_exp),
    eff       = c( eff_esc,  eff_exp),
    dose      = c(dose_esc, dose_exp),
    beta0_muE     = beta0s_m,
    # beta1_muE     = beta1_muE,
    beta2_muE     = beta2s_m,
    beta1_shapeE  = beta1_shapeE,
    beta0_tauE    = beta0s_p,
    #beta1_tauE    = beta1_tauE,
    beta2_tauE    = beta2s_p,
    # tau_1E   = tau_1E,
    beta1_rateE   = beta1_rateE
    # tau_2E    = tau_2E,
    )
	jags_mod <- jags.model('bugs/eff_logistic.bug', 
	  data     = dlist,
    n.chains = 1, 
		quiet    = TRUE)
		
	update(jags_mod, 1e4, quiet = TRUE)

	coda_samples_temp <- coda.samples(model = jags_mod, 
	  variable.names = c('beta_0E', 'beta_1E', 'beta_2E', 'p_eff'), 
		n.iter         = 1e4, 
		quiet          = TRUE)
  coda_samples <- coda_samples_temp[[1]]
  post_sh <- coda_samples[, 'p_eff']

  # beta posterior: no shared parameters
  # dist'n is available directly, but sampling allows for easier posterior 
  # computaion with MEMs
  n <- length(eff_exp)
  # post_ns <- rbeta(1e4, sum(eff_exp) + a_dE, n - sum(eff_exp) + b_dE)
  post_ns <- rbeta(1e4, sum(eff_exp) + 1, n - sum(eff_exp) + 1)
  post <- structure(cbind(post_sh, post_ns),
    colnames = c("sharing", "no_sharing"))
  post
}

tox_post_fun <- function(optimal_dose, tox_esc, tox_exp, dose_esc, dose_exp) {
  require(rjags)
  
  # jags model: posterior with shared parameters
  dlist <- list(optimal_dose = optimal_dose,
    n      = length(tox_esc) + length(tox_exp),
    tox    = c( tox_esc,  tox_exp),
    dose   = c(dose_esc, dose_exp),
    beta0_muT  = beta0t_m,
    # beta1_muT  = beta1_muT,
    beta1_shapeT = beta1_shapeT,
    beta0_tauT = beta0t_p,
    # beta1_tauT = beta1_tauT
    #tau_1T = tau_1T
    beta1_rateT = beta1_rateT
    )
	jags_mod <- jags.model('bugs/tox_logistic.bug', 
	  data     = dlist,
    n.chains = 1, 
		quiet    = TRUE)
		
	update(jags_mod, 1e4, quiet = TRUE)

	coda_samples_temp <- coda.samples(model = jags_mod, 
	  variable.names = c('beta_0T', 'beta_1T','p_tox'), 
		n.iter         = 1e4, 
		quiet          = TRUE)
  coda_samples <- coda_samples_temp[[1]]
  post_sh <- coda_samples[, 'p_tox']

  # beta posterior: no shared parameters
  # dist'n is available directly, but sampling allows for easier posterior 
  # computaion with MEMs
  n <- length(tox_exp)
  # post_ns <- rbeta(1e4, sum(tox_exp) + a_dT, n - sum(tox_exp) + b_dT)
  post_ns <- rbeta(1e4, sum(tox_exp) + 1, n - sum(tox_exp) + 1)
  post <- structure(cbind(post_sh, post_ns),
    colnames = c("sharing", "no_sharing"))
  post
}

quiet <- function(fun) {
  function(...) {
    capture.output(out <- fun(...))
    out
  }
}

failwith <- function (default = NULL, f, quiet = TRUE) 
{
    function(...) {
        out <- default
        try(out <- f(...), silent = quiet)
        out
    }
}

# function to simulate data ----

simulate_data <- function(scenario, data_set_number, keep_dose0 = FALSE) {
  
  # tox and eff probabilities
  ps <- true.prob(scenario)$ps
  pt <- true.prob(scenario)$pt
  
  # simulate escalation cohort
  esl <- SurvEffTox_sim(ps, pt, keep_dose0)

  
  optimal_dose <- esl$optimal_dose  
  
  if(optimal_dose == 0) {
    dec1 <- dec2 <- dec3 <- NA
  } else {
    dec1 <- DEC(optimal_dose, ps, pt, data_gen_type = "no_effect")
    dec2 <- DEC(optimal_dose, ps, pt, data_gen_type = "hiTox_loEff")
    dec3 <- DEC(optimal_dose, ps, pt, data_gen_type = "ITH")
  }
  
  # o1 <- estimate(esl, DEC(optimal_dose, ps, pt, "no_effect"))
  # o2 <- estimate(esl, DEC(optimal_dose, ps, pt, "hiTox_loEff"))
  # o3 <- estimate(esl, DEC(optimal_dose, ps, pt, "ITH"))
  
  # oo <- do.call(rbind, list(o1, o2, o3))
  
  # out <- list(
  #   seed          = seed, 
  #   scenario      = scenario, 
  #   data_gen_type = data_gen_type,
  #   optimal_dose  = optimal_dose,
  #   prob_eff      = prob_eff,
  #   prob_tox      = prob_tox)
  
  # out <- data.frame(seed, scenario, oo)
  
  out <- list(scenario = scenario, 
    data_set_number = data_set_number, 
    esl             = esl, 
    dec1            = dec1, 
    dec2            = dec2,
    dec3            = dec3,
    optimal_dose    = optimal_dose)
  
  out
  
}

# function to apply estimators ----
apply_estimators <- function(data_set, include.ci = FALSE) {
  o1 <- estimate(data_set$esl, data_set$dec1, include.ci = include.ci)
  o2 <- estimate(data_set$esl, data_set$dec2, include.ci = include.ci)
  o3 <- estimate(data_set$esl, data_set$dec3, include.ci = include.ci)
  oo <- rbind(o1, o2, o3)
  out <- data.frame(
    data_set_num = data_set$data_set_number,
    scenario     = data_set$scenario,
    oo)
  out
}

apply_mem <- function(data_set, prior_weight) {
  o1 <- estimate_mem(data_set$esl, data_set$dec1, prior_weight)
  o2 <- estimate_mem(data_set$esl, data_set$dec2, prior_weight)
  o3 <- estimate_mem(data_set$esl, data_set$dec3, prior_weight)
  oo <- rbind(o1, o2, o3)
  out <- data.frame(
    data_set_num = data_set$data_set_number,
    scenario     = data_set$scenario,
    prior_weight = prior_weight,
    oo)
  out
}

# for debugging:
# seed <- 1
# scenario <- 1

# debug(sim_and_estimate)
# debug(estimate)
# debug(DEC)
# debug(RaoBlackwell)
# foo <- quiet(sim_and_estimate)(seed, scenario, data_gen_type)
# foo <- sim_and_estimate(seed, scenario)
