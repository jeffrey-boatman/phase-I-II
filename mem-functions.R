# mean parameters  
beta0_muT <- -3
beta1_muT <-  1
beta0_muE <- -1
beta1_muE <-  1
beta2_muE <-  0
  
# sd parameters
beta0_sdT <- 3
beta1_sdT <- 2
beta0_sdE <- 3
beta1_sdE <- 2
beta2_sdE <- 0.5

# precision parameters
beta0_tauT <- 1 / (beta0_sdT ^ 2)
beta1_tauT <- 1 / (beta1_sdT ^ 2)
beta0_tauE <- 1 / (beta0_sdE ^ 2)
beta1_tauE <- 1 / (beta1_sdE ^ 2)
beta2_tauE <- 1 / (beta2_sdE ^ 2)
  
a_dE <- 1
b_dE <- 1
a_dT <- 1
b_dT <- 1

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

	set.seed(1234)
	
  rand_beta_0E <- beta0_muE + 1 / (sqrt(beta0_tauE)) * rnorm(1e4) 
	rand_beta_1E <- beta1_muE + 1 / (sqrt(beta1_tauE)) * rnorm(1e4)
	rand_beta_2E <- beta2_muE + 1 / (sqrt(beta2_tauE)) * rnorm(1e4) 
	
	rand_beta_0T <- beta0_muT + 1 / (sqrt(beta0_tauT)) * rnorm(1e4) 
	rand_beta_1T <- beta1_muT + 1 / (sqrt(beta1_tauT)) * rnorm(1e4) 
	
  rand_p_eff <- rbeta(1e4, a_dE, b_dE)
  rand_p_tox <- rbeta(1e4, a_dT, b_dT)
  
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
    p_tox   = mean(tox_post),
    tox_lwr = unname(quantile(tox_post, 0.025)),
    tox_upr = unname(quantile(tox_post, 0.975)),
    p_eff   = mean(eff_post),
    eff_lwr = unname(quantile(eff_post, 0.025)),
    eff_upr = unname(quantile(eff_post, 0.975)))
	
  out
}

eff_post_fun <- function(optimal_dose, eff_esc, eff_exp, dose_esc, dose_exp) {
  require(rjags)
  
  # jags model: posterior with shared parameters
  dlist <- list(optimal_dose = optimal_dose,
    n         = length(eff_esc) + length(eff_exp),
    eff       = c( eff_esc,  eff_exp),
    dose      = c(dose_esc, dose_exp),
    beta0_muE     = beta0_muE,
    beta1_muE     = beta1_muE,
    beta2_muE     = beta2_muE,
    # beta1_shapeE  = beta1_shapeE,
    beta0_tauE    = beta0_tauE,
    beta1_tauE    = beta1_tauE,
    beta2_tauE    = beta2_tauE
    # tau_1E   = tau_1E,
    # beta1_rateE   = beta1_rateE,
    # tau_2E    = tau_2E,
    )
	jags_mod <- jags.model('eff_logistic.bug', 
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
  post_ns <- rbeta(1e4, sum(eff_exp) + a_dE, n - sum(eff_exp) + b_dE)
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
    beta0_muT  = beta0_muT,
    beta1_muT  = beta1_muT,
    # beta1_shapeT = beta1_shapeT,
    beta0_tauT = beta0_tauT,
    beta1_tauT = beta1_tauT
    #tau_1T = tau_1T
    # beta1_rateT = beta1_rateT
    )
	jags_mod <- jags.model('tox_logistic.bug', 
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
  post_ns <- rbeta(1e4, sum(tox_exp) + a_dT, n - sum(tox_exp) + b_dT)
  post <- structure(cbind(post_sh, post_ns),
    colnames = c("sharing", "no_sharing"))
  post
}

