
# Purpose: Produce latex table of summary statistics for each combining esimates by survival and toxicity



################ Arrange Simulation Results by Dose (1-4) in Matrix Form ################

# PURPOSE: Obtained organized summary statistics by dose for each estimator 

# INPUT: 
## output: Survival and toxicity estimates from EstiResults Function
## name: Name the esimtator desired
## prob.wt: Vector of probability of selecting a dose used as probability weights
## s.truth: True survival probabilities for each dose
## t.truth: True toxicity probabilities for each dose
## ITH: Indictor for inter-trial heterogeneity present ("yes") or not ("no")
## n_sim: Number of simulations

# OUTPUT: 
## ps: list of average, SD, bias, MSE for each dose and a weighted average for survival
## pt: list of average, SD, bias, MSE for each dose and a weighted average for toxcity

SimEval <- function(output, name = c("exp", "pool", "rao", "beta", "esl", "mb"), prob.wt, s.truth, t.truth, ITH = "no", n_sim = 1000){
	
	if(name == "exp") {s.index = 2; t.index = (s.index + 3)}
	if(name == "pool") {s.index = 8; t.index = (s.index + 3)}
	if(name == "rao") {s.index = 14; t.index = (s.index + 1)}
	if(name == "beta") {s.index = 16; t.index = (s.index + 3)}
	if(name == "esl") {s.index = 22; t.index = (s.index + 3)}
	if(name == "mb") {s.index = 28; t.index = (s.index + 3)}
	
	if (ITH == "no"){
		# Pulls matrix for a particular estimator by optimal dose (i.e 1 - 4)
		s.mat <- ByDoseMat(output, s.index, n_sim)
		t.mat <- ByDoseMat(output, t.index, n_sim)
		
		# Evaluate the statistics by obtaining the average, SD, bias, and MSE for efficacy and toxicity by dose
		s.res <- StatProp(s.mat, s.truth, prob.wt)
		t.res <- StatProp(t.mat, t.truth, prob.wt)
	} else if (ITH == "yes") {
		# Pulls matrix for a particular estimator by optimal dose (i.e 1 - 4)
		s.mat <- ByDoseMat(output, s.index, n_sim)
		t.mat <- ByDoseMat(output, t.index, n_sim)
		
		# Evaluate the statistics by obtaining the average, SD, bias, and MSE for efficacy and toxicity by dose
		s.res <- JStatProp(s.mat, s.truth, prob.wt)
		t.res <- JStatProp(t.mat, t.truth, prob.wt)
				
	}

	return(list(ps = s.res, pt = t.res))	
	
}


## PURPOSE: Take simulated data and organize by dose (Function used by SimEval)
ByDoseMat <- function(output, index, n_sim = 1000){
	
	# Create an empty mat to hold simulation results
	mat <- matrix(NA, ncol = 4, nrow = n_sim, byrow=T)
	
	# The First Column in 'output' is the optimal dose (See ESTIMATOR function)
	for (i in 1:n_sim) {
		if (output[1,i] == 1){
		mat[i,1] <- output[index, i]
		} else if (output[1,i] == 2) {
		mat[i,2] <- output[index, i]
		} else if (output[1,i] == 3) {
		mat[i,3] <- output[index, i]
		} else  if(output[1,i] ==4) {
		mat[i,4] <- output[index, i]
		} }
	
	return(mat)
	
}

## PURPOSE: Take organized data and finds the average, SD, bias, and MSE; along with weighted averages of those (Function used by SimEval)
StatProp <- function(mat, truth, prob.wt){
	
	# Average, sd, bias, and MSE of estimates
	avg <- apply(mat, 2, mean, na.rm=TRUE)
	std <- apply(mat, 2, sd, na.rm=TRUE)
	bias <- avg - truth
	mse <- std^2 + bias^2
	
	# Obtain weighted average, sd, bias, and MSE
	w.avg <- sum(avg * prob.wt, na.rm=T)
	w.std <- sum(std * prob.wt, na.rm=T)
	w.bias <- sum(abs(bias) * prob.wt, na.rm=T)
	w.mse <- sum(mse * prob.wt, na.rm=T)
	
	# Consolidate dose-level (average, sd, bias, and mse) with weighted values
	new.avg <- c(avg, w.avg)
	new.std <- c(std, w.std)
	new.bias <- c(bias, w.bias)
	new.mse <- c(mse, w.mse)
	
	return(list(avg = new.avg, std = new.std, bias = new.bias, mse = new.mse))
}

## PURPOSE: Take organized jittered data and finds the average, SD, bias, and MSE; along with weighted averages of those (Function used by SimEval)
JStatProp <- function(mat, truth, prob.wt){
	
	# Find bias and Bias^2 for each trial 
	trial_bias <- mat - truth
	trial_mse <- (trial_bias)^2
	
	# Bias and MSE of estimates
	bias <- apply(trial_bias, 2, mean, na.rm = T)
	mse <- apply(trial_mse, 2, mean, na.rm = T)
	

	# Obtain weighted average, sd, bias, and MSE
	w.bias <- sum(abs(bias) * prob.wt, na.rm=T)
	w.mse <- sum(mse * prob.wt, na.rm=T)
	
	# Consolidate dose-level (bias, and mse) with weighted values
	new.bias <- c(bias, w.bias)
	new.mse <- c(mse, w.mse)
	
	return(list(bias = new.bias, mse = new.mse))
}

## PURPOSE: Takes simulated data from EffToxSim to find probability selection for each dose
ProbWt <- function(sim_results){
	
	results <- sapply(1:n_sim, function(i) sim_results[[i]]$dat)
	prob.selc <- c(mean(results[1,] == 1), mean(results[1,] == 2), mean(results[1,] == 3), mean(results[1,] == 4))
	return(prob.selc)

}

## PURPOSE: Gives true probability values for efficacy and toxicity
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

#################################################################################


################ Final Results Matrix Table ################

# PURPOSE: Using SimEval, ResultsTable produces three tables: Simulations Subject Results Table, Survival Results Table, and Toxicity Results Table.

# INPUT: 
## output: Survival and toxicity estimates from EstiResults Function
## prob.wt: Vector of probability of selecting a dose used as probability weights
## s.truth: True survival probabilities for each dose
## t.truth: True toxicity probabilities for each dose
## ITH: Indictor for inter-trial heterogeneity present ("yes") or not ("no")
## n_sim: Number of simulations

# OUTPUT: 
## results_table: Simulation Subject Results Table
## Eff_table: Survival Results Table Organized by estimator
## Tox_table: Toxicity Results Table Organized by estimator


ResultsTable <- function(output, ps, pt, prob.wt, results, ITH = "no", n_sim = 1000){
	
	# Results Matrix for holding results
	results_table <- matrix(NA, nrow = 5, ncol = 7)
	
	if (ITH == "no"){
		
		mb_est <- SimEval(output, "mb", prob.wt, s.truth = ps, t.truth = pt, n_sim = n_sim)
		esl_est <- SimEval(output, "esl", prob.wt, s.truth = ps, t.truth = pt, n_sim = n_sim)
		exp_est <- SimEval(output, "exp", prob.wt, s.truth = ps, t.truth = pt, n_sim = n_sim)
		pool_est <- SimEval(output, "pool", prob.wt, s.truth = ps, t.truth = pt, n_sim = n_sim)
		rao_est <- SimEval(output, "rao", prob.wt, s.truth = ps, t.truth = pt, n_sim = n_sim)
		beta_est <- SimEval(output, "beta", prob.wt, s.truth = ps, t.truth = pt, n_sim = n_sim)
		
		results_table[1, 2:5] <- pt
		results_table[2, 2:5] <- ps
		results_table[3, 2:5] <- round(calc_dist(ps, pt), 2)
		
		results_table[4,] <- c(round(c(mean(results[1,] == 0), mean(results[1,] == 1), mean(results[1,] == 2), mean(results[1,] == 3), mean(results[1,] == 4)), 2),NA , round(mean(results[6,]), 1))
		results_table[5,2:5] <- round(c(mean(results[2,]), mean(results[3,]), mean(results[4,]), mean(results[5,])), 2)
		est_name <- c("Model-Based Estimates", "Escalation Only","Expansion Only", "Pooled Estimator", "Rao-Blackwell Estimator", "Beta Commensurate Estimator")
		names_row <- c(sapply(est_name, function(name) c(paste(name),"Mean","Bias", "Standard Deviation", "Mean Square Error")))
		
		# Efficacy Table
		Eff_table <- matrix(NA, nrow = 30, ncol = 5)
		
		# Model-based for Efficacy: mean, bias, SD, and MSE
		Eff_table[2, 1:5] <- round(mb_est$ps$avg, 4)
		Eff_table[3, 1:5] <- round(mb_est$ps$bias, 4)
		Eff_table[4, 1:5] <- round(mb_est$ps$std, 4)
		Eff_table[5, 1:5] <- round(mb_est$ps$mse, 4)
	
		# Dose Escalation Cohort for Efficacy: mean, bias, SD, and MSE
		Eff_table[7, 1:5] <- round(esl_est$ps$avg, 4)
		Eff_table[8, 1:5] <- round(esl_est$ps$bias, 4)
		Eff_table[9, 1:5] <- round(esl_est$ps$std, 4)
		Eff_table[10, 1:5] <- round(esl_est$ps$mse, 4)
		
		# Dose Expansion Cohort for Efficacy: mean, bias, SD, and MSE
		Eff_table[12, 1:5] <- round(exp_est$ps$avg, 4)
		Eff_table[13, 1:5] <- round(exp_est$ps$bias, 4)
		Eff_table[14, 1:5] <- round(exp_est$ps$std, 4)
		Eff_table[15, 1:5] <- round(exp_est$ps$mse, 4)
		
		# Pooled data for Efficacy: mean, bias, SD, and MSE
		Eff_table[17, 1:5] <- round(pool_est$ps$avg, 4)
		Eff_table[18, 1:5] <- round(pool_est$ps$bias, 4)
		Eff_table[19, 1:5] <- round(pool_est$ps$std, 4)
		Eff_table[20, 1:5] <- round(pool_est$ps$mse, 4)
		
		# Rao Blackwell for Efficacy: mean, bias, SD, and MSE
		Eff_table[22, 1:5] <- round(rao_est$ps$avg, 4)
		Eff_table[23, 1:5] <- round(rao_est$ps$bias, 4)
		Eff_table[24, 1:5] <- round(rao_est$ps$std, 4)
		Eff_table[25, 1:5] <- round(rao_est$ps$mse, 4)
		
		# Beta Commensurate Prior for Efficacy: mean, bias, SD, and MSE
		Eff_table[27, 1:5] <- round(beta_est$ps$avg, 4)
		Eff_table[28, 1:5] <- round(beta_est$ps$bias, 4)
		Eff_table[29, 1:5] <- round(beta_est$ps$std, 4)
		Eff_table[30, 1:5] <- round(beta_est$ps$mse, 4)
		
		rownames(Eff_table) <- names_row
		
		#################################################
		
		# Toxicity Table
		Tox_table <- matrix(NA, nrow = 30, ncol = 5)
		
		# Model-based for Toxicity: mean, bias, SD, and MSE
		Tox_table[2, 1:5] <- round(mb_est$pt$avg, 4)
		Tox_table[3, 1:5] <- round(mb_est$pt$bias, 4)
		Tox_table[4, 1:5] <- round(mb_est$pt$std, 4)
		Tox_table[5, 1:5] <- round(mb_est$pt$mse, 4)
	
		# Dose Escalation Cohort for Toxicity: mean, bias, SD, and MSE
		Tox_table[7, 1:5] <- round(esl_est$pt$avg, 4)
		Tox_table[8, 1:5] <- round(esl_est$pt$bias, 4)
		Tox_table[9, 1:5] <- round(esl_est$pt$std, 4)
		Tox_table[10, 1:5] <- round(esl_est$pt$mse, 4)
		
		# Dose Expansion Cohort for Toxicity: mean, bias, SD, and MSE
		Tox_table[12, 1:5] <- round(exp_est$pt$avg, 4)
		Tox_table[13, 1:5] <- round(exp_est$pt$bias, 4)
		Tox_table[14, 1:5] <- round(exp_est$pt$std, 4)
		Tox_table[15, 1:5] <- round(exp_est$pt$mse, 4)
		
		# Pooled data for Toxicity: mean, bias, SD, and MSE
		Tox_table[17, 1:5] <- round(pool_est$pt$avg, 4)
		Tox_table[18, 1:5] <- round(pool_est$pt$bias, 4)
		Tox_table[19, 1:5] <- round(pool_est$pt$std, 4)
		Tox_table[20, 1:5] <- round(pool_est$pt$mse, 4)
		
		# Rao Blackwell for Toxicity: mean, bias, SD, and MSE
		Tox_table[22, 1:5] <- round(rao_est$pt$avg, 4)
		Tox_table[23, 1:5] <- round(rao_est$pt$bias, 4)
		Tox_table[24, 1:5] <- round(rao_est$pt$std, 4)
		Tox_table[25, 1:5] <- round(rao_est$pt$mse, 4)
		
		# Beta Commensurate Prior for Toxicity: mean, bias, SD, and MSE
		Tox_table[27, 1:5] <- round(beta_est$pt$avg, 4)
		Tox_table[28, 1:5] <- round(beta_est$pt$bias, 4)
		Tox_table[29, 1:5] <- round(beta_est$pt$std, 4)
		Tox_table[30, 1:5] <- round(beta_est$pt$mse, 4)
		
		rownames(Tox_table) <- names_row
		
	}
	
	if (ITH == "yes"){
		
		results_table <- NA
		
		mb_est <- SimEval(output, "mb", prob.wt, s.truth = ps, t.truth = pt, ITH)
		esl_est <- SimEval(output, "esl", prob.wt, s.truth = ps, t.truth = pt, ITH)
		exp_est <- SimEval(output, "exp", prob.wt, s.truth = ps, t.truth = pt, ITH)
		pool_est <- SimEval(output, "pool", prob.wt, s.truth = ps, t.truth = pt, ITH)
		rao_est <- SimEval(output, "rao", prob.wt, s.truth = ps, t.truth = pt, ITH)
		beta_est <- SimEval(output, "beta", prob.wt, s.truth = ps, t.truth = pt, ITH)
		
		est_name <- c("Model-Based Estimates", "Escalation Only","Expansion Only", "Pooled Estimator", "Rao-Blackwell Estimator", "Beta Commensurate Estimator")
		
		names_row <- c(sapply(est_name, function(name) c(paste(name),"Bias", "Mean Square Error")))
		
		# Efficacy Table
		Eff_table <- matrix(NA, nrow = 18, ncol = 5)
		
		# Model-based for Efficacy: mean, bias, SD, and MSE
		Eff_table[2, 1:5] <- round(mb_est$ps$bias, 4)
		Eff_table[3, 1:5] <- round(mb_est$ps$mse, 4)
	
		# Dose Escalation Cohort for Efficacy: mean, bias, SD, and MSE
		Eff_table[5, 1:5] <- round(esl_est$ps$bias, 4)
		Eff_table[6, 1:5] <- round(esl_est$ps$mse, 4)
		
		# Dose Expansion Cohort for Efficacy: mean, bias, SD, and MSE
		Eff_table[8, 1:5] <- round(exp_est$ps$bias, 4)
		Eff_table[9, 1:5] <- round(exp_est$ps$mse, 4)
		
		# Pooled data for Efficacy: mean, bias, SD, and MSE
		Eff_table[11, 1:5] <- round(pool_est$ps$bias, 4)
		Eff_table[12, 1:5] <- round(pool_est$ps$mse, 4)
		
		# Rao Blackwell for Efficacy: mean, bias, SD, and MSE
		Eff_table[14, 1:5] <- round(rao_est$ps$bias, 4)
		Eff_table[15, 1:5] <- round(rao_est$ps$mse, 4)
		
		# Beta Commensurate Prior for Efficacy: mean, bias, SD, and MSE
		Eff_table[17, 1:5] <- round(beta_est$ps$bias, 4)
		Eff_table[18, 1:5] <- round(beta_est$ps$mse, 4)
		
		rownames(Eff_table) <- names_row
		
		#################################################
		
		# Toxicity Table
		Tox_table <- matrix(NA, nrow = 18, ncol = 5)
		
		# Model-based for Toxicity: mean, bias, SD, and MSE
		Tox_table[2, 1:5] <- round(mb_est$pt$bias, 4)
		Tox_table[3, 1:5] <- round(mb_est$pt$mse, 4)
	
		# Dose Escalation Cohort for Toxicity: mean, bias, SD, and MSE
		Tox_table[5, 1:5] <- round(esl_est$pt$bias, 4)
		Tox_table[6, 1:5] <- round(esl_est$pt$mse, 4)
		
		# Dose Expansion Cohort for Toxicity: mean, bias, SD, and MSE
		Tox_table[8, 1:5] <- round(exp_est$pt$bias, 4)
		Tox_table[9, 1:5] <- round(exp_est$pt$mse, 4)
		
		# Pooled data for Toxicity: mean, bias, SD, and MSE
		Tox_table[11, 1:5] <- round(pool_est$pt$bias, 4)
		Tox_table[12, 1:5] <- round(pool_est$pt$mse, 4)
		
		# Rao Blackwell for Toxicity: mean, bias, SD, and MSE
		Tox_table[14, 1:5] <- round(rao_est$pt$bias, 4)
		Tox_table[15, 1:5] <- round(rao_est$pt$mse, 4)
		
		# Beta Commensurate Prior for Toxicity: mean, bias, SD, and MSE
		Tox_table[17, 1:5] <- round(beta_est$pt$bias, 4)
		Tox_table[18, 1:5] <- round(beta_est$pt$mse, 4)
		
		rownames(Tox_table) <- names_row
		
	}
	
	return(list(results_table = results_table, Eff_table = Eff_table, Tox_table = Tox_table))
}




###### Producing latex form of the Tables ######
library(boot)
n_sim <- 1000

#output <- cbind(sample(1:4, 100, replace =T), replicate(32,runif(100)))

# Scenario 1 Simulation Results Table
load('~/Desktop/EffToxSIM/Simulated Data/sim_results1.RData')
load('~/Desktop/EffToxSIM/results/output1.RData')
ps <- true.prob(case = 1)$ps
pt <- true.prob(case = 1)$pt
prob.wt1 <- ProbWt(sim_results1)
ps.lwr <- sapply(1:n_sim, function(i) sim_results1[[i]]$dec$ps.lwr)
pt.lwr <- sapply(1:n_sim, function(i) sim_results1[[i]]$dec$pt.lwr)
dat <- sapply(1:n_sim, function(i) sim_results1[[i]]$dat)
S1_table <- ResultsTable(output1$output1, ps, pt, prob.wt1, dat)
S1j_table <- ResultsTable(output1$output2, ps.lwr, pt.lwr, prob.wt1, dat, ITH = "yes")
write.table(S1_table$results_table, "results_tableS1.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S1_table$Eff_table, "Eff_tableS1.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S1_table$Tox_table, "Tox_tableS1.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S1j_table$Eff_table, "Eff_tableS1j.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S1j_table$Tox_table, "Tox_tableS1j.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)


# Scenario 2 Simulation Results Table
load('~/Desktop/EffToxSIM/Simulated Data/sim_results2.RData')
load('~/Desktop/EffToxSIM/results/output2.RData')
ps <- true.prob(case = 2)$ps
pt <- true.prob(case = 2)$pt
prob.wt2 <- ProbWt(sim_results2)
ps.lwr <- sapply(1:n_sim, function(i) sim_results2[[i]]$dec$ps.lwr)
pt.lwr <- sapply(1:n_sim, function(i) sim_results2[[i]]$dec$pt.lwr)
dat <- sapply(1:n_sim, function(i) sim_results2[[i]]$dat)
S2_table <- ResultsTable(output2$output1, ps, pt, prob.wt2, dat)
S2j_table <- ResultsTable(output2$output2, ps.lwr, pt.lwr, prob.wt2, dat, ITH = "yes")
write.table(S2_table$results_table, "results_tableS2.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S2_table$Eff_table, "Eff_tableS2.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S2_table$Tox_table, "Tox_tableS2.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S2j_table$Eff_table, "Eff_tableS2j.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S2j_table$Tox_table, "Tox_tableS2j.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)



# Scenario 3 Simulation Results Table
load('~/Desktop/EffToxSIM/Simulated Data/sim_results3.RData')
load('~/Desktop/EffToxSIM/results/output3.RData')
ps <- true.prob(case = 3)$ps
pt <- true.prob(case = 3)$pt
prob.wt3 <- ProbWt(sim_results3)
ps.lwr <- sapply(1:n_sim, function(i) sim_results3[[i]]$dec$ps.lwr)
pt.lwr <- sapply(1:n_sim, function(i) sim_results3[[i]]$dec$pt.lwr)
dat <- sapply(1:n_sim, function(i) sim_results3[[i]]$dat)
S3_table <- ResultsTable(output3$output1, ps, pt, prob.wt3, dat)
S3j_table <- ResultsTable(output3$output2, ps.lwr, pt.lwr, prob.wt3, dat, ITH = "yes")

write.table(S3_table$results_table, "results_tableS3.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S3_table$Eff_table, "Eff_tableS3.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S3_table$Tox_table, "Tox_tableS3.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S3j_table$Eff_table, "Eff_tableS3j.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S3j_table$Tox_table, "Tox_tableS3j.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)




# Scenario 4 Simulation Results Table
load('~/Desktop/EffToxSIM/Simulated Data/sim_results4.RData')
load('~/Desktop/EffToxSIM/results/output4.RData')
ps <- true.prob(case = 4)$ps
pt <- true.prob(case = 4)$pt
prob.wt4 <- ProbWt(sim_results4)
ps.lwr <- sapply(1:n_sim, function(i) sim_results4[[i]]$dec$ps.lwr)
pt.lwr <- sapply(1:n_sim, function(i) sim_results4[[i]]$dec$pt.lwr)
dat <- sapply(1:n_sim, function(i) sim_results4[[i]]$dat)
S4_table <- ResultsTable(output4$output1, ps, pt, prob.wt1, dat)
S4j_table <- ResultsTable(output4$output2, ps.lwr, pt.lwr, prob.wt4, dat, ITH = "yes")

write.table(S4_table$results_table, "results_tableS4.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S4_table$Eff_table, "Eff_tableS4.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S4_table$Tox_table, "Tox_tableS4.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S4j_table$Eff_table, "Eff_tableS4j.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)
write.table(S4j_table$Tox_table, "Tox_tableS4j.txt", quote = F, sep = " & ", eol = "\\\\\n", na = "", row.names = TRUE, col.names = FALSE)





