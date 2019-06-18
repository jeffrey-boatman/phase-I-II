library(Hmisc) # For binconf (i.e. wilson binomial confidence intervals)
library(rjags)
library(parallel) 
# library(boot)

source("functions.R")

###number of simulations###

n_sim   <- 1000
n_cores <- 12

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

###translate beta1_muT and beta1_sigT into gamma parameters###

beta1_shapeT <- (beta1t_m^2)/(beta1t_s^2)
beta1_rateT <- (beta1t_m)/(beta1t_s^2)

beta1_shapeE <- (beta1s_m^2)/(beta1s_s^2)
beta1_rateE <- (beta1s_m)/(beta1s_s^2)



###translating hyperparameters for SD to precision###

beta0t_p <- 1/(beta0t_s^2)
# beta1t_p <- 1/(beta1t_s^2)

beta0s_p <- 1/(beta0s_s^2)
# beta1s_p <- 1/(beta1s_s^2)
beta2s_p <- 1/(beta2s_s^2)

###calculate p using contour propossed by Thall and Cook###

find_p <- function(p) {
  (((1 - ps_star) / (1 - ps_min)) ^ p + (pt_star / pt_max) ^ p - 1) ^ 2 
}

p <- optimize(find_p, c(-100,100))$minimum

load("./RData/data_sets.RData")

# debug(apply_estimators)
# debug(estimate)
# debug(Pool)
# debug(RaoBlackwell)
# debug(BetaComP)
# debug(calc_MEM)
# apply_estimators(data_sets[[1]], include.ci = FALSE)

host       <- system2("hostname", stdout = TRUE)
servers    <- c("carbon", "cesium", "chromium", "potassium", "silicon")
hosts      <- paste0(servers, ".ccbr.umn.edu")
n_serv     <- length(hosts)
which_serv <- match(host, hosts)

RNGkind("L'Ecuyer-CMRG")
set.seed(1)

# out1 <- mclapply(data_sets[1:2],
#   FUN      = quiet(apply_estimators),
#   mc.cores = 2)
# out2 <- mclapply(data_sets[1:2],
#   FUN      = quiet(apply_estimators),
#   mc.cores = 2)

# rbind(out1[[1]], out2[[1]])


# args <- expand.grid(
#   seed     = seq_len(n_sim),
#   scenario = 1:4
# )

n_sim <- length(data_sets)

# which_serv_args <- rep(seq_len(n_serv), nrow(args) / n_serv)
# args <- args[which_serv_args == which_serv, ]

which_sims <- which(rep(1:n_serv, n_sim / n_serv) == which_serv)

# seed     <- args$seed
# scenario <- args$scenario

# sim_out <- mcmapply(FUN = failwith(NA, quiet(sim_and_estimate)),
#   seed     = seed,
#   scenario = scenario,
#   mc.cores = n_cores,
#   SIMPLIFY = FALSE)

# apply_estimators(data_sets[[47]])

sim_out <- mclapply(data_sets[which_sims],
  FUN        = failwith(NA, quiet(apply_estimators)),
  include.ci = FALSE,
  mc.cores   = n_cores)


sim_out <- do.call(rbind, sim_out)

write.table(sim_out,
  file      = paste0("./output/out_", which_serv, ".txt"),
  quote     = FALSE,
  row.names = FALSE)
