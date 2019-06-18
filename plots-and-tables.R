options(digits = 3)
library(tidyverse)
library(xtable)

outfiles <- paste0("./output/out_", 1:5, ".txt")

rrl <- lapply(outfiles, read.table, header = TRUE)

rr <- as_tibble(do.call(rbind, rrl))

setdiff(1:1000, rr$data_set_num)

# remove the model-based and escalation-only estimators
rr <- rr %>% select(-contains("mb"), -contains("esl"))

# format the data-generation names
old_names <- c("hiTox_loEff", "ITH", "no_effect")
new_names <- c("HTLE", "ITH", "NE")
rr <- rr %>% 
  mutate(data_gen_type = new_names[match(data_gen_type, old_names)])

# selection probabilities
sel_prob <- as_tibble(prop.table(table(rr$scenario, rr$optimal_dose), 1)) %>%
  rename(scenario = Var1, dose = Var2, prob = n) %>%
  arrange(scenario, dose) %>%
  mutate(
    scenario = as.integer(scenario),
    dose = as.integer(dose))

# bias and mse, toxicity ----
tox_by_dose <- rr %>%
  select(data_set_num, data_gen_type, scenario, optimal_dose, prob_tox, contains("tox_est")) %>%
  gather(contains("tox_est"), key = estimator, value = estimate) %>%
  separate(estimator, into = c("estimator", "junk"), "_tox_") %>%
  select(-junk) %>%
  mutate(error = estimate - prob_tox) %>%
  rename(dose = optimal_dose) %>%
  group_by(estimator, data_gen_type, scenario, dose) %>%
  summarize(
    SE   = sd(estimate   , na.rm = TRUE),
    mean = mean(estimate , na.rm = TRUE),
    bias = mean(error    , na.rm = TRUE),
    mse  = mean(error ^ 2, na.rm = TRUE)
  )

tox <- tox_by_dose %>%
  left_join(sel_prob) %>%
  group_by(estimator, data_gen_type, scenario) %>%
  summarize(
    bias = weighted.mean(bias, prob),
    mse  = weighted.mean(mse , prob)
  )

# format estimator and names
old_names <- c("cp", "dec", "mem", "pool", "rb")
new_names <- c("CP", "DEC", "MEMs", "Pooled", "UMVUE")

tox <- tox %>% 
  ungroup() %>%
  mutate(
    estimator = new_names[match(estimator, old_names)],
    scenario  = paste0("Scenario ", scenario)
  )

# tox$estimator <- factor(tox$estimator, 
#   levels  = c("DEC", "Pooled", "UMVUE", "CP", "MEMs"),
#   ordered = TRUE)
# tox$estimator <- as.character(tox$estimator)

ord <- c("DEC", "Pooled", "UMVUE", "CP", "MEMs")

font_size <- 12.5
types <- c("NE", "ITH", "HTLE")
types_title <- c(
  "No Inter-Trial Heterogeneity",
  "Random Inter-Trial Heterogeneity",
  "Random Lower Efficacy, Higher Toxicity"
)

# ~ toxicity plots ~ ----
for(type in types) {
  pdf(paste0("./plots/tox-", type, ".pdf"), height = 8, width = 12)
  temp_title <- paste0("Toxicity, ", types_title[match(type, types)])
  p1 <- tox %>% filter(data_gen_type == type) %>%
    ggplot(mapping = aes(x = estimator, y  = abs(bias), fill = estimator)) +  
    geom_bar(stat = "identity") +
    ylab("Absolute Bias") +
    xlab("Estimator") +
    labs(fill = "Estimator", title = temp_title) +
    facet_wrap(~ scenario, nrow = 1) +
    scale_x_discrete(limits = ord) +
    scale_fill_discrete(limits = ord) +
    theme_classic() +
    theme(text =element_text(size = font_size))
  p2 <- tox %>% filter(data_gen_type == type) %>%
    ggplot(mapping = aes(x = estimator, y  = mse, fill = estimator)) +  
    geom_bar(stat = "identity") +
    ylab("MSE") +
    xlab("Estimator") +
    labs(fill = "Estimator") +
    facet_wrap(~ scenario, nrow = 1) +
    scale_x_discrete(limits = ord) +
    scale_fill_discrete(limits = ord) +
    theme_classic() +
    theme(text =element_text(size = font_size))
  gridExtra::grid.arrange(p1, p2, nrow = 2)
  dev.off()
}

# bias and mse, efficacy ----
eff_by_dose <- rr %>%
  select(data_set_num, data_gen_type, scenario, optimal_dose, prob_eff, contains("eff_est")) %>%
  gather(contains("eff_est"), key = estimator, value = estimate) %>%
  separate(estimator, into = c("estimator", "junk"), "_eff_") %>%
  select(-junk) %>%
  mutate(error = estimate - prob_eff) %>%
  rename(dose = optimal_dose) %>%
  group_by(estimator, data_gen_type, scenario, dose) %>%
  summarize(
    SE   = sd(estimate   , na.rm = TRUE),
    mean = mean(estimate , na.rm = TRUE),
    bias = mean(error    , na.rm = TRUE),
    mse  = mean(error ^ 2, na.rm = TRUE)
  )

eff <- eff_by_dose %>%
  left_join(sel_prob) %>%
  group_by(estimator, data_gen_type, scenario) %>%
  summarize(
    bias = weighted.mean(bias, prob),
    mse  = weighted.mean(mse , prob)
  )

# format estimator and names
old_names <- c("cp", "dec", "mem", "pool", "rb")
new_names <- c("CP", "DEC", "MEMs", "Pooled", "UMVUE")

eff <- eff %>% 
  ungroup() %>%
  mutate(
    estimator = new_names[match(estimator, old_names)],
    scenario  = paste0("Scenario ", scenario)
  )

# eff$estimator <- factor(eff$estimator, 
#   levels  = c("DEC", "Pooled", "UMVUE", "CP", "MEMs"),
#   ordered = TRUE)

for(type in c("NE", "ITH", "HTLE")) {
  pdf(paste0("./plots/eff-", type, ".pdf"), height = 8, width = 12)
  temp_title <- paste0("Efficacy, ", types_title[match(type, types)])
  p1 <- eff %>% filter(data_gen_type == type) %>%
    ggplot(mapping = aes(x = estimator, y  = abs(bias), fill = estimator)) +  
    geom_bar(stat = "identity") +
    ylab("Absolute Bias") +
    xlab("Estimator") +
    labs(fill = "Estimator", title = temp_title) +
    facet_wrap(~ scenario, nrow = 1) +
    scale_x_discrete(limits = ord) +
    scale_fill_discrete(limits = ord) +
    theme_classic() +
    theme(text =element_text(size = font_size))
  p2 <- eff %>% filter(data_gen_type == type) %>%
    ggplot(mapping = aes(x = estimator, y  = mse, fill = estimator)) +  
    geom_bar(stat = "identity") +
    ylab("MSE") +
    xlab("Estimator") +
    labs(fill = "Estimator") +
    facet_wrap(~ scenario, nrow = 1) +
    scale_x_discrete(limits = ord) +
    scale_fill_discrete(limits = ord) +
    theme_classic() +
    theme(text =element_text(size = font_size))
  gridExtra::grid.arrange(p1, p2, nrow = 2)
  dev.off()
}

# ~ latex tables ~ ----

# ~~ selection probabilities ~~ ----
# <- 
sel_prob_table <- prop.table(table(rr$scenario, rr$optimal_dose), 1)
# sel_prob_table <- round(sel_prob_table, 3)
sel_prob_table <- apply(sel_prob_table, c(1, 2), function(x) sprintf("%.3f", x))
which_row <- 1:4
which_col <- c(3, 1, 4, 2)
bf <- function(x) paste0("\\bf{", x, "}")
for(i in 1:4) {
  sel_prob_table[which_row[i], which_col[i]] <- bf(sel_prob_table[which_row[i], which_col[i]])
}
rownames(sel_prob_table) <- paste0("Scenario ", 1:4)
colnames(sel_prob_table) <- paste0("Dose", 1:4) 

print(xtable(sel_prob_table),
  file                   = "./tables/sim-selection-probs.tex",
  only.contents          = TRUE,
  sanitize.text.function = identity,
  include.rownames       = TRUE)



old_names <- c("cp", "dec", "mem", "pool", "rb")
new_names <- c("CP", "DEC", "MEMs", "Pooled", "UMVUE")

# ~~ efficacy ~~ ----
eff_by_dose <- eff_by_dose %>% 
  ungroup() %>%
    mutate(
    estimator = new_names[match(estimator, old_names)]
    )

eff_by_dose$estimator <- factor(eff_by_dose$estimator, 
  levels  = c("DEC", "Pooled", "UMVUE", "CP", "MEMs"),
  ordered = TRUE)


estimators <- levels(eff_by_dose$estimator)
data_gens <- na.omit(unique(eff_by_dose$data_gen_type)) # !!!! ---- na.omit?
for(cur_data_gen in data_gens) {
  for(cur_scenario in 1:4) {
    out <- NULL
    for(cur_estimator in estimators) {
      #cur_scenario  <- 1
      # cur_data_gen  <- "NE"
      # cur_estimator <- "CP"
      cur_sel_prob  <- sel_prob %>% 
        filter(scenario == cur_scenario) %>%
        arrange(dose) %>%
        select(prob)
      mcur_sel_prob <- t(matrix(cur_sel_prob[,,drop = TRUE], 4, 4))

      eff_cur <- eff_by_dose %>%
        filter(scenario == cur_scenario, data_gen_type == cur_data_gen, estimator == cur_estimator) %>%
        arrange(dose) %>%
        select(mean, bias, SE, mse)
      teff_cur <- t(eff_cur)
      wa <- rowSums(teff_cur * mcur_sel_prob, na.rm = TRUE)
      teff_cur <- cbind(teff_cur, wa)
      teff_cur <- rbind(NA, teff_cur)
      teff_cur <- as.data.frame(teff_cur)

      rnames <- c(cur_estimator, paste0("\\qquad ", c("Mean", "Bias", "SE", "MSE")))

      teff_cur <- cbind(rnames, teff_cur)

      colnames(teff_cur) <- c("Estimator", paste0("Dose ", 1:4), "Weighted Avg.")
      # rownames(teff_cur) <- c(cur_estimator, paste0("\\qquad ", c("Mean", "Bias", "SE", "MSE")))

      out <- rbind(out, teff_cur)
    }

    outfile <- paste0("./tables/eff-", cur_data_gen, "-scenario-", cur_scenario, ".tex")
    hlines <- c(0, 0:(length(estimators) - 1) * 5 + 1, 1:length(estimators) * 5)
    print(xtable(out, digits = c(0, 0, 3, 3, 3, 3, 3)),
      file                   = outfile,
      only.contents          = TRUE,
      sanitize.text.function = identity,
      include.rownames       = FALSE,
      hline.after            = hlines,
      NA.string              = "-")
    cat("Generating Table for Efficacy, trial type: ", cur_data_gen, ", scenario: ", cur_scenario, "\n")
  }
}

# ~~ toxicity ~~ ----
tox_by_dose <- tox_by_dose %>% 
  ungroup() %>%
    mutate(
    estimator = new_names[match(estimator, old_names)]
    )

tox_by_dose$estimator <- factor(tox_by_dose$estimator, 
  levels  = c("DEC", "Pooled", "UMVUE", "CP", "MEMs"),
  ordered = TRUE)


estimators <- levels(tox_by_dose$estimator)
data_gens <- na.omit(unique(tox_by_dose$data_gen_type)) # !!!! ---- na.omit?
for(cur_data_gen in data_gens) {
  for(cur_scenario in 1:4) {
    out <- NULL
    for(cur_estimator in estimators) {
      #cur_scenario  <- 1
      # cur_data_gen  <- "NE"
      # cur_estimator <- "CP"
      cur_sel_prob  <- sel_prob %>% 
        filter(scenario == cur_scenario) %>%
        arrange(dose) %>%
        select(prob)
      mcur_sel_prob <- t(matrix(cur_sel_prob[,,drop = TRUE], 4, 4))

      tox_cur <- tox_by_dose %>%
        filter(scenario == cur_scenario, data_gen_type == cur_data_gen, estimator == cur_estimator) %>%
        arrange(dose) %>%
        select(mean, bias, SE, mse)
      ttox_cur <- t(tox_cur)
      wa <- rowSums(ttox_cur * mcur_sel_prob)
      ttox_cur <- cbind(ttox_cur, wa)
      ttox_cur <- rbind(NA, ttox_cur)
      ttox_cur <- as.data.frame(ttox_cur)

      rnames <- c(cur_estimator, paste0("\\qquad ", c("Mean", "Bias", "SE", "MSE")))

      ttox_cur <- cbind(rnames, ttox_cur)

      colnames(ttox_cur) <- c("Estimator", paste0("Dose ", 1:4), "Weighted Avg.")
      # rownames(ttox_cur) <- c(cur_estimator, paste0("\\qquad ", c("Mean", "Bias", "SE", "MSE")))

      out <- rbind(out, ttox_cur)
    }

    outfile <- paste0("./tables/tox-", cur_data_gen, "-scenario-", cur_scenario, ".tex")
    hlines <- c(0, 0:(length(estimators) - 1) * 5 + 1, 1:length(estimators) * 5)
    print(xtable(out, digits = c(0, 0, 3, 3, 3, 3, 3)),
      file                   = outfile,
      only.contents          = TRUE,
      sanitize.text.function = identity,
      include.rownames       = FALSE,
      hline.after            = hlines)
    cat("Generating Table for Toxicity, trial type: ", cur_data_gen, ", scenario: ", cur_scenario, "\n")
  }
}
