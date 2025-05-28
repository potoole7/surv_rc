#### Fit PO and Cox to skewed PWE data (w/ PO treatment effect) ####

#### libs ####

library(survival)
library(survminer)
library(timereg)
library(dplyr)
library(ggplot2)
library(parallel)
library(flexsurv)
library(pammtools)
library(patchwork)
source("00_source.R")

#### Functions ####

# TODO Change name, very verbose (also do with var names)
rpwexp_inv <- function(
    haz0, times, beta = NULL) {
  # censoring_time = NULL, censoring_rate = NULL) {
  # Validate inputs
  m <- length(times)
  if (length(haz0) != m) {
    stop("haz0 and times must have the same length")
  }

  # Generate event time (t_event)
  u <- runif(1)
  if (!is.null(beta) && beta != 0) {
    u0 <- u / (u + exp(beta) * (1 - u))
  } else {
    u0 <- u
  }
  eta <- -log(u0)

  # Precompute cumulative hazards
  H <- numeric(m + 1)
  H[1] <- 0.0
  for (i in 1:m) {
    prev_time <- if (i == 1) 0.0 else times[i - 1]
    delta <- times[i] - prev_time
    H[i + 1] <- H[i] + haz0[i] * delta
  }

  # Find interval and compute t_event
  j <- findInterval(eta, H, rightmost.closed = FALSE)
  if (j <= m) {
    prev_time <- if (j == 1) 0.0 else times[j - 1]
    t_event <- prev_time + (eta - H[j]) / haz0[j]
  } else {
    t_event <- times[m] + (eta - H[m + 1]) / haz0[m]
  }

  # Generate censoring time (t_censor)
  # if (!is.null(censoring_time)) {
  #   t_censor <- censoring_time # Fixed censoring time
  # } else if (!is.null(censoring_rate)) {
  #   t_censor <- rexp(1, rate = censoring_rate) # Random exponential censoring
  # } else {
  #   t_censor <- Inf # No censoring
  # }
  #
  # # Compute observed time and status
  # observed_time <- min(t_event, t_censor)
  # status <- as.integer(t_event <= t_censor)
  observed_time <- t_event

  # if (n == 1) {
  #   return(data.frame(time = observed_time))
  #   # return(list(time = observed_time, status = status))
  # } else {
  #   return(rbind(
  #     data.frame(time = observed_time),
  #     rpwexp_inv(haz0, times, beta, n - 1)
  #   ))
  # }
  return(observed_time)
}

# generate dataset with
rpwexp_q_po <- \(n, times, quantiles, beta, x) {
  # only binomial x supported for now
  stopifnot(all(x %in% c(0, 1)))
  stopifnot(length(x) == n) # x must be same length as n!

  haz0 <- haz(times, quantiles)
  n0 <- sum(x == 0)
  df_baseline <- data.frame(
    time = replicate(n0, rpwexp_inv(haz0, times))
  )

  # Treatment group (beta = -0.5) with fixed censoring
  # df_treatment <- rpwexp_inv(
  #   haz0, times,
  #   beta = beta1, n = 200
  # )
  n1 <- sum(x == 1)
  df_treatment <- data.frame(
    time = replicate(n1, rpwexp_inv(haz0, times, beta = beta1))
  )

  # Combine baseline and treatment data
  # TODO Need to combine and order like in `x`!
  data_po <- rbind(
    cbind(df_baseline, group = 0),
    cbind(df_treatment, group = 1)
  )
  return(data_po)
}

#### metadata ####

# set.seed(123)
times <- c(1, 1.1, 2.5)
quantiles <- c(0.25, 0.5, 0.75)
beta1 <- 0.5
cores <- parallel::detectCores() - 1

#### Plot example ####

data_po_plt <- rpwexp_q_po(
  n = 10000,
  times = times,
  quantiles = quantiles,
  beta = beta1,
  x = rbinom(10000, 1, 0.5) # Random treatment assignment
)

data_po_plt <- data_po_plt |>
  mutate(event = rbinom(n(), 1, 0.8))

# Fit Kaplan-Meier curves
fit_po <- survfit(Surv(time, event) ~ group, data = data_po_plt)

# Plot
p2 <- ggsurvplot(
  fit_po,
  data = data_po_plt,
  risk.table = FALSE,
  palette = ggsci::pal_nejm()(4)[3:4],
  theme = evc::evc_theme(),
  legend.labs = c("Control", "Treatment"),
  lwd = 1,
  alpha = 0.6,
  legend.title = "Group",
  xlim = c(0, 11),
  # add custom x axis breaks
  break.time.by = 1
)
p2

ggsave(p2$plot, filename = "report/plots/02_sim_weird_dist_data_po.png", width = 8, height = 5)


#### Example fit ####

set.seed(123)
data <- rpwexp_q_po(
  n = 200,
  times = times,
  quantiles = c(0.25, 0.5, 0.75),
  beta = beta1,
  x = rbinom(200, 1, 0.5) # Random treatment assignment
) |>
  mutate(
    event = rbinom(n(), 1, 0.8) # censoring indicator
  )

coef(coxph(Surv(time, event) ~ group, data = data))

# try fit PO model
out <- timereg::prop.odds(
  Event(time, event) ~ group,
  data = data
)
out$gamma # close enough!
out$var.gamma[[1]]

# repeat 500 times
n_sims <- 500
df_lst <- lapply(seq_len(n_sims), \(i) {
  print(paste("Simulation", i, "of", n_sims))
  rpwexp_q_po(
    n = 200,
    times = times,
    quantiles = quantiles,
    beta = beta1,
    x = rbinom(200, 1, 0.5) # Random treatment assignment
  ) |>
    mutate(
      event = rbinom(n(), 1, 0.8) # censoring indicator
    )
})

# fit PO model to each simulated dataset
# out_lst <- mclapply(df_lst, \(df) {
out_lst <- mclapply(seq_along(df_lst), \(i) {
  system(sprintf("echo %s", paste("Simulation", i, "of", n_sims)))
  timereg::prop.odds(
    Event(time, event) ~ group,
    # data = df
    data = df_lst[[i]]
  )
}, mc.cores = cores)

coefs_po <- sapply(out_lst, function(x) -x$gamma)

# plot, get mean, bias and variance
hist_fun(coefs_po)
mean(coefs_po)
formatC(mean(coefs_po - beta1), format = "e")
# variance
formatC(var(coefs_po), format = "e")

# try fit Cox
coefs_coxph <- unlist(mclapply(df_lst, \(x) {
  # coef(coxph(Surv(time) ~ group, data = x))
  unname(coef(coxph(Surv(time, event) ~ group, data = x)))
}, mc.cores = cores))

hist_fun(coefs_coxph, beta = NULL)
mean(coefs_coxph)
formatC(mean(coefs_coxph - beta1), format = "e")
# variance
formatC(var(coefs_coxph), format = "e") # TODO Compare to theory

#### Grid search ####

# grid of values
# time_vals <- c("1, 1.1", "1, 1.1, 2.5", "1, 1.1, 7", "1, 1.1, 11")
# n_sims <- 500
n_sims <- 1000
# time_vals <- c("1, 1.1, 1.5", "1, 1.1, 2.5", "1, 1.1, 7", "1, 1.1, 11")
time_vals <- c(
  "1, 1.1, 2.5",
  "1, 1.5, 2.5",
  "1, 1.05, 2.5",
  "1, 1.01, 2.5"
)
n_vals <- c(100, 200, 500)
# censor_prop <- c(0.2, 0.5, 0.8)
censor_prop <- 0.2
# treat_prop <- c(0.2, 0.5, 0.8)
treat_prop <- 0.5

grid <- tidyr::crossing(
  "time"   = time_vals,
  "n"      = n_vals,
  "censor" = censor_prop,
  "treat"  = treat_prop
) |>
  tidyr::separate(
    time,
    into = c("time1", "time2", "time3"), sep = ",\\s*", remove = FALSE
  ) %>%
  mutate(across(c("time1", "time2", "time3"), as.numeric))

set.seed(123)
# run simulation for each row in grid
sim_grid_lst <- lapply(seq_len(nrow(grid)), \(i) {
  x <- grid[i, ]
  times <- c(x$time1, x$time2, x$time3)
  n <- x$n

  mclapply(seq_len(n_sims), \(i) {
    # simulate data
    rpwexp_q_ph(
      n = n,
      times = times, # Quantile times (25th, 50th, 75th)
      quantiles = quantiles,
      beta = beta1,
      x = rbinom(n, 1, x$treat) # Random treatment assignment
    ) |>
      mutate(
        event = rbinom(n(), 1, x$censor), # censoring indicator
      )
  })
})

# test
prop.odds(Event(time, event) ~ group, data = sim_grid_lst[[1]][[1]])$gamma[[1]]

# fit PO
# po_res <- mclapply(sim_grid_lst, \(x) {
#   # fit Cox to each, extract coefficient
#   sapply(x, \(y) prop.odds(Event(time, event) ~ group, data = y)$gamma[[1]])
# }, mc.cores = cores)
po_res <- lapply(seq_along(sim_grid_lst), \(i) {
  print(sprintf("Row %s of %s", i, nrow(grid)))
  # fit Cox to each, extract coefficient
  mclapply(seq_along(sim_grid_lst[[i]]), \(j) {
    system(sprintf("echo %s", paste("Simulation", j, "of", length(sim_grid_lst[[i]]))))
    prop.odds(Event(time, event) ~ group, data = sim_grid_lst[[i]][[j]])$gamma[[1]]
  }, mc.cores = cores)
})

# join to grid
grid_res_po_df <- bind_rows(lapply(seq_along(po_res), \(i) {
  data.frame(grid[i, ], coefs = unlist(po_res[[i]]))
})) |>
  filter(abs(coefs) < 5) |>
  group_by(time, n, censor, treat) |>
  mutate(
    mean_coefs = mean(coefs),
    bias = abs(mean(coefs - beta1)),
    var = var(coefs)
  )

pgrid_po <- grid_res_po_df |>
  ggplot(aes(x = coefs)) +
  geom_histogram(colour = "black", fill = "grey") +
  labs(
    x = latex2exp::TeX("$\\beta$"),
    y = "Count"
  ) +
  # geom_vline(xintercept = mean(grid_res_ph_df$coefs), col = "blue", lty = 1, lwd = 1.2) +
  geom_vline(aes(xintercept = mean_coefs), col = "blue", lty = 1, lwd = 1.2) +
  geom_vline(xintercept = beta1, col = "red", lty = 2, lwd = 1.2) +
  # add bias and variance as text
  geom_text(
    aes(x = 1.8, y = 300, label = paste0("Bias: ", round(bias, 3)))
  ) +
  geom_text(
    # aes(x = 2, y = 250, label = paste0("s: ", round(var, 3))),
    # label with s^2
    aes(x = 1.8, y = 250, label = paste0(latex2exp::TeX("Var: "), round(var, 3)))
  ) +
  # add bias as annotation
  # facet_wrap(n ~ time, scales = "free") +
  facet_grid(n ~ time) +
  evc::evc_theme()
pgrid_po

ggsave(plot = pgrid_po, filename = "report/plots/03_sim_res_po_grid.png", width = 10, height = 8)

grid_res_po_df |>
  distinct(time, n, mean_coefs, bias, var)

# fit Cox
cox_res <- mclapply(sim_grid_lst, \(x) {
  # fit Cox to each, extract coefficient
  sapply(x, \(y) unname(coef(coxph(Surv(time, event) ~ group, data = y))))
}, mc.cores = cores)

# join to grid
grid_res_ph_df <- bind_rows(lapply(seq_along(cox_res), \(i) {
  data.frame(grid[i, ], coefs = cox_res[[i]])
})) |>
  filter(abs(coefs) < 5) |>
  group_by(time, n, censor, treat) |>
  mutate(mean_coefs = mean(coefs))

pgrid_ph <- grid_res_ph_df |>
  ggplot(aes(x = coefs)) +
  geom_histogram(colour = "black", fill = "grey") +
  labs(
    x = latex2exp::TeX("$\\beta$"),
    y = "Count"
  ) +
  # geom_vline(xintercept = mean(grid_res_ph_df$coefs), col = "blue", lty = 1, lwd = 1.2) +
  geom_vline(aes(xintercept = mean_coefs), col = "blue", lty = 1, lwd = 1.2) +
  # geom_vline(xintercept = beta1, col = "red", lty = 2, lwd = 1.2) +
  geom_text(
    aes(x = 1.8, y = 250, label = paste0(latex2exp::TeX("Mean: "), round(mean_coefs, 3)))
  ) +
  facet_grid(n ~ time) +
  evc::evc_theme()
pgrid_ph

ggsave(plot = pgrid_ph, filename = "report/plots/03_sim_res_ph_grid_po_dat.png", width = 10, height = 8)

grid_res_ph_df |>
  distinct(time, n, mean_coefs)

# Investigate how well this agrees with theory

# fit model
set.seed(123)
fit <- coxph(Surv(time, event) ~ group, data = sim_grid_lst[[1]][[1]])
sandwich_mat <- vcov(fit) %*% 
  matrix(norm(coxph.detail(fit)$score, type = "2")) %*% 
  vcov(fit)
# calculate confidence interval
coef(fit)["group"]
coef(fit)["group"] + c(-1, 1) * qnorm(0.975) * sqrt(sandwich_mat[[1]])
set.seed()123


