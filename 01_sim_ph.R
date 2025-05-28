#### Fit Cox and PEM to skewed PWE data ####

# Add simulations for PEM model (PH and PO) (done)

# For report:
# Fix PO model (or data simulation), not working (done)
# Plot of "weird" simulation data
# TODO Simulation results for Cox and PEM under PH assumption
# TODO Vary times for quantiles, as well as n
# TODO Simulation results for Cox, PO and PEM under PH assumption
# TODO Vary t and n
# TODO Estimate theoretical covariance matrix using what Karim asked for!

#### Libs ####

library(dplyr)
library(ggplot2)
library(pammtools)
library(flexsurv)
library(parallel)
library(splines)
library(nphsim)
library(timereg)
source("code/source.R")

#### Metadata ####

# times at which we observe S(t1) = 0.25, S(t2) = 0.5, S(t3) = 0.75
# t1 and t2 are close so that we have a weird shape!
times <- c(1, 1.1, 2)
quantiles <- c(0.25, 0.5, 0.75) # quantiles of S(t)
beta1 <- 0.5 # treatment effect

cores <- parallel::detectCores() - 1

# plot colours
colours <- ggsci::pal_nejm()(3)

#### Plotting simulated data ####

# plot KM curve for different time points
plot_times <- list(
  c(1, 1.1, 2.5), # original
  c(1, 1.5, 2.5), # spread out
  # c(1, 1.1, 5), # spread last time (doesn't seem to do anything??)
  c(1, 1.05, 2.5), # squish first two
  c(1, 1.01, 2.5) # squish even more
)

set.seed(123)
plot_dat <- bind_rows(lapply(plot_times, \(x) {
  haz0 <- haz(x, quantiles)
  # Simulate data for each time point
  # rpwexp_q_ph(
  #   n = 200,
  #   times = x,
  #   quantiles = c(0.25, 0.5, 0.75),
  #   beta = beta1,
  #   x = rbinom(200, 1, 0.5) # Random treatment assignment
  # ) |>
  #   mutate(
  #     event = rbinom(n(), 1, 0.8), # censoring indicator
  #     times = paste(x, collapse = ", ")
  #   )
  data.frame(
    "time"  = nphsim::rpwexp(n = 10000, rate = haz0, intervals = x),
    "times" = paste(x, collapse = ", ")
  )
}))

fit <- survfit(
  Surv(time) ~ times,
  data = plot_dat
)

p0 <- ggsurvplot(
  fit,
  data = plot_dat,
  risk.table = FALSE,
  palette = ggsci::pal_nejm()(length(plot_times)), # use NEJM palette for colours,
  theme = evc::evc_theme(),
  legend.labs = c(
    c("1, 1.01, 2.5"),
    c("1, 1.05, 2.5"),
    c("1, 1.1, 2.5"),
    # c("1, 1.1, 5"),
    c("1, 1.5, 2.5")
  ),
  lwd = 1.5,
  alpha = 0.7,
  legend.title = "Time points",
  xlim = c(0, 6),
  # add custom x axis breaks
  break.time.by = 2
)
p0

ggsave(plot = p0$plot, "report/plots/00_sim_weird_dist_data.png", width = 8, height = 5)

# Simulate PH data and plot
set.seed(123)
data_ph <- rpwexp_q_ph(
  n = 100000,
  times = times,
  quantiles = c(0.25, 0.5, 0.75),
  beta = beta1,
  x = rbinom(100000, 1, 0.5) # Random treatment assignment
) |>
  # mutate(
  #   event = rbinom(n(), 1, 0.8) # censoring indicator
  # ) |>
  identity()

fit_ph <- survfit(
  Surv(time) ~ group,
  data = data_ph
)

p1 <- ggsurvplot(
  fit_ph,
  data = data_ph,
  risk.table = FALSE,
  palette = ggsci::pal_nejm()(2),
  theme = evc::evc_theme(),
  legend.labs = c("Treatment", "Control"),
  lwd = 1.5,
  alpha = 0.6,
  legend.title = "Group",
  xlim = c(0, 3),
  # add custom x axis breaks
  break.time.by = 1
)

saveRDS(p1, file = "report/plots/01_sim_weird_dist_data_ph.rds")
ggsave(p1$plot, filename = "report/plots/01_sim_weird_dist_data_ph.png", width = 8, height = 5)


#### Example (PH) ####

# Simulate data with beta1 = 0.5
set.seed(123)
data <- rpwexp_q_ph(
  n = 200,
  times = times,
  quantiles = c(0.25, 0.5, 0.75),
  beta = beta1,
  x = rbinom(200, 1, 0.5) # Random treatment assignment
) |>
  mutate(
    event = rbinom(n(), 1, 0.8) # censoring indicator
  )

# plot data
# p1 <- data |>
#   # mutate(group = factor(group)) |>
#   mutate(group = factor(group, labels = c("Control", "Treatment"))) |>
#   ggplot(aes(x = time, colour = group)) +
#   geom_density(alpha = 0.9) +
#   geom_vline(
#     xintercept = times,
#     linetype   = "dashed",
#     colour     = "grey",
#     size       = 1
#   ) +
#   labs(
#     x = "Time",
#     y = "Density",
#     colour = "Group"
#   ) +
#   ggsci::scale_colour_nejm() +
#   # in legend, add fill to boxes
#   guides(colour = guide_legend(override.aes = list(fill = c("Control" = colours[1], "Treatment" = colours[2])))) +
#   scale_x_continuous(breaks = seq(0, 3, by = 0.5)) +
#   evc::evc_theme()
# p1
# ggsave(plot = p1, "report/plots/01_sim_weird_dist_data_ph.png", width = 8, height = 5)

# not bad!
coef(coxph(Surv(time, event) ~ group, data = data))

# try Poisson GLM as test
data_ped <- data |>
  as_ped(Surv(time, event) ~ group, cut = times)
pam <- mgcv::gam(
  ped_status ~ interval + group,
  data = data_ped,
  family = poisson(),
  offset = offset
)
coef(pam)["group"]
# Works as well!

# see how PEM does
pem_fit <- fit_pem(data, times)
data_ped <- pem_fit$data
coef(pem_fit$fit)["group"]
# excellent also! Little bit slower, how would it do with more/less times?

#### Simulation study ####

# generate simulated datasets
n_sims <- 500
set.seed(123)
sim_data_lst <- mclapply(seq_len(n_sims), \(i){
  rpwexp_q_ph(
    n = 200,
    times = times, # Quantile times (25th, 50th, 75th)
    quantiles = quantiles,
    beta = beta1,
    x = rbinom(200, 1, 0.5) # Random treatment assignment
  ) |>
    mutate(
      event = rbinom(n(), 1, 0.8), # censoring indicator
    )
}, mc.cores = cores)

# fit Cox to each, extract coefficient
coefs_coxph <- unlist(mclapply(sim_data_lst, \(x) {
  # coef(coxph(Surv(time) ~ group, data = x))
  unname(coef(coxph(Surv(time, event) ~ group, data = x)))
}, mc.cores = cores))

# plot histogram with true value
hist(coefs_coxph)
abline(v = beta1, col = "red", lty = 2, lwd = 2)

p2 <- hist_fun(coefs_coxph)

ggsave(plot = p2, "report/plots/02_sim_res_ph_cox.png", width = 8, height = 5)

# bias
mean(coefs_coxph)
formatC(mean(coefs_coxph - beta1), format = "e")
# variance
formatC(var(coefs_coxph), format = "e")

# now try for PEM
coefs_pem <- unlist(mclapply(sim_data_lst, \(x) {
  # coefs_pem <- unlist(lapply(seq_along(sim_data_lst), \(i) {
  unname(coef(fit_pem(x, times))["group"])
}, mc.cores = cores))

p3 <- hist_fun(coefs_pem)

ggsave(plot = p3, filename = "report/plots/02a_sim_res_ph_pem.png", width = 8, height = 5)

hist(coefs_pem)
abline(v = beta1, col = "red", lty = 2, lwd = 2)
mean(coefs_pem)
formatC(mean(coefs_pem - beta1), format = "e") # worse by a factor of 10 (0.0001 vs 0.007)
formatC(var(coefs_pem), format = "e") # similar

# ## Further test examples ##
#
# # TODO Spread out times, see if bias reduces
# # TODO If not, change binomial split (50/50, 80/20, etc)
# # TODO See changes due to difference in n, t
#
# # spread out t
# # times <- c(1, 1.1, 2.5) # original
# times2 <- c(1, 1.1, 12)
# set.seed(123)
# data <- rpwexp_q_ph(
#   n = 1000,
#   times = times2,
#   quantiles = c(0.25, 0.5, 0.75),
#   beta = beta1,
#   x = rbinom(1000, 1, 0.5) # Random treatment assignment
# ) |>
#   mutate(
#     event = rbinom(n(), 1, 0.8)
#   )
#
# # plot data
# data |>
#   ggplot(aes(x = time, colour = factor(group))) +
#   geom_density() +
#   labs(
#     x = "Time",
#     y = "Density",
#     colour = "Group"
#   ) +
#   evc::evc_theme()
#
# # fit
# coef(coxph(Surv(time, event) ~ group, data = data))

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

  mclapply(seq_len(n_sims), \(j) {
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

cox_res <- mclapply(sim_grid_lst, \(x) {
  # fit Cox to each, extract coefficient
  sapply(x, \(y) unname(coef(coxph(Surv(time, event) ~ group, data = y))))
}, mc.cores = cores)

# join to grid
grid_res_ph_df <- bind_rows(lapply(seq_along(cox_res), \(i) {
  data.frame(grid[i, ], coefs = cox_res[[i]])
}))

grid_res_ph_df <- grid_res_ph_df |>
  filter(abs(coefs) < 5) |>
  group_by(time, n, censor, treat) |>
  mutate(
    mean_coefs = mean(coefs),
    bias = abs(mean(coefs - beta1)),
    var = var(coefs)
  )

pgrid_ph <- grid_res_ph_df |>
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
pgrid_ph

ggsave(plot = pgrid_ph, filename = "report/plots/02b_sim_res_ph_cox_grid.png", width = 10, height = 8)

grid_res_ph_df |>
  distinct(time, n, mean_coefs, bias, var)

# TODO Run PEM on each of these also (will be slow)
pem_res <- mclapply(seq_along(sim_grid_lst), \(i) {
  # fit pem to each, extract coefficient
  # sapply(x, \(y) unname(coef(pempem(Surv(time, event) ~ group, data = y))))
  print(sprintf("Row %s of %s", i, nrow(grid)))
  sapply(seq_along(sim_grid_lst[[i]]), \(j) {
    system(sprintf("echo %s", paste("Simulation", j, "of", length(sim_grid_lst[[i]]))))
    fit <- fit_pem(data = sim_grid_lst[[i]][[j]], times = times)
    unname(coef(fit)["group"])
  })
}, mc.cores = cores)

saveRDS(pem_res, file = "pem_res.RDS")

# join to grid
grid_res_pem_df <- bind_rows(lapply(seq_along(pem_res), \(i) {
  data.frame(grid[i, ], coefs = pem_res[[i]])
}))

grid_res_pem_df <- grid_res_pem_df |>
  filter(abs(coefs) < 5) |>
  group_by(time, n, censor, treat) |>
  mutate(
    mean_coefs = mean(coefs),
    bias = abs(mean(coefs - beta1)),
    var = var(coefs)
  )

pgrid_pem <- grid_res_pem_df |>
  ggplot(aes(x = coefs)) +
  geom_histogram(colour = "black", fill = "grey") +
  labs(
    x = latex2exp::TeX("$\\beta$"),
    y = "Count"
  ) +
  # geom_vline(xintercept = mean(grid_res_pem_df$coefs), col = "blue", lty = 1, lwd = 1.2) +
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
pgrid_pem

ggsave(plot = pgrid_pem, filename = "report/plots/02b_sim_res_pem_cox_grid.png", width = 10, height = 8)

# #### Proportional Odds test case ####
#
# # simulate data under PO model
# set.seed(123)
# data_po <- rpwexp_q_po(
#   n = 1000,
#   times = times,
#   quantiles = c(0.25, 0.5, 0.75),
#   beta = beta1,
#   x = rbinom(1000, 1, 0.5) # Random treatment assignment
# )
#
# # add right censoring
# data_po <- data_po %>%
#   mutate(
#     event = rbinom(n(), 1, 0.8), # censoring indicator
#   )
#
# # plot data
# data_po |>
#   ggplot(aes(x = time, colour = factor(group))) +
#   geom_density() +
#   labs(
#     x = "Time",
#     y = "Density",
#     colour = "Group"
#   ) +
#   evc::evc_theme()
#
# # try proportional odds
# # TODO Fix, why is this not working?? Or is it just a fluke?
# po <- timereg::prop.odds(Event(time, event) ~ group, data = data_po)
# coef(po)
#
# # not bad either, even though we have PO model!
# coef(coxph(Surv(time, event) ~ group, data = data_po))
#
# # see how PEM does
# coef(fit_pem(data_po, times))["group"] # does better! (For this example)
#
#
# #### PO simulations ####
#
# # generate simulated datasets
# n_sims <- 500
# set.seed(123)
# sim_data_po_lst <- mclapply(seq_len(n_sims), \(i){
#   rpwexp_q_po(
#     # n = 1000,
#     n = 200,
#     times = times,
#     quantiles = quantiles,
#     beta = beta1,
#     x = rbinom(1000, 1, 0.5)
#   ) |>
#     mutate(
#       event = rbinom(n(), 1, 0.8)
#     )
# }, mc.cores = cores)
#
# # under PO
# coefs_po <- unlist(mclapply(sim_data_po_lst, \(x) {
#   coef(timereg::prop.odds(Event(time, event) ~ group, data = x))
#   # coef(timereg::prop.odds(Event(time, cens.code = 1) ~ group, data = x))
#   # coef(timereg::prop.odds(Event(time) ~ group, data = x))
# }, mc.cores = cores))
#
# hist(coefs_po)
# abline(v = beta1, col = "red", lty = 2, lwd = 2)
#
# # fit Cox, evaluate performance
# coefs_cox_po <- unlist(mclapply(sim_data_po_lst, \(x) {
#   coef(coxph(Surv(time, event) ~ group, data = x))
# }, mc.cores = cores))
#
# hist(coefs_cox_po)
# abline(v = beta1, col = "red", lty = 2, lwd = 2)
#
# mean(coefs_cox_po - beta1) # bias is larger by factor of 100
# var(coefs_cox_po) # variance also larger by factor of 10
#
# # now try for PEM
# coefs_pem_po <- unlist(mclapply(sim_data_po_lst, \(x) {
#   coef(fit_pem(x, times))["group"]
# }, mc.cores = cores))
#
# hist(coefs_pem_po)
# abline(v = beta1, col = "red", lty = 2, lwd = 2)
# # extremely precise! ~0.0002 => better than Cox by factor of 100
# mean(coefs_pem_po - beta1)
# var(coefs_pem_po) # However, similar variance to Cox here
