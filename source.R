#### Functions for report ####

pot_theme <- \() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      legend.position = "bottom",
      plot.title = ggplot2::element_text(size = 16, hjust = 0.5),
      axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(
        size = 14,
        face = "bold"
      ), legend.text = ggplot2::element_text(size = 12),
      strip.text = ggplot2::element_text(size = 13, face = "bold"),
      strip.background = ggplot2::element_rect(fill = NA, colour = "black"),
      plot.tag = ggplot2::element_text(size = 16, face = "bold"),
      panel.background = ggplot2::element_rect(fill = NA, colour = "black")
    )
}

#### Functions ####

# function to plot histogram of parameter estimates from models
hist_fun <- \(coefs, beta = 0.5) {
  x <- coefs
  if (!is.data.frame(x)) {
    x <- data.frame("coef" = x)
  }
  stopifnot(is.data.frame(x) && names(x) == "coef")
  x |>
    ggplot(aes(x = coef)) +
    geom_histogram(colour = "black", fill = "#D3D3D3") +
    labs(
      x = latex2exp::TeX("$\\beta$"),
      y = "Count"
    ) +
    geom_vline(xintercept = mean(x$coef), col = "blue", lty = 1, lwd = 1.2) +
    geom_vline(xintercept = beta, col = "red", lty = 2, lwd = 1.2) +
    pot_theme() +
    NULL
}

# function to plot density with quantiles
plot_q <- \(x, quantiles = c(0.25, 0.5, 0.75), ...) {
  plot(density(x), ...)
  if (!is.null(quantiles)) {
    qs <- quantile(x, probs = quantiles)
    for (q in qs) {
      abline(v = q, lty = 3, col = "red")
    }
  }
}

# function to calculate piecewise exponential hazards
haz <- \(times, quantiles) {
  stopifnot(length(times) == length(quantiles))
  stopifnot(quantiles == sort(quantiles))

  # need quantiles of survival function
  surv_q <- 1 - quantiles
  haz <- numeric(length(times))
  for (i in seq_along(times)) {
    if (i == 1) {
      haz[i] <- -log(1 - quantiles[i]) / times[i]
    } else {
      haz[i] <- -log(surv_q[i] / surv_q[i - 1]) /
        (times[i] - times[i - 1])
    }
  }
  return(haz)
}

# function to simulate from piecewise exponential distribution (under PH mod)
rpwexp_q_ph <- \(
  n, # number of observations
  times, # breakpoints
  quantiles, # quantiles of CDF (reversed for S)
  beta = 0, # coefficient for X
  x = NULL
) {
  # Calculate baseline hazard rates
  haz_0 <- haz(times, quantiles)

  # only binomial x supported for now
  stopifnot(all(x %in% c(0, 1)))
  stopifnot(length(x) == n) # x must be same length as n!

  # Generate survival times for each group
  t <- numeric(n)

  # Use baseline hazard for x = 0
  n_x0 <- sum(x == 0)
  if (n_x0 > 0) {
    t[x == 0] <- nphsim::rpwexp(
      n         = n_x0,
      rate      = haz_0,
      intervals = times
    )
  }

  # use PH specification of hazard for x = 1
  n_x1 <- sum(x == 1)
  if (n_x1 > 0) {
    # TODO Or is it minus beta??
    haz_1 <- haz_0 * exp(beta) # PH model
    # t[x == 1] <- VirtualPop::r.pw_exp(
    #   n = n_x1,
    #   breakpoints = times,
    #   rates = haz_1
    # )
    t[x == 1] <- nphsim::rpwexp(
      n         = n_x1,
      rate      = haz_1,
      intervals = times
    )
  }

  # Return data frame with times and covariates
  data.frame(time = t, group = x)
}

# fit PEM model using flexsurv
fit_pem <- \(
  ped_form = survival::Surv(time, event) ~ group, # formula for ped trans
  flex_form = survival::Surv( # formula for flexsurv model
    tstart, tend, ped_status
  ) ~ interval + group + offset,
  data, # data in time-to-event form
  times # cutpoints for trans
) {
  # convert to ped object
  data_ped <- pammtools::as_ped(data = data, formula = ped_form, cut = times)

  # initial parameter estimation (gives starting values)
  tmp <- flexsurv::flexsurvreg(
    formula = flex_form,
    data    = data_ped,
    dist    = "exp"
  )

  # pull start values, ensure rate > 0 & set offset coefficient is 1
  init_vals <- coef(tmp)
  init_vals["rate"] <- exp(init_vals["rate"])
  init_vals["offset"] <- 1

  # Fit while fixing offset parameter at 1
  offset_idx <- which(names(init_vals) == "offset")
  pem_fit <- flexsurv::flexsurvreg(
    formula   = flex_form,
    data      = data_ped,
    dist      = "exp",
    inits     = init_vals,
    fixedpars = offset_idx
  )

  pem_fit$data <- data_ped # add ped data to fit object
  return(pem_fit)
}

# function to plot density with quantiles
plot_q <- \(x, quantiles = c(0.25, 0.5, 0.75), ...) {
  plot(density(x), ...)
  if (!is.null(quantiles)) {
    qs <- quantile(x, probs = quantiles)
    for (q in qs) {
      abline(v = q, lty = 3, col = "red")
    }
  }
}


#### sourced functions from Sam ####

# TODO Add hll, hg, Hg, Hlleg, Hll
hll <- function(t) {
  1 / (1 + t)
}

hg <- function(t) {
  exp(t)
}

Hg <- function(t) {
  exp(t) - 1
}

Hlleg <- function(t, psi) {
  log(1 + psi * (exp(t) - 1))
}


Hll <- function(t) {
  log(1 + t)
}

#### Weibull PO ####

hwbPO <- function(t, q1, median, psi) {
  shape <- (log(log(4 / 3)) - log(log(2))) / (log(q1 / median))
  scale <- median / ((log(2))^(1 / shape))
  hll(psi * Hg((t / scale)^shape)) * psi * hg((t / scale)^shape) * (shape / scale) * ((t / scale)^(shape - 1))
}

HwbPO <- function(t, q1, median, psi) {
  shape <- (log(log(4 / 3)) - log(log(2))) / (log(q1 / median))
  scale <- median / ((log(2))^(1 / shape))
  Hlleg((t / scale)^shape, psi)
}

hwbPO <- function(t, q1, median, psi) {
  shape <- (log(log(4 / 3)) - log(log(2))) / (log(q1 / median))
  scale <- median / ((log(2))^(1 / shape))
  hll(psi * Hg((t / scale)^shape)) * psi * hg((t / scale)^shape) * (shape / scale) * ((t / scale)^(shape - 1))
}

HwbPO <- function(t, q1, median, psi) {
  shape <- (log(log(4 / 3)) - log(log(2))) / (log(q1 / median))
  scale <- median / ((log(2))^(1 / shape))
  Hlleg((t / scale)^shape, psi)
}

#### Log-logistic PO for flexsurv models ####

hllgPO <- function(t, q1, median, psi) {
  scale <- median
  shape <- log(3) / log(median / q1)
  hll(psi * (t / scale)^shape) * shape * psi * (t^(shape - 1)) / (scale^shape)
}

HllgPO <- function(t, q1, median, psi) {
  scale <- median
  shape <- log(3) / log(median / q1)
  Hll(psi * (t / scale)^shape)
}

Hinv_llgPO <- function(z, q1, median, psi) {
  scale <- median
  shape <- log(3) / log(median / q1)
  scale * (Hg(z) / psi)^(1 / shape)
}


#### Gompertz PO ####

hgompPO <- function(t, q1, median, psi) {
  shape <- log(log(1 + log(4 / 3)) / log(1 + log(2))) / log(q1 / median)
  scale <- median / ((log(1 + log(2)))^(1 / shape))
  psi * hll(psi * Hg(Hg((t / scale)^shape))) *
    hg(Hg((t / scale)^shape)) *
    hg((t / scale)^shape) *
    (shape / scale) * ((t / scale)^(shape - 1))
}

HgompPO <- function(t, q1, median, psi) {
  shape <- log(log(1 + log(4 / 3)) / log(1 + log(2))) / log(q1 / median)
  scale <- median / ((log(1 + log(2)))^(1 / shape))
  Hlleg(exp((t / scale)^shape) - 1, psi)
}

Hinv_gompPO <- function(z, q1, median, psi) {
  shape <- log(log(1 + log(4 / 3)) / log(1 + log(2))) / log(q1 / median)
  scale <- median / ((log(1 + log(2)))^(1 / shape))
  scale * Hll(Hlleg(z, 1 / psi))^(1 / shape)
}

#### Custom flexsurv distributions ####

custom.llgPO <-
  list(
    name = "llgPO",
    pars = c("psi", "q1", "median"),
    location = "psi",
    transforms = c(log, log, log),
    inv.transforms = c(exp, exp, exp),
    inits = function(t) {
      c(1, 1, 1)
    }
  )

custom.wbPO <-
  list(
    name = "wbPO",
    pars = c("psi", "q1", "median"),
    location = "psi",
    transforms = c(log, log, log),
    inv.transforms = c(exp, exp, exp),
    inits = function(t) {
      c(1, 1, 1)
    }
  )

custom.gompPO <-
  list(
    name = "gompPO",
    pars = c("psi", "q1", "median"),
    location = "psi",
    transforms = c(log, log, log),
    inv.transforms = c(exp, exp, exp),
    inits = function(t) {
      c(1, 1, 1)
    }
  )
