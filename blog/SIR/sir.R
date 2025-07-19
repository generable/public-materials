library(ggplot2)
library(deSolve)
library(tidyr)
theme_set(theme_minimal())

# S: susceptible, I: infected, R: recovered
SIR <- function(t, state, pars) {
  with(as.list(c(state, pars)), {
    N <- S + I + R
    dS_dt <- -b * I / N * S
    dI_dt <-  b * I / N * S - g * I
    dR_dt <-  g * I
    return(list(c(dS_dt, dI_dt, dR_dt)))
  })
}

n <- 100
n_pop <- 763
pars <- c(b = 2, g = 0.5)
state <- c(S = n_pop - 1, I = 1, R = 0)
times <- seq(1, 15, length = n)
sol <-
  ode(
    y = state,
    times = times,
    func = SIR,
    parms = pars,
    method = "ode45"
  )
sol <- as.data.frame(sol)
sol_long <- sol |>
  tidyr::pivot_longer(-time, names_to = "state", values_to = "value")

sol_long |>
  ggplot(aes(time, value, color = state)) +
  geom_line() +
  guides(color = guide_legend(title = NULL)) +
  scale_color_discrete(labels = c("Infected", "Recovered", "Susceptible")) +
  xlab("Days") + ylab("Number of people") +
  ggtitle("Basic SIR model", subtitle = "3-state ODE with beta = 2,
          gamma = 0.5, R0 = 4")

set.seed(1234)
y <- rpois(n, sol$I)

library(cmdstanr)
library(posterior)

data <- list(n_obs = n, n_pop = n_pop, y = y, t0 = 0, ts = times)
sir_mod <- cmdstan_model("sir.stan")
sir_fit <- sir_mod$sample(
  data = data, 
  seed = 1234,
  chains = 4, 
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500
)

theta <- as_draws_rvars(sir_fit$draws(variables = "theta"))
summarise_draws(theta, default_convergence_measures())
summarise_draws(theta, ~ quantile(.x, probs = c(0.05, 0.5, 0.95)))

y_hat <- as_draws_rvars(sir_fit$draws(variables = "y_hat"))

# helper function used for plotting the results
make_df <- function(d, interval = .90, t, obs, pop) {
  S <- d$y_hat[, 1] # Susceptible draws, not used
  I <- d$y_hat[, 2] # Infected draws
  R <- d$y_hat[, 3] # Recovered draws, not used
  
  # compute the uncertainty interval
  low_quant <- (1 - interval) / 2
  high_quant <- interval + low_quant
  low <- apply(I, 2, quantile, probs = low_quant) * pop
  high <- apply(I, 2, quantile, probs = high_quant) * pop
  
  return(tibble(low, high, times = t, obs))
}

d <- make_df(d = y_hat, interval = 0.90, times, obs = sol$I, pop = n_pop)
ggplot(aes(times, obs), data = d) +
  geom_ribbon(aes(ymin = low, ymax = high), fill = "grey70") +
  geom_line(color = "red", linewidth = 0.3) +
  xlab("Days") + ylab("Infections (90% Uncertainty)") +
  ggtitle("SIR Estimation",
          subtitle = "Red curve is the true rate")

pct_train <- 0.30
n_train <- floor(n * pct_train)
n_pred <- n - n_train
times_pred <- times[(n_train + 1):n]
y_train <- y[1:n_train]
data <- list(n_obs = n_train, n_pred = n_pred, 
             n_pop = n_pop, y = y_train, 
             t0 = 0, ts = times[1:n_train], ts_pred = times_pred)


sir_pred_mod <- cmdstan_model("sir-pred.stan")
sir_pred_fit <- sir_pred_mod$sample(
  data = data, 
  seed = 1234,
  chains = 4, 
  parallel_chains = 4,
  iter_warmup = 500,
  iter_sampling = 500,
  adapt_delta = 0.98
)

y_hat <- as_draws_rvars(sir_pred_fit$draws(variables = "y_hat"))
y_hat_pred <- as_draws_rvars(sir_pred_fit$draws(variables = "y_hat_pred"))

d_train <- make_df(y_hat, 0.90, data$ts, sol$I[1:n_train], n_pop)
d_pred <- make_df(y_hat_pred, 0.90, data$ts_pred, sol$I[(n_train + 1):n], n_pop)
d <- dplyr::bind_rows(d_train, d_pred) 

ggplot(aes(times, obs), data = d) +
  geom_ribbon(aes(ymin = low, ymax = high), 
              fill = "grey70", alpha = 1/2) +
  geom_line(color = "red", linewidth = 0.3) +
  geom_vline(xintercept = d$times[n_train], linetype = "dotdash") +
  xlab("Days") + ylab("Infections (90% Uncertainty)") +
  ggtitle("SIR prediction trained on 5 days of data",
          subtitle = "Red curve is the true rate")
