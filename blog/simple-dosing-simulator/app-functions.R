library(tidyverse)
logit  <- function(p) log(p/(1-p))
invlogit <- function(x) 1/(1+exp(-x))

# Hill/Emax term in [0, 1): grows with dose, saturates at 1 as dose -> inf
hill_term <- function(dose, ed50, hill = 1) (dose^hill) / (ed50^hill + dose^hill)

# Probability models ----------------------------------------------------
# Monotone ↑ bleed with dose; monotone ↓ stroke with dose (on logit scale)
p_bleed <- function(dose,
                    p0 = 0.03,           # baseline bleed risk at dose≈0
                    emax = log(2),       # max increase on log-odds (≈2x odds)
                    ed50 = 2,            # dose for half-max effect
                    hill = 1.5) {
  invlogit( logit(p0) + emax * hill_term(dose, ed50, hill) )
}

p_stroke <- function(dose,
                     p0 = 0.08,          # baseline stroke risk at dose≈0
                     emax = log(2),      # max DECREASE on log-odds
                     ed50 = 2.5,
                     hill = 1.2) {
  invlogit( logit(p0) - emax * hill_term(dose, ed50, hill) )
}

sim_risks <- function(p0_bleed = 0.03, p0_stroke = 0.08, 
                      doses = c(0, 0.5, 1, 2, 3, 5),
                      weights = c('Bleed' = 0.4, 'Stroke' = 0.8)) {
  grid <- data.frame(dose = seq(0, max(doses), length.out = 200))
  grid$p_bleed  <- p_bleed(grid$dose,  p0 = p0_bleed, emax = log(3),  ed50 = 2,   hill = 1.4)
  grid$p_stroke <- p_stroke(grid$dose, p0 = p0_stroke, emax = log(2.5), ed50 = 2.5, hill = 1.1)
  
  p_risks <- grid |>
    gather(event, probability, starts_with('p_')) |>
    mutate(event = str_to_title(str_replace_all(str_remove(event, '^p_'), '_', ' ')),
           dose = 10*dose) |>
    rename(Event = event) |>
    ggplot(aes(x = dose, y = probability, colour = Event, linetype = Event, group = Event)) +
    geom_line(size = 1.2) +
    scale_x_continuous('Dose (mg)') +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_manual(values = c('Bleed' = 'skyblue', 'Stroke' = 'purple')) +
    scale_linetype_manual(values = c('Bleed' = 1, 'Stroke' = 2)) +
    theme_minimal() +
    labs(y = "Probability", title = "Risk of bleed (solid blue, ↑) and stroke (dashed purple, ↓) vs Dose")
  
  p_risks2 <- grid |>
    gather(event, probability, starts_with('p_')) |>
    mutate(event = str_to_title(str_replace_all(str_remove(event, '^p_'), '_', ' ')),
           dose = 10*dose) |>
    rename(Event = event) |>
    group_by(dose) |>
    mutate(cumprob = sum(probability)) |>
    ungroup() |>
    ggplot(aes(x = dose, y = probability, colour = Event, group = Event)) +
    geom_line(size = 1.2) +
    geom_line(aes(y = cumprob, colour = 'Composite'),
              size = 1.2) +
    scale_x_continuous('Dose (mg)') +
    scale_y_continuous(labels = scales::percent) +
    scale_colour_manual(values = c('Bleed' = 'skyblue', 'Stroke' = 'purple', 'Composite' = 'orange')) +
    theme_minimal() +
    labs(y = "Probability", title = "Risk of bleed (solid blue, ↑) and stroke (dashed purple, ↓) vs Dose")
  
  p_composite <- grid |>
    mutate(probability = p_stroke + p_bleed,
           dose = 10*dose) |>
    ggplot(aes(x = dose, y = probability, colour = 'Composite', group = 'Composite')) +
    geom_line(size = 1.2) +
    scale_x_continuous('Dose (mg)') +
    scale_y_continuous(labels = scales::percent) +
    theme_minimal() +
    labs(y = "Probability of Bleed or Stroke", title = "Risk of composite endpoint (any event) vs Dose", colour = 'Event')
  
  event_weights <- list('weight' = weights) |>
    as.data.frame() |>
    rownames_to_column('event')
  
  p_utility <- grid |>
    gather(event, probability, starts_with('p_')) |>
    mutate(event = str_to_title(str_replace_all(str_remove(event, '^p_'), '_', ' ')),
           dose = 10*dose) |>
    left_join(event_weights) |>
    group_by(dose) |>
    summarise(harmfulness = sum(probability * weight)) |>
    ungroup() |>
    ggplot(aes(x = dose, y = harmfulness)) + 
    geom_line(size = 1.2) +
    theme_minimal() +
    scale_x_continuous('Dose (mg)') +
    scale_y_continuous('Harmfulness', labels = scales::percent) +
    labs(y = "Probability", title = "Risk of harm (weighted event rate) vs Dose", colour = 'Event')
  list(risks = p_risks, risks2 = p_risks2, composite = p_composite, utility = p_utility)
}