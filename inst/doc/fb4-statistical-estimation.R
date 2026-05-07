## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 5,
  message   = FALSE,
  warning   = FALSE
)
library(fb4package)

## ----bio-setup----------------------------------------------------------------
data(fish4_parameters)

chinook_db <- fish4_parameters[["Oncorhynchus tshawytscha"]]
stage      <- names(chinook_db$life_stages)[1]
sp_params  <- chinook_db$life_stages[[stage]]
sp_info    <- chinook_db$species_info

days <- 1:100
temp_data <- data.frame(
  Day         = days,
  Temperature = round(10 + 5 * sin(2 * pi * (days - 60) / 365), 2)
)
diet_data <- list(
  proportions = data.frame(Day = days, Alewife = 0.7, Shrimp = 0.3),
  prey_names  = c("Alewife", "Shrimp"),
  energies    = data.frame(Day = days, Alewife = 4900, Shrimp = 3200)
)

bio <- Bioenergetic(
  species_params      = sp_params,
  species_info        = sp_info,
  environmental_data  = list(temperature = temp_data),
  diet_data           = diet_data,
  simulation_settings = list(initial_weight = 1800, duration = 100)
)
bio$species_params$predator$ED_ini <- 4500
bio$species_params$predator$ED_end <- 5500

## ----obs-weights--------------------------------------------------------------
set.seed(42)
obs_weights <- rnorm(30, mean = 3000, sd = 80)

## ----binary-search------------------------------------------------------------
res_bs <- run_fb4(bio,
                  strategy  = "binary_search",
                  fit_to    = "Weight",
                  fit_value = 3000)

cat("p estimate  :", round(res_bs$summary$p_value, 4), "\n")
cat("Final weight:", round(res_bs$summary$final_weight, 1), "g\n")

## ----mle, cache=TRUE----------------------------------------------------------
res_mle <- run_fb4(bio,
                   strategy         = "mle",
                   fit_to           = "Weight",
                   observed_weights = obs_weights)

cat("p estimate :", round(res_mle$summary$p_value, 4), "\n")
cat("SE         :", round(res_mle$summary$p_se,    4), "\n")
cat("Converged  :", res_mle$summary$converged, "\n")

## ----bootstrap, cache=TRUE----------------------------------------------------
res_boot <- run_fb4(bio,
                    strategy         = "bootstrap",
                    fit_to           = "Weight",
                    observed_weights = obs_weights,
                    n_bootstrap      = 100,
                    upper            = 1,
                    parallel         = FALSE)

cat("p mean :", round(res_boot$summary$p_mean, 4), "\n")
cat("p SD   :", round(res_boot$summary$p_sd,   4), "\n")

## ----mr-ref-p-----------------------------------------------------------------
res_ref <- run_fb4(bio,
                   strategy  = "binary_search",
                   fit_to    = "Weight",
                   fit_value = 2500)
ref_p <- res_ref$summary$p_value
cat("Reference p (target 2 500 g):", round(ref_p, 4), "\n")

## ----mr-simulate--------------------------------------------------------------
set.seed(123)
n_fish <- 20

# Cohort: individual initial weights varying ± 15 % around 1 800 g
initial_weights <- round(runif(n_fish, min = 1800 * 0.85, max = 1800 * 1.15))

# Each fish has its own feeding level drawn from a population distribution
# centred on the reference p estimated above (sigma_p = 0.06)
true_p <- pmax(0.3, pmin(1.5, rnorm(n_fish, mean = ref_p, sd = 0.06)))

# Simulate final weights: run individual direct simulations + measurement noise
final_weights <- vapply(seq_len(n_fish), function(i) {
  bio_i <- bio
  bio_i$simulation_settings$initial_weight <- initial_weights[i]
  res_i <- run_fb4(bio_i, strategy = "direct", p_value = true_p[i])
  round(res_i$summary$final_weight * rnorm(1, mean = 1, sd = 0.02), 1)
}, numeric(1))

mr_data <- data.frame(
  individual_id  = sprintf("fish_%02d", seq_len(n_fish)),
  initial_weight = initial_weights,
  final_weight   = final_weights
)

head(mr_data, 6)

## ----hierarchical, cache=TRUE, dependson=c("mr-ref-p", "mr-simulate")---------
res_hier <- run_fb4(
  x                = bio,
  strategy         = "hierarchical",
  backend          = "tmb",
  fit_to           = "Weight",
  observed_weights = mr_data
)

cat(sprintf("Population mean p  (mu_p)   : %.4f\n", res_hier$summary$mu_p_estimate))
cat(sprintf("Population SD  p  (sigma_p) : %.4f\n", res_hier$summary$sigma_p_estimate))
cat(sprintf("Number of individuals       : %d\n",   res_hier$summary$n_individuals))
cat(sprintf("Converged                   : %s\n",   res_hier$summary$converged))

## ----compare-123, echo=FALSE--------------------------------------------------
comp123 <- data.frame(
  Strategy    = c("binary_search", "mle", "bootstrap"),
  p_estimate  = round(c(res_bs$summary$p_value,
                         res_mle$summary$p_value,
                         res_boot$summary$p_mean), 4),
  Uncertainty = round(c(NA,
                         res_mle$summary$p_se,
                         res_boot$summary$p_sd), 4),
  Type        = c("none (point estimate)",
                  "SE (asymptotic)",
                  "SD (non-parametric)"),
  stringsAsFactors = FALSE
)
knitr::kable(comp123,
             col.names = c("Strategy", "p estimate", "SE / SD",
                           "Uncertainty type"),
             align  = c("l", "r", "r", "l"),
             caption = "Strategies 1–3 fitted to the same data (Chinook salmon, 100-day simulation, target weight ≈ 3 000 g).")

## ----compare-hier, echo=FALSE-------------------------------------------------
comp_hier <- data.frame(
  Parameter   = c("mu_p (population mean)",
                  "sigma_p (population SD)",
                  "n_individuals",
                  "converged"),
  Value       = c(round(res_hier$summary$mu_p_estimate,    4),
                  round(res_hier$summary$sigma_p_estimate, 4),
                  res_hier$summary$n_individuals,
                  as.character(res_hier$summary$converged)),
  stringsAsFactors = FALSE
)
knitr::kable(comp_hier,
             col.names = c("Parameter", "Estimated value"),
             align     = c("l", "r"),
             caption   = paste("Hierarchical model results for the simulated",
                               "mark-recapture cohort (n =",
                               res_hier$summary$n_individuals,
                               "fish, true mu_p =", round(ref_p, 3), ")."))

