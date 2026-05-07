## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse  = TRUE,
  comment   = "#>",
  fig.width = 7,
  fig.height = 5,
  message   = FALSE,
  warning   = FALSE
)
library(fb4package)

## ----species-db---------------------------------------------------------------
data(fish4_parameters)
cat("Species in database:", length(fish4_parameters), "\n")
head(names(fish4_parameters), 8)

## ----chinook-params-----------------------------------------------------------
chinook_db  <- fish4_parameters[["Oncorhynchus tshawytscha"]]
stage       <- names(chinook_db$life_stages)[1]
cat("Life stage used:", stage, "\n")

sp_params   <- chinook_db$life_stages[[stage]]
sp_info     <- chinook_db$species_info

# Key consumption parameters
cat("CEQ:", sp_params$consumption$CEQ,
    "  CA:", sp_params$consumption$CA,
    "  CB:", sp_params$consumption$CB, "\n")

## ----env-diet-----------------------------------------------------------------
days <- 1:100

# Seasonal temperature (°C)
temp_data <- data.frame(
  Day         = days,
  Temperature = round(10 + 5 * sin(2 * pi * (days - 60) / 365), 2)
)

# Two-prey diet (constant proportions, J/g energy densities)
diet_data <- list(
  proportions = data.frame(Day = days, Alewife = 0.7, Shrimp = 0.3),
  prey_names  = c("Alewife", "Shrimp"),
  energies    = data.frame(Day = days, Alewife = 4900, Shrimp = 3200)
)

head(temp_data, 3)

## ----bio-object---------------------------------------------------------------
bio <- Bioenergetic(
  species_params      = sp_params,
  species_info        = sp_info,
  environmental_data  = list(temperature = temp_data),
  diet_data           = diet_data,
  simulation_settings = list(initial_weight = 1800, duration = 100)
)

# Required when PREDEDEQ = 2 (energy density as function of weight)
bio$species_params$predator$ED_ini <- 4500
bio$species_params$predator$ED_end <- 5500

print(bio)

## ----direct-------------------------------------------------------------------
res_direct <- run_fb4(bio, strategy = "direct", p_value = 0.8)
print(res_direct)

## ----binary-search------------------------------------------------------------
res_bs <- run_fb4(bio,
                  strategy  = "binary_search",
                  fit_to    = "Weight",
                  fit_value = 3000)

cat("Estimated p-value :", round(res_bs$summary$p_value, 4), "\n")
cat("Final weight      :", round(res_bs$summary$final_weight, 1), "g\n")
cat("Total consumption :", round(res_bs$summary$total_consumption_g, 1), "g\n")

## ----obs-weights-def, include=FALSE-------------------------------------------
set.seed(42)
obs_weights <- rnorm(30, mean = 3000, sd = 80)

## ----mle, cache=TRUE, dependson="obs-weights-def"-----------------------------
set.seed(42)
obs_weights <- rnorm(30, mean = 3000, sd = 80)   # 30 observed final weights

res_mle <- run_fb4(bio,
                   strategy         = "mle",
                   fit_to           = "Weight",
                   observed_weights = obs_weights)

cat("p-value (MLE)  :", round(res_mle$summary$p_value, 4), "\n")
cat("SE             :", round(res_mle$summary$p_se,    4), "\n")
cat("Converged      :", res_mle$summary$converged, "\n")

## ----bootstrap, cache=TRUE, dependson=c("obs-weights-def", "mle")-------------
set.seed(42)
res_boot <- run_fb4(bio,
                    strategy         = "bootstrap",
                    fit_to           = "Weight",
                    observed_weights = obs_weights,
                    upper            = 1,
                    n_bootstrap      = 100,
                    parallel         = FALSE)

cat("p mean (bootstrap):", round(res_boot$summary$p_mean, 4), "\n")
cat("p SD              :", round(res_boot$summary$p_sd,   4), "\n")

## ----compare, echo=FALSE------------------------------------------------------
comparison <- data.frame(
  Strategy        = c("binary_search", "mle", "bootstrap"),
  p_estimate      = round(c(res_bs$summary$p_value,
                             res_mle$summary$p_value,
                             res_boot$summary$p_mean), 4),
  SE_or_SD        = round(c(NA, res_mle$summary$p_se, res_boot$summary$p_sd), 4),
  Uncertainty     = c("none (point estimate)", "SE (asymptotic)", "SD (non-parametric)"),
  stringsAsFactors = FALSE
)
knitr::kable(comparison,
             col.names = c("Strategy", "p estimate", "SE / SD", "Uncertainty type"),
             align     = c("l", "r", "r", "l"),
             caption   = "All strategies use the same bio object (Chinook salmon, 100-day simulation). MLE and bootstrap fit to obs_weights; binary search fits to a target weight.")

## ----analysis-----------------------------------------------------------------
growth   <- analyze_growth_patterns(res_bs)
feeding  <- analyze_feeding_performance(res_bs)
energy   <- analyze_energy_budget(res_bs)

cat(sprintf("Final weight              : %.1f g\n",
            growth$final_weight$estimate))
cat(sprintf("Specific growth rate      : %.4f g/g/day\n",
            growth$specific_growth_rate$estimate))
cat(sprintf("Total consumption         : %.1f g\n",
            feeding$total_consumption$estimate))
cat(sprintf("Gross conversion efficiency: %.3f\n",
            feeding$gross_conversion_efficiency$estimate))

## ----plot-growth, fig.cap="Daily growth trajectory from binary search fit."----
plot(res_bs, type = "growth")

## ----plot-energy, fig.cap="Daily energy budget partitioning."-----------------
plot(res_bs, type = "energy")

## ----plot-dashboard, fig.cap="Full dashboard (growth, consumption, temperature, energy)."----
plot(res_bs, type = "dashboard")

## ----sensitivity, fig.cap="Temperature and p-value sensitivity grid.", message=FALSE, warning=FALSE----
plot(bio,
     type         = "sensitivity",
     temperatures = seq(6, 20, by = 2),
     p_values     = seq(0.3, 1.0, by = 0.1))

