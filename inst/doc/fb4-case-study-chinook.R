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

## ----species-params-----------------------------------------------------------
data(fish4_parameters)

chinook_db <- fish4_parameters[["Oncorhynchus tshawytscha"]]
stage      <- if ("juvenile" %in% names(chinook_db$life_stages)) {
                "juvenile"
              } else {
                names(chinook_db$life_stages)[1]
              }

sp_params <- chinook_db$life_stages[[stage]]
sp_info   <- chinook_db$species_info
sp_info$life_stage <- stage

cat("Life stage  :", stage, "\n")
cat("CEQ (consumption equation):", sp_params$consumption$CEQ, "\n")
cat("REQ (respiration equation):", sp_params$respiration$REQ, "\n")

## ----temperature--------------------------------------------------------------
set.seed(42)
days      <- 1:180
base_temp <- 7 + 7 * sin(pi * (days - 20) / 180)   # peak ~14 °C at day 110
temp_vec  <- pmax(2, base_temp + rnorm(180, 0, 0.4))

temp_data <- data.frame(
  Day         = days,
  Temperature = round(temp_vec, 2)
)

cat(sprintf("Temperature range : %.1f – %.1f °C\n",
            min(temp_data$Temperature), max(temp_data$Temperature)))
cat(sprintf("Mean temperature  : %.1f °C\n", mean(temp_data$Temperature)))

## ----diet---------------------------------------------------------------------
alewife <- pmax(0, 0.55 + 0.25 * sin(pi * (days - 30) / 180))
shrimp  <- pmax(0, 0.28 - 0.10 * sin(pi * (days - 30) / 180))
inverts <- pmax(0, 1 - alewife - shrimp)
total   <- alewife + shrimp + inverts

diet_props <- data.frame(
  Day     = days,
  Alewife = round(alewife / total, 4),
  Shrimp  = round(shrimp  / total, 4),
  Inverts = round(inverts / total, 4)
)

prey_energy <- data.frame(
  Day     = days,
  Alewife = 4900,   # J/g (wet weight)
  Shrimp  = 3200,
  Inverts = 2600
)

cat("Diet proportions (first 3 days):\n")
print(head(diet_props, 3))

## ----bio-object---------------------------------------------------------------
bio_chinook <- Bioenergetic(
  species_params     = sp_params,
  species_info       = sp_info,
  environmental_data = list(temperature = temp_data),
  diet_data          = list(
    proportions = diet_props,
    prey_names  = c("Alewife", "Shrimp", "Inverts"),
    energies    = prey_energy
  ),
  simulation_settings = list(initial_weight = 5, duration = 180)
)

# Predator energy density: linear interpolation from 4 200 to 5 000 J/g
# as the fish accumulate lipids through summer (PREDEDEQ = 1 from DB)
bio_chinook$species_params$predator$ED_ini <- 4200
bio_chinook$species_params$predator$ED_end <- 5000

print(bio_chinook)

## ----plot-setup, fig.cap="Model setup dashboard: environmental and diet data coverage."----
plot(bio_chinook, type = "dashboard")

## ----plot-temperature, fig.cap="Seasonal temperature profile used in the simulation."----
plot(bio_chinook, type = "temperature", colors = "red")

## ----plot-diet, fig.cap="Daily diet composition over the 180-day season."-----
plot(bio_chinook, type = "diet", colors = "green")

## ----binary-search------------------------------------------------------------
res_bs <- run_fb4(
  x         = bio_chinook,
  fit_to    = "Weight",
  fit_value = 40,
  strategy  = "binary_search",
  verbose   = FALSE
)

cat(sprintf("Estimated p-value     : %.4f\n", res_bs$summary$p_value))
cat(sprintf("Final weight          : %.1f g\n", res_bs$summary$final_weight))
cat(sprintf("Total consumption     : %.1f g\n", res_bs$summary$total_consumption_g))
cat(sprintf("Simulation converged  : %s\n",     res_bs$summary$converged))

## ----plot-growth, fig.cap="Daily weight trajectory from binary search fit."----
plot(res_bs, type = "growth")

## ----plot-energy, fig.cap="Daily energy budget partitioning (J/day)."---------
plot(res_bs, type = "energy")

## ----plot-dashboard, fig.cap="Simulation dashboard: growth, consumption, temperature, and energy."----
plot(res_bs, type = "dashboard")

## ----direct-------------------------------------------------------------------
res_direct <- run_fb4(
  x         = bio_chinook,
  fit_to    = "p_value",
  fit_value = 0.75,
  strategy  = "direct_p_value",
  verbose   = FALSE
)

cat(sprintf("Final weight at p = 0.75 : %.1f g\n", res_direct$summary$final_weight))
cat(sprintf("Total consumption        : %.1f g\n", res_direct$summary$total_consumption_g))

## ----bootstrap, cache=TRUE----------------------------------------------------
set.seed(123)
n_obs           <- 25
final_wt_true   <- res_bs$summary$final_weight
obs_weights     <- rnorm(n_obs, mean = final_wt_true, sd = final_wt_true * 0.08)

res_boot <- run_fb4(
  x                = bio_chinook,
  fit_to           = "Weight",
  observed_weights = obs_weights,
  strategy         = "bootstrap",
  n_bootstrap      = 100,
  upper            = 1,
  parallel         = FALSE,
  confidence_level = 0.95,
  verbose          = FALSE
)

cat(sprintf("p mean (bootstrap) : %.4f\n", res_boot$summary$p_mean))
cat(sprintf("p SD               : %.4f\n", res_boot$summary$p_sd))
cat(sprintf("95%% CI             : [%.4f, %.4f]\n",
            res_boot$method_data$confidence_intervals$p_ci_lower,
            res_boot$method_data$confidence_intervals$p_ci_upper))

## ----plot-boot, cache=TRUE, fig.cap="Bootstrap distribution of estimated p-values with 95% CI."----
plot(res_boot, type = "uncertainty")

## ----analysis-----------------------------------------------------------------
growth_stats   <- analyze_growth_patterns(res_bs)
feeding_stats  <- analyze_feeding_performance(res_bs)
energy_budget  <- analyze_energy_budget(res_bs)

# Growth metrics
cat("=== Growth ===\n")
cat(sprintf("  Final weight            : %.1f g\n",
            growth_stats$final_weight$estimate))
cat(sprintf("  Total growth            : %.1f g\n",
            growth_stats$total_growth$estimate))
cat(sprintf("  Specific growth rate    : %.4f g/g/day\n",
            growth_stats$specific_growth_rate$estimate))

# Feeding performance
cat("\n=== Feeding performance ===\n")
cat(sprintf("  Total consumption       : %.1f g\n",
            feeding_stats$total_consumption$estimate))
cat(sprintf("  Gross conv. efficiency  : %.3f\n",
            feeding_stats$gross_conversion_efficiency$estimate))

