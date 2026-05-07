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

## ----base-model---------------------------------------------------------------
# Species parameters (Stewart & Ibarra 1991)
sp_params_chinook <- list(
  consumption = list(
    CEQ = 2, CA = 0.303, CB = -0.275,
    CQ = 3.0, CTO = 15.0, CTM = 20.0,
    CTL = 24.0, CK1 = 0.1, CK4 = 0.13
  ),
  respiration = list(
    REQ = 1, RA = 0.00264, RB = -0.217,
    RQ = 0.06818, RTO = 0.0234, RTM = 0.0,
    RTL = 25.0, RK1 = 1.0, RK4 = 0.13, RK5 = 0.0
  ),
  activity  = list(ACT = 1.0, BACT = 0.0),
  sda       = list(SDA = 0.172),
  egestion  = list(EGEQ = 1, FA = 0.212, FB = -0.222, FG = 0.631),
  excretion = list(EXEQ = 1, UA = 0.0314, UB = 0.58, UG = -0.299),
  predator  = list(
    PREDEDEQ = 2,
    Alpha1 = 4500, Beta1 = 6.0,
    Alpha2 = 5500, Beta2 = 2.0,
    Cutoff = 250
  )
)

days <- 1:120

# Baseline temperature: cool Pacific NW summer
temp_base <- data.frame(
  Day         = days,
  Temperature = round(7 + 6 * sin(pi * (days - 20) / 120), 2)
)

# Diet
diet_props  <- data.frame(Day = days, Alewife = 0.65, Shrimp = 0.35)
prey_energy <- data.frame(Day = days, Alewife = 4900, Shrimp = 3200)

# Build baseline object
make_bio <- function(temp_df, sp_params, initial_weight = 5) {
  bio <- Bioenergetic(
    species_params     = sp_params,
    species_info       = list(scientific_name = "Oncorhynchus tshawytscha",
                               common_name = "Chinook salmon",
                               life_stage  = "juvenile"),
    environmental_data = list(temperature = temp_df),
    diet_data          = list(
      proportions = diet_props,
      prey_names  = c("Alewife", "Shrimp"),
      energies    = prey_energy
    ),
    simulation_settings = list(initial_weight = initial_weight, duration = 120)
  )
  bio
}

bio_base <- make_bio(temp_base, sp_params_chinook)

cat(sprintf("Baseline temperature range : %.1f – %.1f °C\n",
            min(temp_base$Temperature), max(temp_base$Temperature)))

## ----manual-sensitivity-------------------------------------------------------
offsets   <- c(-3, -2, -1, 0, 1, 2, 3, 4)
target_wt <- 40

results_sens <- lapply(offsets, function(delta) {
  temp_shifted <- temp_base
  temp_shifted$Temperature <- pmax(1, temp_base$Temperature + delta)

  bio_s <- make_bio(temp_shifted, sp_params_chinook)
  res   <- tryCatch(
    run_fb4(bio_s,
            fit_to    = "Weight",
            fit_value = target_wt,
            strategy  = "binary_search",
            verbose   = FALSE),
    error = function(e) NULL
  )

  if (is.null(res)) {
    return(data.frame(offset = delta, p_value = NA, final_weight = NA))
  }
  data.frame(
    offset       = delta,
    p_value      = round(res$summary$p_value,       4),
    final_weight = round(res$summary$final_weight,   1)
  )
})

sens_df <- do.call(rbind, results_sens)

knitr::kable(
  sens_df,
  caption    = paste0("p-value required to achieve ", target_wt,
                      " g across temperature offsets"),
  col.names  = c("Temp. offset (°C)", "Estimated p", "Final weight (g)")
)

## ----plot-sensitivity-manual, fig.cap="Estimated p-value as a function of temperature offset."----
plot(sens_df$offset, sens_df$p_value,
     type = "b", pch = 19, lwd = 2, col = "steelblue",
     xlab = "Temperature offset (°C)",
     ylab = "Estimated p-value",
     main = "Feeding level required to achieve target weight")
abline(v = 0, lty = 2, col = "gray50")
grid()

## ----sensitivity-plot, fig.cap="Built-in temperature sensitivity plot. Each point shows the p required to achieve the target weight at a given temperature offset."----
plot(bio_base,
     type      = "sensitivity",
     temp_offsets = c(-3, -2, -1, 0, 1, 2, 3, 4),
     fit_to    = "Weight",
     fit_value = 40,
     colors    = "blue")

## ----climate-scenarios--------------------------------------------------------
p_fixed   <- 0.65
scenarios <- c("Baseline" = 0, "+1 °C" = 1, "+2 °C" = 2, "+4 °C" = 4)

scenario_results <- lapply(names(scenarios), function(sc) {
  delta <- scenarios[[sc]]
  temp_s <- temp_base
  temp_s$Temperature <- pmax(1, temp_base$Temperature + delta)

  bio_s <- make_bio(temp_s, sp_params_chinook)
  res   <- run_fb4(bio_s,
                   fit_to    = "p_value",
                   fit_value = p_fixed,
                   strategy  = "direct_p_value",
                   verbose   = FALSE)

  if (is.null(res)) return(NULL)

  data.frame(
    Scenario    = sc,
    Delta_T     = delta,
    Final_wt_g  = round(res$summary$final_weight,       1),
    Consumption = round(res$summary$total_consumption_g, 1)
  )
})

scenario_df <- do.call(rbind, scenario_results)

knitr::kable(
  scenario_df,
  caption   = paste0("Projected outcomes at p = ", p_fixed,
                     " under four temperature scenarios"),
  col.names = c("Scenario", "Δ T (°C)", "Final weight (g)", "Total consumption (g)")
)

## ----plot-climate, fig.cap="Projected final weight under climate change scenarios (p fixed at 0.65)."----
bar_cols <- c("steelblue", "gold", "darkorange", "firebrick")
bp <- barplot(scenario_df$Final_wt_g,
              names.arg = scenario_df$Scenario,
              col       = bar_cols,
              ylab      = "Final weight (g)",
              main      = paste("Growth at p =", p_fixed,
                                "under warming scenarios"),
              ylim      = c(0, max(scenario_df$Final_wt_g) * 1.25),
              border    = NA)

text(bp, scenario_df$Final_wt_g + 0.5,
     labels = paste0(scenario_df$Final_wt_g, " g"),
     cex = 0.85, font = 2)
grid(nx = NA, ny = NULL, col = "gray85")

## ----thermal-window-----------------------------------------------------------
# Sweep a fine temperature grid and record required p
temp_grid   <- seq(2, 24, by = 1)
p_at_target <- sapply(temp_grid, function(t_const) {
  temp_const <- data.frame(Day = days, Temperature = t_const)
  bio_t <- make_bio(temp_const, sp_params_chinook)
  tryCatch({
    res <- run_fb4(bio_t,
                   fit_to    = "Weight",
                   fit_value = target_wt,
                   strategy  = "binary_search",
                   verbose   = FALSE)
    res$summary$p_value
  }, error = function(e) NA_real_)
})

perf_df <- data.frame(Temperature = temp_grid, p_required = p_at_target)
feasible <- perf_df[!is.na(perf_df$p_required) & perf_df$p_required <= 1, ]

cat(sprintf("Feasible temperature range for %.0f g in %d days : %.0f – %.0f °C\n",
            target_wt, max(days),
            min(feasible$Temperature), max(feasible$Temperature)))

## ----plot-thermal-window, fig.cap="Thermal performance window: p required to achieve 40 g over 120 days at constant temperatures. The shaded region shows the feasible window (p ≤ 1)."----
plot(perf_df$Temperature, perf_df$p_required,
     type = "l", lwd = 2, col = "steelblue",
     xlab = "Constant water temperature (°C)",
     ylab = "p required to achieve 40 g",
     main = "Thermal performance window — juvenile Chinook",
     ylim = c(0, max(perf_df$p_required, na.rm = TRUE) * 1.1))

# Shade feasible zone
polygon(c(min(feasible$Temperature), feasible$Temperature,
          max(feasible$Temperature)),
        c(0, feasible$p_required, 0),
        col = adjustcolor("steelblue", alpha.f = 0.2),
        border = NA)

abline(h = 1, lty = 2, col = "firebrick", lwd = 1.5)
text(max(temp_grid) - 2, 1.02, "p = 1 (maximum ration)", col = "firebrick", cex = 0.8)
grid()

