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

## ----load-db------------------------------------------------------------------
data(fish4_parameters)

n_entries <- length(fish4_parameters)
cat("Total entries in database :", n_entries, "\n")
cat("First 12 entries:\n")
print(head(names(fish4_parameters), 12))

## ----db-structure-------------------------------------------------------------
# Inspect a single entry
chinook <- fish4_parameters[["Oncorhynchus tshawytscha"]]

cat("=== Species info ===\n")
print(chinook$species_info)

cat("\n=== Life stages available ===\n")
print(names(chinook$life_stages))

cat("\n=== Sources ===\n")
print(chinook$sources)

## ----search-name--------------------------------------------------------------
# Find all Oncorhynchus species
onco_keys <- grep("Oncorhynchus", names(fish4_parameters),
                  value = TRUE, ignore.case = TRUE)
cat("Oncorhynchus species in database:\n")
print(onco_keys)

## ----search-common------------------------------------------------------------
# Extract all common names and search
common_names <- sapply(fish4_parameters, function(x) {
  x$species_info$common_name %||% NA_character_
})

trout_idx <- grep("trout|salmon|char", common_names,
                  value = FALSE, ignore.case = TRUE)
cat("Salmonid entries (trout / salmon / char):\n")
print(common_names[trout_idx])

## ----life-stages--------------------------------------------------------------
# Summary of life stages per species
stage_summary <- sapply(fish4_parameters, function(x) {
  paste(names(x$life_stages), collapse = ", ")
})

# Show species that have more than one life stage
multi_stage <- stage_summary[sapply(strsplit(stage_summary, ", "), length) > 1]
cat("Species with multiple life stages (first 8):\n")
print(head(multi_stage, 8))

## ----compare-params-----------------------------------------------------------
# Extract CA, CB, CEQ for all entries (first life stage)
param_table <- do.call(rbind, lapply(names(fish4_parameters), function(sp) {
  entry <- fish4_parameters[[sp]]
  stage <- names(entry$life_stages)[1]
  cons  <- entry$life_stages[[stage]]$consumption

  data.frame(
    Species    = sp,
    Stage      = stage,
    CEQ        = cons$CEQ  %||% NA,
    CA         = cons$CA   %||% NA,
    CB         = cons$CB   %||% NA,
    CTO        = cons$CTO  %||% NA,
    CTM        = cons$CTM  %||% NA,
    stringsAsFactors = FALSE
  )
}))

# Show top 10 by CA (highest maximum ration)
param_table_sorted <- param_table[order(-param_table$CA, na.last = TRUE), ]
knitr::kable(
  head(param_table_sorted[, c("Species", "Stage", "CEQ", "CA", "CB", "CTO", "CTM")], 10),
  caption = "Top 10 entries by CA (maximum consumption coefficient)",
  digits  = 4
)

## ----temp-optima-plot, fig.cap="Distribution of consumption temperature optima (CTO) across all database entries."----
cto_vals <- na.omit(param_table$CTO)

hist(cto_vals,
     breaks = 15,
     main   = "Consumption temperature optima (CTO)",
     xlab   = "CTO (°C)",
     ylab   = "Number of entries",
     col    = "steelblue",
     border = "white")
abline(v = median(cto_vals), col = "tomato", lwd = 2, lty = 2)
legend("topright", legend = paste("Median =", round(median(cto_vals), 1), "°C"),
       col = "tomato", lty = 2, lwd = 2, bty = "n")

## ----prededeq-dist------------------------------------------------------------
prededeq_vals <- sapply(names(fish4_parameters), function(sp) {
  entry <- fish4_parameters[[sp]]
  stage <- names(entry$life_stages)[1]
  entry$life_stages[[stage]]$predator$PREDEDEQ %||% NA
})

cat("Distribution of PREDEDEQ across database:\n")
print(table(prededeq_vals, useNA = "ifany"))

## ----from-db------------------------------------------------------------------
# Retrieve brown trout if available, otherwise fall back to first salmonid
target_sp <- grep("Salmo trutta|brown trout",
                  names(fish4_parameters),
                  value = TRUE, ignore.case = TRUE)

if (length(target_sp) == 0) {
  target_sp <- onco_keys[1]   # Fall back to first Oncorhynchus
}
target_sp <- target_sp[1]
cat("Using species:", target_sp, "\n")

entry  <- fish4_parameters[[target_sp]]
stage  <- names(entry$life_stages)[1]
params <- entry$life_stages[[stage]]
info   <- entry$species_info
info$life_stage <- stage

# Minimal 60-day setup
days_short <- 1:60
bio_db <- Bioenergetic(
  species_params     = params,
  species_info       = info,
  environmental_data = list(
    temperature = data.frame(Day = days_short,
                             Temperature = 10 + 3 * sin(pi * days_short / 60))
  ),
  diet_data = list(
    proportions = data.frame(Day = days_short, Invertebrates = 1),
    prey_names  = "Invertebrates",
    energies    = data.frame(Day = days_short, Invertebrates = 2800)
  ),
  simulation_settings = list(initial_weight = 50, duration = 60)
)

# Supply energy density bounds for PREDEDEQ = 1
bio_db$species_params$predator$ED_ini <- 4000
bio_db$species_params$predator$ED_end <- 4500

res_db <- run_fb4(bio_db,
                  fit_to    = "Weight",
                  fit_value = 80,
                  strategy  = "binary_search",
                  verbose   = FALSE)

cat(sprintf("Estimated p : %.4f  |  Final weight : %.1f g\n",
            res_db$summary$p_value,
            res_db$summary$final_weight))

