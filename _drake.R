# This file serves the r_*() functions (e.g. r_make()) documented at
# https://ropenscilabs.github.io/drake-manual/projects.html#safer-interactivity # nolint
# and
# https://docs.ropensci.org/drake/reference/r_make.html

# Load your packages and supporting functions into your session.
# If you use supporting scripts like the ones below,
# you will need to supply them yourself. Examples:
# https://github.com/wlandau/drake-examples/tree/master/main/R
source("R/packages.R")  # Loads packages, e.g. library(drake).
source("R/analyses.R")  # Defines custom functions for analyses.
source("R/figures.R")   # Defines custom functions for figures.
source("R/tables.R")    # Defines custom functions for tables.
source("R/plan.R")      # Creates drake plan.

# _drake.R must end with a call to drake_config().
# The arguments to drake_config() are basically the same as those to make().
drake::drake_config(plan, verbose = 1L, lock_envir = FALSE)
