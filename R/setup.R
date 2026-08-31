# ---------------------------------------------------------------------------
# setup.R
# Single entry point sourced by every analysis script:
#   source("R/setup.R")
# It loads the configuration and all shared helpers, in dependency order.
# ---------------------------------------------------------------------------

source("R/config.R")
source("R/utils.R")
source("R/plots.R")
source("R/plot_volcano.R")
