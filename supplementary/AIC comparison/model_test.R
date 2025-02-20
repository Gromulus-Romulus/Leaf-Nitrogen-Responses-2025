# Following numerical experimentation
# tested to find best fitting empirical model
# using AIC metrics.
#
# Load best model identified in spreadsheet.
#
# @author Nathan D. Malamud
# @date 2025-02-01
# 

# Libraries ----

# Import model from RDS ----

# Best working model based on AIC metrics with random slopes
#   - grt^3 ~ 1/LDMC + sqrt(LMA) + 1/CHL + (1/LDMC + sqrt(LMA) + 1/CHL | species)
#
# Caution: this is an empirical approach based on AIC metrics. Consultation with literature
# is needed in order to assess biological significance of variable transforms.
#

