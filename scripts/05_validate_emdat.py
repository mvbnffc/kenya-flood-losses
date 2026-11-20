"""
05_validate_emdat.py: compare simulated people affected with EM-DAT.

Reads:  config.yaml (validation.*), files.event_losses, files.year_losses, EM-DAT export
Writes: figures/validation_emdat.png, and prints the comparison

Method (Olcese et al. 2024):
- EM-DAT: Kenya floods in validation.emdat_years, people affected per event and per year.
  Note: EM-DAT is incomplete for small events, and affected is not the same as exposed;
  treat the comparison as order-of-magnitude.
- Draw validation.n_samples windows of the same number of years from the simulated years.
  In each, keep the largest N events, where N is the number of EM-DAT events, and build
  an annual exceedance curve.
- Plot the spread of simulated curves (e.g. 5-95%) with EM-DAT points, for mode full.
- Print where the EM-DAT years fall within the simulated spread.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
