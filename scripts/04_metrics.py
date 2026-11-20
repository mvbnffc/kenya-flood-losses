"""
04_metrics.py: headline metrics and consistency checks.

Reads:  config.yaml, files.zone_aal_maps, files.event_losses, files.year_losses, files.zones
Writes: files.metrics, figures/metrics_*.png

Method:
- AAL from the event set (mean of year totals) vs AAL from the maps (02), per mode.
  They should agree within a few percent for mode "full"; print the difference.
- Per-event exceedance curve (largest event per year, OEP) and annual exceedance curve
  (year totals, AEP) at RP 2, 5, 10, 25, 50, 100, 200, 250, 500, 1000, per mode, with the
  gap between modes.
- Coverage: share of AAL and of exposed people at RP 100 by zone rule (own / basin /
  fallback).
- Top 10 zones by AAL, with their site and rule.
- Figures: exceedance curves for both modes; map of AAL by zone.
- Print candidate events for 06: for each event_maps.oep_rps, the three simulated events
  whose national exposure (full dependence) is closest to the OEP value; and the largest
  event in each event-set region.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
