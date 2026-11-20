"""
02_zone_exposure.py: exposure per zone at each mapped return period, and AAL from the maps.

Reads:  config.yaml (hazard.*, exposure.*, zones.*), hazard maps, exposure raster, files.zones
Writes: files.zone_exposure (hybas_id, rp, exposed), files.zone_aal_maps (hybas_id, aal)

Method:
- Work on the hazard grid. Reproject the exposure raster onto it with sum-preserving
  resampling, so totals are unchanged. Rasterise the zone polygons (hybas_id) onto the
  same grid. Process in windows if memory needs it.
- For each return period in hazard.rps: exposed = exposure in cells with depth above
  hazard.depth_threshold_m, summed by zone.
- Check that exposure does not decrease with return period in any zone; print zones that
  break this (map artefacts) rather than fixing them.
- AAL from the maps per zone: integrate exposed over annual exceedance probability
  (1/rp) with the trapezoidal rule, with zero below the lowest mapped return period
  (hazard.below_lowest_rp) and the top map held beyond the highest. Also print the AAL
  that the lowest bin would add if exposure at the lowest RP extended to RP = 2, as a
  sensitivity.
- Print national exposed totals per return period and national AAL.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
