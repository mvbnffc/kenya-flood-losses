"""
01_zones.py: assign every HydroBASINS unit in Kenya to a modelled river site.

Reads:  config.yaml (event_set.*, zones.*), HydroBASINS units, event-set sites and site QC
Writes: files.zones (one row per unit), files.zones_map

Method:
- Units: HydroBASINS polygons that intersect the Kenya boundary. Keep HYBAS_ID, MAIN_BAS,
  NEXT_DOWN, UP_AREA, SUB_AREA.
- Modelled sites: event-set sites passing event_set.qc_max_zero_share and
  event_set.qc_min_mean_flow_m3s (the same filter the event set used).
- Each modelled site sits in one unit (point in polygon).
- Rule "own": the unit contains a modelled site (if more than one, the one with the
  largest upstream area).
- Rule "basin": otherwise, follow NEXT_DOWN downstream to the first unit with a modelled
  site, within the same MAIN_BAS and within zones.max_site_distance_km. This is the
  "small tributary borrows its river's return period" rule.
- Rule "upstream": otherwise, if the unit lies downstream of one or more modelled sites in
  the same MAIN_BAS (on their NEXT_DOWN path), take the one with the largest upstream area,
  within zones.max_site_distance_km. Covers reaches below the last site (lower Tana, Nzoia).
- Rule "fallback": no site found that way. Assign the nearest modelled site in the same
  region, and flag it. A unit's region is the region of the nearest modelled site in its
  MAIN_BAS; if its MAIN_BAS has no modelled site, the region of the nearest modelled site.
  Site regions come from the event-set config (simulation.regions: outlet_id -> region).
  03 also runs these units as independent draws, to bracket the effect.
- Columns: hybas_id, site_id, rule, distance_km, unit_area_km2, up_area_km2.
- Print counts and areas by rule. Map units coloured by rule, with site points.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
