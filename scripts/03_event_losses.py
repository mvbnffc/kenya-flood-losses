"""
03_event_losses.py: exposure for every simulated event and year.

Reads:  config.yaml (event_set.*, events.*), event_rp and events_sim (event set),
        files.zones, files.zone_exposure
Writes: files.event_losses (mode, event_id, sim_year, exposed),
        files.year_losses (mode, sim_year, exposed, n_events), including zero years

Method:
- For each event, each zone takes its site's rp_years. Sites missing from event_rp for that
  event are below events.rp_floor (2) and give zero exposure.
- Exposure at a zone's return period: interpolate the zone's exposure-RP curve
  (events.interpolation, linear in log RP), zero below the lowest mapped RP, top map held
  above the highest.
- Dependence modes (events.dependence_modes):
  - full: every zone uses its assigned site's return period.
  - independent_fallback: "fallback" zones draw their own return period in each event,
    p ~ Uniform(0, 1), rp = 1 / (1 - exp(-lambda * (1 - p))), seeded with events.seed;
    all other zones as in full.
- Year totals sum events per sim_year; years with no events are included with zero,
  for all event_set.n_years years.
- Print events, years, mean annual exposure per mode.
"""
import yaml

cfg = yaml.safe_load(open("config.yaml"))

# To be written in an agent session.
