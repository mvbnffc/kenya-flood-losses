"""
06_event_maps.py: depth footprints for a few selected simulated events.

Reads:  config.yaml (event_maps.*, hazard.*, zones.kenya_boundary, event_set.event_rp,
        event_set.events_sim, files.zones, files.event_losses), hazard maps, zones,
        event_rp, events_sim, event losses from 03
Writes: one depth GeoTIFF per event (event_maps.out_dir), one PNG per event
        (event_maps.figures_dir), event_maps.summary (one row per event)

Purpose: illustrative scenario maps for reports and workshops. No other script reads these.

Method:
- Events: the event_ids in event_maps.events. Fail if the list is empty or an id is unknown.
- Rasterise the zones (HYBAS_ID) and the Kenya boundary onto the hazard grid once. Cells
  outside a zone or outside Kenya are nodata.
- Each zone takes its site's rp_years for the event from event_rp. A zone whose site has no
  row for the event (below the event-set RP floor) is dry.
- Depth per cell: linear in log(return period) between the two mapped return periods
  (hazard.rps) that bracket the zone's event RP. Below the lowest mapped RP: dry. Above the
  highest: the highest map. Hazard nodata means dry. Multiply raw values by hazard.depth_to_m.
- Hold at most two hazard maps in memory: loop over consecutive pairs in hazard.rps and fill
  the cells of the zones whose event RP falls in that interval.
- Per event, write a row: event_id, simulated year, zones wet, flooded area in km2
  (depth > hazard.depth_threshold_m; cell area from latitude), people flooded from the
  footprint, and the same event's exposure from 03 (full dependence). Print both exposure
  figures.
- PNG: depth in metres over a grey Kenya outline, thin zone outlines, title with event_id,
  simulated year and people flooded. Caption: "Illustrative scenario. Depths can jump at
  zone boundaries."
- GeoTIFFs: float32, metres, compressed, nodata set.

Licence: depth GeoTIFFs built from Fathom stay internal until the licence terms for depth
footprints are confirmed.
"""