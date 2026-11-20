"""
01_zones.py: assign every HydroBASINS unit in Kenya to a modelled river site.

Reads:  config.yaml (event_set.*, zones.*), HydroBASINS units, event-set sites and site QC
Writes: files.zones (one row per unit), files.zones_map, files.site_placement (one row per site)

Method:
- Units: HydroBASINS polygons that intersect the Kenya boundary. Keep HYBAS_ID, MAIN_BAS,
  NEXT_DOWN, UP_AREA, SUB_AREA.
- Modelled sites: event-set sites passing event_set.qc_max_zero_share and
  event_set.qc_min_mean_flow_m3s (the same filter the event set used).
- Each modelled site sits in one unit: point in polygon, unless site upstream area / UP_AREA is
  outside zones.up_area_ratio_bounds; then the unit within zones.site_snap_km (km to degrees at
  the equator) with the smallest |log(site upstream area / UP_AREA)|.
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
- Distances: straight-line (haversine) from the unit centroid to the site.
- Columns: hybas_id, main_bas, site_id, region, rule, distance_km, unit_area_km2, up_area_km2.
- Print counts, areas and distance_km by rule; the modelled sites that lose "own" to a
  larger site in the same unit (unused downstream); and per site whether it was snapped, the
  distance moved and the final upstream area / UP_AREA, listing those still outside the bounds.
- Map units coloured by rule, with the Kenya boundary and site points.
"""
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

cfg = yaml.safe_load(open("config.yaml"))
es, zc = cfg["event_set"], cfg["zones"]
max_km = zc["max_site_distance_km"]
RULES = ["own", "basin", "upstream", "fallback"]
R_KM = 6371.0
km_per_deg = 2 * np.pi * R_KM / 360  # at the equator


def haversine_km(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, (lon1, lat1, lon2, lat2))
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * R_KM * np.arcsin(np.sqrt(a))


def downstream(h):
    """Units below unit h on its NEXT_DOWN path, while the path stays in the Kenya unit set."""
    path = []
    h = next_down[h]
    while h in next_down:
        path.append(h)
        h = next_down[h]
    return path


# Units: bbox filter on read (the Africa file is large), then keep those touching Kenya
bnd = gpd.read_file(zc["kenya_boundary"])
units = gpd.read_file(zc["hydrobasins"], bbox=tuple(bnd.total_bounds))
assert bnd.crs.equals(units.crs, ignore_axis_order=True), (bnd.crs, units.crs)
units = units[units.intersects(bnd.union_all())].reset_index(drop=True)
units = units[["HYBAS_ID", "MAIN_BAS", "NEXT_DOWN", "UP_AREA", "SUB_AREA", "geometry"]]
next_down = dict(zip(units.HYBAS_ID, units.NEXT_DOWN))
print(f"{len(units)} HydroBASINS units intersect Kenya")

# Modelled sites, with region and containing unit
regions = yaml.safe_load(open(es["config"]))["simulation"]["regions"]
qc = pd.read_parquet(es["site_qc"])[["site_id", "zero_share", "mean_flow_m3s"]]
sites = pd.read_parquet(es["sites"]).merge(qc, on="site_id", validate="1:1")
sites = sites[(sites.zero_share <= es["qc_max_zero_share"]) & (sites.mean_flow_m3s >= es["qc_min_mean_flow_m3s"])]
sites["region"] = sites.outlet_id.map(regions)
assert sites.region.notna().all(), sites[sites.region.isna()]
pts = gpd.GeoDataFrame(sites, geometry=gpd.points_from_xy(sites.lon, sites.lat), crs=units.crs)
pts = gpd.sjoin(pts, units[["HYBAS_ID", "geometry"]], predicate="intersects", how="left")
assert pts.HYBAS_ID.notna().all(), pts[pts.HYBAS_ID.isna()][["site_id", "name", "lon", "lat"]]
assert pts.site_id.is_unique, pts[pts.site_id.duplicated(keep=False)][["site_id", "HYBAS_ID"]]
pts = pts.drop(columns="index_right").astype({"HYBAS_ID": "int64"}).reset_index(drop=True)
print(f"{len(pts)} modelled sites")

# Snap sites whose unit's UP_AREA (measured at the unit outlet) is far from the site's
# upstream area to the best-matching unit within site_snap_km
lo, hi = zc["up_area_ratio_bounds"]
pts["ratio_pip"] = pts.upstream_area_km2 / pts.HYBAS_ID.map(dict(zip(units.HYBAS_ID, units.UP_AREA)))
pts["snapped"], pts["moved_km"] = False, 0.0
bad = np.flatnonzero((pts.ratio_pip < lo) | (pts.ratio_pip > hi))
qi, ui = units.sindex.query(pts.geometry.values[bad], predicate="dwithin", distance=zc["site_snap_km"] / km_per_deg)
for k, i in enumerate(bad):
    cand = ui[qi == k]
    j = cand[np.argmin(np.abs(np.log(pts.upstream_area_km2[i] / units.UP_AREA.values[cand])))]
    if units.HYBAS_ID[j] != pts.HYBAS_ID[i]:
        pts.loc[i, ["HYBAS_ID", "snapped"]] = units.HYBAS_ID[j], True
        pts.loc[i, "moved_km"] = units.geometry[j].distance(pts.geometry[i]) * km_per_deg
pts = pts.merge(units[["HYBAS_ID", "MAIN_BAS", "UP_AREA"]], on="HYBAS_ID", validate="m:1")
pts["ratio"] = pts.upstream_area_km2 / pts.UP_AREA
cols = ["site_id", "name", "lon", "lat", "region", "HYBAS_ID", "upstream_area_km2", "UP_AREA",
        "ratio_pip", "snapped", "moved_km", "ratio"]
Path(cfg["files"]["site_placement"]).parent.mkdir(parents=True, exist_ok=True)
pts[["site_id", "snapped", "moved_km", "ratio_pip", "ratio", "HYBAS_ID"]].rename(
    columns={"ratio": "ratio_final", "HYBAS_ID": "hybas_id"}).to_csv(cfg["files"]["site_placement"], index=False)
print(f"\n{len(bad)} sites outside {lo}-{hi} by point in polygon; {pts.snapped.sum()} snapped.")
print("Per site (ratio_pip: point in polygon; ratio: final):")
print(pts.sort_values("ratio")[cols].to_string(index=False, float_format="%.2f"))
out = pts[(pts.ratio < lo) | (pts.ratio > hi)]
print(f"\n{len(out)} sites still outside {lo}-{hi} after snapping:")
print(out.sort_values("ratio")[cols].to_string(index=False, float_format="%.2f"))

# Rule "own": one site per unit, the largest upstream area; the rest are dropped
pts = pts.sort_values("upstream_area_km2", ascending=False)
win = pts.drop_duplicates("HYBAS_ID").reset_index(drop=True)
lost = pts[pts.HYBAS_ID.duplicated()].merge(
    win[["HYBAS_ID", "site_id", "upstream_area_km2"]], on="HYBAS_ID", suffixes=("", "_kept"))
print(f"\n{len(lost)} modelled sites lose 'own' to a larger site in the same unit (unused downstream):")
print(lost[["site_id", "name", "upstream_area_km2", "HYBAS_ID", "site_id_kept", "upstream_area_km2_kept"]]
      .to_string(index=False, float_format="%.0f"))

cent = units.centroid  # geographic-CRS centroids; adequate for level-10 units
D = haversine_km(cent.x.values[:, None], cent.y.values[:, None], win.lon.values[None, :], win.lat.values[None, :])
site_in = dict(zip(win.HYBAS_ID, win.index))
n = len(units)
site_j = np.full(n, -1)
rule = np.full(n, "", dtype=object)
for h, j in site_in.items():
    i = units.index[units.HYBAS_ID == h][0]
    site_j[i], rule[i] = j, "own"

# Rule "basin": first site downstream, if within max distance
for i, h in enumerate(units.HYBAS_ID):
    if rule[i]:
        continue
    first = next((d for d in downstream(h) if d in site_in), None)
    if first is not None and D[i, site_in[first]] <= max_km:
        site_j[i], rule[i] = site_in[first], "basin"

# Rule "upstream": among sites whose downstream path passes through the unit, within max
# distance, the largest upstream area
above = {}
for h, j in site_in.items():
    for d in downstream(h):
        above.setdefault(d, []).append(j)
for i, h in enumerate(units.HYBAS_ID):
    if rule[i]:
        continue
    near = [j for j in above.get(h, []) if D[i, j] <= max_km]
    if near:
        site_j[i], rule[i] = max(near, key=lambda j: win.upstream_area_km2[j]), "upstream"

# Rule "fallback": nearest site in the unit's region
for i in np.flatnonzero(rule == ""):
    cand = np.flatnonzero(win.MAIN_BAS.values == units.MAIN_BAS[i])
    if len(cand) == 0:
        cand = np.arange(len(win))
    region = win.region[cand[np.argmin(D[i, cand])]]
    cand = np.flatnonzero(win.region.values == region)
    site_j[i], rule[i] = cand[np.argmin(D[i, cand])], "fallback"
assert (site_j >= 0).all()

zones = pd.DataFrame({
    "hybas_id": units.HYBAS_ID, "main_bas": units.MAIN_BAS,
    "site_id": win.site_id.values[site_j], "region": win.region.values[site_j],
    "rule": rule, "distance_km": D[np.arange(n), site_j],
    "unit_area_km2": units.SUB_AREA, "up_area_km2": units.UP_AREA,
})
Path(cfg["files"]["zones"]).parent.mkdir(parents=True, exist_ok=True)
zones.to_parquet(cfg["files"]["zones"], index=False)

summary = zones.groupby("rule").agg(units=("hybas_id", "size"), area_km2=("unit_area_km2", "sum")).reindex(RULES)
summary["area_share"] = summary.area_km2 / summary.area_km2.sum()
print("\nUnits and area by rule (unit areas are whole HydroBASINS units, including parts outside Kenya):")
print(summary.to_string(float_format="%.3f"))
print("\ndistance_km by rule:")
print(zones.groupby("rule").distance_km.describe(percentiles=[0.5, 0.9]).reindex(RULES).to_string(float_format="%.1f"))
print("\nUnits by region and rule:")
print(pd.crosstab(zones.region, zones.rule)[RULES].to_string())

# Map
fig, ax = plt.subplots(figsize=(9, 10))
gpd.GeoDataFrame(zones, geometry=units.geometry.values, crs=units.crs).plot(
    column="rule", categorical=True, categories=RULES, legend=True, ax=ax, edgecolor="white", linewidth=0.1)
bnd.boundary.plot(ax=ax, color="black", linewidth=0.8)
ax.scatter(win.lon, win.lat, s=10, c="black", zorder=3)
ax.set_title(f"Zones by rule ({len(win)} sites used, max distance {max_km} km)")
ax.set_axis_off()
Path(cfg["files"]["zones_map"]).parent.mkdir(parents=True, exist_ok=True)
fig.savefig(cfg["files"]["zones_map"], dpi=200, bbox_inches="tight")
print(f"\nWrote {cfg['files']['zones']} and {cfg['files']['zones_map']}")
