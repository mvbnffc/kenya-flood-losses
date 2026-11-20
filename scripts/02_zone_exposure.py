"""
02_zone_exposure.py: exposure per zone at each mapped return period, and AAL from the maps.

Reads:  config.yaml (hazard.*, exposure.*, zones.hydrobasins, zones.kenya_boundary, files.zones),
        hazard maps, exposure raster, HydroBASINS units, zones table
Writes: files.zone_exposure (hybas_id, rp, exposed), files.zone_aal_maps (hybas_id, aal)

Method:
- Grids: the exposure raster and the hazard maps share one grid. Fail if shape, transform or
  CRS differ. No resampling.
- Hazard: raw value times hazard.depth_to_m. Hazard nodata means dry (depth 0), not missing.
- Exposure: nodata and negative values set to 0. Print the national total.
- Zones: HydroBASINS polygons joined to files.zones on hybas_id, rasterised onto the grid.
  Count only cells inside the Kenya boundary. Print the exposure inside Kenya that falls in
  no zone; it should be about zero.
- Read whole rasters, one hazard map at a time. For each return period in hazard.rps:
  exposed = exposure in Kenya cells with depth above hazard.depth_threshold_m, summed by
  zone with np.bincount.
- Permanent water: cells with depth >= hazard.permanent_water_depth_m in the lowest
  return-period map are lakes (e.g. Turkana); exclude them at every return period and
  print the exposure excluded.
- Check that exposure does not decrease with return period in any zone; print zones that
  break this (map artefacts) rather than fixing them.
- AAL from the maps per zone: integrate exposed over annual exceedance probability (1/rp)
  with the trapezoidal rule, zero below the lowest mapped return period
  (hazard.below_lowest_rp) and the top map held beyond the highest. As a sensitivity, print
  the AAL that would be added if exposure at the lowest RP extended to AEP = 1.
- Print national exposed per return period, national AAL, and exposed at hazard.summary_rp
  and AAL by zone rule (own, basin, upstream, fallback) and by region, with shares.
"""
import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import yaml
from rasterio.features import rasterize

cfg = yaml.safe_load(open("config.yaml"))
hz = cfg["hazard"]
rps = hz["rps"]
assert rps == sorted(rps), rps
assert hz["below_lowest_rp"] == "zero", hz["below_lowest_rp"]
assert hz["summary_rp"] in rps, hz["summary_rp"]
paths = [hz["path_pattern"].format(rp=rp) for rp in rps]

# Grids
with rasterio.open(cfg["exposure"]["path"]) as r:
    grid = (r.shape, r.transform, r.crs)
    pop = r.read(1, masked=True).filled(0).astype("float64")
for p in paths:
    with rasterio.open(p) as r:
        assert (r.shape, r.transform, r.crs) == grid, (p, (r.shape, r.transform, r.crs), grid)
assert np.isfinite(pop).all(), "non-finite exposure values"
pop[pop < 0] = 0
print(f"National total exposure ({cfg['exposure']['kind']}, {cfg['exposure']['year']}): {pop.sum():,.0f}")

# Zones and Kenya mask on the grid (cell-centre rule)
zones = pd.read_parquet(cfg["files"]["zones"])
bnd = gpd.read_file(cfg["zones"]["kenya_boundary"])
units = gpd.read_file(cfg["zones"]["hydrobasins"], bbox=tuple(bnd.total_bounds), columns=["HYBAS_ID"])
units = zones[["hybas_id"]].merge(units, left_on="hybas_id", right_on="HYBAS_ID", how="left", validate="1:1")
assert units.geometry.notna().all(), units[units.geometry.isna()].hybas_id
assert bnd.crs.equals(grid[2], ignore_axis_order=True), (bnd.crs, grid[2])
nz = len(zones)
zone = rasterize(zip(units.geometry, range(1, nz + 1)), out_shape=grid[0], transform=grid[1], fill=0, dtype="int32")
kenya = rasterize(((g, 1) for g in bnd.geometry), out_shape=grid[0], transform=grid[1], fill=0, dtype="uint8")
in_kenya = (kenya == 1) & (pop > 0)
print(f"Exposure inside Kenya: {pop[in_kenya].sum():,.0f}")
print(f"Exposure inside Kenya in no zone: {pop[in_kenya & (zone == 0)].sum():,.0f}")

# Flat arrays of exposed cells in Kenya and in a zone; free the full grids
idx = np.flatnonzero(in_kenya & (zone > 0))
z = zone.ravel()[idx] - 1
w = pop.ravel()[idx]
del pop, zone, kenya, in_kenya

# Exposure per return period, one map at a time (lowest first: it defines permanent water)
exposed = np.zeros((nz, len(rps)))
for k, (rp, p) in enumerate(zip(rps, paths)):
    with rasterio.open(p) as r:
        depth = r.read(1, masked=True).filled(0).ravel()[idx] * hz["depth_to_m"]
    assert np.isfinite(depth).all(), p
    if k == 0:
        lake = depth >= hz["permanent_water_depth_m"]
        print(f"Exposure excluded as permanent water (depth >= {hz['permanent_water_depth_m']} m at RP{rp}): "
              f"{w[lake].sum():,.0f} in {lake.sum():,} cells")
        idx, z, w, depth = idx[~lake], z[~lake], w[~lake], depth[~lake]
    exposed[:, k] = np.bincount(z, weights=w * (depth > hz["depth_threshold_m"]), minlength=nz)

pd.DataFrame({
    "hybas_id": np.repeat(zones.hybas_id.values, len(rps)),
    "rp": np.tile(rps, nz),
    "exposed": exposed.ravel(),
}).to_parquet(cfg["files"]["zone_exposure"], index=False)

# Monotonicity check
drop = np.diff(exposed, axis=1)
bad = (drop < 0).any(axis=1)
print(f"\n{bad.sum()} of {nz} zones have exposure decreasing with return period somewhere.")
worst = pd.DataFrame({"hybas_id": zones.hybas_id, "rule": zones.rule, "region": zones.region,
                      "largest_drop": -drop.min(axis=1), "at_rp": np.array(rps[1:])[drop.argmin(axis=1)],
                      "n_drops": (drop < 0).sum(axis=1)})[bad]
print("Five largest drops (at_rp: the return period where exposure is lower than at the one before):")
print(worst.nlargest(5, "largest_drop").to_string(index=False, float_format="%.0f"))

# AAL from the maps
aep = 1 / np.array(rps, dtype=float)
aal = (0.5 * (exposed[:, :-1] + exposed[:, 1:]) * (aep[:-1] - aep[1:])).sum(axis=1) + exposed[:, -1] * aep[-1]
extra = exposed[:, 0] * (1 - aep[0])
pd.DataFrame({"hybas_id": zones.hybas_id, "aal": aal}).to_parquet(cfg["files"]["zone_aal_maps"], index=False)

print("\nNational exposed per return period:")
print(pd.Series(exposed.sum(axis=0), index=pd.Index(rps, name="rp"), name="exposed").to_string(float_format="{:,.0f}".format))
print(f"\nNational AAL from the maps: {aal.sum():,.0f}")
print(f"Sensitivity: holding RP{rps[0]} exposure down to AEP 1 would add {extra.sum():,.0f} "
      f"({extra.sum() / aal.sum():.0%} of AAL)")

s = zones[["rule", "region"]].assign(**{f"exposed_rp{hz['summary_rp']}": exposed[:, rps.index(hz["summary_rp"])],
                                        "aal": aal})
for by, order in [("rule", ["own", "basin", "upstream", "fallback"]), ("region", None)]:
    t = s.drop(columns="region" if by == "rule" else "rule").groupby(by).sum()
    t = t.reindex(order) if order else t.sort_values("aal", ascending=False)
    for c in list(t.columns):
        t[f"{c}_share"] = t[c] / t[c].sum()
    print(f"\nBy zone {by}:")
    print(t.to_string(formatters={c: ("{:.1%}" if c.endswith("share") else "{:,.0f}").format for c in t.columns}))
