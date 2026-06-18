# Algorithm Summary

A concise explanation of how `weather_daily_avg` turns hourly GFS-Wave/WAVEWATCH III snapshots
into daily-averaged grids. For build instructions, CLI usage, and exhaustive derivations with
worked examples, see [`README.md`](README.md). For a specific run's file manifest, see
[`average_weather_description.md`](average_weather_description.md).

---

## 1. What it does

For each calendar day, the pipeline reduces up to 8 hourly snapshots (00/03/06/09/12/15/18/21
UTC) into one representative daily grid, per field, per grid cell. It does **not** use the
arithmetic mean. Added resistance from waves and aerodynamic drag from wind both scale with the
*square* of height/speed (`R ∝ H²`, `F ∝ U²`), so a day with one calm half and one rough half is
physically worse than a day that's uniformly "average" — averaging the raw values would
understate that. Every "average" in this pipeline is instead an **energy-weighted average**,
chosen so the representative daily value preserves the total resistance-impulse the ship would
actually experience over the day, not just the mean magnitude.

## 2. Grids

The input archive is **not one uniform grid** — it's two:

| Grid | Fields | Shape (rows × cols) | Resolution |
|---|---|---|---|
| Wave grid | `sigwh`, `wsh`, `wsp`, `wsd`, `pwh`, `pwd` | 621 × 1440 | 0.25° |
| Wind grid | `was`, `wad` | 721 × 1440 | 0.25° |

Both share the same column convention: `lon = col * 0.25`, 0° to 359.75°E.

### The row 0 = south convention

This is the part most likely to surprise someone integrating with this pipeline: **row 0 is the
south end of both grids, and latitude increases northward as the row index increases.** Most
gridded earth-science datasets do the opposite (row 0 = north pole, decreasing southward), and
that was the originally-documented assumption here too — it was wrong, for both grids, and had
to be corrected after the data didn't match it.

```
wind grid:  lat = -90.0 + row * 0.25     row ∈ [0, 720]   ->  lat ∈ [-90, 90]  (full pole-to-pole)
wave grid:  lat = -75.0 + row * 0.25     row ∈ [0, 620]   ->  lat ∈ [-75, 80]  (poles truncated)
```

The wave grid has 100 fewer rows than the wind grid — it's missing the southernmost ~15° and
northernmost ~10° of latitude, consistent with WAVEWATCH III excluding regions of persistent
polar ice where open-water waves don't form. The wind grid has exactly the row count needed for
a full pole-to-pole span at 0.25° resolution (`(180/0.25)+1 = 721`), so it isn't truncated.

Neither `.npy` file carries coordinate metadata, so this had to be determined empirically rather
than assumed:
- **Wave grid:** correlating the real NaN (land/ice) pattern in `sigwh` against rasterized
  coastlines, across multiple days and every plausible row offset/order. A single sharp peak
  (~92.6% agreement) identified the south-first, north-increasing convention; verified visually
  (continents align with coastlines exactly at that offset).
- **Wind grid:** has no NaN/land mask to correlate against, so it was confirmed two other ways —
  visual landmark displacement matching the wave grid's, and independently, winter-storm wind
  maxima in `was` landing at the latitude a real North Pacific/Atlantic cyclone should occupy
  only under the south-first convention (an implausible Southern Ocean latitude under the
  opposite one).

If you're geo-referencing either grid directly, use the formulas above — don't assume the
textbook "row 0 = north" convention.

### Direction encoding

`wsd`, `wad`, and `pwd` are stored on disk as `int8` 16-point-compass codes (22.5° steps), not
raw degrees: `wsd`/`pwd` are 0-indexed (codes 0–15, code 0 = North), `wad` is 1-indexed (codes
1–16, code 1 = North). These are decoded to degrees before any averaging.

## 3. The per-day averaging math

Let `N_valid` be the count of non-NaN hourly values at a cell on a given day. A cell needs at
least 2 valid hourly timesteps to produce a non-NaN output (configurable via
`--min-valid-hours`); NaN always means land/ice or insufficient data, never a real zero.

**Height and speed fields** (`sigwh`, `wsh`, `was`) use the root-mean-square (energy) average:

```
H_eff = sqrt( (1 / N_valid) * sum_t( H_t^2 ) )      -- over valid hours only
```

This is the representative constant value whose total `H²`-impulse over the day equals the
actual sum of `H_t²` — i.e. it reproduces the same total resistance/drag the real, varying day
would have produced.

**Period** (`wsp`) uses an energy-flux-weighted mean, weighted by `wsh²` (its own paired wind-sea
height, not the combined `sigwh`):

```
T_eff = sum_t( wsh_t^2 * wsp_t ) / sum_t( wsh_t^2 )
```

A calm hour (`wsh = 0`) contributes zero weight rather than dragging the average toward whatever
arbitrary period value happened to be recorded during that calm — period is only meaningful when
there's energy behind it.

**Direction fields** (`wsd`, `wad`, `pwd`) use an energy-weighted circular mean (plain
angle-averaging breaks down on a circular quantity — the mean of 359° and 1° must be ~0°, not
180°). Each is weighted by its own physically-paired intensity field, not a generic combined
proxy:

```
u = sum_t( W_t^2 * sin(theta_t) )
v = sum_t( W_t^2 * cos(theta_t) )
theta_eff = atan2(u, v) mod 360
```

| Direction field | Weight `W` |
|---|---|
| `wsd` | `wsh` (wind-sea height) |
| `wad` | `was` (wind speed) |
| `pwd` | `pwh` (peak wave height — loaded only for this; not itself an output) |

**Swell residual** (`swell_residual`) isn't averaged independently — it's derived from the
already-daily-averaged `sigwh`/`wsh`, using spectral energy superposition:

```
H_swell = sqrt( max( sigwh_daily^2 - wsh_daily^2, 0 ) )
```

The `max(..., 0)` clamps floating-point noise (independent model runs can produce `wsh` very
slightly larger than `sigwh`) to zero rather than letting it produce a spurious NaN from a
negative square root.

## 4. Output

8 fields per day, all written as `float64` regardless of input dtype: `sigwh`, `wsh`, `wsp`,
`wsd`, `was`, `wad`, `pwd`, `swell_residual`. See README.md Section 3 for units and grid shape
per field, and Section 4 for full worked-example derivations.
