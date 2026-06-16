# Slide 1 — Why Energy-Weighted Averaging?

**Title:** The Math Behind the Average

**Bullets:**
- Wave resistance scales as R ∝ H²; wind drag scales as F ∝ U² — both nonlinear.
- A day that's half calm and half rough is worse than a uniformly "average" day.
- A plain arithmetic mean understates that; the energy (RMS) average does not.
- Example: heights [0m, 4m] average to 2.0m, but the energy average is 2.83m.

**Speaker notes:**
Because resistance and drag scale with the square of height/speed, the representative daily value needs to preserve the total H² (or U²) over the day, not just the average magnitude. The energy average — the root-mean-square, sqrt(mean(H²)) — does exactly that: it's the constant value that would produce the same total resistance-impulse as the real, varying day. By the QM-AM inequality, this RMS value is always ≥ the arithmetic mean, and equal only when conditions never vary at all. Worked example: a day with hourly heights [0m, 4m] has an arithmetic mean of 2.0m, but an energy average of sqrt((0² + 4²)/2) = 2.83m. The implied resistance ratio is (2.83/2.0)² = 2.0 — using the arithmetic mean would understate the day's resistance by 100% in this example. That error would propagate directly into the routing heuristic's edge weights, so every field in this pipeline uses the energy-weighted form instead.

---

# Slide 2 — How the Daily Average Works

**Title:** From Hourly Snapshots to One Daily Value

**Bullets:**
- Wave resistance and wind drag scale with the *square* of height/speed, not linearly.
- So each day uses an energy-weighted average, not a simple mean.
- Heights and wind speed: root-mean-square (RMS) of the hourly values.
- Wave and wind directions: a circular mean, weighted by each direction's own paired height/speed.

**Speaker notes:**
Added resistance from waves scales as R ∝ H², and wind drag scales as F ∝ U². That means a day that's half calm and half rough is physically worse than a uniformly "average" day — a plain arithmetic mean would understate that. So every field here uses an energy-weighted average instead: it's the constant value that reproduces the same total resistance-impulse the ship actually experienced. For heights and wind speed (sigwh, wsh, was), that's RMS: sqrt(mean(H²)). For period (wsp), it's weighted by wsh² so a calm hour contributes zero weight rather than skewing the period toward a meaningless reading. For direction fields (wsd, wad, pwd), plain angle-averaging breaks down (mean of 359° and 1° should be ~0°, not 180°), so we use a circular mean via atan2 of weighted sin/cos sums — and each direction is weighted by its *own* paired intensity (wsh² for wave direction, was² for wind direction, pwh² for peak-wave direction), not a generic combined energy number, to avoid one component's direction bleeding into another's.

---

# Slide 3 — What Gets Produced

**Title:** Grids and Daily Output Artifacts

**Bullets:**
- Two source grids, both 0.25° resolution: a wave grid and a wind grid.
- Row 0 is the south end for both grids — latitude increases northward.
- 8 fields produced per day: wave height, wind-sea height/period/direction, wind speed/direction, peak-wave direction, and swell residual.
- Land/ice and insufficient data stay NaN — never a fake zero.

**Speaker notes:**
The wave grid (sigwh, wsh, wsp, wsd, pwh, pwd) is 621×1440; the wind grid (was, wad) is 721×1440 — they're not the same grid, and the wave grid is missing the southernmost ~15° and northernmost ~10° (consistent with WAVEWATCH III excluding persistent polar ice). Both grids store row 0 at the south end with latitude increasing northward — the opposite of the typical "row 0 = north pole" convention, discovered empirically (by correlating land/ice patterns against real coastlines, and cross-checking storm-track latitudes) since the source files carry no coordinate metadata. The 8 daily output fields are: sigwh (combined wave height), wsh (wind-sea height), wsp (wind-sea period), wsd (wind-sea direction), was (wind speed), wad (wind direction), pwd (peak-wave direction), and swell_residual (non-local swell energy, derived from sigwh and wsh). A grid cell only gets a real value if at least 2 of the day's 8 hourly snapshots were valid there; otherwise it's NaN — and NaN always means missing/land/ice, never an actual zero reading.
