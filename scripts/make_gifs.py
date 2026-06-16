#!/usr/bin/env python3
"""Renders one animated GIF per daily-averaged output field, on a real map.

Reads the <output-dir>/<YYYY>/<MM>/<DD>/{sigwh,wsh,wsp,wsd,was,wad,pwd,swell_residual}.npy
files produced by weather_daily_avg and animates them over the available days,
with coastlines and land shading for geographic context.

Grid geo-referencing (see README.md Section 2 / src/config.hpp):
  - wave grid (sigwh, wsh, wsp, wsd, pwd, swell_residual): 621 rows, lat = -75 + row*0.25
  - wind grid (was, wad): 721 rows, lat = -90 + row*0.25
  - both grids: row 0 is the SOUTH end, latitude increases northward with row;
    1440 cols, lon = col*0.25 (0 to 359.75 E).

Usage:
    python3 scripts/make_gifs.py --output-dir test_output --out-dir gifs
"""

import argparse
import sys
from pathlib import Path

import numpy as np

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as exc:
    sys.exit(
        "Missing a plotting dependency. Install with:\n"
        "  pip install -r scripts/requirements.txt\n"
        f"Original error: {exc}"
    )

WAVE_GRID_ROWS = 621
WIND_GRID_ROWS = 721
WAVE_GRID_LAT_AT_ROW0 = -75.0  # south end; lat increases northward with row
WIND_GRID_LAT_AT_ROW0 = -90.0  # south pole; lat increases northward with row
GRID_COLS = 1440
RESOLUTION_DEG = 0.25

# field name -> (units, colormap, fixed (vmin, vmax) or None for data-driven)
FIELDS = {
    "sigwh": ("m", "viridis", None),
    "wsh": ("m", "viridis", None),
    "wsp": ("s", "viridis", None),
    "wsd": ("deg, CW from N", "hsv", (0.0, 360.0)),  # circular field
    "was": ("m/s", "viridis", None),
    "wad": ("deg, CW from N", "hsv", (0.0, 360.0)),  # circular field
    "pwd": ("deg, CW from N", "hsv", (0.0, 360.0)),  # circular field
    "swell_residual": ("m", "viridis", None),
}

WIND_GRID_FIELDS = {"was", "wad"}


def grid_coords(field: str):
    lon = np.arange(GRID_COLS) * RESOLUTION_DEG
    if field in WIND_GRID_FIELDS:
        lat = WIND_GRID_LAT_AT_ROW0 + np.arange(WIND_GRID_ROWS) * RESOLUTION_DEG
    else:
        lat = WAVE_GRID_LAT_AT_ROW0 + np.arange(WAVE_GRID_ROWS) * RESOLUTION_DEG
    return lon, lat


def discover_days(output_dir: Path):
    """Returns sorted (date_str, day_dir) pairs for every YYYY/MM/DD directory present."""
    days = []
    for year_dir in sorted(output_dir.glob("[0-9][0-9][0-9][0-9]")):
        for month_dir in sorted(year_dir.glob("[0-9][0-9]")):
            for day_dir in sorted(month_dir.glob("[0-9][0-9]")):
                date_str = f"{year_dir.name}-{month_dir.name}-{day_dir.name}"
                days.append((date_str, day_dir))
    return days


def load_frames(field: str, days):
    frames = []
    for date_str, day_dir in days:
        npy_path = day_dir / f"{field}.npy"
        if npy_path.exists():
            frames.append((date_str, np.load(npy_path)))
    return frames


def value_range(field: str, frames):
    _units, _cmap, fixed_range = FIELDS[field]
    if fixed_range is not None:
        return fixed_range
    all_vals = np.concatenate([f.ravel() for _, f in frames])
    finite = all_vals[np.isfinite(all_vals)]
    return 0.0, float(finite.max()) if finite.size else 1.0


def render_gif(field: str, days, out_path: Path, fps: float, use_map: bool,
               figsize=(11, 5), dpi=100):
    units, cmap, _ = FIELDS[field]
    frames = load_frames(field, days)
    if not frames:
        print(f"  [{field}] no data found across {len(days)} day(s), skipping")
        return
    vmin, vmax = value_range(field, frames)

    from matplotlib.animation import FuncAnimation, PillowWriter

    fig = plt.figure(figsize=figsize, dpi=dpi)

    if use_map:
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        lon, lat = grid_coords(field)
        ax = plt.axes(projection=ccrs.PlateCarree())
        ax.set_global()
        ax.add_feature(cfeature.LAND, facecolor="lightgray", zorder=0)
        ax.coastlines(resolution="50m", color="black", linewidth=0.5, zorder=2)
        mesh = ax.pcolormesh(
            lon, lat, frames[0][1],
            transform=ccrs.PlateCarree(), shading="auto",
            cmap=cmap, vmin=vmin, vmax=vmax, zorder=1,
        )
    else:
        # Raw mode: no lat/lon transform, no coastlines -- just the array
        # exactly as stored on disk (row 0 at the top, as numpy/imshow would
        # naturally show it), to verify orientation independent of any
        # geo-referencing assumption.
        ax = plt.axes()
        ax.set_xlabel("col (0-based)")
        ax.set_ylabel("row (0-based)")
        mesh = ax.imshow(frames[0][1], origin="upper", cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")

    fig.colorbar(mesh, ax=ax, orientation="vertical", shrink=0.7, label=f"{field} [{units}]")
    title = ax.set_title(f"{field} -- {frames[0][0]}")

    def update(i):
        date_str, data = frames[i]
        if use_map:
            mesh.set_array(data.ravel())
        else:
            mesh.set_data(data)
        title.set_text(f"{field} -- {date_str}")
        return mesh, title

    anim = FuncAnimation(fig, update, frames=len(frames), blit=False)
    anim.save(out_path, writer=PillowWriter(fps=fps))
    plt.close(fig)
    print(f"  [{field}] wrote {out_path} ({len(frames)} frame(s))")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default="test_output", help="root of daily-averaged output (default: test_output)")
    parser.add_argument("--out-dir", default="gifs", help="directory to write GIFs into (default: gifs)")
    parser.add_argument("--fps", type=float, default=2.0, help="frames per second (default: 2)")
    parser.add_argument("--fields", nargs="+", choices=list(FIELDS), default=list(FIELDS), help="subset of fields to render")
    parser.add_argument("--step", type=int, default=1,
                        help="use every Nth day (e.g. 7 for weekly; default 1 = every day)")
    parser.add_argument("--dpi", type=int, default=100,
                        help="output resolution in dots-per-inch (default: 100)")
    parser.add_argument("--figsize", type=float, nargs=2, metavar=("W", "H"), default=(11, 5),
                        help="figure width and height in inches (default: 11 5)")
    parser.add_argument(
        "--no-map", action="store_true",
        help="skip the cartopy map/coastlines and render the raw array as-is "
             "(row 0 at top, no lat/lon transform) -- use this to verify "
             "orientation independent of any geo-referencing assumption",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    days = discover_days(output_dir)
    if not days:
        sys.exit(f"No <YYYY>/<MM>/<DD> directories found under {output_dir}")

    days = days[::args.step]
    print(f"Rendering {len(days)} day(s) (step={args.step}): {days[0][0]} .. {days[-1][0]}")

    figsize = tuple(args.figsize)
    for field in args.fields:
        render_gif(field, days, out_dir / f"{field}.gif", args.fps,
                   use_map=not args.no_map, figsize=figsize, dpi=args.dpi)


if __name__ == "__main__":
    main()
