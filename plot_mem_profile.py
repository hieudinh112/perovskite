#!/usr/bin/env python
"""Plot qcpbc.exe RSS over time for one THC-AO-K run, with the run's phases shaded.

Usage: plot_mem_profile.py RUN_DIR OUT.png [--limit GiB] [--model GiB]

RUN_DIR is the run directory made by run_qcpbc_fasrc.sbatch: mem_samples.txt ("epoch rss_kb=N" every ~10 s) and
the job's .out. Phase boundaries come from the .out timers, anchored at the end of the metric pinv, which is the
first sample after the global peak below 40% of it (every THC-AO-K run so far peaks inside the pinv).
"""
import argparse
import datetime as dt
import glob
import os
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.dates as mdates
import matplotlib.pyplot as plt

SURFACE, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS, BAND = "#e1e0d9", "#c3c2b7", "#f0efec"
SERIES, CRITICAL = "#2a78d6", "#d03b3b"


def timer(out, label):
    m = re.search(re.escape(label) + r"\s+wall\s+([0-9.]+) s", out)
    return float(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run_dir")
    ap.add_argument("png")
    ap.add_argument("--limit", type=float, help="job memory request (GiB)")
    ap.add_argument("--model", type=float, help="predicted pinv peak (GiB)")
    a = ap.parse_args()

    t, g = [], []
    for line in open(os.path.join(a.run_dir, "mem_samples.txt")):
        f = line.split()
        if len(f) == 2 and f[1].startswith("rss_kb="):
            t.append(float(f[0]))
            g.append(int(f[1][7:]) / 1048576)
    out_path = sorted(glob.glob(os.path.join(a.run_dir, "*.out")))[0]
    out = open(out_path).read()
    nk = 1 + max(int(x) for x in re.findall(r"== (\d+)-\w\w k-point nmo", out))
    nbsf = int(re.search(r"There are (\d+) basis functions in the central unit cell", out).group(1))
    nisdf = int(re.search(r"(?:Loaded|Saved) (\d+) ISDF interpolation points", out).group(1))
    done = "Total job time" in out

    ipk = max(range(len(g)), key=g.__getitem__)
    ipk0 = next(i for i in range(len(g)) if g[i] >= g[ipk] - 0.05)  # first sample at the peak value
    ipk1 = max(i for i in range(len(g)) if g[i] >= g[ipk] - 0.05)  # last
    iend = next(i for i in range(ipk, len(g)) if g[i] < 0.4 * g[ipk])
    t_pinv1 = t[iend]
    t_pinv0 = t_pinv1 - timer(out, "S1 metric pinv")
    t_m0 = t_pinv0 - timer(out, "S0 Gram build")
    bands = []
    sel = timer(out, "ISDF point selection:")
    if sel and sel > 60:
        t_sel1 = t_m0 - (timer(out, "THCvec interp-point collocation:") or 0)
        bands += [("setup", t[0], t_sel1 - sel), ("ISDF point selection", t_sel1 - sel, t_m0)]
    else:
        bands += [("setup + Gram", t[0], t_pinv0)]
    bands.append((f"metric pinv: {nk} SVDs at once", t_pinv0, t_pinv1))
    mwall = timer(out, "M-matrix build:")
    if mwall:
        bands += [("block loop", t_pinv1, t_m0 + mwall), ("SCF", t_m0 + mwall, t[-1])]
    else:
        bands.append(("block loop (running)", t_pinv1, t[-1]))

    T = [dt.datetime.fromtimestamp(x) for x in t]
    D = lambda x: dt.datetime.fromtimestamp(x)
    ymax = max(a.limit or 0, g[ipk]) * 1.12

    fig, (ax, az) = plt.subplots(1, 2, figsize=(13, 5.2), dpi=160, gridspec_kw=dict(width_ratios=[2.2, 1]),
                                 facecolor=SURFACE)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.80, bottom=0.11, wspace=0.12)
    name = os.path.basename(os.path.normpath(a.run_dir))
    fig.text(0.065, 0.945, f"Memory peaks inside the metric pinv: {g[ipk]:.1f} GiB", fontsize=14,
             fontweight="semibold", color=INK)
    through = "finished" if done else f"running, data through {T[-1]:%m-%d %H:%M}"
    fig.text(0.065, 0.895, f"{name}  ·  qcpbc.exe RSS, 10 s samples  ·  {nk} k-points, nbsf {nbsf}, "
             f"nisdf {nisdf}  ·  {through}", fontsize=9.5, color=INK2)

    for axx in (ax, az):
        axx.set_facecolor(SURFACE)
        for i, (lab, b0, b1) in enumerate(bands):
            if i % 2 == 1:
                axx.axvspan(D(b0), D(b1), color=BAND, lw=0, zorder=0)
        axx.fill_between(T, g, color=SERIES, alpha=0.10, lw=0, zorder=1)
        axx.plot(T, g, color=SERIES, lw=2, solid_joinstyle="round", solid_capstyle="round", zorder=3)
        if a.limit:
            axx.axhline(a.limit, color=CRITICAL, lw=1, zorder=2)
        if a.model:
            axx.axhline(a.model, color=MUTED, lw=1, zorder=2)
        axx.plot([T[ipk0]], [g[ipk0]], "o", ms=8, color=SERIES, mec=SURFACE, mew=2, zorder=4)
        axx.set_ylim(0, ymax)
        axx.grid(axis="y", color=GRID, lw=1)
        axx.set_axisbelow(True)
        for s in ("top", "right"):
            axx.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            axx.spines[s].set_color(AXIS)
        axx.tick_params(colors=MUTED, labelsize=9, length=3)
        axx.tick_params(axis="x", labelcolor=INK2)
        axx.tick_params(axis="y", labelcolor=INK2)
    ax.set_ylabel("RSS (GiB)", color=INK2, fontsize=10)
    az.tick_params(axis="y", labelleft=False)

    # Whole run: band names along the top, one value label per phase.
    ax.set_xlim(T[0], T[-1])
    span = t[-1] - t[0]
    for lab, b0, b1 in bands:
        if (b1 - b0) / span > 0.07:
            ax.text(D((b0 + b1) / 2), ymax * 0.965, lab, ha="center", va="top", fontsize=9, color=INK2)
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    i_mid = min(range(len(t)), key=lambda i: abs(t[i] - (t_pinv0 + 0.35 * (t_pinv1 - t_pinv0))))
    ax.text(T[i_mid], g[i_mid] + ymax * 0.025, f"{g[i_mid]:.0f} GiB", ha="center", va="bottom", fontsize=9,
            color=INK)
    ax.text(T[ipk0], g[ipk0] + ymax * 0.03, f"peak {g[ipk]:.1f} GiB\n{T[ipk0]:%H:%M}-{T[ipk1]:%H:%M}",
            ha="right", va="bottom", fontsize=9, color=INK)
    for lab, b0, b1 in bands:
        ga = [y for x, y in zip(t, g) if b0 + 60 < x <= b1]
        if b0 >= t_pinv1 and ga:
            ax.text(D((b0 + b1) / 2), max(ga) + ymax * 0.025, f"{min(ga):.0f}-{max(ga):.0f} GiB", ha="center",
                    va="bottom", fontsize=9, color=INK)
    if a.limit:
        ax.text(T[-1] - (T[-1] - T[0]) * 0.01, a.limit - ymax * 0.012, f"job request {a.limit:.0f} GiB",
                ha="right", va="top", fontsize=9, color=INK2)
    if a.model:
        ax.text(T[0] + (T[-1] - T[0]) * 0.01, a.model + ymax * 0.012, f"memory model {a.model:.0f} GiB",
                ha="left", va="bottom", fontsize=9, color=INK2)

    # Zoom on the end of the pinv, where the peak forms.
    z0, z1 = t[ipk0] - 2700, t_pinv1 + 720
    az.set_xlim(D(z0), D(z1))
    az.xaxis.set_major_locator(mdates.MinuteLocator(byminute=[0, 30]))
    az.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    az.set_title("Zoom: end of the pinv", fontsize=10, color=INK2, loc="left")
    az.text(T[ipk0], g[ipk0] + ymax * 0.03, f"{g[ipk]:.1f} GiB", ha="center", va="bottom", fontsize=9,
            color=INK)
    az.text(D(t_pinv1 - 60), 0.4 * g[ipk], f"pinv ends\n{D(t_pinv1):%H:%M}", ha="right", va="center",
            fontsize=9, color=INK)
    az.axvline(D(t_pinv1), color=AXIS, lw=1, zorder=2)

    fig.savefig(a.png, facecolor=SURFACE)
    print(f"wrote {a.png}")
    print(f"peak {g[ipk]:.1f} GiB {T[ipk0]:%H:%M:%S}-{T[ipk1]:%H:%M:%S}; pinv {D(t_pinv0):%H:%M:%S}-"
          f"{D(t_pinv1):%H:%M:%S}; M build starts {D(t_m0):%H:%M:%S}")


if __name__ == "__main__":
    main()
