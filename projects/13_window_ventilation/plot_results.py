#!/usr/bin/env python3
"""Plot comparison.csv as images/results.svg (standard library only).

Run after analyze.py. Left: share of room air that came from outside (tracer).
Right: volume-averaged room temperature.
"""
import csv
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
COLOURS = {"wide_open": "#d9622b", "tilted": "#2a7fb8"}
LABELS = {"wide_open": "Wide open", "tilted": "Tilted 10°"}
PANEL_W, PANEL_H, LEFT, TOP, GAP = 400, 260, 60, 50, 90


def load(path):
    series = {}
    with path.open() as stream:
        for row in csv.DictReader(stream):
            series.setdefault(row["scenario"], []).append(
                (float(row["time_s"]), float(row["fresh_air_percent"]),
                 float(row["mean_temperature_C"])))
    return series


def ticks(lo, hi, count=5):
    step = (hi - lo) / count
    return [lo + i * step for i in range(count + 1)]


def panel(x0, title, series, column, y_range, y_label, y_ticks=5):
    t_end = max(rows[-1][0] for rows in series.values())
    y_lo, y_hi = y_range

    def sx(t):
        return x0 + PANEL_W * t / t_end

    def sy(v):
        return TOP + PANEL_H * (1 - (v - y_lo) / (y_hi - y_lo))

    out = [f'<text x="{x0}" y="{TOP - 16}" class="title">{title}</text>',
           f'<rect x="{x0}" y="{TOP}" width="{PANEL_W}" height="{PANEL_H}" class="frame"/>']
    for v in ticks(y_lo, y_hi, y_ticks):
        out.append(f'<line x1="{x0}" x2="{x0 + PANEL_W}" y1="{sy(v):.1f}" y2="{sy(v):.1f}" class="grid"/>')
        out.append(f'<text x="{x0 - 8}" y="{sy(v) + 4:.1f}" class="tick" text-anchor="end">{v:g}</text>')
    for t in ticks(0, t_end):
        out.append(f'<text x="{sx(t):.1f}" y="{TOP + PANEL_H + 18}" class="tick" '
                   f'text-anchor="middle">{t:g}</text>')
    out.append(f'<text x="{x0 + PANEL_W / 2}" y="{TOP + PANEL_H + 38}" class="axis" '
               'text-anchor="middle">time (s)</text>')
    out.append(f'<text transform="translate({x0 - 42},{TOP + PANEL_H / 2}) rotate(-90)" '
               f'class="axis" text-anchor="middle">{y_label}</text>')
    for scenario, rows in series.items():
        points = " ".join(f"{sx(r[0]):.1f},{sy(r[column]):.1f}" for r in [(0, 0, 20)] + rows)
        out.append(f'<polyline points="{points}" fill="none" stroke="{COLOURS[scenario]}" '
                   'stroke-width="2.5"/>')
        last = rows[-1]
        out.append(f'<text x="{sx(last[0]) - 4:.1f}" y="{sy(last[column]) - 8:.1f}" '
                   f'class="value" text-anchor="end" fill="{COLOURS[scenario]}">'
                   f'{last[column]:.1f}</text>')
    return out


def main():
    source = ROOT / "comparison.csv"
    if not source.is_file():
        print("Missing comparison.csv; run analyze.py first.", file=sys.stderr)
        return 2
    series = load(source)
    width = 2 * PANEL_W + GAP + LEFT + 30
    height = TOP + PANEL_H + 90
    temperatures = [r[2] for rows in series.values() for r in rows] + [20]
    t_lo = 4 * int(min(temperatures) // 4)  # 4 °C grid steps
    body = panel(LEFT, "Room air that came from outside", series, 1, (0, 100), "fresh air (%)")
    body += panel(LEFT + PANEL_W + GAP, "Mean room temperature", series, 2, (t_lo, 20), "°C",
                  (20 - t_lo) // 4)
    legend_x = LEFT
    for scenario in series:
        body.append(f'<line x1="{legend_x}" x2="{legend_x + 24}" y1="{height - 18}" '
                    f'y2="{height - 18}" stroke="{COLOURS[scenario]}" stroke-width="3"/>')
        body.append(f'<text x="{legend_x + 30}" y="{height - 13}" class="axis">'
                    f'{LABELS.get(scenario, scenario)}</text>')
        legend_x += 140
    svg = "\n".join([
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        f'<rect width="{width}" height="{height}" fill="#f7fafc"/>',
        "<style>text{font-family:Arial,sans-serif;fill:#1e3445}.title{font-size:16px;"
        "font-weight:bold}.tick{font-size:12px;fill:#4a5d6b}.axis{font-size:13px}"
        ".value{font-size:13px;font-weight:bold}.frame{fill:#fff;stroke:#9fb0bc}"
        ".grid{stroke:#e3e9ee}</style>",
        *body, "</svg>", ""])
    target = ROOT / "images" / "results.svg"
    target.write_text(svg)
    print(f"Wrote {target.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
