#!/usr/bin/env python3
"""Compare computed room-air replacement, temperature, and exchange rate."""
import csv
import json
import math
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parent
SCENARIOS = ("wide_open", "tilted")


def read_table(case, name, filename, columns):
    files = list((case / "postProcessing" / name).glob(f"*/{filename}"))
    if len(files) != 1:
        raise ValueError(f"{case.name}: expected one {name} history; found {len(files)}")
    text = files[0].read_text()
    rows = []
    for line in text.splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = [float(value) for value in line.split()]
        if len(values) != columns or not all(map(math.isfinite, values)):
            raise ValueError(f"invalid data in {files[0]}: {line}")
        if rows and values[0] <= rows[-1][0]:
            raise ValueError(f"non-increasing times in {files[0]}")
        rows.append(values)
    if not rows:
        raise ValueError(f"no data in {files[0]}")
    return text, rows


def history(case):
    info = json.loads((case / "case_info.json").read_text())
    header, stats = read_table(case, "roomStats", "volFieldValue.dat", 3)
    _, flux = read_table(case, "windowFlux", "surfaceFieldValue.dat", 2)
    volume = re.search(r"^#\s*Volume\s*:\s*(\S+)", header, re.M | re.I)
    if not volume:
        raise ValueError(f"{case.name}: missing room volume in statistics header")
    volume = float(volume.group(1))
    if not math.isfinite(volume) or volume <= 0:
        raise ValueError(f"{case.name}: invalid room volume")
    if len(stats) != len(flux):
        raise ValueError(f"{case.name}: room and window histories have different lengths")
    result = []
    for (time, original, temperature), (flux_time, sum_abs_phi) in zip(stats, flux):
        if not math.isclose(time, flux_time, abs_tol=1e-7, rel_tol=1e-8):
            raise ValueError(f"{case.name}: unsynchronised room and window histories")
        if sum_abs_phi < 0:
            raise ValueError(f"{case.name}: negative absolute flux")
        # With an incompressible room and one opening, inflow = outflow.
        exchange = 0.5 * sum_abs_phi
        result.append(dict(scenario=case.name, time_s=time,
                           original_air_fraction=original,
                           fresh_air_percent=100 * (1 - original),
                           mean_temperature_C=temperature - 273.15,
                           exchange_m3_s=exchange, exchange_ACH=3600 * exchange / volume,
                           # Air inventory only: radiator input is not integrated.
                           air_sensible_energy_decrease_J=1.2 * 1005 * volume *
                           (info["room_temperature"] - temperature)))
    return info, volume, result


def main():
    output = ROOT / "comparison.csv"
    histories = [history(ROOT / "cases" / scenario)[2] for scenario in SCENARIOS]
    times = [[row["time_s"] for row in rows] for rows in histories]
    if len(times[0]) != len(times[1]) or any(
        not math.isclose(a, b, abs_tol=1e-7, rel_tol=1e-8) for a, b in zip(*times)
    ):
        raise ValueError("the two cases must cover the same output times")
    with output.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(histories[0][0]))
        writer.writeheader()
        for rows in histories:
            writer.writerows(rows)
    for rows in histories:
        row = rows[-1]
        print(f"{row['scenario']}: t={row['time_s']:g} s, "
              f"fresh air={row['fresh_air_percent']:.2f}%, "
              f"mean T={row['mean_temperature_C']:.2f} C, "
              f"instantaneous exchange={row['exchange_ACH']:.2f} ACH")
    print(f"Wrote {output.name}")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        print(f"Cannot analyse: {error}. Run ./Allrun or ./run_docker.sh first.", file=sys.stderr)
        sys.exit(2)
