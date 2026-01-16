#!/usr/bin/env python3
"""
Plot shell center distances versus shell index from one or more
"*shells_statistics.dat" files produced by identify_shells.py.

Usage examples:
  python plot_shell_centers.py nanoonion_shells_6shells_statistics.dat
  python plot_shell_centers.py stats1.dat stats2.dat --output shell_centers.png
"""

import argparse
import os
from typing import List, Tuple, Dict

import matplotlib.pyplot as plt


def parse_shell_statistics_file(path: str) -> List[Tuple[int, float]]:
    """Parse a *shells_statistics.dat file and return list of (shell_id, center_distance).

    The file is expected to contain header lines starting with '#', followed by
    whitespace-separated columns where the first two are Shell_ID and Center_Distance(Angstrom).
    """
    shell_id_to_center: Dict[int, float] = {}
    with open(path, "r") as f:
        for line in f:
            if not line.strip() or line.lstrip().startswith('#'):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                shell_id = int(parts[0])
                center_distance = float(parts[1])
                shell_id_to_center[shell_id] = center_distance
            except ValueError:
                # Skip malformed rows
                continue

    # Sort by shell_id to ensure correct ordering on x-axis
    sorted_items = sorted(shell_id_to_center.items(), key=lambda kv: kv[0])
    return sorted_items


def build_label_from_filename(path: str) -> str:
    """Create a concise label from the file name."""
    base = os.path.basename(path)
    name, _ = os.path.splitext(base)
    return name


def plot_shell_centers(series: List[Tuple[str, List[Tuple[int, float]]]], output: str = None) -> str:
    """Plot shell centers for one or more series.

    series: list of (label, [(shell_id, center_distance), ...])
    output: optional output file path. If not provided, a name is derived from first label.

    Returns the output file path used.
    """
    fig, ax = plt.subplots(1, 1, figsize=(10, 6))

    for label, items in series:
        shell_indices = [sid for sid, _ in items]
        center_distances = [cd for _, cd in items]
        ax.plot(shell_indices, center_distances, marker='o', linewidth=2, label=label)

    ax.set_xlabel('Shell Index (n)')
    ax.set_ylabel('Shell Center Distance (Angstrom)')
    ax.set_title('Shell Center Distances vs Shell Index')
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()

    if output is None:
        default_label = series[0][0] if series else 'shell_centers'
        output = f"{default_label}_centers_plot.png"

    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close()
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description='Plot shell center distances vs shell index from *shells_statistics.dat files')
    parser.add_argument('stats_files', nargs='+', help='One or more *shells_statistics.dat files')
    parser.add_argument('--output', help='Output PNG file path (default derived from first file)')

    args = parser.parse_args()

    series: List[Tuple[str, List[Tuple[int, float]]]] = []
    for path in args.stats_files:
        if not os.path.exists(path):
            print(f"Warning: {path} not found; skipping")
            continue
        data = parse_shell_statistics_file(path)
        if not data:
            print(f"Warning: {path} contained no usable rows; skipping")
            continue
        label = build_label_from_filename(path)
        series.append((label, data))

    if not series:
        print("Error: No valid data to plot.")
        return

    output_path = plot_shell_centers(series, args.output)
    print(f"Shell centers plot saved: {output_path}")


if __name__ == '__main__':
    main()


