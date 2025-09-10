#!/usr/bin/env python3
"""
plot_hysteresis.py

Usage:
  python3 scripts/plot_hysteresis.py --data output/data_preprocessed.csv \
      --model output/model_fit.csv --poly output/poly_fit.csv --metrics output/metrics.json --which both --save output/plot.png

This script:
 - plots preprocessed data (H vs M) as points
 - overlays parametric model (if provided) as green line
 - overlays polynomial up (red) and down (blue) from poly CSV (if provided)
 - annotates metrics for data, parametric, polynomial (if metrics.json contains them)
 - optionally saves PNG with --save <path>
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import json
import os
import sys

def load_csv(path):
    try:
        return pd.read_csv(path)
    except Exception as e:
        print(f"Error reading {path}: {e}", file=sys.stderr)
        return None

def annotate_metrics(ax, metrics, x0, color, label_prefix):
    txt = (f"{label_prefix}:\n"
           f"Ms = {metrics['Ms']:.4g}\n"
           f"Mr = {metrics['Mr']:.4g}\n"
           f"Hc = {metrics['Hc']:.4g}\n"
           f"Area= {metrics['Area']:.4g}")
    ax.text(x0, 0.98, txt, transform=ax.transAxes,
            fontsize=9, verticalalignment='top',
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'),
            color=color)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", required=True, help="CSV: H,M")
    parser.add_argument("--model", help="Parametric model CSV: H,ModelM")
    parser.add_argument("--poly", help="Polynomial CSV: H,Mup,Mdown")
    parser.add_argument("--metrics", required=True, help="metrics JSON from C++")
    parser.add_argument("--which", default="both", choices=["data","model","both"],
                        help="what to plot")
    parser.add_argument("--save", help="If provided, save the figure to this path")
    args = parser.parse_args()

    data = load_csv(args.data)
    model = load_csv(args.model) if args.model else None
    poly = load_csv(args.poly) if args.poly else None

    # load metrics JSON
    try:
        with open(args.metrics, "r") as fh:
            metrics = json.load(fh)
    except Exception as e:
        print(f"Error reading metrics JSON: {e}", file=sys.stderr)
        metrics = {}

    fig, ax = plt.subplots(figsize=(7,7))

    # raw data scatter
    if data is not None and args.which in ("data","both"):
        if "H" in data.columns and "M" in data.columns:
            ax.scatter(data["H"], data["M"], s=8, alpha=0.6, label="Data (preprocessed)")
        else:
            print("Data CSV missing H/M columns", file=sys.stderr)

    # parametric model line
    if model is not None and args.which in ("model","both"):
        if "H" in model.columns and "ModelM" in model.columns:
            ax.plot(model["H"], model["ModelM"], '-', linewidth=1.5, label="Parametric model", color="green")
        else:
            print("Model CSV missing H/ModelM columns", file=sys.stderr)

    # polynomial up/down
    if poly is not None and args.which in ("model","both"):
        # poly CSV format expected: H,Mup,Mdown
        if {"H","Mup","Mdown"}.issubset(set(poly.columns)):
            ax.plot(poly["H"], poly["Mup"], '-', linewidth=1.2, label="Poly up", color="red")
            ax.plot(poly["H"], poly["Mdown"], '-', linewidth=1.2, label="Poly down", color="blue")
            # fill between (visualize loop area)
            try:
                ax.fill_between(poly["H"], poly["Mup"], poly["Mdown"], color="gray", alpha=0.15)
            except Exception:
                pass
        else:
            print("Poly CSV missing expected columns H,Mup,Mdown", file=sys.stderr)

    ax.set_xlabel("H (Oe)")
    ax.set_ylabel("M (normalized)")
    ax.grid(True)
    ax.legend()

    # annotations: arrange left, center, right
    left_x = 0.02; mid_x = 0.35; right_x = 0.68
    if "data" in metrics:
        annotate_metrics(ax, metrics["data"], left_x, "black", "Data")
    if "parametric" in metrics:
        annotate_metrics(ax, metrics["parametric"], mid_x, "green", "Parametric")
    if "polynomial" in metrics:
        annotate_metrics(ax, metrics["polynomial"], right_x, "red", "Polynomial")

    plt.tight_layout()

    if args.save:
        outdir = os.path.dirname(args.save)
        if outdir and not os.path.exists(outdir):
            os.makedirs(outdir, exist_ok=True)
        plt.savefig(args.save, dpi=200)
        print(f"Saved figure to {args.save}")
    else:
        plt.show()

if __name__ == "__main__":
    main()