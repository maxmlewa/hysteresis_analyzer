#!/usr/bin/env python3
"""
Plot the hysteresis loop data and model fits

Usage: python3 plot_hysteresis.py -- data output/data_preprocessed.csv \
                                  -- model output/model_fit.csv \
                                  -- metrics output/metrics.json \
                                  -- which both \
                                  -- save output/ loop_plot.png

Options:
--data    CSV with columns: H,M (preprocessed data)
--model   CSV with columns: H,ModelM (parametric model sampled)
--metrics JSON file containing derived metrics for 'data' and 'model'
--which  one of: data | model | both    (default: both)
--save   path to save PNG (if omitted, interactive window is shown)
--dpi    DPI for saved PNG (default 150)
--marker-size size of scatter markers (default 4)
--linewidth model line width (default 1.5)
"""
import argparse
import json
import os
import sys

import matplotlib.pyplot as plt
import numpy as np

def read_csv(path):
    try:
        import csv
        H = []
        M = []
        with open(path, 'r') as f:
            rdr = csv.reader(f)
            header = next(rdr)

            #figure column indices
            idxH = 0
            idxM = 1

            if any(h.lower().startswith('h') for h in header): #trying to detecct headers
                for i,h in enumerate(header):
                    if h.strip().lower().startswith('h'): idxH = i
                    if h.strip().lower().startswith('m'): idxM = i

            else: #then the header is numeric and should be treated as raw data
                f.seek(0)
                rdr = csv.reader(f)

            for row in rdr:
                if len(row) < max(idxH, idxM) + 1: continue
                try:
                    H.append(float(row[idxH]))
                    M.append(float(row[idxM]))
                except:
                    continue
        return np.array(H), np.array(M)
    except Exception as e:
        print("Error reading the CSV:", e)
        sys.exit(1)



def annotate_metrics(ax, x, y, label_prefix = "Data", color = None):
    #x, y are dicts with the keys as the parameters: Ms, Mr, Hc, Area
    if x is None:
        return

    txt = (f"{label_prefix}:\n"
           f"Ms = {x.get('Ms', np.nan):.4g}:\n"
           f"Mr = {x.get('Mr', np.nan):.4g}:\n"
           f"Hc = {x.get('Hc', np.nan):.4g}:\n"
           f"Area = {x.get('Area', np.nan):.4g}:\n")

    #place the annotation to the top-left/ change according to the plot forms

    ax.txt(0.02, 0.98 if color is None else 0.98, txt,
           transform= ax.transAxes, fontsize = 9,
           verticalalignments = 'top', 
           bbox = dict(facecolor='white', alpha=0.8, edgecolor='none'))
    
def main():
    parser = argparse.ArgumentParser(description="Plot hysteresis data & model")
    parser.add_argument("--data", type=str, help="CSV: H,M (preprocessed data)")
    parser.add_argument("--model", type=str, help="CSV: H,ModelM (model samples)")
    parser.add_argument("--metrics", type=str, help="JSON metrics file")
    parser.add_argument("--which", type=str, choices=['data','model','both'], default='both',
                        help="Which series to display")
    parser.add_argument("--save", type=str, default=None, help="Save PNG to path")
    parser.add_argument("--dpi", type=int, default=150)
    parser.add_argument("--marker-size", type=float, default=4.0)
    parser.add_argument("--linewidth", type=float, default=1.5)
    args = parser.parse_args()

    H_data = M_data = None
    H_model = M_model = None
    metrics = {}

    if args.data:
        H_data, M_data = read_csv(args.data)
    if args.model:
        H_model, M_model = read_csv(args.model)
    if args.metrics and os.path.exists(args.metrics):
        try:
            with open(args.metrics, 'r') as f:
                metrics = json.load(f)
        except Exception as e:
            print("Warning: could not read metrics JSON:", e)

    fig, ax = plt.subplots(figsize=(8,6))
    ax.set_xlabel("Magnetizing Field H (Oe)")
    ax.set_ylabel("Magnetization M (normalized)")
    ax.grid(True, linestyle=':', linewidth=0.6)

    plotted = False
    if args.which in ('data','both') and H_data is not None:
        ax.scatter(H_data, M_data, s=args.marker_size, label="Preprocessed data", alpha=0.9)
        plotted = True
    if args.which in ('model','both') and H_model is not None:
        # sort model by H for clean line
        order = np.argsort(H_model)
        ax.plot(H_model[order], M_model[order], linewidth=args.linewidth, label="Parametric model")
        plotted = True

    # zero axes lines
    ax.axhline(0, color='k', linewidth=0.6, linestyle='-')
    ax.axvline(0, color='k', linewidth=0.6, linestyle='-')

    # annotate metrics if available
    data_metrics = metrics.get('data', None)
    model_metrics = metrics.get('model', None)
    annotate_metrics(ax, data_metrics, None, label_prefix="Data")
    annotate_metrics(ax, model_metrics, None, label_prefix="Model")

    if not plotted:
        print("No data plotted. Provide --data and/or --model files.")
        return

    ax.legend(loc='best')
    fig.tight_layout()

    if args.save:
        plt.savefig(args.save, dpi=args.dpi)
        print("Saved plot to:", args.save)
    else:
        plt.show()

if __name__ == "__main__":
    main()