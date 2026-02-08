#!/usr/bin/env python3
"""
charge_balance_vhh_windows.py

Streaming analysis of charge "balancing" across IMGT-position windows for VHH datasets.

Works with your current schema:
- IMGT columns present: typically 66-104, 118-128, plus 42,49,50,52 (may vary)
- Computes per-window net charge, then correlation matrices:
    1) overall raw correlation
    2) overall residualized on covariates (default: total_charge + len_cdr3 + len_v)
    3) within-family pooled raw correlation (controls for family)
    4) within-family pooled residual correlation (controls for family + covariates)

Outputs:
- windows.tsv
- corr_raw.overall.tsv
- corr_resid.overall.tsv
- corr_raw.within_family_pooled.tsv
- corr_resid.within_family_pooled.tsv
- family_counts.tsv (optional but useful)

Example:
  python3 charge_balance_vhh_windows.py \
    --inputs data/databases/annotated/vhh_full_annotated_v4/*dedup.csv \
    --outdir out_charge_bal \
    --window-size 8 \
    --require-valid
"""

from __future__ import annotations

import argparse
import os
import glob
import time
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

import numpy as np
import pandas as pd


AA_CHARGE = {"K": 1.0, "R": 1.0, "D": -1.0, "E": -1.0}


def aa_to_charge(s: pd.Series, his_charge: float = 0.0) -> np.ndarray:
    """
    Convert a column of single-letter amino acids to charge values.
    Robust to pd.NA by coercing NA -> "" prior to comparisons.
    """
    arr = s.astype("string").fillna("").str.upper().to_numpy()

    out = np.zeros((len(arr),), dtype=np.float32)
    out += (arr == "K").astype(np.float32)
    out += (arr == "R").astype(np.float32)
    out -= (arr == "D").astype(np.float32)
    out -= (arr == "E").astype(np.float32)

    if his_charge != 0.0:
        out += (arr == "H").astype(np.float32) * np.float32(his_charge)

    return out


def is_trueish(series: pd.Series) -> np.ndarray:
    if series.dtype == bool:
        return series.to_numpy()
    vv = series.astype("string").str.lower()
    return (vv == "true") | (vv == "1") | (vv == "t") | (vv == "yes")


def detect_imgt_positions_from_header(path: str) -> List[int]:
    cols = list(pd.read_csv(path, nrows=0).columns)
    pos = []
    for c in cols:
        if c.startswith("imgt_"):
            try:
                pos.append(int(c.split("_", 1)[1]))
            except ValueError:
                pass
    return sorted(set(pos))


def intersect_imgt_positions(paths: List[str]) -> List[int]:
    """Intersection of numeric IMGT positions across all input files (safe)."""
    sets = []
    for p in paths:
        sets.append(set(detect_imgt_positions_from_header(p)))
    inter = set.intersection(*sets) if sets else set()
    return sorted(inter)


def group_into_runs(sorted_positions: List[int]) -> List[List[int]]:
    """Group positions into consecutive runs: e.g., [42],[49,50],[66..104],[118..128]."""
    if not sorted_positions:
        return []
    runs = [[sorted_positions[0]]]
    for x in sorted_positions[1:]:
        if x == runs[-1][-1] + 1:
            runs[-1].append(x)
        else:
            runs.append([x])
    return runs


def build_windows_from_runs(runs: List[List[int]], window_size: int) -> List[Tuple[int, int, List[int]]]:
    """
    For each run, bin into fixed windows of length window_size in numeric space,
    but only include positions actually present.
    """
    windows: List[Tuple[int, int, List[int]]] = []
    for run in runs:
        start = run[0]
        end = run[-1]
        run_set = set(run)
        a = start
        while a <= end:
            b = min(end, a + window_size - 1)
            poslist = [p for p in range(a, b + 1) if p in run_set]
            if poslist:
                windows.append((a, b, poslist))
            a = b + 1
    return windows


@dataclass
class Accumulator:
    W: int
    covariates: List[str]   # without intercept
    p: int                 # with intercept
    n: int = 0
    sumY: Optional[np.ndarray] = None  # (W,)
    Syy: Optional[np.ndarray] = None   # (W,W)
    Sxx: Optional[np.ndarray] = None   # (p,p)
    Sxy: Optional[np.ndarray] = None   # (p,W)

    def __post_init__(self):
        self.sumY = np.zeros((self.W,), dtype=np.float64)
        self.Syy = np.zeros((self.W, self.W), dtype=np.float64)
        self.Sxx = np.zeros((self.p, self.p), dtype=np.float64)
        self.Sxy = np.zeros((self.p, self.W), dtype=np.float64)

    def update(self, X: np.ndarray, Y: np.ndarray):
        self.n += Y.shape[0]
        self.sumY += Y.sum(axis=0)
        self.Syy += Y.T @ Y
        self.Sxx += X.T @ X
        self.Sxy += X.T @ Y

    def cov_raw(self) -> np.ndarray:
        if self.n < 2:
            return np.full((self.W, self.W), np.nan)
        return (self.Syy - np.outer(self.sumY, self.sumY) / self.n) / (self.n - 1)

    def corr_from_cov(self, cov: np.ndarray) -> np.ndarray:
        var = np.diag(cov)
        denom = np.sqrt(np.outer(var, var))
        with np.errstate(divide="ignore", invalid="ignore"):
            corr = cov / denom
        np.fill_diagonal(corr, 1.0)
        return corr

    def corr_raw(self) -> np.ndarray:
        return self.corr_from_cov(self.cov_raw())

    def cov_resid(self) -> np.ndarray:
        if self.n <= self.p:
            return np.full((self.W, self.W), np.nan)
        try:
            invSxx = np.linalg.inv(self.Sxx)
        except np.linalg.LinAlgError:
            return np.full((self.W, self.W), np.nan)
        Sres = self.Syy - (self.Sxy.T @ invSxx @ self.Sxy)
        df = self.n - self.p
        return Sres / df

    def corr_resid(self) -> np.ndarray:
        return self.corr_from_cov(self.cov_resid())


def build_X(chunk: pd.DataFrame, covariates: List[str], total_charge: np.ndarray) -> np.ndarray:
    cols = [np.ones((len(chunk),), dtype=np.float64)]
    for c in covariates:
        if c == "total_charge":
            cols.append(total_charge.astype(np.float64))
        else:
            if c not in chunk.columns:
                raise ValueError(f"Covariate '{c}' not found in columns.")
            cols.append(pd.to_numeric(chunk[c], errors="coerce").to_numpy(dtype=np.float64))
    X = np.vstack(cols).T
    return X


def write_windows(path: str, windows: List[Tuple[int, int, List[int]]]):
    rows = []
    for i, (a, b, pos) in enumerate(windows):
        rows.append({
            "window_index": i,
            "imgt_start": a,
            "imgt_end": b,
            "n_pos": len(pos),
            "positions": ",".join(map(str, pos))
        })
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def write_matrix(path: str, labels: List[str], mat: np.ndarray):
    pd.DataFrame(mat, index=labels, columns=labels).to_csv(path, sep="\t", float_format="%.6g")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs="+", required=True, help="Input CSVs or globs")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--chunk-size", type=int, default=250_000)
    ap.add_argument("--window-size", type=int, default=8)
    ap.add_argument("--require-valid", action="store_true")
    ap.add_argument("--his-charge", type=float, default=0.0)
    ap.add_argument("--covariates", nargs="*", default=["total_charge", "len_cdr3", "len_v"])
    ap.add_argument("--min-family-count", type=int, default=5000,
                    help="Families with fewer rows than this are ignored in within-family pooled stats.")

    # Progress/debug
    ap.add_argument("--log-every", type=int, default=10,
                    help="Print progress every N chunks total across all files (0 disables).")
    ap.add_argument("--max-chunks", type=int, default=0,
                    help="If >0, stop after this many chunks total (debug).")

    args = ap.parse_args()

    # Expand globs
    paths: List[str] = []
    for x in args.inputs:
        paths.extend(glob.glob(x))
    paths = sorted(set(paths))
    if not paths:
        raise SystemExit("No input files matched.")

    os.makedirs(args.outdir, exist_ok=True)

    kept_positions = intersect_imgt_positions(paths)
    if not kept_positions:
        raise SystemExit("No numeric imgt_# columns found across all inputs.")

    runs = group_into_runs(kept_positions)
    windows = build_windows_from_runs(runs, window_size=args.window_size)
    if not windows:
        raise SystemExit("No windows formed (unexpected).")

    pos_cols = [f"imgt_{p}" for p in kept_positions]
    pos_to_idx = {p: i for i, p in enumerate(kept_positions)}
    W = len(windows)
    labels = [f"{a}-{b}" for (a, b, _) in windows]

    write_windows(os.path.join(args.outdir, "windows.tsv"), windows)

    p = 1 + len(args.covariates)
    overall = Accumulator(W=W, covariates=args.covariates, p=p)

    fam_acc: Dict[str, Accumulator] = {}
    fam_counts: Dict[str, int] = {}

    t0 = time.time()
    print(f"[start] files={len(paths)} windows={W} positions={len(kept_positions)} "
          f"chunk_size={args.chunk_size} window_size={args.window_size} covariates={args.covariates}",
          flush=True)

    chunks_total = 0
    rows_total = 0

    for path in paths:
        print(f"[file] {os.path.basename(path)}", flush=True)

        for chunk in pd.read_csv(path, chunksize=args.chunk_size, low_memory=False):
            if args.require_valid and "valid" in chunk.columns:
                chunk = chunk[is_trueish(chunk["valid"])]
            if len(chunk) == 0:
                continue

            chunks_total += 1
            rows_total += len(chunk)

            # Debug stop
            if args.max_chunks and chunks_total >= args.max_chunks:
                print(f"[stop] reached --max-chunks={args.max_chunks}", flush=True)
                break

            # Progress log
            if args.log_every and (chunks_total % args.log_every == 0):
                dt = time.time() - t0
                rate = rows_total / dt if dt > 0 else 0.0
                print(f"[prog] chunks={chunks_total} rows={rows_total:,} elapsed={dt:,.1f}s "
                      f"rate={rate:,.0f} rows/s (current={os.path.basename(path)})",
                      flush=True)

            missing = [c for c in pos_cols if c not in chunk.columns]
            if missing:
                raise SystemExit(f"{path} missing IMGT columns (intersection should prevent this?): {missing[:10]}")

            n = len(chunk)
            P = len(pos_cols)
            Q = np.zeros((n, P), dtype=np.float32)
            for j, c in enumerate(pos_cols):
                Q[:, j] = aa_to_charge(chunk[c], his_charge=args.his_charge)

            Y = np.zeros((n, W), dtype=np.float32)
            for w_idx, (_, _, poslist) in enumerate(windows):
                idxs = [pos_to_idx[p] for p in poslist]
                Y[:, w_idx] = Q[:, idxs].sum(axis=1)

            total_charge = Q.sum(axis=1).astype(np.float64)
            X = build_X(chunk, args.covariates, total_charge=total_charge)

            ok = np.isfinite(X).all(axis=1)
            if not np.all(ok):
                X = X[ok]
                Y = Y[ok]
                if Y.shape[0] == 0:
                    continue
                chunk = chunk.loc[chunk.index[ok]]

            overall.update(X=X, Y=Y.astype(np.float64))

            if "family" in chunk.columns:
                fam = chunk["family"].astype("string").fillna("NA").to_numpy()
                for f in np.unique(fam):
                    m = (fam == f)
                    nf = int(m.sum())
                    if nf == 0:
                        continue
                    if f not in fam_acc:
                        fam_acc[f] = Accumulator(W=W, covariates=args.covariates, p=p)
                        fam_counts[f] = 0
                    fam_acc[f].update(X=X[m], Y=Y[m].astype(np.float64))
                    fam_counts[f] += nf

        if args.max_chunks and chunks_total >= args.max_chunks:
            break

    write_matrix(os.path.join(args.outdir, "corr_raw.overall.tsv"), labels, overall.corr_raw())
    write_matrix(os.path.join(args.outdir, "corr_resid.overall.tsv"), labels, overall.corr_resid())

    if fam_acc:
        kept = [(f, fam_counts[f]) for f in fam_acc.keys() if fam_counts[f] >= args.min_family_count]
        kept.sort(key=lambda x: x[1], reverse=True)

        S_within = np.zeros((W, W), dtype=np.float64)
        N_total = 0
        n_fam = 0
        for f, nf in kept:
            acc = fam_acc[f]
            if acc.n < 2:
                continue
            S_within += (acc.Syy - np.outer(acc.sumY, acc.sumY) / acc.n)
            N_total += acc.n
            n_fam += 1

        if N_total > n_fam + 1:
            cov_within = S_within / (N_total - n_fam)
            corr_within = overall.corr_from_cov(cov_within)
            write_matrix(os.path.join(args.outdir, "corr_raw.within_family_pooled.tsv"), labels, corr_within)

        Sres_within = np.zeros((W, W), dtype=np.float64)
        df_total = 0
        for f, nf in kept:
            acc = fam_acc[f]
            if acc.n <= acc.p:
                continue
            try:
                invSxx = np.linalg.inv(acc.Sxx)
            except np.linalg.LinAlgError:
                continue
            Sres_f = acc.Syy - (acc.Sxy.T @ invSxx @ acc.Sxy)
            df_f = acc.n - acc.p
            Sres_within += Sres_f
            df_total += df_f

        if df_total > 2:
            cov_res_within = Sres_within / df_total
            corr_res_within = overall.corr_from_cov(cov_res_within)
            write_matrix(os.path.join(args.outdir, "corr_resid.within_family_pooled.tsv"), labels, corr_res_within)

        pd.DataFrame(
            [{"family": f, "n": fam_counts[f]} for f, _ in sorted(fam_counts.items(), key=lambda x: x[1], reverse=True)]
        ).to_csv(os.path.join(args.outdir, "family_counts.tsv"), sep="\t", index=False)

    dt = time.time() - t0
    rate = rows_total / dt if dt > 0 else 0.0
    print(f"[done] outdir={args.outdir} rows={rows_total:,} chunks={chunks_total} "
          f"elapsed={dt:,.1f}s rate={rate:,.0f} rows/s",
          flush=True)
    print(f"[done] IMGT positions used: {len(kept_positions)} in {len(runs)} runs; windows={W}; rows_accumulated={overall.n}",
          flush=True)


if __name__ == "__main__":
    main()
