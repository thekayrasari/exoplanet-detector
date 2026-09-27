#!/usr/bin/env python3
"""
exoplanet_pipeline.py
=====================
Unified, single-file Automated Exoplanet Discovery Pipeline.

Merges three previously-separate scripts into one file with zero extra
file dependencies beyond standard scientific Python packages:
  - the core BLS discovery/vetting engine
  - the AFK daily batch scheduler/distributor
  - the Tkinter desktop GUI

No arguments -> launches the GUI. With arguments -> CLI subcommands:

    python exoplanet_pipeline.py                      # launch GUI
    python exoplanet_pipeline.py gui                   # launch GUI (explicit)
    python exoplanet_pipeline.py scan --star "TOI-700"
    python exoplanet_pipeline.py scan --targets sample_targets.txt --outdir scan_results
    python exoplanet_pipeline.py afk --time-limit 60 --source TIC --workers 3
    python exoplanet_pipeline.py self-test

Install: pip install lightkurve numpy matplotlib scipy pandas astroquery
"""

# ==========================================================================
# IMPORTS
# ==========================================================================
import os
import sys
import csv
import json
import time
import math
import argparse
import functools
import threading
import subprocess
import urllib.request as _urllib_request
from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # headless backend; GUI uses subprocesses, never imports pyplot directly
import matplotlib.pyplot as plt

print = functools.partial(print, flush=True)

try:
    import lightkurve as lk
except ImportError:
    lk = None  # only truly required by `scan`/`afk`/`self-test`; GUI can still launch to show the warning

import warnings
warnings.filterwarnings('ignore')


# ==========================================================================
# CONSTANTS
# ==========================================================================
R_SUN_TO_R_EARTH = 109.076

CSV_HEADERS = [
    "Target", "TIC_ID", "Mission", "Sectors_Found", "Planet_Label",
    "Period_days", "Epoch_BJD", "Depth_percent", "Duration_hours", "SNR",
    "Odd_Even_Ratio", "Stellar_Radius_Rsun", "Stellar_Radius_Source",
    "Planet_Radius_Rearth", "Planet_Type", "Multi_Planet_Count",
    "Exofop_Status", "Vetting_Status", "Plot_Path", "Discovery_Timestamp",
]

STATE_FILENAME = "scanned_state.json"
REQUIRED_PACKAGES = ["lightkurve", "pandas", "numpy", "matplotlib", "scipy", "astroquery"]

_OFFLINE_FALLBACK_TARGETS = [
    "TIC 25155310", "TIC 42446519", "TIC 44713209", "TIC 88846397", "TIC 149603524",
    "TIC 144432041", "TIC 201505066", "TIC 207265987", "TIC 260919380", "TIC 261136674",
    "TIC 269701147", "TIC 278956474", "TIC 282350145", "TIC 300914493", "TIC 347436069",
    "TIC 350757263", "TIC 55525572", "TIC 99709840", "TIC 150428135", "TIC 182091106",
    "TIC 219990199", "TIC 231663901", "TIC 251456914", "TIC 272018806", "TIC 280198812",
    "TIC 298435447", "TIC 307210830", "TIC 308538095",
]

THIS_SCRIPT = str(Path(__file__).resolve())


# ==========================================================================
# CROSS-PROCESS FILE LOCK
# ==========================================================================
class SimpleFileLock:
    """Minimal cross-platform, cross-process lock via an atomic O_CREAT|O_EXCL
    lockfile. Protects the catalog CSV/JSON both across threads in one process
    and across the parallel `scan` subprocesses spawned by `afk`."""

    def __init__(self, protected_path, timeout=30.0, poll_interval=0.05):
        self.lock_path = str(protected_path) + ".lock"
        self.timeout = timeout
        self.poll_interval = poll_interval
        self._fd = None

    def __enter__(self):
        start = time.time()
        while True:
            try:
                self._fd = os.open(self.lock_path, os.O_CREAT | os.O_EXCL | os.O_RDWR)
                return self
            except FileExistsError:
                if time.time() - start > self.timeout:
                    try:
                        os.remove(self.lock_path)  # reclaim a stale lock from a crashed process
                    except OSError:
                        pass
                    start = time.time()
                    continue
                time.sleep(self.poll_interval)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fd is not None:
            try:
                os.close(self._fd)
            except OSError:
                pass
        try:
            os.remove(self.lock_path)
        except OSError:
            pass
        return False


# ==========================================================================
# CORE DISCOVERY & VETTING ENGINE
# ==========================================================================
class AutomatedExoplanetFinder:
    def __init__(self, output_dir="scan_results", min_period=0.5, max_period=30.0,
                 snr_threshold=6.0, max_sectors=3, plot_mode="interesting",
                 toi_cache_hours=24.0, bls_target_grid_points=100000):
        self.output_dir = Path(output_dir)
        self.plots_dir = self.output_dir / "plots"
        self.plots_dir.mkdir(parents=True, exist_ok=True)

        self.min_period = min_period
        self.max_period = max_period
        self.snr_threshold = snr_threshold
        self.max_sectors = max_sectors
        self.plot_mode = plot_mode  # "all" | "interesting" | "pass_only"
        self.toi_cache_hours = toi_cache_hours
        # Target POINT COUNT rather than a raw astropy `frequency_factor` - see
        # _adaptive_frequency_factor() for why the raw parameter is dangerous
        # to hardcode (its scaling direction is counter-intuitive).
        self.bls_target_grid_points = bls_target_grid_points

        self.csv_path = self.output_dir / "candidate_catalog.csv"
        self.json_path = self.output_dir / "candidate_summary.json"

        self.exofop_toi_df = None
        self._seen_signals = set()  # in-memory duplicate-protection index

        self._load_exofop_catalog()
        self._init_catalog_files()
        self._load_existing_catalog_keys()

    # ------------------------------------------------------------------
    # Catalog bootstrap / duplicate protection
    # ------------------------------------------------------------------
    def _load_exofop_catalog(self):
        """Pre-fetch or load a cached NASA ExoFOP TOI catalog. The cache expires
        after `toi_cache_hours` (default 24h) instead of being reused forever."""
        try:
            cache_file = self.output_dir / "exofop_toi_cache.csv"
            cache_is_fresh = (
                cache_file.exists()
                and os.path.getsize(cache_file) > 1000
                and (time.time() - cache_file.stat().st_mtime) < self.toi_cache_hours * 3600.0
            )
            if cache_is_fresh:
                self.exofop_toi_df = pd.read_csv(cache_file)
                print(f"📋 Using cached ExoFOP TOI catalog ({cache_file.name}, "
                      f"{(time.time() - cache_file.stat().st_mtime) / 3600.0:.1f}h old).")
            else:
                url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
                df = pd.read_csv(url)
                if df is not None and len(df) > 0:
                    df.to_csv(cache_file, index=False)
                    self.exofop_toi_df = df
                    print(f"📋 Refreshed ExoFOP TOI catalog cache ({len(df)} rows).")
        except Exception as e:
            try:
                cache_file = self.output_dir / "exofop_toi_cache.csv"
                if cache_file.exists():
                    self.exofop_toi_df = pd.read_csv(cache_file)
                    print(f"⚠️ ExoFOP live refresh failed ({e}); using stale cached copy.")
                else:
                    self.exofop_toi_df = None
                    print(f"⚠️ ExoFOP catalog unavailable ({e}); cross-referencing disabled.")
            except Exception:
                self.exofop_toi_df = None

    def _init_catalog_files(self):
        if not self.csv_path.exists():
            with open(self.csv_path, "w", newline="", encoding="utf-8") as f:
                csv.writer(f).writerow(CSV_HEADERS)

    @staticmethod
    def _make_dedup_key(target_name, period, epoch, planet_label):
        """Epoch is normalized into [0, period) so the same physical transit
        reported with an epoch offset by whole cycles (common across different
        multi-sector baselines) is still recognized as the same signal."""
        period = float(period)
        epoch_norm = float(epoch) % period if period else float(epoch)
        return (str(target_name).strip(), round(period, 3), round(epoch_norm, 3), str(planet_label))

    def _load_existing_catalog_keys(self):
        if not self.csv_path.exists():
            return
        try:
            df = pd.read_csv(self.csv_path)
            required = {"Target", "Period_days", "Epoch_BJD"}
            if not required.issubset(df.columns):
                archive_path = self.output_dir / f"candidate_catalog_legacy_{int(time.time())}.csv"
                self.csv_path.rename(archive_path)
                print(f"⚠️ Existing catalog uses an older schema. Archived to "
                      f"'{archive_path.name}' and started a fresh catalog.")
                self._init_catalog_files()
                return

            labels = df["Planet_Label"] if "Planet_Label" in df.columns else pd.Series(["b"] * len(df))
            loaded = 0
            for target, period, epoch, label in zip(df["Target"], df["Period_days"], df["Epoch_BJD"], labels):
                try:
                    self._seen_signals.add(self._make_dedup_key(target, period, epoch, label))
                    loaded += 1
                except (TypeError, ValueError):
                    continue
            if loaded:
                print(f"📇 Loaded {loaded} existing catalog entries for duplicate protection.")
        except Exception as e:
            print(f"⚠️ Could not pre-load existing catalog for dedup ({e}); "
                  f"duplicate protection will only cover this session.")

    # ------------------------------------------------------------------
    # TIC resolution + real stellar radius
    # ------------------------------------------------------------------
    def resolve_target_catalog_info(self, target_name):
        """Resolve a target's real TIC ID and stellar radius from the TESS
        Input Catalog via astroquery (not from the lightkurve search-result
        table, which has no stellar parameters)."""
        info = {"tic_id": None, "r_star": 1.0, "r_star_source": "DEFAULT_ASSUMED_SOLAR"}
        try:
            from astroquery.mast import Catalogs
            digits = "".join(c for c in target_name if c.isdigit())
            if target_name.strip().upper().startswith("TIC") and digits:
                table = Catalogs.query_criteria(catalog="TIC", ID=int(digits))
            else:
                table = Catalogs.query_object(target_name, radius=0.01, catalog="TIC")

            if table is not None and len(table) > 0:
                row = table[0]
                if "ID" in table.colnames and row["ID"] is not None:
                    info["tic_id"] = int(row["ID"])
                if "rad" in table.colnames:
                    rad_val = row["rad"]
                    if rad_val is not None and not (isinstance(rad_val, float) and np.isnan(rad_val)) and rad_val > 0:
                        info["r_star"] = float(rad_val)
                        info["r_star_source"] = "TIC_CATALOG"
        except Exception as e:
            print(f"⚠️ TIC catalog lookup failed for '{target_name}': {e}. Using R*=1.0 Rsun default.")

        if info["r_star_source"] == "TIC_CATALOG":
            print(f"⭐ Host Star Radius (TIC {info['tic_id']}): {info['r_star']:.2f} R_sun")
        else:
            print("⭐ Host Star Radius: defaulted to 1.00 R_sun (TIC catalog lookup unavailable)")
        return info

    def calculate_planet_radius(self, depth_percent, r_star=1.0):
        depth_fraction = max(0.0, depth_percent / 100.0)
        r_planet = r_star * np.sqrt(depth_fraction) * R_SUN_TO_R_EARTH

        if r_planet < 1.25:
            p_type = "Earth-sized (<1.25 Re)"
        elif r_planet < 2.0:
            p_type = "Super-Earth (1.25-2.0 Re)"
        elif r_planet < 4.0:
            p_type = "Sub-Neptune (2.0-4.0 Re)"
        elif r_planet < 8.0:
            p_type = "Neptune-like (4.0-8.0 Re)"
        else:
            p_type = "Gas Giant / Jovian (>8.0 Re)"
        return r_planet, p_type

    def cross_reference_exofop(self, tic_id, detected_period):
        if tic_id is None:
            return "UNLISTED_TARGET"
        if self.exofop_toi_df is None or len(self.exofop_toi_df) == 0:
            return "EXOFOP_UNAVAILABLE"
        try:
            df = self.exofop_toi_df
            tic_col = "TIC ID" if "TIC ID" in df.columns else ("TIC" if "TIC" in df.columns else None)
            period_col = "Period (days)" if "Period (days)" in df.columns else ("Period" if "Period" in df.columns else None)
            toi_col = "TOI" if "TOI" in df.columns else None

            if tic_col and period_col:
                matches = df[df[tic_col] == tic_id]
                if len(matches) > 0:
                    for _, row in matches.iterrows():
                        known_p = float(row[period_col]) if pd.notnull(row[period_col]) else None
                        toi_name = f"TOI-{row[toi_col]}" if toi_col and pd.notnull(row[toi_col]) else "Known TOI"
                        if known_p and known_p > 0 and abs(detected_period - known_p) / known_p < 0.04:
                            return f"KNOWN_TOI ({toi_name}, P={known_p:.2f}d)"
                        elif known_p and known_p > 0 and (
                            abs(detected_period - 2 * known_p) / (2 * known_p) < 0.04
                            or abs(detected_period - known_p / 2) / (known_p / 2) < 0.04
                        ):
                            return f"KNOWN_TOI_HARMONIC ({toi_name}, P={known_p:.2f}d)"
                    return "NEW_PLANET_IN_KNOWN_SYSTEM"
            return "NEW_UNREPORTED_CANDIDATE"
        except Exception:
            return "CROSS_REF_ERROR"

    # ------------------------------------------------------------------
    # Data acquisition
    # ------------------------------------------------------------------
    def fetch_and_preprocess(self, target_name, mission="TESS"):
        print(f"\n{'='*70}")
        print(f"📡 Fetching data for: {target_name} ({mission})")
        print(f"{'='*70}")

        catalog_info = self.resolve_target_catalog_info(target_name)

        try:
            search_result = lk.search_lightcurve(target_name, mission=mission, author="SPOC")
            if len(search_result) == 0:
                search_result = lk.search_lightcurve(target_name, mission=mission)
            if len(search_result) == 0:
                print(f"❌ No light curves found for '{target_name}'.")
                return None, 0, catalog_info

            max_sec = getattr(self, 'max_sectors', 3)
            if len(search_result) > max_sec:
                search_result = search_result[-max_sec:]

            print(f"✓ Found {len(search_result)} sector observation(s). Downloading...")

            downloaded_lcs = []
            for i, item in enumerate(search_result, 1):
                try:
                    print(f"  ⬇️ Downloading sector {i}/{len(search_result)}...")
                    single_lc = item.download()
                    if single_lc is not None:
                        downloaded_lcs.append(single_lc)
                except Exception as e:
                    print(f"  ⚠️ Skipping sector {i} due to download error: {e}")
                    continue

            if not downloaded_lcs:
                print(f"❌ Could not download valid light curves for '{target_name}'.")
                return None, 0, catalog_info

            lc_collection = lk.LightCurveCollection(downloaded_lcs)
            lc_stitched = lc_collection.stitch().remove_nans()

            window_length = min(1001, len(lc_stitched) // 5)
            if window_length % 2 == 0:
                window_length += 1
            if window_length < 21:
                window_length = 21

            lc_flat = lc_stitched.flatten(window_length=window_length)
            lc_clean = lc_flat.remove_outliers(sigma=3.5)

            print(f"✓ Stitched & cleaned: {len(lc_clean)} data points across {len(search_result)} sector(s).")
            return lc_clean, len(search_result), catalog_info

        except Exception as e:
            print(f"❌ Error downloading/processing {target_name}: {e}")
            return None, 0, catalog_info

    # ------------------------------------------------------------------
    # BLS search
    # ------------------------------------------------------------------
    @staticmethod
    def _adaptive_frequency_factor(min_p, max_p, time_baseline_days, target_points, min_duration=0.05):
        """Pick astropy's BLS `frequency_factor` to hit a target grid size.
        Per astropy's `BoxLeastSquares.autoperiod` docs, df = frequency_factor *
        min(duration) / (time_baseline)**2, so a LARGER factor makes the grid
        COARSER. We solve for the factor that gives ~target_points grid points,
        so resolution stays good regardless of observing baseline length."""
        if time_baseline_days <= 0 or min_p <= 0 or max_p <= min_p:
            return 10.0
        denom = (1.0 / min_p) - (1.0 / max_p)
        if denom <= 0:
            return 10.0
        ideal_factor = denom * (time_baseline_days ** 2) / (min_duration * target_points)
        return max(1.0, ideal_factor)

    def run_bls_pass(self, lc, min_period=None, max_period=None):
        min_p = min_period if min_period else self.min_period
        time_span = lc.time.max().value - lc.time.min().value
        max_p = max_period if max_period else min(self.max_period, time_span / 2.0)
        if max_p <= min_p:
            max_p = min_p + 2.0

        freq_factor = self._adaptive_frequency_factor(min_p, max_p, time_span, self.bls_target_grid_points)
        bls = lc.to_periodogram(method='bls', minimum_period=min_p, maximum_period=max_p, frequency_factor=freq_factor)

        best_period = float(bls.period_at_max_power.value)
        best_epoch = float(bls.transit_time_at_max_power.value)
        best_duration = float(bls.duration_at_max_power.value)
        best_depth = float(bls.depth_at_max_power.value)

        in_transit_mask = bls.get_transit_mask(period=best_period, transit_time=best_epoch, duration=best_duration)
        out_of_transit_flux = lc.flux[~in_transit_mask].value
        out_of_transit_std = float(np.nanstd(out_of_transit_flux)) if len(out_of_transit_flux) > 0 else 1.0
        snr = float(best_depth / out_of_transit_std) if out_of_transit_std > 0 else 0.0

        return {
            "period": best_period, "epoch": best_epoch, "duration": best_duration * 24.0,
            "depth_percent": best_depth * 100.0, "snr": snr,
            "bls_object": bls, "in_transit_mask": in_transit_mask,
        }

    def run_multi_planet_search(self, lc, r_star=1.0):
        """2-Pass BLS search: primary planet, then mask and search for a
        secondary sibling planet. Both signals are returned so the caller can
        persist both (the original bug this pipeline had was discarding pass2)."""
        print("🔍 Pass 1: Running primary planet BLS search...")
        pass1_results = self.run_bls_pass(lc)
        r_p1, p_type1 = self.calculate_planet_radius(pass1_results["depth_percent"], r_star)
        pass1_results["planet_radius_rearth"], pass1_results["planet_type"] = r_p1, p_type1

        print(f"  - Primary Period: {pass1_results['period']:.4f} days")
        print(f"  - Depth:          {pass1_results['depth_percent']:.3f}% (Rp = {r_p1:.2f} R_earth -> {p_type1})")
        print(f"  - SNR:            {pass1_results['snr']:.2f}")

        pass2_results = None
        multi_planet_count = 1

        if pass1_results["snr"] >= 4.0:
            print("🪐 Pass 2: Masking primary transit to search for secondary sibling planets...")
            mask1 = pass1_results["in_transit_mask"]
            lc_residual = lc[~mask1]

            if len(lc_residual) > 1000:
                try:
                    pass2_raw = self.run_bls_pass(lc_residual)
                    period_diff = abs(pass2_raw["period"] - pass1_results["period"]) / pass1_results["period"]
                    if pass2_raw["snr"] >= 5.0 and period_diff > 0.06:
                        r_p2, p_type2 = self.calculate_planet_radius(pass2_raw["depth_percent"], r_star)
                        pass2_raw["planet_radius_rearth"], pass2_raw["planet_type"] = r_p2, p_type2
                        pass2_results = pass2_raw
                        multi_planet_count = 2
                        print("  ✨ MULTI-PLANET DETECTED! Secondary Planet Candidate:")
                        print(f"     - Period: {pass2_raw['period']:.4f} days, Depth: {pass2_raw['depth_percent']:.3f}% (Rp = {r_p2:.2f} R_earth)")
                except Exception as e:
                    print(f"  (Pass 2 search skipped: {e})")

        return pass1_results, pass2_results, multi_planet_count

    # ------------------------------------------------------------------
    # Vetting
    # ------------------------------------------------------------------
    def vet_candidate(self, lc, bls_results):
        print("🛡️ Running Automated False-Positive Vetting Engine...")

        period = bls_results["period"]
        epoch = bls_results["epoch"]
        duration_days = bls_results["duration"] / 24.0
        depth_pct = bls_results["depth_percent"]
        snr = bls_results["snr"]

        if snr < self.snr_threshold:
            print(f"  ⚠️ Low SNR ({snr:.2f} < {self.snr_threshold}).")
            return "WEAK_SIGNAL", 1.0, {"odd_depth": 0.0, "even_depth": 0.0, "ratio": 1.0}

        if depth_pct > 4.0:
            print(f"  ⚠️ Excessive depth ({depth_pct:.2f}% > 4.0%). Flagged as Eclipsing Binary / Starspot.")
            return "FALSE_POSITIVE_EB", 1.0, {"odd_depth": depth_pct, "even_depth": depth_pct, "ratio": 1.0}

        time_vals = lc.time.value
        flux_vals = lc.flux.value

        transit_indices = np.round((time_vals - epoch) / period)
        phase = ((time_vals - epoch + 0.5 * period) % period) - 0.5 * period
        half_dur = duration_days / 2.0

        in_transit = np.abs(phase) <= half_dur
        odd_mask = in_transit & (np.abs(transit_indices) % 2 == 1)
        even_mask = in_transit & (np.abs(transit_indices) % 2 == 0)
        out_transit = np.abs(phase) > (half_dur * 2.0)

        out_baseline = np.median(flux_vals[out_transit]) if np.any(out_transit) else 1.0
        odd_depth = (out_baseline - np.median(flux_vals[odd_mask])) * 100.0 if np.any(odd_mask) else depth_pct
        even_depth = (out_baseline - np.median(flux_vals[even_mask])) * 100.0 if np.any(even_mask) else depth_pct

        # A vanished even-transit depth is itself a classic EB signature - don't
        # force ratio=1.0 ("looks fine") near that boundary. Only treat the
        # ratio as clean when BOTH depths are near the noise floor; the pass/
        # fail decision below uses an absolute-difference threshold that scales
        # with transit depth so shallow real mismatches aren't missed.
        depth_floor = 0.02  # percent
        diff = odd_depth - even_depth
        if abs(odd_depth) < depth_floor and abs(even_depth) < depth_floor:
            ratio = 1.0
        else:
            denom = even_depth if abs(even_depth) >= depth_floor else math.copysign(depth_floor, even_depth or 1.0)
            ratio = float(np.clip(odd_depth / denom, -50.0, 50.0))

        vetting_details = {"odd_depth": float(odd_depth), "even_depth": float(even_depth), "ratio": float(ratio)}

        mismatch_threshold = max(0.15, 0.30 * depth_pct)
        if abs(ratio - 1.0) > 0.35 or abs(diff) > mismatch_threshold:
            print(f"  ⚠️ Odd/Even depth mismatch (Odd: {odd_depth:.3f}%, Even: {even_depth:.3f}%, Ratio: {ratio:.2f}). Flagged as Eclipsing Binary.")
            return "FALSE_POSITIVE_EB", ratio, vetting_details

        print(f"  ✅ Vetting PASSED! Planetary candidate verified (SNR={snr:.1f}, Odd/Even={ratio:.2f}).")
        return "PASS_CANDIDATE", ratio, vetting_details

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------
    def _should_generate_plot(self, vetting_status):
        if self.plot_mode == "all":
            return True
        if self.plot_mode == "pass_only":
            return vetting_status == "PASS_CANDIDATE"
        return vetting_status in ("PASS_CANDIDATE", "FALSE_POSITIVE_EB")

    def generate_diagnostic_plot(self, target_name, lc, bls_results, vetting_status, vetting_details,
                                  r_star, exofop_status, pass2_results=None, planet_label="b"):
        fig = plt.figure(figsize=(15, 11))
        r_p = bls_results["planet_radius_rearth"]
        p_type = bls_results["planet_type"]

        header_title = (
            f"Exoplanet Candidate Report: {target_name} [{planet_label}] [{vetting_status}]\n"
            f"Planet Size: {r_p:.2f} R_⊕ ({p_type}) | Host Star: {r_star:.2f} R_☉ | Status: {exofop_status}"
        )
        fig.suptitle(header_title, fontsize=13, fontweight='bold', y=0.98)

        period = bls_results["period"]
        epoch = bls_results["epoch"]
        bls = bls_results["bls_object"]

        gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.2], hspace=0.35, wspace=0.25)

        ax1 = fig.add_subplot(gs[0, :])
        ax1.plot(lc.time.value, lc.flux.value, 'k.', alpha=0.25, markersize=1.5, label="Detrended Flux")
        ax1.set_title("1. Cleaned & Stitched TESS Light Curve", fontsize=11, fontweight='bold')
        ax1.set_xlabel("Time (BJD)"); ax1.set_ylabel("Normalized Flux"); ax1.grid(True, alpha=0.3)

        ax2 = fig.add_subplot(gs[1, 0])
        ax2.plot(bls.period.value, bls.power.value, 'b-', lw=1.2, label="BLS Power")
        ax2.axvline(period, color='r', linestyle='--', label=f"Candidate P = {period:.4f} d")
        if pass2_results:
            ax2.axvline(pass2_results["period"], color='green', linestyle=':', label=f"Sibling P2 = {pass2_results['period']:.4f} d")
        ax2.set_title("2. BLS Periodogram Power Spectrum", fontsize=11, fontweight='bold')
        ax2.set_xlabel("Orbital Period (days)"); ax2.set_ylabel("BLS Power")
        ax2.legend(loc='upper right', fontsize=9); ax2.grid(True, alpha=0.3)

        ax3 = fig.add_subplot(gs[1, 1])
        folded = lc.fold(period=period, epoch_time=epoch)
        ax3.plot(folded.phase.value, folded.flux.value, '.', color='gray', alpha=0.2, markersize=1.5)
        binned = folded.bin(time_bin_size=period / 80.0)
        ax3.plot(binned.phase.value, binned.flux.value, 'r.', markersize=5, label="Binned Flux")
        ax3.set_title(f"3. Phase Folded (P = {period:.4f} d, Epoch = {epoch:.2f})", fontsize=11, fontweight='bold')
        ax3.set_xlabel("Phase (days)"); ax3.set_ylabel("Normalized Flux")
        ax3.legend(loc='lower right', fontsize=9); ax3.grid(True, alpha=0.3)

        ax4 = fig.add_subplot(gs[2, 0])
        zoom_window = bls_results["duration"] / 24.0 * 3.0
        mask = np.abs(folded.phase.value) <= zoom_window
        ax4.plot(folded.phase.value[mask] * 24.0, folded.flux.value[mask], '.', color='lightgray', alpha=0.5, markersize=3)
        binned_zoom = folded[mask].bin(time_bin_size=period / 200.0)
        ax4.plot(binned_zoom.phase.value * 24.0, binned_zoom.flux.value, 'ro-', markersize=4, lw=1.5, label="Binned Transit Profile")
        ax4.set_title(f"4. Zoomed Transit Profile (Depth: {bls_results['depth_percent']:.3f}%, Rp: {r_p:.2f} R_⊕)", fontsize=11, fontweight='bold')
        ax4.set_xlabel("Hours from Mid-Transit"); ax4.set_ylabel("Normalized Flux")
        ax4.legend(loc='lower right', fontsize=9); ax4.grid(True, alpha=0.3)

        ax5 = fig.add_subplot(gs[2, 1])
        time_vals = lc.time.value
        phase = ((time_vals - epoch + 0.5 * period) % period) - 0.5 * period
        transit_indices = np.round((time_vals - epoch) / period)
        odd = (np.abs(transit_indices) % 2 == 1) & (np.abs(phase) <= zoom_window)
        even = (np.abs(transit_indices) % 2 == 0) & (np.abs(phase) <= zoom_window)

        ax5.plot(phase[odd] * 24.0, lc.flux.value[odd], 'b.', alpha=0.35, label=f"Odd Transits (Depth: {vetting_details['odd_depth']:.2f}%)")
        ax5.plot(phase[even] * 24.0, lc.flux.value[even], 'g.', alpha=0.35, label=f"Even Transits (Depth: {vetting_details['even_depth']:.2f}%)")
        ax5.set_title(f"5. False-Positive Vetting (Odd/Even Ratio: {vetting_details.get('ratio', 1.0):.2f})", fontsize=11, fontweight='bold')
        ax5.set_xlabel("Hours from Mid-Transit"); ax5.set_ylabel("Normalized Flux")
        ax5.legend(loc='lower right', fontsize=9); ax5.grid(True, alpha=0.3)

        safe_name = target_name.replace(" ", "_").replace("/", "_")
        plot_file = self.plots_dir / f"{safe_name}_{planet_label}_report.png"
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        plt.close(fig)

        print(f"📊 Publication-grade plot saved: {plot_file}")
        return str(plot_file)

    # ------------------------------------------------------------------
    # Catalog writing
    # ------------------------------------------------------------------
    def log_results(self, target_name, tic_id, sectors_found, bls_results, vetting_status, ratio,
                     r_star, r_star_source, exofop_status, multi_planet_count, plot_path, planet_label="b"):
        period = bls_results["period"]
        epoch = bls_results["epoch"]
        key = self._make_dedup_key(target_name, period, epoch, planet_label)

        with SimpleFileLock(self.csv_path):
            if key in self._seen_signals:
                print(f"⚠️ Signal '{planet_label}' for {target_name} (P={period:.4f}d) already in catalog. Skipping duplicate.")
                return
            self._seen_signals.add(key)

            r_p = bls_results["planet_radius_rearth"]
            p_type = bls_results["planet_type"]
            timestamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

            row = {
                "Target": target_name, "TIC_ID": tic_id if tic_id is not None else "",
                "Mission": "TESS", "Sectors_Found": sectors_found, "Planet_Label": planet_label,
                "Period_days": f"{period:.6f}", "Epoch_BJD": f"{epoch:.6f}",
                "Depth_percent": f"{bls_results['depth_percent']:.4f}",
                "Duration_hours": f"{bls_results['duration']:.2f}", "SNR": f"{bls_results['snr']:.2f}",
                "Odd_Even_Ratio": f"{ratio:.2f}", "Stellar_Radius_Rsun": f"{r_star:.2f}",
                "Stellar_Radius_Source": r_star_source, "Planet_Radius_Rearth": f"{r_p:.2f}",
                "Planet_Type": p_type, "Multi_Planet_Count": multi_planet_count,
                "Exofop_Status": exofop_status, "Vetting_Status": vetting_status,
                "Plot_Path": plot_path or "N/A", "Discovery_Timestamp": timestamp,
            }

            with open(self.csv_path, "a", newline="", encoding="utf-8") as f:
                csv.DictWriter(f, fieldnames=CSV_HEADERS).writerow(row)

            catalog = {}
            if self.json_path.exists():
                try:
                    with open(self.json_path, "r", encoding="utf-8") as f:
                        catalog = json.load(f)
                except Exception:
                    catalog = {}

            entry = dict(row)
            entry["multi_planet_count"] = multi_planet_count
            existing_list = catalog.get(target_name, [])
            if not isinstance(existing_list, list):
                existing_list = [existing_list]  # migrate legacy single-object format
            existing_list.append(entry)
            catalog[target_name] = existing_list

            with open(self.json_path, "w", encoding="utf-8") as f:
                json.dump(catalog, f, indent=2)

    # ------------------------------------------------------------------
    # Orchestration
    # ------------------------------------------------------------------
    def process_target(self, target_name):
        try:
            lc, sectors_found, catalog_info = self.fetch_and_preprocess(target_name)
            if lc is None:
                return "FAILED"

            r_star = catalog_info["r_star"]
            r_star_source = catalog_info["r_star_source"]
            tic_id = catalog_info["tic_id"]

            pass1_results, pass2_results, multi_planet_count = self.run_multi_planet_search(lc, r_star)
            if pass1_results is None:
                return "FAILED"

            signals = [("b", pass1_results)]
            if pass2_results is not None:
                signals.append(("c", pass2_results))

            logged_statuses = []
            for label, bls_res in signals:
                exofop_status = self.cross_reference_exofop(tic_id, bls_res["period"])
                if label == "b":
                    print(f"🛰️ NASA ExoFOP Cross-Reference: {exofop_status}")

                vetting_status, ratio, vetting_details = self.vet_candidate(lc, bls_res)

                plot_path = None
                if self._should_generate_plot(vetting_status):
                    plot_path = self.generate_diagnostic_plot(
                        target_name, lc, bls_res, vetting_status, vetting_details,
                        r_star, exofop_status, pass2_results if label == "b" else None,
                        planet_label=label,
                    )

                self.log_results(
                    target_name, tic_id, sectors_found, bls_res, vetting_status, ratio,
                    r_star, r_star_source, exofop_status, multi_planet_count, plot_path,
                    planet_label=label,
                )
                logged_statuses.append(vetting_status)

            return logged_statuses[0] if logged_statuses else "FAILED"

        except Exception as e:
            print(f"❌ Error processing target '{target_name}': {e}")
            return "FAILED"


# ==========================================================================
# SELF-TEST SUITE (zero network access required)
# ==========================================================================
def make_synthetic_lightcurve(period, epoch, duration_days, depth, n_points=20000,
                               time_span=27.0, noise_std=0.0015, seed=42,
                               even_depth_override=None):
    """Synthetic LightCurve with an injected box transit. If even_depth_override
    is set, even-numbered transits get that depth instead (simulates an
    eclipsing binary's classic odd/even asymmetry)."""
    rng = np.random.default_rng(seed)
    t = np.linspace(0, time_span, n_points)
    phase = ((t - epoch + 0.5 * period) % period) - 0.5 * period
    transit_idx = np.round((t - epoch) / period)
    in_transit = np.abs(phase) <= (duration_days / 2.0)

    flux = np.ones(n_points)
    if even_depth_override is not None:
        odd_mask = in_transit & (np.abs(transit_idx) % 2 == 1)
        even_mask = in_transit & (np.abs(transit_idx) % 2 == 0)
        flux[odd_mask] -= depth
        flux[even_mask] -= even_depth_override
    else:
        flux[in_transit] -= depth

    flux += rng.normal(0, noise_std, n_points)
    return lk.LightCurve(time=t, flux=flux)


def run_self_test_suite():
    """Runs all offline regression checks. Returns True iff everything passed."""
    import shutil
    import tempfile

    scratch = Path(tempfile.mkdtemp(prefix="exoplanet_pipeline_selftest_"))
    results = {}

    def new_finder(subdir):
        return AutomatedExoplanetFinder(output_dir=str(scratch / subdir), snr_threshold=6.0, plot_mode="pass_only")

    # --- Test 1: BLS period & depth recovery ---
    print("\n[TEST 1] BLS period & depth recovery on a clean injected transit...")
    finder = new_finder("t1")
    lc = make_synthetic_lightcurve(period=3.14, epoch=1.0, duration_days=0.1, depth=0.01)
    bls = finder.run_bls_pass(lc, min_period=0.5, max_period=10.0)
    r_p, p_type = finder.calculate_planet_radius(bls["depth_percent"], r_star=1.0)
    bls["planet_radius_rearth"], bls["planet_type"] = r_p, p_type
    period_ok = abs(bls["period"] - 3.14) / 3.14 < 0.01
    depth_ok = abs(bls["depth_percent"] - 1.0) < 0.3
    vstatus, ratio, _ = finder.vet_candidate(lc, bls)
    print(f"    period={bls['period']:.4f}d (want ~3.14, {'OK' if period_ok else 'FAIL'})")
    print(f"    depth={bls['depth_percent']:.3f}% (want ~1.0, {'OK' if depth_ok else 'FAIL'})")
    print(f"    vetting={vstatus} ({'OK' if vstatus == 'PASS_CANDIDATE' else 'FAIL'})")
    results["BLS period/depth recovery"] = period_ok and depth_ok and vstatus == "PASS_CANDIDATE"

    # --- Test 2: duplicate protection ---
    print("\n[TEST 2] Duplicate catalog protection...")
    finder = new_finder("t2")
    lc = make_synthetic_lightcurve(period=5.0, epoch=2.0, duration_days=0.12, depth=0.008)
    bls = finder.run_bls_pass(lc, min_period=0.5, max_period=15.0)
    r_p, p_type = finder.calculate_planet_radius(bls["depth_percent"], r_star=1.0)
    bls["planet_radius_rearth"], bls["planet_type"] = r_p, p_type
    vstatus, ratio, _ = finder.vet_candidate(lc, bls)
    for _ in range(3):
        finder.log_results("TEST-DEDUP-STAR", 999999, 1, bls, vstatus, ratio, 1.0, "TIC_CATALOG",
                            "NEW_UNREPORTED_CANDIDATE", 1, None, planet_label="b")
    df = pd.read_csv(finder.csv_path)
    n_rows = len(df[df["Target"] == "TEST-DEDUP-STAR"])
    print(f"    {n_rows} row(s) after 3 identical log_results() calls (want 1, {'OK' if n_rows == 1 else 'FAIL'})")
    results["Duplicate protection"] = n_rows == 1

    # --- Test 3: multi-planet logging ---
    print("\n[TEST 3] Multi-planet secondary signal is actually persisted...")
    finder = new_finder("t3")
    lc1 = make_synthetic_lightcurve(period=4.2, epoch=1.5, duration_days=0.1, depth=0.012)
    bls1 = finder.run_bls_pass(lc1, min_period=0.5, max_period=15.0)
    r_p1, p_type1 = finder.calculate_planet_radius(bls1["depth_percent"], r_star=1.0)
    bls1["planet_radius_rearth"], bls1["planet_type"] = r_p1, p_type1
    v1, r1, _ = finder.vet_candidate(lc1, bls1)
    finder.log_results("TEST-MULTI-STAR", 888888, 1, bls1, v1, r1, 1.0, "TIC_CATALOG",
                        "NEW_UNREPORTED_CANDIDATE", 2, None, planet_label="b")
    lc2 = make_synthetic_lightcurve(period=9.7, epoch=3.0, duration_days=0.08, depth=0.004)
    bls2 = finder.run_bls_pass(lc2, min_period=0.5, max_period=20.0)
    r_p2, p_type2 = finder.calculate_planet_radius(bls2["depth_percent"], r_star=1.0)
    bls2["planet_radius_rearth"], bls2["planet_type"] = r_p2, p_type2
    v2, r2, _ = finder.vet_candidate(lc2, bls2)
    finder.log_results("TEST-MULTI-STAR", 888888, 1, bls2, v2, r2, 1.0, "TIC_CATALOG",
                        "NEW_UNREPORTED_CANDIDATE", 2, None, planet_label="c")
    df = pd.read_csv(finder.csv_path)
    rows = df[df["Target"] == "TEST-MULTI-STAR"]
    labels = set(rows["Planet_Label"])
    ok3 = len(rows) == 2 and labels == {"b", "c"}
    print(f"    {len(rows)} row(s), labels={labels} (want 2, {{'b','c'}}, {'OK' if ok3 else 'FAIL'})")
    results["Multi-planet logging"] = ok3

    # --- Test 4: odd/even EB detection regression ---
    print("\n[TEST 4] Odd/even EB detection when even-transit depth is ~0...")
    finder = new_finder("t4")
    lc = make_synthetic_lightcurve(period=3.5, epoch=1.0, duration_days=0.1, depth=0.006,
                                    even_depth_override=0.0001, noise_std=0.0008)
    # Hand vet_candidate the TRUE injected period/epoch/duration directly rather
    # than routing through run_bls_pass, which would lock onto the sub-harmonic
    # of only the deep transits (a separate, real phenomenon) and defeat this
    # specific unit test of the vetting math.
    bls = {"period": 3.5, "epoch": 1.0, "duration": 0.1 * 24.0, "depth_percent": 0.3, "snr": 7.5}
    r_p, p_type = finder.calculate_planet_radius(bls["depth_percent"], r_star=1.0)
    bls["planet_radius_rearth"], bls["planet_type"] = r_p, p_type
    vstatus, ratio, details = finder.vet_candidate(lc, bls)
    ok4 = vstatus == "FALSE_POSITIVE_EB"
    print(f"    odd={details['odd_depth']:.3f}%, even={details['even_depth']:.3f}%, ratio={ratio:.2f}")
    print(f"    vetting_status={vstatus} (want FALSE_POSITIVE_EB, {'OK' if ok4 else 'FAIL'})")
    results["Odd/even EB detection (regression)"] = ok4

    print(f"\n{'='*60}\nSELF-TEST SUMMARY\n{'='*60}")
    all_ok = True
    for name, ok in results.items():
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        all_ok = all_ok and ok

    shutil.rmtree(scratch, ignore_errors=True)
    print(f"\n{'✅ ALL TESTS PASSED' if all_ok else '❌ SOME TESTS FAILED'}")
    return all_ok


# ==========================================================================
# AFK / DISTRIBUTION LOGIC
# ==========================================================================
def fetch_exofop_toi_list():
    print("🌐 Auto-fetching latest candidate target list from NASA ExoFOP...")
    url = "https://exofop.ipac.caltech.edu/tess/download_toi.php?sort=toi&output=csv"
    try:
        df = pd.read_csv(url)
        col = "TIC ID" if "TIC ID" in df.columns else ("TIC" if "TIC" in df.columns else None)
        if col:
            tic_ids = [f"TIC {int(tic)}" for tic in df[col].dropna().astype(float).unique()]
            print(f"✓ Successfully fetched {len(tic_ids)} candidate targets from ExoFOP.")
            return tic_ids, df
    except Exception as e:
        print(f"⚠️ ExoFOP live fetch warning: {e}. Falling back to default target list.")
    return ["TOI-700", "LHS 1140", "L 98-59", "TIC 251456914", "TOI-1338"], None


def _filter_by_magnitude(tic_ids, max_mag, chunk_size=200):
    from astroquery.mast import Catalogs
    keep = []
    for i in range(0, len(tic_ids), chunk_size):
        batch = tic_ids[i:i + chunk_size]
        try:
            table = Catalogs.query_criteria(catalog="TIC", ID=batch)
            for row in table:
                tmag = row["Tmag"] if "Tmag" in table.colnames else None
                if tmag is not None and not (isinstance(tmag, float) and tmag != tmag) and tmag <= max_mag:
                    keep.append(int(row["ID"]))
        except Exception as e:
            print(f"  ⚠️ Magnitude filter batch failed ({e}); keeping batch unfiltered.")
            keep.extend(batch)
    return keep


def fetch_mast_tic_targets(count=100, max_mag=None, sector=None, exclude_known_tois=True, exofop_df=None):
    """Real MAST query for observed TIC targets, honestly-labeled fallback if
    the live query fails (e.g. no network access to MAST)."""
    print(f"🌐 Querying MAST for up to {count} real observed TESS TIC targets...")
    try:
        from astroquery.mast import Observations
        query_kwargs = dict(obs_collection="TESS", dataproduct_type="timeseries")
        if sector is not None:
            query_kwargs["sequence_number"] = sector
        obs_table = Observations.query_criteria(**query_kwargs)
        if obs_table is None or len(obs_table) == 0:
            raise RuntimeError("MAST returned no timeseries observations for the given criteria")

        raw_names = [str(n) for n in obs_table["target_name"] if n not in (None, "")]
        tic_ids = sorted({int(n) for n in raw_names if n.isdigit()})
        if not tic_ids:
            raise RuntimeError("MAST results contained no parseable numeric TIC target names")
        print(f"  ✓ MAST returned {len(tic_ids)} unique real observed TIC IDs.")

        if max_mag is not None:
            before = len(tic_ids)
            tic_ids = _filter_by_magnitude(tic_ids, max_mag)
            print(f"  ✓ Magnitude filter (Tmag <= {max_mag}): {before} -> {len(tic_ids)} targets.")

        if exclude_known_tois and exofop_df is not None and len(exofop_df) > 0:
            tic_col = "TIC ID" if "TIC ID" in exofop_df.columns else ("TIC" if "TIC" in exofop_df.columns else None)
            if tic_col:
                known = set(int(x) for x in exofop_df[tic_col].dropna().astype(float).unique())
                before = len(tic_ids)
                tic_ids = [t for t in tic_ids if t not in known]
                print(f"  ✓ Excluded {before - len(tic_ids)} already-known TOIs; {len(tic_ids)} unreported remain.")

        if not tic_ids:
            raise RuntimeError("No targets survived filtering")
        return [f"TIC {t}" for t in tic_ids[:count]], "MAST_LIVE"

    except Exception as e:
        print(f"⚠️ Live MAST query unavailable ({e}).")
        print("   Falling back to the built-in list of well-known observed targets (NOT a live MAST result).")
        return _OFFLINE_FALLBACK_TARGETS[:count], "OFFLINE_FALLBACK"


def load_scan_state(outdir):
    path = Path(outdir) / STATE_FILENAME
    if path.exists():
        try:
            with open(path, "r", encoding="utf-8") as f:
                return json.load(f)
        except Exception:
            return {}
    return {}


def save_scan_state(outdir, state):
    path = Path(outdir) / STATE_FILENAME
    try:
        with open(path, "w", encoding="utf-8") as f:
            json.dump(state, f, indent=2)
    except Exception as e:
        print(f"⚠️ Could not save scan state ({e}).")


def filter_unseen_targets(targets, state, cooldown_days, force_rescan):
    if force_rescan or cooldown_days <= 0:
        return targets
    cutoff = time.time() - cooldown_days * 86400.0
    fresh, skipped = [], 0
    for t in targets:
        last = state.get(t, {}).get("last_scanned_epoch", 0)
        (fresh if last < cutoff else []).append(t) if last < cutoff else None
        if last < cutoff:
            pass
        else:
            skipped += 1
    fresh = [t for t in targets if state.get(t, {}).get("last_scanned_epoch", 0) < cutoff]
    skipped = len(targets) - len(fresh)
    if skipped:
        print(f"⏭️  Skipping {skipped} target(s) already scanned within the last "
              f"{cooldown_days:.0f} day(s). Use --force-rescan to override.")
    return fresh


def send_webhook(url, payload):
    try:
        data = json.dumps(payload).encode("utf-8")
        req = _urllib_request.Request(url, data=data, headers={"Content-Type": "application/json"})
        _urllib_request.urlopen(req, timeout=10)
    except Exception as e:
        print(f"⚠️ Webhook notification failed: {e}")


def scan_one_target(target, args):
    """Runs one target as an isolated subprocess (this same script, `scan`
    subcommand) so a real, enforced timeout can kill a hung download - a
    thread-based timeout cannot do that in Python."""
    cmd = [
        sys.executable, "-u", THIS_SCRIPT, "scan",
        "--star", target, "--outdir", args.outdir,
        "--min-period", str(args.min_period), "--max-period", str(args.max_period),
        "--snr-threshold", str(args.snr_threshold), "--max-sectors", str(args.max_sectors),
        "--plot-mode", args.plot_mode,
    ]
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=args.per_target_timeout)
        stdout = result.stdout or ""
    except subprocess.TimeoutExpired:
        print(f"⏰ {target} exceeded --per-target-timeout ({args.per_target_timeout}s); killed.")
        return target, "FAILED"
    except Exception as e:
        print(f"❌ {target} subprocess raised: {e}")
        return target, "FAILED"

    status = "FAILED"
    for line in stdout.splitlines():
        if line.startswith("RESULT_STATUS::"):
            parts = line.split("::")
            if len(parts) == 3 and parts[1] == target:
                status = parts[2]
    return target, status


def _handle_afk_result(future, target, summary_counts, candidates_found, state, args):
    try:
        _, status = future.result()
    except Exception as e:
        print(f"❌ {target} raised while collecting result: {e}")
        status = "FAILED"
    if status not in summary_counts:
        status = "FAILED"
    summary_counts[status] += 1
    state[target] = {"last_scanned_epoch": time.time(), "last_status": status}

    if status == "PASS_CANDIDATE":
        candidates_found.append(target)
        print(f"🎉 PASS_CANDIDATE: {target}")
        if args.webhook_url:
            send_webhook(args.webhook_url, {
                "target": target, "status": status,
                "timestamp": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            })
    else:
        print(f"   [{status}] {target}")


def run_afk_session(args):
    """The full AFK batch-scan orchestration (parallel workers, time budget,
    resumable state, summary report)."""
    exofop_df = None
    if args.source == "TOI":
        targets, exofop_df = fetch_exofop_toi_list()
        target_source_label = "EXOFOP_TOI"
    else:
        try:
            _, exofop_df = fetch_exofop_toi_list()
        except Exception:
            exofop_df = None
        targets, target_source_label = fetch_mast_tic_targets(
            count=args.max_targets, max_mag=args.max_mag, sector=args.sector,
            exclude_known_tois=args.exclude_known_tois, exofop_df=exofop_df,
        )

    if not targets:
        print("❌ No targets found to scan.")
        sys.exit(1)
    targets = targets[:args.max_targets]

    state = load_scan_state(args.outdir)
    targets = filter_unseen_targets(targets, state, args.rescan_after_days, args.force_rescan)

    if not targets:
        print("\n💡 Nothing new to scan - every candidate target was already scanned within the "
              f"cooldown window ({args.rescan_after_days:.0f} days). Use --force-rescan to override.")
        return

    print(f"\n{'='*70}\n🤖 STARTING AFK AUTOMATION SESSION\n{'='*70}")
    print(f"⏱️  Time Limit:      {args.time_limit} minutes")
    print(f"🎯 Target Source:   {args.source} ({target_source_label})")
    print(f"🎯 Targets Queued:  {len(targets)}")
    print(f"⚙️  Parallel Workers: {args.workers}")
    print(f"⏰ Per-Target Timeout: {args.per_target_timeout}s")
    print(f"📁 Output Dir:      {args.outdir}\n{'='*70}\n")

    start_time = time.time()
    time_limit_seconds = args.time_limit * 60.0
    summary_counts = {"PASS_CANDIDATE": 0, "FALSE_POSITIVE_EB": 0, "WEAK_SIGNAL": 0, "FAILED": 0}
    candidates_found = []
    completed = 0

    with ThreadPoolExecutor(max_workers=args.workers) as executor:
        futures = {}
        target_iter = iter(targets)

        def submit_next():
            try:
                t = next(target_iter)
            except StopIteration:
                return False
            futures[executor.submit(scan_one_target, t, args)] = t
            return True

        for _ in range(args.workers):
            if not submit_next():
                break

        while futures:
            elapsed = time.time() - start_time
            if elapsed >= time_limit_seconds:
                print(f"\n⏰ Time limit of {args.time_limit} minutes reached! "
                      f"Letting {len(futures)} in-flight target(s) finish, submitting no more.")
                for future in as_completed(list(futures.keys())):
                    target = futures.pop(future)
                    _handle_afk_result(future, target, summary_counts, candidates_found, state, args)
                    completed += 1
                break

            for future in as_completed(list(futures.keys()), timeout=max(1, time_limit_seconds - elapsed)):
                target = futures.pop(future)
                _handle_afk_result(future, target, summary_counts, candidates_found, state, args)
                completed += 1
                submit_next()
                break

    save_scan_state(args.outdir, state)

    total_elapsed_mins = (time.time() - start_time) / 60.0
    print(f"\n{'='*70}\n🏆 AFK SCANNING SESSION COMPLETED in {total_elapsed_mins:.1f} minutes!\n{'='*70}")
    print(f"  - Total Stars Scanned:      {completed}")
    print(f"  - Candidates Discovered:    {summary_counts['PASS_CANDIDATE']} 🎉")
    print(f"  - False Positives (EB):     {summary_counts['FALSE_POSITIVE_EB']}")
    print(f"  - Weak Signals / Noise:     {summary_counts['WEAK_SIGNAL']}")
    print(f"  - Failed / No Data:         {summary_counts['FAILED']}")

    if candidates_found:
        print(f"\n🌟 NEW PASS CANDIDATES DISCOVERED ({len(candidates_found)}):")
        try:
            cat_df = pd.read_csv(Path(args.outdir) / "candidate_catalog.csv")
            for cand in candidates_found:
                match = cat_df[(cat_df["Target"] == cand) & (cat_df["Vetting_Status"] == "PASS_CANDIDATE")]
                if len(match) > 0:
                    row = match.iloc[-1]
                    print(f"   👉 {cand}: Period={row.get('Period_days', 'N/A')}d, "
                          f"Size={row.get('Planet_Radius_Rearth', 'N/A')} R_earth "
                          f"({row.get('Planet_Type', 'N/A')}), ExoFOP={row.get('Exofop_Status', 'N/A')}")
                else:
                    print(f"   👉 {cand}")
        except Exception:
            for cand in candidates_found:
                print(f"   👉 {cand}")
    else:
        print("\n💡 No strong new candidates passed vetting in this session. Keep scanning. 🔭")

    print(f"\n📁 Catalogs exported to: {Path(args.outdir) / 'candidate_catalog.csv'}")
    print(f"🖼️  Diagnostic plots in:  {Path(args.outdir) / 'plots'}")
    print(f"🗒️  Scan state saved to:  {Path(args.outdir) / STATE_FILENAME}\n")


# ==========================================================================
# GUI
# ==========================================================================
BG = "#0d1117"; BG2 = "#161b22"; ACCENT = "#58a6ff"; ACCENT2 = "#3fb950"
WARN = "#f85149"; YELLOW = "#e3b341"; FG = "#c9d1d9"; FG_DIM = "#8b949e"
FONT = ("Segoe UI", 10); FONT_SM = ("Segoe UI", 9); FONT_MONO = ("Consolas", 9)


def open_path_cross_platform(path):
    path = str(path)
    if sys.platform.startswith("win"):
        os.startfile(path)  # noqa: only reached on Windows
    elif sys.platform == "darwin":
        subprocess.run(["open", path], check=False)
    else:
        subprocess.run(["xdg-open", path], check=False)


def launch_gui():
    import tkinter as tk
    from tkinter import ttk, scrolledtext, filedialog, messagebox

    class ExoplanetHunterApp(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("🪐 Exoplanet Hunter")
            self.geometry("1050x720")
            self.minsize(900, 620)
            self.configure(bg=BG)

            self._scan_process = None
            self._scan_thread = None
            self._running = False
            self._sort_state = {}
            self._all_rows = []
            self._row_plot_paths = {}

            self._build_ui()
            self._check_dependencies()

        # ---- build ----
        def _build_ui(self):
            hdr = tk.Frame(self, bg=BG, pady=10)
            hdr.pack(fill=tk.X, padx=18)
            tk.Label(hdr, text="🪐 Exoplanet Hunter", font=("Segoe UI", 16, "bold"), bg=BG, fg=ACCENT).pack(side=tk.LEFT)
            tk.Label(hdr, text="Automated TESS Discovery & Vetting Pipeline", font=FONT, bg=BG, fg=FG_DIM).pack(side=tk.LEFT, padx=14, pady=4)
            self._status_label = tk.Label(hdr, text="● Ready", font=FONT_SM, bg=BG, fg=ACCENT2)
            self._status_label.pack(side=tk.RIGHT, padx=4)
            ttk.Separator(self, orient=tk.HORIZONTAL).pack(fill=tk.X, padx=14)

            body = tk.Frame(self, bg=BG)
            body.pack(fill=tk.BOTH, expand=True, padx=14, pady=8)
            left_outer = tk.Frame(body, bg=BG2, relief=tk.FLAT, bd=0)
            left_outer.pack(side=tk.LEFT, fill=tk.Y, ipadx=4)
            left_outer.pack_propagate(False)
            left_outer.configure(width=310)
            right_panel = tk.Frame(body, bg=BG, padx=10, pady=0)
            right_panel.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

            self._build_config_panel(left_outer)
            self._build_log_panel(right_panel)
            self._build_status_bar()

        def _build_config_panel(self, outer_parent):
            canvas = tk.Canvas(outer_parent, bg=BG2, highlightthickness=0)
            vsb = ttk.Scrollbar(outer_parent, orient="vertical", command=canvas.yview)
            canvas.configure(yscrollcommand=vsb.set)
            vsb.pack(side=tk.RIGHT, fill=tk.Y)
            canvas.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

            parent = tk.Frame(canvas, bg=BG2, padx=14, pady=12)
            window_id = canvas.create_window((0, 0), window=parent, anchor="nw")
            parent.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
            canvas.bind("<Configure>", lambda e: canvas.itemconfigure(window_id, width=e.width))
            canvas.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(int(-1 * (e.delta / 120)), "units"))

            tk.Label(parent, text="TARGET", font=("Segoe UI", 9, "bold"), bg=BG2, fg=FG_DIM).pack(anchor=tk.W, pady=(4, 2))
            self._mode = tk.StringVar(value="single")
            modes = tk.Frame(parent, bg=BG2); modes.pack(fill=tk.X, pady=2)
            for label, val in [("Single Star", "single"), ("Target File", "file"), ("AFK Auto", "afk")]:
                tk.Radiobutton(modes, text=label, variable=self._mode, value=val, command=self._on_mode_change,
                               bg=BG2, fg=FG, selectcolor=BG, activebackground=BG2, font=FONT_SM).pack(side=tk.LEFT, padx=2)

            self._star_frame = tk.Frame(parent, bg=BG2)
            self._star_frame.pack(fill=tk.X, pady=3)
            tk.Label(self._star_frame, text="Star / TIC ID:", font=FONT_SM, bg=BG2, fg=FG_DIM).pack(anchor=tk.W)
            self._star_var = tk.StringVar(value="TOI-700")
            tk.Entry(self._star_frame, textvariable=self._star_var, font=FONT_MONO, bg="#21262d", fg=FG,
                     insertbackground=FG, relief=tk.FLAT, bd=4).pack(fill=tk.X, pady=2)

            self._file_frame = tk.Frame(parent, bg=BG2)
            tk.Label(self._file_frame, text="Target File:", font=FONT_SM, bg=BG2, fg=FG_DIM).pack(anchor=tk.W)
            file_row = tk.Frame(self._file_frame, bg=BG2); file_row.pack(fill=tk.X)
            self._file_var = tk.StringVar(value="sample_targets.txt")
            tk.Entry(file_row, textvariable=self._file_var, font=FONT_MONO, bg="#21262d", fg=FG,
                     insertbackground=FG, relief=tk.FLAT, bd=4).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(0, 4))
            tk.Button(file_row, text="Browse", font=FONT_SM, bg="#21262d", fg=ACCENT, relief=tk.FLAT,
                      command=self._browse_file).pack(side=tk.RIGHT)

            self._afk_frame = tk.Frame(parent, bg=BG2)
            tk.Label(self._afk_frame, text="AFK Data Source:", font=FONT_SM, bg=BG2, fg=FG_DIM).pack(anchor=tk.W)
            self._afk_source = tk.StringVar(value="TIC")
            src_row = tk.Frame(self._afk_frame, bg=BG2); src_row.pack(fill=tk.X, pady=(0, 4))
            for val in ["TIC", "TOI"]:
                tk.Radiobutton(src_row, text=val, variable=self._afk_source, value=val, bg=BG2, fg=FG,
                               selectcolor=BG, activebackground=BG2, font=FONT_SM).pack(side=tk.LEFT, padx=4)

            afk_fields = [("Time Limit (mins)", "time_limit", "60"), ("Max Targets", "max_targets", "100"),
                          ("Parallel Workers", "workers", "3"), ("Rescan Cooldown (days)", "rescan_days", "14"),
                          ("Max Magnitude (opt.)", "max_mag", "")]
            self._afk_vars = {}
            for label, key, default in afk_fields:
                row = tk.Frame(self._afk_frame, bg=BG2); row.pack(fill=tk.X, pady=2)
                tk.Label(row, text=label + ":", font=FONT_SM, bg=BG2, fg=FG_DIM, width=19, anchor=tk.W).pack(side=tk.LEFT)
                var = tk.StringVar(value=default)
                tk.Entry(row, textvariable=var, font=FONT_MONO, bg="#21262d", fg=FG, insertbackground=FG,
                         relief=tk.FLAT, bd=4, width=8).pack(side=tk.LEFT, padx=4)
                self._afk_vars[key] = var

            self._force_rescan_var = tk.BooleanVar(value=False)
            tk.Checkbutton(self._afk_frame, text="Force rescan (ignore cooldown)", variable=self._force_rescan_var,
                           bg=BG2, fg=FG, selectcolor=BG, activebackground=BG2, font=FONT_SM).pack(anchor=tk.W, pady=(2, 4))

            webhook_row = tk.Frame(self._afk_frame, bg=BG2); webhook_row.pack(fill=tk.X, pady=2)
            tk.Label(webhook_row, text="Webhook URL (opt.):", font=FONT_SM, bg=BG2, fg=FG_DIM, anchor=tk.W).pack(anchor=tk.W)
            self._webhook_var = tk.StringVar(value="")
            tk.Entry(webhook_row, textvariable=self._webhook_var, font=FONT_MONO, bg="#21262d", fg=FG,
                     insertbackground=FG, relief=tk.FLAT, bd=4).pack(fill=tk.X, pady=2)

            ttk.Separator(parent, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
            tk.Label(parent, text="PARAMETERS", font=("Segoe UI", 9, "bold"), bg=BG2, fg=FG_DIM).pack(anchor=tk.W, pady=(0, 4))
            params = [("Max Sectors", "max_sectors", "2"), ("Min Period (days)", "min_period", "0.5"),
                      ("Max Period (days)", "max_period", "30.0"), ("SNR Threshold", "snr_threshold", "6.0")]
            self._params = {}
            for label, key, default in params:
                row = tk.Frame(parent, bg=BG2); row.pack(fill=tk.X, pady=2)
                tk.Label(row, text=label + ":", font=FONT_SM, bg=BG2, fg=FG_DIM, width=17, anchor=tk.W).pack(side=tk.LEFT)
                var = tk.StringVar(value=default)
                tk.Entry(row, textvariable=var, font=FONT_MONO, bg="#21262d", fg=FG, insertbackground=FG,
                         relief=tk.FLAT, bd=4, width=8).pack(side=tk.LEFT, padx=4)
                self._params[key] = var

            plot_row = tk.Frame(parent, bg=BG2); plot_row.pack(fill=tk.X, pady=2)
            tk.Label(plot_row, text="Plot Mode:", font=FONT_SM, bg=BG2, fg=FG_DIM, width=17, anchor=tk.W).pack(side=tk.LEFT)
            self._plot_mode = tk.StringVar(value="interesting")
            ttk.Combobox(plot_row, textvariable=self._plot_mode, values=["all", "interesting", "pass_only"],
                         state="readonly", width=11, font=FONT_SM).pack(side=tk.LEFT, padx=4)

            ttk.Separator(parent, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
            tk.Label(parent, text="OUTPUT", font=("Segoe UI", 9, "bold"), bg=BG2, fg=FG_DIM).pack(anchor=tk.W, pady=(0, 4))
            out_row = tk.Frame(parent, bg=BG2); out_row.pack(fill=tk.X, pady=2)
            tk.Label(out_row, text="Output Dir:", font=FONT_SM, bg=BG2, fg=FG_DIM, width=10, anchor=tk.W).pack(side=tk.LEFT)
            self._outdir_var = tk.StringVar(value="scan_results")
            tk.Entry(out_row, textvariable=self._outdir_var, font=FONT_MONO, bg="#21262d", fg=FG, insertbackground=FG,
                     relief=tk.FLAT, bd=4).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=(4, 4))
            tk.Button(out_row, text="…", font=FONT_SM, bg="#21262d", fg=ACCENT, relief=tk.FLAT,
                      command=self._browse_outdir, width=3).pack(side=tk.RIGHT)

            ttk.Separator(parent, orient=tk.HORIZONTAL).pack(fill=tk.X, pady=10)
            self._run_btn = tk.Button(parent, text="▶  START SCAN", font=("Segoe UI", 11, "bold"), bg=ACCENT2,
                                       fg="#0d1117", relief=tk.FLAT, command=self._toggle_scan, padx=8, pady=6)
            self._run_btn.pack(fill=tk.X, pady=2)
            tk.Button(parent, text="🗂 Open Results Folder", font=FONT_SM, bg="#21262d", fg=ACCENT, relief=tk.FLAT,
                      command=self._open_results).pack(fill=tk.X, pady=2)
            self._on_mode_change()

        def _build_log_panel(self, parent):
            top = tk.Frame(parent, bg=BG); top.pack(fill=tk.BOTH, expand=True)
            self._notebook = ttk.Notebook(top); self._notebook.pack(fill=tk.BOTH, expand=True)
            style = ttk.Style()
            style.configure("TNotebook", background=BG, borderwidth=0)
            style.configure("TNotebook.Tab", background=BG2, foreground=FG, padding=(10, 4))
            style.map("TNotebook.Tab", background=[("selected", BG)], foreground=[("selected", ACCENT)])

            log_tab = tk.Frame(self._notebook, bg=BG)
            self._notebook.add(log_tab, text="  📡 Live Log  ")
            self._log_box = scrolledtext.ScrolledText(log_tab, font=FONT_MONO, bg="#0d1117", fg=FG, relief=tk.FLAT,
                                                       state=tk.DISABLED, wrap=tk.WORD, padx=8, pady=6, insertbackground=FG)
            self._log_box.pack(fill=tk.BOTH, expand=True)
            self._log_box.tag_config("ok", foreground=ACCENT2)
            self._log_box.tag_config("warn", foreground=WARN)
            self._log_box.tag_config("accent", foreground=ACCENT)
            self._log_box.tag_config("yellow", foreground=YELLOW)
            self._log_box.tag_config("dim", foreground=FG_DIM)

            res_tab = tk.Frame(self._notebook, bg=BG)
            self._notebook.add(res_tab, text="  🔭 Candidates  ")
            self._build_results_table(res_tab)

        def _build_results_table(self, parent):
            search_row = tk.Frame(parent, bg=BG); search_row.pack(fill=tk.X, padx=4, pady=(4, 0))
            tk.Label(search_row, text="🔎", font=FONT_SM, bg=BG, fg=FG_DIM).pack(side=tk.LEFT)
            self._search_var = tk.StringVar(value="")
            self._search_var.trace_add("write", lambda *_: self._apply_filter())
            tk.Entry(search_row, textvariable=self._search_var, font=FONT_MONO, bg="#21262d", fg=FG,
                     insertbackground=FG, relief=tk.FLAT, bd=4).pack(side=tk.LEFT, fill=tk.X, expand=True, padx=6)
            tk.Label(search_row, text="Double-click a row to open its diagnostic plot", font=FONT_SM, bg=BG, fg=FG_DIM).pack(side=tk.RIGHT)

            cols = ("Target", "TIC ID", "Label", "Period (days)", "Depth %", "Rp (R⊕)", "Type", "ExoFOP", "Status", "SNR")
            self._tree = ttk.Treeview(parent, columns=cols, show="headings", selectmode="browse")
            self._tree.pack(fill=tk.BOTH, expand=True, padx=4, pady=4)
            self._tree.bind("<Double-1>", self._on_tree_double_click)

            style = ttk.Style()
            try:
                style.theme_use("clam")
            except tk.TclError:
                pass
            style.configure("Treeview", background=BG2, fieldbackground=BG2, foreground=FG, rowheight=24, font=FONT_SM)
            style.configure("Treeview.Heading", background="#21262d", foreground=ACCENT, font=("Segoe UI", 9, "bold"))
            style.map("Treeview", background=[("selected", "#264f78")])

            widths = [150, 85, 55, 95, 70, 80, 130, 170, 110, 55]
            for col, w in zip(cols, widths):
                self._tree.heading(col, text=col, command=lambda c=col: self._sort_col(c))
                self._tree.column(col, width=w, anchor=tk.CENTER if col != "Target" else tk.W)

            self._tree.tag_configure("pass", foreground=ACCENT2)
            self._tree.tag_configure("new", foreground=YELLOW)
            self._tree.tag_configure("fail", foreground=WARN)
            self._tree.tag_configure("weak", foreground=FG_DIM)

            scroll = ttk.Scrollbar(parent, orient=tk.HORIZONTAL, command=self._tree.xview)
            self._tree.configure(xscrollcommand=scroll.set)
            scroll.pack(fill=tk.X, padx=4)

            clear_row = tk.Frame(parent, bg=BG); clear_row.pack(fill=tk.X, padx=4, pady=4)
            tk.Button(clear_row, text="🔄 Refresh", font=FONT_SM, bg="#21262d", fg=ACCENT, relief=tk.FLAT,
                      command=self._refresh_results).pack(side=tk.LEFT, padx=4)
            tk.Button(clear_row, text="🗑 Clear Table", font=FONT_SM, bg="#21262d", fg=WARN, relief=tk.FLAT,
                      command=self._clear_results).pack(side=tk.LEFT, padx=4)

        def _build_status_bar(self):
            bar = tk.Frame(self, bg="#21262d", height=28); bar.pack(fill=tk.X, side=tk.BOTTOM)
            self._progress = ttk.Progressbar(bar, mode="indeterminate", length=180)
            self._progress.pack(side=tk.RIGHT, padx=10, pady=4)
            self._bottom_label = tk.Label(bar, text="Ready to scan.", font=FONT_SM, bg="#21262d", fg=FG_DIM)
            self._bottom_label.pack(side=tk.LEFT, padx=10)

        # ---- logic ----
        def _on_mode_change(self):
            mode = self._mode.get()
            self._star_frame.pack_forget(); self._file_frame.pack_forget(); self._afk_frame.pack_forget()
            if mode == "single":
                self._star_frame.pack(fill=tk.X, pady=3)
            elif mode == "file":
                self._file_frame.pack(fill=tk.X, pady=3)
            else:
                self._afk_frame.pack(fill=tk.X, pady=3)

        def _browse_file(self):
            path = filedialog.askopenfilename(title="Select Target List",
                                               filetypes=[("Text files", "*.txt"), ("CSV files", "*.csv"), ("All", "*.*")])
            if path:
                self._file_var.set(path)

        def _browse_outdir(self):
            path = filedialog.askdirectory(title="Select Output Directory")
            if path:
                self._outdir_var.set(path)

        def _validate_inputs(self):
            errors = []

            def check_float(label, value):
                try:
                    float(value)
                except ValueError:
                    errors.append(f"{label} must be a number (got '{value}')")

            def check_int(label, value):
                try:
                    int(value)
                except ValueError:
                    errors.append(f"{label} must be a whole number (got '{value}')")

            check_int("Max Sectors", self._params["max_sectors"].get())
            check_float("Min Period", self._params["min_period"].get())
            check_float("Max Period", self._params["max_period"].get())
            check_float("SNR Threshold", self._params["snr_threshold"].get())

            mode = self._mode.get()
            if mode == "single" and not self._star_var.get().strip():
                errors.append("Star / TIC ID cannot be empty")
            elif mode == "file":
                if not self._file_var.get().strip():
                    errors.append("Target file path cannot be empty")
                elif not Path(self._file_var.get()).exists():
                    errors.append(f"Target file not found: {self._file_var.get()}")
            elif mode == "afk":
                check_int("Time Limit", self._afk_vars["time_limit"].get())
                check_int("Max Targets", self._afk_vars["max_targets"].get())
                check_int("Parallel Workers", self._afk_vars["workers"].get())
                check_float("Rescan Cooldown", self._afk_vars["rescan_days"].get())
                max_mag = self._afk_vars["max_mag"].get().strip()
                if max_mag:
                    check_float("Max Magnitude", max_mag)
            return errors

        def _build_command(self):
            mode = self._mode.get()
            max_s = self._params["max_sectors"].get()
            min_p = self._params["min_period"].get()
            max_p = self._params["max_period"].get()
            snr = self._params["snr_threshold"].get()
            outdir = self._outdir_var.get()
            plot_mode = self._plot_mode.get()
            base = [sys.executable, "-u", THIS_SCRIPT]

            if mode == "afk":
                cmd = base + [
                    "afk", "--source", self._afk_source.get(),
                    "--time-limit", self._afk_vars["time_limit"].get(),
                    "--max-targets", self._afk_vars["max_targets"].get(),
                    "--workers", self._afk_vars["workers"].get(),
                    "--rescan-after-days", self._afk_vars["rescan_days"].get(),
                    "--outdir", outdir, "--max-sectors", max_s, "--min-period", min_p,
                    "--max-period", max_p, "--snr-threshold", snr, "--plot-mode", plot_mode,
                ]
                if self._force_rescan_var.get():
                    cmd.append("--force-rescan")
                max_mag = self._afk_vars["max_mag"].get().strip()
                if max_mag:
                    cmd += ["--max-mag", max_mag]
                webhook = self._webhook_var.get().strip()
                if webhook:
                    cmd += ["--webhook-url", webhook]
                return cmd
            elif mode == "single":
                return base + ["scan", "--star", self._star_var.get(), "--outdir", outdir,
                               "--max-sectors", max_s, "--min-period", min_p, "--max-period", max_p,
                               "--snr-threshold", snr, "--plot-mode", plot_mode]
            else:
                return base + ["scan", "--targets", self._file_var.get(), "--outdir", outdir,
                               "--max-sectors", max_s, "--min-period", min_p, "--max-period", max_p,
                               "--snr-threshold", snr, "--plot-mode", plot_mode]

        def _toggle_scan(self):
            self._stop_scan() if self._running else self._start_scan()

        def _start_scan(self):
            errors = self._validate_inputs()
            if errors:
                messagebox.showerror("Fix these before starting", "\n".join(f"• {e}" for e in errors))
                return
            self._running = True
            self._run_btn.configure(text="⏹  STOP SCAN", bg=WARN, fg="white")
            self._status_label.configure(text="● Running", fg=YELLOW)
            self._progress.start(12)
            self._clear_log()
            self._log("🚀 Starting scan session...\n", "accent")
            self._bottom_label.configure(text="Scanning...")
            cmd = self._build_command()
            cwd = str(Path(__file__).parent)
            self._scan_thread = threading.Thread(target=self._run_subprocess, args=(cmd, cwd), daemon=True)
            self._scan_thread.start()

        def _stop_scan(self):
            if self._scan_process:
                try:
                    self._scan_process.terminate()
                except Exception:
                    pass
            self._on_scan_done(aborted=True)

        def _run_subprocess(self, cmd, cwd):
            try:
                self._scan_process = subprocess.Popen([str(c) for c in cmd], cwd=cwd, stdout=subprocess.PIPE,
                                                        stderr=subprocess.STDOUT, text=True, bufsize=1,
                                                        encoding="utf-8", errors="replace")
                for line in self._scan_process.stdout:
                    self.after(0, self._append_log_line, line)
                self._scan_process.wait()
            except Exception as e:
                self.after(0, self._log, f"\n❌ Process error: {e}\n", "warn")
            finally:
                self.after(0, self._on_scan_done)

        def _on_scan_done(self, aborted=False):
            self._running = False
            self._scan_process = None
            self._run_btn.configure(text="▶  START SCAN", bg=ACCENT2, fg="#0d1117")
            self._progress.stop()
            if aborted:
                self._log("\n⏹ Scan stopped by user.\n", "warn")
                self._status_label.configure(text="● Stopped", fg=WARN)
                self._bottom_label.configure(text="Scan stopped.")
            else:
                self._log("\n✅ Scan session complete!\n", "ok")
                self._status_label.configure(text="● Idle", fg=ACCENT2)
                self._bottom_label.configure(text="Scan complete. Check Candidates tab.")
                self._refresh_results()
                self._notebook.select(1)

        def _append_log_line(self, line):
            stripped = line.rstrip()
            if "❌" in stripped or "Error" in stripped or "error" in stripped:
                tag = "warn"
            elif "✅" in stripped or "✓" in stripped or "PASS" in stripped:
                tag = "ok"
            elif "🔍" in stripped or "📡" in stripped or "⬇️" in stripped or "🌐" in stripped:
                tag = "accent"
            elif "⚠️" in stripped or "WEAK" in stripped or "LOW" in stripped:
                tag = "yellow"
            elif "====" in stripped or "COMPLETE" in stripped or "BATCH" in stripped:
                tag = "accent"
            else:
                tag = None
            self._log(line, tag)

        def _log(self, text, tag=None):
            self._log_box.configure(state=tk.NORMAL)
            self._log_box.insert(tk.END, text, tag) if tag else self._log_box.insert(tk.END, text)
            self._log_box.see(tk.END)
            self._log_box.configure(state=tk.DISABLED)

        def _clear_log(self):
            self._log_box.configure(state=tk.NORMAL)
            self._log_box.delete("1.0", tk.END)
            self._log_box.configure(state=tk.DISABLED)

        def _refresh_results(self):
            outdir = Path(self._outdir_var.get())
            csv_path = outdir / "candidate_catalog.csv"
            self._all_rows = []
            if not csv_path.exists():
                self._clear_results()
                return
            try:
                with open(csv_path, newline="", encoding="utf-8") as f:
                    for row in csv.DictReader(f):
                        self._all_rows.append(row)
            except Exception as e:
                self._log(f"\n⚠️ Could not load catalog: {e}\n", "warn")
                return
            self._apply_filter()

        def _apply_filter(self):
            query = self._search_var.get().strip().lower()
            self._clear_results()
            for row in self._all_rows:
                target = row.get("Target", ""); tic_id = row.get("TIC_ID", "")
                label = row.get("Planet_Label", ""); period = row.get("Period_days", "")
                depth = row.get("Depth_percent", ""); rp = row.get("Planet_Radius_Rearth", "")
                ptype = row.get("Planet_Type", ""); exofop = row.get("Exofop_Status", "")
                status = row.get("Vetting_Status", ""); snr = row.get("SNR", "")
                plot_path = row.get("Plot_Path", "")

                if query:
                    haystack = " ".join([target, tic_id, ptype, exofop, status]).lower()
                    if query not in haystack:
                        continue

                if status == "PASS_CANDIDATE" and "NEW_UNREPORTED" in exofop:
                    tag = "new"
                elif status == "PASS_CANDIDATE":
                    tag = "pass"
                elif "FALSE_POSITIVE" in status:
                    tag = "fail"
                else:
                    tag = "weak"

                item_id = self._tree.insert("", tk.END, values=(target, tic_id, label, period, depth, rp, ptype,
                                                                 exofop, status, snr), tags=(tag,))
                self._row_plot_paths[item_id] = plot_path

        def _clear_results(self):
            for row in self._tree.get_children():
                self._tree.delete(row)
            self._row_plot_paths = {}

        def _sort_col(self, col):
            ascending = not self._sort_state.get(col, False)
            self._sort_state = {col: ascending}
            items = [(self._tree.set(k, col), k) for k in self._tree.get_children("")]
            try:
                items.sort(key=lambda x: float(x[0]), reverse=not ascending)
            except ValueError:
                items.sort(key=lambda x: x[0].lower(), reverse=not ascending)
            for idx, (_, k) in enumerate(items):
                self._tree.move(k, "", idx)
            for c in self._tree["columns"]:
                self._tree.heading(c, text=c.split(" ▲")[0].split(" ▼")[0])
            self._tree.heading(col, text=col + (" ▲" if ascending else " ▼"))

        def _on_tree_double_click(self, event):
            item_id = self._tree.identify_row(event.y)
            if not item_id:
                return
            plot_path = self._row_plot_paths.get(item_id)
            if not plot_path or plot_path == "N/A":
                messagebox.showinfo("No Plot", "No diagnostic plot was generated for this row (depends on --plot-mode).")
                return
            p = Path(plot_path)
            if not p.exists():
                messagebox.showinfo("Not Found", f"Plot file not found on disk:\n{plot_path}")
                return
            try:
                open_path_cross_platform(str(p.resolve()))
            except Exception as e:
                messagebox.showerror("Could Not Open Plot", str(e))

        def _open_results(self):
            path = str(Path(self._outdir_var.get()).resolve())
            if not os.path.exists(path):
                messagebox.showinfo("Not Found", f"Output directory does not exist yet:\n{path}")
                return
            try:
                open_path_cross_platform(path)
            except Exception as e:
                messagebox.showerror("Could Not Open Folder", str(e))

        def _check_dependencies(self):
            def _check_packages():
                import_stmt = "import " + ", ".join(REQUIRED_PACKAGES)
                try:
                    result = subprocess.run([sys.executable, "-c", import_stmt], capture_output=True, text=True, timeout=30)
                    if result.returncode != 0:
                        err_lines = (result.stderr or "").strip().splitlines()
                        err_msg = err_lines[-1] if err_lines else "Unknown import error"
                        self.after(0, lambda: messagebox.showwarning(
                            "Missing Python Packages",
                            f"One or more required packages could not be imported:\n\n{err_msg}\n\n"
                            "Install them with:\npip install lightkurve numpy matplotlib scipy pandas astroquery"))
                except Exception:
                    pass
            threading.Thread(target=_check_packages, daemon=True).start()

    app = ExoplanetHunterApp()
    app.mainloop()


# ==========================================================================
# CLI
# ==========================================================================
def build_argparser():
    parser = argparse.ArgumentParser(
        prog="exoplanet_pipeline.py",
        description="Unified Automated Exoplanet Finder & Vetting Pipeline (engine + AFK scheduler + GUI, single file).",
    )
    sub = parser.add_subparsers(dest="command")

    p_scan = sub.add_parser("scan", help="Scan a single star or a target list file")
    p_scan.add_argument("--star", type=str, help="Single target star name or TIC ID")
    p_scan.add_argument("--targets", type=str, help="Path to a text file of targets (one per line)")
    p_scan.add_argument("--outdir", type=str, default="scan_results")
    p_scan.add_argument("--min-period", type=float, default=0.5)
    p_scan.add_argument("--max-period", type=float, default=30.0)
    p_scan.add_argument("--snr-threshold", type=float, default=6.0)
    p_scan.add_argument("--max-sectors", type=int, default=3)
    p_scan.add_argument("--plot-mode", type=str, choices=["all", "interesting", "pass_only"], default="interesting")
    p_scan.add_argument("--toi-cache-hours", type=float, default=24.0)
    p_scan.add_argument("--bls-target-points", type=int, default=100000)

    p_afk = sub.add_parser("afk", help="Run an AFK batch scan session across many targets")
    p_afk.add_argument("--source", type=str, choices=["TOI", "TIC"], default="TIC")
    p_afk.add_argument("--time-limit", type=int, default=60)
    p_afk.add_argument("--max-targets", type=int, default=100)
    p_afk.add_argument("--max-sectors", type=int, default=2)
    p_afk.add_argument("--min-period", type=float, default=0.5)
    p_afk.add_argument("--max-period", type=float, default=30.0)
    p_afk.add_argument("--snr-threshold", type=float, default=6.0)
    p_afk.add_argument("--plot-mode", type=str, choices=["all", "interesting", "pass_only"], default="interesting")
    p_afk.add_argument("--max-mag", type=float, default=None)
    p_afk.add_argument("--sector", type=int, default=None)
    p_afk.add_argument("--exclude-known-tois", action="store_true", default=True)
    p_afk.add_argument("--include-known-tois", dest="exclude_known_tois", action="store_false")
    p_afk.add_argument("--outdir", type=str, default="daily_afk_results")
    p_afk.add_argument("--workers", type=int, default=3)
    p_afk.add_argument("--per-target-timeout", type=int, default=600)
    p_afk.add_argument("--rescan-after-days", type=float, default=14.0)
    p_afk.add_argument("--force-rescan", action="store_true")
    p_afk.add_argument("--webhook-url", type=str, default=None)

    sub.add_parser("self-test", help="Run the offline regression suite (no network needed)")
    sub.add_parser("gui", help="Launch the graphical interface (default with no arguments)")

    return parser


def cmd_scan(args):
    targets_to_scan = []
    if args.star:
        targets_to_scan.append(args.star)
    elif args.targets:
        if os.path.exists(args.targets):
            with open(args.targets, "r", encoding="utf-8") as f:
                targets_to_scan = [line.strip() for line in f if line.strip() and not line.startswith("#")]
        else:
            print(f"❌ Targets file '{args.targets}' not found.")
            sys.exit(1)
    else:
        print("💡 No target specified. Defaulting to example candidate: TOI-700")
        targets_to_scan.append("TOI-700")

    finder = AutomatedExoplanetFinder(
        output_dir=args.outdir, min_period=args.min_period, max_period=args.max_period,
        snr_threshold=args.snr_threshold, max_sectors=args.max_sectors, plot_mode=args.plot_mode,
        toi_cache_hours=args.toi_cache_hours, bls_target_grid_points=args.bls_target_points,
    )

    print(f"\n🚀 Starting Professional Exoplanet Finder on {len(targets_to_scan)} target(s)...")
    summary_counts = {"PASS_CANDIDATE": 0, "FALSE_POSITIVE_EB": 0, "WEAK_SIGNAL": 0, "FAILED": 0}

    for idx, target in enumerate(targets_to_scan, 1):
        print(f"\n[{idx}/{len(targets_to_scan)}] Processing {target}...")
        status = finder.process_target(target)
        summary_counts[status if status in summary_counts else "FAILED"] += 1
        # Machine-readable marker: run_afk_session()/scan_one_target() greps
        # this out of captured subprocess stdout.
        print(f"RESULT_STATUS::{target}::{status}")

    print(f"\n{'='*70}\n✨ BATCH SCAN COMPLETE!\n{'='*70}")
    print(f"  - Candidates Discovered:    {summary_counts['PASS_CANDIDATE']}")
    print(f"  - False Positives (EB):     {summary_counts['FALSE_POSITIVE_EB']}")
    print(f"  - Weak Signals / Noise:    {summary_counts['WEAK_SIGNAL']}")
    print(f"  - Failed / No Data:         {summary_counts['FAILED']}")
    print(f"\n📁 Catalogs exported to: {finder.csv_path}")
    print(f"🖼️ Diagnostic plots in:  {finder.plots_dir}\n")


def main():
    if len(sys.argv) == 1:
        launch_gui()
        return

    parser = build_argparser()
    args = parser.parse_args()

    if args.command in (None, "gui"):
        launch_gui()
    elif args.command == "scan":
        if lk is None:
            print("❌ lightkurve is not installed. Run: pip install lightkurve numpy matplotlib scipy pandas astroquery")
            sys.exit(1)
        cmd_scan(args)
    elif args.command == "afk":
        run_afk_session(args)
    elif args.command == "self-test":
        if lk is None:
            print("❌ lightkurve is not installed. Run: pip install lightkurve numpy matplotlib scipy pandas astroquery")
            sys.exit(1)
        sys.exit(0 if run_self_test_suite() else 1)


if __name__ == "__main__":
    main()
