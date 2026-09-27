<div align="center">

# 🪐 Exoplanet Detector

**Automated TESS Exoplanet Discovery & Vetting Pipeline**

[![Python](https://img.shields.io/badge/Python-3.9%2B-blue?logo=python&logoColor=white)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![TESS](https://img.shields.io/badge/Data-NASA%20TESS-orange?logo=nasa)](https://tess.mit.edu/)
[![ExoFOP](https://img.shields.io/badge/Catalog-ExoFOP%20TOI-blueviolet)](https://exofop.ipac.caltech.edu/tess/)
[![Lightkurve](https://img.shields.io/badge/Powered%20by-Lightkurve-red)](https://docs.lightkurve.org/)

*Find new exoplanet candidates from NASA TESS data — with a GUI, a CLI, and fully automated AFK scanning. Built for citizen scientists and researchers alike.*

</div>

---

## ✨ Features

| Feature | Description |
|---|---|
| **BLS Periodogram** | Box Least Squares transit search with adaptive frequency grid (no crash on large datasets) |
| **Multi-planet search** | 2-pass transit masking to detect sibling planets in the same system |
| **Physical planet sizes** | Calculates planet radius in Earth radii from stellar parameters |
| **ExoFOP cross-reference** | Automatically tags detections as `KNOWN_TOI` or `NEW_UNREPORTED_CANDIDATE` |
| **False-positive vetting** | Odd/even transit depth test to filter eclipsing binaries |
| **5-panel diagnostic plots** | Publication-style PNG reports per candidate |
| **AFK batch mode** | Time-limited sessions scanning 100s of real TESS targets overnight |
| **Real MAST queries** | Fetches only stars that TESS actually observed (no fake sequential IDs) |
| **GUI included** | Dark-themed Tkinter desktop app — no extra install needed |
| **CSV + JSON export** | Full catalog of all results with physical parameters |

---

## 📦 Installation

### Requirements

- Python 3.9 or newer
- Internet access (to download TESS data from MAST)

### Install dependencies

```bash
pip install lightkurve numpy matplotlib scipy pandas astroquery
```

### Clone the repository

```bash
git clone https://github.com/thekayrasari/exoplanet-detector.git
cd exoplanet-detector
```

---

## 🚀 Quick Start

### 1. Launch the GUI (easiest)

```bash
python exoplanet_pipeline.py
```

### 2. Scan a single star (CLI)

```bash
python exoplanet_pipeline.py scan --star "TOI-700"
python exoplanet_pipeline.py scan --star "TIC 149603524"
```

### 3. Scan a list of targets

Create `targets.txt` (one per line, `#` for comments):

```
# My favourite candidates
TOI-700
LHS 1140
TIC 149603524
```

Then run:

```bash
python exoplanet_pipeline.py scan --targets targets.txt --outdir results
```

### 4. AFK overnight session

```bash
# Scan up to 100 real TESS TIC stars over 60 minutes (3 parallel workers)
python exoplanet_pipeline.py afk --source TIC --time-limit 60 --workers 3

# Scan ExoFOP TOI candidates instead
python exoplanet_pipeline.py afk --source TOI --time-limit 120 --max-targets 200
```

### 5. Self-test (offline, no network needed)

```bash
python exoplanet_pipeline.py self-test
```

---

## ⚙️ CLI Reference

### `scan` subcommand

| Argument | Default | Description |
|---|---|---|
| `--star TEXT` | — | Single target (star name or TIC ID) |
| `--targets FILE` | — | Text file with one target per line |
| `--outdir DIR` | `scan_results` | Output directory |
| `--min-period DAYS` | `0.5` | Minimum orbital period to search |
| `--max-period DAYS` | `30.0` | Maximum orbital period to search |
| `--snr-threshold` | `6.0` | Minimum SNR to declare a candidate |
| `--max-sectors INT` | `3` | Max TESS sectors to download per target |
| `--plot-mode` | `interesting` | `all` / `interesting` / `pass_only` |
| `--bls-target-points` | `100000` | BLS frequency grid size |

### `afk` subcommand

| Argument | Default | Description |
|---|---|---|
| `--source` | `TIC` | `TIC` (MAST live query) or `TOI` (ExoFOP) |
| `--time-limit MINS` | `60` | Hard stop after N minutes |
| `--max-targets INT` | `100` | Max targets to attempt |
| `--workers INT` | `3` | Parallel scan workers |
| `--sector INT` | — | Restrict to a specific TESS sector |
| `--max-mag FLOAT` | — | Only stars brighter than this TESS magnitude |
| `--exclude-known-tois` | on | Skip already-catalogued TOIs (TIC mode) |
| `--rescan-after-days` | `14` | Cooldown before re-scanning the same target |
| `--force-rescan` | — | Ignore cooldown |
| `--webhook-url URL` | — | POST JSON results here when done |
| `--outdir DIR` | `daily_afk_results` | Output directory |
| `--per-target-timeout` | `600` | Seconds before a stuck target is killed |

---

## 📁 Output Structure

```
scan_results/
├── candidate_catalog.csv      # All results with physical parameters
├── candidate_summary.json     # Same data in JSON format
├── exofop_toi_cache.csv       # Cached ExoFOP TOI catalog (24h TTL)
├── scanned_state.json         # AFK session resume state
└── plots/
    ├── TIC_149603524_report.png
    └── ...
```

### `candidate_catalog.csv` columns

| Column | Description |
|---|---|
| `Target` | Star identifier |
| `Period_days` | Best-fit orbital period |
| `Depth_percent` | Transit depth (%) |
| `Duration_hours` | Transit duration |
| `Planet_Radius_Rearth` | Planet radius in Earth radii |
| `Planet_Type` | Rocky / Super-Earth / Sub-Neptune / Neptune-like / Gas Giant / Hot Jupiter |
| `SNR` | Signal-to-Noise ratio |
| `Odd_Even_Ratio` | >1.2 flags potential eclipsing binary |
| `Exofop_Status` | `KNOWN_TOI` / `NEW_UNREPORTED_CANDIDATE` / `EXOFOP_UNAVAILABLE` |
| `Vetting_Status` | `PASS_CANDIDATE` / `WEAK_SIGNAL` / `FALSE_POSITIVE_EB` |
| `Multi_Planet_Count` | Number of planet signals found in this system |

---

## 🔭 How It Works

```
1. Fetch TESS light curve(s) from MAST  (SPOC pipeline prioritised)
         ↓
2. Stitch & clean multi-sector data  (sigma-clip, normalise)
         ↓
3. BLS periodogram  (adaptive grid — no grid explosion)
         ↓
4. Pass 1 planet extracted & masked from light curve
         ↓
5. Pass 2 BLS on residuals → sibling planet search
         ↓
6. Vetting: SNR check + odd/even depth test (EB filter)
         ↓
7. Cross-reference NASA ExoFOP TOI catalog
         ↓
8. Generate 5-panel diagnostic plot
         ↓
9. Append to CSV / JSON catalog
```

---

## 🧠 TOI vs TIC — Which Should I Use?

| Source | What it is | Best for |
|---|---|---|
| **TOI** | Pre-flagged TESS Objects of Interest | Testing your setup; confirming known candidates |
| **TIC** | All TESS-observed stars | **Genuine discovery** — most have never been individually analysed |

> Use `--source TIC` for AFK sessions to maximise your chance of finding something new. The pipeline queries MAST for stars where TESS actually produced a timeseries, so you never waste time on targets with no data.

---

## 🤝 Contributing

Contributions are welcome! Please read [CONTRIBUTING.md](CONTRIBUTING.md) first.

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-new-feature`
3. Commit your changes: `git commit -m 'Add some feature'`
4. Push the branch: `git push origin feature/my-new-feature`
5. Open a Pull Request

---

## 📄 License

This project is licensed under the **MIT License** — see [LICENSE](LICENSE) for details.

---

## 🙏 Acknowledgements

- **[NASA TESS Mission](https://tess.mit.edu/)** — for the photometric data
- **[Lightkurve](https://docs.lightkurve.org/)** — Python library powering all TESS data access and BLS analysis
- **[NASA ExoFOP](https://exofop.ipac.caltech.edu/tess/)** — TOI cross-reference catalog
- **[MAST Archive](https://mast.stsci.edu/)** — for hosting TESS data
- **[Planet Hunters TESS](https://www.zooniverse.org/projects/nora-dot-eisner/planet-hunters-tess)** — the citizen science project that inspired this tool

---

<div align="center">
Made with ❤️ for the citizen science astronomy community
</div>
