# Changelog

All notable changes to this project are documented here.

---

## [1.0.0] - 2026-09-27

### Added
- Unified single-file pipeline (`exoplanet_pipeline.py`) merging core engine, AFK scheduler, and GUI
- `scan` CLI subcommand for single-star and batch-file scanning
- `afk` CLI subcommand for time-limited automated overnight sessions
- `self-test` subcommand for offline regression testing
- Dark-themed Tkinter GUI with live log, sortable results table, and plot viewer
- Real MAST timeseries query for TIC targets (replaces fake sequential ID generation)
- Multi-planet 2-pass BLS transit masking
- Physical planet radius calculation in Earth radii from TIC stellar parameters
- NASA ExoFOP TOI live cross-referencing (`KNOWN_TOI` vs `NEW_UNREPORTED_CANDIDATE`)
- Odd/even transit depth false-positive vetting (eclipsing binary filter)
- 5-panel publication-grade diagnostic PNG plots
- CSV + JSON candidate catalog export
- Cross-process file lock for safe parallel worker catalog writes
- AFK session state persistence (`scanned_state.json`) with configurable cooldown
- Webhook URL support for posting results to external services
- Adaptive BLS frequency factor (prevents grid explosion on multi-sector data)
- Per-target subprocess timeout to prevent stuck workers killing a session
- `.gitignore` excluding scan outputs, TESS cache, and bytecode

### Fixed
- `bls.get_in_transit_mask` → correct API `bls.get_transit_mask`
- `args.time-limit` argparse dash → underscore access bug
- BLS grid explosion (`ValueError: period contains 100M+ points`) via adaptive `frequency_factor`
- FITS cache corruption handled by per-sector download loop with fault isolation
- ExoFOP catalog download using `pd.read_csv(url)` instead of fragile `urllib`

---

## [0.x] — Pre-release development (private)

Earlier iterations included separate `auto_exoplanet_finder.py`, `run_daily_scan.py`,
and `gui_app.py` scripts that were merged and superseded by v1.0.0.
