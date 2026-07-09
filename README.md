# CMB Acoustic Peak Extraction

Feature extraction on the Planck cosmic microwave background TT power spectrum. The script reads a spectrum file, smooths it with a Savitzky-Golay filter, locates the first three acoustic peaks by prominence and spacing, converts each peak location to an angular scale via theta = 180 / ell, and writes the measurements to CSV with figures.

![Detected acoustic peaks](outputs/cmb_power_spectrum_peaks.png)

Measured on the included Planck TT spectrum: peaks at ell = 224, 514, and 822, first peak angular scale 0.80 degrees, mean peak spacing 299 in ell. The scope is peak measurement and derived scales, not cosmological parameter fitting.

## Run

```bash
python3 -m pip install numpy pandas matplotlib scipy
python3 src/main.py
```

The script looks for a spectrum file in `data/*.txt`, `data/*.csv`, then the repo root, and exits with an error if none is found. It ships with `data/planck_tt_power_spectrum.txt` (Planck TT, columns ell, D_ell, lower error, upper error). Any whitespace or comma delimited file with ell and D_ell columns works; reported uncertainty columns are picked up when present.

## Outputs

Written to `outputs/`:

- `peak_table.csv`: peak number, ell, D_ell, smoothed D_ell, prominence, angular scale, spacing from previous peak, and signal-to-uncertainty ratio when the source file reports errors
- `spectrum_summary.csv`: source file, peak count, first peak location and angular scale, mean peak spacing
- Figures: linear spectrum with labeled peaks, log-log spectrum, smoothed spectrum with detections, uncertainty by multipole

## Method notes

Peak detection runs on the smoothed spectrum with a prominence threshold of 0.20 standard deviations and a minimum separation of 120 in ell, then keeps the three strongest detections ordered by ell. The first peak near one degree traces the sound horizon at recombination; the spacing of the following peaks reflects the harmonic structure of the same acoustic oscillations. Smoothing window and threshold choices shift detected locations by a few multipoles, which sets the effective precision of the measurement.
