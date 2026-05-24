# CMB Power Spectrum Feature Analysis

This project analyzes the cosmic microwave background temperature TT power spectrum using Python. It reads a spectrum file, plots D_ell against multipole moment ell, marks the strongest acoustic peaks, and saves the measured peak locations to CSV.

The main goal is feature extraction, not cosmological parameter fitting. The script focuses on where the first few acoustic peaks appear and how those locations relate to angular scale through the approximation theta ≈ 180 / ell.

## What the script does

- Loads a CMB TT power spectrum file
- Smooths the spectrum with a Savitzky-Golay filter
- Detects the first three acoustic peaks using spacing and prominence thresholds
- Estimates angular scale from each peak location
- Saves peak measurements and a summary table
- Plots the spectrum, smoothed spectrum, log-scaled spectrum, and reported uncertainty when uncertainty columns are present

## Data

The analysis is set up for the Planck TT power spectrum text file:

    data/planck_tt_power_spectrum.txt

If the file is not present, the script stops instead of generating fake data. This keeps the outputs tied to the input spectrum used in the project.

## Run

Install the dependencies:

    python3 -m pip install numpy pandas matplotlib scipy

Run the analysis:

    python3 src/main.py

The script writes its figures and CSV files into `outputs/`.

## Outputs

Main figures:

- `outputs/cmb_power_spectrum_peaks.png`
- `outputs/cmb_power_spectrum_log.png`
- `outputs/cmb_smoothed_peak_detection.png`
- `outputs/cmb_uncertainty_by_ell.png`

Tables:

- `outputs/peak_table.csv`
- `outputs/spectrum_summary.csv`

The peak table includes the detected peak number, ell location, D_ell value, smoothed D_ell value, peak prominence, angular scale, and spacing from the previous detected peak. If reported uncertainty columns are present in the source file, the table also includes the mean D_ell error and the peak signal-to-uncertainty ratio.

## Scientific context

The first acoustic peak in the CMB temperature spectrum appears near ell ≈ 220, which corresponds to an angular scale near one degree. This scale is important because it reflects the apparent size of sound-horizon structure at recombination.

The second and third peaks appear at higher ell values. Their locations show smaller angular features in the early-universe plasma. This project only extracts peak positions and simple derived values. It does not fit a cosmological model.

## Notes

Peak detection depends on smoothing, prominence, and spacing choices. Those settings are useful for identifying the dominant acoustic peaks, but they are not a replacement for a full CMB likelihood analysis.
