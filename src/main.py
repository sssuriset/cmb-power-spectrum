import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks, savgol_filter

def make_sample_data(path):
    ell = np.arange(2, 2501)

    peak1 = 5600 * np.exp(-0.5 * ((ell - 220) / 65) ** 2)
    peak2 = 2600 * np.exp(-0.5 * ((ell - 540) / 85) ** 2)
    peak3 = 2300 * np.exp(-0.5 * ((ell - 800) / 95) ** 2)
    damping = np.exp(-ell / 2400)

    baseline = 400 + 900 * np.exp(-ell / 900)
    noise = np.random.default_rng(4).normal(0, 110, size=len(ell))

    dl = (baseline + peak1 + peak2 + peak3) * damping + noise
    dl = np.clip(dl, 0, None)

    data = pd.DataFrame({
        "ell": ell,
        "D_ell": dl
    })

    data.to_csv(path, index=False)


def find_data_file():
    patterns = [
        "data/*.txt",
        "data/*.csv",
        "*.txt",
        "*.csv"
    ]

    files = []

    for pattern in patterns:
        files.extend(glob.glob(pattern))

    blocked_words = [
        "README",
        "requirements",
        "peak_table",
        "summary",
        "output",
        "result"
    ]

    files = [
        f for f in files
        if not any(word.lower() in os.path.basename(f).lower() for word in blocked_words)
    ]

    if files:
        return files[0]

    os.makedirs("data", exist_ok=True)
    sample_path = "data/sample_cmb_power_spectrum.csv"
    make_sample_data(sample_path)

    return sample_path


def load_spectrum():
    path = find_data_file()

    try:
        data = pd.read_csv(path, comment="#")
    except Exception:
        data = pd.read_csv(path, delim_whitespace=True, comment="#", header=None)

    # If pandas read the whole row as one column, retry as whitespace-separated.
    if len(data.columns) == 1:
        data = pd.read_csv(path, delim_whitespace=True, comment="#", header=None)

    # Normalize column names.
    data.columns = [str(col).strip() for col in data.columns]

    lower_cols = {str(col).lower(): col for col in data.columns}

    ell_candidates = ["ell", "l", "multipole", "#ell"]
    dl_candidates = ["d_ell", "dl", "d_l", "power", "tt", "cl", "dell"]

    ell_col = None
    dl_col = None

    for candidate in ell_candidates:
        if candidate in lower_cols:
            ell_col = lower_cols[candidate]
            break

    for candidate in dl_candidates:
        if candidate in lower_cols:
            dl_col = lower_cols[candidate]
            break

    if ell_col is not None and dl_col is not None:
        clean = data[[ell_col, dl_col]].copy()
        clean.columns = ["ell", "D_ell"]
    elif len(data.columns) >= 2:
        numeric = data.apply(pd.to_numeric, errors="coerce")
        numeric = numeric.dropna(axis=1, how="all")

        if len(numeric.columns) < 2:
            raise ValueError(
                f"Could not identify numeric ell and D_ell columns in {path}."
            )

        clean = numeric.iloc[:, :2].copy()
        clean.columns = ["ell", "D_ell"]
    else:
        raise ValueError(f"Could not identify ell and D_ell columns in {path}.")

    clean["ell"] = pd.to_numeric(clean["ell"], errors="coerce")
    clean["D_ell"] = pd.to_numeric(clean["D_ell"], errors="coerce")

    clean = clean.dropna()
    clean = clean.sort_values("ell")
    clean = clean[(clean["ell"] > 0) & (clean["D_ell"] > 0)]

    if clean.empty:
        raise ValueError(f"No valid positive ell and D_ell values found in {path}.")

    return clean, path


def smooth_spectrum(dl):
    window = 61

    if len(dl) < window:
        window = len(dl) // 2 * 2 - 1

    if window < 7:
        return dl

    return savgol_filter(dl, window_length=window, polyorder=3)


def acoustic_scale_degrees(ell):
    return 180 / ell


def detect_peaks(data):
    ell = data["ell"].to_numpy()
    dl = data["D_ell"].to_numpy()

    smoothed = smooth_spectrum(dl)

    peaks, properties = find_peaks(
        smoothed,
        prominence=np.std(smoothed) * 0.20,
        distance=120
    )

    peak_table = pd.DataFrame({
        "peak_number": np.arange(1, len(peaks) + 1),
        "ell": ell[peaks],
        "D_ell": dl[peaks],
        "smoothed_D_ell": smoothed[peaks],
        "prominence": properties["prominences"],
        "angular_scale_degrees": acoustic_scale_degrees(ell[peaks])
    })

    peak_table = peak_table.sort_values("smoothed_D_ell", ascending=False)
    peak_table = peak_table.head(3)
    peak_table = peak_table.sort_values("ell").reset_index(drop=True)
    peak_table["peak_number"] = np.arange(1, len(peak_table) + 1)
    peak_table["spacing_from_previous_ell"] = peak_table["ell"].diff()

    return peak_table, smoothed


def save_linear_plot(data, peak_table):
    plt.figure(figsize=(10, 5))
    plt.plot(data["ell"], data["D_ell"], linewidth=1, label="CMB power spectrum")
    plt.scatter(peak_table["ell"], peak_table["D_ell"], s=45, label="Detected acoustic peaks")

    for _, row in peak_table.iterrows():
        plt.annotate(
            f"Peak {int(row['peak_number'])}\nell={int(row['ell'])}",
            (row["ell"], row["D_ell"]),
            textcoords="offset points",
            xytext=(8, 8)
        )

    plt.xlabel("Multipole moment ell")
    plt.ylabel("D_ell")
    plt.title("CMB Temperature Power Spectrum with Acoustic Peaks")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/cmb_power_spectrum_peaks.png", dpi=300)
    plt.close()


def save_log_plot(data, peak_table):
    plt.figure(figsize=(10, 5))
    plt.plot(data["ell"], data["D_ell"], linewidth=1, label="CMB power spectrum")
    plt.scatter(peak_table["ell"], peak_table["D_ell"], s=45, label="Detected acoustic peaks")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Multipole moment ell")
    plt.ylabel("D_ell")
    plt.title("CMB Power Spectrum on Log Scale")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/cmb_power_spectrum_log.png", dpi=300)
    plt.close()


def save_smoothed_plot(data, smoothed, peak_table):
    plt.figure(figsize=(10, 5))
    plt.plot(data["ell"], data["D_ell"], alpha=0.45, linewidth=1, label="Raw spectrum")
    plt.plot(data["ell"], smoothed, linewidth=2, label="Smoothed spectrum")
    plt.scatter(peak_table["ell"], peak_table["smoothed_D_ell"], s=45, label="Detected peaks")
    plt.xlabel("Multipole moment ell")
    plt.ylabel("D_ell")
    plt.title("Smoothed CMB Spectrum and Peak Detection")
    plt.legend()
    plt.tight_layout()
    plt.savefig("outputs/cmb_smoothed_peak_detection.png", dpi=300)
    plt.close()


def main():
    os.makedirs("outputs", exist_ok=True)

    data, source_path = load_spectrum()
    peak_table, smoothed = detect_peaks(data)

    peak_table.to_csv("outputs/peak_table.csv", index=False)

    if len(peak_table) > 1:
        mean_spacing = peak_table["spacing_from_previous_ell"].dropna().mean()
    else:
        mean_spacing = np.nan

    summary = pd.DataFrame({
        "metric": [
            "source_file",
            "number_of_detected_peaks",
            "first_peak_ell",
            "first_peak_angular_scale_degrees",
            "mean_peak_spacing_ell"
        ],
        "value": [
            source_path,
            len(peak_table),
            peak_table.iloc[0]["ell"] if len(peak_table) else np.nan,
            peak_table.iloc[0]["angular_scale_degrees"] if len(peak_table) else np.nan,
            mean_spacing
        ]
    })

    summary.to_csv("outputs/spectrum_summary.csv", index=False)

    save_linear_plot(data, peak_table)
    save_log_plot(data, peak_table)
    save_smoothed_plot(data, smoothed, peak_table)

    print("Loaded data from:", source_path)
    print("\nDetected acoustic peaks:")
    print(peak_table)
    print("\nSaved outputs/peak_table.csv")
    print("Saved outputs/spectrum_summary.csv")


if __name__ == "__main__":
    main()
