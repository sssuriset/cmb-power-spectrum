import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

# load real Planck data
data = np.loadtxt("COM_PowerSpect_CMB-TT-full_R3.01.txt")

l = data[:, 0]
Dl = data[:, 1]
err_low = data[:, 2]
err_high = data[:, 3]

# find significant peaks
peaks, properties = find_peaks(Dl, distance=150, prominence=500)

# keep only peaks in the acoustic-peak region
major_peaks = peaks[l[peaks] > 150]

# take the first three in left-to-right order
top_peaks = major_peaks[:3]

for i, p in enumerate(top_peaks, start=1):
    print(f"Peak {i} at l ≈ {int(l[p])}")

# main plot
plt.figure(figsize=(8,5))

plt.errorbar(
    l, Dl,
    yerr=[err_low, err_high],
    fmt='o',
    markersize=3,
    alpha=0.4,
    label="Planck TT Data"
)

p0 = top_peaks[0]
plt.scatter(l[p0], Dl[p0], color="red", s=60, label=f"Peak 1 (l ≈ {int(l[p0])})")
plt.text(l[p0] + 20, Dl[p0], f"{int(l[p0])}", fontsize=8)

for i, p in enumerate(top_peaks[1:], start=2):
    plt.scatter(l[p], Dl[p], color="orange", s=40, label=f"Peak {i} (l ≈ {int(l[p])})")
    plt.text(l[p] + 20, Dl[p], f"{int(l[p])}", fontsize=8)

plt.xlabel("Multipole moment l")
plt.ylabel("D_l (μK^2)")
plt.title("CMB Temperature Angular Power Spectrum (Planck Data)")
plt.xlim(0, 2000)
plt.legend()

plt.savefig("power_spectrum.png", dpi=300)
plt.show()

# log-scale plot
plt.figure()

plt.errorbar(
    l, Dl,
    yerr=[err_low, err_high],
    fmt='o',
    markersize=3,
    alpha=0.4
)

plt.xscale("log")
plt.yscale("log")
plt.xlim(10, 2000)

plt.xlabel("Multipole moment l")
plt.ylabel("D_l (μK^2)")
plt.title("CMB Power Spectrum (Log Scale)")

plt.savefig("power_spectrum_log.png", dpi=300)
plt.show()