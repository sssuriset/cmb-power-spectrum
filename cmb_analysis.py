import numpy as np
import matplotlib.pyplot as plt

# load real Planck data
data = np.loadtxt("COM_PowerSpect_CMB-TT-full_R3.01.txt")

l = data[:, 0]
Dl = data[:, 1]
err_low = data[:, 2]
err_high = data[:, 3]

# find peak
peak_index = np.argmax(Dl)

# plot
plt.figure(figsize=(8,5))

plt.errorbar(
    l, Dl,
    yerr=[-err_low, err_high],
    fmt='o',
    markersize=3,
    label="Planck TT Data"
)

plt.scatter(
    l[peak_index],
    Dl[peak_index],
    color="red",
    label=f"Peak l ≈ {int(l[peak_index])}"
)

plt.xlabel("Multipole moment l")
plt.ylabel("D_l (μK^2)")
plt.title("CMB Temperature Angular Power Spectrum (Planck Data)")
plt.xlim(0, 2000)
plt.legend()

plt.savefig("power_spectrum.png", dpi=300)
plt.show()