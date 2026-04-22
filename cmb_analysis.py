import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# load data
df = pd.read_csv("cmb_data.csv")

ell = df["ell"].to_numpy()
Dl = df["Dl"].to_numpy()

# find first major peak
peak_index = np.argmax(Dl)
peak_ell = ell[peak_index]
peak_value = Dl[peak_index]

print(f"First major visible peak near l = {peak_ell}")
print(f"Peak value = {peak_value:.2f}")

# plot spectrum
plt.figure()
plt.plot(ell, Dl, marker="o", label="CMB Power Spectrum")
plt.scatter([peak_ell], [peak_value], color="red", label=f"Peak near l = {peak_ell}")
plt.xlabel("Multipole moment l")
plt.ylabel("D_l (μK^2)")
plt.title("CMB Temperature Angular Power Spectrum")
plt.legend()
plt.savefig("power_spectrum.png", dpi=300)
plt.show()