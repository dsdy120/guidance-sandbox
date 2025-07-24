import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("eci_orbit.csv")
plt.plot(df["x"], df["y"])
plt.axis("equal")
plt.xlabel("X (km)")
plt.ylabel("Y (km)")
plt.title("ECI Orbit (XY projection)")
plt.grid()
plt.show()