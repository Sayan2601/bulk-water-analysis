import pandas as pd
import matplotlib.pyplot as plt

csv_file = "tcf-15ps-p1-p2-p3.csv"

df = pd.read_csv(csv_file)

plt.figure(figsize=(8, 5))

plt.plot(df["time"], df["p1"], label=r"$C_1(t)$")
plt.plot(df["time"], df["p2"], label=r"$C_2(t)$")
plt.plot(df["time"], df["p3"], label=r"$C_3(t)$")

plt.xlabel("Time (ps)")
plt.ylabel("Orientational Correlation")
plt.title("Orientational Correlation Function")

plt.legend()
plt.tight_layout()
plt.savefig("p1_p2_p3_correlation.png", dpi=300)
plt.show()
