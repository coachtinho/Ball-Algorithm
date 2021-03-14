import pandas as pd
import matplotlib.pyplot as plt
import sys

argv = sys.argv

if len(argv) != 3:
    print("Usage: python visualizer.py input_file output_file")

data = pd.read_csv(argv[1])

x = data["x"].values
y = data["y"].values

plt.figure()
plt.title("Points")
plt.xlabel("x")
plt.ylabel("y")
plt.scatter(x, y)
plt.savefig(f"{argv[2]}")
