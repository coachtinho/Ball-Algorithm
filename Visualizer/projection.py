import pandas as pd
import matplotlib.pyplot as plt
import sys

argv = sys.argv

if len(argv) != 3:
    print("Usage: python projection.py input_file output_file")

data = pd.read_csv(argv[1])

x = data["x"].values
y = data["y"].values
proj_x = data["projx"].values
proj_y = data["projy"].values

plt.figure()
plt.title("Projections")
plt.xlabel("x")
plt.ylabel("y")
for i in range(1, len(x)):
    plt.scatter(x[i], y[i], c="blue")
    plt.scatter(proj_x[i], proj_y[i], c="red")
    plt.plot([x[i], proj_x[i]], [y[i], proj_y[i]], color="red")
plt.scatter(x[0], y[0], c="blue")
plt.scatter(proj_x[0], proj_y[0], c="blue")
plt.plot([x[0], proj_x[0]], [y[0], proj_y[0]], color="red")

plt.savefig(f"{argv[2]}")
