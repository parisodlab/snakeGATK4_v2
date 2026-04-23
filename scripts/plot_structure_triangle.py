import sys
import pandas as pd
import matplotlib.pyplot as plt
import math

qfile, metafile, outfile = sys.argv[1:4]

# STRUCTURE q-file: usually one row per individual, columns are Q1..QK
q = pd.read_csv(qfile, delim_whitespace=True, header=None)
if q.shape[1] != 3:
    raise ValueError("Triangle plot requires K=3")

q.columns = ["Q1", "Q2", "Q3"]

meta = pd.read_csv(metafile, sep=None, engine="python")
# assume same order as structure input; if not, merge by sample name after adding names explicitly

# barycentric to Cartesian
x = q["Q2"] + 0.5 * q["Q3"]
y = (math.sqrt(3) / 2.0) * q["Q3"]

plt.figure(figsize=(6, 6))
plt.scatter(x, y, s=25)

# triangle edges
plt.plot([0, 1], [0, 0], "k-")
plt.plot([0, 0.5], [0, math.sqrt(3)/2], "k-")
plt.plot([1, 0.5], [0, math.sqrt(3)/2], "k-")

plt.text(-0.03, -0.03, "Cluster 1")
plt.text(1.01, -0.03, "Cluster 2")
plt.text(0.48, math.sqrt(3)/2 + 0.03, "Cluster 3")

plt.axis("off")
plt.tight_layout()
plt.savefig(outfile)