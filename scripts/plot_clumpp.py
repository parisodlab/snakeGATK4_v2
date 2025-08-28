import pandas as pd
import matplotlib.pyplot as plt
import re
import sys

clumpp_file = sys.argv[1]  # e.g., "results/clumpp/output_K2.txt"
structure_file = sys.argv[2]  # e.g., "data/structure_input.str"
outfile = sys.argv[3]  # e.g., "results/plots/clumpp_K2.pdf"
K = int(sys.argv[4])  # number of clusters

# --- Step 1: Extract unique individual names (grouping by consecutive identical names) ---
with open(structure_file) as f:
    lines = f.readlines()

names = []
last_name = None
for line in lines:
    if line.strip() == "":
        continue
    name = line.strip().split()[0]
    if name != last_name:
        names.append(name)
        last_name = name

print(f"Detected {len(names)} unique individuals from STRUCTURE input")

# --- Step 2: Parse CLUMPP output ---
records = []
with open(clumpp_file) as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        m = re.match(r"^\s*(\d+).*?\((\d+)\).*?:\s*(.+)$", line)
        if not m:
            continue
        idx = int(m.group(1))
        missing = int(m.group(2))
        props = m.group(3).split()
        if len(props) != K:
            raise ValueError(f"Line {idx}: expected {K} proportions, got {len(props)}")
        props = list(map(float, props))
        try:
            name = names[idx - 1]  # CLUMPP uses 1-based index
        except IndexError:
            raise ValueError(
                f"Index {idx} out of range. STRUCTURE input only has {len(names)} unique names."
            )
        records.append(
            {
                "Index": idx,
                "Name": name,
                "Missing": missing,
                **{f"Cluster_{i+1}": props[i] for i in range(K)},
            }
        )

df = pd.DataFrame(records).sort_values("Index").reset_index(drop=True)

# --- Step 3: Plot ---

fig, ax = plt.subplots(figsize=(max(10, len(df) / 4), 4))
bottom = [0] * len(df)
for i in range(K):
    ax.bar(df["Name"], df[f"Cluster_{i+1}"], bottom=bottom, label=f"Cluster {i+1}")
    bottom = [bottom[j] + df[f"Cluster_{i+1}"].iloc[j] for j in range(len(df))]

# Missingness percentage as text above bars
for i, miss in enumerate(df["Missing"]):
    ax.text(i, 1.02, f"{miss:.0f}%", ha="center", va="bottom", fontsize=6, rotation=90)

ax.set_xticks(range(len(df)))
ax.set_xticklabels(df["Name"], rotation=90, fontsize=6)
ax.set_ylabel("Ancestry proportion")
ax.set_xlabel("Individual")
ax.set_title(f"CLUMPP K={K}")

# Background shading for missingness
norm = plt.Normalize(df["Missing"].min(), df["Missing"].max())
cmap = plt.cm.Reds
for i, miss in enumerate(df["Missing"]):
    ax.axvspan(i - 0.5, i + 0.5, color=cmap(norm(miss)), alpha=0.2)

cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)
cbar.set_label("Missingness (%)")

ax.set_ylim(0, 1.1)
ax.legend(ncol=2, fontsize="small")
plt.tight_layout()
plt.savefig(outfile)
