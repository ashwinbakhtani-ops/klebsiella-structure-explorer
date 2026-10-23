import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.DSSP import residue_max_acc

model = PDBParser(QUIET=True).get_structure("ndm1", "structures/pdb3spu.ent")[0]
chain = next(iter(model))
for r in [r for r in chain if r.id[0] != " "]:
    chain.detach_child(r.id)

ShrakeRupley().compute(chain, level="R")
MAXACC = residue_max_acc["Sander"]

rows = []
for r in chain:
    if "CA" not in r:
        continue
    name = r.get_resname()
    if name not in MAXACC:
        continue
    rel = min(r.sasa / MAXACC[name], 1.0)
    b = float(np.mean([a.get_bfactor() for a in r]))
    rows.append((r.id[1], r.get_resname(), rel, b))

num = np.array([x[0] for x in rows])
rel = np.array([x[2] for x in rows])
bf = np.array([x[3] for x in rows])
z = (bf - bf.mean()) / bf.std()

for n in (88, 154):
    i = int(np.where(num == n)[0][0])
    print(f"Residue {n} {rows[i][1]}: relative SASA {rel[i]:.2f}, "
          f"B-factor {bf[i]:.1f} (z = {z[i]:+.2f}, "
          f"more flexible than {100 * (bf < bf[i]).mean():.0f}% of residues)")
print(f"Median relative SASA: {np.median(rel):.2f}, mean B-factor: {bf.mean():.1f}")

LIG = [120, 122, 124, 189, 208, 250]
GREY, BLUE, RED = "#c9c9c9", "#1f77b4", "#d62728"

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6.5), sharex=True,
                               constrained_layout=True)
for ax, y, lab in ((ax1, rel, "Relative solvent accessibility"),
                   (ax2, z, "B-factor (z-score)")):
    ax.plot(num, y, color=GREY, lw=1.2, zorder=1)
    for l in LIG:
        ax.axvline(l, color=BLUE, alpha=0.25, lw=1)
    for n, name in ((88, "V88L"), (154, "M154L")):
        i = int(np.where(num == n)[0][0])
        ax.scatter(n, y[i], s=90, color=RED, zorder=3)
        ax.annotate(name, (n, y[i]), xytext=(n + 6, y[i] + (0.08 if ax is ax1 else 0.5)),
                    color=RED, fontsize=11, fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.9))
    ax.set_ylabel(lab)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
ax1.axhline(0.2, color="black", ls=":", lw=0.8)
ax1.text(num.max(), 0.21, "buried < 0.2", ha="right", fontsize=8)
ax2.axhline(0, color="black", ls=":", lw=0.8)
ax2.set_xlabel("Residue number (NDM-1)")
fig.suptitle("Burial and flexibility of NDM-5 substitution sites in NDM-1 (PDB 3SPU)",
             fontsize=12)
fig.text(0.01, -0.01, "Blue lines: Zn-coordinating residues", fontsize=8, color=BLUE)
fig.savefig("figures/fig_burial_flexibility.png", dpi=220, bbox_inches="tight")
print("saved figures/fig_burial_flexibility.png")
