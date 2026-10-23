import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio.PDB import PDBParser

model = PDBParser(QUIET=True).get_structure("ndm1", "structures/pdb3spu.ent")[0]
chain = next(iter(model))
coords = np.array([a.coord for a in chain.get_atoms()])
zn = [a for a in model.get_atoms() if a.element == "ZN"]
zn_A = [z for z in zn if np.linalg.norm(coords - z.coord, axis=1).min() < 3.0]

sites = []
for z in zn_A:
    lig = []
    for a in model.get_atoms():
        if a.element in ("N", "O", "S") and np.linalg.norm(a.coord - z.coord) < 2.6:
            lig.append((a, a.get_parent()))
    sites.append((z, lig))

def label(r):
    if r.get_resname() in ("HOH", "WAT"):
        return "H2O/OH"
    return f"{r.get_resname().title()}{r.id[1]}"

for k, (z, lig) in enumerate(sites, 1):
    desc = ", ".join(f"{label(r)}:{a.get_id()} {np.linalg.norm(a.coord - z.coord):.2f}A" for a, r in lig)
    print(f"Zn{k}: {desc}")

pts = [z.coord for z, _ in sites]
for _, lig in sites:
    pts += [a.coord for a, _ in lig]
pts = np.array(pts)
c = pts.mean(0)
P = np.linalg.svd(pts - c)[2][:2]
proj = lambda x: (np.asarray(x) - c) @ P.T

COL = {"N": "#1f77b4", "O": "#d62728", "S": "#e0a800"}
fig, ax = plt.subplots(figsize=(7.5, 6))
drawn = set()
for k, (z, lig) in enumerate(sites, 1):
    zp = proj(z.coord)
    for a, r in lig:
        ap = proj(a.coord)
        ax.plot([zp[0], ap[0]], [zp[1], ap[1]], "--", color="#666", lw=1.3, zorder=1)
        mid = (zp + ap) / 2
        d = np.linalg.norm(a.coord - z.coord)
        ax.text(mid[0], mid[1], f"{d:.1f}", fontsize=7, color="#444", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.85), zorder=2)
        if id(a) not in drawn:
            drawn.add(id(a))
            ax.scatter(*ap, s=300, color=COL[a.element], edgecolor="white", zorder=3)
            v = ap - zp
            v = v / (np.linalg.norm(v) + 1e-9)
            t = ap + 0.55 * v
            ax.text(t[0], t[1], label(r), fontsize=10, ha="center", va="center", fontweight="bold")
for k, (z, _) in enumerate(sites, 1):
    zp = proj(z.coord)
    ax.scatter(*zp, s=900, color="#9a9a9a", edgecolor="black", zorder=4)
    ax.text(zp[0], zp[1], f"Zn{k}", ha="center", va="center", fontsize=9, fontweight="bold", zorder=5)

handles = [Line2D([], [], marker="o", ls="", color=COL[e], markersize=9, label=n)
           for e, n in (("N", "Nitrogen"), ("O", "Oxygen"), ("S", "Sulfur"))]
ax.legend(handles=handles, fontsize=8, frameon=False, loc="lower right")
ax.margins(0.2)
ax.set_aspect("equal")
ax.axis("off")
ax.set_title("NDM-1 di-zinc active site (PDB 3SPU, chain A)", fontsize=12)
fig.text(0.5, 0.01, "2D projection of 3D coordinates. Dashed lines: contacts < 2.6 A, distances in A.",
         ha="center", fontsize=8, color="#555")
fig.savefig("figures/fig_zinc_site.png", dpi=220, bbox_inches="tight")
print("saved figures/fig_zinc_site.png")
