import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio.PDB import PDBParser, PPBuilder

parser = PDBParser(QUIET=True)

STRUCTS = {
    "NDM-1": ("3SPU", "Class B (metallo-beta-lactamase)"),
    "KPC-2": ("2OV5", "Class A (serine)"),
    "OXA-48": ("3HBR", "Class D (serine)"),
    "CTX-M-15": ("4HBT", "Class A (serine)"),
}

def ca_coords(pid):
    model = parser.get_structure(pid, f"structures/pdb{pid.lower()}.ent")[0]
    chain = next(iter(model))
    out = []
    for pp in PPBuilder().build_peptides(chain):
        for r in pp:
            if "CA" in r:
                out.append(r["CA"].coord)
    return np.array(out)

fig, axes = plt.subplots(2, 2, figsize=(10, 9), constrained_layout=True)
for ax, (name, (pid, cls)) in zip(axes.ravel(), STRUCTS.items()):
    xyz = ca_coords(pid)
    d = np.linalg.norm(xyz[:, None] - xyz[None], axis=-1)
    im = ax.imshow(d, cmap="magma_r", vmin=0, vmax=30)
    ax.set_title(f"{name}\n{cls}", fontsize=11)
for ax in axes[1]:
    ax.set_xlabel("Residue")
for ax in axes[:, 0]:
    ax.set_ylabel("Residue")
fig.colorbar(im, ax=axes.ravel().tolist(), shrink=0.8, label="CA-CA distance (Å)")
fig.suptitle("Residue distance maps of South Asian AMR enzymes", fontsize=14)
fig.savefig("figures/fig_distance_maps.png", dpi=220, bbox_inches="tight")

model = parser.get_structure("ndm1", "structures/pdb3spu.ent")[0]
chain = next(iter(model))
coords = np.array([a.coord for a in chain.get_atoms()])
zn = [a for a in model.get_atoms() if a.element == "ZN"]
zn_A = [z for z in zn if np.linalg.norm(coords - z.coord, axis=1).min() < 3.0]

for i, z in enumerate(zn_A, 1):
    lig = sorted({(r.id[1], r.get_resname()) for r in chain if r.id[0] == " "
                  for a in r if np.linalg.norm(a.coord - z.coord) < 2.6})
    print(f"Zn site {i}:", ", ".join(f"{n}{num}" for num, n in lig))

rows = []
for r in chain:
    if r.id[0] != " " or "CA" not in r:
        continue
    d = min(np.linalg.norm(a.coord - z.coord) for a in r for z in zn_A)
    rows.append((r.id[1], r.get_resname(), d))

for num in (88, 154):
    hit = [x for x in rows if x[0] == num][0]
    print(f"Residue {num}: {hit[1]}, {hit[2]:.1f} A from nearest Zn")

LIG = {120, 122, 124, 189, 208, 250}
GREY, BLUE, RED = "#c9c9c9", "#1f77b4", "#d62728"
other = [r for r in rows if r[0] not in LIG and r[0] not in (88, 154)]
ligs = [r for r in rows if r[0] in LIG]

fig, ax = plt.subplots(figsize=(10, 4.5))
ax.axhspan(0, 3, color=BLUE, alpha=0.10)
ax.scatter([r[0] for r in other], [r[2] for r in other], s=10, color=GREY, label="Other residues")
ax.scatter([r[0] for r in ligs], [r[2] for r in ligs], s=45, color=BLUE,
           label="Zn ligands (His120/122/189, Asp124, Cys208, His250)")
for num, label in ((88, "V88L"), (154, "M154L")):
    d = [r[2] for r in rows if r[0] == num][0]
    ax.scatter(num, d, s=130, color=RED, zorder=3)
    ax.annotate(f"{label}\n{d:.1f} Å", (num, d), xytext=(num + 12, d + 5),
                fontsize=11, color=RED, fontweight="bold", zorder=4,
                bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="none", alpha=0.9),
                arrowprops=dict(arrowstyle="-", color=RED))
ax.text(46, 0.7, "first coordination shell (< 3 Å)", fontsize=8, color=BLUE, ha="left")
ax.set_xlabel("Residue number (NDM-1)")
ax.set_ylabel("Distance to nearest Zn (Å)")
ax.set_title("NDM-5 substitutions (V88L, M154L) mapped onto NDM-1 (PDB 3SPU)", fontsize=12)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.legend(fontsize=8, frameon=False, loc="upper right")
fig.savefig("figures/fig_mutation_map.png", dpi=220, bbox_inches="tight")
print("saved figures/fig_distance_maps.png and figures/fig_mutation_map.png")
