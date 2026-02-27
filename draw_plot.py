import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots(figsize=(12, 6))

colors_upper = ["#f3b5bf", "#d67f8f", "#ab495c", "#d67f8f", "#f3b5bf"]  
colors_lower = ["#f9d2e2", "#f4b3ce", "#f194b8", "#f4b3ce", "#f9d2e2"] 

bins = [70, 90, 110, 140, 170, 200]
labels_lower = [
    "CR2",
    "VR2",
    "CR0",
    "VR2",
    "CR2"
]
labels_upper = [
    "CR1",
    "VR1",
    "SR",
    "VR1",
    "CR1"
]

y_cut = 0.95  # NN score cut

# --- draw shaded rectangles below and above the cut ---
for i in range(len(bins)-1):
    width = bins[i+1] - bins[i]

    # lower region (0 -> y_cut)
    ax.add_patch(patches.Rectangle(
        (bins[i], 0), width, y_cut,
        facecolor=colors_lower[i],
        edgecolor='none', alpha=0.8
    ))
    x_center = bins[i] + width / 2
    y_center = y_cut / 2
    ax.text(x_center, y_center, labels_lower[i],
            ha='center', va='center',
            fontsize=15, fontweight='bold', color='black')

    # upper region (y_cut -> 1)
    ax.add_patch(patches.Rectangle(
        (bins[i], y_cut), width, 1 - y_cut,
        facecolor=colors_upper[i],
        edgecolor='none', alpha=0.8
    ))
    y_center_upper = y_cut + (1 - y_cut) / 2
    ax.text(x_center, y_center_upper, labels_upper[i],
            ha='center', va='center',
            fontsize=15, fontweight='bold', color='black')

# --- cut line ---
ax.axhline(y=y_cut, color="black", linestyle="--", linewidth=1.5)
ax.text(215, y_cut + 0.02, f"NN cut = {y_cut}", ha="right", va="bottom", fontsize=12)

# axis settings
ax.set_xticks(bins)
ax.set_xticklabels([str(b) for b in bins])

ax.set_xlim(50, 220)
ax.set_ylim(0, 1)

ax.set_xlabel("m(H) [GeV]")
ax.set_ylabel("NN Score")

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()
