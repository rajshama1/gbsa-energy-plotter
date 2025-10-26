import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import re
from matplotlib import rcParams

# Use Unicode-friendly font
rcParams["font.family"] = ["Liberation Sans", "DejaVu Sans", "sans-serif"]

# Load your data
df = pd.read_csv("example_gbsa.csv")

# Extract means and SDs from "value (SD)" strings
def extract_mean_sd(series):
    # Match: number (number)
    extracted = series.str.extract(r'([-+]?\d*\.\d+|\d+)\s*\(([-+]?\d*\.\d+|\d+)\)')
    mean = extracted[0].astype(float)
    sd = extracted[1].astype(float)
    return mean, sd

# Columns to process
columns_to_extract = [
    "ENTROPY (-TS)", "VDWAALS", "EEL", "EGB", "ESURF",
    "GGAS", "GSOLV", "TOTAL", "Gbinding"
]
new_col_names = [
    "ENTROPY (-T\u0394S)", "\u0394VDWAALS", "\u0394EEL", "\u0394EGB", "\u0394ESURF",
    "\u0394GGAS", "\u0394GSOLV", "\u0394TOTAL", "\u0394G_binding"
]

# Create cleaned DataFrame
df_clean = pd.DataFrame()
df_clean["Complex"] = df["Complex"]
for old_col, new_col in zip(columns_to_extract, new_col_names):
    mean, sd = extract_mean_sd(df[old_col])
    df_clean[new_col] = mean
    df_clean[new_col + "_SD"] = sd  # keep SD

# Melt for seaborn
df_melted = df_clean.melt(
    id_vars="Complex", 
    value_vars=new_col_names, 
    var_name="Energy Term", 
    value_name="Energy (kcal/mol)"
)

# Add SD column for error bars
sd_data = df_clean.melt(
    id_vars="Complex",
    value_vars=[c + "_SD" for c in new_col_names],
    var_name="Energy Term_SD",
    value_name="SD"
)
# Match Energy Terms correctly
sd_data["Energy Term"] = sd_data["Energy Term_SD"].str.replace("_SD", "", regex=False)
df_melted["SD"] = sd_data["SD"]

# Custom palette
custom_palette = {
    'ENTROPY (-T\u0394S)': '#1f77b4',
    '\u0394VDWAALS': '#ff7f0e',
    '\u0394EEL': '#2ca02c',
    '\u0394EGB': '#d62728',
    '\u0394ESURF': '#9467bd',
    '\u0394GGAS': '#8c564b',
    '\u0394GSOLV': '#e377c2',
    '\u0394TOTAL': 'black',
    '\u0394G_binding': '#7f7f7f'
}

# Plot
plt.figure(figsize=(16, 8), dpi=600)
sns.set(style="whitegrid")
ax = plt.gca()

# Barplot with error bars
bars = sns.barplot(
    data=df_melted,
    x="Complex",
    y="Energy (kcal/mol)",
    hue="Energy Term",
    palette=custom_palette,
    edgecolor="black",
    linewidth=0.8,
    ci=None
)

for bar in ax.patches:
    bar.set_width(bar.get_width() * 0.9)  # shrink width to 70%

# Add error bars + annotations
for i, bar in enumerate(ax.patches):
    height = bar.get_height()
    x = bar.get_x() + bar.get_width() / 2
    term = df_melted.iloc[i]["Energy Term"]
    sd = df_melted.iloc[i]["SD"]

    if term != "\u0394G_binding":  # skip SD for ΔG_binding
        if np.isfinite(height) and np.isfinite(sd):
            ax.errorbar(x, height, yerr=sd, color='black', capsize=3, fmt='none', lw=1)
# Always annotate (vertical placement)
    if np.isfinite(height):
        ax.text(
            x,
            height + (2 if height >= 0 else -2),  # adjust offset
            f"{height:.2f}",
            ha='center',
            va='bottom' if height >= 0 else 'top',
            fontsize=14,
            fontweight="bold",
            rotation=90,          # 🔹 Vertical text
            color="black"
        )
"""
    # Always annotate
    if np.isfinite(height):
        if term == "\u0394G_binding":
            ax.text(
                x, height + (1 if height >= 0 else -1),
                f"{height:.2f}",
                ha='center', va='bottom' if height >= 0 else 'top',
                fontsize=14, fontweight="bold", color="black"  # bigger, bold
            )
        else:
            ax.text(
                x, height + (1 if height >= 0 else -1),
                f"{height:.2f}",
                ha='center', va='bottom' if height >= 0 else 'top',
                fontsize=14, fontweight="bold"
            )
"""
# Background
compounds = df_clean["Complex"].tolist()
background_colors = ['lightyellow', 'lavender']
for i in range(len(compounds)):
    ax.axvspan(i - 0.5, i + 0.5, facecolor=background_colors[i % 2], alpha=0.2, zorder=0)
    
for spine in ax.spines.values():
    spine.set_edgecolor("black")
    spine.set_linewidth(1.2)   # control thickness
    
ax.tick_params(color="black", labelcolor="black", width=1.2)

# Labels and formatting
plt.ylabel("Energy (kcal/mol)", fontsize=22, fontweight="bold")

ax.set_xlabel("")

plt.xticks(
    ticks=range(len(df_clean["Complex"])),
    labels=df_clean["Complex"],
    rotation=0,
    fontsize=22,
    ha='center',
    fontweight="bold"
)

import matplotlib.ticker as ticker

# === Y-axis settings ===
plt.yticks(np.arange(-100, 60, 20), fontsize=22)  # major ticks every 10 units
ax.set_ylim(-100, 60)
ax.set_xlim(-0.5, len(compounds) - 0.5)
ax.axhline(0, color='black', linewidth=1.2)

# === Gridline control ===
# Add minor ticks every 5 units
ax.yaxis.set_minor_locator(ticker.MultipleLocator(5))

# Draw only Y gridlines
ax.grid(axis='y', which='major', alpha=0.2, linewidth=1.2)   # major (every 10)
ax.grid(axis='y', which='minor', alpha=0.2, linestyle='--', linewidth=0.8)  # minor (every 5)


leg = plt.legend(
    bbox_to_anchor=(1.02, 1.0),
    loc="upper left",
    fontsize=16,
    title="Energy Term",
    title_fontsize=18,
    frameon=True
)

# Now you can style the frame
leg.get_frame().set_edgecolor("black")
leg.get_frame().set_linewidth(1.2)

ax.spines["bottom"].set_linewidth(1.2)
ax.spines["left"].set_linewidth(1.2)
ax.tick_params(width=1.2)
#ax.grid(axis='y', alpha=0.2)
#ax.grid(alpha=0.2)

plt.tight_layout()
plt.savefig("MMGBSA_plot_with_SD.tif", format="tiff", dpi=600)
plt.show()
