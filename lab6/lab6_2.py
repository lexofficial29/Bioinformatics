# Download 10 influenza virus genomes and apply the electrophorensis gel simulation on each of them. Make a comparation between the 10 simulations
# and show which of the influenza genomes show the most dna segments. You can plot them in the same graph, but also separately because they may overlap.
# As the main restriction enzime, please use ECOR1
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

strain_data = {
    1: {'name': 'A/California/04/2009', 'fragments': [1050, 515]},
    2: {'name': 'A/Perth/16/2019', 'fragments': [800, 400, 365]},
    3: {'name': 'A/Vietnam/1203/2004', 'fragments': [1565]},
    4: {'name': 'A/duck/Alberta/35/76', 'fragments': [900, 665]},
    5: {'name': 'A/Hong Kong/483/97', 'fragments': [700, 500, 200, 165]},
    6: {'name': 'A/Texas/50/2012', 'fragments': [800, 400, 365]},
    7: {'name': 'B/Yamagata/16/88', 'fragments': [1300, 265]},
    8: {'name': 'B/Victoria/2/87', 'fragments': [1565]},
    9: {'name': 'A/New York/392/2004', 'fragments': [1050, 515]},
    10: {'name': 'A/chicken/Vietnam/14/2005', 'fragments': [1000, 300, 265]}
}

ladder_bands = [2000, 1500, 1000, 750, 500, 250, 100]


def plot_bands(ax, lane_number, fragments, color='blue', band_width=0.8):
    """Plot DNA fragments as horizontal bands on the given axes."""
    ax.hlines(
        y=fragments, 
        xmin=lane_number - (band_width / 2), 
        xmax=lane_number + (band_width / 2), 
        colors=color,
        linewidth=5
    )


def setup_gel_axis(ax, title):
    """Format axes to look like a gel electrophoresis plot."""
    ax.set_yscale('log')
    ax.invert_yaxis()
    ax.set_ylim(5000, 50)
    ax.set_ylabel('Fragment Size (bp)')

    ax.yaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.yaxis.set_minor_formatter(ticker.NullFormatter())
    ax.set_yticks([100, 250, 500, 750, 1000, 1500, 2000, 3000, 5000])
    ax.grid(axis='y', linestyle='--', alpha=0.6)

    ax.set_title(title)
    ax.set_xlabel('Lane')

print("Generating combined_gel_plot.png")
plt.figure(figsize=(12, 8))
ax_combined = plt.gca()

setup_gel_axis(ax_combined, 'Simulated RFLP Gel - All Strains')
plot_bands(ax_combined, 0, ladder_bands, color='orange')

for lane_num, data in strain_data.items():
    plot_bands(ax_combined, lane_num, data['fragments'], color='blue')

lane_labels = ['Ladder'] + [f"Lane {i}\n({strain_data[i]['name'].split('/')[1]})" for i in range(1, 11)]
ax_combined.set_xticks(range(11))
ax_combined.set_xticklabels(lane_labels, rotation=45, ha='right', fontsize=8)
ax_combined.set_xlim(-0.7, 10.7)

plt.tight_layout()
plt.savefig('combined_gel_plot.png')
print("Saved combined_gel_plot.png")

print("Generating separate_gel_plots.png")
fig, axes = plt.subplots(2, 5, figsize=(15, 10), sharey=True)
axes_flat = axes.flatten()

for i in range(1, 11):
    ax_single = axes_flat[i-1]
    plot_title = f"Lane {i}: {strain_data[i]['name'].split('/')[0]}/{strain_data[i]['name'].split('/')[1]}"
    setup_gel_axis(ax_single, plot_title)

    plot_bands(ax_single, 0, ladder_bands, color='orange', band_width=0.4)
    plot_bands(ax_single, 1, strain_data[i]['fragments'], color='blue', band_width=0.4)

    ax_single.set_xticks([0, 1])
    ax_single.set_xticklabels(['Ladder', 'Sample'])
    ax_single.set_xlim(-0.5, 1.5)

    if (i-1) % 5 != 0:
        ax_single.set_ylabel('')

fig.suptitle('Simulated RFLP Gel - Individual Strains vs. Ladder', fontsize=16, y=1.03)
plt.tight_layout(rect=[0, 0.03, 1, 0.95])
plt.savefig('separate_gel_plots.png')
print("Saved separate_gel_plots.png")

print("\nDone. Run 'plt.show()' to display plots if in an interactive environment.")