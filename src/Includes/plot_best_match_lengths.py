import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


def generate_histogram(
    tsv_path="minimap2_best_match_lengths.tsv",
    output_png="minimap2_best_match_lengths_updated.png",
):
    # Load the verified dataset
    df = pd.read_csv(tsv_path, sep="\t")

    # Calculate summary metrics
    unmapped_count = (df["Best_Match_Length"] == 0).sum()
    median_val = df["Best_Match_Length"].median()
    mean_val = df["Best_Match_Length"].mean()

    # Configure plot styling
    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")

    # Plot histogram
    sns.histplot(
        df["Best_Match_Length"],
        bins=30,
        color="skyblue",
        edgecolor="black",
        kde=False,
    )

    # Add title and axis labels
    plt.title(
        "Minimap2: Best Match Alignment Lengths (Including Unmapped Reads)",
        fontsize=14,
        fontweight="bold",
    )
    plt.xlabel("Alignment Length (bp)", fontsize=12)
    plt.ylabel("Query Count", fontsize=12)

    # Add vertical lines for Median and Mean
    plt.axvline(
        median_val,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Median: {median_val:.1f} bp",
    )
    plt.axvline(
        mean_val,
        color="green",
        linestyle=":",
        linewidth=1.5,
        label=f"Mean: {mean_val:.1f} bp",
    )

    # Annotate the unmapped reads peak at 0 bp
    plt.annotate(
        f"Unmapped (0 bp): {unmapped_count:,}",
        xy=(0, unmapped_count),
        xytext=(20, unmapped_count * 0.85),
        arrowprops=dict(facecolor="black", shrink=0.05, width=1, headwidth=6),
        fontsize=11,
        fontweight="bold",
        color="darkred",
    )

    plt.legend(fontsize=11)
    plt.tight_layout()

    # Save high-resolution PNG
    plt.savefig(output_png, dpi=300)
    plt.close()
    print(f"Successfully generated and saved plot to {output_png}")


if __name__ == "__main__":
    generate_histogram()