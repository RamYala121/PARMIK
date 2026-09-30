import os
import matplotlib.pyplot as plt
import pandas as pd
import pysam
import seaborn as sns


def parse_sam_best_matches(sam_file):
    """Parses a SAM file and extracts the alignment length for the best match per query."""
    best_matches = {}

    print(f"Parsing {sam_file}...")
    samfile = pysam.AlignmentFile(sam_file, "r")

    for read in samfile.fetch(until_eof=True):
        # Skip unmapped reads
        if read.is_unmapped:
            continue

        query_id = read.query_name
        alignment_len = read.query_alignment_length  # Length of the aligned part of the query

        # Extract alignment score (AS tag) if present, otherwise use mapq or alignment length
        score = read.get_tag("AS") if read.has_tag("AS") else read.mapping_quality

        # Keep the record with the highest alignment score/mapping quality for each query
        if query_id not in best_matches or score > best_matches[query_id]["score"]:
            best_matches[query_id] = {
                "score": score,
                "alignment_length": alignment_len,
            }

    samfile.close()

    df = pd.DataFrame.from_dict(best_matches, orient="index")
    print(f"Extracted best matches for {len(df)} mapped queries.")
    return df


def plot_histogram(df, output_img="minimap2_best_match_lengths.png"):
    """Generates and saves a histogram of best match alignment lengths."""
    if df.empty:
        print("No mapped reads found in the SAM file.")
        return

    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")

    # Draw histogram with KDE curve
    sns.histplot(
        df["alignment_length"],
        kde=True,
        bins=30,
        color="skyblue",
        edgecolor="black",
    )

    plt.title(
        "Minimap2: Distribution of Best Match Alignment Lengths per Query",
        fontsize=14,
        fontweight="bold",
    )
    plt.xlabel("Alignment Length (bp)", fontsize=12)
    plt.ylabel("Query Count", fontsize=12)

    # Annotate summary statistics on the plot
    median_len = df["alignment_length"].median()
    mean_len = df["alignment_length"].mean()
    plt.axvline(
        median_len,
        color="red",
        linestyle="--",
        linewidth=1.5,
        label=f"Median: {median_len:.1f} bp",
    )
    plt.axvline(
        mean_len,
        color="green",
        linestyle=":",
        linewidth=1.5,
        label=f"Mean: {mean_len:.1f} bp",
    )

    plt.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(output_img, dpi=300)
    print(f"Histogram successfully saved to {output_img}")


if __name__ == "__main__":
    sam_filename = "minimap2_output.sam"

    if not os.path.exists(sam_filename):
        print(f"Error: {sam_filename} does not exist in the current directory.")
    else:
        df_results = parse_sam_best_matches(sam_filename)
        plot_histogram(df_results)