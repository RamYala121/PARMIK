import pandas as pd
import pysam


def extract_best_match_lengths(sam_file, output_tsv):
    best_matches = {}

    # Open and parse SAM file
    samfile = pysam.AlignmentFile(sam_file, "r")

    for read in samfile.fetch(until_eof=True):
        if read.is_unmapped:
            continue

        query_id = read.query_name
        alignment_len = read.query_alignment_length

        # Get alignment score (AS tag) or fall back to mapping quality
        score = read.get_tag("AS") if read.has_tag("AS") else read.mapping_quality

        # Keep the match with the highest score per query read
        if query_id not in best_matches or score > best_matches[query_id]["score"]:
            best_matches[query_id] = {
                "score": score,
                "alignment_length": alignment_len,
            }

    samfile.close()

    # Convert to DataFrame and write to TSV
    df = pd.DataFrame(
        [
            (q_id, data["alignment_length"])
            for q_id, data in best_matches.items()
        ],
        columns=["Query_ID", "Best_Match_Length"],
    )

    df.to_csv(output_tsv, sep="\t", index=False)
    print(f"Saved raw data to {output_tsv}")


if __name__ == "__main__":
    extract_best_match_lengths(
        "minimap2_output.sam", "minimap2_best_match_lengths.tsv"
    )