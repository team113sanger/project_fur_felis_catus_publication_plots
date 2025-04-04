#!/usr/bin/env python
import argparse
import pandas as pd
import re
import os

from fur2cosmic.utils.logging_utils import setup_logging, get_package_logger

LOGGER = get_package_logger()


def create_synthetic_data(
    input_path,
    output_path,
    genes=["TP53"],
    variant_types=["missense_variant", "synonymous_variant"],
):
    """
    Create synthetic COSMIC CMC data by:
    1. Filtering by specified genes and variant types
    2. Replacing COSMIC IDs with synthetic ones

    Args:
        input_path (str): Path to the original COSMIC CMC file
        output_path (str): Path to save the synthetic file
        genes (list): List of HUGO gene symbols to filter by (GENE_SYMBOL column)
        variant_types (list): List of variant types to filter by (MUTATION_DESCRIPTION column)
    """
    # Read the input file with low_memory=False to handle mixed data types
    LOGGER.info(f"Reading input file: {input_path}")
    df = pd.read_csv(input_path, sep="\t", low_memory=False)
    LOGGER.info(f"Original dataset has {len(df)} rows")

    # Filter by genes and variant types
    gene_mask = df["GENE_SYMBOL"].isin(genes)
    variant_mask = df["MUTATION_DESCRIPTION"].isin(variant_types)
    df = df[gene_mask & variant_mask]
    LOGGER.info(f"After filtering, dataset has {len(df)} rows")

    # Reset the index to use it as row-index (starting from 0)
    df = df.reset_index(drop=True)

    # List of columns to modify
    cosmic_columns = [
        "COSMIC_SAMPLE_ID",
        "COSMIC_GENE_ID",
        "COSMIC_PHENOTYPE_ID",
        "COSMIC_STUDY_ID",
        "GENOMIC_MUTATION_ID",
        "LEGACY_MUTATION_ID",
        "MUTATION_ID",
    ]

    # Anonymize SAMPLE_NAME
    if "SAMPLE_NAME" in df.columns:
        LOGGER.info("Anonymizing SAMPLE_NAME column")
        # Convert index to string before concatenation to avoid dtype warnings
        df["SAMPLE_NAME"] = "ANON_" + df.index.astype(str)

    # Create synthetic IDs for COSMIC columns
    for column in cosmic_columns:
        if column not in df.columns:
            LOGGER.warning(f"Column '{column}' not found in the input file")
            continue

        LOGGER.info(f"Processing column: {column}")

        # Convert column to string to ensure consistent handling
        # This prevents dtype warnings when working with mixed types
        if not pd.api.types.is_string_dtype(df[column]):
            df[column] = df[column].astype(str)

        # Skip NA values
        mask = ~df[column].isna()

        # Create a function to handle each pattern
        def transform_cosmic_id(x):
            match = re.match(r"COS([A-Za-z])", str(x))
            if match:
                letter = match.group(1)
                return f"SYNTH_{letter}_"
            return ""

        # For values that aren't NA, transform them using string operations
        if mask.any():
            # First convert matching patterns to SYNTH_{letter}_
            df.loc[mask, column] = df.loc[mask, column].apply(transform_cosmic_id)

            # Then concatenate with the index (ensure both sides are strings)
            df.loc[mask, column] = df.loc[mask, column] + df.loc[mask].index.astype(str)

            # For values that don't match the COS pattern (resulting in empty strings),
            # just use the index - explicitly convert to string
            empty_mask = (df[column] == "") & mask
            df.loc[empty_mask, column] = df.loc[empty_mask].index.astype(str)

    # Save the synthetic data
    LOGGER.info(f"Saving synthetic data to {output_path}")
    df.to_csv(output_path, sep="\t", index=False)
    LOGGER.info(f"Generated {len(df)} rows of synthetic data")
    return


def main():
    parser = argparse.ArgumentParser(description="Create synthetic COSMIC CMC data")
    parser.add_argument("input_path", help="Path to the original COSMIC CMC file")
    parser.add_argument("output_path", help="Path to save the synthetic file")
    parser.add_argument(
        "--genes",
        nargs="+",
        default=["TP53"],
        help="HUGO gene symbols to filter by (default: TP53)",
    )
    parser.add_argument(
        "--variant_types",
        nargs="+",
        default=["missense_variant", "synonymous_variant"],
        help="Variant types to filter by (default: missense_variant, synonymous_variant)",
    )

    args = parser.parse_args()

    # Check if the input file exists
    if not os.path.exists(args.input_path):
        LOGGER.error(f"Input file '{args.input_path}' not found")
        return

    # Create the output directory if it doesn't exist
    output_dir = os.path.dirname(args.output_path)
    if output_dir and not os.path.exists(output_dir):
        LOGGER.info(f"Creating output directory: {output_dir}")
        os.makedirs(output_dir)

    create_synthetic_data(
        args.input_path, args.output_path, args.genes, args.variant_types
    )
    LOGGER.info("Synthetic data creation completed successfully")


if __name__ == "__main__":
    setup_logging()
    main()
