from pathlib import Path
import argparse
import sys

import fur2cosmic.main
from fur2cosmic import constants
from fur2cosmic.utils.logging_utils import setup_logging, get_package_logger

LOGGER = get_package_logger()


def setup_argparser() -> argparse.ArgumentParser:
    """
    Set up and return an argument parser for the script.

    Returns:
        argparse.ArgumentParser: Configured argument parser instance
    """
    parser = argparse.ArgumentParser(
        description="Map feline mutations to human and compare with COSMIC data.",
        prog=constants.PROGRAM_NAME,
    )
    parser.add_argument(
        "gene_symbol",
        type=str,
        metavar="GENE-SYMBOL",
        help="HUFO gene symbol to analyze, e.g., TP53",
    )
    parser.add_argument(
        "cosmic_file",
        type=Path,
        metavar="CMC-FILE",
        help="Path to a COSMIC Cancer Mutation Census (CMC) data file",
    )
    parser.add_argument(
        "maf_files",
        nargs="+",
        metavar="MAF-FILE",
        help="One or more MAF files containing feline mutations, e.g., study1.maf study2.maf",
        type=Path,
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        type=Path,
        default=Path.cwd(),
        dest="output_dir",
        help="Directory to save the output files (default: current directory)",
        required=False,
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=600,
        help="DPI for saving plots (default: 600)",
        required=False,
    )
    parser.add_argument(
        "--skip-plots", action="store_true", help="Skip generating comparison plots"
    )
    parser.add_argument(
        "--use-json",
        action="store_true",
        help="Use packaged JSON data for pre-calculated nucleotide alignments not present in Ensembl API",
    )
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {fur2cosmic.__version__}",
        help="Show the version of the package",
    )

    return parser


def cli():
    setup_logging()
    parser = setup_argparser()
    args = parser.parse_args()

    # Check if the input files exist
    inputs = [args.cosmic_file] + args.maf_files
    for input_file in inputs:
        if not input_file.exists():
            LOGGER.error(f"Input file '{input_file}' not found")
            sys.exit(1)

    # Check if the output directory exists, create it if not
    if not args.output_dir.exists():
        LOGGER.info(f"Creating output directory: {args.output_dir}")
        args.output_dir.mkdir(parents=True, exist_ok=True)

    # Run the main function from the main module
    dataframe, figure = fur2cosmic.main.main(
        args.gene_symbol,
        args.cosmic_file,
        args.maf_files,
        skip_plots=args.skip_plots,
        use_json=args.use_json,
    )

    # Save the output files in the specified directory
    output_file_path = args.output_dir / f"{args.gene_symbol}.fur2cosmic.tsv"
    LOGGER.info(f"Saving data to {output_file_path}")
    dataframe.to_csv(output_file_path, index=False, sep="\t")

    if figure is not None:
        fig_pdf_file = args.output_dir / f"{args.gene_symbol}.fur2cosmic.pdf"
        fig_png_file = args.output_dir / f"{args.gene_symbol}.fur2cosmic.png"

        LOGGER.info(f"Saving plot to {fig_pdf_file}")
        figure.savefig(fig_pdf_file, format="pdf", bbox_inches="tight", dpi=args.dpi)

        LOGGER.info(f"Saving plot to {fig_png_file}")
        figure.savefig(fig_png_file, format="png", bbox_inches="tight", dpi=args.dpi)

    LOGGER.info("Done.")
    return
