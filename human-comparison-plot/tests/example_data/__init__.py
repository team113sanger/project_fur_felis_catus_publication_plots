import importlib.resources
from pathlib import Path


def get_synthetic_COSMIC_CMC_taster_tar_gz() -> Path:
    """
    Retrieve the path to the synthetic COSMIC CMC file, filtered for TP53 mutations only.

    This file is used for testing purposes and contains synthetic and anonymized
    data based on the original COSMIC CMC file.

    The file was generated by running while in a virtual environment:
    python tests/create_synthetic_cmc.py Cosmic_MutantCensus_v101_GRCh38.tsv Synthetic_CMC_taster.TP53.tsv --genes TP53
    gzip -9 Synthetic_CMC_taster.TP53.tsv
    mv Synthetic_CMC_taster.TP53.tsv.gz tests/example_data/

    You will need to acquire the original COSMIC CMC file from the COSMIC website (this will
    require registration):
    https://cancer.sanger.ac.uk/cmc/home
    """
    # Use importlib.resources to access the package data
    import tests.example_data

    file_name = "Synthetic_CMC_taster.TP53.tsv.gz"
    tsvgz_as_traversible = importlib.resources.files(tests.example_data) / file_name
    tsvgz_as_path = Path(str(tsvgz_as_traversible)).resolve()
    if not tsvgz_as_path.exists():
        raise FileNotFoundError(
            f"TSV file {file_name} not found in package data directory."
        )
    return tsvgz_as_path


def get_truncacted_maf_file() -> Path:
    """
    Retrieve the path to an example MAF file that only has CATD252a mutations.
    """
    # Use importlib.resources to access the package data
    import tests.example_data

    maf_file_name = "example.CATD252a.maf"
    maf_as_traversible = importlib.resources.files(tests.example_data) / maf_file_name
    maf_as_path = Path(str(maf_as_traversible)).resolve()
    if not maf_as_path.exists():
        raise FileNotFoundError(
            f"MAF file {maf_file_name} not found in package data directory."
        )
    return maf_as_path


def get_example_maf_gz_files() -> list[Path]:
    """
    Retrieve the paths to non-truncated, compressed MAF files.
    """
    # Use importlib.resources to access the package data
    import tests.example_data

    maf_file_names = [
        "example_1.maf.gz",
        "example_2.maf.gz",
    ]
    maf_as_paths = []
    for maf_file_name in maf_file_names:
        maf_as_traversible = (
            importlib.resources.files(tests.example_data) / maf_file_name
        )
        maf_as_path = Path(str(maf_as_traversible)).resolve()
        if not maf_as_path.exists():
            raise FileNotFoundError(
                f"MAF.GZ file {maf_file_name} not found in package data directory."
            )
        maf_as_paths.append(maf_as_path)
    return maf_as_paths
