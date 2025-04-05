# This is where you may define fixtures that are shared across multiple test files.
#
# See this FAQ: https://docs.pytest.org/en/8.2.x/how-to/index.html
# Or read this section of the docs: https://docs.pytest.org/en/8.2.x/how-to/fixtures.html#scope-sharing-fixtures-across-classes-modules-packages-or-session
import typing as t
from pathlib import Path
import os
import shutil
import gzip

import pytest

import tests.example_data


# CONSTANTS
ENVVAR_MAF_DIR = "TEST_MAF_DIR"


# FIXTURES


@pytest.fixture
def synthetic_cmc_file(tmp_path: Path) -> t.Generator[Path, None, None]:
    cmc_tar_gz = tests.example_data.get_synthetic_COSMIC_CMC_taster_tar_gz()
    temp_dir = tmp_path / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Decompress the tsv.gz file into a tsv file
    extracted_file = temp_dir / "Synthetic_CMC.TP53.tsv"

    with gzip.open(cmc_tar_gz, "rb") as f_in:
        with open(extracted_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    if not extracted_file.exists():
        raise FileNotFoundError(f"Extracted file {extracted_file} not found.")

    yield extracted_file


@pytest.fixture
def maf_files_from_env() -> t.Generator[t.List[Path], None, None]:
    # Check if the environment variable is set or raise an error
    maf_dir = os.getenv(ENVVAR_MAF_DIR, "")

    # Check if the directory specified by the environment variable exists or skip the test
    conditions = [
        (
            lambda: maf_dir is None,
            f"Environment variable {ENVVAR_MAF_DIR!r} is not set. Please set it to the directory containing MAF files.",
        ),
        (
            lambda: not os.path.exists(maf_dir),
            f"Directory {maf_dir!r} does not exist.",
        ),
        (
            lambda: not os.path.isdir(maf_dir),
            f"Path {maf_dir!r} is not a directory.",
        ),
        (
            lambda: not os.access(maf_dir, os.R_OK),
            f"Directory {maf_dir!r} is not readable.",
        ),
        (
            lambda: not os.access(maf_dir, os.X_OK),
            f"Directory {maf_dir!r} is not executable.",
        ),
        (
            lambda: len(list(Path(maf_dir).glob("*.maf"))) == 0,
            f"No MAF files found in {maf_dir!r}.",
        ),
    ]
    for condition, message in conditions:
        if condition():
            pytest.skip(f"Skipping test: {message}")

    maf_files = list(Path(maf_dir).glob("*.maf"))

    yield maf_files


@pytest.fixture
def maf_files_from_example_data(
    tmp_path: Path,
) -> t.Generator[t.List[Path], None, None]:
    maf_gz_files = tests.example_data.get_example_maf_gz_files()

    # Decompress the gz files into a temporary directory
    temp_dir = tmp_path / "temp"
    temp_dir.mkdir(parents=True, exist_ok=True)
    for maf_gz_file in maf_gz_files:
        with gzip.open(maf_gz_file, "rb") as f_in:
            extracted_file = temp_dir / maf_gz_file.name.replace(".gz", "")
            with open(extracted_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        if not extracted_file.exists():
            raise FileNotFoundError(f"Extracted file {extracted_file} not found.")

    # Collect all MAF files in the temporary directory
    maf_as_paths = list(temp_dir.glob("*.maf"))
    if not maf_as_paths:
        raise FileNotFoundError(f"No MAF files found in {temp_dir}.")

    yield maf_as_paths


@pytest.fixture
def CATD252a_maf_file() -> t.Generator[Path, None, None]:
    """
    Retrieve the path to an example MAF file that only has CATD252a mutations.
    """
    # Use importlib.resources to access the package data
    import tests.example_data

    maf_file = tests.example_data.get_truncacted_maf_file()
    if "CATD252a" not in maf_file.name:
        raise ValueError(
            f"File {maf_file} does not contain CATD252a mutations. Is this the correct file?"
        )
    yield maf_file
