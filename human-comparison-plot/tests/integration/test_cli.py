import typing as t
import subprocess
import shlex
import shutil
from pathlib import Path

import fur2cosmic
from fur2cosmic import constants

MODULE_NAME = fur2cosmic.__name__
PROGRAM_NAME = constants.PROGRAM_NAME

# HELPERS


def get_subprocess_message(subproces_result: subprocess.CompletedProcess) -> str:
    indent = " " * 2

    msg = (
        f"Error running CLI command. "
        f"{indent}Command: {subproces_result.args}\n"
        f"{indent}Return code: {subproces_result.returncode}\n"
        f"{indent}Stdout: {subproces_result.stdout!r}\n"
        f"{indent}Stderr: {subproces_result.stderr!r}"
    )
    return msg


# TESTS


def test_cli(
    synthetic_cmc_file: Path, maf_files_from_env: t.List[Path], tmp_path: Path
):
    # Given
    gene_symbol = "TP53"
    cosmic_cmc_file = synthetic_cmc_file
    maf_files = maf_files_from_env
    output_dir = tmp_path / "output"

    # Usage: fur2cosmic [-h] [-o OUTPUT_DIR] [--skip-plots] [--use-json] [--version] GENE-SYMBOL CMC-FILE MAF-FILE [MAF-FILE ...]
    # So we need to pass the gene symbol, cosmic file and maf files
    cmd = f"{PROGRAM_NAME} {gene_symbol} {cosmic_cmc_file} {' '.join(map(str, maf_files))} -o {output_dir}"

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg

    # Finally
    # Check if the output directory exists and contains a tsv, a png and a pdf file
    assert output_dir.exists()
    assert output_dir.is_dir()
    assert len(list(output_dir.glob("*.tsv"))) == 1
    assert len(list(output_dir.glob("*.png"))) == 1
    assert len(list(output_dir.glob("*.pdf"))) == 1


def test_cli__on_path():
    # When
    result = shutil.which(PROGRAM_NAME)

    assert result is not None, f"{PROGRAM_NAME} is not in PATH, has the name changed?"


def test_cli__version():
    # Given
    cmd = f"{PROGRAM_NAME} --version"
    expected_version = fur2cosmic.__version__

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout
    assert expected_version in subproces_result.stdout


def test_cli__help():
    # Given
    cmd = f"{PROGRAM_NAME} --help"

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout


def test_python_dash_m__version():
    # Given
    cmd = f"python -m {MODULE_NAME} --version"
    expected_version = fur2cosmic.__version__

    # When
    subproces_result = subprocess.run(shlex.split(cmd), capture_output=True, text=True)

    # Then
    errmsg = get_subprocess_message(subproces_result)
    assert subproces_result.returncode == 0, errmsg
    assert PROGRAM_NAME in subproces_result.stdout
    assert expected_version in subproces_result.stdout
