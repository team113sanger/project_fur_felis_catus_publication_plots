import typing as t
from pathlib import Path

from fur2cosmic import main
import fur2cosmic.constants as C


# TESTS


def test_main__env_mafs(
    synthetic_cmc_file: Path, maf_files_from_env: t.List[Path]
) -> None:
    # Given
    gene_symbol = "TP53"
    cosmic_file = synthetic_cmc_file
    skip_plots = False
    use_json = False

    # When
    dataframe, figure = main.main(
        gene_symbol=gene_symbol,
        cosmic_file=cosmic_file,
        maf_files=maf_files_from_env,
        skip_plots=skip_plots,
        use_json=use_json,
    )

    # Then
    assert dataframe is not None
    assert figure is not None

    assert len(dataframe) > 0
    assert len(dataframe.columns) > 0


def test_main__example_mafs(
    synthetic_cmc_file: Path, maf_files_from_example_data: t.List[Path]
) -> None:
    # Given
    gene_symbol = "TP53"
    cosmic_file = synthetic_cmc_file
    skip_plots = False
    use_json = False

    # When
    dataframe, figure = main.main(
        gene_symbol=gene_symbol,
        cosmic_file=cosmic_file,
        maf_files=maf_files_from_example_data,
        skip_plots=skip_plots,
        use_json=use_json,
    )

    # Then
    assert dataframe is not None
    assert figure is not None

    assert len(dataframe) > 0
    assert len(dataframe.columns) > 0


def test_main__with_one_specific_maf(
    synthetic_cmc_file: Path, CATD252a_maf_file: Path
) -> None:
    # Given
    gene_symbol = "TP53"
    cosmic_file = synthetic_cmc_file
    skip_plots = False
    use_json = False

    expected_patient_id = "CATD252a"
    expected_feline_nuc_mutation = "c.988T>G"
    expected_humanized_nuc_mutation = "c.700T>G"

    # When
    dataframe, figure = main.main(
        gene_symbol=gene_symbol,
        cosmic_file=cosmic_file,
        maf_files=[CATD252a_maf_file],
        skip_plots=skip_plots,
        use_json=use_json,
    )

    # Then
    assert dataframe is not None
    assert figure is not None

    assert len(dataframe) > 0
    assert len(dataframe.columns) > 0

    # Finally
    # We create a mask with the patient id and feline mutation then we check if
    # the humanized mutation is the same
    col_patient = C.COL_PATIENT_ID
    col_feline_mut = C.COL_OUTPUT__DNA_CHANGE_ORIGINAL
    col_humanized_mut = C.COL_OUTPUT__DNA_CHANGE_AS_HUMAN
    mask = (dataframe[col_patient] == expected_patient_id) & (
        dataframe[col_feline_mut] == expected_feline_nuc_mutation
    )
    assert len(dataframe[mask]) == 1
    assert (
        dataframe[mask][col_humanized_mut].values[0] == expected_humanized_nuc_mutation
    )
