import importlib.resources
from pathlib import Path


def get_non_ensembly_homology_json() -> Path:
    """
    Retrieve the path to the non-ensembly homology JSON file.

    In some situations Ensembl doesn't recognise a gene from two species as
    homologous, as they do not appear in the same Ensembl gene tree.

    However for CDKN2A, our feline baitset uses the human CDKN2A gene with
    homology previously established. This JSON corresponds to clustal2
    alignments of the CDKN2A gene from human and feline.
    """
    # Use importlib.resources to access the package data
    import fur2cosmic.data

    file_name = "gapped_cdna.json"
    json_as_traversible = importlib.resources.files(fur2cosmic.data) / file_name
    json_as_path = Path(str(json_as_traversible)).resolve()
    if not json_as_path.exists():
        raise FileNotFoundError(
            f"JSON file {file_name} not found in package data directory."
        )
    return json_as_path
