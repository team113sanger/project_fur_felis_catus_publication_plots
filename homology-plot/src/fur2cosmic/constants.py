PROGRAM_NAME = "fur2cosmic"

COL_OUTPUT__FELINE_NUC_POS = "Feline_Nuc_Pos"
COL_OUTPUT__HUMAN_NUC_POS = "Human_Nuc_Pos"
COL_PATIENT_ID = "Patient_ID"
COL_OUTPUT__DNA_CHANGE_ORIGINAL = "Original HGVS DNA Change"
COL_OUTPUT__PROTEIN_CHANGE_ORIGINAL = "Original HGVS Protein Change"
COL_OUTPUT__DNA_CHANGE_AS_HUMAN = "Humanized HGVS DNA Change"
COL_OUTPUT__PROTEIN_CHANGE_AS_HUMAN = "Humanized HGVS Protein Change"
COL_OUTPUT__STUDY = "Study"
COL_OUTPUT__TUMOUR_TYPE = "Tumour_Type"
ENSEMBL_SERVER = "https://rest.ensembl.org"

# -----------------------------
# Study ID to Tumour Type Mapping
# -----------------------------

STUDY_ID_TO_TUMOUR_TYPE_DICT = {
    "6555_2711": "LUCA",
    "6711_2820": "cSCC",
    "6712_2822": "oSCC",
    "6713_2821": "cMCT",
    "6841_2964": "MEN",
    "6864_2965": "PANC",
    "6945_3142": "CHOL",
    "6973_2987": "OSA",
    "6982_3135": "LYM",
    "6990_3065": "MAM",
    "7040_3064": "BCC",
    "7097_3073": "CRC",
    "7098_3140": "GLIO",
}


# Use the keys from STUDY_ID_TO_TUMOUR_TYPE_DICT to build a fixed global color mapping


def _set_global_study_colors():
    import matplotlib as mpl
    import matplotlib.pyplot as plt

    all_study_ids = sorted(list(STUDY_ID_TO_TUMOUR_TYPE_DICT.keys()), reverse=True)
    num_colors = len(all_study_ids)
    # Get the colormap
    cmap = plt.colormaps["tab20"]
    # Create a ListedColormap with exactly num_colors evenly spaced colors from the original colormap
    colors = [
        cmap(i / (num_colors - 1)) if num_colors > 1 else cmap(0)
        for i in range(num_colors)
    ]
    _discrete_cmap = mpl.colors.ListedColormap(colors)

    return {study: _discrete_cmap(i) for i, study in enumerate(all_study_ids)}


GLOBAL_STUDY_COLORS = _set_global_study_colors()
