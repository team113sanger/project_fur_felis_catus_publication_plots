# The domestic cat oncogenome - individual analysis scripts

This repository serves as sub-collection of individual analysis scripts & processes that could not
be classified in the same way as the parent repository [domestic-cat-oncogenome (github)](https://github.com/team113sanger/The-pan-cancer-oncogenomic-landscape-of-domestic-cats-reveals-shared-alterations-with-human-cancer).

These analyses are part of the manuscript:
> **_Francis, B., Ludwig, L. et al 2025 - The domestic cat oncogenome_**


## Code organisation

Plots are organised by subdirectory, each containing R or Python code as well
documentation, dependencies, results and data files.

| Subdirectory | Description | Language |
|--------------|-------------|----------|
| `human-comparison-plot`   | Feline-to-human comparative oncogenomics [see the README.md](human-comparison-plot/README.md) | Python |
| `cna-plot` | Copy number alterations from CNVKit segmentation data [see the README.md](cna-plot/README.md) | Python |
| `germline-plot`   | Germline variant calling across feline samples [see the README.md](germline-plot/README.md) | R |
| `oncoprint-plot`   | Fur CNVKit derived oncoprints across feline samples [see the README.md](oncoprint-plot/README.md) | Python |
| `figure_1b`   | Summary comparison of mutations across cancer types [see the README.md](figure_1b/README.md) | R |
| `penetrance-plot` | Visualizing genome-wide DNA copy number alteration penetrance [see the README.md](penetrance-plot/README.md) | Python |
| `tmb-plot` | Tumour mutational burden across feline samples [see the README.md](tmb-plot/README.md) | R |
| `mutational-signature-plot` | Mutational signature categorization and analysis across feline samples [see the README.md](mutational-signature-plot/README.md) | Python |

## Maintenance script
`git_subclone.py` - Clone a git repository at a specific tag or commit hash into
this repo. While `git submodules` facilitate a similar and superior workflow it
is not ideal for mirroring repositories from the internal GitLab to Github nor
when creating a DOI with Zenodo. This script instead creates a single git
entity.

For provenance and versioning of each subdirectory, inspect the `REPO_SERIES` list inside `git_subclone.py`.