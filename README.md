# Feline Fur Plots

A collection of code corresponding to different publication plots.

## Code organisation

Plots are organised by subdirectory, each containing R or Python code as well
documentation, results and data files.

| Subdirectory | Description | Language |
|--------------|-------------|----------|
| `homology-plot`   | Feline-to-human comparative oncogenomics [see the README.md](homology-plot/README.md) | Python |
| `germline-plot`   | Germline variant calling across feline samples [see the README.md](germline-plot/README.md) | R |

## TODO

- [x] Create a `homology` repo
- [x] Include `homology` repo via `git_subclone.py`
- [x] Include `germline` repo via `git_subclone.py`
- [ ] Add data files to `homology` subdirectory / repo
- [ ] ? Add plots to `homology` subdirectory / repo
- [ ] Either host an `oncoplots` repo or include code directly in this repo
- [ ] Either host an `cnv_plots` repo or include code directly in this repo
- [ ] ? Add some code to handle data file duplicates

## Scripts
`git_subclone.py` - Clone a git repository at a specific tag or commit hash into this repo. While `git submodules` facilitate a similar and superior workflow it is not ideal for mirroring
repositories from Gitlab to Github nor when creating a DOI with Zenodo. This script instead creates a single git entity.