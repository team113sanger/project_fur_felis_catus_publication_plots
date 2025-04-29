# Feline Fur Plots

A collection of code corresponding to different publication plots.

## Code organisation

Plots are organised by subdirectory, each containing R or Python code as well
documentation, results and data files.

| Subdirectory | Description | Language |
|--------------|-------------|----------|
| `human-comparison-plot`   | Feline-to-human comparative oncogenomics [see the README.md](human-comparison-plot/README.md) | Python |
| `cna-plot` | Copy number alterations from CNVKit segmentation data [see the README.md](cna-plot/README.md) | Python |
| `germline-plot`   | Germline variant calling across feline samples [see the README.md](germline-plot/README.md) | R |
| `oncoprint-plot`   | Fur CNVKit derived oncoprints across feline samples [see the README.md](oncoprint-plot/README.md) | Python |
| `figure_1b`   | Summary comparison of mutations across cancer types [see the README.md](figure_1b/README.md) | R |

## TODO

- [x] Create a `human-comparison-plot` repo
- [x] Include `human-comparison-plot` repo via `git_subclone.py`
- [x] Include `germline-plot` repo via `git_subclone.py`
- [x] Add data files to `human-comparison-plot` subdirectory / repo
- [x] Add plots to `human-comparison-plot` subdirectory / repo
- [x] Include `oncoplots` repo
    - [ ] Oncoplots is a WIP
- [x] Code for figure 1b
- [x] Host an `cnv_plots` repo
- [x] Include `cnv_plots` repo via `git_subclone.py`


## Scripts
`git_subclone.py` - Clone a git repository at a specific tag or commit hash into this repo. While `git submodules` facilitate a similar and superior workflow it is not ideal for mirroring
repositories from Gitlab to Github nor when creating a DOI with Zenodo. This script instead creates a single git entity.