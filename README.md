# modelling
Jan 27, 2025

[![](https://www.repostatus.org/badges/latest/wip.svg)](https://www.repostatus.org/#wip)
[![lint](https://github.com/Bioforest-project/modelling/workflows/lint/badge.svg)](https://github.com/Bioforest-project/modelling/actions?query=workflow%3Alint)

**modelling** is a sub-project of the
[**BioForest**](https://github.com/Bioforest-project) project aimed at
species related data (correct names, functional traits, phylogeny) as
part of the data preparation for data preparation within the project.

## Usage

All **modelling** analyses rely on the quarto documents (`files.qmd`)
that can be run with R and associated environment defined with
[renv](#0).

## Project

**modelling** includes:

- Analyse of the data with associated [documentation and
  figures](https://bioforest-project.github.io/modelling/):
  - Reproductive analyses in `files.qmd`
  - Resulting pages in `docs/`
  - Document structure definition in `_quarto.yml`
- Automatic per site reports for PIs `analyses/`
- Data in `data/` with:
  - All raw data in `raw_data/`
  - All derived data in `derived_sata/`
- Stan Bayesian models in `models/`
- Stan Bayesian models outputs in `chains/`
- R scripts with funtions in `r/`
- Intermediary files in `outputs/`
- Figures in `figures/`
- R environment definition with
  [renv](https://rstudio.github.io/renv/articles/renv.html) in `renv/`
  and `renv/lock`
- R files (`.Rbuildignore` , `.Rdata` , `.Rprofile` , `.Rhistory`,
  `.lintr`)
- Git and GitHub files (`.gitignore` , `.github/`)
- Project documentation (`README.qmd` , `README.md` , `NEWS.md`,
  `LICENSE` )

## Contribution

You can contribute to the project by forking the repository on github
and cloning the fork to your machine using several options, including
GitHub desktop GUI. Further informations on contribution are detailed in
the online document:
<https://bioforest-project.github.io/data_preparation/98_contributing.html>.

## Help

Please preferentially create an issue on GitHub for any questions, bugs
or help needed regarding **modelling**:
<a href="https://github.com/Bioforest-project/species/issues"
class="uri">https://github.com/Bioforest-project/modelling/issues</a> .
You may however reach us by mail with people from the core group (see
below).

## Core group

- Sylvain Schmitt (sylvain.schmitt@cirad.fr)
- Camille Piponiot-Laroche (camille.piponiot-laroche@cirad.fr)
- Bruno Hérault (bruno.herault@cirad.fr)

The whole group consist of participants to the [Bioforest
project](https://www.fondationbiodiversite.fr/la-frb-en-action/programmes-et-projets/le-cesab/bioforest/).

![](https://www.fondationbiodiversite.fr/wp-content/uploads/2023/10/bioforest-ws1_web.jpeg)
