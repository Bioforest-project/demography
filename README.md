
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Demography <img src="figures/compendium-sticker.jpg" align="right" style="float:right; height:120px;"/>

<!-- badges: start -->
<!-- badges: end -->
<p align="left">
• <a href="#overview">Overview</a><br> •
<a href="#features">Features</a><br> •
<a href="#content">Content</a><br> •
<a href="#installation">Installation</a><br> •
<a href="#usage">Usage</a><br> • <a href="#citation">Citation</a><br> •
<a href="#contributing">Contributing</a><br> •
<a href="#acknowledgments">Acknowledgments</a><br> •
<a href="#references">References</a><br>
</p>

## Overview

**Demography** is a subproject of the Bioforest project aimed at
studying the interactions between tree biodiversity, forest dynamics and
climate in managed tropical forests with a pan-tropical approach. The
subproject **Demagraphy** aims to investigate the role of demographic
processes in explaining post-logging changes in biodiversity and
composition. Its main research questions are:

- What is the contribution of demographic processes (growth,
  recruitment, mortality) to post-logging changes in biomass,
  biodiversity and composition?

- What is the link between biomass and biodiversity recovery processes?

- How do abiotic factors and pre-logging biodiversity shape demographic
  processes of recovery?

The core group is composed of Camille Piponiot, Géraldine Derroire,
Euridice Honorio, Dwinda Putri, and Irié Zo-Bi.

## Content

This repository is structured as follow:

- [`DESCRIPTION`](https://github.com/Bioforest-project/Demography/tree/main/DESCRIPTION):
  contains project metadata (authors, date, dependencies, etc.)

- [`make.R`](https://github.com/Bioforest-project/Demography/tree/main/make.R):
  main R script to run the entire project

- [`R/`](https://github.com/Bioforest-project/Demography/tree/main/R):
  contains R functions developed especially for this project

- [`analyses/`](https://github.com/Bioforest-project/Demography/tree/main/analyses):
  contains scripts used for data analysis and for creating figures and
  outputs

- [`data/`](https://github.com/Bioforest-project/Demography/tree/main/data):
  contains additional data used for the analysis

- [`figures/`](https://github.com/Bioforest-project/Demography/tree/main/figures):
  contains figures

- [`outputs/`](https://github.com/Bioforest-project/Demography/tree/main/outputs):
  contains additional outputs such as tables

## Installation

To install this compendium:

- [Fork](https://docs.github.com/en/get-started/quickstart/contributing-to-projects)
  this repository using the GitHub interface.
- [Clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository)
  your fork using `git clone fork-url` (replace `fork-url` by the URL of
  your fork). Alternatively, open [RStudio
  IDE](https://posit.co/products/open-source/rstudio/) and create a New
  Project from Version Control.

## Usage

Launch the
[`make.R`](https://github.com/cpiponiot/Demography/tree/main/make.R)
file with:

    source("make.R")

**Notes**

- All required packages listed in the `DESCRIPTION` file will be
  installed (if necessary)
- All required packages and R functions will be loaded
- Some analyses listed in the `make.R` might take time

## Citation

## Contributing

All types of contributions are encouraged and valued. For more
information, check out our [Contributor
Guidelines](https://github.com/cpiponiot/Demography/blob/main/CONTRIBUTING.md).

Please note that this project is released with a [Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.

## Acknowledgments

## References
