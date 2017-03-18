# Welcome to `Deleterious Alleles`
This is a research repo for our recent project "**incomplete dominance of deleterious alleles contribute substantially to trait variation and heterosis in maize**". The manuscript could be found via [bioRxiv]().

## Introduction
In this study, we take advantage of the genetic and genomic tools available in maize to investigate the contribution of deleterious alleles to phenotypic variation. And then we use the deleterious annotation to inform our genomic prediction model to improve the prediction accuracy for phenotypic traits and heterosis.

## Architecture about this Repo
This project contains ~400 commits. A `largedata` directory was intentionally ignored by `github` because of large size of the files. To guide the visitors having a better sense about the project, here we briefly introduce the functions or sepecific purposes of the directories. The layout of directories is based on the idea from [ProjectTemplate](http://projecttemplate.net/architecture.html). 

1. **cache**: Here we store intermediate data sets that are generated during a preprocessing step.
2. **data**: Here we store our raw data of small size. Data of large size, i.e. > 100M, store in a `largedata` folder that has been ignored using `gitignore`.
3. **doc**: Documentation codes (i.e. Rmd files) for generating the figures.
4. **graphs**: Graphs produced during the study.
5. **lib**: Some functions for our work.
6. **munge**: Here we store some preprocessing or data munging codes.
7. **profilling**: Analysis scripts for the project. It contains some sub-directories.

## Figures
Rmd code to generate Figures in the paper. The `doc/` fold is in `gitignore` at the moment. It will be released soon.

1. **Figure 1**: Panels of this figure can be produced with code [`doc/Fig1.Rmd`](https://github.com/yangjl/pvpDiallel/blob/master/doc/Fig1.Rmd). And the [Figure1](https://github.com/yangjl/pvpDiallel/blob/master/doc/Fig1.pdf) in a pdf version.
2. **Figure 2**: Panels of this figure can be produced with code [`doc/Fig2.Rmd`](https://github.com/yangjl/pvpDiallel/blob/master/doc/Fig2.Rmd). And the [Figure2](https://github.com/yangjl/pvpDiallel/blob/master/doc/Fig2.pdf) in a pdf version.
3. **Figure 3**: Panels of this figure can be produced with code [`doc/Fig3.Rmd`](https://github.com/yangjl/pvpDiallel/blob/master/doc/Fig3.Rmd). And the [Figure3](https://github.com/yangjl/pvpDiallel/blob/master/doc/Fig3.pdf) in a pdf version.

## License
This repo is free and open source for research usage, licensed under [GPLv2](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html).

