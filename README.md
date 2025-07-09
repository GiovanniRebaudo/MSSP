# Multivariate Species Sampling Process (MSSP)

R codes for inference and prediction under multivariate species sampling process (MSSP).

**Authors**: [Beatrice Franzolini](https://beatricefranzolini.github.io) and [Giovanni Rebaudo](https://giovannirebaudo.github.io)

#### Overview 
This repository is associated with the article [Franzolini, Lijoi, Prünster and Rebaudo (2025) **Multivariate species sampling models.**  arXiv:2503.24004]

The key contribution of the paper is outlined below.
 
> Species sampling models provide a general framework for random discrete distributions and are tailored for exchangeable data. However, they fall short when used to model heterogeneous data collected from related sources or distinct experimental conditions. To address this limitation, partial exchangeability serves as the ideal probabilistic invariance condition. While numerous models exist for partially exchangeable observations, a unifying framework, similar to species sampling models, is currently absent. In this paper, we introduce multivariate species sampling models, which are a general class of models characterized by their partially exchangeable partition probability function. These models encompass existing nonparametric models for partial exchangeable data, thereby highlighting their core distributional properties and induced learning mechanisms. Our results enable an in-depth comprehension of the induced dependence structure as well as facilitate the development of new models.

This repository provides codes to replicate the results in Franzolini, Lijoi, Prünster and Rebaudo (2025) **Multivariate species sampling models.**  arXiv:2503.24004.**

More precisely, we provide the `R` code to implement inference and prediction under some relevant multivariate species sampling process (MSSP).

The repository contains the following:

1. `MAB_main_simul.R` code to reproduce the synthetic-data analysis results in the article;
2. `MAB_main_trees.R` code to reproduce the real-data analysis results in the article;
3. `Functions-and-Extra` functions and code needed to perform the analysis in `MAB_main_simul.R` and `MAB_main_trees.R`;
4. `Data-and-Results` folder with data and results of the analyses.

#### Questions or bugs
For bug reporting purposes, e-mail [Beatrice Franzolini](https://beatricefranzolini.github.io) (franzolini@pm.me) and [Giovanni Rebaudo](https://giovannirebaudo.github.io) (giovanni.rebaudo@unito.it).

#### Citation
Please cite the following publication if you use this repository in your research: [Franzolini, Lijoi, Prünster and Rebaudo (2025) **Multivariate species sampling models.**  arXiv:2503.24004]




