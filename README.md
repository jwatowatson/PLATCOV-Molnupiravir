# Molnupiravir versus Paxlovid
 
 This repository provides the code and data underlying the molnupiravir versus ritonavir-boosted nirmatrelvir (paxlovid) versus no study drug comparison in the PLATCOV trial. The clinical paper is available as a preprint: "Molnupiravir versus ritonavir-boosted nirmatrelvir; a randomised controlled adaptive trial comparison of antiviral effects in early symptomatic COVID-19 (PLATCOV)"
 
 
## Overview

The primary three-way analysis is given in the file *Molnupiravir_analysis.qmd*. This constructs the mITT population, makes a baseline characteristics table, sets up the data for the analysis, and then plots the outputs. It also looks at fever clearance and time to symptom resolution. All models are run in *stan* and the stan code is provided in the folder *Stan_models*.

The main result is Figure in the paper:


![](/Molnupiravir_analysis_files/figure-html/Figure_main-1.png "Comparison of viral clearance in the three randomised arms").
