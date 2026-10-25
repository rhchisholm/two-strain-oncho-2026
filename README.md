# two-strain-oncho-2026

This repository contains the model code, analysis scripts, and processed data associated with the manuscript:

**Modelling the impact of mass drug administration with ivermectin or moxidectin on human onchocerciasis when there is sub-optimal response of parasites to treatment**

Rebecca H. Chisholm, Shilian Xu, Kwadwo K Frempong, Himal Shrestha, Joseph Kwadwo Larbi Opare, Odame D. Asiedu, Ernest Mensah, Makedonka Mitreva, Warwick N Grant, Shannon M Hedtke

## Overview

This study uses a mathematical model of human onchocerciasis transmission to investigate how sub-optimal response (SOR) of *Onchocerca volvulus* to treatment may affect the impact of mass drug administration (MDA).

The model distinguishes between:

- **GR worms** — parasites with a good response to ivermectin; and
- **SOR worms** — parasites with a sub-optimal response to ivermectin.

The analyses examine the impact of annual and biannual MDA with ivermectin (MDAi) or moxidectin (MDAm), including how outcomes depend on baseline endemic transmission, treatment frequency, the presence of SOR worms, and assumptions used to calibrate the model.

The model is calibrated to epidemiological and MDAi coverage data from Asubende, Ghana.

## Repository structure

The repository contains code for:

- calibration of the variable-response model in which GR and SOR worms can coexist;
- calibration of the GR-only model;
- Latin hypercube sampling (LHS) stability and calibration diagnostics;
- sensitivity analyses examining alternative calibration assumptions;
- simulation of annual and biannual MDAi and MDAm intervention scenarios;
- generation of manuscript and supplementary figures; and
- generation of supplementary tables.

Detailed README files within the relevant analysis folders provide additional information about individual scripts and their required inputs and outputs.

## Reproducibility

The full calibration and intervention analyses generate large simulation files containing individual model trajectories. These files are not included in the repository because of their size.

Instead, the repository contains compact processed datasets sufficient to reproduce the manuscript figures and tables without rerunning the computationally intensive model simulations.

Two workflows are therefore provided.

### Reproducing the manuscript figures and tables

The manuscript figures and tables can be reproduced directly from the compact processed data included in the repository.

The full workflow also includes the model calibration and intervention simulation scripts. These analyses can be used to regenerate the large intermediate simulation files that are omitted from the repository.

## Model calibration

The variable-response model is calibrated in two stages.

### Stage 1

The three parameters governing inheritance of the parasite response phenotype are sampled using Latin hypercube sampling. Parameter combinations are retained if they produce a plausible pre-MDA proportion of SOR worms and satisfy the specified inheritance-parameter ordering.

### Stage 2

The retained Stage 1 parameter combinations are combined with LHS samples of parameters governing the response of worms to ivermectin.

The primary calibration uses the 20-year annual MDAi program in Asubende and retains parameter combinations consistent with:

- final mf prevalence of 10–12%; and
- an increase in SOR worm proportion greater than 0.10 during the annual MDAi program.

Additional constraints are applied to the inheritance parameters following Stage 2 calibration.

The GR-only model is calibrated separately to the final mf prevalence following the same annual MDAi program.

## Calibration sensitivity analyses

Four sensitivity analyses examine alternative assumptions used during model calibration:

1. **Moderate increase** — an increase in SOR worm proportion of 0.05–0.10.
2. **Weak increase** — an increase in SOR worm proportion of 0–0.05.
3. **No offspring-phenotype selection** — offspring-phenotype selection is removed, with no calibration criterion imposed on the change in SOR worm proportion.
4. **Wider prevalence** — the accepted final mf prevalence range is widened from 10–12% to 8–14%.

The sensitivity analyses are used to assess the robustness of model outcomes to uncertainty in the calibration assumptions.

## Intervention scenarios

The calibrated models are used to simulate MDA programs under different:

- baseline endemic mf prevalence levels;
- treatment frequencies;
- treatment drugs (ivermectin or moxidectin); and
- assumptions about the presence or absence of SOR worms.

For comparisons of annual and biannual MDA with the same number of treatment rounds, 20 MDA rounds correspond to 20 years of annual treatment and 10 years of biannual treatment.

Outcomes include final mf prevalence, time to the provisional operational threshold for treatment interruption and commencement of surveillance (pOTTIS), and changes in the proportion of SOR worms.

## Terminology

Some historical code objects and filenames use the term `posteriors`.

These objects contain parameter combinations retained because they satisfy specified biological and epidemiological calibration criteria. They are **not Bayesian posterior samples**.

Accordingly, the manuscript and current documentation refer to these as **retained parameter distributions** or **retained parameter combinations**.

Reported 2.5th–97.5th percentile intervals describe the distributions of retained parameter values or model outcomes and should not be interpreted as Bayesian credible intervals.

## Software

The analyses were performed in **R**.

The scripts use standard R packages including packages for:

- Latin hypercube sampling;
- data manipulation;
- plotting; and
- combining figure panels.

Individual scripts load the packages required for the corresponding analysis.

## Contact

For questions about the model or repository, please contact the corresponding authors:

**Rebecca H. Chisholm**  
La Trobe University

**Shannon M. Hedtke**  
La Trobe University