# SCPathology_publication
The Single-Cell Pathology Landscape of Breast Cancer Reveals Intra-Tumor Cellular Organization and Novel Subgroups

This repository contains all code used to produce the results and figures of the publication "The Single-Cell Pathology Landscape of Breast Cancer". All data, including tiff images, masks, single-cell and patient data are available on Zenodo (10.5281/zenodo.3518284).

Matlab:
Image and other early analysis steps were performed using Matlab. Since the single-cell data was extracted using histoCAT, the Matlab scripts assume a data structure as in a loaded histoCAT session. Saved histoCAT sessions can be downloaded from Zenodo (10.5281/zenodo.3518284).

R:
Downstream analysis was performed using R pipelines. The R analysis is divided into one notebook for the analysis of the first TMA of 281 patients from University Hospital Basel and a second one for comparison and analysis of the second multi-core cohort from Univerity Hospital Zurich. All input data required to reproduce the figures of this publication are available on Zenodo(10.5281/zenodo.3518284). The BaselTMA and ZurichTMA folders contain the input data for the respective R pipelines.
