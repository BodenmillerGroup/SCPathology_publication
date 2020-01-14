[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3518284.svg)](https://doi.org/10.5281/zenodo.3518284)
# SCPathology_publication
The Single-Cell Pathology Landscape of Breast Cancer Reveals Intra-Tumor Cellular Organization and Novel Subgroups

This repository contains all code used to produce the results and figures of the publication "The Single-Cell Pathology Landscape of Breast Cancer". All data, including tiff images, masks, single-cell and patient data are available on Zenodo (10.5281/zenodo.3518284).

## Matlab scripts:
Image and other early analysis steps were performed using Matlab. Since the single-cell data was extracted using histoCAT, the Matlab scripts assume a data structure as in a loaded histoCAT session. Saved histoCAT sessions can be downloaded from Zenodo (10.5281/zenodo.3518284).

## R scripts:
Downstream analysis was performed using R pipelines. The R analysis is divided into one notebook for the analysis of the first TMA of 281 patients from University Hospital Basel and a second one for comparison and analysis of the second multi-core cohort from Univerity Hospital Zurich. All input data required to reproduce the figures of this publication are available on Zenodo(10.5281/zenodo.3518284). The BaselTMA and ZurichTMA folders contain the input data for the respective R pipelines.

## Data organization on Zenodo:
OMEandSingleCellMasks.zip contains the ome-tiff stacks and the single-cell masks.
TumorStroma_masks.zip contains masks for tumor and stromal regions.
SingleCell_and_Metadata.zip contains the single-cell and patient data as well as all other input data for the R pipelines provided here.

| Where to find:                            | Subpath                                                     |
| ----------------------------------------- | ----------------------------------------------------------- |
| Patient and core metadata BaselTMA        | SingleCell_and_Metadata/BaselTMA/Basel_PatientMetadata.csv  |
| Patient and core metadata ZurichTMA       | SingleCell_and_Metadata/ZurichTMA/Zuri_PatientMetadata.csv  |
| Single-cell data BaselTMA                 | SingleCell_and_Metadata/BaselTMA/SC_dat.csv                 |
| Single-cell data ZurichTMA                | SingleCell_and_Metadata/ZurichTMA/SC_dat.csv                |
| Single-cell segmentation masks both TMAs  | OMEandSingleCellMasks/Basel_Zuri_masks/                     |
| Image tiffs both TMAs                     | OMEandSingleCellMasks/ome/                                  |
| Antibody panel                            | SingleCell_and_Metadata/Basel_Zuri_StainingPanel            |

### Important notes when working with the data provided on Zenodo: 
- The single-cell data provided for downstream R analysis is already spillover corrected.
- The single-cell masks that were generated using CellProfiler do not always have sequential single-cell labels. Every now and then an ID is skipped. This can cause issues in histoCAT and therefore the single cells are automatically relabelled sequentially during loading into histoCAT. We exported the single-cell data for downstream R analysis from histoCAT and therefore the single-cell labels are the newly assigned sequential ones and match the labels in the histoCAT sessions. However, the original mask files that are also provided here still contain the original labels from CellProfiler. Therefore, for matching the single-cell data provided here directly to the masks (e.g. for visualization of single-cell features on the image outside of histoCAT), the single-cell labels in the mask need to be relabelled sequentially first.
