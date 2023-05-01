# Omix Insight R

## Introduction
OmixInsightR is a package to analyze and visualize multi-omics data. It can be used to analyze and visualize the data with cBioportal data format. It contains the functions to analyze and visualize the data, such as the mutation, copy number variation, mRNA expression, microRNA expression, protein expression, DNA methylation, and clinical data.

## Installation
```{r}
# Install devtools if necessary
install.packages("devtools")

# Install OmixInsightR from github
devtools::install_github("biominer-lab/omix-insight-r")
```

## Usage
```{r}
library(OmixInsightR)
```

## For Developers

### How to organize the code

1. import-* functions in the R folder
    import-* function is used to import the data from the cBioportal data format. The data format is described in the [cBioportal website](https://docs.cbioportal.org/file-formats/#introduction). The import-* function should return an object of the class OmixInsightData. The OmixInsightData class is defined in the R/omix-insight-data.R file. The OmixInsightData class contains the following slots:
    - study: a data frame containing the study information
    - metadata: a list containing the metadata information. all metadata files from cBioportal data are stored in the list. The name of the list is the name of the metadata file. The content of the list is a data frame containing the metadata information.
    - caseList: a data frame containing the case list information (Row: sample ID, Column: data type; 1: included, 0: excluded; The first column is the sample ID, and the second column is the patient ID, and the rest columns are the data types)
    - clinicalData: a list which contains two data frames. The first data frame is the clinical data, and the second data frame is the clinical data dictionary which contains annotations from the first four rows of the clinical data file.
    - sampleData: a list which contains two data frames. The first data frame is the sample data, and the second data frame is the sample data dictionary which contains annotations from the first four rows of the sample data file.
    - cnaData: a data frame containing the copy number variation data (Row: gene ID, Column: sample ID, Value: copy number variation value)
    - mutationData: a data frame is same as the MAF file from cBioportal data.
    - mrnaData: a list which contains metadata and data. The metadata is a data frame containing the metadata information. The data is a data frame containing the mRNA expression data (Row: gene ID, Column: sample ID, Value: mRNA expression value)
    - mirnaData: a list which contains metadata and data. The metadata is a data frame containing the metadata information. The data is a data frame containing the microRNA expression data (Row: gene ID, Column: sample ID, Value: microRNA expression value)


2. export-* functions in the R folder
3. plot-* functions in the R folder
4. transform-* functions in the R folder