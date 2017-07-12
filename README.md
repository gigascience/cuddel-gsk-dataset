# cuddel-gsk-dataset

## Background

This project involves reproducing the analysis of a metabolomics dataset 
using R and enabling its science to be as open as possible.

The dataset was generated using UPLC-MS from samples of plasma taken from 
healthy volunteers undergoing a longitudinal study to determine the human 
metabolic profiles associated with food intake and exercise. Details of the 
sample analysis by mass spectrometry and the processing of the subsequent 
data can be found in  [thesis.md](./docs/thesis.md).

Originally, the computational analysis of the mass spectrometry data was done
using a combination of Matlab and [Taverna](https://taverna.incubator.apache.org)
workflows. Since Matlab requires a paid licence for its use and it is not 
known whether the workflows are executable, the decision was made to 
reproduce the complete analysis in [R](https://cran.r-project.org).

After the metabolomics data analysis has been reproduced in R, a number of 
tools will be used to replicate the results and make the entire methodology 
to be accessible from an open science manner. In addition to this GitHub 
repository, these tools include:
 
* [protocols.io](https://www.protocols.io)
* [Docker](https://www.docker.com)
* [ISA](http://isa-tools.org)
* [Common Workflow Language](https://github.com/common-workflow-language/common-workflow-language)
* [Galaxy](https://galaxyproject.org)
* [CodeOcean](https://codeocean.com)
* [MetaboLights](http://www.ebi.ac.uk/metabolights/)
* [GigaDB](http://gigadb.org)

  
## Current status

As of 12th July 2017, reproduction of the metabolomics analysis using R is 
still underway. See [gsk.R](./R/gsk.R) for details.