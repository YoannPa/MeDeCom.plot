# MeDeCom.plot - Enhanced plots for MeDeCom results

![GitHub repo size](https://img.shields.io/github/repo-size/YoannPa/MeDeCom.plot)
![GitHub issues](https://img.shields.io/github/issues-raw/YoannPa/MeDeCom.plot)
![GitHub closed issues](https://img.shields.io/github/issues-closed-raw/YoannPa/MeDeCom.plot)  

_**MeDeCom.plot** provides MeDeCom user with additionnal plots that they can use to display results from a MeDeComSet._  

**Author: PAGEAUD Y.<sup>1</sup>**  
**How to cite:** _Pageaud Y. et al., MeDeCom.plot - Enhanced plots for MeDeCom results._  

![GitHub R package version](https://img.shields.io/github/r-package/v/YoannPa/MeDeCom.plot?label=Package%20version&logo=RStudio&logoColor=white&style=for-the-badge)  
<img src="https://img.shields.io/static/v1?label=compatibility&message=4.2.0&color=blue&logo=R&logoColor=white&style=for-the-badge" />  
![GitHub last commit](https://img.shields.io/github/last-commit/YoannPa/MeDeCom.plot?logo=git&style=for-the-badge)  
![GitHub](https://img.shields.io/github/license/YoannPa/MeDeCom.plot?color=brightgreen&style=for-the-badge)  


## Prerequisites (install MeDeCom, RnBeads, BiocompR, and Bioconductor dependencies)
In R execute the following command:
```R
devtools::install_github("CompEpigen/MeDeCom")
devtools::install_github("epigen/RnBeads")
devtools::install_github("YoannPa/BiocompR@7db9e0d7b4a653bd7338bd7b2df8a91ae7d65955") # MeDeCom.plot was developed under this version of BioCompR.
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(pkgs = c("AnnotationDbi", "Biobase", "BiocGenerics", "Biostrings", "FDb.InfiniumMethylation.hg19", "GenomeInfoDb", "GenomicFeatures", "GenomicRanges", "illuminaio", "methylumi", "minfi", "org.Hs.eg.db", "SummarizedExperiment", "TxDb.Hsapiens.UCSC.hg19.knownGene", "XVector"))
```

## Install MeDeCom.plot
In R execute the following command:
```R
devtools::install_github("YoannPa/MeDeCom.plot")
```

## Content
Currently the package MeDeCom.plot contains **40 exported functions**:

* `cve_plot3d()` - Plots MeDeComSet C.V. error following kappa and lambda values tested.  
* `LMCprop_heatmap()` - Plots a heatmap of LMCs proportions in all methylation samples.  
* `LMCs_dendrogram()` - Plots a UPGMA dendrogram clustering of LMCs with reference methylomes.  
* `mds_plot()` - Plots a MDS of LMCs and samples.  
* `rmse_plot3d()` - 	Plots MeDeComSet R.M.S. error following kappa and lambda values tested.  

## Problems ? / I need help !
For any questions **Not related to bugs or development** please check the section "**Known Issues**" available below. If the issue you experience is not adressed in the known issues you can write me at [y.pageaud@dkfz.de](y.pageaud@dkfz.de).

## Technical questions / Development / Feature request
If you encounters issues or if a feature you would expect is not available in a MeDeCom.plot function, please check if an existing issue adresses your point [here](https://github.com/YoannPa/MeDeCom.plot/issues/). If not, create a [new issue here](https://github.com/YoannPa/MeDeCom.plot/issues/new).  

## References
⚠️ **Work in progress !**  

1. [_Lutsik P. et al., MeDeCom: discovery and quantification of latent components of heterogeneous methylomes. Genome Biol. 2017 Mar 24;18(1):55. doi: 10.1186/s13059-017-1182-6. PMID: 28340624; PMCID: PMC5366155._](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1182-6)  
2. [_Pageaud Y. et al., BiocompR - Advanced visualizations for data comparison._](https://github.com/YoannPa/BiocompR)  
3. [_Package ‘lattice’_](https://cran.r-project.org/web/packages/lattice/lattice.pdf)  
4. [_3D Graph - Lattice Package_](https://myrcodes.blogspot.com/2015/10/3d-graph-lattice-package.html)  

## Licence
MeDeCom.plot is currently under the GPL-3.0 licence.  

