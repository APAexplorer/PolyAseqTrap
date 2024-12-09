# PolyAseqTrap
A universal tool for genome-wide identification and quantification of polyadenylation sites from different 3′ end sequencing data

About
====================
PolyAseqTrap is an R package designed to identify and quantify polyA sites from various 3′ sequencing datasets (e.g., DRS, PAT-seq, PASC-seq, 3′READS). It utilizes a polyA read prioritization strategy with detailed post-inspection to minimize false polyA site calls and accurately determine their precise locations. Notably, PolyAseqTrap incorporates a transferable, cross-species deep learning model to resolve the persistent challenge of internal priming. Furthermore, it includes a weighted density peak clustering method that considers the microheterogeneity of polyadenylation across species to define polyA site clusters (PACs). The package also provides extensive tools for annotation, validation, and visualization of polyA sites, delivering well-structured reports to facilitate seamless analysis.

* The PolyAseqTrap package consists of six main modules.

<img src="https://github.com/APAexplorer/PolyAseqTrap/blob/main/img/schema.png" alt="schema" width="1000"/>


Installing scNPF
=============
Mandatory 
---------

* R (>3.1). [R 3.5](https://www.r-project.org/) is recommended.

Required R Packages
---------
* [igraph](https://cran.r-project.org/web/packages/igraph/index.html), [WGCNA](https://cran.r-project.org/web/packages/WGCNA/index.html), [foreach](https://cran.r-project.org/web/packages/foreach/index.html), [doParallel](https://cran.r-project.org/web/packages/doParallel/index.html), [Matrix](https://cran.r-project.org/web/packages/Matrix/index.html), [plyr](https://cran.r-project.org/web/packages/plyr/index.html),  

Suggested R Packages
---------
* [gplots](https://cran.r-project.org/web/packages/gplots/index.html), [RColorBrewer](https://cran.r-project.org/web/packages/RColorBrewer/index.html),  

Installation
---------
* Install the R package using the following commands on the R console:
```
install.packages("devtools")
library(devtools)
install_github("BMILAB/scNPF")
library(scNPF)
```


