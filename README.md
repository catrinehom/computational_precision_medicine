Optimizing models for subtype prediction in Lung Adenocarcinoma
====================

Description
--------------------

This is a project on optimizing models for subtype prediction in Lung Adenocarcinoma

This has been done by: 

- Varying the number of genes in the subtype signature

- Using three different models: 
  * K Nearest Neighbors (KNN) (where K=17)
  * Distance-To-Centroid (DTC)
  * single sample Gene Set Enrichment Analysis (ssGSEA)


Data
--------------------


The data contains 2 tables: 
- FPKM data set
- Clinical data set

The data was originally found on [UCSC Xeno](https://xenabrowser.net/datapages/?cohort=GDC%20TCGA%20Lung%20Adenocarcinoma%20(LUAD)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443):

Dependencies
--------------------
- [R](https://cran.r-project.org/bin/windows/base/) >= 3.6.3, including packages:
  * tidyverse
  * GSVA


Installation
--------------------
The following code download this data analysis pipeline, when run in the terminal:

```
git clone https://github.com/catrinehom/computational_precision_medicine.git
```


Contributors
--------------------

* Catrine Høm (catrinehom)
* Julie Zimmermann (julzim)
* Oriade Simpson (s172084)
* Paolo Federico (Drivahah)

