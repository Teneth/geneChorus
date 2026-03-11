An R based method for for determining the underlying correlative networks of genes in single cell data. Uses progressive rounds of correlation at separate quality levels to distinguish neighborhoods, the final version of which is termed Visualization Of eXpression groups (VOX), in contrast to TF-node based models. All genes will be classified using this method but any gene with <0.1 R<sup>2</sup> to its home network can be filtered as noise. Also generates gene lists for ontology and numerous data displays.


Here are the dependencies, everything should be on CRAN. This can easily work without Seurat if you splice in your dataframe and cell metadata manually, but the function is not currently configured for those inputs.
```
library(Seurat)
library(ggplot2)
library(dplyr)
library(randomcoloR)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
```


How to run the workflow:

```
source("2026.02.11_geneCHORUS_v12.6_function_script_part1.R")

gene_network_table <- CreateGeneChorusReport(Object_data = your_rds,    ### Provide a seurat object with normalized data in @assays$RNA$data
                                      output.dir.location = "/u/project/location",
                                      cluster_search = "test",      ## "fast" or "test", test will try to guess the right amount of clusters wanted and may be slow, fast will use an approximate cell number scalar
                                      Exp.Name = "ExpName",       ## Just a title for the data
                                   noise_floor=0.1){      ##Minimum floor for genes to be paired into network. Filteration step, raw networks reported before this is triggered 
                                        

```
