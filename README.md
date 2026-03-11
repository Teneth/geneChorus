A method for for determining the underlying correlative networks of genes in single cell data. Uses progressive rounds of correlation at separate quality levels to distinguish neighborhoods, the final version of which is termed Visualization Of eXpression groups (VOX), in contrast to TF-node based models. All genes will be classified using this method but any gene with <0.1 R2 to its home network can be filtered as noise. Also generates gene lists for ontology and numerous data displays.



```
gene_network_table <- CreateGeneChorusReport(Object_data = your_rds,
                                      output.dir.location = "/u/project/location",
                                      cluster_search = "test", ## "fast" or "test", test will try to guess the right amount of clusters wanted, fast will use a tested cell number scalar
                                      Exp.Name = "ExpName",## Just a title for the data)
                                   noise_floor=0.1){   ##Minimum floor for genes to be paired into network. Filteration step, raw networks reported beore this is triggered 
                                        

```
