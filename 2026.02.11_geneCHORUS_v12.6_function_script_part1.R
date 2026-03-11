# echo "R CMD BATCH /u/project/kp1/jlangerm/Projects/IGVF/R_Analysis/Analysis/2025.09.11_YS3v6_geneCHORUS_v12.4_script_part1.R" | qsub -cwd -V -N gc_ys3_pt11 -l h_rt=18:00:00,h_data=4G -pe shared 16



library(Seurat)
library(ggplot2)
library(dplyr)
library(randomcoloR)
library(reshape2)
library(pheatmap)
library(RColorBrewer)
# library(ggsankey)


# #### parameters
# 
# your_rds <-  readRDS("/u/project/kp1/jlangerm/Projects/IGVF/R_Analysis/Analysis/2025-02-24_YS3v5_ChorusX_v11.5_pt1/IGVF.YS3v5.Merge.object.obj")
# 
# ###subsample for testing
# 
# your_rds <- your_rds[,1:10000]
# 
# 
# 
# Object_data <- your_rds
# cluster_search = "test" ## "fast" or "test"
# Exp.Name <- "YS3v7"
# output.dir.location <- "/u/project/kp1/jlangerm/Projects/IGVF/R_Analysis/Analysis/"
# 
# 
# 








# markers<-read.table("/u/project/kp1/jlangerm/Projects/IGVF/R_Analysis/Analysis/2025-05-07_YS3v6_atac_gene_link3/YS3v6.metadata.with.TrajScore.txt",
#                     sep="\t")


CreateGeneChorusReport <- function(Object_data = your_rds,
                                      output.dir.location = "/u/project/location",
                                      cluster_search = "test", ## "fast" or "test", test will try to guess the right amount of clusters wanted, fast will use a tested cell number scalar
                                      Exp.Name = "YS3v7",## Just a title for the data)
                                   noise_floor=0.1){   ##Minimum floor for genes to be paired into network. Filteration step, raw networks reported beore this is triggered 
                                        


####





##### Data preparation

Dir.Name <- paste0(Exp.Name,"_geneCHORUS_v0.12.5")
dir.create(paste0(output.dir.location ,Sys.Date(),"_", Dir.Name))
outdir <- c(paste0(output.dir.location,Sys.Date(),"_", Dir.Name,"/"))





##


##check for normalization
if(nrow(Object_data@assays$RNA$data)<1){
  Object_data <- NormalizeData(Object_data, normalization.method = "LogNormalize", scale.factor = 10000)
}
Object_data <- Object_data[rowSums(Object_data@assays$RNA$data)>0,]

## Create cell design file
markers <- as.data.frame(Object_data@meta.data)

##Check for umap coords for display
if(is.null(Object_data@reductions$pca)==T){
if( is.null(Object_data@reductions$umap)==T){
  # print("PCA and UMAP found")
# } else {
  print("UMAP and PCA not found, running UMAP derivation")
  Object_data = FindVariableFeatures(Object_data, verbose = F)
  all.genes <- rownames(Object_data)
  Object_data <- ScaleData(Object_data,  vars.to.regress = c("nFeature_RNA"))
  
  # IGVF_Merge <- ScaleData(IGVF_Merge, features = all.genes, vars.to.regress = c("nFeature_RNA"))
  
  
  
  Object_data <- RunPCA(Object_data, features = VariableFeatures(object = Object_data))
  Object_data <- FindNeighbors(Object_data, dims = 1:20)
  Object_data <- FindClusters(Object_data, resolution = 0.5)
  Object_data = RunUMAP(Object_data, dims = 1:20, verbose = T)
  
}
}

umap.meta <- as.data.frame(Object_data@reductions$umap@cell.embeddings)

names(umap.meta)<- c("X","Y")
markers$X <- umap.meta$X
markers$Y <- umap.meta$Y


 



pca.data <- as.data.frame((as.data.frame(Object_data@reductions$pca@cell.embeddings)))
# pca.data<- as.data.frame(t(pca.data))

###Figure out good center for cutoff
###This is the per cluster square sum operation

# sum(as.matrix(dist(df)^2)) / (2 * nrow(df))
### It turns out running a dist calc on thousands of cells in hundreds of clusters will take weeks
# c(1000,10000,25000,50000)

if(cluster_search=="fast"){

  approximate_cells_per_cluster <- sqrt(nrow(markers))/2
  
  ###Floor for noisey cell mixing
  if(approximate_cells_per_cluster<25){
    approximate_cells_per_cluster<-25
  }
  
  Center <- round(nrow(markers)/approximate_cells_per_cluster)
  
  markers$Cluster <- as.factor(kmeans(pca.data, iter.max=30, nstart=5, centers=Center )$cluster)
  
}



## Get a gestimate of K ranges I should work in 

if(cluster_search!="fast"){

centers.to.test.vector <- c(5,10,25)
centers.to.test.vector2 <- c(round(nrow(markers)/10000),round(nrow(markers)/1000),round(nrow(markers)/500), round(nrow(markers)/100), round(nrow(markers)/50), round(nrow(markers)/25), round(nrow(markers)/10)  )
centers.to.test.vector2<- centers.to.test.vector2[centers.to.test.vector2>25]

centers.to.test.vector<- c(centers.to.test.vector,centers.to.test.vector2)

Center<-5

# centers.to.test.vector <- c(5,10,25)

center.status.df <- NULL
for(Center in centers.to.test.vector){
  print(Center)
# selected.center.number<-round(nrow(markers)/num.of.cells)
  markers$Cluster <- as.factor(kmeans(pca.data, iter.max=30, nstart=5, centers= Center)$cluster)
  
  # wcss.df <- NULL
  # for(Cluster in unique(markers$Cluster)){
  #   print(Cluster)
  #   tmp.df <- Object_data@assays$RNA$data[, colnames(Object_data) %in% rownames(markers[markers$Cluster==Cluster,])]
  #   
  #   tmp.cluster.num<-sum(as.matrix(dist(tmp.df)^2)) / (2 * nrow(tmp.df))
  #   wcss.df <- c(wcss.df, tmp.cluster.num)
  #   }
  # 
  
  
  ##Going to have to use the mean correlation metric instead
  
  
  
  # cutoff.k.merge <- 0.99
  # markers$Cluster <- 0
  # markers$Cluster <- as.factor(kmeans(pca.data, iter.max=30, nstart=5, centers= selected.center.number))$cluster
  tag.table <- cbind.data.frame(rn=row.names(Object_data@assays$RNA$data) )
  # markers2 <- markers[markers$Timepoint==Day,]
  for( Sample in sort(unique(markers$Cluster))  ){
    if(length(row.names(markers[markers$Cluster %in% Sample,]) )>1){
      tag.temp <-cbind.data.frame(rn=row.names(Object_data),
                                  Sample=rowMeans(Object_data@assays$RNA$data[ ,colnames(Object_data) %in% row.names(markers[markers$Cluster %in% Sample,]) ]) )
      tag.table <- cbind.data.frame(tag.table, tag.temp$Sample)
      names(tag.table)[ncol(tag.table)] <- paste(Sample)} else{
        
        tag.temp <- Object_data@assays$RNA$data[ ,colnames(Object_data) %in% row.names(markers[markers$Cluster %in% Sample,]) ]
        tag.table$temp <- as.data.frame(tag.temp)
        colnames(tag.table)[ncol(tag.table)] <- Sample
      }
  }
  row.names(tag.table) <- tag.table$rn
  tag.table$rn <- NULL
  temp.corr <- cor(tag.table)
  diag(temp.corr) <- 0
  
  temp.corr2 <-temp.corr
  temp.corr2[temp.corr2==0]<-1
  
    tmp.center.status.df <- cbind.data.frame( Centers=Center,MeanCorr= mean(temp.corr),
                                              MaxMean= mean(apply(temp.corr,1,max)),
                                              MinMean= mean(apply(temp.corr2,1,min)))
                                              
                                              
    center.status.df<-rbind(center.status.df, tmp.center.status.df)

    
    ###Check to see if both maxmean and meancorr have decreased
    if(nrow(center.status.df)==1){next}
    if(center.status.df$MeanCorr[nrow(center.status.df)]< center.status.df$MeanCorr[nrow(center.status.df)-1]){
      print("Mean decreased")
      if(center.status.df$MaxMean[nrow(center.status.df)]< center.status.df$MaxMean[nrow(center.status.df)-1]){
        print("Max decreased, centers found at ")
        print(Center)
        break
      }
    }
    
}


}
# 
# ##Result 70k
# Centers  MeanCorr   MaxMean   MinMean
# 1       5 0.6319897 0.8960535 0.7229668
# 2      10 0.7195102 0.9190210 0.6892255
# 3      25 0.7678204 0.9457801 0.6424741
# 4      77 0.7652521 0.9680790 0.5529364
# 5     155 0.7654970 0.9712410 0.5245414
# 6     773 0.7424232 0.9661157 0.4826080
# 7    1546 0.7216635 0.9503943 0.3508444
# 8    3093 0.6863784 0.9220941 0.2639837
# 9    7732 0.6029551 0.8559628 0.1773481
# 
# #### Result 10k
# > center.status.df
# Centers  MeanCorr   MaxMean   MinMean
# 1       5 0.7303566 0.9459544 0.8728301
# 2      10 0.7828445 0.9434230 0.7831126
# 3      25 0.8157722 0.9606455 0.7132621
# 4     100 0.8028425 0.9650941 0.6184452
# 5     200 0.7903634 0.9562368 0.5742333
# 6     400 0.7664313 0.9378375 0.5054699
# 7    1000 0.6953511 0.8862752 0.3276466



### What does this mean??
#####One interpretation of this is that the return tail is when we start to see dropout effect
##Max 155 is ~155 clusters of 450 cells each
### Largest amount is 8k clusters with 10 cells each
## Previously I ran this with 100 cells clusters I think

markers$Cluster <- paste0("C",markers$Cluster)

length(markers$Cluster[!duplicated(markers$Cluster)])
markers$meanX <- 0 ; markers$meanY <- 0
for(Cluster in unique(markers$Cluster)){
  markers[markers$Cluster==Cluster,]$meanX <- mean(markers[markers$Cluster==Cluster,]$X)
  markers[markers$Cluster==Cluster,]$meanY <- mean(markers[markers$Cluster==Cluster,]$Y)
}


###This can be exceded by the pallete collection

palette.vector <- distinctColorPalette(length(unique(markers$Cluster)))
# 
plot1<-ggplot(markers)+
  geom_point(aes(X,Y, color=(Cluster)))+
  geom_text(aes(meanX,meanY, label=Cluster))+
  # scale_colour_discrete(type =palette.vector)+
  labs( x="UMAP 1", y="UMAP 2", color="Clusters")+
  theme_classic()+
  theme(legend.position = "")
ggsave(plot=plot1,filename=paste0(outdir,  Exp.Name, ".OverKClusters.labeled_nice.UMAP.png"),
       height=10, width=15)

write.table(markers, paste0(outdir,Exp.Name,".processed_and_Clustered_for_part2.metadata.txt"),
            sep="\t", quote=F)




####
## Start VOX workflow

print("Starting geneChorus network discovery")

##Does this actually require this? I don't think so. looks like I run the correlation after simplifying the data
## Call it a lot though. might be a good edit to stop as.data.frame bottlenecks but Im not sure how this behaves on massive data anyway
normalized.data <-as.data.frame(Object_data@assays$RNA$data)




###Have to rerun tag table now with new cluster names
tag.table <- NULL
tag.table$rn <- row.names(normalized.data)
for( Sample1 in sort(unique(markers$Cluster))  ){
  print(Sample1)
  
  if(length(row.names(markers[markers$Cluster %in% Sample1,]) )>1){
    
    tag.temp <-cbind.data.frame(rn=row.names(normalized.data),
                                Sample1=rowMeans(normalized.data[ ,colnames(normalized.data) %in% row.names(markers[markers$Cluster %in% Sample1,]) ]) )
    tag.table <- cbind.data.frame(tag.table, tag.temp$Sample1)
    names(tag.table)[ncol(tag.table)] <- paste(Sample1)} else{
      
      tag.temp <- normalized.data[ ,colnames(normalized.data) %in% row.names(markers[markers$Cluster %in% Sample1,]) ]
      tag.table$temp <- tag.temp
      names(tag.table)[ncol(tag.table)] <- paste(Sample1)
    }
}
row.names(tag.table) <- tag.table$rn
tag.table$rn <- NULL

new.sample.df<-as.data.frame(t(t(tag.table)/colSums(tag.table)))
new.sample.df$rn <- NULL


###Size factoring high minimums to prepare for correlation

new.sample.df2 <-new.sample.df
new.sample.df2[new.sample.df2==0]<-1000
checker.min <- cbind.data.frame(Cell=colnames(new.sample.df2), ColMins=apply(new.sample.df2,2,min))

##size factor of new.sample.df2
checker.min$SizeFactor <- checker.min$ColMins/mean(checker.min$ColMins)
checker.min$AdjustFactor <- checker.min$SizeFactor
checker.min[checker.min$SizeFactor<=1,]$AdjustFactor <- 1
sized.sample.df <-as.data.frame(t( t(new.sample.df)/checker.min$AdjustFactor))
new.sample.df<-sized.sample.df



##########
##Builds variance tag
gene.info <- cbind.data.frame(Gene=row.names(new.sample.df), Mean=rowMeans(new.sample.df))
gene.info$Stdev <- apply(new.sample.df,1, FUN=sd)
gene.info$Variance <- gene.info$Stdev/gene.info$Mean
gene.info$Bin.Mean <- cut(gene.info$Mean,5, labels=F)

gene.info2 <- gene.info
tissue.names <-names(new.sample.df)
for(Tissue in tissue.names){
  gene.info2$input <- new.sample.df[,Tissue]
  # gene.info2$input <- new.sample.df[!rownames(new.sample.df) %in% rownames(noise.store),Tissue]
  names(gene.info2)[ncol(gene.info2)] <- Tissue
}


scaled.tag<-(gene.info2[,c(6:ncol(gene.info2))]-gene.info2$Mean)/gene.info2$Stdev
scaled.tag <- scaled.tag[!row.names(scaled.tag) %in% (row.names(scaled.tag[is.na(scaled.tag),])) ,]


scaled.bin.tag <- scaled.tag
scaled.bin.tag<- scaled.bin.tag*scaled.bin.tag*(scaled.bin.tag/abs(scaled.bin.tag))








##Tests to observe scaled gene data
# test.df <-scaled.tag[row.names(scaled.tag) %in% c("MUC5B","BPIFB1", "XBP1", "LTF","CFB",  "FCGBP","LPO" , "USP32","RTN2","SMCO4","SOD2"),]
# cor(t (test.df))
# test.df <-scaled.bin.tag[row.names(scaled.bin.tag) %in% c("MUC5B","BPIFB1", "XBP1", "LTF","CFB",  "FCGBP","LPO" , "USP32","RTN2","SMCO4","SOD2"),]
# cor(t (test.df))


###Round 1 start

print("Starting A Pass")

variance.tag <-gene.info

gene.cor <- cor( t(scaled.bin.tag) )
diag(gene.cor)<-0

tester.melt <- melt(gene.cor)


tester.melt<- tester.melt[tester.melt$value>0.1,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]

# ###Choose an initial cutoff
# plot1<-ggplot(tester.melt)+
#   theme_classic()+
#   geom_density(aes(value))+
#   geom_vline(xintercept = c(  mean(tester.melt$value)+sd(tester.melt$value) ))
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Vox.firstround.intercorr2_all.cutoff.distrubtion.png"),
#        height=3, width=5)

tester.melt<- tester.melt[tester.melt$value>0.7,]


variance.tag$CG1<- paste0("A", 1:nrow(variance.tag))
variance.tag$CG2 <- variance.tag$CG1
# Transcode="RCN3"

for(Transcode in tester.melt$Var1 ){
  # print(paste0("Converting - ",Transcode))
  change.target  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
  change.target2 <- variance.tag[variance.tag$Gene==change.target,]$CG1
  # print(paste0("changing to " , change.target2))
  variance.tag[variance.tag$Gene==Transcode,]$CG2 <- as.character(change.target2)

}


variance.tag$Transcode <- variance.tag$CG2
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag$CG2)){
  # print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    # print("next")
    }
  dictionary.spread<-Transcode
  tester.df <- variance.tag[,c("CG2","CG1")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$CG2 %in% dictionary.spread,]$CG2)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$CG1 %in% dictionary.spread,]$CG2)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$CG2 %in% dictionary.spread,]$CG1)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$CG1 %in% dictionary.spread,]$CG1)
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>5){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag[variance.tag$CG2 %in% dictionary.spread,]$Transcode <- paste0("AA", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag$CG2<- variance.tag$Transcode
# print(paste0("Now have ", length(unique(variance.tag$CG2)), " transcodes compared to ", length(unique(variance.tag$CG1))) )




corr.round.x.df <- cbind.data.frame(Corr.Cutoff=c(0.95,0.92,0.89,0.85),
                                    GrabColumn=c("CG2","CG3","CG4","CG5"),
                                    NewColumn=c("CG3","CG4","CG5","CG6"),
                                    GroupName=c("AB","AC","AD","AE") )

for(Corr.Cutoff in corr.round.x.df$Corr.Cutoff){
  # print(Corr.Cutoff)
  
  GrabColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GrabColumn
  NewColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$NewColumn
  GroupName <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GroupName
  
  trans.means <- cbind.data.frame(rn= names(new.sample.df))
  ##Recalculate transcode means 
  for(Transcode in unique(variance.tag[,GrabColumn]) ){
    top.genes<- variance.tag[ variance.tag[,GrabColumn] %in%  Transcode,]$Gene
    trans.means$Gene <- colMeans(new.sample.df[top.genes,])
    names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  }
  row.names(trans.means) <- trans.means$rn
  trans.means$rn <- NULL
  exp.level.tester <- as.data.frame(t(trans.means))
  exp.level.tester<-as.data.frame(t(t(exp.level.tester)/colSums(exp.level.tester)))
  
  group.info <- cbind.data.frame(Gene=row.names(exp.level.tester), Mean=rowMeans(exp.level.tester))
  group.info$Stdev <- apply(exp.level.tester,1, FUN=sd)
  group.scale<-(exp.level.tester-group.info$Mean)/group.info$Stdev
  # scaled.tag <- scaled.tag[!row.names(scaled.tag) %in% (row.names(scaled.tag[is.na(scaled.tag),])) ,]
  
  #########squared scale tag
  group.scaled.bin.tag <- group.scale
  group.scaled.bin.tag<- group.scaled.bin.tag*group.scaled.bin.tag*(group.scaled.bin.tag/abs(group.scaled.bin.tag))
  
  ###Just using exp.level.tester here
  gene.cor <- cor( t(group.scaled.bin.tag) )
  diag(gene.cor)<-0
  

  
  tester.melt <- melt(gene.cor)
  tester.melt<- tester.melt[tester.melt$value>0.5,]
  tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]

  tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]

  variance.tag$NewColumn <- variance.tag[,GrabColumn]
  names(variance.tag)[ncol(variance.tag)] <- NewColumn

  tester.melt<-tester.melt[tester.melt$value>Corr.Cutoff,]
  

  for(Transcode in tester.melt$Var1 ){
    change.target2  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
    variance.tag[variance.tag[,GrabColumn] %in% Transcode,NewColumn] <- as.character(change.target2)
  }
  
  variance.tag$Transcode <- variance.tag[,NewColumn]
  dictionary.summary <- NULL
  tb.counter <-1
  for(Transcode in unique(variance.tag[,NewColumn])){
    if(Transcode %in% dictionary.summary==T){next
    }
    dictionary.spread<-Transcode
    tester.df <- variance.tag[,c(GrabColumn, NewColumn)]
    n=1
    repeat{
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,NewColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,NewColumn])
      
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread<-unique(dictionary.spread)
      n=n+1
      if(n>3){break}
    }
    dictionary.spread<-unique(dictionary.spread)
    variance.tag[variance.tag[,NewColumn] %in% dictionary.spread,]$Transcode <- paste0(GroupName, tb.counter)
    tb.counter <- tb.counter+1
    dictionary.summary<- c(dictionary.summary, dictionary.spread)
  }
  variance.tag[,NewColumn]<- variance.tag$Transcode
  # print(paste0("Now have ", length(unique(variance.tag[,NewColumn])), " transcodes compared to ", length(unique(variance.tag[,GrabColumn]))) )
}




##############
### Flag Core patterns 

transcode.table <- as.data.frame(table(variance.tag$CG6))
transcode.table2 <- transcode.table[transcode.table$Freq>9,]
transcode.table2$GroupMean <- 0
for(Transcode in transcode.table2$Var1){
  
  genes.of.interest <- variance.tag[variance.tag$CG6==Transcode,]$Gene
  transcode.table2[transcode.table2$Var1==Transcode,]$GroupMean <- mean(variance.tag[variance.tag$Gene %in% genes.of.interest,]$Mean)
  
}

transcode.table2 <- transcode.table2[transcode.table2$GroupMean > mean(variance.tag$Mean)/10,]
variance.tag$FirstPass <- "B1"
variance.tag[variance.tag$CG6 %in% transcode.table2$Var1,]$FirstPass <- variance.tag[variance.tag$CG6 %in% transcode.table2$Var1,]$CG6
# write.table(variance.tag, paste0(outdir, Exp.Name, ".vox_corr_table.step1.gene.table.txt"),
#             sep="\t")










############################
### Level B1 merges

print("Starting B Pass")
variance.tag2 <- variance.tag[variance.tag$FirstPass=="B1",]


############################
### Level B1 merges


variance.tag2 <- variance.tag[variance.tag$FirstPass=="B1",]


# scaled.tag

p2.scaled.bin.tag <- scaled.tag[row.names(scaled.tag) %in% variance.tag2$Gene,]
p2.scaled.bin.tag<- p2.scaled.bin.tag*p2.scaled.bin.tag*(p2.scaled.bin.tag/abs(p2.scaled.bin.tag))
gene.cor <- cor( t(p2.scaled.bin.tag) )
diag(gene.cor)<-0

tester.melt <- melt(gene.cor)
tester.melt<- tester.melt[tester.melt$value>0.1,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
# plot1<-ggplot(tester.melt)+
#   theme_classic()+
#   geom_density(aes(value))+
#   geom_vline(xintercept = c(  mean(tester.melt$value)+sd(tester.melt$value) ), color="red")+
#   geom_vline(xintercept=(0.6))
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Vox.secondpass.intercorr_all.cutoff.distrubtion.png"),
#        height=3, width=5)

tester.melt<- tester.melt[tester.melt$value>0.6,]


variance.tag2$P2CG1 <- variance.tag2$Gene

for(Transcode in tester.melt$Var1 ){
  # print(paste0("Converting - ",Transcode))

  change.target  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
  change.target2 <- variance.tag2[variance.tag2$Gene==change.target,]$Gene
  
  # print(paste0("changing to " , change.target2))
  variance.tag2[variance.tag2$Gene==Transcode,]$P2CG1 <- as.character(change.target2)
}


variance.tag2$Transcode <- variance.tag2$P2CG1
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$P2CG1)){
  # print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    # print("next")
    }
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("P2CG1","Gene")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P2CG1 %in% dictionary.spread,]$P2CG1)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$P2CG1)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P2CG1 %in% dictionary.spread,]$Gene)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$Gene)
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$P2CG1 %in% dictionary.spread,]$Transcode <- paste0("BA", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$P2CG1<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$P2CG1)), " transcodes compared to ", length(unique(variance.tag2$Gene))) )






#######Compression pass


GrabColumn="P2CG1"
NewColumn <- "P2CG2"
GroupName <-"BB"


corr.round.x.df <- cbind.data.frame(Corr.Cutoff=c(0.85,0.8,0.75,0.7),
                                    GrabColumn=c("P2CG1","P2CG2","P2CG3","P2CG4"),
                                    NewColumn=c("P2CG2","P2CG3","P2CG4","P2CG5"),
                                    GroupName=c("BB","BC","BD","BE") )



for(Corr.Cutoff in corr.round.x.df$Corr.Cutoff){
  # print(Corr.Cutoff)
  
  GrabColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GrabColumn
  NewColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$NewColumn
  GroupName <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GroupName
  
  
  trans.means <- cbind.data.frame(rn= names(new.sample.df))
  ##Recalculate transcode means 
  for(Transcode in unique(variance.tag2[,GrabColumn]) ){
    
    # print(Transcode)
    top.genes<- variance.tag2[ variance.tag2[,GrabColumn] %in%  Transcode,]$Gene
    trans.means$Gene <- colMeans(new.sample.df[top.genes,])
    names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  }
  row.names(trans.means) <- trans.means$rn
  trans.means$rn <- NULL
  exp.level.tester <- as.data.frame(t(trans.means))
  exp.level.tester<-as.data.frame(t(t(exp.level.tester)/colSums(exp.level.tester)))
  
  
  
  group.info <- cbind.data.frame(Gene=row.names(exp.level.tester), Mean=rowMeans(exp.level.tester))
  group.info$Stdev <- apply(exp.level.tester,1, FUN=sd)
  # gene.info$Variance <- gene.info$Stdev/gene.info$Mean
  # gene.info$Bin.Mean <- cut(gene.info$Mean,5, labels=F)
  
  group.scale<-(exp.level.tester-group.info$Mean)/group.info$Stdev
  # scaled.tag <- scaled.tag[!row.names(scaled.tag) %in% (row.names(scaled.tag[is.na(scaled.tag),])) ,]
  
  #########squared scale tag
  group.scaled.bin.tag <- group.scale
  group.scaled.bin.tag<- group.scaled.bin.tag*group.scaled.bin.tag*(group.scaled.bin.tag/abs(group.scaled.bin.tag))
  
  
  ###Just using exp.level.tester here
  gene.cor <- cor( t(group.scaled.bin.tag) )
  diag(gene.cor)<-0
  
  
  

  tester.melt <- melt(gene.cor)
  tester.melt<- tester.melt[tester.melt$value>0.5,]
  tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
  
  tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]

  
  tester.melt<-tester.melt[tester.melt$value>Corr.Cutoff,]
  
  
  
  variance.tag2$NewColumn <- variance.tag2[,GrabColumn]
  names(variance.tag2)[ncol(variance.tag2)] <- NewColumn
  

  for(Transcode in tester.melt$Var1 ){
    change.target2  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
    variance.tag2[variance.tag2[,GrabColumn] %in% Transcode,NewColumn] <- as.character(change.target2)
  }
  
  
  
  
  variance.tag2$Transcode <- variance.tag2[,NewColumn]
  dictionary.summary <- NULL
  tb.counter <-1
  for(Transcode in unique(variance.tag2[,NewColumn])){
    # print(Transcode)
    if(Transcode %in% dictionary.summary==T){next
      # print("next")
    }
    dictionary.spread<-Transcode
    tester.df <- variance.tag2[,c(GrabColumn, NewColumn)]
    n=1
    repeat{
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,NewColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,NewColumn])
      
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread<-unique(dictionary.spread)
      n=n+1
      if(n>3){break}
    }
    dictionary.spread<-unique(dictionary.spread)
    variance.tag2[variance.tag2[,NewColumn] %in% dictionary.spread,]$Transcode <- paste0(GroupName, tb.counter)
    tb.counter <- tb.counter+1
    dictionary.summary<- c(dictionary.summary, dictionary.spread)
  }
  
  
  variance.tag2[,NewColumn]<- variance.tag2$Transcode
  # print(paste0("Now have ", length(unique(variance.tag2[,NewColumn])), " transcodes compared to ", length(unique(variance.tag2[,GrabColumn]))) )
  
  
  
}


# 
# write.table(variance.tag2, paste0(outdir, Exp.Name, ".vox_corr_table.step2.gene.table.txt"),
#             sep="\t")
# 







#######################################
### Round 3!!!!

##############
### Flag Core patterns 


transcode.table <- as.data.frame(table(variance.tag2$P2CG5))
transcode.table2 <- transcode.table[transcode.table$Freq>9,]
transcode.table2$GroupMean <- 0
for(Transcode in transcode.table2$Var1){
  genes.of.interest <- variance.tag2[variance.tag2$P2CG5==Transcode,]$Gene
  transcode.table2[transcode.table2$Var1==Transcode,]$GroupMean <- mean(variance.tag2[variance.tag2$Gene %in% genes.of.interest,]$Mean)
  
}
transcode.table2 <- transcode.table2[transcode.table2$GroupMean > mean(variance.tag$Mean)/10,]

second.pass.genes <- variance.tag2[variance.tag2$P2CG5 %in% transcode.table2$Var1,]$Gene

variance.tag$SecondPass <- "C1"
variance.tag[variance.tag$FirstPass!="B1",]$SecondPass <- variance.tag[variance.tag$FirstPass!="B1",]$FirstPass
variance.tag[variance.tag$Gene %in% second.pass.genes,]$SecondPass <- variance.tag2[variance.tag2$P2CG5 %in% transcode.table2$Var1,]$P2CG5

# 
# write.table(variance.tag, paste0(outdir, Exp.Name, ".vox_corr_table.step2_all.gene.table.txt"),
#             sep="\t")
# 




############################
### Level C1 merges
print("Starting C Pass")

variance.tag2 <- variance.tag[variance.tag$SecondPass=="C1",]


# scaled.tag

p2.scaled.bin.tag <- scaled.tag[row.names(scaled.tag) %in% variance.tag2$Gene,]
p2.scaled.bin.tag<- p2.scaled.bin.tag*p2.scaled.bin.tag*(p2.scaled.bin.tag/abs(p2.scaled.bin.tag))
gene.cor <- cor( t(p2.scaled.bin.tag) )
diag(gene.cor)<-0

tester.melt <- melt(gene.cor)
tester.melt<- tester.melt[tester.melt$value>0.1,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
# plot1<-ggplot(tester.melt)+
#   theme_classic()+
#   geom_density(aes(value))+
#   geom_vline(xintercept = c(  mean(tester.melt$value)+sd(tester.melt$value) ), color="red")+
#   geom_vline(xintercept=(0.55))
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Vox.thirdpass.intercorr_all.cutoff.distrubtion.png"),
#        height=3, width=5)

tester.melt<- tester.melt[tester.melt$value>0.55,]


variance.tag2$P3CG1 <- variance.tag2$Gene


for(Transcode in tester.melt$Var1 ){
  # print(paste0("Converting - ",Transcode))
  change.target  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
  change.target2 <- variance.tag2[variance.tag2$Gene==change.target,]$Gene
  
  # print(paste0("changing to " , change.target2))
  variance.tag2[variance.tag2$Gene==Transcode,]$P3CG1 <- as.character(change.target2)
}


variance.tag2$Transcode <- variance.tag2$P3CG1
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$P3CG1)){
  # print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    # print("next")
    }
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("P3CG1","Gene")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P3CG1 %in% dictionary.spread,]$P3CG1)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$P3CG1)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P3CG1 %in% dictionary.spread,]$Gene)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$Gene)
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$P3CG1 %in% dictionary.spread,]$Transcode <- paste0("CA", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$P3CG1<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$P3CG1)), " transcodes compared to ", length(unique(variance.tag2$Gene))) )






#######Compression pass


GrabColumn="P2CG1"
NewColumn <- "P2CG2"
GroupName <-"BB"


corr.round.x.df <- cbind.data.frame(Corr.Cutoff=c(0.75,0.7,0.65,0.6),
                                    GrabColumn=c("P3CG1","P3CG2","P3CG3","P3CG4"),
                                    NewColumn=c("P3CG2","P3CG3","P3CG4","P3CG5"),
                                    GroupName=c("CB","CC","CD","CE") )



for(Corr.Cutoff in corr.round.x.df$Corr.Cutoff){
  # print(Corr.Cutoff)
  
  GrabColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GrabColumn
  NewColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$NewColumn
  GroupName <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GroupName
  
  
  trans.means <- cbind.data.frame(rn= names(new.sample.df))
  ##Recalculate transcode means 
  for(Transcode in unique(variance.tag2[,GrabColumn]) ){
    
    # print(Transcode)
    top.genes<- variance.tag2[ variance.tag2[,GrabColumn] %in%  Transcode,]$Gene
    trans.means$Gene <- colMeans(new.sample.df[top.genes,])
    names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  }
  row.names(trans.means) <- trans.means$rn
  trans.means$rn <- NULL
  exp.level.tester <- as.data.frame(t(trans.means))
  exp.level.tester<-as.data.frame(t(t(exp.level.tester)/colSums(exp.level.tester)))
  
  
  
  group.info <- cbind.data.frame(Gene=row.names(exp.level.tester), Mean=rowMeans(exp.level.tester))
  group.info$Stdev <- apply(exp.level.tester,1, FUN=sd)
  # gene.info$Variance <- gene.info$Stdev/gene.info$Mean
  # gene.info$Bin.Mean <- cut(gene.info$Mean,5, labels=F)
  
  group.scale<-(exp.level.tester-group.info$Mean)/group.info$Stdev
  # scaled.tag <- scaled.tag[!row.names(scaled.tag) %in% (row.names(scaled.tag[is.na(scaled.tag),])) ,]
  
  #########squared scale tag
  group.scaled.bin.tag <- group.scale
  group.scaled.bin.tag<- group.scaled.bin.tag*group.scaled.bin.tag*(group.scaled.bin.tag/abs(group.scaled.bin.tag))
  
  
  ###Just using exp.level.tester here
  gene.cor <- cor( t(group.scaled.bin.tag) )
  diag(gene.cor)<-0
  
  
  

  tester.melt <- melt(gene.cor)
  tester.melt<- tester.melt[tester.melt$value>0.5,]
  tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
  tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
  
  tester.melt<-tester.melt[tester.melt$value>Corr.Cutoff,]
  
  variance.tag2$NewColumn <- variance.tag2[,GrabColumn]
  names(variance.tag2)[ncol(variance.tag2)] <- NewColumn
  

  for(Transcode in tester.melt$Var1 ){
    change.target2  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
    variance.tag2[variance.tag2[,GrabColumn] %in% Transcode,NewColumn] <- as.character(change.target2)
  }
  
  
  
  
  variance.tag2$Transcode <- variance.tag2[,NewColumn]
  dictionary.summary <- NULL
  tb.counter <-1
  for(Transcode in unique(variance.tag2[,NewColumn])){
    # print(Transcode)
    if(Transcode %in% dictionary.summary==T){next
      # print("next")
    }
    dictionary.spread<-Transcode
    tester.df <- variance.tag2[,c(GrabColumn, NewColumn)]
    n=1
    repeat{
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,NewColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,NewColumn])
      
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread<-unique(dictionary.spread)
      n=n+1
      if(n>3){break}
    }
    dictionary.spread<-unique(dictionary.spread)
    variance.tag2[variance.tag2[,NewColumn] %in% dictionary.spread,]$Transcode <- paste0(GroupName, tb.counter)
    tb.counter <- tb.counter+1
    dictionary.summary<- c(dictionary.summary, dictionary.spread)
  }
  
  
  variance.tag2[,NewColumn]<- variance.tag2$Transcode
  # print(paste0("Now have ", length(unique(variance.tag2[,NewColumn])), " transcodes compared to ", length(unique(variance.tag2[,GrabColumn]))) )
  
  
  
}


# 
# write.table(variance.tag2, paste0(outdir, Exp.Name, ".vox_corr_table.step3.gene.table.txt"),
#             sep="\t")
# 



##############
### Flag Core patterns 4


transcode.table <- as.data.frame(table(variance.tag2$P3CG5))
transcode.table2 <- transcode.table[transcode.table$Freq>9,]
transcode.table2$GroupMean <- 0
for(Transcode in transcode.table2$Var1){
  genes.of.interest <- variance.tag2[variance.tag2$P3CG5==Transcode,]$Gene
  transcode.table2[transcode.table2$Var1==Transcode,]$GroupMean <- mean(variance.tag2[variance.tag2$Gene %in% genes.of.interest,]$Mean)
  
}
transcode.table2 <- transcode.table2[transcode.table2$GroupMean > mean(variance.tag$Mean)/10,]

second.pass.genes <- variance.tag2[variance.tag2$P3CG5 %in% transcode.table2$Var1,]$Gene

variance.tag$ThirdPass <- "D1"

variance.tag[variance.tag$SecondPass!="C1",]$ThirdPass <- variance.tag[variance.tag$SecondPass!="C1",]$SecondPass
variance.tag[variance.tag$Gene %in% second.pass.genes,]$ThirdPass <- variance.tag2[variance.tag2$P3CG5 %in% transcode.table2$Var1,]$P3CG5


# write.table(variance.tag, paste0(outdir, Exp.Name, ".vox_corr_table.step3_all.gene.table.txt"),
#             sep="\t")
# 


############################
### Level D1 merges
print("Starting D Pass")

variance.tag2 <- variance.tag[variance.tag$ThirdPass=="D1",]


# scaled.tag

p2.scaled.bin.tag <- scaled.tag[row.names(scaled.tag) %in% variance.tag2$Gene,]
p2.scaled.bin.tag<- p2.scaled.bin.tag*p2.scaled.bin.tag*(p2.scaled.bin.tag/abs(p2.scaled.bin.tag))
gene.cor <- cor( t(p2.scaled.bin.tag) )
diag(gene.cor)<-0

tester.melt <- melt(gene.cor)
tester.melt<- tester.melt[tester.melt$value>0.1,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
# plot1<-ggplot(tester.melt)+
#   theme_classic()+
#   geom_density(aes(value))+
#   geom_vline(xintercept = c(  mean(tester.melt$value)+sd(tester.melt$value) ), color="red")+
#   geom_vline(xintercept=(0.55))
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Vox.fourthpass.intercorr_all.cutoff.distrubtion.png"),
#        height=3, width=5)

tester.melt<- tester.melt[tester.melt$value>0.50,]


variance.tag2$P4CG1 <- variance.tag2$Gene


for(Transcode in tester.melt$Var1 ){
  # print(paste0("Converting - ",Transcode))

  change.target  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
  change.target2 <- variance.tag2[variance.tag2$Gene==change.target,]$Gene
  
  # print(paste0("changing to " , change.target2))
  variance.tag2[variance.tag2$Gene==Transcode,]$P4CG1 <- as.character(change.target2)
}


variance.tag2$Transcode <- variance.tag2$P4CG1
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$P4CG1)){
  # print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    # print("next")
    }
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("P4CG1","Gene")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P4CG1 %in% dictionary.spread,]$P4CG1)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$P4CG1)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P4CG1 %in% dictionary.spread,]$Gene)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$Gene)
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$P4CG1 %in% dictionary.spread,]$Transcode <- paste0("DA", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$P4CG1<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$P4CG1)), " transcodes compared to ", length(unique(variance.tag2$Gene))) )






#######Compression pass


GrabColumn="P2CG1"
NewColumn <- "P2CG2"
GroupName <-"BB"


corr.round.x.df <- cbind.data.frame(Corr.Cutoff=c(0.70,0.65,0.60,0.55),
                                    GrabColumn=c("P4CG1","P4CG2","P4CG3","P4CG4"),
                                    NewColumn=c("P4CG2","P4CG3","P4CG4","P4CG5"),
                                    GroupName=c("DB","DC","DD","DE") )



for(Corr.Cutoff in corr.round.x.df$Corr.Cutoff){
  # print(Corr.Cutoff)
  
  GrabColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GrabColumn
  NewColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$NewColumn
  GroupName <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GroupName
  
  
  trans.means <- cbind.data.frame(rn= names(new.sample.df))
  ##Recalculate transcode means 
  for(Transcode in unique(variance.tag2[,GrabColumn]) ){
    
    # print(Transcode)
    top.genes<- variance.tag2[ variance.tag2[,GrabColumn] %in%  Transcode,]$Gene
    trans.means$Gene <- colMeans(new.sample.df[top.genes,])
    names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  }
  row.names(trans.means) <- trans.means$rn
  trans.means$rn <- NULL
  exp.level.tester <- as.data.frame(t(trans.means))
  exp.level.tester<-as.data.frame(t(t(exp.level.tester)/colSums(exp.level.tester)))
  
  
  
  group.info <- cbind.data.frame(Gene=row.names(exp.level.tester), Mean=rowMeans(exp.level.tester))
  group.info$Stdev <- apply(exp.level.tester,1, FUN=sd)
  # gene.info$Variance <- gene.info$Stdev/gene.info$Mean
  # gene.info$Bin.Mean <- cut(gene.info$Mean,5, labels=F)
  
  group.scale<-(exp.level.tester-group.info$Mean)/group.info$Stdev
  # scaled.tag <- scaled.tag[!row.names(scaled.tag) %in% (row.names(scaled.tag[is.na(scaled.tag),])) ,]
  
  #########squared scale tag
  group.scaled.bin.tag <- group.scale
  group.scaled.bin.tag<- group.scaled.bin.tag*group.scaled.bin.tag*(group.scaled.bin.tag/abs(group.scaled.bin.tag))
  
  
  ###Just using exp.level.tester here
  gene.cor <- cor( t(group.scaled.bin.tag) )
  diag(gene.cor)<-0
  

  
  tester.melt <- melt(gene.cor)
  tester.melt<- tester.melt[tester.melt$value>0.5,]
  tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
  
  tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
  
  tester.melt<-tester.melt[tester.melt$value>Corr.Cutoff,]
  
  variance.tag2$NewColumn <- variance.tag2[,GrabColumn]
  names(variance.tag2)[ncol(variance.tag2)] <- NewColumn

  for(Transcode in tester.melt$Var1 ){

    change.target2  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
    variance.tag2[variance.tag2[,GrabColumn] %in% Transcode,NewColumn] <- as.character(change.target2)
  }
  
  
  
  
  variance.tag2$Transcode <- variance.tag2[,NewColumn]
  dictionary.summary <- NULL
  tb.counter <-1
  for(Transcode in unique(variance.tag2[,NewColumn])){
    # print(Transcode)
    if(Transcode %in% dictionary.summary==T){next
      # print("next")
    }
    dictionary.spread<-Transcode
    tester.df <- variance.tag2[,c(GrabColumn, NewColumn)]
    n=1
    repeat{
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,NewColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,NewColumn])
      
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread<-unique(dictionary.spread)
      n=n+1
      if(n>3){break}
    }
    dictionary.spread<-unique(dictionary.spread)
    variance.tag2[variance.tag2[,NewColumn] %in% dictionary.spread,]$Transcode <- paste0(GroupName, tb.counter)
    tb.counter <- tb.counter+1
    dictionary.summary<- c(dictionary.summary, dictionary.spread)
  }
  
  
  variance.tag2[,NewColumn]<- variance.tag2$Transcode
  # print(paste0("Now have ", length(unique(variance.tag2[,NewColumn])), " transcodes compared to ", length(unique(variance.tag2[,GrabColumn]))) )
  
  
  
}



# write.table(variance.tag2, paste0(outdir, Exp.Name, ".vox_corr_table.step4.gene.table.txt"),
#             sep="\t")




##############
### Flag Core patterns 5


transcode.table <- as.data.frame(table(variance.tag2$P4CG5))
transcode.table2 <- transcode.table[transcode.table$Freq>9,]
transcode.table2$GroupMean <- 0
for(Transcode in transcode.table2$Var1){
  genes.of.interest <- variance.tag2[variance.tag2$P4CG5==Transcode,]$Gene
  transcode.table2[transcode.table2$Var1==Transcode,]$GroupMean <- mean(variance.tag2[variance.tag2$Gene %in% genes.of.interest,]$Mean)
  
}
transcode.table2 <- transcode.table2[transcode.table2$GroupMean > mean(variance.tag$Mean)/10,]

second.pass.genes <- variance.tag2[variance.tag2$P4CG5 %in% transcode.table2$Var1,]$Gene

variance.tag$FourthPass <- "E1"

variance.tag[variance.tag$ThirdPass!="D1",]$FourthPass <- variance.tag[variance.tag$ThirdPass!="D1",]$ThirdPass
variance.tag[variance.tag$Gene %in% second.pass.genes,]$FourthPass <- variance.tag2[variance.tag2$P4CG5 %in% transcode.table2$Var1,]$P4CG5


# write.table(variance.tag, paste0(outdir, Exp.Name, ".vox_corr_table.step3_all.gene.table.txt"),
#             sep="\t")

############################
### Level E1 merges

variance.tag2 <- variance.tag[variance.tag$FourthPass=="E1",]

p2.scaled.bin.tag <- scaled.tag[row.names(scaled.tag) %in% variance.tag2$Gene,]
p2.scaled.bin.tag<- p2.scaled.bin.tag*p2.scaled.bin.tag*(p2.scaled.bin.tag/abs(p2.scaled.bin.tag))
gene.cor <- cor( t(p2.scaled.bin.tag) )
diag(gene.cor)<-0

tester.melt <- melt(gene.cor)
tester.melt<- tester.melt[tester.melt$value>0.1,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
# plot1<-ggplot(tester.melt)+
#   theme_classic()+
#   geom_density(aes(value))+
#   geom_vline(xintercept = c(  mean(tester.melt$value)+sd(tester.melt$value) ), color="red")+
#   geom_vline(xintercept=(0.45))
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Vox.fifthpass.intercorr_all.cutoff.distrubtion.png"),
#        height=3, width=5)

##change
tester.melt<- tester.melt[tester.melt$value>0.45,]
##change
variance.tag2$P5CG1 <- variance.tag2$Gene


for(Transcode in tester.melt$Var1 ){
  # print(paste0("Converting - ",Transcode))
  change.target  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
  change.target2 <- variance.tag2[variance.tag2$Gene==change.target,]$Gene
  # print(paste0("changing to " , change.target2))
  variance.tag2[variance.tag2$Gene==Transcode,]$P5CG1 <- as.character(change.target2)
}


variance.tag2$Transcode <- variance.tag2$P5CG1
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$P5CG1)){
  # print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    # print("next")
    }
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("P5CG1","Gene")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P5CG1 %in% dictionary.spread,]$P5CG1)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$P5CG1)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P5CG1 %in% dictionary.spread,]$Gene)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$Gene)
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$P5CG1 %in% dictionary.spread,]$Transcode <- paste0("EA", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$P5CG1<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$P5CG1)), " transcodes compared to ", length(unique(variance.tag2$Gene))) )






#######Compression pass


GrabColumn="P2CG1"
NewColumn <- "P2CG2"
GroupName <-"BB"


corr.round.x.df <- cbind.data.frame(Corr.Cutoff=c(0.65,0.60,0.55,0.50),
                                    GrabColumn=c("P5CG1","P5CG2","P5CG3","P5CG4"),
                                    NewColumn=c("P5CG2","P5CG3","P5CG4","P5CG5"),
                                    GroupName=c("EB","EC","ED","EE") )



for(Corr.Cutoff in corr.round.x.df$Corr.Cutoff){
  # print(Corr.Cutoff)
  
  GrabColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GrabColumn
  NewColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$NewColumn
  GroupName <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GroupName
  
  
  trans.means <- cbind.data.frame(rn= names(new.sample.df))
  ##Recalculate transcode means 
  for(Transcode in unique(variance.tag2[,GrabColumn]) ){
    
    # print(Transcode)
    top.genes<- variance.tag2[ variance.tag2[,GrabColumn] %in%  Transcode,]$Gene
    trans.means$Gene <- colMeans(new.sample.df[top.genes,])
    names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  }
  row.names(trans.means) <- trans.means$rn
  trans.means$rn <- NULL
  exp.level.tester <- as.data.frame(t(trans.means))
  exp.level.tester<-as.data.frame(t(t(exp.level.tester)/colSums(exp.level.tester)))
  
  group.info <- cbind.data.frame(Gene=row.names(exp.level.tester), Mean=rowMeans(exp.level.tester))
  group.info$Stdev <- apply(exp.level.tester,1, FUN=sd)
  group.scale<-(exp.level.tester-group.info$Mean)/group.info$Stdev
  
  #########squared scale tag
  group.scaled.bin.tag <- group.scale
  group.scaled.bin.tag<- group.scaled.bin.tag*group.scaled.bin.tag*(group.scaled.bin.tag/abs(group.scaled.bin.tag))
  
  
  ###Just using exp.level.tester here
  gene.cor <- cor( t(group.scaled.bin.tag) )
  diag(gene.cor)<-0
  

  tester.melt <- melt(gene.cor)
  tester.melt<- tester.melt[tester.melt$value>0.4,]
  tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
  
  tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
  
  tester.melt<-tester.melt[tester.melt$value>Corr.Cutoff,]
  
  variance.tag2$NewColumn <- variance.tag2[,GrabColumn]
  names(variance.tag2)[ncol(variance.tag2)] <- NewColumn
  

  for(Transcode in tester.melt$Var1 ){
    change.target2  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
    variance.tag2[variance.tag2[,GrabColumn] %in% Transcode,NewColumn] <- as.character(change.target2)
  }
  
  
  variance.tag2$Transcode <- variance.tag2[,NewColumn]
  dictionary.summary <- NULL
  tb.counter <-1
  for(Transcode in unique(variance.tag2[,NewColumn])){
    # print(Transcode)
    if(Transcode %in% dictionary.summary==T){next
      # print("next")
    }
    dictionary.spread<-Transcode
    tester.df <- variance.tag2[,c(GrabColumn, NewColumn)]
    n=1
    repeat{
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,NewColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,NewColumn])
      
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread<-unique(dictionary.spread)
      n=n+1
      if(n>3){break}
    }
    dictionary.spread<-unique(dictionary.spread)
    variance.tag2[variance.tag2[,NewColumn] %in% dictionary.spread,]$Transcode <- paste0(GroupName, tb.counter)
    tb.counter <- tb.counter+1
    dictionary.summary<- c(dictionary.summary, dictionary.spread)
  }
  
  
  variance.tag2[,NewColumn]<- variance.tag2$Transcode
  # print(paste0("Now have ", length(unique(variance.tag2[,NewColumn])), " transcodes compared to ", length(unique(variance.tag2[,GrabColumn]))) )
  
}



# write.table(variance.tag2, paste0(outdir, Exp.Name, ".vox_corr_table.step5.gene.table.txt"),
#             sep="\t")




##############
### Flag Core patterns 5


transcode.table <- as.data.frame(table(variance.tag2$P5CG5))
transcode.table2 <- transcode.table[transcode.table$Freq>9,]
transcode.table2$GroupMean <- 0
for(Transcode in transcode.table2$Var1){
  genes.of.interest <- variance.tag2[variance.tag2$P5CG5==Transcode,]$Gene
  transcode.table2[transcode.table2$Var1==Transcode,]$GroupMean <- mean(variance.tag2[variance.tag2$Gene %in% genes.of.interest,]$Mean)
  
}
transcode.table2 <- transcode.table2[transcode.table2$GroupMean > mean(variance.tag$Mean)/10,]

second.pass.genes <- variance.tag2[variance.tag2$P5CG5 %in% transcode.table2$Var1,]$Gene

variance.tag$Pass5 <- "F1"

variance.tag[variance.tag$FourthPass!="E1",]$Pass5 <- variance.tag[variance.tag$FourthPass!="E1",]$FourthPass
variance.tag[variance.tag$Gene %in% second.pass.genes,]$Pass5 <- variance.tag2[variance.tag2$P5CG5 %in% transcode.table2$Var1,]$P5CG5


# write.table(variance.tag, paste0(outdir, Exp.Name, ".vox_corr_table.step5_all.gene.table.txt"),
#             sep="\t")




############################
### Level F1 merges - no holds barred

variance.tag2 <- variance.tag[variance.tag$Pass5=="F1",]

p2.scaled.bin.tag <- scaled.tag[row.names(scaled.tag) %in% variance.tag2$Gene,]
p2.scaled.bin.tag<- p2.scaled.bin.tag*p2.scaled.bin.tag*(p2.scaled.bin.tag/abs(p2.scaled.bin.tag))
gene.cor <- cor( t(p2.scaled.bin.tag) )
diag(gene.cor)<-0

tester.melt <- melt(gene.cor)
tester.melt<- tester.melt[tester.melt$value>0.1,]
tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
# plot1<-ggplot(tester.melt)+
#   theme_classic()+
#   geom_density(aes(value))+
#   geom_vline(xintercept = c(  mean(tester.melt$value)+sd(tester.melt$value) ), color="red")+
#   geom_vline(xintercept=(0.2))
# ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Vox.pass6.intercorr_all.cutoff.distrubtion.png"),
#        height=3, width=5)

##change
tester.melt<- tester.melt[tester.melt$value>0.2,]
##change
variance.tag2$P6CG1 <- variance.tag2$Gene


for(Transcode in tester.melt$Var1 ){
  # print(paste0("Converting - ",Transcode))
  change.target  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
  change.target2 <- variance.tag2[variance.tag2$Gene==change.target,]$Gene
  # print(paste0("changing to " , change.target2))
  variance.tag2[variance.tag2$Gene==Transcode,]$P6CG1 <- as.character(change.target2)
}


variance.tag2$Transcode <- variance.tag2$P6CG1
dictionary.summary <- NULL
tb.counter <-1
for(Transcode in unique(variance.tag2$P6CG1)){
  # print(Transcode)
  if(Transcode %in% dictionary.summary==T){next
    # print("next")
    }
  dictionary.spread<-Transcode
  tester.df <- variance.tag2[,c("P6CG1","Gene")]
  n=1
  repeat{
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P6CG1 %in% dictionary.spread,]$P6CG1)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$P6CG1)
    
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$P6CG1 %in% dictionary.spread,]$Gene)
    dictionary.spread <- c(dictionary.spread,tester.df[tester.df$Gene %in% dictionary.spread,]$Gene)
    dictionary.spread<-unique(dictionary.spread)
    n=n+1
    if(n>3){break}
  }
  dictionary.spread<-unique(dictionary.spread)
  variance.tag2[variance.tag2$P6CG1 %in% dictionary.spread,]$Transcode <- paste0("FA", tb.counter)
  tb.counter <- tb.counter+1
  dictionary.summary<- c(dictionary.summary, dictionary.spread)
}
variance.tag2$P6CG1<- variance.tag2$Transcode
print(paste0("Now have ", length(unique(variance.tag2$P6CG1)), " transcodes compared to ", length(unique(variance.tag2$Gene))) )






#######Compression pass


GrabColumn="P2CG1"
NewColumn <- "P2CG2"
GroupName <-"BB"


corr.round.x.df <- cbind.data.frame(Corr.Cutoff=c(0.6,0.40,0.3,0.20),
                                    GrabColumn=c("P6CG1","P6CG2","P6CG3","P6CG4"),
                                    NewColumn=c("P6CG2","P6CG3","P6CG4","P6CG5"),
                                    GroupName=c("FB","FC","FD","FE") )



for(Corr.Cutoff in corr.round.x.df$Corr.Cutoff){
  # print(Corr.Cutoff)
  
  GrabColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GrabColumn
  NewColumn <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$NewColumn
  GroupName <- corr.round.x.df[corr.round.x.df$Corr.Cutoff==Corr.Cutoff,]$GroupName
  
  
  trans.means <- cbind.data.frame(rn= names(new.sample.df))
  ##Recalculate transcode means 
  for(Transcode in unique(variance.tag2[,GrabColumn]) ){
    
    # print(Transcode)
    top.genes<- variance.tag2[ variance.tag2[,GrabColumn] %in%  Transcode,]$Gene
    trans.means$Gene <- colMeans(new.sample.df[top.genes,])
    names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
  }
  row.names(trans.means) <- trans.means$rn
  trans.means$rn <- NULL
  exp.level.tester <- as.data.frame(t(trans.means))
  exp.level.tester<-as.data.frame(t(t(exp.level.tester)/colSums(exp.level.tester)))
  
  group.info <- cbind.data.frame(Gene=row.names(exp.level.tester), Mean=rowMeans(exp.level.tester))
  group.info$Stdev <- apply(exp.level.tester,1, FUN=sd)
  group.scale<-(exp.level.tester-group.info$Mean)/group.info$Stdev
  
  #########squared scale tag
  group.scaled.bin.tag <- group.scale
  group.scaled.bin.tag<- group.scaled.bin.tag*group.scaled.bin.tag*(group.scaled.bin.tag/abs(group.scaled.bin.tag))
  
  
  ###Just using exp.level.tester here
  gene.cor <- cor( t(group.scaled.bin.tag) )
  diag(gene.cor)<-0
  

  tester.melt <- melt(gene.cor)
  tester.melt<- tester.melt[tester.melt$value>0.2,]
  tester.melt<- tester.melt[with(tester.melt, order(Var1, -value)),]
  
  tester.melt<- tester.melt[!duplicated(tester.melt$Var1),]
  
  tester.melt<-tester.melt[tester.melt$value>Corr.Cutoff,]
  
  variance.tag2$NewColumn <- variance.tag2[,GrabColumn]
  names(variance.tag2)[ncol(variance.tag2)] <- NewColumn
  

  for(Transcode in tester.melt$Var1 ){
    change.target2  <- tester.melt[tester.melt$Var1==Transcode,]$Var2
    variance.tag2[variance.tag2[,GrabColumn] %in% Transcode,NewColumn] <- as.character(change.target2)
  }
  
  
  variance.tag2$Transcode <- variance.tag2[,NewColumn]
  dictionary.summary <- NULL
  tb.counter <-1
  for(Transcode in unique(variance.tag2[,NewColumn])){
    # print(Transcode)
    if(Transcode %in% dictionary.summary==T){next
      # print("next")
    }
    dictionary.spread<-Transcode
    tester.df <- variance.tag2[,c(GrabColumn, NewColumn)]
    n=1
    repeat{
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,NewColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,NewColumn])
      
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,NewColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread <- c(dictionary.spread,tester.df[tester.df[,GrabColumn] %in% dictionary.spread,GrabColumn])
      dictionary.spread<-unique(dictionary.spread)
      n=n+1
      if(n>3){break}
    }
    dictionary.spread<-unique(dictionary.spread)
    variance.tag2[variance.tag2[,NewColumn] %in% dictionary.spread,]$Transcode <- paste0(GroupName, tb.counter)
    tb.counter <- tb.counter+1
    dictionary.summary<- c(dictionary.summary, dictionary.spread)
  }
  
  
  variance.tag2[,NewColumn]<- variance.tag2$Transcode
  # print(paste0("Now have ", length(unique(variance.tag2[,NewColumn])), " transcodes compared to ", length(unique(variance.tag2[,GrabColumn]))) )
  
}


# 
# write.table(variance.tag2, paste0(outdir, Exp.Name, ".vox_corr_table.step5.gene.table.txt"),
#             sep="\t")




##############
### Flag Core patterns 5


transcode.table <- as.data.frame(table(variance.tag2$P6CG5))
transcode.table2 <- transcode.table[transcode.table$Freq>5,]
transcode.table2$GroupMean <- 0
for(Transcode in transcode.table2$Var1){
  genes.of.interest <- variance.tag2[variance.tag2$P6CG5==Transcode,]$Gene
  transcode.table2[transcode.table2$Var1==Transcode,]$GroupMean <- mean(variance.tag2[variance.tag2$Gene %in% genes.of.interest,]$Mean)
  
}
# transcode.table2 <- transcode.table2[transcode.table2$GroupMean > mean(variance.tag$Mean)/10,]

second.pass.genes <- variance.tag2[variance.tag2$P6CG5 %in% transcode.table2$Var1,]$Gene

variance.tag$Pass6 <- "G1"

variance.tag[variance.tag$Pass5!="F1",]$Pass6 <- variance.tag[variance.tag$Pass5!="F1",]$Pass5
variance.tag[variance.tag$Gene %in% second.pass.genes,]$Pass6 <- variance.tag2[variance.tag2$P6CG5 %in% transcode.table2$Var1,]$P6CG5

# 
# write.table(variance.tag, paste0(outdir, Exp.Name, ".vox_corr_table.step6_all.gene.table.txt"),
#             sep="\t", quote=F)




print(paste0("Now have ", length(unique(variance.tag[,"Pass6"])), " transcodes compared to ", length(unique(variance.tag[,"FirstPass"])), " starting") )










#############
## Blow up low intracorr groups and rechair genes

Transcode="AE3"

master.record.table <- NULL
for(Transcode in unique(variance.tag$Pass6)){
  # print(Transcode)
  genes.of.interest <- variance.tag[variance.tag$Pass6==Transcode,]$Gene
  
  tmp.norm.data <- normalized.data[row.names(normalized.data) %in% genes.of.interest,]
  
  
  
  tmp.cor <- cor(t(tmp.norm.data))
  
  tmp.record.table <- cbind.data.frame(Transcode=Transcode, Corr.Mean= mean (mean(tmp.cor)), Exp.Max=max(colMeans(tmp.norm.data)) )
  master.record.table<- rbind(master.record.table, tmp.record.table)
}

mean(master.record.table$Corr.Mean)


plot1 <- ggplot(master.record.table)+
  geom_point(aes(Corr.Mean,Exp.Max))+
  theme_classic()+
  xlim(c(0,0.2))+
  ylim(c(0,2))
ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Pass6.Innercorr.vs.ExpMax.scatter.png"),
       height=4, width=5)

plot1 <- ggplot(master.record.table)+
  geom_density(aes(Corr.Mean), fill="black")+
  geom_vline(xintercept = c(mean(master.record.table$Corr.Mean)), color="red")+
  geom_vline(xintercept = c(mean(master.record.table$Corr.Mean)+sd(master.record.table$Corr.Mean)), color="purple")+
  theme_classic()
# xlim(c(0,0.2))+
# ylim(c(0,2))
ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Pass6.Innercorr_perGene.Distribution.png"),
       height=4, width=10)



## get first letter for quality report
master.record.table$Tier <- sapply(strsplit(master.record.table$Transcode,""),"[[",1)


master.record.table2 <- master.record.table[with(master.record.table, order(Corr.Mean)),] 

head(master.record.table2,60)

master.record.table2$Qual <- "Low"
master.record.table2[master.record.table2$Corr.Mean>0.07,]$Qual <- "High"

dim(master.record.table[master.record.table<0.07,])

plot1 <- ggplot(master.record.table2)+
  geom_bar(aes(Tier,fill=Qual), stat="count")+
  # geom_vline(xintercept = c(mean(master.record.table$Corr.Mean)), color="red")+
  # geom_vline(xintercept = c(mean(master.record.table$Corr.Mean)+sd(master.record.table$Corr.Mean)), color="purple")+
  theme_classic()
# xlim(c(0,0.2))+
# ylim(c(0,2))
ggsave(plot=plot1, paste0(outdir, Exp.Name, ".Pass6.Innercorr_Qual.per.TranscodeTier.Bargraph.png"),
       height=4, width=10)





genes.of.interest<- c("MUC5B","KRT5","LTF","DNAH11","LYZ","TFF1","BPIFB1", "XBP1", "LTF","CFB",  "FCGBP")
master.record.table[master.record.table$Transcode %in% variance.tag[variance.tag$Gene %in% genes.of.interest,]$Pass6,]


###
fuzzy.transcodes <- master.record.table[master.record.table$Corr.Mean<0.1,]$Transcode

dim(variance.tag[variance.tag$Pass6 %in% fuzzy.transcodes,])



### Rechair 

print("Final retest of genes to see if they prefer a different emergent pattern")



#####################
### Rechair
GrabColumn <- "Pass6"
NewColumn <- "Group1"
GroupName <- "GA"

variance.tag$Rechair <- variance.tag$Pass6

variance.tag2 <- variance.tag[!variance.tag$Pass6 %in% fuzzy.transcodes,]


trans.means <- cbind.data.frame(rn= names(new.sample.df))
##Recalculate transcode means
for(Transcode in unique(variance.tag2[,GrabColumn]) ){
  
  # print(Transcode)
  top.genes<- variance.tag2[ variance.tag2[,GrabColumn] %in%  Transcode,]$Gene
  trans.means$Gene <- colMeans(new.sample.df[top.genes,])
  names(trans.means)[ncol(trans.means)] <- paste0(Transcode)
}
row.names(trans.means) <- trans.means$rn
trans.means$rn <- NULL
exp.level.tester <- as.data.frame(t(trans.means))
exp.level.tester<-as.data.frame(t(t(exp.level.tester)/colSums(exp.level.tester)))


group.info <- cbind.data.frame(Gene=row.names(exp.level.tester), Mean=rowMeans(exp.level.tester))
group.info$Stdev <- apply(exp.level.tester,1, FUN=sd)
group.scale<-(exp.level.tester-group.info$Mean)/group.info$Stdev
# squared.group <- group.scale*group.scale *(group.scale/abs(group.scale))


Group="AE395"
##Testing with VoxGroup will have to rewrite to Transcode
for(Group in unique(variance.tag[variance.tag$Pass6 %in% fuzzy.transcodes,]$Rechair)){
  # print(Group)
  genes.of.interest <- variance.tag[variance.tag$Rechair %in% Group,]$Gene
  
  ##For full gene
  # tmp.segment <- normalized.data[row.names(normalized.data) %in% genes.of.interest,]
  
  ##For gene summary (WAY FASTER)
  # tmp.segment <- (new.sample.df[genes.of.interest,])
  # rechair.corr<-cor(t(tmp.segment),t(exp.level.tester))
  
  tmp.segment <- (scaled.tag[genes.of.interest,])
  rechair.corr<-cor(t(tmp.segment),t(group.scale))
  
  # print(table(colnames(rechair.corr)[apply(rechair.corr, 1, which.max)]))
  rechair.melt <- melt(rechair.corr)
  rechair.melt<- rechair.melt[with(rechair.melt, order(-value)),]
  rechair.melt <- rechair.melt[!duplicated(rechair.melt$Var1),]
  
  
  ##Too many crappy merges out there installing floor
  rechair.melt <- rechair.melt[rechair.melt$value>0.25,]
  
  ###For low ex[ressed genes, basicallyanything works, for high you need a higher corr floor
  #### for now just doing 0.25
  
  ##
  genes.of.interest <- genes.of.interest[genes.of.interest %in% rechair.melt$Var1]
  
  variance.tag[variance.tag$Gene %in% genes.of.interest,]$Rechair <- as.character(rechair.melt[match(variance.tag[variance.tag$Gene %in% genes.of.interest,]$Gene,rechair.melt$Var1),]$Var2)
}

nrow(variance.tag[variance.tag$Rechair!=variance.tag$Pass6,])




dim(variance.tag[!variance.tag$Rechair %in% fuzzy.transcodes,])

table(variance.tag[variance.tag$Rechair %in% fuzzy.transcodes,]$Rechair)

###Stick the rest in G1 for now??????
###Will make a the n1-n4 thing with a broad kmeans distribution


####
###Check if there is any that didnt find a home and classify as noise

if(length(variance.tag[variance.tag$Rechair %in% fuzzy.transcodes,]$Rechair)>0){

ending.cluster.num<- 5

markers$Cluster_mid <- as.factor(kmeans(markers[,c("X","Y")], iter.max=30, nstart=5, centers=ending.cluster.num)$cluster) 
markers$Cluster_mid <- paste0("C",markers$Cluster_mid)



# markers$Cluster <- markers$NewCluster
markers$meanX2 <- 0 ; markers$meanY2 <- 0
for(Cluster in unique(markers$Cluster_mid)){
  markers[markers$Cluster_mid==Cluster,]$meanX2 <- mean(markers[markers$Cluster_mid==Cluster,]$X)
  markers[markers$Cluster_mid==Cluster,]$meanY2 <- mean(markers[markers$Cluster_mid==Cluster,]$Y)
}

# markers$NextCluster <- NULL

palette.vector <- distinctColorPalette(length(unique(markers$Cluster_mid)))
# 
plot1<-ggplot(markers)+
  geom_point(aes(X,Y, color=(Cluster_mid)))+
  geom_text(aes(meanX2,meanY2, label=Cluster_mid), size=5)+
  scale_colour_discrete(type =palette.vector)+
  labs( x="UMAP 1", y="UMAP 2", color="Clusters")+
  theme_classic()+
  theme(legend.position = "")
# plot1
ggsave(plot=plot1,filename=paste0(outdir,  Exp.Name, ".MidClustersForNoise.", ending.cluster.num,"_labeled_nice.UMAP.png"),
       height=12, width=14)




###Have to rerun tag table now with new cluster names
tag.table2 <- NULL
tag.table2$rn <- row.names(normalized.data)
for( Sample1 in sort(unique(markers$Cluster_mid))  ){
  # print(Sample1)
  
  if(length(row.names(markers[markers$Cluster_mid %in% Sample1,]) )>1){
    
    tag.temp <-cbind.data.frame(rn=row.names(normalized.data),
                                Sample1=rowMeans(normalized.data[ ,colnames(normalized.data) %in% row.names(markers[markers$Cluster_mid %in% Sample1,]) ]) )
    tag.table2 <- cbind.data.frame(tag.table2, tag.temp$Sample1)
    names(tag.table2)[ncol(tag.table2)] <- paste(Sample1)} else{
      
      tag.temp <- normalized.data[ ,colnames(normalized.data) %in% row.names(markers[markers$Cluster_mid %in% Sample1,]) ]
      tag.table2$temp <- tag.temp
      names(tag.table2)[ncol(tag.table2)] <- paste(Sample1)
    }
}
row.names(tag.table2) <- tag.table2$rn
tag.table2$rn <- NULL



noise.table <-cbind.data.frame(Gene=variance.tag[variance.tag$Rechair %in% fuzzy.transcodes,]$Gene,
                               Cluster=1)


head(tag.table2[noise.table$Gene, ])

noise.table$MaxCluster<-colnames(tag.table2)[apply(tag.table2[noise.table$Gene, ],1,which.max)]

noise.table$Cluster <- paste0("G", noise.table$MaxCluster)

variance.tag[variance.tag$Rechair %in% fuzzy.transcodes,]$Rechair <- noise.table[match(variance.tag[variance.tag$Rechair %in% fuzzy.transcodes,]$Gene, noise.table$Gene),]$Cluster

}


length(unique(variance.tag$Rechair))

###We can stop here
## 144 after these two operations

vox.tag<- variance.tag

vox.tag$Transcode <- vox.tag$Rechair
###Can do a similar operation for supergroup








###Final Organization

##Depending on what was picked have to rerun above fingerprint? Just paste here....
new.sample.group.data <- cbind.data.frame(Sample=names(new.sample.df))
Group<- "G1"
for(Group in sort(unique(vox.tag$Transcode))){
  # print(Group)
  gene.list <- (vox.tag[vox.tag$Transcode==Group,]$Gene)
  new.sample.group.data <- cbind.data.frame( new.sample.group.data ,Score=colSums(new.sample.df[gene.list,]) ) 
  names(new.sample.group.data)[ncol(new.sample.group.data)]<- paste0(Group)
}
new.sample.group.data$Sample <- NULL
new.sample.group.data <-((new.sample.group.data)/rowSums(new.sample.group.data))
rowSums(new.sample.group.data)

enrich.data <- NULL
enrich.data <- cbind.data.frame(Group=names(new.sample.group.data))
for(Sample in names(new.sample.df)){
  group.mean.tester <- new.sample.group.data
  group.mean.tester$Sample <- NULL
  diff.groups.tester<-(group.mean.tester[Sample,]+0.01)/(colMeans(group.mean.tester[row.names(group.mean.tester)!=Sample,])+0.01)
  diff.groups.tester<-unlist(diff.groups.tester)
  enrich.data<-cbind(enrich.data, Data=diff.groups.tester)
  names(enrich.data)[ncol(enrich.data)] <- Sample
  
}

row.names(enrich.data) <- enrich.data$Group
enrich.data$Group<-NULL

ml3.table<-cbind.data.frame(Group=names(apply(enrich.data,1, max)),
                            Max=unlist(apply(enrich.data,1, max)),
                            Min=unlist(apply(enrich.data,1, min)))
ml3.table$Max_Up <- (ml3.table$Max-1)
ml3.table$Min_Up <- (1-ml3.table$Min)
ml3.table$AvgEnrich <- (ml3.table$Max_Up+ ml3.table$Min_Up)/2
ml3.table$Enrich_Bin <- cut(ml3.table$AvgEnrich, labels=F, breaks=20)


ml3.table <- ml3.table[with(ml3.table, order(-Enrich_Bin)),]

# height=4, width=5)
plot1<-ggplot(ml3.table)+
  geom_bar(aes(Group, AvgEnrich), stat="identity",fill="black")+
  theme_classic()+
  geom_hline(yintercept = median(ml3.table$AvgEnrich), color="red")+
  theme(axis.text.x=element_blank())
ggsave(plot=plot1,paste0(outdir, Exp.Name, ".VoxGroup.Enrichment.data.bargraph.png"),
       height = 3, width=4)



pheatmap(enrich.data/rowSums(enrich.data),
         filename=paste0(outdir, Exp.Name, ".VoxEnrich.Data.heatmap.rowNorm.png"),
         height=8, width=12)






### Method to order VOX groups by seurat clustering using MAX expression

##Rename clusters by order post agreement
markers.temp <- markers 
markers.temp$XY <- (markers.temp$X)^2+(markers.temp$Y)
markers.temp<- markers.temp[with(markers.temp, order(-markers.temp$XY)),]
n=1
markers$NewCluster <- "C_0"
for(NewCluster in unique(markers.temp$seurat_cluster)){
  markers[markers$seurat_cluster==NewCluster,]$NewCluster <- paste0("C_", n)
  n=n+1
}


markers$Cluster_mid<- markers$NewCluster



###Have to rerun tag table now with new cluster names
tag.table2 <- NULL
tag.table2$rn <- row.names(normalized.data)
for( Sample1 in sort(unique(markers$Cluster_mid))  ){
  # print(Sample1)
  
  if(length(row.names(markers[markers$Cluster_mid %in% Sample1,]) )>1){
    
    tag.temp <-cbind.data.frame(rn=row.names(normalized.data),
                                Sample1=rowMeans(normalized.data[ ,colnames(normalized.data) %in% row.names(markers[markers$Cluster_mid %in% Sample1,]) ]) )
    tag.table2 <- cbind.data.frame(tag.table2, tag.temp$Sample1)
    names(tag.table2)[ncol(tag.table2)] <- paste(Sample1)} else{
      
      tag.temp <- normalized.data[ ,colnames(normalized.data) %in% row.names(markers[markers$Cluster_mid %in% Sample1,]) ]
      tag.table2$temp <- tag.temp
      names(tag.table2)[ncol(tag.table2)] <- paste(Sample1)
    }
}
row.names(tag.table2) <- tag.table2$rn
tag.table2$rn <- NULL

louvain.group.data <- cbind.data.frame(Sample=names(tag.table2))
Group<- "G1"
for(Group in sort(unique(vox.tag$Transcode))){
  # print(Group)
  gene.list <- (vox.tag[vox.tag$Transcode==Group,]$Gene)
  louvain.group.data <- cbind.data.frame( louvain.group.data ,Score=colSums(tag.table2[gene.list,]) ) 
  names(louvain.group.data)[ncol(louvain.group.data)]<- paste0(Group)
}
louvain.group.data$Sample <- NULL
louvain.group.data <-((louvain.group.data)/rowSums(louvain.group.data))
rowSums(louvain.group.data)




sg.table <-cbind.data.frame(Group=colnames(louvain.group.data),
                            Cluster=1)


# head(tag.table2[noise.table$Gene, ])

sg.table$MaxCluster<-rownames(louvain.group.data)[apply(louvain.group.data[, ],2,which.max)]

####Should probably come back and rewrite this to be a max specificity or max differential term
sg.table$MaxExp <- apply(louvain.group.data,2,max)


order.of.clusters <- paste0("C_", c(1:length(unique(markers$seurat_clusters))))
sg.table$MaxCluster <- factor(sg.table$MaxCluster, levels=c(order.of.clusters,  paste0("C_",(length(unique(markers$seurat_clusters))+1) ) ))

if(length(grep("G",sg.table$Group))>0 ){
sg.table[grepl("G", sg.table$Group),]$MaxCluster<- paste0("C_",(length(unique(markers$seurat_clusters))+1) )
}

sg.table<-sg.table[with(sg.table, order(MaxCluster, -MaxExp)),]
sg.table$Group_New <- paste0("V" ,c(1:nrow(sg.table)) )

if(length(grep("G",sg.table$Group))>0 ){
sg.table[grepl("G", sg.table$Group),]$Group_New <- paste0("Z.",sg.table[grepl("G", sg.table$Group),]$Group_New)
}

vox.tag$VoxGroup <- sg.table[match(vox.tag$Transcode,sg.table$Group),]$Group_New

length(unique(vox.tag$VoxGroup))
table(vox.tag$VoxGroup)


vox.tag$VoxGroup<- factor(vox.tag$VoxGroup, levels=unique(sg.table$Group_New))
vox.tag<-vox.tag[with(vox.tag, order(vox.tag$VoxGroup)),]

vox.tag$MetaLayerEnrichment <- ml3.table[match(vox.tag$Transcode, ml3.table$Group),]$AvgEnrich

vox.tag$MetaLayerCorrelation <- 0

for(Group in unique(vox.tag$VoxGroup)){
  # print(Group)
  genes.of.interest<- vox.tag[vox.tag$VoxGroup==Group,]$Gene
  norm.data.chunk <- normalized.data[row.names(normalized.data) %in% genes.of.interest,]
  
  score.piece<-colMeans(normalized.data[row.names(normalized.data) %in% genes.of.interest, ])
  score.cor <- cor(t(norm.data.chunk), score.piece)
  
  
  final.gene.tester<- cbind.data.frame(Gene=row.names(score.cor),Specificity=0, Score_Corr =as.numeric(score.cor))
  
  vox.tag[vox.tag$VoxGroup==Group,]$MetaLayerCorrelation <- final.gene.tester[match(vox.tag[vox.tag$VoxGroup==Group,]$Gene,final.gene.tester$Gene),]$Score_Corr
  
}

vox.tag[is.na(vox.tag)]<-0






vox.tag<-vox.tag[with(vox.tag, order(vox.tag$VoxGroup, -MetaLayerCorrelation)),]





write.table(vox.tag, paste0(outdir, Exp.Name,"_geneCHORUS_v12.6_raw.VoxGroups.txt"),
            sep="\t", quote=F)







# write.table(vox.tag, paste0(outdir, Exp.Name,"_geneCHORUS_v12.VoxGroups.txt"),
#             sep="\t", quote=F)
# 
# 

vox.tag$VoxGroup<- factor(vox.tag$VoxGroup, levels=c(unique(sg.table$Group_New),"Noise"))



vox.tag[vox.tag$MetaLayerCorrelation<noise_floor,]$VoxGroup <- "Noise"
network.df <- as.data.frame(table(vox.tag$VoxGroup))
if(nrow(network.df[network.df$freq<10,])>0){
  vox.tag[vox.tag$VoxGroup %in% network.df[network.df$Freq<10,]$Var1,]$VoxGroup <- "Noise"
}

# 
# write.table(vox.tag, paste0(outdir, Exp.Name,"_geneCHORUS_v12.VoxGroups.txt"),
#             sep="\t", quote=F)


#################
## Data Quality Segment, plot and grade networks

print("Testing data quality and outputting QA plots")

################################
###Plot group expression on UMAP


dir.create(paste0(outdir, "VoxGroup.Expression/"))
outdir2<- paste0(outdir, "VoxGroup.Expression/")

# variance.tag$ML1_Index <- paste0(variance.tag$MetaLayer2,"-",variance.tag$MetaLayer1)

transcode.table <- as.data.frame(table(vox.tag$VoxGroup))

transcode.table<- transcode.table[transcode.table$Var1!="Noise",]
# n=1
# Exp.Name<- "YS3.Filt"

Group="V10"

expression.level.table <- NULL
for(Group in unique(transcode.table$Var1)){
  # print(Group)
  plot.data <- markers
  genes.of.interest <- vox.tag[vox.tag$VoxGroup==Group,]$Gene
  plot.data$Score <- colMeans(normalized.data[row.names(normalized.data) %in% genes.of.interest,colnames(normalized.data) %in% row.names(plot.data)])
  
  amount.of.genes<-transcode.table[transcode.table$Var1==Group,]$Freq
  
  plot.data <- plot.data[with(plot.data, order(Score)),]
  # ggplot(plot.data)+
  #   # facet_wrap(~Type)+
  #   geom_point(data=markers[,c("X","Y")], aes(X,Y), color="grey55")+
  #   geom_point( aes(X,Y, color=Score) )+
  #   theme_classic()+
  #   labs(color=Group)+
  #   theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  #   ggtitle(paste0("VOX ",Group, " with ", amount.of.genes," genes"))+
  #   scale_colour_gradientn(colours=c("lightcyan3", "gold", "orange","orangered","purple"))
  # ggsave(filename=paste0(outdir2,  Exp.Name, ".Vox.avgExp.VoxGroup_",Group,".UMAP.png"),
  #        height=4, width=6)
  
  
  
  
  ##Intercorrelation bit
  
  temp.norm.data <- normalized.data[row.names(normalized.data) %in% genes.of.interest,]
  temp.corr.data <- cor(t(temp.norm.data))
  diag(temp.corr.data) <- 0
  

  pheatmap(temp.corr.data,
           show_rownames = F, show_colnames = F,
           filename=paste0(outdir2,  Exp.Name, ".Vox.Intercorr.VoxGroup_",Group,".UMAP.png"),
           height=6, width=7)
  
  ##This reports the scaled maximum per vox 
  # length(plot.data$Score[plot.data$Score> mean(sort(plot.data$Score, decreasing=TRUE)[1:10] ) /2 ])
  # 
  tmp.exp.table <- cbind.data.frame(Vox=Group, MetaSD=  sd(plot.data$Score), 
                                    MetaBeatsMax=length(plot.data$Score[plot.data$Score> mean(sort(plot.data$Score, decreasing=TRUE)[1:10] ) /2 ]),
                                    MetaIntercorr=mean(temp.corr.data) )
  expression.level.table <- rbind(expression.level.table, tmp.exp.table)
  
}












expression.level.table$Vox <- factor(expression.level.table$Vox, levels=unique(vox.tag$VoxGroup))
plot1<-ggplot(expression.level.table)+
  geom_bar(aes(Vox, MetaSD,fill=MetaSD), stat="identity")+
  scale_fill_gradientn(colors=c("grey55", "orange","orangered","purple"))+
  theme_classic()+
  scale_y_continuous( expand = c(0, 0)) +
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(plot=plot1, filename=paste0(outdir, Exp.Name, ".VoxGroup.MetaSD.Summary.bargraph.png"),
       height=4, width=25)

# expression.level.table$Vox <- factor(expression.level.table$Vox, levels=unique(vox.tag$VoxGroup))
plot1<-ggplot(expression.level.table)+
  geom_bar(aes(Vox, MetaBeatsMax,fill=MetaBeatsMax), stat="identity")+
  scale_fill_gradientn(colors=c("skyblue2", "orange4","orangered4","purple3"))+
  theme_classic()+
  labs(y="Cells over MAX/2")+
  scale_y_continuous( expand = c(0, 0)) +
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(plot=plot1, filename=paste0(outdir, Exp.Name, ".VoxGroup.MetaBeatsMax.Summary.bargraph.png"),
       height=4, width=25)

plot1<-ggplot(expression.level.table)+
  geom_bar(aes(Vox, MetaIntercorr,fill=MetaIntercorr), stat="identity")+
  scale_fill_gradientn(colors=c("black", "darkgreen","green3","orange"))+
  theme_classic()+
  scale_y_continuous( expand = c(0, 0)) +
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(plot=plot1, filename=paste0(outdir, Exp.Name, ".VoxGroup.MetaIntercorr.Summary.bargraph.png"),
       height=4, width=25)
##This is it, this is awesome


metalayer.table <- vox.tag[!duplicated(vox.tag$VoxGroup),c("VoxGroup","MetaLayerEnrichment")]
metalayer.table$VoxGroup <- factor(metalayer.table$VoxGroup, levels=unique(vox.tag$VoxGroup))

plot1<-ggplot(metalayer.table)+
  geom_bar(aes(VoxGroup, MetaLayerEnrichment,fill=MetaLayerEnrichment), stat="identity")+
  scale_fill_gradientn(colors=c("grey55", "skyblue","orange","red2"))+
  theme_classic()+
  labs(y="Enrichment")+
  scale_y_continuous( expand = c(0, 0)) +
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(plot=plot1, filename=paste0(outdir, Exp.Name, ".VoxGroup.Enrichment.Summary.bargraph.png"),
       height=4, width=25)



expression.level.table$MetaEnrichment <- metalayer.table[match(expression.level.table$Vox,metalayer.table$VoxGroup),]$MetaLayerEnrichment




################################
###Plot simulcount on UMAP


dir.create(paste0(outdir, "VoxGroup.Simulcount/"))
outdir2<- paste0(outdir, "VoxGroup.Simulcount/")

# variance.tag$ML1_Index <- paste0(variance.tag$MetaLayer2,"-",variance.tag$MetaLayer1)

transcode.table <- as.data.frame(table(vox.tag$VoxGroup))
# n=1
# Exp.Name<- "YS3.Filt"


Group="V140"

cohesive.table<-NULL
for(Group in unique(vox.tag$VoxGroup)){
  # print(Group)
  plot.data <- markers
  genes.of.interest <- vox.tag[vox.tag$VoxGroup==Group,]$Gene
  
  temp.data<-(new.sample.df[row.names(new.sample.df) %in% genes.of.interest, ])
  temp.data[temp.data>0]<-1
  
  score.table <- cbind.data.frame(Score=colSums(temp.data), Cluster=colnames(new.sample.df))
  plot.data$Score <- score.table[match(plot.data$Cluster,score.table$Cluster),]$Score
  amount.of.genes<-transcode.table[transcode.table$Var1==Group,]$Freq
  
  plot.data <- plot.data[with(plot.data, order(Score)),]
  plot1<-ggplot(plot.data)+
    # facet_wrap(~Type)+
    geom_point(data=markers[,c("X","Y")], aes(X,Y), color="grey55")+
    geom_point( aes(X,Y, color=Score) )+
    theme_classic()+
    labs(color=Group)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle(paste0("VOX ",Group, " Simulcount with ", amount.of.genes," genes"))+
    scale_colour_gradientn(colours=c("grey55", "blue", "purple","orange","cyan"))
  ggsave(plot=plot1,filename=paste0(outdir2,  Exp.Name, ".Vox.avgExp.VoxGroup_",Group,".UMAP.png"),
         height=4, width=6)
  
  
  #this takes the average of the top 10 simulcount cells 
  # mean(sort((colSums(temp.data)), decreasing=TRUE)[1:10])
  # tmp.cohesive.table <- cbind.data.frame(Vox=Group, MetaCohesiveness=  mean(sort((colSums(temp.data)), decreasing=TRUE)[1:10])/amount.of.genes   )
  ### Max cohesiveness wasnt that great anyway, need to come up with new simulcount metric
  ### For me it is the distance from the max or min
  
  ##Calculate Max and min, then subtract both and take the lesser sum
  
  
  score.table$MinDiff <- score.table$Score - min(score.table$Score)
  score.table$MaxDiff <- max(score.table$Score) - score.table$Score
  score.table$CohesiveDistance <- apply(score.table[,c("MinDiff","MaxDiff")], 1, min)/max(score.table$Score)
  
  
  tmp.cohesive.table <- cbind.data.frame(Vox=Group, MetaCohesiveness= 0.5-mean(score.table$CohesiveDistance))
  cohesive.table <- rbind(cohesive.table, tmp.cohesive.table)
}



cohesive.table$Vox <- factor(cohesive.table$Vox, levels=unique(vox.tag$VoxGroup))
plot1<-ggplot(cohesive.table)+
  geom_bar(aes(Vox, MetaCohesiveness,fill=MetaCohesiveness), stat="identity")+
  scale_fill_gradientn(colors=c("grey55", "blue", "purple","orange","cyan"))+
  theme_classic()+
  scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  theme(axis.text.x=element_text(angle=90, hjust=1))
ggsave(plot=plot1, filename=paste0(outdir, Exp.Name, ".VoxGroup.Simulcount.Cohesiveness.Summary.bargraph.png"),
       height=4, width=25)

# vox.tag$MetaCohesiveness <- cohesive.table[match(vox.tag$VoxGroup,cohesive.table$Vox),]$MetaCohesiveness



####Integrated Quality Plot
cohesive.table<- cohesive.table[cohesive.table$Vox!="Noise",]

norm.level.table<-cbind.data.frame(Vox=expression.level.table$Vox,
                                   MetaIntercorr=expression.level.table$MetaIntercorr/max(expression.level.table$MetaIntercorr),
                                   MetaSD=expression.level.table$MetaSD/(max(expression.level.table$MetaSD)),
                                   MetaEnrichment=expression.level.table$MetaEnrichment/max(expression.level.table$MetaEnrichment),
                                   MetaCohesiveness= cohesive.table$MetaCohesiveness/max(cohesive.table$MetaCohesiveness))

# MetaBeatsMax MetaIntercorr MetaEnrichment

heatmap.data <- norm.level.table
row.names(heatmap.data) <- heatmap.data$Vox
heatmap.data$Vox<- NULL


heatmap.colors = colorRampPalette(c("blue", "lightblue4", "lightblue1","darkseagreen1","yellow","yellow","orange","orange", "orangered","orangered"))(50)



pheatmap(t(heatmap.data),
         color=heatmap.colors,
         cluster_cols=F, cluster_rows=F,
         filename=paste0(outdir, Exp.Name, ".VoxGroup.Quality.Metrics.Summary.heatmap.png"),
         height=2.5, width=25)

###Cutoff version for classification

heatmap.data2 <-heatmap.data

# MetaIntercorr     MetaSD MetaEnrichment MetaCohesiveness

cutoff.val<-mean(heatmap.data2$MetaIntercorr)
heatmap.data2$MetaIntercorr[heatmap.data2$MetaIntercorr>=cutoff.val/2] <- 1
heatmap.data2$MetaIntercorr[heatmap.data2$MetaIntercorr<cutoff.val/2 & heatmap.data2$MetaIntercorr>cutoff.val/3] <- 0.5
heatmap.data2$MetaIntercorr[heatmap.data2$MetaIntercorr<=cutoff.val/3] <- 0

cutoff.val<-mean(heatmap.data2$MetaSD)
heatmap.data2$MetaSD[heatmap.data2$MetaSD>=cutoff.val/2] <- 1
heatmap.data2$MetaSD[heatmap.data2$MetaSD<cutoff.val/2 & heatmap.data2$MetaSD>cutoff.val/3] <- 0.5
heatmap.data2$MetaSD[heatmap.data2$MetaSD<=cutoff.val/3] <- 0

cutoff.val<-mean(heatmap.data2$MetaEnrichment)
heatmap.data2$MetaEnrichment[heatmap.data2$MetaEnrichment>=cutoff.val/2] <- 1
heatmap.data2$MetaEnrichment[heatmap.data2$MetaEnrichment<cutoff.val/2 & heatmap.data2$MetaEnrichment>cutoff.val/3] <- 0.5
heatmap.data2$MetaEnrichment[heatmap.data2$MetaEnrichment<=cutoff.val/3] <- 0

cutoff.val<-mean(heatmap.data2$MetaCohesiveness)
heatmap.data2$MetaCohesiveness[heatmap.data2$MetaCohesiveness>=cutoff.val/1.2] <- 1
heatmap.data2$MetaCohesiveness[heatmap.data2$MetaCohesiveness<cutoff.val/1.2 & heatmap.data2$MetaCohesiveness>cutoff.val/1.8] <- 0.5
heatmap.data2$MetaCohesiveness[heatmap.data2$MetaCohesiveness<=cutoff.val/1.8] <- 0


heatmap.colors<-c("white","grey65","black")

pheatmap(t(heatmap.data2),
         color=heatmap.colors,
         cluster_cols=F, cluster_rows=F,
         filename=paste0(outdir, Exp.Name, ".VoxGroup.Quality.Metrics.Summary_Binary.heatmap.png"),
         height=2.5, width=25)

### So here at the end we can call Network Quality based on these binaries


naming.df <- heatmap.data2

naming.df[naming.df==1]<-"A"
naming.df[naming.df==0.5]<-"B"
naming.df[naming.df==0]<-"C"

# MetaIntercorr     MetaSD MetaEnrichment MetaCohesiveness

norm.level.table$QualScore <- paste0(naming.df$MetaIntercorr,".",naming.df$MetaSD,".",naming.df$MetaEnrichment,".",naming.df$MetaCohesiveness)


write.table(norm.level.table, paste0(outdir, Exp.Name, ".VoxGroup.Quality.Table.txt"),
            sep="\t", quote=F)



###Thought about reordering dynamically based on quality
##  but really similar grouping type seems more useful

###Okay but we can drop the really bad groups


cutoff.table <- cbind.data.frame(Group=rownames(heatmap.data2), QualityScore=rowSums(heatmap.data2))
cutoff.table$Qual <- 1
cutoff.table[cutoff.table$QualityScore < 2,]$Qual <- 0

posthoc.noise.groups <- cutoff.table[cutoff.table$Qual==0,]$Group


# vox.tag$VoxGroup<- factor(vox.tag$VoxGroup, levels=c(unique(sg.table$Group_New), paste0("Z.",posthoc.noise.groups))  )

# vox.tag[vox.tag$VoxGroup %in% posthoc.noise.groups,]$VoxGroup <- paste0( "Z.",vox.tag[vox.tag$VoxGroup %in% posthoc.noise.groups,]$VoxGroup)


vox.tag<-vox.tag[with(vox.tag, order(vox.tag$VoxGroup, -MetaLayerCorrelation)),]

vox.tag$QualityScore <- norm.level.table[match( vox.tag$VoxGroup,norm.level.table$Vox),]$QualScore

# 
# write.table(vox.tag, paste0(outdir, Exp.Name,"_geneCHORUS_v12.6_raw.VoxGroups.txt"),
#             sep="\t", quote=F)


print("networks written")

# 
# vox.tag[vox.tag$MetaLayerCorrelation<0.1,]$VoxGroup <- "Noise"
# network.df <- as.data.frame(table(vox.tag$VoxGroup))
# if(nrow(network.df[network.df$freq<10,])>0){
#   vox.tag[vox.tag$VoxGroup %in% network.df[network.df$Freq<10,]$Var1,]$VoxGroup <- "Noise"
# }
# 

write.table(vox.tag, paste0(outdir, Exp.Name,"_geneCHORUS_v12.6_filtered.VoxGroups.txt"),
            sep="\t", quote=F)







################################
###Plot group expression on UMAP

print("plotting Correlation Group Expression")

dir.create(paste0(outdir, "VoxGroup.Expression/"))
outdir2<- paste0(outdir, "VoxGroup.Expression/")

transcode.table <- as.data.frame(table(vox.tag$VoxGroup))

for(Group in unique(vox.tag$VoxGroup)){
  # print(Group)
  plot.data <- markers
  genes.of.interest <- vox.tag[vox.tag$VoxGroup==Group,]$Gene
  plot.data$Score <- colMeans(normalized.data[row.names(normalized.data) %in% genes.of.interest,colnames(normalized.data) %in% row.names(plot.data)])
  
  amount.of.genes<-transcode.table[transcode.table$Var1==Group,]$Freq
  
  plot.data <- plot.data[with(plot.data, order(Score)),]
  ggplot(plot.data)+
    # facet_wrap(~Type)+
    geom_point(data=markers[,c("X","Y")], aes(X,Y), color="grey55")+
    geom_point( aes(X,Y, color=Score) )+
    theme_classic()+
    labs(color=Group)+
    theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
    ggtitle(paste0("VOX ",Group, " with ", amount.of.genes," genes"))+
    scale_colour_gradientn(colours=c("lightcyan3", "gold", "orange","orangered","purple"))
  ggsave(filename=paste0(outdir2,  Exp.Name, ".Vox.avgExp.VoxGroup_",Group,".UMAP.png"),
         height=4, width=6)
  
}













########################################
######## Create Outputs for metascape - Keeps order

print("printing outputs for metascape gene ontology")

##All Genes
Program.Table3 <- vox.tag
Program.Table3$FinalSet <- Program.Table3$VoxGroup
Program.Table3$Genes <- Program.Table3$Gene
Network.Source<- "VOX_Groups_All"


if(length(unique(Program.Table3$FinalSet))>28){
  print("Too many networks, dividing into parts")
  Cut.Amount <- c(1:ceiling( length(unique(Program.Table3$FinalSet))/28))
  Take.Amount <- 1
  Network.Cut.Table <- cbind.data.frame(Network=unique(Program.Table3$FinalSet),Cut.Identity=0)
  for(Cutn in c(1: max(Cut.Amount))){
    top.amount <- min( c(Take.Amount+round(nrow(Network.Cut.Table)/Cut.Amount), nrow(Network.Cut.Table) ) )
    Network.Cut.Table[c(Take.Amount:top.amount),]$Cut.Identity<-  Cutn
    Take.Amount <- Take.Amount+ round(nrow(Network.Cut.Table)/ max(Cut.Amount))+1
  }
  Program.Table3$Cut.Identity <- Network.Cut.Table[match(Program.Table3$FinalSet, Network.Cut.Table$Network),]$Cut.Identity
  for(Cut.Identity in unique(Program.Table3$Cut.Identity)){
    print(paste0("part ", Cut.Identity))
    Metascape.Multi <- NULL
    for(Set in unique(Program.Table3[Program.Table3$Cut.Identity==Cut.Identity,]$FinalSet)){
      tmp.Prog.Table <- Program.Table3[Program.Table3$FinalSet==Set,]
      Metascape.tmp<-cbind.data.frame( Name=paste0(Set),Genes=paste(as.character(tmp.Prog.Table$Genes), collapse=",")  )
      Metascape.Multi <- rbind(Metascape.Multi, Metascape.tmp)
    }
    print("done parsing")
    write.table(Metascape.Multi, paste0(outdir, Exp.Name,".",Network.Source,".Part",Cut.Identity,".Metanetworks.MultiGO.List.txt"), sep="\t", col.names =F, row.names = F,quote=F)
    
  }
}else{
  ##One loop version
  print("Network Numbers Okay")
  Metascape.Multi <- NULL
  for(Set in unique(Program.Table3$Group)){
    Program.Table.temp <- Program.Table3[Program.Table3$Group==Set,]
    Metascape.tmp<-cbind.data.frame( Name=paste0(Set),Genes=paste(as.character(Program.Table.temp$Genes), collapse=",")  )
    Metascape.Multi <- rbind(Metascape.Multi, Metascape.tmp)
  }
  print("done full")
  write.table(Metascape.Multi, paste0(outdir, Exp.Name,".VoxGroups.MultiGO.List.txt"), sep="\t", col.names =F, row.names = F,quote=F)
}




vox.tag

}


# install.packages("remotes")
# remotes::install_github("davidsjoberg/ggsankey")
# 

# install.packages("ggplot2")

# install.packages("dplyr")


# 
# sankey.df <- vox.tag[,c("FirstPass","SecondPass","ThirdPass","FourthPass","Pass5","Pass6", "Rechair", "VoxGroup")] %>%
#   make_long(FirstPass, SecondPass, ThirdPass, FourthPass, Pass5,Pass6,Rechair,VoxGroup)
# 
# 
# 
# plot1<- ggplot(sankey.df, aes(x = x, 
#                               next_x = next_x, 
#                               node = node, 
#                               next_node = next_node,
#                               fill = factor(node))) +
#   geom_sankey() +
#   theme_sankey(base_size = 16)+
#   theme(legend.position="")
# ggsave(plot=plot1,filename=paste0(outdir,  Exp.Name, ".Vox.Sankey.png"),
#        height=5, width=10)
# 
# 
# 
# 
# 
# 
# 
#   vox.tag
# }
# 














