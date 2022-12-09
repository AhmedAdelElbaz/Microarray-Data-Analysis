###########################################
#  - Multi Omics Course                   # 
#  - Gene Expression Profiling - Cisplatin#
#  - treated and untreated samples        #
#  - GSE77515 Microarray Analysis         #    
#  - 2022- 12-09                          #
#  - Copyright: @Ahmed Elbaz              #
###########################################
# R version 4.2.0 (2022-04-22)
# Platform: x86_64-apple-darwin17.0 (64-bit)
# Running under: macOS Monterey 12.6
# 
# Matrix products: default
# LAPACK: /Library/Frameworks/R.framework/Versions/4.2/Resources/lib/libRlapack.dylib
# 
# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

####load libraries ########  
library(readxl)
library(readr)
library(ggfortify)
library(matrixStats)
library(ComplexHeatmap)
library(genefilter)
library(multtest)
library(Biobase)
library(affy)
library(affyPLM)
library(pd.hugene.1.0.st.v1)
library(hugene10stv1cdf)
BiocManager::install('hugene10sttranscriptcluster.db')
library(hugene10sttranscriptcluster.db)
library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)

# Loading data and Create AffyObject and Create the Phenotable: --------------------------------------

celFilesDirectory="D:/MAIN/STUDY/PROGRAMMING/R/MICROARRAY DATA SETS ANALYSIS/GSE77515"

cels = list.files(celFilesDirectory, pattern = "CEL")
  cels
  rm(cels)

affyData <- ReadAffy(celfile.path=celFilesDirectory)
  affyData

phenotable <- read_excel(paste0(celFilesDirectory,"/phenotable.xlsx"))
  View(phenotable)
  
##############Exception Here################
samples <- c()
for (sample in sampleNames(affyData)){
  name <- strsplit(sample, "_")[[1]][1]
  samples <- c(samples , name)
}
rm(name)
sampleNames(affyData) <- samples
#############Exception Here################# 
  class(affyData)
  sampleNames(affyData)
  featureNames(affyData)
  head(featureNames(affyData))
  annotation(affyData)
  dim(affyData)
  dim(exprs(affyData))
# see how the RAW expression look like without processing : notice the large values, exprs:to extract the expression matrix
  head(Biobase::exprs(affyData))
  dim(Biobase::exprs(affyData))

exp.raw= Biobase::exprs(affyData)

#Save:
write.csv(exp.raw , paste0(celFilesDirectory,"exp.raw.csv"))

#Convert into Tidy Data ....Takes 5 Minutes
exp.raw_Tidy <-  as_tibble(exp.raw) %>% gather(samples , expression)%>% merge(phenotable ,  by = 'samples' , all.x = T)
  names(exp.raw_Tidy)
  head(exp.raw_Tidy)


exp.raw_Tidy_gg <- ggplot2::ggplot(exp.raw_Tidy)

#Histogram
exp.raw_Tidy_gg + 
  geom_density(aes(x = log(expression) , 
                   group = samples , 
                   fill = treated) , 
               adjust = 1.2 , alpha = .4) + 
  theme_classic() +
  xlab('Flouresence Dye Intensity') +
  ylab('Density')

#Boxplot
exp.raw_Tidy_gg + 
  geom_boxplot(aes(x = samples, 
                   y = log(expression) , 
                   fill = samples)) + 
  theme(axis.text.x = element_text(angle = 90))

#Prinipal Component Analysis
PCA <- prcomp(t(exp.raw),
              center = T,
              scale = T)
  dim(PCA$x)

#Plot PCA
autoplot(PCA,
         data = phenotable,
         colour = 'sample.type2',
         size = 2.5,
         #frame = T,
         #frame.type = 'norm'
         )



# Create eset Object from affyobject --------------------------------------

# threestep (background correction, normalization, summarization from probes to probesets)
eset = threestep(affyData,
                 background.method = "IdealMM",
                 normalize.method = "quantile",
                 summary.method = "average.log")
  head(eset@assayData$exprs)

#### if Generates NAs go to expresso
eset <- expresso(affyData, 
                 bgcorrect.method="rma", 
                 normalize.method="constant",
                 pmcorrect.method="pmonly",
                 summary.method="avgdiff")
  head(eset@assayData$exprs)
  eset@assayData$exprs[1:10,1:2]
  log2(eset@assayData$exprs[1:10,1:2])
# the espresso function doesn't do log transformation.
#don't forget to do it youorself... .check the ranges 
#range(exprs(eset2))
exprs(eset)=log2(exprs(eset))
  head(eset@assayData$exprs)
  eset@assayData$exprs[1:10,1:2]

exp=as.matrix(eset@assayData$exprs)
  head(exp)

#Saving

#write.exprs(eset, "preprocessed.RDATA")
write.csv(exp, "exp.csv")

#Convert into Tidy Data ....Takes 5 Minutes
exp_Tidy <-  as_tibble(exp) %>% gather(samples , expression)%>% merge(phenotable ,  by = 'samples' , all.x = T)
  names(exp_Tidy)
  head(exp_Tidy)

exp_Tidy_gg <- ggplot2::ggplot(exp_Tidy)

#Histogram
exp_Tidy_gg + 
  geom_density(aes(x = log(expression) , 
                   group = samples , 
                   fill = treated) , 
               adjust = 1.2 , alpha = .4) + 
  theme_classic() +
  xlab('log2 Expression') +
  ylab('Density')

#Boxplot
exp_Tidy_gg + 
  geom_boxplot(aes(x = samples, 
                   y = log(expression) , 
                   fill = samples)) + 
  theme(axis.text.x = element_text(angle = 90))

#Prinipal Component Analysis
PCA.exp <- prcomp(t(exp),
              center = T,
              scale = T)
  dim(PCA.exp$x)

#Plot PCA
autoplot(PCA.exp,
         data = phenotable,
         colour = 'sample.type2',
         size = 2.5,
         #frame = T,
         #frame.type = 'norm'
)

exp_df=cbind(probe_id=rownames(exp),exp)
#### mapping the probe ids into gene symbole
hugene10sttranscriptcluster()
mapper = hugene10sttranscriptclusterSYMBOL
#mapper <- read_delim("mapper.txt","\t", escape_double = FALSE, trim_ws = TRUE)
#names(mapper)
map.df = as.data.frame(mapper)
  head(map.df)
# merge the two data frames to have the symbole annotation in the data object
exp_df=merge(exp_df,map.df,by="probe_id",all.x=T)
  names(exp_df)
# do i need the probe id again?  no, then let's drop it
exp_df=exp_df[,-1]
# remove nulls : some probes were not mapped to any gene symbol, to remove the rows has no gene symbol
exp_df=exp_df[ ! is.na(exp_df$symbol),]
# check duplciation of of gene symbols?  
x=duplicated(exp_df$symbol)  
  sum(x)

### yes .. why ? probesets?  solutions : aggregation
exp.data=exp_df[-dim(exp_df)[2]]
exp.data=apply(exp.data,2, as.numeric)

####Aggregate
exp.data.agg= aggregate(exp.data, list(exp_df$symbol),FUN=mean)
  names(exp.data.agg)
rownames(exp.data.agg)=exp.data.agg$Group.1
exp.data.agg=exp.data.agg[- 1]

# save the final object in a RDATA  object
#phenotable$treated= unlist(sapply(strsplit(phenotable$Sample.Type, "_"),c)[2,])
exp=exp.data.agg
save(exp,phenotable, file="GSE77515.RDATA")
###################   do the PCA and exploratory analysis  ############
#Convert into Tidy Data ....Takes 5 Minutes
exp_Tidy <-  as_tibble(exp) %>% gather(samples , expression)%>% merge(phenotable ,  by = 'samples' , all.x = T)
  names(exp_Tidy)
  head(exp_Tidy)

exp_Tidy_gg <- ggplot2::ggplot(exp_Tidy)

#Histogram
exp_Tidy_gg + 
  geom_density(aes(x = log(expression) , 
                   group = samples , 
                   fill = treated) , 
               adjust = 1.2 , alpha = .4) + 
  theme_classic() +
  xlab('log2 Expression') +
  ylab('Density')

#Boxplot
exp_Tidy_gg + 
  geom_boxplot(aes(x = samples, 
                   y = log(expression) , 
                   fill = samples)) + 
  theme(axis.text.x = element_text(angle = 90))

#Prinipal Component Analysis
PCA.exp <- prcomp(t(exp),
                  center = T,
                  scale = T)
dim(PCA.exp$x)

#Plot PCA
autoplot(PCA.exp,
         data = phenotable,
         colour = 'sample.type2',
         size = 2.5,
         #frame = T,
         #frame.type = 'norm'
)


###################   DO differential expression analysis ############
groups=unique(phenotable$sample.type2)
group1=groups[1]
group2=groups[2]
group3=groups[3]

group1.columns=  phenotable[ phenotable$sample.type2 ==group1,]$samples
group2.columns=  phenotable[ phenotable$sample.type2 ==group2,]$samples
group3.columns=  phenotable[ phenotable$sample.type2 ==group3,]$samples

#DEG2 Cis24 and Ctrl
exp=exp[,c(group1.columns,group3.columns)]
## calculating LFC,1 is rows
lfc=apply(exp,1, function(x)  mean(x[group3.columns]) -mean(x[group1.columns]))
res=as.data.frame(lfc)

## calcualting p values
f=factor( c( rep(1, length(group1.columns)) , rep(2, length(group3.columns)) ))
pval=rowttests(as.matrix(exp),f)$p.value

correctPvalueandReturnAll<-function(tt.pval,method)
{
  mt=mt.rawp2adjp(tt.pval,proc=method)
  adjp=mt$adjp[order(mt$index),]
  return(adjp[,2])
}
adj.pval=correctPvalueandReturnAll(pval,"BH")



res=cbind(lfc,pval,adj.pval)
#####  selection criteria for identifying DEGS
res.degs=res[abs(lfc) > log2(2)   & adj.pval <0.05,]  ##### identify DEGs based on both  LFC and the significane level 
res.degs=res.degs[order(res.degs[,3]) ,]
degs.genes= rownames(res.degs)

exp.degs=exp[degs.genes,]
#export them for further analysis in DAVID quote is "" on colnames andrownames
write.table(degs.genes,file = "DEGs_Cis24_Ctrl.txt",row.names = F,col.names = F,quote = F)#to get txt file

#bar Plot 
group1.3.genes <- exp[degs.genes[1:10],c(group1.columns,group3.columns)]
group1.3.genes <- mutate(group1.3.genes, genes = row.names(group1.3.genes))
group1.3.genes <- gather(group1.3.genes,samples,expr,-genes)
group1.3.genes <- merge(group1.3.genes, phenotable, by = 'samples', all.x = T)
ggplot(group1.3.genes , aes(x = genes , y = expr , fill = sample.type2)) + geom_bar(stat="identity" , position = 'dodge') + 
  theme(axis.text.x = element_text(angle = 90))
#Volcano Plot
library(EnhancedVolcano)
res1 <- as.data.frame(res)
EnhancedVolcano(res1,
                lab = row.names(res1),
                x = 'lfc',
                y =  'adj.pval',
                title = 'GSE77515-DGEs',
                #selectLab = c('VCAM1','KCTD12','ADAM12'),
                ylab = bquote(~-log[10]~ 'adj. P'),
                xlab = bquote(~Log[2]~ 'fold change'),
                boxedLabels = TRUE,
                drawConnectors = TRUE,
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 1.0,
                labSize = 6.0,
                colAlpha = 1,
                cutoffLineType = 'blank',
                cutoffLineCol = 'black',
                cutoffLineWidth = 0.8,
                gridlines.major = FALSE,
                gridlines.minor = FALSE) #+ 
  #coord_flip()


#try at home as an assignment:
#https://bioconductor.org/packages/devel/bioc/vignettes/EnhancedVolcano/inst/doc/EnhancedVolcano.html
#https://biocorecrg.github.io/CRG_RIntroduction/volcano-plots.html
#https://samdsblog.netlify.app/post/visualizing-volcano-plots-in-r/



### creating a heatmap for the top 100 DEG genes ## use always ComplexHeatMap
# 1- get the expression profiles of the top 100 degs only
n_genes = 100
topDEGs <- exp.degs[1:n_genes , ]

#Color Functions:
##Main Color Function
col_fun <- colorRamp2(seq(min(topDEGs), max(topDEGs), length = 6), rev(heat.colors(6)))

##Sample.type function
sample.type_fun <- list(sample.type = c(group1 = "#FAB50C" , 
                                        #group2 = "#C58D03"  , 
                                        group3 = "#513C07"))
#Annotations
sample.type <- c(rep(group1, length(group1.columns)), rep(group3, length(group1.columns)))
bottom_annotaion = HeatmapAnnotation(sample.type = sample.type , col = sample.type_fun)

ComplexHeatmap::Heatmap(as.matrix(topDEGs),
                        col = col_fun,
                        name = "Abs. Expr.",
                        na_col = "black",
                        #rect_gp = gpar(col = "white", lwd = 2),
                        row_dend_side = "right",
                        column_dend_side = "bottom",
                        column_dend_height = unit(2, "cm"),
                        row_dend_width = unit(2, "cm"),
                        row_names_gp = gpar(fontsize = 4.2),
                        column_names_gp = gpar(fontsize = 8),
                        show_column_dend = TRUE,
                        show_row_dend = TRUE,
                        #top_annotation = top_annotation,
                        #bottom_annotation = bottom_annotaion,
                        row_dend_gp = gpar(col='#A31002'),
                        column_dend_gp = gpar(col='#A31002')
)



exp.genes <- as.character(unlist(row.names(exp)))

degs.genes[stringr::str_detect(degs.genes, "IL")]
exp.genes[stringr::str_detect(exp.genes, "HOXC")]
res["FAM193A",]

#Prinipal Component Analysis
PCA.exp <- prcomp(t(topDEGs),
                  center = T,
                  scale = T)
  dim(PCA.exp$x)

#Plot PCA
autoplot(PCA.exp,
         data = phenotable[phenotable$samples == group1.columns | phenotable$samples ==  group3.columns,],
         colour = 'sample.type2',
         size = 2.5,
         #frame = T,
         #frame.type = 'norm'
)

### creating  3D PCA

### creating  a schematic figure for study design
