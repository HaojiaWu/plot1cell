#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Generate the files required by CellPhoneDB
#'
#' This function can generate the files used as input into CellPhoneDB analysis. 
#' Two files will be created: normalized count matrix with the designated group
#' and the metadata to descript each cell.
#'
#' @param seu_obj A complete Seurat object
#' @param group The group selected for analysis
#' @return Output two txt files for CellPhoneDB
#' @export
creat_cellphonedb_file <- function(seu_obj, group){
  counts <- as.data.frame(
    as.matrix(
      seu_obj@assays$RNA@data)
  )
  mouse_genes<-rownames(counts)
  human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                   values = mouse_genes , mart = mouse, attributesL = c("hgnc_symbol",'ensembl_gene_id'), 
                   martL = human, uniqueRows=F)
  genesV2$Gene.stable.ID<-as.character(genesV2$Gene.stable.ID)
  genesV2$MGI.symbol<-as.character(genesV2$MGI.symbol)
  genesV2<-genesV2[!duplicated(genesV2$Gene.stable.ID),]
  genesV2<-genesV2[!duplicated(genesV2$MGI.symbol),]
  counts<-counts[genesV2$MGI.symbol,]
  rownames(counts)<-genesV2$Gene.stable.ID
  metadata <- data.frame(Cell = colnames(counts),
                         cell_type = as.character(seu_obj@active.ident)
  )
  metadata$Cell<-gsub('\\-','\\.',metadata$Cell)
  counts<-na.omit(counts)
  counts$Gene<-rownames(counts)
  counts<-counts[,c(ncol(counts), 1:(ncol(counts)-1))]
  colnames(counts)<-gsub('\\-','\\.',colnames(counts))
  
  write.table(counts,
              file =paste0(group, '_counts.txt'),
              quote = F,
              col.names = T,
              row.names = F,
              sep = '\t')
  
  write.table(metadata,
              file = paste0(group, '_metadata.txt'),
              quote = F,
              col.names = T,
              row.names = F,
              sep = '\t')
  
}

#' Generate the files required by pySCENIC
#'
#' This function can generate the file used as input into pySCENIC analysis. 
#' A loom file will be created.
#'
#' @param seu_obj A complete Seurat object
#' @param celltypes The cell types being selected for analysis
#' @param min.gene Cutoff to filter the low expression genes
#' @return A loom file
#' @export
create_pyscenic_file <- function (seu_obj, celltypes, min.gene = 1) {
  seu_sub<-subset(seu_obj, idents=celltypes)
  sub_count<-seu_sub@assays$RNA@counts
  sub_count<-sub_count[rowSums(sub_count) > min.gene,]
  seu_sub<-subset(seu_sub, features=rownames(sub_count))
  SaveH5Seurat(seu_sub, filename = "seu_sub.h5Seurat")
  Convert("seu_sub.h5Seurat", dest = "h5ad")
  cat("Please run the following lines in python.\n
      import scanpy as sc \n
      import numpy as np \n
      import loompy as lp \n
      adata = sc.read(\"seu_sub.h5ad\") \n
      row_attrs = { \n
      \"Gene\": np.array(adata.var_names), \n
      } \n
      col_attrs = { \n
     \"CellID\": np.array(adata.obs_names), \n
     \"nGene\": np.array(np.sum(adata.X.transpose()>0, axis=0)).flatten(), \n
      \"nUMI\": np.array(np.sum(adata.X.transpose(),axis=0)).flatten(), \n
      } \n
      lp.create(\"seu_sub_filtered.loom\",adata.X.transpose(),row_attrs, \n
          col_attrs) \n
  ")
}


#' Compute and plot correlations between two datasets
#'
#' This function can compute the correlations on the samples from two count matrix.
#' It will use Pearson method as default.
#'
#' @param data1 The first data frame with genes in rows and samples in columns
#' @param data2 The second data frame with genes in rows and samples in columns
#' @param ngenes Number of top variable genes used for the computation
#' @param method.use Default is "pearson". Other methods include kendall" or "spearman". See ?cor.
#' @param color.use Color palette for the heatmap plot.
#' @return A heatmap plot
#' @export
run_correlation<-function(
  data1, 
  data2, 
  ngenes=2000, 
  method.use='pearson', 
  color.use=NULL
  ){
  genes1<-rownames(data1)
  genes2<-rownames(data2)
  comm_genes<-intersect(genes1, genes1)
  data1<-data1[comm_genes,]
  data2<-data2[comm_genes,]
  merged<-cbind(data1, data2)
  datExpr<-as.matrix(merged)
  Pvars <- rowVars(datExpr)
  select <- order(Pvars, decreasing = TRUE)[seq_len(min(ngenes, length(Pvars)))]
  clust_select<-datExpr[select,]
  pearsonplot<-cor(clust_select, use="complete.obs", method=method.use)
  cell1<-length(colnames(data1))
  cell2<-length(colnames(data2))
  pearsonplot<-pearsonplot[1:cell1, (1+cell1):ncol(pearsonplot)]
  pearsonplot<-round (pearsonplot, 2)
  pearsonplot<-as.matrix(pearsonplot)
  colors = color.use
  if(is.null(colors)){
    colors<-colorRampPalette(c("lemonchiffon1","white","darkred"))(255)
  }
  pearsonplot[pearsonplot<0]<-0
  pheatmap(pearsonplot, col=colors, display_numbers = F, cluster_cols = F, cluster_rows = F, fontsize = 8)
}

#' A function to convert ensembl ID to gene name
#'
#' This function is for replace the ensembl ID with actual gene names
#'
#' @param count Count data input
#' @param species value can be "mouse" or "human"
#' @return A new count data with external gene names in the rows
#' @export

convert_geneid <- function(
  count, 
  species
  ){
  if (species == "mouse") {
    mart = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  }
  else if (species == "human") {
    mart = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
  }
  else {
    print("Please specify a species!")
  }
  results <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                   mart = mart)
  results$external_gene_name[results$external_gene_name==""]<-NA
  results<-na.omit(results)
  all1 <- data.frame(as.matrix(count))
  all1$gene <- rownames(count)
  annotate2 <- results[which(results$ensembl_gene_id %in% all1$gene), 
  ]
  all1 <- all1[annotate2$ensembl_gene_id, ]
  all1 <- cbind(all1, annotate2)
  all1 <- all1[!duplicated(all1$external_gene_name), ]
  rownames(all1) <- all1$external_gene_name
  new_count <- all1[, -c((ncol(all1) - 2):ncol(all1))]
  new_count <- new_count[-nrow(new_count), ]
  new_count
}

#' A function to process the h5 file from CellBender and generate QC plots
#'
#' This function works only with the h5 file output from CellBender. It takes the 
#' count table, clusters the cells, removes the doublets (DoubletFinder), and 
#' recluster the cells, and finally output a Seurat object
#'
#' @param cellbender_h5 The h5 file from CellBender output
#' @param sampleID ID given to the run
#' @param out_dir Output directory
#' @param species It can be "human" or "mouse"
#' @param type It can be "cell" or "nucleus". If your data is from snRNA-seq, use "type = nucleus". Otherwise, use "type = cell".
#' @return QC plots and seurat objects before and after QC
#' @export

data_processing <- function(
  cellbender_h5, 
  sampleID, 
  out_dir, 
  species,
  type
  ){
  new_df<-Read10X_h5(cellbender_h5)
  colnames(new_df)<-paste(sampleID, colnames(new_df), sep="_")
  ##### 1.initial QC and clustering ####
  Sample4<- CreateSeuratObject(counts = new_df, project = sampleID, min.cells = 3, min.features = 200)
  metadata<-Sample4@meta.data
  png(paste0(out_dir,sampleID, "_counts_histogram_before_clean.png"), units="in", width=6, height=4, res=300)
  hist(
    metadata$nFeature_RNA,
    breaks = 1000
  )
  dev.off()
  if(species=="human"){
    Sample4[["percent.mt"]] <- PercentageFeatureSet(Sample4, pattern = "^MT-")
  } else {
    Sample4[["percent.mt"]] <- PercentageFeatureSet(Sample4, pattern = "^mt-")
  }
  png(paste0(out_dir,sampleID, "_counts_QC_before_clean.png"), units="in", width=8, height=4, res=300)
  print(VlnPlot(Sample4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
  dev.off()
  if(type=="nucleus"){
    Sample4<- subset(Sample4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 5)
  } else {
    Sample4<- subset(Sample4, subset = nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 50)
  }
  Sample4<- NormalizeData(Sample4)
  Sample4<- FindVariableFeatures(Sample4, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(Sample4)
  Sample4<- ScaleData(Sample4, features = all.genes)
  Sample4<- RunPCA(Sample4, features = VariableFeatures(object = Sample4), npcs = 50)
  png(paste0(out_dir,sampleID, "_PCA_elbowplot_before_clean.png"), units="in", width=6, height=4, res=300)
  print(ElbowPlot(Sample4, ndims = 50))
  dev.off()
  Sample4<- FindNeighbors(Sample4, dims = 1:30)
  Sample4<- FindClusters(Sample4, resolution = 0.5)
  Sample4<- RunUMAP(Sample4, dims = 1:30)
  Sample4<- CreateSeuratObject(counts = Sample4@assays$RNA@counts, project = sampleID, min.cells = 3, min.features = 200)
  metadata<-Sample4@meta.data
  png(paste0(out_dir,sampleID, "_counts_histogram_before_clean.png"), units="in", width=6, height=4, res=300)
  hist(
    metadata$nFeature_RNA,
    breaks = 1000
  )
  dev.off()
  if(species=="human"){
    Sample4[["percent.mt"]] <- PercentageFeatureSet(Sample4, pattern = "^MT-")
  } else {
    Sample4[["percent.mt"]] <- PercentageFeatureSet(Sample4, pattern = "^mt-")
  }
  png(paste0(out_dir,sampleID, "_counts_QC_before_clean.png"), units="in", width=8, height=4, res=300)
  print(VlnPlot(Sample4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
  dev.off()
  Sample4<- NormalizeData(Sample4)
  Sample4<- FindVariableFeatures(Sample4, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(Sample4)
  Sample4<- ScaleData(Sample4, features = all.genes)
  Sample4<- RunPCA(Sample4, features = VariableFeatures(object = Sample4), npcs = 50)
  png(paste0(out_dir,sampleID, "_PCA_elbowplot_before_clean.png"), units="in", width=6, height=4, res=300)
  print(ElbowPlot(Sample4, ndims = 50))
  dev.off()
  Sample4<- FindNeighbors(Sample4, dims = 1:30)
  Sample4<- FindClusters(Sample4, resolution = 0.5)
  Sample4<- RunUMAP(Sample4, dims = 1:30)
  Sample4<-RunTSNE(Sample4, dims = 1:30, perplexity=50)
  png(paste0(out_dir,sampleID, "_UMAP_before_clean.png"), units="in", width=6, height=4, res=300)
  print(DimPlot(Sample4, reduction = "umap"))
  dev.off()
  png(paste0(out_dir,sampleID, "_tSNE_before_clean.png"), units="in", width=6, height=4, res=300)
  print(DimPlot(Sample4, reduction = "tsne"))
  dev.off()
  genes_select<-c('Nphs1','Slc5a12',"Slc7a12",'Slc12a1',"Slc12a3","Emcn",'Fhl2',"Tnc","Aqp2","Ptprc","Slc26a4",'Kit')
  if(species=="human"){
  genes_select<-toupper(genes_select)
  }
  png(paste0(out_dir,sampleID, "_markers_before_clean.png"), units="in", width=15, height=15, res=300)
  print(FeaturePlot(Sample4, features = genes_select, reduction = 'tsne'))
  dev.off()
  png(paste0(out_dir,sampleID, "_nGene_tSNE_before_clean.png"), units="in", width=6, height=4, res=300)
  print(FeaturePlot(Sample4, features = "nFeature_RNA", reduction = 'tsne'))
  dev.off()
  saveRDS(Sample4, file = paste0(out_dir,sampleID,'_seurat.rds'))
  Sys.time()
  print(paste(sampleID,"initial clustering done!", sep = ":"))
  
  ### 2.doubleFinder to remove doublets ###
  doublet_rate<-0.018+8.4e-06*ncol(Sample4)  ##This function is based on the data from the Demuxlet paper##
  ndoublets<-doublet_rate*ncol(Sample4)
  Sample4<-doubletFinder_v3(Sample4,PCs = 1:30, nExp = ndoublets, pK = 0.09)
  png(paste0(out_dir,sampleID, "_doublets_tSNE.png"), units="in", width=6, height=4, res=300)
  print(TSNEPlot(Sample4, group.by=paste0('DF.classifications_0.25_0.09_', ndoublets)))
  dev.off()
  Sample4<-SetIdent(Sample4, value = paste0('DF.classifications_0.25_0.09_', ndoublets))
  Sample4_clean<-subset(Sample4, idents = 'Singlet')
  Sys.time()
  print(paste(sampleID,"doubletFinder done!", sep = ":"))
  
  ### 3.recluster the singlets ####
  Sample4_clean <- CreateSeuratObject(counts = Sample4_clean@assays$RNA@counts, project = sampleID, min.cells = 3, min.features = 200)
  Sample4_clean
  metadata<-Sample4_clean@meta.data
  png(paste0(out_dir,sampleID, "_counts_histogram_after_clean.png"), units="in", width=6, height=4, res=300)
  hist(
    metadata$nFeature_RNA,
    breaks = 1000
  )
  dev.off()
  Sample4_clean[["percent.mt"]] <- PercentageFeatureSet(Sample4_clean, pattern = "^mt-")
  png(paste0(out_dir,sampleID, "_counts_QC_after_clean.png"), units="in", width=8, height=4, res=300)
  print(VlnPlot(Sample4_clean, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1))
  dev.off()
  Sample4_clean <- NormalizeData(Sample4_clean)
  Sample4_clean <- FindVariableFeatures(Sample4_clean, selection.method = "vst", nfeatures = 2000)
  all.genes <- rownames(Sample4_clean)
  Sample4_clean <- ScaleData(Sample4_clean, features = all.genes)
  Sample4_clean <- RunPCA(Sample4_clean, features = VariableFeatures(object = Sample4_clean), npcs = 50)
  png(paste0(out_dir,sampleID, "_PCA_elbowplot_after_clean.png"), units="in", width=6, height=4, res=300)
  print(ElbowPlot(Sample4_clean, ndims = 50))
  dev.off()
  Sample4_clean <- FindNeighbors(Sample4_clean, dims = 1:25)
  Sample4_clean <- FindClusters(Sample4_clean, resolution = 0.6)
  Sample4_clean <- RunUMAP(Sample4_clean, dims = 1:25, min.dist = 0.6, spread = 3)
  print(paste(sampleID,"final clustering done!", sep = ":"))
  png(paste0(out_dir,sampleID, "_UMAP_after_clean.png"), units="in", width=6, height=4, res=300)
  print(DimPlot(Sample4_clean, reduction = "umap", label = T))
  dev.off()
  png(paste0(out_dir,sampleID, "_markers_after_clean.png"), units="in", width=15, height=15, res=300)
  print(FeaturePlot(Sample4_clean, features = genes_select))
  dev.off()
  png(paste0(out_dir,sampleID, "_nGene_UMAP_after_clean.png"), units="in", width=6, height=4, res=300)
  print(FeaturePlot(Sample4_clean, features = 'nFeature_RNA'))
  dev.off()
  saveRDS(Sample4_clean, file=paste0(out_dir,sampleID,'_seurat_clean.rds'))
  print(paste(sampleID,"all done!", sep = ":"))
  Sys.time()
}


#' A function to extract gene counts for ploting
#'
#' This function is a modified Seurat::FetchData function to extract gene 
#' counts and the associated meta data for ploting. It returns a dataframe
#' with the requested information from the Seurat object.
#'
#' @param seu_obj A finished Seurat Object with cell type annotation in the active.ident slot
#' @param features Gene names to extract expression data
#' @param cell.types The cell types to be inspected. By default, it will incorporate all cell types.
#' @param data.type The data slot to be accessed. By default, the "data" slot will be used.
#' @param meta.groups The colnames in the meta.data slot you want to include.
#' @return A data frame with the requested info.
#' @export
#' 
extract_gene_count <- function(
  seu_obj, 
  features, 
  cell.types=NULL, 
  data.type="data", 
  meta.groups=NULL
 ){
  if(is.null(cell.types)){
    cell.types=levels(seu_obj)
  }
  seu_obj@meta.data$celltype<-as.character(seu_obj@active.ident)
  if(is.null(meta.groups)){
    meta.groups=colnames(seu_obj@meta.data)
  }
  if(!is.null(cell.types)){
    new_seu<-subset(seu_obj, idents=cell.types)
  }
  feature_count<-Seurat::FetchData(new_seu, slot = data.type, vars = c(features,meta.groups,"celltype"))
  umap_data<-data.frame(new_seu[["umap"]]@cell.embeddings)
  feature_count$UMAP1<-umap_data$UMAP_1
  feature_count$UMAP2<-umap_data$UMAP_2
  feature_count
}


#' A function to make gene name first letter capital
#'
#' The function is modified from this thread: https://stackoverflow.com/questions/18509527/first-letter-to-upper-case/18509816
#'
#' @param gene Gene name
#' @export
#' 

firstup <- function(
  gene
  ){
  x <- tolower(gene)
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#' A function to change the strip background color in ggplot
#' @param ggplt_obj A ggplot object
#' @param type Strip on the "top" or "right" side only or "both" sides
#' @param strip.color A color vector
#' @export
#' 
change_strip_background <- function(
  ggplt_obj, 
  type = "top",
  strip.color=NULL
  ){
  g <- ggplot_gtable(ggplot_build(ggplt_obj))
  if(type == "top"){
    strip_both <- which(grepl('strip-t', g$layout$name))
    fills<-strip.color
    if(is.null(fills)){
    fills<- scales::hue_pal(l=90)(length(strip_both))
    }
  } else if(type=="right"){
    strip_both <- which(grepl('strip-r', g$layout$name))
    fills<-strip.color
    if(is.null(fills)){
      fills<- scales::hue_pal(l=90)(length(strip_both))
    }
  } else {
    strip_t <- which(grepl('strip-t', g$layout$name))
    strip_r <- which(grepl('strip-r', g$layout$name))
    strip_both<-c(strip_t, strip_r)
    fills<-strip.color
    if(is.null(fills)){
      fills<- c(scales::hue_pal(l=90)(length(strip_t)),scales::hue_pal(l=90)(length(strip_r)))
    }
  }
  k <- 1
  for (i in strip_both) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
    k <- k+1
  }
  g
}

#' A function to generate an example dataset for demo
#' @export
Install.example<-function(){
  getGEOSuppFiles(GEO = "GSE139107")
  setwd("GSE139107/")
  all_files<-list.files(pattern = "dge")
  all_ct<-list()
  for(i in 1:length(all_files)){
    ct_data<-read.delim(all_files[i])
    ct_data<-Matrix(as.matrix(ct_data), sparse = T)
    all_ct[[i]]<-ct_data
    print(all_files[i])
  }
  all.count<-RowMergeSparseMatrices(all_ct[[1]], all_ct[-1])
  meta.data<-read.delim('GSE139107_MouseIRI.metadata.txt.gz')
  all.count<-all.count[,rownames(meta.data)]
  iri <- CreateSeuratObject(counts = all.count, min.cells = 0, min.features = 0, meta.data = meta.data)
  iri.list <- SplitObject(iri, split.by = "orig.ident")
  iri.list <- lapply(X = iri.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
  })
  features <- SelectIntegrationFeatures(object.list = iri.list)
  iri.list <- lapply(X = iri.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
  })
  anchors <- FindIntegrationAnchors(object.list = iri.list, reference = c(1, 2), reduction = "rpca",
                                    dims = 1:50)
  iri.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
  iri.integrated <- ScaleData(iri.integrated, verbose = FALSE)
  iri.integrated <- RunPCA(iri.integrated, verbose = FALSE)
  iri.integrated <- RunUMAP(iri.integrated, dims = 1:25, min.dist = 0.2)
  iri.integrated<-SetIdent(iri.integrated, value = 'celltype')
  levels(iri.integrated)<-c("PTS1"  ,  "PTS2"  ,  "PTS3", "NewPT1", "NewPT2",
                            "DTL-ATL", "MTAL", "CTAL1" ,  "CTAL2","MD","DCT" ,
                            "DCT-CNT","CNT","PC1","PC2","ICA","ICB","Uro","Pod","PEC",
                            "EC1","EC2","Fib","Per","Mø","Tcell")
  levels(iri.integrated@meta.data$Group)<-c("Control","4hours","12hours","2days","14days","6weeks" )
  DefaultAssay(iri.integrated)<-"RNA"
  setwd("../")
  unlink("GSE139107/",recursive=TRUE)
  iri.integrated
}

#' A function to order genes upregulated from the first to the last column
#' @param df A data frames with genes in row and samples in column
#' @export
order_gene_up<-function(df){
  min.col <- function(m, ...) max.col(-m, ...)
  df$celltype<-colnames(df)[min.col(df,ties.method="first")]
  df$celltype<-factor(df$celltype, levels = names(df))
  df<-df[order(df$celltype),]
  gene.names<-list()
  for (i in names(df)){
    aa<-df[df$celltype==i,]
    aa<-aa[order(aa[,i],decreasing = T),]
    gene.names[[i]]<-rownames(aa)
  }
  gene.names<-as.character(unlist(gene.names))
  return(gene.names)
}

#' A function to order genes downregulated from the first to the last column
#' @param df A data frames with genes in row and samples in column
#' @export
order_gene_down<-function(df){
  df$celltype<-colnames(df)[max.col(df,ties.method="first")]
  df$celltype<-factor(df$celltype, levels = names(df))
  df<-df[order(df$celltype),]
  gene.names<-list()
  for (i in names(df)){
    aa<-df[df$celltype==i,]
    aa<-aa[order(aa[,i],decreasing = F),]
    gene.names[[i]]<-rownames(aa)
  }
  gene.names<-as.character(unlist(gene.names))
  return(gene.names)
}

