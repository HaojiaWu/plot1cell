#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Plot gene expression across groups using ComplexHeatmap
#'
#' This function is for identifying the group-specific genes in a selected celltype
#' and plot the expression of those genes in heatmap.
#'
#' @param seu_obj A complete Seurat object
#' @param celltype Cell types selected for gene plot.
#' @param group The group to show on x axis. One of the column names in meta.data.
#' @param gene_highlight Gene names showing on the rows. Default: all genes
#' @param logfc Fold change to select the genes
#' @param return_marker If TRUE, a list of specific gene will be returned.
#' @param col_fun Heatmap color key.
#' @return A ComplexHeatmap object or/and a gene list
#' @export

complex_heatmap_unique<-function(
  seu_obj,
  celltype,
  group,
  gene_highlight=NULL,
  logfc=0.5,
  return_marker=FALSE,
  col_fun = colorRamp2(c(-2, -1, 0, 1, 2), rev(c("#BF0080", "#CE6EAE", "#dddddd", "#6EAE6E", "#008000")))
){
cell1<-subset(seu_obj, idents=celltype)
cell1<-SetIdent(cell1, value = group)
group_levels<-levels(seu_obj@meta.data[,group])
if (is.null(group_levels)){
  seu_obj@meta.data[,group] <-factor(seu_obj@meta.data[,group], levels = names(table(seu_obj@meta.data[,group])))
  group_levels<-levels(seu_obj@meta.data[,group])
  }
levels(cell1)<-group_levels
cell1_avg<-AverageExpression(cell1, verbose = F, return.seurat = F, assays = "RNA")
cell1_avg<-cell1_avg$RNA
cell1_avg<-data.frame(cell1_avg)
all_markers<-FindAllMarkers(cell1, min.pct = 0.1, logfc.threshold = logfc,verbose = F)
all_markers1<-all_markers[all_markers$avg_log2FC>0,]
all_markers2<-all_markers[all_markers$avg_log2FC<0,]
unique_list1<-list()
for(i in 1:length(group_levels)){
  group1<-all_markers1[all_markers1$cluster==group_levels[i],]
  group2<-all_markers1[all_markers1$cluster!=group_levels[i],]
  group1_unique<-setdiff(group1$gene, intersect(group1$gene, group2$gene))
  unique_list1[[i]]<-group1_unique
}
unique_list2<-list()
for(i in 1:length(group_levels)){
  group1<-all_markers2[all_markers2$cluster==group_levels[i],]
  group2<-all_markers2[all_markers2$cluster!=group_levels[i],]
  group1_unique<-setdiff(group1$gene, intersect(group1$gene, group2$gene))
  unique_list2[[i]]<-group1_unique
}
unique_list<-c(unique_list1,unique_list2)
unique_genes<-as.character(unlist(unique_list))
data_plot<-cell1_avg[unique_genes,]
gene_num<-c()
for(i in 1:length(unique_list)){
  gene_num[i]<-length(unique_list[[i]])
}
unique_list<-unique_list[which(gene_num!=0)]
col_split<-group_levels
col_split<-factor(col_split, levels = group_levels)
gene_groups<-c(paste0(group_levels, "_up"),paste0(group_levels, "_down"))
gene_groups<-gene_groups[which(gene_num!=0)]
row_split<-list()
for(i in 1:length(gene_groups)){
  row_split[[i]]<-rep(gene_groups[i], length(unique_list[[i]]))
}
row_split<-as.character(unlist(row_split))
row_split<-factor(row_split, levels = gene_groups)
term = list()
for(i in 1:length(gene_groups)){
  txt1<-as.character(unlist(strsplit(gene_groups[i], split = "_")))
  term[[i]]<-data.frame(txt=txt1, index=c(1,2))
}
names(term) = gene_groups
data_plot<-t(scale(t(data_plot)))
label1<-gene_highlight
if(is.null(label1)){
  label1=rownames(data_plot)
}
ht_opt$message = FALSE
ht<- Heatmap(data_plot,name = "mat", cluster_rows = F, cluster_columns = F, 
             column_title = NULL,col = col_fun, row_title = NULL,
             cluster_row_slices = FALSE, cluster_column_slices = FALSE,
             column_split  = col_split, row_split = row_split,
             column_gap = unit(0, "mm"), row_gap = unit(0, "mm"), 
             heatmap_legend_param = list(direction = "horizontal",title = celltype)
) +
  rowAnnotation(link = anno_mark(at = match(label1,unique_genes), which = 'row',
                                 labels = label1, 
                                 labels_gp = gpar(fontsize = 10), padding = unit(1, "mm"))) +
  rowAnnotation(wc = anno_word_cloud(row_split, term, add_new_line = TRUE,
                                     value_range = c(1, 2), fontsize_range = c(12, 12)) ) 
ht=draw(ht,heatmap_legend_side = "top")
down_groups<-sum(grepl("down", gene_groups))
gaps<-c(1:length(group_levels), 1:length(group_levels))
gaps<-gaps[which(gene_num!=0)]
for(i in 1:length(gene_groups)){
  if(i<=down_groups){
    decorate_heatmap_body("mat", row_slice = i, column_slice = 1, code = {
      grid.rect(unit(gaps[i]-1, "npc"), unit(1, "npc"), 
                width = 1 * unit(1, "npc"),
                height = 1 * unit(1, "npc"),
                gp = gpar(lwd = 2, lty = 1, fill=NA, col='black'), just = c("left", "top") 
      )
    }
    )
  } else {
    decorate_heatmap_body("mat", row_slice = i, column_slice = 1, code = {
      grid.rect(unit(gaps[i]-1, "npc"), unit(1, "npc"), 
                width = 1 * unit(1, "npc"),
                height = 1 * unit(1, "npc"),
                gp = gpar(lwd = 2, lty = 1, fill=NA, col='blue'), just = c("left", "top") 
      )
    }
    )
  }
}
if(return_marker){
  names(unique_list)<-gene_groups
  return(unique_list)
}
}
