#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' UpSet plot visualize the number of unique and shared DEGs across group 
#'
#' This function takes Seurat object as input and visualize the genes that 
#' are unique to a particular group or shared by multiple group. 
#' @param seu_obj A complete Seurat object.
#' @param celltype The cell type to analyze.
#' @param group Group factor in meta data.
#' @param logfc Log fold change to select the DEGs
#' @param min_size Minimal number of observations in an intersection for it to be included
#' @return An UpSet plot
#' @export

complex_upset_plot<-function(
  seu_obj, 
  celltype, 
  group,
  logfc = 0.5,
  min_size = 1
){
  cell1<-subset(seu_obj, idents=celltype)
  cell1<-SetIdent(cell1, value = group)
  group_levels<-levels(seu_obj@meta.data[,group])
  if (is.null(group_levels)){
    seu_obj@meta.data[,group] <-factor(seu_obj@meta.data[,group], levels = names(table(seu_obj@meta.data[,group])))
    group_levels<-levels(seu_obj@meta.data[,group])
  }
  levels(cell1)<-group_levels
  all_markers<-FindAllMarkers(cell1, min.pct = 0.1, logfc.threshold = logfc,verbose = F)
  all_markers1<-all_markers[all_markers$avg_log2FC>0,]
  all_markers2<-all_markers[all_markers$avg_log2FC<0,]
  gene_list1<-list()
  for(i in 1:length(group_levels)){
    cluster_marker <- all_markers1[all_markers1$cluster == group_levels[i],]$gene
    cluster_marker <- data.frame("gene" = cluster_marker)
    cluster_marker$cell1 <- 1
    colnames(cluster_marker)[2] <- group_levels[i]
    gene_list1[[i]] <- cluster_marker
  }
  combined_data1 <- purrr::reduce(gene_list1, full_join)
  combined_data1[is.na(combined_data1)] <- 0
  gene_list2<-list()
  for(i in 1:length(group_levels)){
    cluster_marker <- all_markers2[all_markers2$cluster == group_levels[i],]$gene
    cluster_marker <- data.frame("gene" = cluster_marker)
    cluster_marker$cell1 <- 1
    colnames(cluster_marker)[2] <- group_levels[i]
    gene_list2[[i]] <- cluster_marker
  }
  combined_data2 <- purrr::reduce(gene_list2, full_join)
  combined_data2[is.na(combined_data2)] <- 0
  combined_data1$Direction<-"Upregulated"
  combined_data2$Direction<-"Downregulated"
  combined_data<-rbind(combined_data1, combined_data2)
  all_genes<-combined_data$gene
  gene_count<-data.frame(table(all_markers$gene))
  colnames(gene_count)[1]<-"gene"
  combined_data<-merge(combined_data, gene_count, by='gene')
  combined_data$type<-ifelse(combined_data$Freq==1, "Unique","Shared")
  main_bar_color<-hue_pal()(length(group_levels))
  metadata<-data.frame(set=group_levels)
  metadata$color_col<-metadata$set
  upset(combined_data, group_levels, 
        base_annotations=list(
          "Unique or shared DEG"=intersection_size(
            counts=T,
            mapping=aes(fill=Direction),
            width=0.7,
            alpha=0.4
          ) + scale_fill_manual(values= c("blue", "orange"))+
            theme_void()+
            theme(legend.position = "top", legend.title = element_blank())
        ),
        set_sizes=(
          upset_set_size(
            geom=geom_bar(
              aes(fill=type, x=group),
              width=0.7
            ),
            position='right'
          )+ scale_fill_manual(values= c("hotpink",'green'))+theme_void()+
            theme(axis.line.x = element_line(colour = "black"),
                  axis.ticks.x =element_line(size = 0.5, color="black") ,
                  axis.ticks.length = unit(.05, "cm"),
                  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8))
        ),
        width_ratio=0.1,
        stripes=upset_stripes(
          geom = geom_point(size=0.1),
          mapping=aes(color=color_col),
          data=metadata
        ),
        name = celltype,
        min_size = min_size
  )
}