#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' UpSet plot for the cell markers
#'
#' This function takes the marker gene table as input and visualize the genes that 
#' are unique to a particular cell type or shared by multiple cell types. 
#' @param deg_table The dataframe from Seurat::FindAllMarkers function.
#' @param celltypes The cell types of interest
#' @param textscale Scale of the text font size
#' @param pointsize The size of the dots in the matrix
#' @param nintersects The number of intersections showing in the main bar
#' @param alpha The transparency of the bars
#' @return An UpSet plot
#' @export
plot_upset<-function(
  deg_table, 
  celltypes = NULL, 
  nintersects = NULL,
  textscale = 1.5,
  pointsize = 2,
  alpha = 0.3
  ){
    fullcelltypes<-as.character(unique(deg_table$cluster))
    if(is.null(celltypes)){
      celltypes = fullcelltypes
    }
    gene_list<-list()
    for(i in 1:length(fullcelltypes)){
      cluster_marker <- deg_table[deg_table$cluster == fullcelltypes[i],]$gene
      cluster_marker <- data.frame("gene" = cluster_marker)
      cluster_marker$cell1 <- 1
      colnames(cluster_marker)[2] <- fullcelltypes[i]
      gene_list[[i]] <- cluster_marker
    }
  combined_data <- reduce(gene_list, full_join)
  combined_data[is.na(combined_data)] <- 0
  group.color<-scales::hue_pal()(length(celltypes))
  names(group.color)<-fullcelltypes
  group.color<-group.color[celltypes]
  metadata<-data.frame(sets=celltypes)
  metadata$group<-metadata$sets
  combined_data<-combined_data[,celltypes]
  sum_col<-colSums(combined_data)
  sum_col<-sum_col[order(sum_col)]
  group.color2<-rev(as.character(group.color[names(sum_col)]))
  if(is.null(nintersects)){
    nintersects = 2^ncol(combined_data)-1
  }
  upset(combined_data,
        set.metadata = list(data = metadata, 
                            plots = list(list(type = "matrix_rows", 
                                              column = "group", 
                                              colors = group.color, 
                                              alpha = alpha))),
        sets.bar.color = group.color2, 
        point.size = pointsize,
        nsets = length(celltypes),
        nintersects = nintersects,
        text.scale = textscale, 
        main.bar.color=colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(nintersects),
        matrix.color= "black",
        mainbar.y.label="Unique or shared DEG",
        sets.x.label="Total DEG"
        )
}