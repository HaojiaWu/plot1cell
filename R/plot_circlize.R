#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Convert coordinates
#'
#' This function converts the Cartesian coordinates to Polar coordinates. 
#' Input data can be the coordinates from tSNE or UMAP. It outputs a matrix with
#' polar coordinates.
#'
#' @param coord_data Cartesian coordinates from tSNE, UMAP, etc.
#' @param zoom value from c(0,1) to adjust the coordinates.
#' @return A matrix with polar coordinates
#' @export
transform_coordinates <- function(
  coord_data, 
  zoom
  ){
  center_data<-coord_data-mean(c(min(coord_data),max(coord_data)))
  max_data<-max(center_data)
  new_data<-center_data*zoom/max_data
  new_data
}

#' Get metadata from a Seurat object
#'
#' This function extracts the metadata from a Seurat object and transforms the
#' UMAP/tSNE coordinates.
#'
#' @param obj SeuratObject
#' @param reductions reductions methods, e.g."umap" or "tsne".
#' @param color Colors assigned to the cell clusters
#' @param coord_scale value from c(0,1) to adjust the UMAP/tSNE coordinates.
#' @return A metadata dataframe
#' @export
get_metadata <- function(
  obj, 
  reductions = "umap", 
  coord_scale = 0.8, 
  color
  ){
  metadata<-obj@meta.data
  metadata$Cluster<-obj@active.ident
  metadata$dim1<-as.numeric(obj[[reductions]]@cell.embeddings[,1])
  metadata$dim2<-as.numeric(obj[[reductions]]@cell.embeddings[,2])
  metadata$x<-transform_coordinates(metadata$dim1, zoom = coord_scale)
  metadata$y<-transform_coordinates(metadata$dim2, zoom = coord_scale)
  color_df<-data.frame(Cluster=levels(obj), Colors=color)
  cellnames<-rownames(metadata)
  metadata$cells<-rownames(metadata)
  metadata<-merge(metadata, color_df, by='Cluster')
  rownames(metadata)<-metadata$cells
  metadata<-metadata[cellnames,]
  metadata
}

#' Make count matrix for the selected markers
#'
#' This function labels the cells based their expression levels of the selected 
#' marker genes.
#'
#' @param obj SeuratObject
#' @param features Selected marker genes
#' @return A dataframe with cells labeled by marker genes
#' @export
mk_marker_ct <- function(
  obj, 
  features
  ){
  dat <- Seurat::FetchData(obj, vars = features)
  ori_names <- rownames(dat)
  zero_ct <- dat[rowSums(dat)==0,]
  non_zero <- dat[rowSums(dat)!=0,]
  max_genes <- colnames(non_zero)[max.col(non_zero,ties.method="first")]
  non_zero <- data.frame(cells=rownames(non_zero), genes=max_genes)
  zero_ct <- data.frame(cells=rownames(zero_ct), genes='No_expr')
  all_cells <- rbind(non_zero, zero_ct)
  rownames(all_cells) <- all_cells$cells
  all_cells <- all_cells[ori_names,]
  all_cells
}

#' Create a dataframe for color mapping
#'
#' This function assigns a color for each value in a vector
#'
#' @param group Group to be assigned color
#' @return A dataframe with colors assigned to groups
#' @export
mk_color_table <- function(group){
  n=length(group)
  colors=scales::hue_pal()(n)
  color_table <- data.frame(group, colors)
  color_table
}

#' Order the cells from each cluster
#'
#' This function orders the cells from each cluster by giving a value from
#' 1 to max
#' @param dat Data input. 
#' @return An vector with ordered cells
#' @export
cell_order <- function(dat){
  celltypes <- names(table(dat$Cluster))
  new_dat <- list()
  for (i in 1:length(celltypes)){
    dat$Cluster<-as.character(dat$Cluster)
    dat1<-dat[dat$Cluster==celltypes[i],]
    dat1$x_polar<-1:nrow(dat1)
    new_dat[[i]]<-dat1
  }
  new_dat<-do.call('rbind', new_dat)
  new_dat
}

#' Create a segment for each element in a group
#'
#' This function creates a segment for each element within a group
#' @param dat Data input. 
#' @param group The group name
#' @return An vector with ordered cells
#' @export
get_segment <- function(
  dat, 
  group
  ){
  dat<-dat[order(dat[,group],decreasing = F), ]
  rownames(dat)<-1:nrow(dat)
  dat<-dat[!duplicated(dat[,group]),]
  dat_seg<-as.integer(rownames(dat))
  dat_seg
}

#' Prepare circlize data for plotting
#'
#' This function creates a segment for each element within a group
#' @param seu_obj Seurat object 
#' @param scale Scale factor to zoom in our zoom out the tSNE/UMAP proportionally
#' @return A data frame for plotting
#' @export
prepare_circlize_data <- function(
  seu_obj, 
  scale =0.8
  ){
  celltypes<-levels(seu_obj)
  cell_colors <- scales::hue_pal()(length(celltypes))
  data_plot <- get_metadata(seu_obj, color = cell_colors, coord_scale = scale)
  data_plot <- cell_order(data_plot)
  data_plot$x_polar2 <- log10(data_plot$x_polar)
  data_plot
}

#' Generate a circlize plot outside the tSNE/UMAP
#'
#' This function generates a circlize plot outside the tSNE/UMAP
#'
#' @param data_plot Data frame prepared by the prepare_circlize_data function
#' @param do.label Whether to label the clusters
#' @param contour.levels Which contour line to be drawn on the plot. Value: 0-1
#' @param bg.color Canvas background color
#' @param pt.size Point size of the graph
#' @param kde2d.n Number of grid points in each direction. A kde2d parameter
#' @param contour.nlevels Total number of levels in contour
#' @return Return a circlize plot
#' @export
plot_circlize <- function(
  data_plot,
  do.label = T,
  contour.levels = c(0.2,0.3),
  pt.size = 0.5,
  kde2d.n = 1000,
  contour.nlevels = 100,
  bg.color='#F9F2E4'
  ) {
  data_plot %>%
    dplyr::group_by(Cluster) %>%
    summarize(x = median(x = x),y = median(x = y)) -> centers
  z <- MASS::kde2d(data_plot$x, data_plot$y, n=kde2d.n)
  celltypes<-names(table(data_plot$Cluster))
  cell_colors <- scales::hue_pal()(length(celltypes))
  circos.clear()
  par(bg = bg.color)
  circos.par(cell.padding=c(0,0,0,0), track.margin=c(0.01,0),"track.height" = 0.01, gap.degree =c(rep(2, (length(celltypes)-1)),12))
  circos.initialize(sectors =  data_plot$Cluster, x = data_plot$x_polar2)
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, bg.border=NA,panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter,
                CELL_META$cell.ylim[2]+ mm_y(4),
                CELL_META$sector.index,
                cex=0.5, col = 'black', facing = "bending.inside", niceFacing = T)
    circos.axis(labels.cex = 0.3, col = 'black', labels.col =  'black')
  })
  for(i in 1:length(celltypes)){
    dd<-data_plot[data_plot$Cluster==celltypes[i],]
    circos.segments(x0 = min(dd$x_polar2), y0 = 0, x1 = max(dd$x_polar2), y1 = 0, col = cell_colors[i],  lwd=3, sector.index = celltypes[i])
  }
  text(x = 1, y=0.1, labels = "Cluster", cex = 0.4, col = 'black',srt=-90)
  points(data_plot$x,data_plot$y, pch = 19, col = alpha(data_plot$Colors,0.2), cex = pt.size);
  contour(z, drawlabels=F, nlevels= 100, levels = contour.levels,col = '#ae9c76', add=TRUE)
  if(do.label){
  text(centers$x,centers$y, labels=centers$Cluster, cex = 0.8, col = 'black')
  }
}

#' Add tracts to the circlize plot
#'
#' This function allows users to add more tracks into the circlize plot
#' @param data_plot Data for circlize plot 
#' @param group The group for showing on the new track
#' @param colors Color palette to color the group
#' @return A new circlize track adding to the current circlize plot
#' @export
add_tract <- function(data_plot, group, colors = NULL){
  circos.track(data_plot$Cluster, data_plot$x_polar2, y=data_plot$dim2, bg.border=NA)
  celltypes<-names(table(data_plot$Cluster))
  group_names<-names(table(data_plot[,group]))
  if(is.null(colors)){
    col_group = scales::hue_pal()(length(group_names))
  } else {
    col_group = colors
  }
  for(i in 1:length(celltypes)) {
    data_plot_cl<-data_plot[data_plot$Cluster==celltypes[i],]
    dat_seg<-get_segment(data_plot_cl, group = group)
    dat_seg2<-c(dat_seg[-1]-1, nrow(data_plot_cl))
    scale_factor<-max(data_plot_cl$x_polar2)/nrow(data_plot_cl)
    dat_seg<-scale_factor*dat_seg
    dat_seg2<-scale_factor*dat_seg2
    circos.segments(x0 = dat_seg, y0 = 0, x1 = dat_seg2, y1 = 0, col = col_group, sector.index = celltypes[i], lwd=3)
  }
}
