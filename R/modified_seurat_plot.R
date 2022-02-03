#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Modified DotPlot
#'
#' This function is a modified version of the DotPlot function in Seurat
#'
#' @param seu_obj A complete Seurat object
#' @param features A vector of gene names.
#' @return A ggplot object
#' @export
modified_dotplot <- function(
  seu_obj, 
  features, 
  col_palette = NULL,
  scale.by='radius'
){
  levels(seu_obj) <- rev(levels(seu_obj))
  dataplot<-DotPlot(seu_obj, features = features)
  dataplot<-dataplot$data
  dataplot$avg.exp<-scale(dataplot$avg.exp)
  colnames(dataplot)[1:2]<-c('Avg.Exp', 'Pct.Exp')
  if(is.null(col_palette)){
    col_palette = colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
  }
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  
  if(max(dataplot$Pct.Exp)>=20) {
    dotplot<-ggplot(dataplot, aes(x = features.plot, y = id, fill=Avg.Exp)) + 
      geom_tile(fill="white", color="white") +
      geom_point(aes( size =Pct.Exp), shape=21, color='grey80')  +  
      scale_fill_gradientn(colours  =  col_palette)+
      scale_size(range = c(0, 10))+
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.line = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.key = element_rect(colour = NA, fill = NA),
            axis.text = element_text(size = 12),axis.title=element_text(size=8),legend.text=element_text(size=8), 
            legend.title = element_text(size = 10),legend.position="right", legend.margin=margin())+ylab("")+xlab("")+
      guides(size = guide_legend(override.aes = list(color='black')))
  } else {
    dotplot<-ggplot(dataplot, aes(x = features.plot, y = id, fill=Avg.Exp)) +  
      geom_tile(fill="white", color="white") +
      geom_point(aes( size =Pct.Exp), shape=21, color='grey80')  +  
      scale_fill_gradientn(colours  = col_palette)+
      scale.func(range = c(0, 10), limits = c(0, 20)) +
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.line = element_blank(),axis.text.x = element_text(angle = 45, hjust = 1), 
            legend.key = element_rect(colour = NA, fill = NA),
            axis.text = element_text(size = 12),axis.title=element_text(size=8),legend.text=element_text(size=8), 
            legend.title = element_text(size = 10),legend.position="right", legend.margin=margin())+ylab("")+xlab("")+
      guides(size = guide_legend(override.aes = list(color='black')))
  }
  dotplot
}

#' Modified DimPlot
#'
#' This function is a modified version of the DimPlot function in Seurat
#'
#' @param seu_obj A complete Seurat object
#' @param colors Colors to label the clusters
#' @param pt.size Size of the data points
#' @param label.box Whether or not to label the cell type name
#' @param label.size Font size of the labels
#' @param title Main title of the plot
#' @return A ggplot object
#' @export
modified_dimplot <- function(
  seu_obj,
  colors = NULL,
  pt.size = 0.5,
  label.box = T,
  label.size = 6,
  title = "UMAP plot"
  ){
  if(is.null(colors)){
    colors = scales::hue_pal()(length(levels(seu_obj)))
  }
  x0<-min(seu_obj@reductions$umap@cell.embeddings[,1])
  x1<-x0+(max(seu_obj@reductions$umap@cell.embeddings[,1])-min(seu_obj@reductions$umap@cell.embeddings[,1]))/8
  y0<-min(seu_obj@reductions$umap@cell.embeddings[,2])
  y1<-y0+(max(seu_obj@reductions$umap@cell.embeddings[,2])-min(seu_obj@reductions$umap@cell.embeddings[,2]))/8
  DimPlot(seu_obj, label = T, cols = colors,label.box = label.box, 
        label.size = label.size, pt.size = pt.size, raster = F)+NoLegend()+NoAxes()+
  geom_segment(aes(x=x0, xend = x1 , y=y0, yend = y0), size=0.8,
               arrow = arrow(length = unit(0.2,"cm"))) +
  geom_segment(aes(x=x0, xend = x0 , y=y0, yend = y1), size=0.8,
               arrow = arrow(length = unit(0.2,"cm"))) +
  xlab("UMAP_1")+theme(axis.title.x = element_text(hjust = 0.05, size = 12))+
  ylab('UMAP_2')+theme(axis.title.y = element_text(hjust = 0.05, angle = 90, size = 12))+
  ggtitle(title)+theme(plot.title = element_text(hjust = 0.5))
}