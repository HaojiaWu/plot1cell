#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Plot single gene across groups
#'
#' This function can be used for plotting a single gene expression across 
#' different groups in a study with complex group design. It uses the data
#' output from the DotPlot function in the Seurat package and uses ggplot2 
#' to re-create a dotplot for better visualization. In the output graph, y axis
#' is the cluster name and x axis is the group ID.
#'
#' @param seu_obj A complete Seurat object
#' @param feature Gene name. Only one gene is allowed.
#' @param groupby The group to show on x axis. One of the column names in meta.data.
#' @param splitby The group to separate the gene expression. One of the column names in meta.data.
#' @param scale.by Methods to scale the dot size. "radius" or "size"
#' @param strip.color Colors for the strip background
#' @param do.scale Whether or not to scale the dot when percentage expression of the gene is less than 20.
#' @return A ggplot object
#' @export
complex_dotplot_single <- function(
  seu_obj, 
  feature, 
  groupby,
  splitby=NULL,
  strip.color=NULL,
  do.scale=T,
  scale.by='radius'
){
  if (is.null(levels(seu_obj@meta.data[,groupby]))){
    seu_obj@meta.data[,groupby] <-factor(seu_obj@meta.data[,groupby], levels = names(table(seu_obj@meta.data[,groupby])))
  }
  groupby_level<-levels(seu_obj@meta.data[,groupby])
  levels(seu_obj)<-rev(levels(seu_obj))
  celltypes<-levels(seu_obj)
  celltypes<-gsub("_", ".", celltypes)
  seu_obj@meta.data$celltype<-as.character(seu_obj@active.ident)
  seu_obj@meta.data$celltype<-gsub("_", ".", seu_obj@meta.data$celltype)
  seu_obj<-SetIdent(seu_obj, value='celltype')
  levels(seu_obj)<-celltypes
  if(!is.null(splitby)){
    if (is.null(levels(seu_obj@meta.data[,splitby]))){
      seu_obj@meta.data[,splitby] <-factor(seu_obj@meta.data[,splitby], levels = names(table(seu_obj@meta.data[,splitby])))
    }
    splitby_level<-levels(seu_obj@meta.data[,splitby])
    count_df<-extract_gene_count(seu_obj, features = feature, meta.groups = c(groupby,splitby))
    count_df$new_group<-paste(count_df[,groupby], count_df[,"celltype"], count_df[,splitby],sep = "_")
    exp_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){mean(expm1(x))})
    pct_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){length(x[x > 0]) / length(x)})
    colnames(exp_df)[2]<-"avg.exp"
    colnames(pct_df)[2]<-"pct.exp"
    data_plot<-merge(exp_df, pct_df, by='new_group')
    data_plot$groupby <- as.character(lapply(X=strsplit(data_plot$new_group, split = "_"),FUN = function(x){x[[1]]}))
    data_plot$celltype <- as.character(lapply(X=strsplit(data_plot$new_group, split = "_"),FUN = function(x){x[[2]]}))
    data_plot$splitby <- as.character(lapply(X=strsplit(data_plot$new_group, split = "_"),FUN = function(x){x[[3]]}))
    data_plot$groupby <- factor(data_plot$groupby, levels = groupby_level)
    data_plot$splitby <- factor(data_plot$splitby, levels = splitby_level)
    data_plot$celltype <- factor(data_plot$celltype, levels = celltypes)
  } else {
  count_df<-extract_gene_count(seu_obj, features = feature, meta.groups = groupby)
  count_df$new_group<-paste(count_df[,groupby], count_df[,"celltype"],sep = "_")
  exp_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){mean(expm1(x))})
  pct_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){length(x[x > 0]) / length(x)})
  colnames(exp_df)[2]<-"avg.exp"
  colnames(pct_df)[2]<-"pct.exp"
  data_plot<-merge(exp_df, pct_df, by='new_group')
  data_plot$groupby <- as.character(lapply(X=strsplit(data_plot$new_group, split = "_"),FUN = function(x){x[[1]]}))
  data_plot$celltype <- as.character(lapply(X=strsplit(data_plot$new_group, split = "_"),FUN = function(x){x[[2]]}))
  data_plot$groupby <- factor(data_plot$groupby, levels = groupby_level)
  data_plot$celltype <- factor(data_plot$celltype, levels = celltypes)
  }
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  )
  data_plot$pct.exp <- round(100 * data_plot$pct.exp, 2)
  data_plot$avg.exp <- scale(data_plot$avg.exp)
  p<-ggplot(data_plot, aes(y = celltype, x = groupby)) +  
    geom_tile(fill="white", color="white") +
    geom_point(aes( colour=avg.exp, size =pct.exp))  +  
    scale_color_gradientn(colours  =  colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255))+ 
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 16,hjust = 0.5, face = 'bold'),
          axis.text = element_text(size = 12),
          axis.title=element_text(size=8),
          legend.text=element_text(size=8),
          legend.title = element_text(size = 12),
          legend.position="right")+
    ylab("")+xlab("")+ggtitle(feature)
  if(do.scale){
    p = p + scale_size(range = c(0, 10))
  } else {
    if(max(data_plot$pct.exp)>=20){
      p = p + scale_size(range = c(0, 10))
    } else {
      p = p + scale.func(range = c(0, 10), limits = c(0, 20))
    }
  }
  if(!is.null(splitby)){
    p <- p +facet_wrap(~splitby, scales = 'free_x')
    g <- change_strip_background(p, type = 'top', n.color = length(levels(seu_obj)), strip.color = strip.color)
    print(grid.draw(g))
  } else {
    p
  }
}


#' Plot multiple genes across groups
#'
#' This function allows for visualization of multiple genes in multiple groups. 
#' It takes the single gene expression data generated by PlotSingleGeneGroup,
#' concatenate all data, and produces a dotplot graph where the group ID are in
#' x axis, wrapped by cell types, genes are on the y axis.
#'
#' @param seu_obj A complete Seurat object
#' @param features A vector of gene names.
#' @param groupby Group ID Must be one of the column names in the meta.data slot of the Seurat object.
#' @param strip.color Colors for the strip background
#' @return A ggplot object
#' @export
complex_dotplot_multiple <- function(
  seu_obj, 
  features, 
  groupby, 
  strip.color = NULL
  ){
 pb <- progress_bar$new(
   format = "  Ploting [:bar] :percent eta: :eta",
   clear = FALSE, total = length(features), width = 100)
 features=rev(features)
 plot_list<-list()
 for(i in 1:length(features)){
  pp<-invisible(
    complex_dotplot_single(seu_obj = seu_obj, feature = features[i], groupby = groupby)
  )
  pp<-pp$data
  pp$gene <- features[i]
  plot_list[[i]]<-pp
  pb$tick()
  Sys.sleep(1 / length(features))
  }
  all_data<-do.call('rbind', plot_list)
  all_data$celltype <- factor(all_data$celltype, levels = levels(seu_obj))
  p <- invisible(
    ggplot(all_data, aes(x = groupby, y = gene)) +  
    geom_tile(fill="white", color="white") +
    geom_point(aes( colour=avg.exp, size =pct.exp), alpha=0.9)  +  
    scale_color_gradientn(colours  =  grDevices::colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255))+ 
    scale_size(range = c(0, 10))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 16,hjust = 0.5, face = 'bold'),
          axis.text = element_text(size = 12),
          axis.title=element_text(size=8),
          legend.text=element_text(size=8),
          legend.title = element_text(size = 12),
          legend.position="right",
          strip.text = element_text(size = 14,colour = 'black',face = 'bold'))+
    ylab("")+xlab("")+ggtitle('')+
    facet_wrap(~celltype, ncol = length(levels(seu_obj)))
  )
  g <- change_strip_background(p, type = 'top', n.color = length(levels(seu_obj)), strip.color = strip.color)
  print(grid.draw(g))
}


