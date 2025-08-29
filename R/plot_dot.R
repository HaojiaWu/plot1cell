#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Plot single gene across groups
#'
#' This function can be used for plotting a single gene expression across 
#' different groups in a study with complex group design.
#'
#' @param seu_obj A complete Seurat object.
#' @param feature Gene name. Only one gene is allowed.
#' @param celltypes Cell types to be included in the dot plot. Default: all cell types.
#' @param groups The group to show on x axis. One of the column names in meta.data.
#' @param splitby The group to separate the gene expression. One of the column names in meta.data.
#' @param scale.by Methods to scale the dot size. "radius" or "size".
#' @param color.palette Color for gene expression.
#' @param strip.color Colors for the strip background.
#' @param font.size Font size for the labels.
#' @param do.scale Whether or not to scale the dot when percentage expression of the gene is less than 20.
#' @param vmin Clip lower limit for gene expression value (after z-scaling).
#' @param vmax Clip upper limit for gene expression value (after z-scaling).
#' @return A ggplot object
#' @export

complex_dotplot_single <- function(
  seu_obj, 
  feature, 
  celltypes=NULL,
  groups,
  splitby=NULL,
  color.palette = NULL,
  font.size = 12,
  strip.color=NULL,
  do.scale=T,
  scale.by='radius',
  vmin = NULL, 
  vmax = NULL
){
  .clip_vals <- function(x, vmin = NULL, vmax = NULL){
      if(!is.null(vmin)) x[x < vmin] <- vmin
      if(!is.null(vmax)) x[x > vmax] <- vmax
      x
    }
  if(is.null(color.palette)){
    color.palette <- colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
  }
  scale.func <- switch(
    EXPR = scale.by,
    'size' = scale_size,
    'radius' = scale_radius,
    stop("'scale.by' must be either 'size' or 'radius'")
  ) ### This function is from Seurat https://github.com/satijalab/seurat
  if(is.null(celltypes)){
    celltypes<-levels(seu_obj)
  }
  if(length(groups)==1){
    groups_level<-levels(seu_obj@meta.data[,groups])
    if (is.null(groups_level)){
      seu_obj@meta.data[,groups] <-factor(seu_obj@meta.data[,groups], levels = names(table(seu_obj@meta.data[,groups])))
      groups_level<-levels(seu_obj@meta.data[,groups])
    } 
    
    if(!is.null(splitby)){
      if (is.null(levels(seu_obj@meta.data[,splitby]))){
        seu_obj@meta.data[,splitby] <-factor(seu_obj@meta.data[,splitby], levels = names(table(seu_obj@meta.data[,splitby])))
      }
      splitby_level<-levels(seu_obj@meta.data[,splitby])
      count_df<-extract_gene_count(seu_obj, features = feature, cell.types = celltypes, meta.groups = c(groups,splitby))
      count_df$new_group<-paste(count_df[,groups], count_df[,"celltype"], count_df[,splitby],sep = "___")
      exp_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){mean(expm1(x))}) 
      pct_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){length(x[x > 0]) / length(x)}) #This is the same data processing as Seurat
      colnames(exp_df)[2]<-"avg.exp"
      colnames(pct_df)[2]<-"pct.exp"
      data_plot<-merge(exp_df, pct_df, by='new_group')
      data_plot$groups <- as.character(lapply(X=strsplit(data_plot$new_group, split = "___"),FUN = function(x){x[[1]]}))
      data_plot$celltype <- as.character(lapply(X=strsplit(data_plot$new_group, split = "___"),FUN = function(x){x[[2]]}))
      data_plot$splitby <- as.character(lapply(X=strsplit(data_plot$new_group, split = "___"),FUN = function(x){x[[3]]}))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$splitby <- factor(data_plot$splitby, levels = splitby_level)
      data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    } else {
      count_df<-extract_gene_count(seu_obj, features = feature, cell.types = celltypes, meta.groups = groups)
      count_df$new_group<-paste(count_df[,groups], count_df[,"celltype"],sep = "___")
      exp_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){mean(expm1(x))})
      pct_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){length(x[x > 0]) / length(x)})
      colnames(exp_df)[2]<-"avg.exp"
      colnames(pct_df)[2]<-"pct.exp"
      data_plot<-merge(exp_df, pct_df, by='new_group')
      data_plot$groups <- as.character(lapply(X=strsplit(data_plot$new_group, split = "___"),FUN = function(x){x[[1]]}))
      data_plot$celltype <- as.character(lapply(X=strsplit(data_plot$new_group, split = "___"),FUN = function(x){x[[2]]}))
      data_plot$groups <- factor(data_plot$groups, levels = groups_level)
      data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    }
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    data_plot$avg.exp <- .clip_vals(data_plot$avg.exp, vmin, vmax)
    p<-ggplot(data_plot, aes(y = celltype, x = groups)) +  
      geom_tile(fill="white", color="white") +
      geom_point(aes( colour=avg.exp, size =pct.exp))  +  
      scale_color_gradientn(colours  =  color.palette )+ 
      theme(panel.background = element_rect(fill = "white", colour = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1, size = font.size),
            plot.title = element_text(size = (font.size +2), hjust = 0.5, face = 'bold'),
            axis.text = element_text(size = font.size),
            legend.text=element_text(size=(font.size-2)),
            legend.key = element_rect(fill = NA, colour = NA) ,
            legend.title = element_text(size = (font.size)),
            strip.text = element_text( size = font.size),
            legend.position="right")+
      ylab("")+xlab("")+ggtitle(feature)+ guides(
  size = guide_legend(
    override.aes = list(shape = 21, fill = NA, colour = "black") 
  )
)
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
      g <- change_strip_background(p, type = 'top', strip.color = strip.color)
      print(grid.draw(g))
    } else {
      p
    }
  } else {  ### group number greater than 1
    gene_count<-extract_gene_count(seu_obj=seu_obj, features = feature, cell.types = celltypes, meta.groups = c(groups, splitby))
    allgroups<-c(groups,splitby )
    for(i in 1:length(allgroups)){
      if (is.null(levels(seu_obj@meta.data[,allgroups[i]]))){
        seu_obj@meta.data[,allgroups[i]] <-factor(seu_obj@meta.data[,allgroups[i]], levels = names(table(seu_obj@meta.data[,allgroups[i]])))
      }
      group_level<-levels(seu_obj@meta.data[,allgroups[i]])
      gene_count[,allgroups[i]]<-factor(gene_count[,allgroups[i]],
                                        levels = group_level)
    }
    gene_count$celltype<-factor(gene_count$celltype, levels = celltypes)
    all_levels<-list()
    for(i in 1:length(groups)){
      if (is.null(levels(seu_obj@meta.data[,groups[i]]))){
        seu_obj@meta.data[,groups[i]] <-factor(seu_obj@meta.data[,groups[i]], levels = names(table(seu_obj@meta.data[,groups[i]])))
      }
      group_level<-levels(seu_obj@meta.data[,groups[i]])
      all_levels[[i]]<-group_level
    }
    all_levels<-as.character(unlist(all_levels))
    data_plot<-list()
    for(i in 1:length(groups)){
      count_df <- gene_count
      count_df$new_group<-paste(gene_count[,groups[i]], gene_count[,"celltype"],sep = "___")
      exp_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){mean(expm1(x))})
      pct_df<-aggregate(.~new_group, data=count_df[,c('new_group',feature)], FUN=function(x){length(x[x > 0]) / length(x)})
      colnames(exp_df)[2]<-"avg.exp"
      colnames(pct_df)[2]<-"pct.exp"
      df1<-merge(exp_df, pct_df, by='new_group')
      df1$groupID <- groups[i]
      data_plot[[i]] <- df1
    }
    data_plot <- do.call("rbind", data_plot)
    data_plot$groups <- as.character(lapply(X=strsplit(data_plot$new_group, split = "___"),FUN = function(x){x[[1]]}))
    data_plot$celltype <- as.character(lapply(X=strsplit(data_plot$new_group, split = "___"),FUN = function(x){x[[2]]}))
    data_plot$groups <- factor(data_plot$groups, levels = all_levels)
    data_plot$celltype <- factor(data_plot$celltype, levels = rev(celltypes))
    data_plot$groupID <- factor(data_plot$groupID, levels = groups)
    data_plot$pct.exp <- round(100 * data_plot$pct.exp, 2)
    data_plot$avg.exp <- scale(data_plot$avg.exp)
    data_plot$avg.exp <- .clip_vals(data_plot$avg.exp, vmin, vmax)

    if(is.null(splitby)){
      p<-ggplot(data_plot, aes(y = celltype, x = groups)) +  
        geom_tile(fill="white", color="white") +
        geom_point(aes( colour=avg.exp, size =pct.exp))  +  
        scale_color_gradientn(colours  =  color.palette )+ 
        theme(panel.background = element_rect(fill = "white", colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1, size = font.size),
              plot.title = element_text(size = (font.size +2), hjust = 0.5, face = 'bold'),
              axis.text = element_text(size = font.size),
              legend.text=element_text(size=(font.size-2)),
              legend.title = element_text(size = (font.size)),
              legend.key = element_rect(fill = NA, colour = NA) ,
              strip.text = element_text( size = font.size),
              legend.position="right")+
        ylab("")+xlab("")+ggtitle(feature)+facet_wrap(~groupID, scales = 'free_x')+ guides(
  size = guide_legend(
    override.aes = list(shape = 21, fill = NA, colour = "black") 
  )
)
      if(do.scale){
        p = p + scale_size(range = c(0, 10))
      } else {
        if(max(data_plot$pct.exp)>=20){
          p = p + scale_size(range = c(0, 10))
        } else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 20))
        }
      }
      g <- change_strip_background(p, type = 'top',  strip.color = strip.color)
      print(grid::grid.draw(g))
    } else {
      df2<-reshape2::melt(gene_count[,c(groups, splitby)], measure.vars  = groups)
      df2<-df2[!duplicated(df2$value),]
      colnames(df2)[colnames(df2) == "value"]<-"groups"
      data_plot2<-list()
      for(i in 1:length(groups)){
        df3<-data_plot[data_plot$groupID==groups[i],]
        df4<-df2[df2$variable==groups[i],c('groups', splitby[i])]
        colnames(df4)[2]<-"split"
        df5<-merge(df3, df4, by='groups')
        data_plot2[[i]]<-df5
      }
      data_plot2<-do.call("rbind", data_plot2)
      fill_x1<-grDevices::rainbow(length(groups), alpha = 0.5)
      fill_x2<-list()
      for(i in 1:length(splitby)){
        n_col<-unique(gene_count[, splitby[i]])
        fill_x2[[i]]<-scales::hue_pal(l=90)(length(n_col))
      }
      fill_x2<-as.character(unlist(fill_x2))
      fill_x <- c(fill_x1, fill_x2)
      p<-ggplot(data_plot2, aes(y = celltype, x = groups)) +  
        geom_tile(fill="white", color="white") +
        geom_point(aes( colour=avg.exp, size =pct.exp))  +  
        scale_color_gradientn(colours  =  color.palette )+ 
        theme(panel.background = element_rect(fill = "white", colour = "black"),
              axis.text.x = element_text(angle = 45, hjust = 1, size = font.size),
              plot.title = element_text(size = (font.size +2), hjust = 0.5, face = 'bold'),
              axis.text = element_text(size = font.size),
              legend.text=element_text(size=(font.size-2)),
              legend.key = element_rect(fill = NA, colour = NA) ,
              legend.title = element_text(size = (font.size)),
              strip.text = element_text( size = font.size),
              legend.position="right")+
        ylab("")+xlab("")+ggtitle(feature)+
        facet_nested(~ groupID + split, scales = "free_x", 
                     strip = strip_nested( background_x = elem_list_rect(fill = fill_x)))+ guides(
  size = guide_legend(
    override.aes = list(shape = 21, fill = NA, colour = "black") 
  )
)
      if(do.scale){
        p = p + scale_size(range = c(0, 10))
      } else {
        if(max(data_plot$pct.exp)>=20){
          p = p + scale_size(range = c(0, 10))
        } else {
          p = p + scale.func(range = c(0, 10), limits = c(0, 20))
        }
      }
      p
    }
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
#' @param celltypes Cell types to be included in the dot plot. Default: all cell types.
#' @param groups Group ID must be one of the column names in the meta.data slot of the Seurat object.
#' @param color.palette Color for gene expression.
#' @param strip.color Colors for the strip background
#' @param vmin Clip lower limit for gene expression value (after z-scaling).
#' @param vmax Clip upper limit for gene expression value (after z-scaling).
#' @param strip.alpha Set transparency for the strip.color
#' @return A ggplot object
#' @export

complex_dotplot_multiple <- function(
  seu_obj, 
  features, 
  celltypes=NULL,
  groups, 
  color.palette = NULL,
  strip.color = NULL,
  vmin = NULL, 
  vmax = NULL,
  strip.alpha =1
  ){
 pb <- progress_bar$new(
   format = "  Ploting [:bar] :percent eta: :eta",
   clear = FALSE, total = length(features), width = 100)
 plot_list<-list()
 for(i in 1:length(features)){
  pp<-invisible(
    complex_dotplot_single(seu_obj = seu_obj, feature = features[i], groups = groups, celltypes = celltypes, vmin = vmin, vmax = vmax)
  )
  pp<-pp$data
  pp$gene <- features[i]
  plot_list[[i]]<-pp
  pb$tick()
  Sys.sleep(1 / length(features))
  }
  all_data<-do.call('rbind', plot_list)
  all_data$gene<-factor(all_data$gene, levels = rev(features)) 
  all_data$celltype <- factor(all_data$celltype, levels = levels(seu_obj))
  if(is.null(color.palette)){
    color.palette <- colorRampPalette(c('grey80','lemonchiffon1','indianred1','darkred'))(255)
  }
  p <- invisible(
    ggplot(all_data, aes(x = groups, y = gene)) +  
    geom_tile(fill="white", color="white") +
    geom_point(aes( colour=avg.exp, size =pct.exp), alpha=0.9)  +  
    scale_color_gradientn(colours  =  color.palette)+ 
    scale_size(range = c(0, 10))+
    theme(panel.background = element_rect(fill = "white", colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size=16),
          plot.title = element_text(size = 16,hjust = 0.5, face = 'bold'),
          axis.text = element_text(size = 12),
          axis.title=element_text(size=8),
          legend.text=element_text(size=8),
          legend.key = element_rect(fill = NA, colour = NA) ,
          legend.title = element_text(size = 12),
          legend.position="right",
          strip.text = element_text(size = 14,colour = 'black',face = 'bold'))+
    ylab("")+xlab("")+ggtitle('')+
    facet_wrap(~celltype, ncol = length(levels(seu_obj))) + guides(
  size = guide_legend(
    override.aes = list(shape = 21, fill = NA, colour = "black") 
  )
)
  )
  if (!is.null(strip.color)) strip.color <- scales::alpha(strip.color, strip.alpha)
  g <- change_strip_background(p, type = 'top',  strip.color = strip.color)
  print(grid.draw(g))
}

