#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Violin plot for a single gene across groups
#'
#' This function generates violin plot(s) to compare the expression of a single gene across 
#' different groups or cell types. It is designed for visualizing a complicated scenario: 
#' Gene expression on multiple cell types and multiple conditions.
#'
#' @param seu_obj A complete Seurat object.
#' @param feature Gene name. Only one gene is allowed.
#' @param celltypes Cell types of interest. By default, all cell types are included.
#' @param groups Groups selected for plotting. Support multiple groups.
#' @param add.dot Whether or not to add points on the violins.
#' @param font.size Font size for the labels.
#' @param pt.size Point size for the data points on the violin.
#' @param splitby Group to split the gene expression. Only works when length(groups)==1.
#' @param alpha Point transparency. value from 0 to 1.
#' @param strip.color Colors for the strip background.
#' @return A ggplot object
#' @export


complex_vlnplot_single <- function(
  seu_obj,
  feature,
  celltypes=NULL,
  groups,
  add.dot = T,
  font.size=14,
  pt.size=0.1,
  splitby=NULL,
  alpha=0.5,
  strip.color=NULL
){
  if(length(feature)>1){
    stop("Only one gene is allowed in this method. Please use complex_vlnplot_multiple if you want to plot multiple genes.")
  }
  if(is.null(celltypes)){
    celltypes = levels(seu_obj)
  } 
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
  set.seed(seed = 42)
  noise <- rnorm(n = length(x = gene_count[,feature])) / 100000 ## This follows the same data processing for VlnPlot in Seurat
  gene_count[, feature]<-gene_count[,feature]+noise
  gene_count$celltype<-factor(gene_count$celltype, levels = celltypes)
  if (length(groups)==1) {
    p<-ggplot(gene_count, aes_string(x = groups, y = feature, fill = groups)) +
      geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=alpha, color="pink")+
      xlab("") +
      ylab("") +
      ggtitle(feature)+
      theme(panel.background = element_rect(fill = "white",colour = "black"),
            axis.title = element_text(size = font.size), 
            axis.text.x = element_text(size = font.size, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size=(font.size-2)),
            strip.text = element_text( size = font.size),
            legend.title = element_blank(),
            legend.position = 'none',
            plot.title = element_text(size=(font.size),face = "bold", hjust = 0.5))
    if(add.dot){
      p = p + geom_quasirandom(size=pt.size, alpha=0.2)
    }
    if(is.null(splitby)){
      p <- p + facet_wrap(~celltype, ncol = 1, strip.position = "right")
      g <- change_strip_background(p, type = 'right',  strip.color = strip.color)
      print(grid::grid.draw(g))
    } else {
      p<-p + facet_grid(as.formula(paste("celltype","~", splitby)), scales = "free_x") 
      g <- change_strip_background(p, type = 'both',  strip.color = strip.color)
      print(grid::grid.draw(g))
    }
    
  } else {
    if(is.null(splitby)){
    all_levels<-list()
    for(i in 1:length(groups)){
      if (is.null(levels(seu_obj@meta.data[,groups[i]]))){
        seu_obj@meta.data[,groups[i]] <-factor(seu_obj@meta.data[,groups[i]], levels = names(table(seu_obj@meta.data[,groups[i]])))
      }
      group_level<-levels(seu_obj@meta.data[,groups[i]])
      all_levels[[i]]<-group_level
    }
    all_levels<-as.character(unlist(all_levels))
    gene_count<-reshape2::melt(gene_count[,c(feature, groups, "celltype")], measure.vars  = groups)
    gene_count$value<-factor(gene_count$value, levels = all_levels)
    gene_count$celltype<-factor(gene_count$celltype, levels = celltypes)
    p<-ggplot(gene_count, aes_string(x="value", y=feature, fill="value"))+
      geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=alpha, color="pink")+
      xlab("") + ylab("") + ggtitle(feature) +
      theme(panel.background = element_rect(fill = "white",colour = "black"),
            axis.title = element_text(size = font.size), 
            axis.text.x = element_text(size = font.size, angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size=(font.size-2)),
            strip.text = element_text( size = font.size),
            legend.title = element_blank(),
            legend.position = 'none',
            plot.title = element_text(size=(font.size),face = "bold", hjust = 0.5))+
      facet_grid(celltype~variable, scales = 'free_x')
    if(add.dot){
      p = p + geom_quasirandom(size=pt.size, alpha=0.2)
    }
    g <- change_strip_background(p, type = 'both',  strip.color = strip.color)
    print(grid::grid.draw(g))
    } else {
      count_list<-list()
      for(i in 1:length(groups)){
        df1<-gene_count[, c(groups[i],splitby[i],feature, "celltype")]
        colnames(df1)[1:2]<-c("group","split")
        df1$new_group<-groups[i]
        count_list[[i]]<-df1
      }
      new_count<-do.call("rbind", count_list)
      new_count$celltype<-factor(new_count$celltype, levels = celltypes)
      group_level<-list()
      for(i in 1:length(groups)){
        if (is.null(levels(seu_obj@meta.data[,groups[i]]))){
          seu_obj@meta.data[,groups[i]] <-factor(seu_obj@meta.data[,groups[i]], levels = names(table(seu_obj@meta.data[,groups[i]])))
        }
        group_level[[i]]<-levels(seu_obj@meta.data[,groups[i]])
      }
      group_level<-as.character(unlist(group_level))
      new_count$group<-factor(new_count$group, levels=group_level)
      fill_x1<-grDevices::rainbow(length(groups), alpha = 0.5)
      fill_x2<-list()
      for(i in 1:length(splitby)){
        n_col<-unique(gene_count[, splitby[i]])
        fill_x2[[i]]<-scales::hue_pal(l=90)(length(n_col))
      }
      fill_x2<-as.character(unlist(fill_x2))
      fill_x <- c(fill_x1, fill_x2)
      fill_y <- scales::hue_pal(l=90)(length(celltypes))
      p<- ggplot(new_count, aes_string(x="group", y=feature, fill="group"))+
        geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=alpha, color="pink")+
        xlab("") + ylab("") + ggtitle(feature) +
        theme(panel.background = element_rect(fill = "white",colour = "black"),
              axis.title = element_text(size = font.size), 
              axis.text.x = element_text(size = font.size, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size=(font.size-2)),
              strip.text = element_text( size = font.size),
              legend.title = element_blank(),
              legend.position = 'none',
              plot.title = element_text(size=(font.size),face = "bold", hjust = 0.5)) +
        facet_nested(celltype ~ new_group + split, scales = "free_x", 
                     strip = strip_nested( background_x = elem_list_rect(fill = fill_x),
                                           background_y = elem_list_rect(fill = fill_y)))
      if(add.dot){
        p = p + geom_quasirandom(size=pt.size, alpha=0.2)
        print(p)
      } else {
        p
      }
    }
  }
}

#' Violin plot for multiple genes across groups
#'
#' This function generates violin plot(s) to compare the expression of multiple genes across 
#' different groups or cell types. It is designed for visualizing a complicated scenario: 
#' Gene expression of multiple genes on multiple cell types across groups.
#'
#' @param seu_obj A complete Seurat object
#' @param features Gene name. Only one gene is allowed.
#' @param celltypes Cell types of interest. By default, all cell types are included.
#' @param group Only one groupID is allowed.
#' @param add.dot Whether or not to add points on the violins.
#' @param font.size Font size for the labels.
#' @param pt.size Point size for the data points on the violin
#' @param alpha Point transparency. value from 0 to 1.
#' @param strip.color Colors for the strip background
#' @return A ggplot object
#' @export

complex_vlnplot_multiple <- function(
  seu_obj,
  features,
  celltypes=NULL,
  group,
  add.dot = T,
  font.size=12,
  pt.size=0.1,
  alpha=0.01,
  strip.color = NULL
){
  if(length(features)<2){
    stop("At least two genes are required. For single gene violin plot, please use complex_vlnplot_single instead.")
  }
  if(length(group)>1){
    stop("Use violin plot to show multiple genes in multiple group categories across multiple cell types will look too messy. Please use one group ID only.")
  }
  if(is.null(celltypes)){
    celltypes = levels(seu_obj)
  } 
  gene_count<-extract_gene_count(seu_obj=seu_obj, features = features, cell.types = celltypes, meta.groups = group)
  if (is.null(levels(seu_obj@meta.data[,group]))){
    seu_obj@meta.data[,group] <-factor(seu_obj@meta.data[,group], levels = names(table(seu_obj@meta.data[,group])))
  }
  group_level<-levels(seu_obj@meta.data[,group])
  gene_count[,group]<-factor(gene_count[,group],levels = group_level)
  for(i in 1:length(features)){
    set.seed(seed = 42)
    noise <- rnorm(n = length(x = gene_count[,features[i]])) / 100000 ## This follows the same data processing for VlnPlot in Seurat
    gene_count[, features[i]]<-gene_count[,features[i]]+noise
  }
  gene_count$Cell<-rownames(gene_count)
  gene_count <- reshape2::melt(gene_count, id.vars = c("Cell","celltype",group), measure.vars = features,
                               variable.name = "Genes", value.name = "Expr")
  gene_count[, group]<-factor(gene_count[, group], levels = group_level)
  gene_count[, "celltype"]<-factor(gene_count[, "celltype"], levels = celltypes)
  p<-ggplot(gene_count, aes_string(group, "Expr", fill = group)) +
    geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=0.5, color="pink") +
    scale_y_continuous(expand = c(0, 0), position="left", labels = function(x)
      c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
    theme(panel.background = element_rect(fill = "white",colour = "black"),
          axis.title = element_text(size = font.size), 
          axis.text.x = element_text(size = font.size, angle = 45, hjust = 1, vjust = 1),
          axis.text.y = element_text(size=(font.size)),
          strip.text = element_text( size = font.size),
          legend.title = element_blank(),
          legend.position = 'none') +
    facet_grid(celltype~Genes, scales =  'free_x') +
    xlab("") + ylab("")
  if(add.dot){
    p = p + geom_quasirandom(size=pt.size, alpha=alpha)
  }
  g <- change_strip_background(p, type = 'both',  strip.color = strip.color)
  print(grid::grid.draw(g))
}

