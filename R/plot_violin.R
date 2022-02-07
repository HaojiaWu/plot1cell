#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Violin plot for a single gene across groups
#'
#' This function generates violin plot(s) to compare the expression of a single gene across 
#' different groups or cell types. It is designed for visualizing a complicated scenario: 
#' Gene expression on multiple cell types and multiple conditions.
#'
#' @param seu A complete Seurat object
#' @param feature Gene name. Only one gene is allowed.
#' @param cell.types Cell types of interest. By default, all cell types are included.
#' @param groups Groups to split the plot. Support multiple groups.
#' @param add.dot Whether or not to add points on the violins.
#' @param font.size Font size for the labels.
#' @param pt.size Point size for the data points on the violin
#' @return A ggplot object
#' @export

complex_vlnplot_single <- function(
  seu,
  feature,
  cell.types=NULL,
  groups,
  add.dot = T,
  font.size=14,
  pt.size=0.1,
  split.by=NULL
){
  if(is.null(cell.types)){
    cell.types = levels(seu)
  } 
  gene_count<-extract_gene_count(seu=seu, features = feature, cell.types = cell.types, meta.groups = c(groups, split.by))
  set.seed(seed = 42)
  noise <- rnorm(n = length(x = gene_count[,feature])) / 100000
  gene_count[, feature]<-gene_count[,feature]+noise
  if (length(groups)==1) {
    if(length(cell.types)==1){
      p<-ggplot(gene_count, aes_string(x=groups, y=feature, fill=groups))+
        geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=0.5, color="pink")+
        xlab("") + ylab("") + ggtitle(feature) +
        theme(panel.background = element_rect(fill = "white",colour = "black"),
              axis.title = element_text(size = font.size), 
              axis.text.x = element_text(size = font.size, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size=(font.size-2)),
              legend.title = element_blank(),
              legend.position = 'none',
              plot.title = element_text(size=(font.size+2), hjust = 0.5))
      if(add.dot){
        p = p + geom_quasirandom(size=pt.size, alpha=0.2)
      }
      if(!is.null(split.by)){
        p = p + facet_wrap(as.formula(paste("~", split.by)), scales = 'free_x')
      }
      p
    } else {
      if(is.null(split.by)){
        plot_list<-list()
        for(i in 1:length(cell.types)){
          cell_count <- gene_count[gene_count$celltype==cell.types[i],]
          p<-ggplot(cell_count, aes_string(x=groups, y=feature, fill=groups))+
            geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=0.5, color="pink")+
            xlab("") + ylab(cell.types[i]) +
            theme(panel.background = element_rect(fill = "white",colour = "black"),
                  legend.position = "none", 
                  axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank(), 
                  axis.title.y = element_text(size = font.size, angle = 0), 
                  axis.text.y = element_text(size = (font.size-2)),
                  plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm") ) 
          if(add.dot){
            p = p + geom_quasirandom(size=pt.size, alpha=0.2)
          }
          plot_list[[i]]<-p
        }
        plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
          theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = font.size), axis.ticks.x = element_line())
        p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)  + patchwork::plot_annotation(title = feature) & theme(plot.title = element_text(hjust = 0.5, size = (font.size +2)))
        p
      } else {
        p<-ggplot(gene_count, aes_string(x = "group", y = feature, fill = "group")) +
          facet_grid(as.formula(paste("celltype","~","group2")), scales = "free_x") +
          geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=0.5, color="pink")+
          xlab("") +
          ylab(paste(feature,"expression")) +
          theme_bw() +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                strip.text = element_text( size = font.size),
                axis.text.x = element_text(size=(font.size-2), angle = 45, hjust = 1, vjust = 1),
                axis.title.y = element_text(size = font.size),
                legend.position = "none")
        if(add.dot){
          p = p + geom_quasirandom(size=pt.size, alpha=0.2)
        }
        g <- ggplot_gtable(ggplot_build(p))
        strip_t <- which(grepl('strip-t', g$layout$name))
        strip_r <- which(grepl('strip-r', g$layout$name))
        strip_both<-c(strip_t, strip_r)
        ncol <- length(cell.types) + length(names(table(gene_count[,split.by])))
        fills <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(ncol)
        k <- 1
        for (i in strip_both) {
          j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
          g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
          k <- k+1
        }
        print(grid::grid.draw(g))
      }
    }
  } else {
    if(!is.null(split.by)){
      stop("This function does not support spliting multiple groups. Plots will look too messy! Please select one group only in the 'groups' parameter if you want to use 'split.by'.")
    }
    if(length(cell.types)==1){
      gene_count<-melt(gene_count[,c(feature, groups)], measure.vars  = groups)
      p<-ggplot(gene_count, aes_string(x="value", y=feature, fill="value"))+
        geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=0.5, color="pink")+
        xlab("") + ylab("") + ggtitle(feature) +
        theme(panel.background = element_rect(fill = "white",colour = "black"),
              axis.title = element_text(size = font.size), 
              axis.text.x = element_text(size = font.size, angle = 45, hjust = 1, vjust = 1),
              axis.text.y = element_text(size=(font.size-2)),
              strip.background =element_rect(fill="lemonchiffon1"),
              strip.text = element_text( size = font.size),
              legend.title = element_blank(),
              legend.position = 'none',
              plot.title = element_text(size=(font.size+2), hjust = 0.5))+
        facet_wrap(~variable, scales = 'free_x')
      if(add.dot){
        p = p + geom_quasirandom(size=pt.size, alpha=0.2)
      }
      p
    } else {
      plot_list1<-list()
      plot_list2<-list()
      for(i in 1:length(groups)){
        group=groups[i]
        cell_count<-gene_count[,c(feature, group, "celltype")]
        for(j in 1:length(cell.types)){
          cell_count2 <- cell_count[cell_count$celltype==cell.types[j],]
          p<-ggplot(cell_count2, aes_string(x=group, y=feature, fill=group))+
            geom_violin(scale = 'width', adjust = 1, trim = TRUE, size=0.3, alpha=0.5, color="pink")+
            xlab("") + ylab(cell.types[j]) +
            theme(panel.background = element_rect(fill = "white",colour = "black"),
                  legend.position = "none", 
                  axis.text.x = element_blank(), 
                  axis.ticks.x = element_blank(), 
                  axis.title.y = element_text(size = font.size, angle = 0), 
                  axis.text.y = element_text(size = (font.size-2)),
                  plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm") ) 
          if(add.dot){
            p = p + geom_quasirandom(size=pt.size, alpha=0.2)
          }
          plot_list1[[j]]<-p
        }
        plot_list1[[length(plot_list1)]]<- plot_list1[[length(plot_list1)]] +
          theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 1, size = font.size), axis.ticks.x = element_line(), axis.title.x = element_text(size=font.size))+ xlab(group)
        p2<-patchwork::wrap_plots(plotlist = plot_list1, ncol = 1) 
        plot_list2[[i]]<-p2
        plot_list1<-list()
      }
      p<-patchwork::wrap_plots(plotlist = plot_list2) + patchwork::plot_annotation(title = feature) & theme(plot.title = element_text(hjust = 0.5, size = (font.size+2)))
      p
    }
  }
}