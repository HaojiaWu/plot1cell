#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Modified CoveragePlot
#'
#' This function is a modified version of the CoveragePlot in Signac that allows gene name input and group splitting
#'
#' @param sig_obj A complete Signac object
#' @param gene Gene name
#' @param slide_windows_up Upper range of the slide window
#' @param slide_windows_down Lower range of the slide window
#' @param split_by Group ID in metadata to split peaks
#' @return A ggplot object
#' @export
modified_coverageplot <- function(
  sig_obj, 
  gene, 
  slide_windows_up, 
  slide_windows_down,
  split_by
){
  celltypes<-levels(sig_obj)
  DefaultAssay(sig_obj)<-'peaks'
  sig_obj@meta.data$disease<-sig_obj@meta.data[,split_by]
  sig_obj@meta.data$celltype<-sig_obj@active.ident
  sig_obj@meta.data$id2<-paste(sig_obj@meta.data$celltype,sig_obj@meta.data$disease,  sep = "_")
  sig_obj<- SetIdent(sig_obj, value = 'id2')
  p2 <- Signac::CoveragePlot(
    object = sig_obj,
    region = gene,
    annotation = T,
    peaks = T, 
    extend.upstream = slide_windows_up,
    extend.downstream = slide_windows_down
  )
  gene.ranges <- genes(EnsDb.Hsapiens.v86)
  seqlevelsStyle(gene.ranges) <- 'UCSC'
  gene.ranges <- keepStandardChromosomes(gene.ranges, pruning.mode = 'coarse')
  region2 <- GRangesToString(subset(gene.ranges, symbol==gene))
  if(region2=='--') {
    next
  }
  if(length(region2)>1){
    region2<-region2[1]
  }
  region2<-add_sliding_windows(region = region2, up_stream = slide_windows_up, down_stream = slide_windows_down)
  xlabel<-p2$patches$plots[1][[1]]$labels$x
  ylabel<-p2$patches$plots[1][[1]]$labels$y
  data_plot<-p2$patches$plots[1][[1]]$data
  data_plot$celltype<-gsub("_[^_]*$", "", data_plot$group )
  data_plot$disease<-gsub(".*_", "", data_plot$group )
  data_plot$celltype <- factor(data_plot$celltype, levels = celltypes)
  a1<-Signac:::PeakPlot(sig_obj,region = region2)
  a2<-p2$patches$plots[2][[1]]
  group_id<-names(table(data_plot$disease))
  a3 <- ggplot(
    data = data_plot[data_plot$disease==group_id[1],],
    mapping = aes(x = position, y = coverage, fill = celltype)
  ) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    xlab("")+ylab(as_string(ylabel))+ggtitle(group_id[1])+
    theme_browser(legend = FALSE) +
    facet_grid(celltype ~.)+
    theme(
      panel.spacing.y = unit(x = 0, units = "line"),
      plot.title = element_text(size = 14,hjust = 0.5, face = 'bold'),
      strip.text.y = element_text(angle = 0),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank(),
      axis.line.x = element_blank())
  
  a4 <- ggplot(
    data = data_plot[data_plot$disease==group_id[2],],
    mapping = aes(x = position, y = coverage, fill = celltype)
  ) +
    geom_area(stat = "identity") +
    geom_hline(yintercept = 0, size = 0.1) +
    xlab("")+
    ylab("")+
    ggtitle(group_id[2])+
    theme_browser(legend = FALSE) +
    facet_grid(celltype ~.)+
    theme(
      plot.title = element_text(size = 14,hjust = 0.5, face = 'bold'),
      strip.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.line = element_blank())
  a5 <- a2 + ylab("")
  a6 <- a1 + ylab("")
  p3 <- a3 + a2 + a1  + plot_layout(heights = c(6,1,1))
  p4 <- a4 + a5 + a6  + plot_layout(heights = c(6,1,1))
  plot_grid(p3, p4, ncol = 2, rel_widths = c(2.4, 2))
}



#' Add extended windows for coverage plot
#'
#' This function is to add an extended window for ploting
#'
#' @param region Gene region for plotting
#' @param up_stream Upper range of the slide window
#' @param down_stream Lower range of the slide window
#' @return Gene range for plot
#' @export
add_sliding_windows <- function(
  region,
  up_stream, 
  down_stream
){
  reg_list<-strsplit(region, split = "-")
  reg_list<-as.character(unlist(reg_list))
  new_reg<-paste(reg_list[1], as.integer(reg_list[2])-down_stream, as.integer(reg_list[3])+up_stream, sep = "-")
  new_reg
}