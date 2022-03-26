#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Plot cell fractions across groups
#'
#' This function is to show the cell fraction changes across groups.
#'
#' @param seu_obj A complete Seurat object
#' @param celltypes Cell types to be included in the plot. Default: all cell types.
#' @param groupby The group to show on x axis. One of the column names in meta.data.
#' @param show_replicate Whether or not to show the individual replicate on the graph. If TRUE, the replicate column name needs to specify in the argument rep_colname.
#' @param rep_colname The column name for biological replicates in the meta data.
#' @param strip.color Colors for the strip background
#' @return A ggplot object
#' @export

plot_cell_fraction<-function(
  seu_obj, 
  celltypes=NULL, 
  groupby, 
  show_replicate = FALSE, 
  rep_colname = NULL,
  strip.color = NULL
  ){
  meta_data<-seu_obj@meta.data
  meta_data$celltype<-as.character(Idents(seu_obj))
  groupby_level<-levels(seu_obj@meta.data[,groupby])
  if (is.null(groupby_level)){
    seu_obj@meta.data[,groupby] <-factor(seu_obj@meta.data[,groupby], levels = names(table(seu_obj@meta.data[,groupby])))
    groupby_level<-levels(seu_obj@meta.data[,groupby])
  } 
  if(!show_replicate){
    freq_df<-prop.table(table(meta_data[,"celltype"], meta_data[, groupby]), margin = 2)
    freq_df<-data.frame(freq_df)
    colnames(freq_df)[1:2]<-c("celltype","group")
    freq_df$Freq<-freq_df$Freq*100
  } else {
    if(is.null(rep_colname)){
      stop("Please specify the replicate colname in your meta data!")
    } else {
      meta_data$new_group<-paste(meta_data[,"celltype"], meta_data[,groupby], meta_data[, rep_colname], sep = "___")
      freq_df<-data.frame(table(meta_data$new_group))
      freq_df$Var1<-as.character(freq_df$Var1)
      freq_df$celltype <- as.character(lapply(X=strsplit(freq_df$Var1, split = "___"),FUN = function(x){x[[1]]}))
      freq_df$group <- as.character(lapply(X=strsplit(freq_df$Var1, split = "___"),FUN = function(x){x[[2]]}))
      freq_df$replicate <- as.character(lapply(X=strsplit(freq_df$Var1, split = "___"),FUN = function(x){x[[3]]}))
      total_cell<-aggregate(Freq ~ replicate, freq_df, sum)
      colnames(total_cell)[2]<-"Total"
      freq_df <- merge(freq_df, total_cell, by="replicate")
      freq_df$Freq<-freq_df$Freq*100/freq_df$Total
    }
  }
  freq_df$group <- factor(freq_df$group, levels = groupby_level)
  if(is.null(celltypes)){
    celltypes <- levels(seu_obj)
  }
  freq_df <- freq_df[freq_df$celltype %in% celltypes, ]
  freq_df$celltype <- factor(freq_df$celltype, levels = celltypes)
 p <- ggplot(freq_df, aes(group, Freq, fill=group))+
    geom_bar(position = "dodge",  stat = "summary", fun='mean', width = 0.6)+
    stat_summary(fun.data = mean_se, geom = "errorbar", width=0.2, size=1, color='midnightblue', alpha=0.8)+
    ylab('Percentage of cells')+xlab('')+
    scale_fill_manual(values = brewer.pal(8, 'Set2'))+
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          strip.text = element_text(size = 12),axis.title = element_text(size = 14),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 12),
          legend.position = "none", axis.text.y = element_text(size = 12),
          legend.text = element_text(size = 12), legend.title = element_text(size = 12, face = 'bold'), 
          plot.title = element_text(hjust =0.5))
 
 if(show_replicate){
   p <- p +   geom_quasirandom(size=1,width = 0.2, color='midnightblue', alpha=0.8, groupOnX = F)+facet_wrap(~celltype, scales = 'free_y', ncol = length(celltypes))
   g <- change_strip_background(p, type = 'top', n.color = length(celltypes), strip.color = strip.color)
   print(grid.draw(g))
 } else {
   p
 }
}

