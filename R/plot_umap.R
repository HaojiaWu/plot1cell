#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Plot gene expression on umap
#'
#' This function can be used for plotting a single gene or multiple genes expression across 
#' different groups in a seurat featureplot format. 
#'
#' @param seu_obj A complete Seurat object
#' @param features Gene names to be plotted.
#' @param group The group to show on y axis. One of the column names in meta.data.
#' @param select Select the elements within the group to show.
#' @param cols Change the color legend.
#' @param label.size Change the label size.
#' @return A ggplot object
#' @export

complex_featureplot<-function(
  seu_obj, 
  features, 
  group, 
  select=NULL, 
  cols=NULL, 
  label.size=14, 
  order=F
){
  gene_count<-extract_gene_count(seu_obj,features = features, meta.groups = group)
  if (is.null(levels(seu_obj@meta.data[,group]))){
    seu_obj@meta.data[,group] <-factor(seu_obj@meta.data[,group], levels = names(table(seu_obj@meta.data[,group])))
  }
  group_level<-levels(seu_obj@meta.data[,group])
  gene_count[,group]<-factor(gene_count[,group],
                                    levels = group_level)
  if(!is.null(select)){
    gene_count<-gene_count[gene_count[, group] %in% select,]
  }
  colnames(gene_count)[which(colnames(gene_count)==group)]<-"group"
  all_col<-setdiff(colnames(gene_count), features)
  df_list<-list()
  for(i in 1:length(features)){
    df<-gene_count[, c(features[i], all_col)]
    df$gene<-features[i]
    colnames(df)[1]<-"Exp"
    df$Exp<-rescale(df$Exp, to = c(0,5))
    df_list[[i]]<-df
  }
  df_all<-do.call("rbind", df_list)
  if(is.null(cols)){
    cols=colorRampPalette(c('grey90','lemonchiffon1','indianred1','darkred'))(255)
  }
  df_all$gene<-factor(df_all$gene, levels=features)
  if(order){
    df_all$isExpr<-ifelse(df_all$Exp>0, "Yes", "NO")
    p<-ggplot(df_all, aes(UMAP1, UMAP2))+geom_point(color="gray80",size=0.01)+
      geom_point(data = df_all[df_all$isExpr=="Yes",], aes(UMAP1, UMAP2, color=Exp), size=0.01)
  } else {
    p<-ggplot(df_all, aes(UMAP1, UMAP2, color=Exp))+geom_point(size=0.01)
  }
  p +
    scale_color_gradientn(colours  =  cols,
                          na.value = "white", limits=c(quantile(df_all$Exp, 0,na.rm= T), quantile(df_all$Exp, 1,na.rm= T)),
                          breaks = c(quantile(df_all$Exp, 0,na.rm= T), quantile(df_all$Exp, 1,na.rm= T)), labels = c("min","max"))+ 
    theme(panel.background = element_rect(fill = "white", colour = "white"),
          axis.ticks = element_blank(),
          axis.line = element_blank(),
          axis.text = element_blank(),
          axis.title=element_blank(),
          strip.background = element_blank(),
          strip.text = element_text(size=label.size),
          legend.title = element_text(size = 12),
          legend.position="left")+
    facet_grid(group ~ gene)
}