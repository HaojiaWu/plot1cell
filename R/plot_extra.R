#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Functions
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#' Plot the qPCR results
#'
#' This function is to process the Cq file generated by the qPCR machine 
#' in our lab. 
#' @param qPCR_file Path to the Cq file from qPCR machine. Must be in csv format.
#' @param metadata_file Path to the metadata file to assign your samples into groups. Must be in csv format.
#' @param ref_gene The gene used for normalization. e.g. GAPDH. 
#' @param ref_sample The sample used as reference sample. e.g. sample from the control group.
#' @param file_name The output file name
#' @return There are three output files produced from this script. The csv file contains the average quantitative value (2^-ΔΔCt) for each sample (and each gene) after normalized by the reference gene (e.g. GAPDH) and the reference sample (e.g. sample from the control group). This file can be input into Graphpad Prism. The txt file includes all statistics from comparisons of any two given groups. If the run has two groups only, Welch's t test will be performed. Otherwise, one-way ANOVA with post-hoc Tukey's test will be performed. Finally, the tiff file is a boxplot graph to visualize the gene expression across groups. 
#' @export
plot_qpcr<-function(
  qPCR_file, 
  metadata_file, 
  ref_gene, 
  ref_sample, 
  file_name
  ){
  qPCR_data<-read.csv(qPCR_file)
  meta.data<-read.csv(metadata_file)
  meta.data[,'Sample']<-as.character(meta.data[,'Sample'])
  qPCR_data<-qPCR_data[,c('Target','Sample','Cq')]
  qPCR_data[,'Sample']<-as.character(qPCR_data[,'Sample'])
  all.genes<-unique(qPCR_data$Target)
  gene.exist<-ref_gene %in% all.genes
  if (!gene.exist)
    stop("Reference gene ID is not correct! Please double check!")
  all.samples<- unique(qPCR_data$Sample)
  sample.exist<-ref_sample %in% all.samples
  if (!sample.exist)
    stop("The calibrator sample ID is not correct! Please double check!")
  qPCR_data<-qPCR_data %>% na_if("") %>% na.omit
  qPCR_data<-qPCR_data %>% group_by(Sample, Target) %>% summarize(Mean = mean(Cq, na.rm=TRUE))
  genes<-setdiff(all.genes,ref_gene)
  data_rebuild<-dcast(qPCR_data, Sample ~ Target)
  data_rebuild<-data_rebuild[data_rebuild$Sample %in% meta.data$Sample,]
  new.data<-matrix(data = NA, ncol = length(genes), nrow = nrow(data_rebuild))
  for (i in 1: length(genes)){
    new.data[,i]<-data_rebuild[,genes[i]]-data_rebuild[,ref_gene]
  }
  new.data<-data.frame(new.data)
  colnames(new.data)<-genes
  new.data[,'Sample']<-data_rebuild[,'Sample']
  ref.ct<-new.data[new.data[,'Sample']==ref_sample,][-ncol(new.data)]
  ref.ct<-as.numeric(ref.ct)
  new.data2<-new.data[,-ncol(new.data)]
  new.data2<-data.frame(new.data2)
  for (i in 1: length(genes)){
    new.data2[,i]<-new.data[,genes[i]]-ref.ct[i]
  }
  new.data2<-2^(-new.data2)
  colnames(new.data2)<-genes
  new.data2[,'Sample']<-data_rebuild[,'Sample']
  new.data2<-merge(new.data2, meta.data, by = 'Sample')
  write.csv(new.data2, paste0(file_name, '_processed.csv'))
  groups<-unique(as.character(new.data2$Group))
  if (length(groups) != 1) {
    if (length(groups)<=2) {
      test_stat<-lapply(new.data2[c(-1, -ncol(new.data2))], function(x) t.test(x ~ new.data2$Group))
      sink(paste0(file_name, '_stats.txt'))
      print(test_stat)
      sink() 
    } else {
      new.data2$Group <- factor(new.data2$Group, levels= groups)
      formulae <- lapply(colnames(new.data2)[2:(ncol(new.data2)-1)], function(x) as.formula(paste0(x, " ~ Group")))
      res1 <- lapply(formulae, function(x) summary(aov(x, data = new.data2)))
      names(res1) <- format(formulae)
      res2 <- lapply(formulae, function(x) TukeyHSD(aov(x, data = new.data2)))
      names(res2) <- format(formulae)
      res<-c(res1, res2)
      sink(paste0(file_name, '_stats.txt'))
      print(res)
      sink() 
    }
  }
  new.data2<-melt(new.data2)
  group<-unique(as.character(meta.data$Group))
  new.data2[,'Group']<-factor(new.data2[,'Group'], levels = group)
  p1<-ggplot(new.data2, aes(x=Group, y=value))+
    stat_boxplot(geom ='errorbar', width = 0.2) +
    geom_boxplot( aes(Group, value, color=Group), width=0.5) +    
    geom_quasirandom()+
    ylab("Fold change")+xlab("")+ 
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          strip.text = element_text(size = 12),axis.title = element_text(size = 12),
          axis.text.x = element_blank(),legend.position = "bottom", legend.text = element_text(size = 12), 
          legend.title = element_text(size = 12, face = 'bold'))+
    facet_wrap(~variable, scales = 'free_y', ncol = 2)
  num=length(genes)
  if (num ==1){
    tiff(paste0(file_name,'.tiff'), units="in", width=3, height=3, res=300)
    print(p1)
    invisible(dev.off())
  } else if (num ==2){
    tiff(paste0(file_name,'.tiff'), units="in", width=6, height=3, res=300)
    print(p1)
    invisible(dev.off())
  } else if (num>2 & (num %% 2) == 0) {
    tiff(paste0(file_name,'.tiff'), units="in", width=6, height=3*(num%/%2), res=300)
    print(p1)
    invisible(dev.off())
  } else {
    tiff(paste0(file_name,'.tiff'), units="in", width=6, height=3*(1+num%/%2), res=300)
    print(p1)
    invisible(dev.off())
  }
}
