#' SummarizedExperiment objects
#' @param ... arguments in package SummarizedExperiment. See \code{\link[SummarizedExperiment]{assay}}
#' @export
assay <- function(...){SummarizedExperiment::assay(...)}

#' TaxaBarPlot
#' @description Visualization of Wham upstream analysis
#' @importFrom cowplot theme_cowplot
#' @importFrom cowplot plot_grid
#' @importFrom stats aggregate
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr  ggpubr  methods ComplexHeatmap circlize
#' @param object required. The output of WhamBiobakery, WhamEBI or Wham16s
#' @param filter c("all","custom","DE.filter"). Default is "all". "all" shows all samples and top # taxa(indicated in display.number value (default:30)). "custom" shows taxa provided by user and all samples. "DE.filter" shows taxa which pass DE test cutoff(p.cutoff and effect.size.range) and samples which are in the DE.result contrast.
#' @param taxa.level Required, c("k","p","o","c","f","g","s").Collapse taxa at the provided level.
#' @param display.number applied when filter is set as "all". Default is 30.
#' @param custom required when filter is set as "custom". A vector which contain taxa to be shown
#' @param p.cutoff applied when filter is set as "DE.filter". Default is 0.05.
#' @param effect.size.range applied when filter is set as "DE.filter" and aldex.module is set as t.test.Default range is c(0,0)
#' @param annotation optional. Default is the last variable in design. A character vector indicates what column of metadata will be used for annotation bar
#' @param merge_group Default F. Merging samples under the condition.(mean of relative abundance)
#' @param abundance_cutoff When setting filter argument as "all", the bacterium whose the proportion is less than abundance_cutoff will be removed.
#' @param relative_abundance optional. Default is "to_selection". "to_selection" means counts are divided by the total sum of selected taxonomic component in each library. "to_total" means counts are divided by the total sum of all components in each library
#'
#' @return a barplot of taxonomy
#' @export
#'
#' @examples
#' # see the vignette for more details

TaxaBarPlot = function(object,
                       taxa.level = c("k","p","o","c","f","g","s","otu"),
                       filter = c("all","custom","DE.filter"),
                       display.number = 30,
                       custom,
                       p.cutoff = 0.05,
                       effect.size.range,
                       annotation,
                       merge_group = F,
                       abundance_cutoff = 0.001,
                       relative_abundance = c("to_selection","to_total")

){
  ###parameter confirmation
  relative_abundance = match.arg(relative_abundance,choices = c("to_selection","to_total"))

  filter = match.arg(filter,choices = c("all","custom","DE.filter"))

  if(missing(annotation)){
    annotation = all.vars(object@design)[length(all.vars(object@design))]
  }

  ######check the object belongs to whamInput or whamoutput

  if(!is(object,"WhamOutputBiobakery") & !is(object,"WhamOutput16s") & !is(object,"WhamOutputEBI"))
    stop("the input is not one of following classes : 'WhamOutputBiobakery','WhamOutput16s',or 'WhamOutputEBI'")

  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s")){
    count_mat = assay(object)
  }

  if(is(object,"WhamOutputEBI")){
    count_mat = object@taxaCount
  }

  plot_data.frame = count_mat
  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s") | is(object,"WhamOutputEBI")){
    DE.result = object@DE.result
  }

  if(is(DE.result,"WhamResult")){
    contrast = DE.result@contrast
    DE_index = colData(object)[[contrast[1]]]
    numerato_index = which(DE_index==contrast[2])
    denominator_index = which(DE_index==contrast[3])

    DE.result = DE.result@DE.result
    DE.result$name = rownames(DE.result)
    plot_data.frame = plot_data.frame[,c(numerato_index,denominator_index)]
  }

  if(missing(taxa.level) & filter == "DE.filter"){
    taxa.level =  unique(DE.result$Level)
    taxa.level = c("k","p","o","c","f","g","s","otu")[grep( as.character(unique(DE.result$Level)),c("kingdom","phylum","class","order","family","genus","species","OTU "))]
  }

  if(missing(taxa.level) & c(filter == "all"| filter == "custom")){stop("taxa.level is required")
  }

  taxa.level=match.arg(taxa.level,choices = c("k","p","o","c","f","g","s","otu"))
  #####check the taxa
  Taxa_Info = object@TaxaInfo
  taxaLEVEL = TaxaIden_fn(TaxaInfo = Taxa_Info,taxa.level = taxa.level)
  if(relative_abundance == "to_total"){
    plot_data.frame = plot_data.frame[,!(apply(plot_data.frame,2,sum) == 0)]
    plot_data.frame =  t(t(plot_data.frame)/apply(plot_data.frame,2,sum))
  }

  ########filter####################

  if(filter == "DE.filter"){
    #######only show the taxa or feature which passed the DE test
    barindex = grep(taxa.level,c("k","p","c","o","f","g","s","otu"))

    DEindex = grep(unique(DE.result$Level),c("kingdom","phylum","class","order","family","genus","species","OTU"))
    if(DEindex < barindex)
      stop("During the differential test, taxa are collapsed at level which is higher than the one declaimed in taxa.level, please set up higher taxa level or  set filter as other to use all the taxa")

    if(length(grep("we.ep|we.eBH|wi.ep|wi.eBH",names(DE.result)))!=0){
      We.eBH.cutoff = p.cutoff
      if(missing(effect.size.range)){
        effect.size.range = c(0,0)
      }
      if(effect.size.range[1] > effect.size.range[2])
        stop( "lower limit of effect size cutoff has to be samller than or equal to upper limit")
      DE.result = subset(DE.result, we.eBH < We.eBH.cutoff)
      DE.result = DE.result[DE.result$effect <= effect.size.range[1] | DE.result$effect >= effect.size.range[2],]
    }

    if(length(grep("kw.ep|kw.eBH|glm.ep|glm.eBH",names(DE.result)))!=0){
      if(!missing(effect.size.range)){
        stop("Statistic test are based on glm.module, it is not necessary to provide effect size")
      }
      kw.eBH.cutoff = p.cutoff
      DE.result$name = rownames(DE.result)
      DE.result = subset(DE.result, kw.eBH < kw.eBH.cutoff)
    }

    DE.level = as.character(unique(DE.result$Level))

    ###collapse taxa to DE test level
    taxalevel.DE = apply(object@TaxaTable[,1:grep(DE.level,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel.DE = gsub("unassigned;|;unassigned","",taxalevel.DE)
    count_merged.DE = aggregate(plot_data.frame,by = list(taxalevel.DE),sum)
    rownames(count_merged.DE) = count_merged.DE[,1]
    CollTaxaTable.DE = TaxaToDF_fn(count_merged.DE[,1])
    rownames(CollTaxaTable.DE) = count_merged.DE[,1]
    rownames(count_merged.DE) = count_merged.DE[,1]
    count_merged.DE = count_merged.DE[,-1]
    plot_data.frame.DE = count_merged.DE[unique(DE.result$name),]

    ####collapse taxa to provided level
    taxatable.DE = TaxaToDF_fn(row.names(plot_data.frame.DE))
    if(taxaLEVEL == "OTU"){
      taxatable.DE = TaxaToDF_fn(row.names(plot_data.frame.DE),level_table ="otu")
    }
    taxalevel.DE = apply(taxatable.DE[,1:grep(taxaLEVEL,names(taxatable.DE))],1,paste,collapse =";")
    taxalevel.DE = gsub("unassigned;|;unassigned","",taxalevel.DE)
    count_merged = aggregate(plot_data.frame.DE,by = list(taxalevel.DE),sum)
    rownames(count_merged) = count_merged[,1]
    CollTaxaTable = TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) = count_merged[,1]
    count_merged = count_merged[,-1]
    plot_data.frame = count_merged

    ######the last assigned level as name
    name = regmatches(rownames(CollTaxaTable), regexpr("[^;]+$",rownames(CollTaxaTable)))
    unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__"),name)]
    name[unassigned_index] = paste0(name[unassigned_index],";unassigned")

    if(taxaLEVEL == "OTU"){
      names_coll = rownames(CollTaxaTable)
      names_coll_no_otu = gsub( paste0(";","OTU","_",".*$"),"",names_coll)
      names_coll_otu = gsub( paste0("^.*;OTU_"),"",names_coll)
      names_coll_otu = paste0("OTU_",names_coll_otu )
      name = regmatches(names_coll_no_otu, regexpr("[^;]+$",names_coll_no_otu))
      unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__"),name)]
      name[unassigned_index] = paste0(name[unassigned_index],";unassigned")
      name = paste0(name,"[",names_coll_otu,"]")
    }
  }

  if(filter == "all"){
    taxalevel = apply(object@TaxaTable[,1:grep(taxaLEVEL,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel = gsub("unassigned;|;unassigned","",taxalevel)
    count_merged = aggregate(plot_data.frame,by = list(taxalevel),sum)
    rownames(count_merged) = count_merged[,1]
    CollTaxaTable = TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) = count_merged[,1]
    count_merged = count_merged[,-1]
    plot_data.frame = count_merged

    ######the last assigned level as name#####
    name = regmatches(rownames(CollTaxaTable), regexpr("[^;]+$",rownames(CollTaxaTable)))
    unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__|unassigned"),name)]
    name[unassigned_index] = paste0(name[unassigned_index],";unassigned")

    if(taxaLEVEL == "OTU"){
      names_coll = rownames(CollTaxaTable)
      names_coll_no_otu = gsub( paste0(";","OTU","_",".*$"),"",names_coll)
      names_coll_otu = gsub( paste0("^.*;OTU_"),"",names_coll)
      names_coll_otu = paste0("OTU_",names_coll_otu )
      name = regmatches(names_coll_no_otu, regexpr("[^;]+$",names_coll_no_otu))
      unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__"),name)]
      name[unassigned_index] = paste0(name[unassigned_index],";unassigned")
      name = paste0(name,"[",names_coll_otu,"]")
    }

    ######select top # bug#####
    ra_plot_df = t(t(plot_data.frame)/apply(plot_data.frame,2,sum))
    index_abundance = rowMeans(ra_plot_df) > abundance_cutoff
    plot_data.frame = plot_data.frame[index_abundance,]
    index_variance = order(rowVars(as.matrix(plot_data.frame)),decreasing = T)
    plot_data.frame = plot_data.frame[index_variance,]
    if(display.number > nrow(plot_data.frame)){display.number = nrow(plot_data.frame) }
    plot_data.frame = plot_data.frame[1:display.number,]

    name = name[index_abundance]
    name = name[index_variance]
    name = name[1:display.number]
  }

  if(filter == "custom"){
    taxalevel = apply(object@TaxaTable[,1:grep(taxaLEVEL,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel = gsub("unassigned;|;unassigned","",taxalevel)
    count_merged = aggregate(plot_data.frame,by = list(taxalevel),sum)
    rownames(count_merged) = count_merged[,1]
    CollTaxaTable = TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) = count_merged[,1]
    count_merged = count_merged[,-1]
    custom = paste0(taxa.level,"__",custom)
    index_custom = grep(paste0(custom,collapse = "|"),rownames(count_merged))
    plot_data.frame = count_merged[index_custom,]
    ######the last assigned level as name#####
    name = regmatches(rownames(CollTaxaTable), regexpr("[^;]+$",rownames(CollTaxaTable)))
    unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__|unassigned"),name)]
    name[unassigned_index] = paste0(name[unassigned_index],";unassigned")
    name = name[index_custom]

    if(taxaLEVEL == "OTU"){
      names_coll = rownames(CollTaxaTable)
      names_coll_no_otu = gsub( paste0(";","OTU","_",".*$"),"",names_coll)
      names_coll_otu = gsub( paste0("^.*;OTU_"),"",names_coll)
      names_coll_otu = paste0("OTU_",names_coll_otu )
      name = regmatches(names_coll_no_otu, regexpr("[^;]+$",names_coll_no_otu))
      unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__"),name)]
      name[unassigned_index] = paste0(name[unassigned_index],";unassigned")
      name = paste0(name,"[",names_coll_otu,"]")
      name = name[index_custom]
    }
  }

  ####create data.frame and plot of barplot
  if(relative_abundance == "to_selection"){
    plot_data.frame = plot_data.frame[,!(apply(plot_data.frame,2,sum) == 0)]
    plot_data.frame = t(t(plot_data.frame)/apply(plot_data.frame,2,sum))
  }
  barplot.df =  t(plot_data.frame)
  barplot.df = as.data.frame(barplot.df)

  barplot.df = cbind(as.data.frame(colData(object))[colnames(plot_data.frame),annotation,drop=F],barplot.df)
  barplot.df$sampleID = rownames(barplot.df)
  barplot.df = barplot.df[,c(ncol(barplot.df),1:(ncol(barplot.df)-1))]
  if(merge_group == F){
    barplot.df = barplot.df[order(barplot.df[[2]]),]
  }

  if(merge_group == T){
    #barplot.df = barplot.df[!apply(is.nan(as.matrix(barplot.df[,3:ncol(barplot.df)])),1,any),]
    barplot.df = aggregate(barplot.df[,3:ncol(barplot.df)],by = list(barplot.df[[2]]),mean)
    barplot.df$sampleID = barplot.df$Group.1
    barplot.df = barplot.df[,c(ncol(barplot.df),1:c(ncol(barplot.df)-1))]
    names(barplot.df)[2] = c("group")
  }

  ggplot_df = barplot.df  %>% gather("bacteria","RPK",3:(ncol(barplot.df)))
  names(ggplot_df) = c("sampleID","annotation","bacteria","RPK")
  ggplot_df$sampleID = factor(ggplot_df$sampleID,levels = (barplot.df$sampleID))

  if(relative_abundance == "to_total"){
    y_label  = "relative abundance"
  }

  if(relative_abundance == "to_selection"){
    y_label  = "relative abundance to selection"
  }

  p1 = ggplot(ggplot_df,aes(x = sampleID,y = RPK,fill = bacteria))+
    geom_bar(stat = "identity")+
    theme_cowplot(9)+
    rotate_x_text(90)+
    ylab(y_label)+
    theme(legend.key.size = unit(0.4, 'cm'),plot.margin = unit(c(0, 0, 0, 0), "cm"))+scale_fill_discrete(guide = guide_legend(nrow = 12),labels = name)

  p2 = ggplot(ggplot_df,aes(x=sampleID,y=1,fill= annotation))+
    geom_bar(stat = "identity",width = 1)+
    theme_void()+theme(legend.direction= "vertical")+
    theme(legend.key.size = unit(0.4, 'cm'),plot.margin = unit(c(0, 0, 0, 0), "cm"))+scale_fill_discrete(guide = guide_legend(nrow = 6))

  l1 = get_legend(p1)
  l2 = get_legend(p2)
  l2$layout$l = c(1,1)
  l1$layout$l = c(1,1)
  l1$layout$t = c(1,1)
  l2$layout$b = c(1,1)
  l1$layout$b = c(2,1)

  legend <- plot_grid(l2, l1, ncol = 1,align = "v")+theme(plot.margin = unit(c(0, 1, 0, 0.3), "cm"))

  p1 = ggplot(ggplot_df,aes(x = sampleID,y = RPK,fill = bacteria))+
    geom_bar(stat = "identity")+
    theme_cowplot(9)+
    rotate_x_text(90)+
    ylab(y_label)+
    scale_fill_discrete(labels = name)+ theme(legend.position='none',plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"))

  p2 = ggplot(ggplot_df,aes(x=sampleID,y=1,fill= annotation))+
    geom_bar(stat = "identity",width = 1)+
    theme_void()+ theme(legend.position='none',plot.margin = unit(c(0, 0, 0, 0), "cm"))

  plot = plot_grid(p2,p1,align = "v",axis = "tb",ncol = 1,rel_heights = c(0.5, 15))+theme(plot.margin = unit(c(0.3,0,0.3,0.3), "cm"))
  plot_legend = plot_grid(plot, legend, nrow = 1, align = "h",axis = "rl")+ theme(plot.margin = unit(c(0,0,0,0), "cm"))

  print(plot_legend)
}

#' TaxaHeatmap
#' @description Visualization of Wham upstream analysis
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr  ggpubr  methods ComplexHeatmap circlize
#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom grid gpar
#' @param object required. The output of WhamBiobakery, WhamEBI or Wham16s
#' @param taxa.level Required, c("k","p","o","c","f","g","s").Collapse taxa at the provided level.
#' @param filter c("all","custom","DE.filter"). Default is "all". "all" shows all samples. "custom" shows taxa provided by user and all samples. "DE.filter" shows taxa which pass DE test cutoff(p.cutoff and effect.size.range) and samples which are in the DE.result contrast.
#' @param custom required when filter is set as "custom". A vector which contain taxa to be shown
#' @param p.cutoff applied when filter is set as "DE.filter". Default is 0.05.
#' @param effect.size.range applied when filter is set as "DE.filter" and aldex.module is set as t.test.Default range is c(0,0)
#' @param scale default is TRUE. Z score is used when set as T. Expression counts are used when set as F
#' @param column_split Split on columns.TRUE, FALSE or a vector provided by users.Default value is False. Column will be spilt based on annotation.
#' @param annotation optional. Default is the last variable in design. A character vector indicates what column of metadata will be used for annotation bar
#' @param col an argument in ComplexHeatmap::Heatmap: vector of colors if the color mapping is discrete or a color mapping function if the matrix is continuous numbers (should be generated by colorRamp2). If the matrix is continuous, the value can also be a vector of colors so that colors can be interpolated.
#' @param count c("rawCount","WhamTransformation"). Default value is "rawCount". Setting as "rawCount" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#' @param ... arguments in ComplexHeatmap::Heatmap
#' @return a heatmap of taxonomy
#' @export
#'
#' @examples
#' # see the vignette for more details
TaxaHeatmap = function(object,
                       taxa.level = c("k","p","o","c","f","g","s","otu"),
                       filter = c("all","custom","DE.filter"),
                       custom,
                       p.cutoff=0.05,
                       effect.size.range = c(0,0),
                       scale = T,
                       col,
                       column_split =F,
                       annotation,
                       count = c("rawCount","WhamTransformation"),
                       ...){
  ###parameter confirmation
  if(missing(taxa.level))
    stop("taxa.level is required")
  taxa.level=match.arg(taxa.level,choices = c("k","p","o","c","f","g","s","otu"))

  count = match.arg(count,choices = c("rawCount","WhamTransformation"))


  if(missing(filter)){filter = "all"}
  filter = match.arg(filter,choices = c("all","custom","DE.filter"))

  if(missing(annotation)){
    annotation = all.vars(object@design)[length(all.vars(object@design))]
  }

  ######check the object belongs to whamInput or whamoutput

  if(!is(object,"WhamOutputBiobakery") & !is(object,"WhamOutput16s") & !is(object,"WhamOutputEBI"))
    stop("the input is not one of following classes : 'WhamOutputBiobakery','WhamOutput16s',or 'WhamOutputEBI'")

  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s")){
    count_mat = assay(object)
  }

  if(is(object,"WhamOutputEBI")){
    count_mat = object@taxaCount
  }

  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s") | is(object,"WhamOutputEBI")){
    DE.result = object@DE.result
  }

  #####check the taxa
  Taxa_Info = object@TaxaInfo
  taxaLEVEL = TaxaIden_fn(TaxaInfo = Taxa_Info,taxa.level = taxa.level)
  plot_data.frame = count_mat
  if(count == "WhamTransformation"){
    plot_data.frame = t(t(plot_data.frame)/apply(plot_data.frame,2,sum))*1e6
  }

  if(is(DE.result,"WhamResult")){
    contrast = DE.result@contrast
    DE_index = colData(object)[[contrast[1]]]
    numerato_index = which(DE_index==contrast[2])
    denominator_index = which(DE_index==contrast[3])

    DE.result = DE.result@DE.result
    DE.result$name = rownames(DE.result)
    plot_data.frame = plot_data.frame[,c(numerato_index,denominator_index)]
  }

  ########filter####################

  if(filter == "DE.filter"){
    #######only show the taxa or feature which passed the DE test
    barindex = grep(taxa.level,c("k","p","c","o","f","g","s","otu"))

    DEindex = grep(unique(DE.result$Level),c("kingdom","phylum","class","order","family","genus","species","OTU"))
    if(DEindex < barindex)
      stop("During the differential test, taxa are collapsed at level which is higher than the one declaimed in taxa.level, please set up higher taxa level or  set filter as other to use all the taxa")

    if(length(grep("we.ep|we.eBH|wi.ep|wi.eBH",names(DE.result)))!=0){
      We.eBH.cutoff = p.cutoff
      if(missing(effect.size.range)){
        effect.size.range = c(0,0)
      }
      if(effect.size.range[1] > effect.size.range[2])
        stop( "lower limit of effect size cutoff has to be samller than or equal to upper limit")
      DE.result = subset(DE.result, we.eBH < We.eBH.cutoff)
      DE.result = DE.result[DE.result$effect <= effect.size.range[1] | DE.result$effect >= effect.size.range[2],]
    }

    if(length(grep("kw.ep|kw.eBH|glm.ep|glm.eBH",names(DE.result)))!=0){
      if(!missing(effect.size.range)){
        stop("Statistic test are based on glm.module, it is not necessary to provide effect size")
      }
      kw.eBH.cutoff = p.cutoff
      DE.result$name = rownames(DE.result)
      DE.result = subset(DE.result, kw.eBH < kw.eBH.cutoff)
    }

    DE.level = as.character(unique(DE.result$Level))

    ###collapse taxa to DE test level
    taxalevel.DE = apply(object@TaxaTable[,1:grep(DE.level,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel.DE = gsub("unassigned;|;unassigned","",taxalevel.DE)
    count_merged.DE = aggregate(plot_data.frame,by = list(taxalevel.DE),sum)
    rownames(count_merged.DE) = count_merged.DE[,1]
    CollTaxaTable.DE = TaxaToDF_fn(count_merged.DE[,1])
    rownames(CollTaxaTable.DE) = count_merged.DE[,1]
    rownames(count_merged.DE) = count_merged.DE[,1]
    count_merged.DE = count_merged.DE[,-1]
    plot_data.frame.DE = count_merged.DE[unique(DE.result$name),]

    ####collapse taxa to provided level
    taxatable.DE = TaxaToDF_fn(row.names(plot_data.frame.DE))
    if(taxaLEVEL == "OTU"){
      taxatable.DE = TaxaToDF_fn(row.names(plot_data.frame.DE),level_table ="otu")
    }
    taxalevel.DE = apply(taxatable.DE[,1:grep(taxaLEVEL,names(taxatable.DE))],1,paste,collapse =";")
    taxalevel.DE = gsub("unassigned;|;unassigned","",taxalevel.DE)
    count_merged = aggregate(plot_data.frame.DE,by = list(taxalevel.DE),sum)
    rownames(count_merged) = count_merged[,1]
    CollTaxaTable = TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) = count_merged[,1]
    count_merged = count_merged[,-1]
    plot_data.frame = count_merged

    ######the last assigned level as name
    name = regmatches(rownames(CollTaxaTable), regexpr("[^;]+$",rownames(CollTaxaTable)))
    unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__"),name)]
    name[unassigned_index] = paste0(name[unassigned_index],";unassigned")

    if(taxaLEVEL == "OTU"){
      names_coll = rownames(CollTaxaTable)
      names_coll_no_otu = gsub( paste0(";","OTU","_",".*$"),"",names_coll)
      names_coll_otu = gsub( paste0("^.*;OTU_"),"",names_coll)
      names_coll_otu = paste0("OTU_",names_coll_otu )
      name = regmatches(names_coll_no_otu, regexpr("[^;]+$",names_coll_no_otu))
      unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__"),name)]
      name[unassigned_index] = paste0(name[unassigned_index],";unassigned")
      name = paste0(name,"[",names_coll_otu,"]")
    }
  }

  if(filter == "all"){
    taxalevel = apply(object@TaxaTable[,1:grep(taxaLEVEL,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel = gsub("unassigned;|;unassigned","",taxalevel)
    count_merged = aggregate(plot_data.frame,by = list(taxalevel),sum)
    rownames(count_merged) = count_merged[,1]
    CollTaxaTable = TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) = count_merged[,1]
    count_merged = count_merged[,-1]
    plot_data.frame = count_merged
    ######the last assigned level as name#####
    name = regmatches(rownames(CollTaxaTable), regexpr("[^;]+$",rownames(CollTaxaTable)))
    unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__|unassigned"),name)]
    name[unassigned_index] = paste0(name[unassigned_index],";unassigned")

    if(taxaLEVEL == "OTU"){
      names_coll = rownames(CollTaxaTable)
      names_coll_no_otu = gsub( paste0(";","OTU","_",".*$"),"",names_coll)
      names_coll_otu = gsub( paste0("^.*;OTU_"),"",names_coll)
      names_coll_otu = paste0("OTU_",names_coll_otu )
      name = regmatches(names_coll_no_otu, regexpr("[^;]+$",names_coll_no_otu))
      unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__"),name)]
      name[unassigned_index] = paste0(name[unassigned_index],";unassigned")
      name = paste0(name,"[",names_coll_otu,"]")
    }
  }

  if(filter == "custom"){
    taxalevel = apply(object@TaxaTable[,1:grep(taxaLEVEL,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel = gsub("unassigned;|;unassigned","",taxalevel)
    count_merged = aggregate(plot_data.frame,by = list(taxalevel),sum)
    rownames(count_merged) = count_merged[,1]
    CollTaxaTable = TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) = count_merged[,1]
    count_merged = count_merged[,-1]
    custom = paste0(taxa.level,"__",custom)
    index_custom = grep(paste(custom,collapse = "|"),rownames(count_merged))
    plot_data.frame = count_merged[index_custom,]
    ######the last assigned level as name#####
    name = regmatches(rownames(CollTaxaTable), regexpr("[^;]+$",rownames(CollTaxaTable)))
    unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__|unassigned"),name)]
    name[unassigned_index] = paste0(name[unassigned_index],";unassigned")
    name = name[index_custom]

    if(taxaLEVEL == "OTU"){
      names_coll = rownames(CollTaxaTable)
      names_coll_no_otu = gsub( paste0(";","OTU","_",".*$"),"",names_coll)
      names_coll_otu = gsub( paste0("^.*;OTU_"),"",names_coll)
      names_coll_otu = paste0("OTU_",names_coll_otu )
      name = regmatches(names_coll_no_otu, regexpr("[^;]+$",names_coll_no_otu))
      unassigned_index = c(1:length(name))[-grep(paste0(taxa.level,"__"),name)]
      name[unassigned_index] = paste0(name[unassigned_index],";unassigned")
      name = paste0(name,"[",names_coll_otu,"]")
      name = name[index_custom]
    }
  }


  if(!is.logical(column_split)){
    column.split = column_split
  }
  if(is.logical(column_split)){
    if(column_split == T){
      if(missing(annotation)){
        message("annotation is not provided, the last variable in design will be used as annotation")
        annotation = all.vars(object@design)[length(all.vars(object@design))]
      }
      column.split = colData(object)[colnames(plot_data.frame),annotation]
    }
    if(column_split == F){
      column.split = NULL
    }
  }

  if(scale == T){
    plot_data.frame = log2(plot_data.frame+1)
    plot_data.frame = as.matrix(plot_data.frame)
    plot_data.frame = t(scale(t(plot_data.frame)))
    index = !apply(is.na(plot_data.frame),1,any)
    plot_data.frame = plot_data.frame[!apply(is.na(plot_data.frame),1,any),]
    name = name[index]
    if(missing(col)){
      col =  colorRamp2(c(quantile(plot_data.frame,probs = 0.02), 0, quantile(plot_data.frame,probs = 0.98)), c("blue", "white", "red"))
    }
    ht = Heatmap(plot_data.frame,
                 heatmap_legend_param = list(direction = "horizontal",
                                             title = "z-score",
                                             title_position = "topcenter"),
                 col = col,
                 column_split = column.split,
                 row_labels = name,
                 ...)
    draw(ht,heatmap_legend_side = "top")
  }

  if(scale == F){
    plot_data.frame = log10(plot_data.frame+1)
    plot_data.frame = as.matrix(plot_data.frame)
    index = !apply(is.na(plot_data.frame),1,any)
    plot_data.frame = plot_data.frame[!apply(is.na(plot_data.frame),1,any),]
    name = name[index]
    if(missing(col)){
      col =  colorRamp2(c(0, quantile(plot_data.frame,probs = 0.98)), c("white", "red"))
    }
    ht = Heatmap(plot_data.frame,
                 heatmap_legend_param = list(direction = "horizontal",
                                             title = "expression",
                                             title_position = "topcenter"),
                 col = col,
                 column_split = column.split,
                 row_labels = name,
                 ...)
    ht@matrix_legend_param$labels = 10^(ht@matrix_color_mapping@levels)
    draw(ht,heatmap_legend_side = "top")
  }
}


#' FeatureHeatmap
#' @description Visualization of Wham upstream analysis
#' @importFrom stats aggregate
#' @importFrom stats quantile
#' @importFrom grid gpar
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr ggpubr  methods ComplexHeatmap circlize
#' @param object required. The output of WhamBiobakery, WhamEBI or Wham16s
#' @param filter c("all","custom","DE.filter"). Default is "all". "all" shows all samples. "custom" shows taxa provided by user and all samples. "DE.filter" shows taxa which pass DE test cutoff(p.cutoff and effect.size.range) and samples which are in the DE.result contrast.
#' @param custom required when filter is set as "custom". A vector which contain taxa to be shown
#' @param p.cutoff applied when filter is set as "DE.filter". Default is 0.05.
#' @param effect.size.range applied when filter is set as "DE.filter" and aldex.module is set as t.test.Default range is c(0,0)
#' @param scale default is TRUE. Z score is used when set as T. Expression counts are used when set as F
#' @param col an argument in ComplexHeatmap::Heatmap: vector of colors if the color mapping is discrete or a color mapping function if the matrix is continuous numbers (should be generated by colorRamp2). If the matrix is continuous, the value can also be a vector of colors so that colors can be interpolated.
#' @param count c("rawCount","WhamTransformation"). Default value is "rawCount". Setting as "rawCount" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#' @param column_split Split on columns.TRUE, FALSE or a vector provided by users.Default value is False. Column will be spilt based on annotation.
#' @param annotation optional. Default is the last variable in design. A character vector indicates what column of metadata will be used for annotation bar
#' @param ... arguments in ComplexHeatmap::Heatmap
#' @return a heatmap of Feature
#' @export
#'
#' @examples
#' # see the vignette for more details
FeatureHeatmap=function(object,
                        filter = c("all","custom","DE.filter"),
                        custom,
                        p.cutoff=0.05,
                        effect.size.range = c(0,0),
                        scale =T,
                        column_split =F,
                        annotation,
                        col,
                        count = c("rawCount","WhamTransformation"),
                        ...){
  ######check the object belongs to whamInput or whamoutput
  count = match.arg(count,choices = c("rawCount","WhamTransformation"))

  if(missing(filter)){filter = "all"}
  filter = match.arg(filter,choices = c("all","custom","DE.filter"))
  ######check the object belongs to whamInput or whamoutput
  if(!is(object,"WhamOutputBiobakery") & !is(object,"WhamOutput16s") & !is(object,"WhamOutputEBI"))
    stop("the input is not one of following classes : 'WhamOutputBiobakery','WhamOutput16s',or 'WhamOutputEBI'")
  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s")){
    count_mat = assay(object)
  }
  if(is(object,"WhamOutputEBI")){
    count_mat = object@featureCount
  }
  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s") | is(object,"WhamOutputEBI")){
    DE.result = object@DE.result
  }
  plot_data.frame = count_mat
  if(count == "WhamTransformation"){
    plot_data.frame = t(t(plot_data.frame)/apply(plot_data.frame,2,sum))*1e6
  }

  if(is(DE.result,"WhamResult")){
    contrast = DE.result@contrast
    DE_index = colData(object)[[contrast[1]]]
    numerato_index = which(DE_index==contrast[2])
    denominator_index = which(DE_index==contrast[3])
    plot_data.frame  = plot_data.frame[,c(numerato_index,denominator_index)]
    DE.result = DE.result@DE.result
    DE.result$name = rownames(DE.result)
  }
  ################
  plot_data.frame = aggregate(plot_data.frame,by = list(object@FeatureInfo),sum)
  rownames(plot_data.frame) = plot_data.frame[,1]
  plot_data.frame = plot_data.frame[,-1]

  if(filter == "DE.filter"){
    ##### only show the sample appearing in the declaimed contrast

    #######only show the taxa or feature which passed the DE test
    if(length(grep("we.ep|we.eBH|wi.ep|wi.eBH",names(DE.result)))!=0){
      We.eBH.cutoff = p.cutoff
      if(missing(effect.size.range)){ effect.size.range = c(0,0)}
      if(effect.size.range[1] > effect.size.range[2])
        stop( "lower limit of effect size cutoff has to be samller than or equal to upper limit")
      DE.result = subset(DE.result, we.eBH < We.eBH.cutoff)
      DE.result = DE.result[DE.result$effect <= effect.size.range[1] | DE.result$effect >= effect.size.range[2],]
    }
    if(length(grep("kw.ep|kw.eBH|glm.ep|glm.eBH",names(DE.result)))!=0){
      if(!missing(effect.size.range)){
        stop("Statistic test are based on glm.module, it is not necessary to provide effect size")
      }
      kw.eBH.cutoff = p.cutoff
      DE.result$name = rownames(DE.result)
      DE.result = subset(DE.result, kw.eBH < kw.eBH.cutoff)
    }
    plot_data.frame = plot_data.frame[unique(DE.result$name),]
  }
  if(filter == "custom"){
    plot_data.frame = plot_data.frame[custom,]
  }

  if(!is.logical(column_split)){
    column.split = column_split
  }
  if(is.logical(column_split)){
    if(column_split == T){
      if(missing(annotation)){
        message("annotation is not provided, the last variable in design will be used as annotation")
        annotation = all.vars(object@design)[length(all.vars(object@design))]
      }
      column.split = colData(object)[colnames(plot_data.frame),annotation]
    }
    if(column_split == F){
      column.split = NULL
    }
  }


  if(scale == T){
    plot_data.frame = log2(plot_data.frame+1)
    plot_data.frame = as.matrix(plot_data.frame)
    plot_data.frame = t(scale(t(plot_data.frame)))
    plot_data.frame = plot_data.frame[!apply(is.na(plot_data.frame),1,any),]
    if(missing(col)){
      col =  colorRamp2(c(quantile(plot_data.frame,probs = 0.05), 0, quantile(plot_data.frame,probs = 0.95)), c("blue", "white", "red"))
    }
    ht = Heatmap(plot_data.frame,
                 heatmap_legend_param = list(direction = "horizontal",
                                             title = "z-score",
                                             title_position = "topcenter"),
                 col = col,
                 column_split = column.split,
                 ...)
    draw(ht,heatmap_legend_side = "top")
  }

  if(scale == F){
    plot_data.frame = log10(plot_data.frame+1)
    plot_data.frame = as.matrix(plot_data.frame)
    plot_data.frame = plot_data.frame[!apply(is.na(plot_data.frame),1,any),]
    if(missing(col)){
      col =  colorRamp2(c(0, quantile(plot_data.frame,probs = 0.95)), c("white", "red"))
    }
    ht = Heatmap(plot_data.frame,
                 heatmap_legend_param = list(direction = "horizontal",
                                             title = "expression",
                                             title_position = "topcenter"),
                 col = col,
                 column_split = column.split,
                 ...)
    ht@matrix_legend_param$labels = 10^(ht@matrix_color_mapping@levels)
    draw(ht,heatmap_legend_side = "top")
  }
}

###correlaton function
#' feature correlation
#' @description Visualization of Wham upstream analysis
#' @importFrom stats aggregate p.adjust
#' @importFrom grid grid.text
#' @importFrom matrixStats rowVars
#' @importFrom psych corr.test
#' @param object required. The output of WhamBiobakery, WhamEBI or Wham16s
#' @param feature default is NULL. if not provided and UseDE.result is FALSE, all features will be used.
#' @param group default is NULL. if not provided and UseDE.result is FALSE, all samples will be used
#' @param UseDE.result default value is FALSE. FALSE means using the argument feature and argument group. TRUE means using the features and samples from DE.result (after selection of parameters below).
#' @param p.cutoff Keep genes whose we.eHB values are less than cutoff. The default value is 0.05. Use all gene from DE.result
#' @param effect.size.range Keep genes whose effect size values are out of the range.The default value is c(0,0), Use all gene from DE.result. This only work when the aldex.module is set as t.test
#' @param method  the method of correlation. "spearman","pearson", and "kendall"
#' @param adjust the method of multiple test correction. "holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"
#' @param count c("rawCount","WhamTransformation"). Default value is "rawCount". Setting as "rawCount" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#'
#' @return  a correlation heatmap of gene family
#' @export
#'
#' @examples
#'
#' # see the vignette for more details
FeatureCorr <- function(object,
                        feature,
                        group,
                        UseDE.result=F,
                        p.cutoff=0.05,
                        effect.size.range = c(0,0),
                        method = c("pearson","spearman","kendall"),
                        adjust = c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"),
                        count = c("rawCount","WhamTransformation")
){if(!is(object,"WhamOutputBiobakery") & !is(object,"WhamOutput16s") & !is(object,"WhamOutputEBI"))
    stop("the input is not one of following classes : 'WhamOutputBiobakery','WhamOutput16s',or 'WhamOutputEBI'")
  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s")){
    count_mat = assay(object)
  }
  if(is(object,"WhamOutputEBI")){
    count_mat = object@featureCount
  }
  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s") | is(object,"WhamOutputEBI")){
    DE.result = object@DE.result
  }
  if(missing(method)){
    method = "spearman"
  }else{method = match.arg(method,choices = c("pearson","spearman","kendall"))}
  if(missing(adjust)){
    adjust = "BH"
  }else{adjust =  match.arg(adjust,choices = c("holm","hochberg","hommel","bonferroni","BH","BY","fdr","none"))}
  count = match.arg(count,choices = c("rawCount","WhamTransformation"))
  ##########correlation test function buildup#######
  corr_fn = function(count_matrix,method = method,adjust = adjust){
    message(paste0("features:(totally ",sum(rowVars(as.matrix(count_matrix)) == 0),")"," ",paste(rownames(count_matrix)[rowVars(as.matrix(count_matrix)) == 0],collapse = ",")," is/are removed because the standard deviations of these feature(s) are Zero"))
    count_matrix = count_matrix[rowVars(as.matrix(count_matrix)) != 0,]
    res = corr.test(t(count_matrix),method = method,adjust = "none")
    r_value_matrix = res$r
    p_value_matrix = res$p
    p_value_matrix[col(p_value_matrix) == row(p_value_matrix) | upper.tri(p_value_matrix)] = NA
    p_value_matrix = p.adjust( p_value_matrix, method = adjust, n = length(p_value_matrix))
    p_value_matrix = matrix(p_value_matrix, ncol = nrow(count_matrix) , nrow = nrow(count_matrix))
    p_value_matrix[upper.tri(p_value_matrix)] = t(p_value_matrix)[upper.tri(p_value_matrix)]
    diag(p_value_matrix) = rep(1,nrow(p_value_matrix))
    p_value_matrix[p_value_matrix < 0.001] = "***"
    p_value_matrix[p_value_matrix >= 0.001 & p_value_matrix < 0.01] = "**"
    p_value_matrix[p_value_matrix >= 0.01 & p_value_matrix < 0.05] = "*"
    p_value_matrix[p_value_matrix >= 0.05] = ""
    p_value_matrix[is.na(p_value_matrix)] = ""
    ht = Heatmap(r_value_matrix,
                 heatmap_legend_param = list(direction = "horizontal",
                                             title = "r",
                                             title_position = "topcenter"),
                 col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(p_value_matrix[i,j], x, y)})
    res = draw(ht,heatmap_legend_side = "top")
    return(res)
  }
  count_matrix = count_mat
  if(count == "WhamTransformation"){
    count_matrix = t(t(count_matrix)/apply(count_matrix,2,sum))*1e6
  }
  count_matrix = aggregate(count_matrix ,by =list(object@FeatureInfo),sum)
  rownames(count_matrix) =count_matrix[,1]
  count_matrix = count_matrix[,-1]
  if(!missing(feature) & UseDE.result == T){
    stop("When the argument feature is declaimed, UseDE.result can't be set as TRUE")
  }
  if(!missing(group) & UseDE.result == T){
    stop("When the argument group is declaimed, UseDE.result can't be set as TRUE")
  }
  if(missing(feature) & missing(group) & UseDE.result == F){
    message("Arguments 'group' and 'feature' are not declaimed and at the same time 'UseDE.result' are set as FALSE, so all samples and all features will be in process of correlation test (Usually, this setup will cost a long time\n")
    res = corr_fn(count_matrix = count_matrix, method = method, adjust = adjust)
  }
  if(missing(feature) & missing(group) & UseDE.result == T){
    message("Correlation test will be based on the samples and features in the DE.result. (The we.eHB cutoff and effect size cutoff can be adjusted via parameter we.eHB.cutoff and effect.size.range)\n")
    if(is(DE.result,"data.frame")){
      message("All samples will be used for correlation\n")
      if(length(grep("we.ep|we.eBH|wi.ep|wi.eBH",names(DE.result)))!=0){
        We.eBH.cutoff = p.cutoff
        if(missing(effect.size.range)){effect.size.range = c(0,0)}
        if(effect.size.range[1] > effect.size.range[2])
          stop( "lower limit of effect size cutoff has to be samller than or equal to upper limit")
        DE.result = subset(DE.result, we.eBH < We.eBH.cutoff)
        DE.result = DE.result[DE.result$effect <= effect.size.range[1] | DE.result$effect >= effect.size.range[2],]
      }
      if(length(grep("kw.ep|kw.eBH|glm.ep|glm.eBH",names(DE.result)))!=0){
        if(!missing(effect.size.range)){
          stop("Statistic test are based on glm.module, it is not necessary to provide effect size")
        }
        kw.eBH.cutoff = p.cutoff
        DE.result$name = rownames(DE.result)
        DE.result = subset(DE.result, kw.eBH < kw.eBH.cutoff)
      }
      count_matrix = count_matrix[DE.result$name,]
      res = corr_fn(count_matrix = count_matrix, method = method, adjust = adjust)
    }
    if(is(DE.result,"WhamResult")){
      contrast = DE.result@contrast
      DE_index = colData(object)[[contrast[1]]]
      numerato_index = which(DE_index==contrast[2])
      denominator_index = which(DE_index==contrast[3])
      count_matrix  =count_matrix[,c(numerato_index,denominator_index)]
      DE.result = DE.result@DE.result
      DE.result$name = rownames(DE.result)
      if(length(grep("we.ep|we.eBH|wi.ep|wi.eBH",names(DE.result)))!=0){
        We.eBH.cutoff = p.cutoff
        if(missing(effect.size.range)){effect.size.range = c(0,0)}
        if(effect.size.range[1] > effect.size.range[2])
          stop( "lower limit of effect size cutoff has to be samller than or equal to upper limit")
        DE.result = subset(DE.result, we.eBH < We.eBH.cutoff)
        DE.result = DE.result[DE.result$effect <= effect.size.range[1] | DE.result$effect >= effect.size.range[2],]
      }
      if(length(grep("kw.ep|kw.eBH|glm.ep|glm.eBH",names(DE.result)))!=0){
        if(!missing(effect.size.range)){
          stop("Statistic test are based on glm.module, it is not necessary to provide effect size")
        }
        kw.eBH.cutoff = p.cutoff
        DE.result$name = rownames(DE.result)
        DE.result = subset(DE.result, kw.eBH < kw.eBH.cutoff)
      }
      count_matrix = count_matrix[DE.result$name,]
      res = corr_fn(count_matrix = count_matrix, method = method, adjust = adjust)
    }
  }
  if(!missing(feature) & missing(group) & UseDE.result == F){
    message("Correlation test will be in process of calculation via using all samples and declaimed features\n")
    count_matrix = count_matrix[feature,]
    res = corr_fn(count_matrix = count_matrix, method = method, adjust = adjust)
  }
  if(!missing(feature) & !missing(group) & UseDE.result == F){
    message("Correlation test will be in process of calculation via using samples and  declaimed features\n")
    count_matrix = count_matrix[feature,group]
    res = corr_fn(count_matrix = count_matrix, method = method, adjust = adjust)
  }
  if(missing(feature) & !missing(group) & UseDE.result == F){
    message("Correlation test will be in process of calculation via using all feature and declaimed samples\n")
    count_matrix = count_matrix[feature,group]
    res = corr_fn(count_matrix = count_matrix, method = method, adjust = adjust)
  }
}


#' VolcanoPlot
#' @description Visualization of Wham upstream analysis
#' @import EnhancedVolcano
#' @param object required. The output of WhamBiobakery, WhamEBI or Wham16s
#' @param effect c("effect"), Use effect as effectcutoff. Default is "effect"
#' @param p_value c("we.eBH","we.ep","wi.eBH","wi.ep"). Using one of them as p-value.Default is "we.eBH".
#' @param pCutoff Cut-off for statistical significance. A horizontal line will be drawn at -log10(pCutoff). DEFAULT = 0.05. OPTIONAL
#' @param effectCutoff ut-off for absolute effect size. Vertical lines will be drawn at the negative and positive values of FCCutoff. DEFAULT = 1.0. OPTIONAL.
#' @param ... ...
#' @return a volcano plot of statistical result
#' @export
#'
#' @examples
#' # see the vignette for more details
VolcanoPlot = function(object = NULL,
                       effect = c("effect"),
                       p_value = c("we.eBH","we.ep","wi.eBH","wi.ep"),
                       pCutoff = 0.05,
                       effectCutoff = 1,
                       ...
){
  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s") | is(object,"WhamOutputEBI")){
    count_mat = assay(object)
  }else{
    stop("the input is not one of following classes : 'WhamOutputBiobakery','WhamOutput16s',or 'WhamOutputEBI'")
  }
  if(is(object,"WhamOutputBiobakery") | is(object,"WhamOutput16s") | is(object,"WhamOutputEBI")){
    DE.result = object@DE.result
  }
  effect = match.arg(effect,c("effect"))
  p_value = match.arg(p_value,c("we.eBH","we.ep","wi.eBH","wi.ep"))
  if(is(DE.result,"data.frame")){
    if(length(grep("kw.ep|kw.eBH|glm.ep|glm.eBH",names(DE.result)))!=0){
      stop("vocalnoplot can only be used for t.test module")
    }
    Comparison = unique(DE.result$COMPARISON)
    list_res =list()
    for( i in Comparison){
      VolcanoPlot_df = DE.result[DE.result$COMPARISON == i,]
      list_res[[i]] = EnhancedVolcano(VolcanoPlot_df,
                                      lab = VolcanoPlot_df$name,
                                      x = effect ,
                                      y = p_value ,
                                      xlab = effect,
                                      pCutoff = pCutoff,
                                      FCcutoff = effectCutoff,
                                      title = i,
                                      ylim = c(0, max(-log10(VolcanoPlot_df[,p_value]), na.rm=TRUE) *1.1),
                                      legendLabels  = c("NS", "Effect Size", "P", "P & Effect Size")
      )
    }
  }
  if(is(DE.result,"data.frame")){
    return(list_res)
    print(list_res)
  }
  if(is(DE.result,"WhamResult")){
    contrast = DE.result@contrast
    VolcanoPlot_df = DE.result@DE.result
    p1 = EnhancedVolcano(VolcanoPlot_df,
                         lab = rownames(VolcanoPlot_df),
                         x = effect ,
                         y = p_value ,
                         xlab = effect,
                         pCutoff = pCutoff,
                         FCcutoff = effectCutoff,
                         title = paste0(contrast[2]," vs ",contrast[[3]]),
                         ylim = c(0, max(-log10(VolcanoPlot_df[,p_value]), na.rm=TRUE) *1.1),
                         legendLabels  = c("NS", "Effect Size", "P", "P & Effect Size")
    )
  }
  if(is(DE.result,"WhamResult")){
    return(p1)
    print(p1)
  }
}


