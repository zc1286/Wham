

#############################################Wham Biobakery############
#'@importClassesFrom SummarizedExperiment RangedSummarizedExperiment

setClass("WhamInputfromBiobakery",
         contains = "RangedSummarizedExperiment",
         representation = representation(
           se = "RangedSummarizedExperiment",
           design = "ANY",
           TaxaInfo = "character",
           FeatureInfo = "character",
           TaxaTable = "data.frame"))

#' Workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr  ggpubr methods ComplexHeatmap
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @param se a object whose class is the RangedSummarizedExperiment generated from package SummarizedExperiment. Columns of variables indicate sample information in ColData.
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#' @param TaxaInfo a vector indicates the taxanomic information. The length of the vector is equal to the row number of count table.
#' @param FeatureInfo a vector indicates the genomic information. The length of the vector is equal to the row number of count table.
#' @param TaxaTable a data frame whose columns are kindom, phylum, class, order, family, genus and species indicate the taxanomic information. The row number of TaxaTable is equal to the row number of count table
#'
#' @return a WhamInput object contains all input information
#' @export
#'
#' @examples
#' # see the vignette for more details
WhamInputfromBiobakery <- function(se = NULL,design = NULL,TaxaInfo = NULL,FeatureInfo = NULL,TaxaTable = NULL){
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }
  object <- new("WhamInputfromBiobakery",se,design = design, TaxaInfo =TaxaInfo, FeatureInfo = FeatureInfo,TaxaTable = TaxaTable)
  return(object)
}


#' Workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr  ggpubr methods ComplexHeatmap
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @param countData a data frame or a matrix contains the counts generated from the biobakery process. The first three columns of the count table have to be "Acc","Taxa",and "Feature" respectively(case insensitive)
#' @param colData metadata table contains experiment design with at least one column. The row number of coldata has to be equal to the columns of count table excluding the first three columns.
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#'
#' @return a WhamInput object contains countData, ColData, design, Taxainfo, FeatureInfo, TaxaTable
#' @export
#'
#' @examples
#' # see the vignette for more details
WhamFromBiobakery <- function(countData = NULL, colData = NULL,design = NULL){
  #####check input format
  if( !(is.data.frame(countData) | is.matrix(countData)))
    stop("countData has to be a data frame or a matrix")

  if( !(is.data.frame(colData)))
    stop("colData has to be a data frame")

  if(length(grep("\\bacc\\b",ignore.case = T,names(countData[1]))) != 1)
    stop("the first column name has to be ACC (case insensitive)")

  if(length(grep("\\bfeature\\b",ignore.case = T,names(countData[2]))) !=1)
    stop("the second column name has to be feature (case insensitive)")

  if(length(grep("\\btaxa\\b",ignore.case = T,names(countData[3]))) !=1)
    stop("the third column name has to be taxa (case insensitive)")

  if(ncol(countData) != (nrow(colData)+3))
    stop("the column number of countData doesn't match the row number of colData")

  ####check the countdata value
  countdata_subset = countData[,c(4:ncol(countData))]
  countdata_subset = as.matrix(countdata_subset)

  if(any(is.na(countdata_subset )))
    stop("countData can't contain NA, NaN, NULL value")

  if(any(is.nan(countdata_subset )))
    stop("countData can't contain NA, NaN, NULL value")

  if(any(is.null(countdata_subset )))
    stop("countData can't contain NA, NaN, NULL value")

  if(any(!is.finite(countdata_subset)))
    stop("countData can't contain NA, NaN, NULL value")

  if(any(countdata_subset<0))
    stop("value in countData must be non-negative number")

  ####create se object
  se = SummarizedExperiment(assays = list(countdata_subset),colData = colData)

  #####create TaxaTable
  taxainfo = as.character(countData[,3])
  TaxaTable = TaxaToDF_fn(taxainfo)

  ######create wham object
  object <- WhamInputfromBiobakery(se,design = design,TaxaInfo = as.character(countData[,3]),FeatureInfo = as.character(countData[,2]),TaxaTable = TaxaTable)

  if(!all(all.vars(object@design) %in% names(colData)))
    stop("every variable in design must be contained in the column name of metadata")

  return(object)
}

#################

setClass("WhamResult",
         representation = representation(
           DE.result = "data.frame",
           contrast = "character"))


#' workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @param object required. Object generated from WhamFromBiobakery, WhamFromEBI or WhamFrom16s.
#' @param DE required. Type of differential expression. Options are "taxa",or "feature".
#' @param taxa.level required, c("k","p","o","c","f","g","s"). Collapsing the count to provided taxonomic level. options are "k","p","c","o","f","g","s".Only need to declaimed when argument DE is setup as "taxa". Every taxonimic level has to have "k__","p__","c__","o__","f__","g__",or "s__"
#' @param contrast Specifying what comparison to extract from the WhamInput. a vector indicates the variable, numerator,and denominator. e.g. contrast = c("group","eye","mouth"). The argument contrast, ref.contrast, and All.contrast can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param ref.contrast Specifying a reference element in the variable, All other elements in the variable will be processed one vs one comparison (vs reference variable). A vector indicates the variable, denominator. e.g. contrast = c("cytokine","Buffer"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param All.contrast Specifying a variable. Every two elements in the variable will be processed one vs one comparison. a vector indicates the variable, denominator. e.g. contrast = c("cytokine"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param WhamTransformation c("count","WhamTransformation"). Default value is "count". Setting as "count" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#' @param useMC use multicore by default(FALSE). Multi core processing will be attempted with the BiocParallel package, then the parallel package. If neither are installed, serial processing will be used.
#' @param aldex.module c("t.test","anova"), Default value is t.test(Welch's t and Wilcoxon rank test). "anova" will perform glm and Kruskal Wallace tests for one-way ANOVA.
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr ggpubr methods ComplexHeatmap
#' @importFrom stats aggregate
#' @importFrom utils combn
#' @importFrom stats model.matrix
#'
#' @return a dataframe contains statistic results or a WhamResult object contains statistic results and contrast.
#' @return we.ep - Expected P value of Welch's t test
#' @return we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
#' @return wi.ep - Expected P value of Wilcoxon rank test
#' @return wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
#' @return rab.all: a vector containing the median clr value for each feature
#' @return rab.win.conditionA: a vector containing the median clr value for each feature in condition A
#' @return rab.win.conditionB: a vector containing the median clr value for each feature in condition B
#' @return diff.btw: a vector containing the per-feature median difference between condition A and B
#' @return diff.win: a vector containing the per-feature maximum median difference between Dirichlet instances within conditions
#' @return effect:  median effect size: diff.btw / max(diff.win) for all instances
#' @return overlap: proportion of effect size that overlaps 0 (i.e. no effect)
#' @return kw.ep a vector containing the expected p-value of the Kruskal-Wallis test for each feature
#' @return kw.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @return glm.ep a vector containing the expected p-value of the glm ANOVA for each feature
#' @return glm.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @export
#'
#' @examples
#' # see the vignette for more details

WhamDEFromBiobakery <- function(object = NULL,
                                DE = c("feature","taxa"),
                                contrast,
                                ref.contrast,
                                All.contrast,
                                taxa.level =  c("k","p","o","c","f","g","s"),
                                WhamTransformation = c("count","WhamTransformation"),
                                useMC = F,
                                aldex.module = c("t.test","anova")
){
  #######check object class#####
  if((!missing(contrast))&(!missing(ref.contrast))& !missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(contrast))&(!missing(ref.contrast)))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(contrast))&!missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(ref.contrast))& !missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if(!is(object,"WhamInputfromBiobakery"))
    stop("object has to be the class of WhamInputfromBiobakery")

  if(missing(aldex.module)){
    message("aldex.module is not provided,default choice is 't.test'(Welch's t-test and Wilcoxon rank test)")
  }
  aldex.module = match.arg(aldex.module,c("t.test","anova"))

  ######wham Transformation######
  WhamTransformation = match.arg(WhamTransformation,choices = c("count","WhamTransformation"))

  if(WhamTransformation == "WhamTransformation"){
    count_trans = matrix(as.integer(t(t(assay(object))/apply(assay(object),2,sum)) * 1e6),nrow = nrow(assay(object)))
    colnames(count_trans) = colnames(assay(object))
    assay(object) = count_trans
    message("Wham transformation has been initiated\n")
  }

  ######feature or taxa table will be merged and processed by the following analysis.
  if(missing(DE))
    stop("DE is required, choose 'feature' or 'taxa'")

  DE <- match.arg(DE, choices=c("feature","taxa"))
  if(DE == "feature"){
    count_merged = aggregate(assay(object),by = list(object@FeatureInfo),sum)
    rownames(count_merged) = count_merged[,1]
    count_merged = count_merged[,-1]
  }

  if(DE == "taxa"){
    if(missing(taxa.level))
      stop("you have to provide taxa level which you want to collapse at")
    taxa.level=match.arg(taxa.level,choices = c("k","p","o","c","f","g","s"))

    Taxa_Info = object@TaxaInfo
    taxaLEVEL = TaxaIden_fn(TaxaInfo = Taxa_Info,taxa.level = taxa.level)

    taxalevel = apply(object@TaxaTable[,1:grep(taxaLEVEL,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel = gsub("unassigned;|;unassigned","",taxalevel)

    count_merged = aggregate(assay(object),by = list(taxalevel),sum)
    rownames(count_merged) = count_merged[,1]
    CollTaxaTable = TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) = count_merged[,1]
    count_merged = count_merged[,-1]
  }

  if(aldex.module == "t.test" & length(all.vars(object@design)) == 1){
    message("perform the Welch's t and Wilcoxon rank test\n")
    ##########if contrast is not claimed, then use the first vs the last in the var#########

    ####1.default setting####
    if(missing(contrast) & missing(ref.contrast) & missing(All.contrast)){
      message("contrast is not provided,use the first condition vs the last condition in the variable in design\n")
      designVars <- all.vars(object@design)
      lastVarName <- designVars[length(designVars)]
      lastVar <- colData(object)[[lastVarName]]
      lastVar = factor(lastVar)
      if (is.factor(lastVar)) {
        nlvls <- nlevels(lastVar)
        contrast <- c(lastVarName, levels(lastVar)[1], levels(lastVar)[nlvls])
      }
    }

    ##################t.test aldex function#######
    if(missing(ref.contrast) & missing(All.contrast)){
      DE_index = colData(object)[[contrast[1]]]
      numerato_index = which(DE_index==contrast[2])
      denominator_index = which(DE_index==contrast[3])
      count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
      group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                       paste0("BBBBB_",as.character(DE_index[denominator_index])))
      x <- aldex.clr(count_merged_subset,group_subset,useMC = useMC)
      x.tt <- aldex.ttest(x)
      x.effect <- aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE,useMC = useMC)
      colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
      if(DE == "taxa"){
        Level = data.frame(Level = rep(taxaLEVEL,nrow(x.tt)))
        rank_name = CollTaxaTable[rownames(x.tt),taxaLEVEL]
        result_table = cbind(rank_name,Level,x.tt,x.effect)}

      if(DE == "feature"){
        result_table = cbind(x.tt,x.effect)}
    }

    if(!missing(ref.contrast) & missing(All.contrast) & missing(contrast)){
      DE_index = factor(colData(object)[[ref.contrast[1]]])
      denominator_index = which(DE_index == ref.contrast[2])
      result_table = c()
      for(i in levels(DE_index)[levels(DE_index) != ref.contrast[2]]){
        numerato_index = which(DE_index==i)
        count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
        group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                         paste0("BBBBB_",as.character(DE_index[denominator_index])))
        x <- aldex.clr( count_merged_subset,  group_subset,useMC = useMC)
        x.tt <- aldex.ttest(x)
        x.effect <- aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE,useMC = useMC)
        colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
        result_table_unmerge = cbind(x.tt,x.effect)
        result_table_unmerge$COMPARISON = rep(paste0( i , " vs ", ref.contrast[2]),nrow(result_table_unmerge ))
        result_table_unmerge$name = row.names(result_table_unmerge)
        rownames(result_table_unmerge) = c()
        names(result_table_unmerge) = c("we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap","COMPARISON","name" )

        if(DE == "taxa"){
          result_table_unmerge$Level = rep(taxaLEVEL,nrow(result_table_unmerge))
          result_table_unmerge$rank_name = CollTaxaTable[result_table_unmerge$name,taxaLEVEL]
          result_table = rbind(result_table,result_table_unmerge)
        }

        if(DE == "feature"){
          result_table = rbind(result_table,result_table_unmerge)
        }
      }
      if(DE == "taxa"){result_table = result_table[,c("name","rank_name","Level","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]}
      if(DE == "feature"){result_table = result_table[,c("name","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]}
    }

    if(!missing(All.contrast) & missing(ref.contrast) & missing(contrast)){
      DE_index = factor(colData(object)[[All.contrast]])
      combination.matrix = combn(levels(DE_index),2)

      result_table = c()
      for( i in 1:ncol(combination.matrix)){
        contrast = c(All.contrast,combination.matrix[1,i],combination.matrix[2,i])
        numerato_index = which(DE_index==contrast[2])
        denominator_index =which(DE_index== contrast[3])
        count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
        group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                         paste0("BBBBB_",as.character(DE_index[denominator_index])))
        x <- aldex.clr( count_merged_subset,  group_subset,useMC = useMC)
        x.tt <- aldex.ttest(x)
        x.effect <-aldex.effect(x,include.sample.summary=FALSE, verbose=TRUE,useMC = useMC)
        colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
        result_table_unmerge = cbind(x.tt,x.effect)
        result_table_unmerge$COMPARISON = rep(paste0( contrast[2], " vs ", contrast[3]),nrow(result_table_unmerge))
        result_table_unmerge$name = row.names(result_table_unmerge)
        result_table_unmerge$group = rep(All.contrast,nrow(result_table_unmerge))
        rownames(result_table_unmerge) = c()
        names(result_table_unmerge) = c("we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap","COMPARISON","name","group" )

        if(DE == "taxa"){
          result_table_unmerge$Level = rep(taxaLEVEL,nrow(result_table_unmerge))
          result_table_unmerge$rank_name = CollTaxaTable[result_table_unmerge$name,taxaLEVEL]
          result_table = rbind(result_table,result_table_unmerge)
        }

        if(DE == "feature"){ result_table = rbind(result_table,result_table_unmerge)}
      }

      if(DE == "taxa"){result_table = result_table[,c("group","name","rank_name","Level","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]}
      if(DE == "feature"){result_table = result_table[,c("group","name","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]}
    }

    DE.result = result_table

    ######one vs one will return a object called whamresult
    if(missing(ref.contrast) & missing(All.contrast)){
      result = new("WhamResult",DE.result = DE.result, contrast = contrast)
    }

    ######one vs ref/all will return a merged data.frame
    if(!missing(ref.contrast)){
      result = DE.result
    }

    if(!missing(All.contrast)){
      result = DE.result
    }
  }

  ####anova like aldex function####
  if(aldex.module == "anova" & length(all.vars(object@design)) == 1 ){
    message("perform the glm and Kruskal Wallace tests for one-way ANOVA\n")
    designVars <- all.vars(object@design)
    anova.group <- designVars[length(designVars)]
    x <- aldex.clr(count_merged,
                   conds = as.character(as.data.frame(colData(object))[,anova.group]),
                   useMC = useMC)
    x.kw <- aldex.kw(x)
    result =x.kw

    if(DE == "taxa"){
      result$rank_name = CollTaxaTable[rownames(x.kw),taxaLEVEL]
      result$Level = rep(taxaLEVEL,nrow(result))
      result= result[,c("rank_name","Level","kw.ep","kw.eBH","glm.ep","glm.eBH")]
    }
  }

  #####multi-factor design two anova/ancova#####
  if(length(all.vars(object@design)) > 1 ){
    message("multi-factor complex design approaching, contrast shouldn't be provided")
    covariates <- as.data.frame(colData(object))[,all.vars(object@design)]
    mm <- model.matrix(object@design, covariates)

    x <- aldex.clr(count_merged,
                   mm,
                   useMC = useMC)
    glm.test <- aldex.glm(x, mm)
    result = glm.test

    if(DE == "taxa"){
      rank_name = CollTaxaTable[rownames(glm.test),taxaLEVEL]
      Level = rep(taxaLEVEL,nrow(result))
      result = cbind(rank_name,Level,result)
    }
  }

  return(result)
}

##########WhamBiobakery
setClass("WhamOutputBiobakery",
         contains = "RangedSummarizedExperiment",
         representation = representation(se = "RangedSummarizedExperiment",
                                         design = "ANY",
                                         TaxaInfo = "character",
                                         FeatureInfo = "character",
                                         DE.result = "ANY",
                                         TaxaTable = "data.frame" )
)
#' workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr  ggpubr methods ComplexHeatmap
#' @importFrom stats aggregate
#' @param countData a data frame or a matrix contains the counts generated from the biobakery process. The first three columns of the count table have to be "Acc","Taxa",and "Feature" respectively(case insensitive)
#' @param colData metadata table contains experiment design with at least one column. The row number of coldata has to be equal to the columns of count table excluding the first three columns.
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#' @param DE required. Type of differential expression. Options are "taxa",or "feature".
#' @param taxa.level required, c("k","p","o","c","f","g","s"). Collapsing the count to provided taxonomic level. options are "k","p","c","o","f","g","s".Only need to declaimed when argument DE is setup as "taxa". Every taxonimic level has to have "k__","p__","c__","o__","f__","g__",or "s__"
#' @param contrast Specifying what comparison to extract from the WhamInput. a vector indicates the variable, numerator,and denominator. e.g. contrast = c("group","eye","mouth"). The argument contrast, ref.contrast, and All.contrast can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param ref.contrast Specifying a reference element in the variable, All other elements in the variable will be processed one vs one comparison (vs reference variable). A vector indicates the variable, denominator. e.g. contrast = c("cytokine","Buffer"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param All.contrast Specifying a variable. Every two elements in the variable will be processed one vs one comparison. a vector indicates the variable, denominator. e.g. contrast = c("cytokine"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param WhamTransformation c("count","WhamTransformation"). Default value is "count". Setting as "count" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#' @param useMC use multicore by default(FALSE). Multi core processing will be attempted with the BiocParallel package, then the parallel package. If neither are installed, serial processing will be used.
#' @param aldex.module c("t.test","anova"), Default value is t.test(Welch's t and Wilcoxon rank test). "anova" will perform glm and Kruskal Wallace tests for one-way ANOVA
#' @return a dataframe contains statistic results or a WhamResult object contains statistic results and contrast.
#' @return we.ep - Expected P value of Welch's t test
#' @return we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
#' @return wi.ep - Expected P value of Wilcoxon rank test
#' @return wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
#' @return rab.all: a vector containing the median clr value for each feature
#' @return rab.win.conditionA: a vector containing the median clr value for each feature in condition A
#' @return rab.win.conditionB: a vector containing the median clr value for each feature in condition B
#' @return diff.btw: a vector containing the per-feature median difference between condition A and B
#' @return diff.win: a vector containing the per-feature maximum median difference between Dirichlet instances within conditions
#' @return effect:  median effect size: diff.btw / max(diff.win) for all instances
#' @return overlap: proportion of effect size that overlaps 0 (i.e. no effect)
#' @return kw.ep a vector containing the expected p-value of the Kruskal-Wallis test for each feature
#' @return kw.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @return glm.ep a vector containing the expected p-value of the glm ANOVA for each feature
#' @return glm.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @export
#'
#' @examples
#' dir.data <-  system.file("extdata","biobakery_sample_input.tsv",package = "Wham")
#' countData <-  read.delim(dir.data)
#' countData <- countData[1:2000,]
#' dir.meta = system.file("extdata","biobakery_sample_metadata.csv",package = "Wham")
#' metadata <- read.csv(dir.meta,row.names = 1)
#' wham_bbk <-  WhamBiobakery(countData = countData,
#'                            colData = metadata,
#'                            DE = "taxa",
#'                            design = ~ Location,
#'                            taxa.level = "s",
#'                            contrast = c("Location","Stool","Arm")
#'                            )
#'
#' # see the vignette for more details

WhamBiobakery <- function(countData =NULL,
                          colData = NULL,
                          design= NULL,
                          DE = c("taxa","feature"),
                          taxa.level =  c("k","p","o","c","f","g","s"),
                          contrast,
                          ref.contrast,
                          All.contrast,
                          aldex.module = c("t.test","anova"),
                          WhamTransformation = c("count","WhamTransformation"),
                          useMC = F
){
  if(DE == "taxa"){ if(missing(taxa.level))
    stop("you have to declaim taxa level which you want to collapse ")}

  output = WhamFromBiobakery(countData = countData,colData = colData,design = design)
  DE_result = WhamDEFromBiobakery(object = output,
                                  DE = DE,
                                  contrast=contrast,
                                  ref.contrast = ref.contrast,
                                  All.contrast = All.contrast,
                                  WhamTransformation = WhamTransformation,
                                  taxa.level = taxa.level,
                                  aldex.module = aldex.module,
                                  useMC = useMC
  )

  se = SummarizedExperiment(assays <- list(assay(output)),colData = colData)

  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }

  result <- new("WhamOutputBiobakery",se,
                design = output@design,
                TaxaInfo = output@TaxaInfo,
                FeatureInfo = output@FeatureInfo,
                TaxaTable = output@TaxaTable,
                DE.result = DE_result)

  return(result)
}




######################################################################
###################EBI################################################
######################################################################


setClass("WhamInputfromEBI",
         contains = "RangedSummarizedExperiment",
         representation = representation(
           se = "RangedSummarizedExperiment",
           design = "ANY",
           TaxaInfo = "character",
           FeatureInfo = "character",
           TaxaTable = "data.frame",
           featureCount = "matrix",
           taxaCount = "matrix"))

#' Workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr  ggpubr methods ComplexHeatmap
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @param se a object whose class is the RangedSummarizedExperiment generated from package SummarizedExperiment. Columns of variables indicate sample information in ColData.
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#' @param TaxaInfo a vector indicates the taxanomic information. The length of the vector is equal to the row number of count table.
#' @param FeatureInfo a vector indicates the genomic information. The length of the vector is equal to the row number of count table.
#' @param TaxaTable a data frame whose columns are kindom, phylum, class, order, family, genus and species indicate the taxanomic information. The row number of TaxaTable is equal to the row number of count table
#' @param featureCount a data frame or a matrix contains the feature counts generated from the EBI process. The first two columns of the count table have to be "Acc" ,and "Feature" respectively(case insensitive)
#' @param taxaCount a data frame or a matrix contains the taxa counts generated from the EBI process. The first column of the count table have to be "taxa" respectively(case insensitive)
#' @return a WhamInput object contains all input information
#' @export
#' @examples
#' # see the vignette for more details

WhamInputfromEBI <- function(se = NULL,design = NULL,TaxaInfo = NULL,FeatureInfo = NULL,TaxaTable = NULL,featureCount = NULL,taxaCount = NULL){
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }
  object <- new("WhamInputfromEBI",se,design = design, TaxaInfo =TaxaInfo, FeatureInfo = FeatureInfo,TaxaTable = TaxaTable,featureCount = featureCount,taxaCount = taxaCount )
  return(object)
}


#' Workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @param DE required. Type of differential expression. Options are "taxa",or "feature".
#' @param featureCount a data frame or a matrix contains the feature counts generated from the EBI process. The first two columns of the count table have to be "Acc" ,and "Feature" respectively(case insensitive)
#' @param taxaCount a data frame or a matrix contains the taxa counts generated from the EBI process. The first column of the count table have to be "taxa" respectively(case insensitive)
#' @param colData metadata table contains experiment design with at least one column. The row number of coldata has to be equal to the columns of count table excluding the first three columns.
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#' @return a WhamInput object contains taxaCount/featureCount, ColData, design, Taxainfo, FeatureInfo, TaxaTable
#' @export
#' @examples
#' # see the vignette for more details
WhamFromEBI <- function(DE = c("feature","taxa"),featureCount = NULL, taxaCount = NULL,colData = NULL,design = NULL){
  if(missing(DE))
    stop("DE is required. 'feature' or 'taxa'")

  #####check input format
  if( !(is.data.frame(featureCount) | is.matrix(featureCount)))
    stop("featureCount has to be a data frame or a matrix")

  if( !(is.data.frame(taxaCount) | is.matrix(taxaCount)))
    stop("taxaCount has to be a data frame or a matrix")

  if( !(is.data.frame(colData)))
    stop("colData has to be a data frame")

  if(length(grep("\\bacc\\b",ignore.case = T,names(featureCount[1]))) != 1)
    stop("the first column name in featureCount has to be ACC (case insensitive)")

  if(length(grep("\\bfeature\\b",ignore.case = T,names(featureCount[2]))) !=1)
    stop("the second column name in featureCounthas to be feature (case insensitive)")

  if(length(grep("\\btaxa\\b",ignore.case = T,names(taxaCount[1]))) !=1)
    stop("the first column name has in taxaFeature to be taxa (case insensitive)")

  if(ncol(taxaCount) != (nrow(colData)+1))
    stop("the column number of taxaData doesn't match the row number of colData")

  if(ncol(featureCount) != (nrow(colData)+2))
    stop("the column number of featureData doesn't match the row number of colData")

  ####check the featureCount value
  featureCount_subset = featureCount[,c(3:ncol(featureCount))]
  featureCount_subset = as.matrix(featureCount_subset)

  if(any(is.na(featureCount_subset )))
    stop("featureCount can't contain NA, NaN, NULL value")

  if(any(is.nan(featureCount_subset )))
    stop("featureCount can't contain NA, NaN, NULL value")

  if(any(is.null(featureCount_subset )))
    stop("featureCount can't contain NA, NaN, NULL value")

  if(any(!is.finite(featureCount_subset)))
    stop("featureCount can't contain NA, NaN, NULL value")

  if(any(featureCount_subset<0))
    stop("value in featureCount must be non-negative number")

  ####check the taxaCount value

  taxaCount_subset = taxaCount[,c(2:ncol(taxaCount))]
  taxaCount_subset = as.matrix(taxaCount_subset)

  if(any(is.na(taxaCount_subset )))
    stop("taxaCount can't contain NA, NaN, NULL value")

  if(any(is.nan(taxaCount_subset )))
    stop("taxaCount can't contain NA, NaN, NULL value")

  if(any(is.null(taxaCount_subset )))
    stop("taxaCount can't contain NA, NaN, NULL value")

  if(any(!is.finite(taxaCount_subset)))
    stop("taxaCount can't contain NA, NaN, NULL value")

  if(any(taxaCount_subset<0))
    stop("value in taxaCount must be non-negative integer")

  DE = match.arg(DE,choices = c("feature","taxa"))
  if(DE == "feature"){
    count_matrix = featureCount_subset
  }

  if(DE == "taxa"){
    count_matrix = taxaCount_subset
  }

  ####create se object
  se = SummarizedExperiment(assays = list(count_matrix),colData = colData)

  #####create TaxaTable
  taxainfo = as.character(taxaCount[,1])
  TaxaTable = TaxaToDF_fn(taxainfo)

  ####create wham object
  object <- WhamInputfromEBI(se,design = design,TaxaInfo = as.character(taxaCount[,1]),FeatureInfo =       as.character(featureCount[,2]),TaxaTable = TaxaTable,featureCount = featureCount_subset,taxaCount = taxaCount_subset)

  if(!all(all.vars(object@design) %in% names(colData)))
    stop("every variable in design must be contained in the column name of metadata")

  return(object)
}


setClass("WhamResult",
         representation = representation(
           DE.result = "data.frame",
           contrast = "character"))
#' workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr  ggpubr methods ComplexHeatmap
#' @importFrom utils combn
#' @importFrom stats aggregate
#' @importFrom stats model.matrix
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @param object required. Object generated from WhamFromBiobakery, WhamFromEBI or WhamFrom16s.
#' @param DE required. Type of differential expression. Options are "taxa",or "feature".
#' @param taxa.level required, c("k","p","o","c","f","g","s"). Collapsing the count to provided taxonomic level. options are "k","p","c","o","f","g","s".Only need to declaimed when argument DE is setup as "taxa". Every taxonimic level has to have "k__","p__","c__","o__","f__","g__",or "s__"
#' @param contrast Specifying what comparison to extract from the WhamInput. a vector indicates the variable, numerator,and denominator. e.g. contrast = c("group","eye","mouth"). The argument contrast, ref.contrast, and All.contrast can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param ref.contrast Specifying a reference element in the variable, All other elements in the variable will be processed one vs one comparison (vs reference variable). A vector indicates the variable, denominator. e.g. contrast = c("cytokine","Buffer"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param All.contrast Specifying a variable. Every two elements in the variable will be processed one vs one comparison. a vector indicates the variable, denominator. e.g. contrast = c("cytokine"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param WhamTransformation c("count","WhamTransformation"). Default value is "count". Setting as "count" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#' @param useMC use multicore by default(FALSE). Multi core processing will be attempted with the BiocParallel package, then the parallel package. If neither are installed, serial processing will be used.
#' @param aldex.module c("t.test","anova"), Default value is t.test(Welch's t and Wilcoxon rank test). "anova" will perform glm and Kruskal Wallace tests for one-way ANOVA
#' @return a dataframe contains statistic results or a WhamResult object contains statistic results and contrast.
#' @return we.ep - Expected P value of Welch's t test
#' @return we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
#' @return wi.ep - Expected P value of Wilcoxon rank test
#' @return wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
#' @return rab.all: a vector containing the median clr value for each feature
#' @return rab.win.conditionA: a vector containing the median clr value for each feature in condition A
#' @return rab.win.conditionB: a vector containing the median clr value for each feature in condition B
#' @return diff.btw: a vector containing the per-feature median difference between condition A and B
#' @return diff.win: a vector containing the per-feature maximum median difference between Dirichlet instances within conditions
#' @return effect:  median effect size: diff.btw / max(diff.win) for all instances
#' @return overlap: proportion of effect size that overlaps 0 (i.e. no effect)
#' @return kw.ep a vector containing the expected p-value of the Kruskal-Wallis test for each feature
#' @return kw.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @return glm.ep a vector containing the expected p-value of the glm ANOVA for each feature
#' @return glm.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @export
#'
#' @examples
#' # see the vignette for more details
WhamDEFromEBI <- function(object = NULL,
                          DE = c("feature","taxa"),
                          contrast,
                          ref.contrast,
                          All.contrast,
                          taxa.level =  c("k","p","o","c","f","g","s"),
                          WhamTransformation = c("count","WhamTransformation"),
                          useMC = F,
                          aldex.module = c("t.test","anova")
){

  if(!(is(object,"WhamInputfromEBI")))
    stop("object has to be the class of WhamInputfromEBI")

  DE <- match.arg(DE, choices=c("feature","taxa"))

  if(is(object,"WhamInputfromEBI")){
    count_mat = assay(object)
  }

  #######check object class#####
  if((!missing(contrast))&(!missing(ref.contrast))& !missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(contrast))&(!missing(ref.contrast)))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(contrast))&!missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(ref.contrast))& !missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if(missing(aldex.module)){
    message("aldex.module is not provided,default choice is 't.test'(Welch's t-test and Wilcoxon rank test)")
  }
  aldex.module = match.arg(aldex.module,c("t.test","anova"))

  ######wham Transformation######
  WhamTransformation = match.arg(WhamTransformation,choices = c("count","WhamTransformation"))

  if(WhamTransformation == "WhamTransformation"){
    count_trans = matrix(as.integer(t(t(count_mat)/apply(count_mat,2,sum)) * 1e6),nrow = nrow(count_mat))
    colnames(count_trans) = colnames(count_mat)
    count_mat = count_trans
    message("Wham transformation has been initiated\n")
  }

  ######feature or taxa table will be merged and processed by the following analysis.

  DE <- match.arg(DE, choices=c("feature","taxa"))
  if(DE == "feature"){
    count_merged = aggregate(count_mat,by = list(object@FeatureInfo),sum)
    rownames(count_merged) = count_merged[,1]
    count_merged = count_merged[,-1]
  }

  if(DE == "taxa"){
    if(missing(taxa.level))
      stop("you have to declaim taxa level which you want to collapse ")

    taxa.level=match.arg(taxa.level,choices = c("k","p","o","c","f","g","s"))

    #####check taxa
    Taxa_Info = object@TaxaInfo
    taxaLEVEL = TaxaIden_fn(TaxaInfo = Taxa_Info,taxa.level = taxa.level)

    ###collapse taxa
    taxalevel = apply(object@TaxaTable[,1:grep(taxaLEVEL,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel = gsub("unassigned;|;unassigned","",taxalevel)

    count_merged = aggregate(assay(object),by = list(taxalevel),sum)
    rownames(count_merged) = count_merged[,1]
    CollTaxaTable = TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) = count_merged[,1]
    count_merged = count_merged[,-1]
  }



  ###################t.test
  if(aldex.module == "t.test" & length(all.vars(object@design)) == 1){
    message("perform the Welch's t and Wilcoxon rank test\n")
    ##########if contrast is not claimed, then use the first vs the last in the var#########
    if(missing(contrast) & missing(ref.contrast) & missing(All.contrast)){
      message("contrast is not provided,use the first condition vs the last condition in the variable in design\n")

      designVars <- all.vars(object@design)
      lastVarName <- designVars[length(designVars)]
      lastVar <- colData(object)[[lastVarName]]
      lastVar = factor(lastVar)
      if (is.factor(lastVar)) {
        nlvls <- nlevels(lastVar)
        contrast <- c(lastVarName, levels(lastVar)[1], levels(lastVar)[nlvls])
      }
    }

    ##################aldex function#######
    if(missing(ref.contrast) & missing(All.contrast)){
      DE_index = colData(object)[[contrast[1]]]
      numerato_index = which(DE_index==contrast[2])
      denominator_index = which(DE_index==contrast[3])

      count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
      group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                       paste0("BBBBB_",as.character(DE_index[denominator_index])))
      x <- aldex.clr( count_merged_subset,  group_subset,useMC = useMC)
      x.tt <- aldex.ttest(x)
      x.effect <-aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE,useMC = useMC)
      colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
      if(DE == "taxa"){
        Level = data.frame(Level = rep(taxaLEVEL,nrow(x.tt)))
        rank_name =  CollTaxaTable[rownames(x.tt),taxaLEVEL]
        result_table = cbind(rank_name,Level,x.tt,x.effect)}

      if(DE == "feature"){
        result_table = cbind(x.tt,x.effect)}
    }

    if(!missing(ref.contrast) & missing(All.contrast)  & missing(contrast)){
      DE_index = factor(colData(object)[[ref.contrast[1]]])
      denominator_index =which(DE_index== ref.contrast[2])
      result_table = c()
      for(i in levels(DE_index)[levels(DE_index) != ref.contrast[2]]){
        numerato_index = which(DE_index==i)
        count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
        group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                         paste0("BBBBB_",as.character(DE_index[denominator_index])))
        x <- aldex.clr( count_merged_subset,  group_subset,useMC = useMC)
        x.tt <- aldex.ttest(x)
        x.effect <- aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE,useMC = useMC)
        colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
        result_table_unmerge = cbind(x.tt,x.effect)
        result_table_unmerge$COMPARISON = rep(paste0( i , " vs ", ref.contrast[2]),nrow(result_table_unmerge ))
        result_table_unmerge$name = row.names(result_table_unmerge)
        rownames(result_table_unmerge) = c()
        names(result_table_unmerge) = c("we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap","COMPARISON","name" )

        if(DE == "taxa"){
          result_table_unmerge$Level = rep(taxaLEVEL,nrow(result_table_unmerge))
          result_table_unmerge$rank_name = CollTaxaTable[result_table_unmerge$name,taxaLEVEL]
          result_table = rbind(result_table,result_table_unmerge)}

        if(DE == "feature"){ result_table = rbind(result_table,result_table_unmerge)}
      }
      if(DE == "taxa"){result_table = result_table[,c("name","rank_name","Level","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]}
      if(DE == "feature"){result_table = result_table[,c("name","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]}
    }

    if(!missing(All.contrast) & missing(ref.contrast) & missing(contrast)){
      DE_index = factor(colData(object)[[All.contrast]])
      combination.matrix = combn(levels(DE_index),2)
      result_table = c()
      for( i in 1:ncol(combination.matrix)){
        contrast = c(All.contrast,combination.matrix[1,i],combination.matrix[2,i])
        numerato_index = which(DE_index==contrast[2])
        denominator_index =which(DE_index== contrast[3])
        count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
        group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                         paste0("BBBBB_",as.character(DE_index[denominator_index])))
        x <- aldex.clr( count_merged_subset,  group_subset,useMC = useMC)
        x.tt <- aldex.ttest(x)
        x.effect <-aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE,useMC = useMC)
        colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
        result_table_unmerge = cbind(x.tt,x.effect)
        result_table_unmerge$COMPARISON = rep(paste0( contrast[2], " vs ", contrast[3]),nrow(result_table_unmerge))
        result_table_unmerge$name = row.names(result_table_unmerge)
        result_table_unmerge$group = rep(All.contrast,nrow(result_table_unmerge))
        rownames(result_table_unmerge) = c()
        names(result_table_unmerge) = c("we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap","COMPARISON","name","group" )

        if(DE == "taxa"){
          result_table_unmerge$Level = rep(taxaLEVEL,nrow(result_table_unmerge))
          result_table_unmerge$rank_name = CollTaxaTable[result_table_unmerge$name,taxaLEVEL]
          result_table = rbind(result_table,result_table_unmerge)}

        if(DE == "feature"){ result_table = rbind(result_table,result_table_unmerge)}
      }

      if(DE == "taxa"){result_table = result_table[,c("group","name","rank_name","Level","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]}
      if(DE == "feature"){result_table = result_table[,c("group","name","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]}
    }

    DE.result = result_table

    ######one vs one will return a object called whamresult
    if(missing(ref.contrast) & missing(All.contrast)){
      result = new("WhamResult",DE.result = DE.result, contrast = contrast)
    }

    ######one vs ref will return a merged data.frame
    if(!missing(ref.contrast)){
      result = DE.result
    }

    if(!missing(All.contrast)){
      result = DE.result
    }
  }

  if(aldex.module == "anova" & length(all.vars(object@design)) == 1 ){
    message("perform the glm and Kruskal Wallace tests for one-way ANOVA\n")
      designVars <- all.vars(object@design)
      anova.group <- designVars[length(designVars)]
    x <- aldex.clr(count_merged,
                   conds = as.character(as.data.frame(colData(object))[,anova.group]),
                   useMC = useMC)
    x.kw <- aldex.kw(x)

    result <- x.kw

    if(DE == "taxa"){
      result$rank_name <- CollTaxaTable[rownames(x.kw),taxaLEVEL]
      result$Level <- rep(taxaLEVEL,nrow(result))
      result <- result[,c("rank_name","Level","kw.ep","kw.eBH","glm.ep","glm.eBH")]
    }
  }

  if(length(all.vars(object@design)) > 1 ){
    message("multi-factor complex design approaching, contrast shouldn't be provided")
    covariates <- as.data.frame(colData(object))[,all.vars(object@design)]
    mm <- model.matrix(object@design, covariates)

    x <- aldex.clr(count_merged,
                   mm,
                   useMC = useMC)
    glm.test <- aldex.glm(x, mm)
    result <- glm.test

    if(DE == "taxa"){
      rank_name <-  CollTaxaTable[rownames(glm.test),taxaLEVEL]
      Level <-  rep(taxaLEVEL,nrow(result))
      result <-  cbind(rank_name,Level,result)
    }
  }

  return(result)
}














setClass("WhamOutputEBI",
         contains = "RangedSummarizedExperiment",
         representation = representation(se = "RangedSummarizedExperiment",
                                         design = "ANY",
                                         TaxaInfo = "character",
                                         FeatureInfo = "character",
                                         DE.result = "ANY",
                                         TaxaTable = "data.frame",
                                         featureCount = "matrix",
                                         taxaCount = "matrix"))


#' workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @param featureCount a data frame or a matrix contains the feature counts generated from the EBI process. The first two columns of the count table have to be "Acc" ,and "Feature" respectively(case insensitive)
#' @param taxaCount a data frame or a matrix contains the taxa counts generated from the EBI process. The first column of the count table have to be "taxa" respectively(case insensitive)
#' @param colData metadata table contains experiment design with at least one column. The row number of coldata has to be equal to the columns of count table excluding the first three columns.
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#' @param DE required. Type of differential expression. Options are "taxa",or "feature".
#' @param taxa.level required, c("k","p","o","c","f","g","s"). Collapsing the count to provided taxonomic level. options are "k","p","c","o","f","g","s".Only need to declaimed when argument DE is setup as "taxa". Every taxonimic level has to have "k__","p__","c__","o__","f__","g__",or "s__"
#' @param contrast Specifying what comparison to extract from the WhamInput. a vector indicates the variable, numerator,and denominator. e.g. contrast = c("group","eye","mouth"). The argument contrast, ref.contrast, and All.contrast can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param ref.contrast Specifying a reference element in the variable, All other elements in the variable will be processed one vs one comparison (vs reference variable). A vector indicates the variable, denominator. e.g. contrast = c("cytokine","Buffer"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param All.contrast Specifying a variable. Every two elements in the variable will be processed one vs one comparison. a vector indicates the variable, denominator. e.g. contrast = c("cytokine"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param WhamTransformation c("count","WhamTransformation"). Default value is "count". Setting as "count" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#' @param useMC use multicore by default(FALSE). Multi core processing will be attempted with the BiocParallel package, then the parallel package. If neither are installed, serial processing will be used.
#' @param aldex.module c("t.test","anova"), Default value is t.test(Welch's t and Wilcoxon rank test). "anova" will perform glm and Kruskal Wallace tests for one-way ANOVA
#' @return a dataframe contains statistic results or a WhamResult object contains statistic results and contrast.
#' @return we.ep - Expected P value of Welch's t test
#' @return we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
#' @return wi.ep - Expected P value of Wilcoxon rank test
#' @return wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
#' @return rab.all: a vector containing the median clr value for each feature
#' @return rab.win.conditionA: a vector containing the median clr value for each feature in condition A
#' @return rab.win.conditionB: a vector containing the median clr value for each feature in condition B
#' @return diff.btw: a vector containing the per-feature median difference between condition A and B
#' @return diff.win: a vector containing the per-feature maximum median difference between Dirichlet instances within conditions
#' @return effect:  median effect size: diff.btw / max(diff.win) for all instances
#' @return overlap: proportion of effect size that overlaps 0 (i.e. no effect)
#' @return kw.ep a vector containing the expected p-value of the Kruskal-Wallis test for each feature
#' @return kw.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @return glm.ep a vector containing the expected p-value of the glm ANOVA for each feature
#' @return glm.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#'
#' @examples
#' dir.feature <-  system.file("extdata","EBI_feature_input.tsv",package = "Wham")
#' feature_EBI <-  read.delim(dir.feature)[,1:32]
#' dir.taxa <-  system.file("extdata","EBI_taxa_input.tsv",package = "Wham")
#' taxa_EBI <-  read.delim(dir.taxa)[,1:31]
#'
#' ##create fake metadata information
#' metadata = data.frame(group = rep(c("arm","mouth","eye"),c(10,10,10)),
#'                       age = sample(15:50,30,
#'                       replace = TRUE))
#' rownames(metadata) <- colnames(feature_EBI)[3:32]
#'
#' WhamEBI(featureCount = feature_EBI,
#'        taxaCount = taxa_EBI,
#'        colData = metadata,
#'        design = ~ group,
#'        DE ="taxa",
#'        taxa.level = "g")
#'
#' # see the vignette for more details
#'
#' @export
WhamEBI <- function(featureCount = NULL,
                    taxaCount = NULL,
                    colData = NULL,
                    design = NULL,
                    DE = c("taxa","feature"),
                    taxa.level =  c("k","p","o","c","f","g","s"),
                    contrast,
                    ref.contrast,
                    All.contrast,
                    WhamTransformation = c("count","WhamTransformation"),
                    useMC = F,
                    aldex.module = c("t.test","anova")
                    ){
  if(missing(DE)){
    stop("DE is required,'taxa' or 'feature'")
  }
  DE <- match.arg(DE,choices = c("feature","taxa"))
  output <- WhamFromEBI(DE = DE,featureCount = featureCount,taxaCount = taxaCount,colData = colData,design = design)
  DE_result <- WhamDEFromEBI(object = output,
                            DE = DE,
                            contrast=contrast,
                            ref.contrast = ref.contrast,
                            All.contrast = All.contrast,
                            WhamTransformation = WhamTransformation,
                            taxa.level = taxa.level,
                            useMC = useMC,
                            aldex.module = aldex.module

  )

  se <-  SummarizedExperiment(assays <- list(assay(output)),colData = colData)

  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }

  result <- new("WhamOutputEBI",se,
                design = output@design,
                TaxaInfo = output@TaxaInfo,
                FeatureInfo = output@FeatureInfo,
                TaxaTable = output@TaxaTable,
                DE.result = DE_result,
                featureCount = output@featureCount,
                taxaCount = output@taxaCount
  )

  return(result)
}


##################################################################
###############################16s################################
##################################################################

setClass("WhamInputfrom16s",
         contains = "RangedSummarizedExperiment",
         representation = representation(
           se = "RangedSummarizedExperiment",
           design = "ANY",
           TaxaInfo = "character",
           TaxaTable = "data.frame"))



#' Workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr ggpubr methods ComplexHeatmap
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @param se a object whose class is the RangedSummarizedExperiment generated from package SummarizedExperiment. Columns of variables indicate sample information in ColData.
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#' @param TaxaInfo a vector indicates the taxanomic information. The length of the vector is equal to the row number of count table.
#' @param TaxaTable a data frame whose columns are kindom, phylum, class, order, family, genus and species indicate the taxonomic information. The row number of TaxaTable is equal to the row number of count table
#'
#' @return a WhamInput object contains all input information
#' @export
#'
#' @examples
#' # see the vignette for more details
WhamInputfrom16s <- function(se = NULL,design = NULL,TaxaInfo = NULL,TaxaTable = NULL){
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }
  object <- new("WhamInputfrom16s",se,design = design, TaxaInfo =TaxaInfo,TaxaTable = TaxaTable)
  return(object)
}

#' Workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr ggpubr methods ComplexHeatmap
#' @importClassesFrom phyloseq phyloseq
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @param taxaCount a data frame or a matrix contains the counts generated from the 16s process. The first column of the count table have to be "Taxa"(case insensitive)
#' @param colData metadata table contains experiment design with at least one column. The row number of coldata has to be equal to the columns of count table excluding the first three columns.
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#' @param phyloseq_object object from package phyloseq
#' @return a WhamInput object contains countData, ColData, design, Taxainfo, FeatureInfo, TaxaTable
#' @export
#' @examples
#' # see the vignette for more details

WhamFrom16s <- function(taxaCount = NULL,colData = NULL,design = NULL,phyloseq_object = NULL ){
  if(!is.null(phyloseq_object)){
    if(!is(phyloseq_object,"phyloseq"))
      stop("'phyloseq_object' must be a phyloseq object")
  }

  if(!is.null(phyloseq_object)& c(!is.null(taxaCount)|!is.null(colData)))
    stop("phyloseq object has been provided, taxaCount or colData is not needed")

  if(!is.null(phyloseq_object)&is.null(taxaCount)&is.null(colData)){
    message("phyloseq object has been provided" )
    if(is(phyloseq_object,"phyloseq")){
      colData <- data.frame(phyloseq_object@sam_data)
      taxaCount <- phyloseq_object@otu_table@.Data
      taxatable <- phyloseq_object@tax_table@.Data
      taxatable[taxatable == ""] = NA
      taxainfo <- apply(taxatable,1,paste,collapse = ";")
      taxainfo <- gsub(";NA|NA;| ","",taxainfo)
      taxainfo <- gsub("d__","k__",taxainfo)
      taxaCount <- as.data.frame(taxaCount)
      taxaCount$taxa <- taxainfo
      taxaCount <- taxaCount[,c(ncol(taxaCount),1: c(ncol(taxaCount) - 1))]
    }
  }

  #####check input format
  if( !(is.data.frame(taxaCount) | is.matrix(taxaCount)))
    stop("taxaCount has to be a data frame or a matrix")

  if( !(is.data.frame(colData)))
    stop("colData has to be a data frame")

  if(length(grep("\\btaxa\\b",ignore.case = T,names(taxaCount[1]))) !=1)
    stop("the first column name has in taxaFeature to be taxa (case insensitive)")

  if(ncol(taxaCount) != (nrow(colData)+1))
    stop("the column number of taxaData doesn't match the row number of colData")

  ####check the taxaCount value

  taxaCount_subset = taxaCount[,c(2:ncol( taxaCount))]
  taxaCount_subset = as.matrix(taxaCount_subset)

  if(any(is.na(taxaCount_subset )))
    stop("taxaCount can't contain NA, NaN, NULL value")

  if(any(is.nan(taxaCount_subset )))
    stop("taxaCount can't contain NA, NaN, NULL value")

  if(any(is.null(taxaCount_subset )))
    stop("taxaCount can't contain NA, NaN, NULL value")

  if(any(!is.finite(taxaCount_subset)))
    stop("taxaCount can't contain NA, NaN, NULL value")

  if(any(taxaCount_subset<0))
    stop("value in taxaCount must be non-negative number")

  ####create se object
  se = SummarizedExperiment(assays = list(taxaCount_subset),colData = colData)

  #####create TaxaTable

  taxainfo = as.character(taxaCount[,1])
  TaxaTable = TaxaToDF_fn(taxainfo)

  #####create wham object
  object <- WhamInputfrom16s(se,design = design,TaxaInfo = as.character(taxaCount[,1]),TaxaTable = TaxaTable)

  if(!all(all.vars(object@design) %in% names(colData)))
    stop("every variable in design must be contained in the column name of metadata")

  return(object)
}


setClass("WhamResult",
         representation = representation(
           DE.result = "data.frame",
           contrast = "character"))


#' workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @param object required. Object generated from WhamFromBiobakery, WhamFromEBI or WhamFrom16s.
#' @param taxa.level required, c("k","p","o","c","f","g","s"). Collapsing the count to provided taxonomic level. options are "k","p","c","o","f","g","s".Only need to declaimed when argument DE is setup as "taxa". Every taxonimic level has to have "k__","p__","c__","o__","f__","g__",or "s__"
#' @param contrast Specifying what comparison to extract from the WhamInput. a vector indicates the variable, numerator,and denominator. e.g. contrast = c("group","eye","mouth"). The argument contrast, ref.contrast, and All.contrast can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param ref.contrast Specifying a reference element in the variable, All other elements in the variable will be processed one vs one comparison (vs reference variable). A vector indicates the variable, denominator. e.g. contrast = c("cytokine","Buffer"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param All.contrast Specifying a variable. Every two elements in the variable will be processed one vs one comparison. a vector indicates the variable, denominator. e.g. contrast = c("cytokine"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param WhamTransformation c("count","WhamTransformation"). Default value is "count". Setting as "count" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#' @param useMC use multicore by default(FALSE). Multi core processing will be attempted with the BiocParallel package, then the parallel package. If neither are installed, serial processing will be used.
#' @param aldex.module c("t.test","anova"), Default value is t.test(Welch's t and Wilcoxon rank test). "anova" will perform glm and Kruskal Wallace tests for one-way ANOVA
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr  ggpubr methods ComplexHeatmap
#' @importFrom stats aggregate
#' @importFrom utils combn
#' @importFrom stats model.matrix
#' @importClassesFrom phyloseq phyloseq
#'
#' @return a dataframe contains statistic results or a WhamResult object contains statistic results and contrast.
#' @return we.ep - Expected P value of Welch's t test
#' @return we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
#' @return wi.ep - Expected P value of Wilcoxon rank test
#' @return wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
#' @return rab.all: a vector containing the median clr value for each feature
#' @return rab.win.conditionA: a vector containing the median clr value for each feature in condition A
#' @return rab.win.conditionB: a vector containing the median clr value for each feature in condition B
#' @return diff.btw: a vector containing the per-feature median difference between condition A and B
#' @return diff.win: a vector containing the per-feature maximum median difference between Dirichlet instances within conditions
#' @return effect:  median effect size: diff.btw / max(diff.win) for all instances
#' @return overlap: proportion of effect size that overlaps 0 (i.e. no effect)
#' @return kw.ep a vector containing the expected p-value of the Kruskal-Wallis test for each feature
#' @return kw.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @return glm.ep a vector containing the expected p-value of the glm ANOVA for each feature
#' @return glm.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @export
#'
#' @examples
#' # see the vignette for more details
WhamDEFrom16s <- function(object = NULL,
                          contrast,
                          ref.contrast,
                          All.contrast,
                          taxa.level =  c("k","p","o","c","f","g","s","otu"),
                          WhamTransformation = c("count","WhamTransformation"),
                          useMC = F,
                          aldex.module = c("t.test","anova")
){
  if(!(is(object,"WhamInputfrom16s")))
    stop("object has to be the class of WhamInputfrom16s")

  if(is(object,"WhamInputfrom16s")){
    count_mat = assay(object)
  }

  #######check object class#####
  if((!missing(contrast))&(!missing(ref.contrast))& !missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(contrast))&(!missing(ref.contrast)))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(contrast))&!missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")

  if((!missing(ref.contrast))& !missing(All.contrast))
    stop("contrast,All.contrast,or ref.contrast can't be declaimed at the same time.")


  if(missing(aldex.module)){
    message("aldex.module is not provided,default choice is 't.test'(Welch's t-test and Wilcoxon rank test)")
  }
  aldex.module = match.arg(aldex.module,c("t.test","anova"))

  ######wham Transformation######
  WhamTransformation = match.arg(WhamTransformation,choices = c("count","WhamTransformation"))

  if(WhamTransformation == "WhamTransformation"){
    count_trans = matrix(as.integer(t(t(count_mat)/apply(count_mat,2,sum)) * 1e6),nrow = nrow(count_mat))
    colnames(count_trans) = colnames(count_mat)
    count_mat = count_trans
    message("Wham transformation has been initiated\n")
  }

  ######feature or taxa table will be merged and processed by the following analysis.
  if(missing(taxa.level))
    stop("you have to declaim taxa level which you want to collapse ")

  if(taxa.level != "otu"){
    taxa.level <- match.arg(taxa.level,choices = c("k","p","o","c","f","g","s","otu"))
    Taxa_Info <- object@TaxaInfo
    taxaLEVEL <- TaxaIden_fn(TaxaInfo = Taxa_Info,taxa.level = taxa.level)
    taxalevel <- apply(object@TaxaTable[,1:grep(taxaLEVEL,names(object@TaxaTable))],1,paste,collapse =";")
    taxalevel <- gsub("unassigned;|;unassigned","",taxalevel)
    count_merged <- aggregate(assay(object),by = list(taxalevel),sum)
    rownames(count_merged) <- count_merged[,1]
    CollTaxaTable <- TaxaToDF_fn(count_merged[,1])
    rownames(CollTaxaTable) <- count_merged[,1]
    count_merged <- count_merged[,-1]
  }

  if(taxa.level == "otu"){
    taxaLEVEL <- "OTU"
    TaxaTable <- object@TaxaTable
    TaxaTable$OTU = paste0("OTU_",1:nrow(TaxaTable))
    CollTaxaTable = TaxaTable
    taxalevel = apply(TaxaTable[,1:grep(taxaLEVEL,names(TaxaTable))],1,paste,collapse =";")
    taxalevel = gsub("unassigned;|;unassigned","",taxalevel)
    count_merged = count_mat
    rownames(count_merged) = taxalevel
    rownames(CollTaxaTable) = taxalevel
  }

  #########t.test#####
  if(aldex.module == "t.test" & length(all.vars(object@design)) == 1){
    message("perform the Welch's t and Wilcoxon rank test\n")
    ##########if contrast is not claimed, then use the first vs the last in the var#########
    if(missing(contrast) & missing(ref.contrast) & missing(All.contrast)){
      message("contrast is not provided, use the first condition vs the last condition in the variable in design\n")
      designVars <- all.vars(object@design)
      lastVarName <- designVars[length(designVars)]
      lastVar <- colData(object)[[lastVarName]]
      lastVar = factor(lastVar)
      if (is.factor(lastVar)) {
        nlvls <- nlevels(lastVar)
        contrast <- c(lastVarName, levels(lastVar)[1], levels(lastVar)[nlvls])
      }
    }
    ##################aldex function#######
    if( missing(ref.contrast) & missing(All.contrast) ){
      DE_index = colData(object)[[contrast[1]]]
      numerato_index = which(DE_index == contrast[2])
      denominator_index = which(DE_index == contrast[3])
      count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
      group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                       paste0("BBBBB_",as.character(DE_index[denominator_index])))

      x <- aldex.clr( count_merged_subset, group_subset, useMC = useMC)
      x.tt <- aldex.ttest(x)
      x.effect <-aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE, useMC = useMC)
      colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
      Level = data.frame(Level = rep(taxaLEVEL,nrow(x.tt)))
      rank_name = CollTaxaTable[rownames(x.tt),taxaLEVEL]
      result_table = cbind(rank_name,Level,x.tt,x.effect)
    }

    if(!missing(ref.contrast) & missing(All.contrast)  & missing(contrast)){
      DE_index = factor(colData(object)[[ref.contrast[1]]])
      denominator_index =which(DE_index== ref.contrast[2])
      result_table = c()
      for(i in levels(DE_index)[levels(DE_index) != ref.contrast[2]]){
        numerato_index = which(DE_index==i)
        count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
        group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                         paste0("BBBBB_",as.character(DE_index[denominator_index])))
        x <- aldex.clr( count_merged_subset, group_subset, useMC = useMC)
        x.tt <- aldex.ttest(x)
        x.effect <- aldex.effect(x,include.sample.summary=FALSE, verbose=TRUE,useMC = useMC)
        colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
        result_table_unmerge = cbind(x.tt,x.effect)
        result_table_unmerge$COMPARISON = rep(paste0( i , " vs ", ref.contrast[2]),nrow(result_table_unmerge ))
        result_table_unmerge$name = row.names(result_table_unmerge)
        rownames(result_table_unmerge) = c()
        names(result_table_unmerge) = c("we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap","COMPARISON","name" )

        result_table_unmerge$Level = rep(taxaLEVEL,nrow(result_table_unmerge))
        result_table_unmerge$rank_name = CollTaxaTable[result_table_unmerge$name,taxaLEVEL]
        result_table = rbind(result_table,result_table_unmerge)
      }
      result_table = result_table[,c("name","rank_name","Level","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]
    }

    if(!missing(All.contrast) & missing(ref.contrast) & missing(contrast)){
      DE_index = factor(colData(object)[[All.contrast]])
      combination.matrix = combn(levels(DE_index),2)
      result_table = c()
      for( i in 1:ncol(combination.matrix)){
        contrast = c(All.contrast,combination.matrix[1,i],combination.matrix[2,i])
        numerato_index = which(DE_index==contrast[2])
        denominator_index =which(DE_index== contrast[3])
        count_merged_subset = count_merged[,c(numerato_index,denominator_index)]
        group_subset = c(paste0("AAAAA_",as.character(DE_index[numerato_index])),
                         paste0("BBBBB_",as.character(DE_index[denominator_index])))
        x <- aldex.clr( count_merged_subset, group_subset, useMC = useMC)
        x.tt <- aldex.ttest(x)
        x.effect <-aldex.effect(x, include.sample.summary=FALSE, verbose=TRUE,useMC = useMC)
        colnames(x.effect) = gsub("AAAAA_|BBBBB_","",colnames(x.effect ))
        result_table_unmerge = cbind(x.tt,x.effect)
        result_table_unmerge$COMPARISON = rep(paste0( contrast[2], " vs ", contrast[3]),nrow(result_table_unmerge))
        result_table_unmerge$name = row.names(result_table_unmerge)
        result_table_unmerge$group = rep(All.contrast,nrow(result_table_unmerge))
        rownames(result_table_unmerge) = c()
        names(result_table_unmerge) = c("we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap","COMPARISON","name","group" )
        result_table_unmerge$Level = rep(taxaLEVEL,nrow(result_table_unmerge))
        result_table_unmerge$rank_name = CollTaxaTable[result_table_unmerge$name,taxaLEVEL]
        result_table = rbind(result_table,result_table_unmerge)
        result_table = result_table[,c("group","rank_name","name","Level","COMPARISON","we.ep","we.eBH","wi.ep","wi.eBH","rab.all","rab.win.numerator","rab.win.denominator","diff.btw","diff.win","effect","overlap" )]
      }
    }

    DE.result = result_table
    ######one vs one will return a object called whamresult
    if(missing(ref.contrast) & missing(All.contrast)){
      result = new("WhamResult",DE.result = DE.result, contrast = contrast)
    }

    ######one vs ref will return a merged data.frame
    if(!missing(ref.contrast)){
      result = DE.result
    }

    if(!missing(All.contrast)){
      result = DE.result
    }
  }
  if(aldex.module == "anova" & length(all.vars(object@design)) == 1 ){
    message("perform the glm and Kruskal Wallace tests for one-way ANOVA\n")
      designVars <- all.vars(object@design)
      anova.group <- designVars[length(designVars)]
    x <- aldex.clr(count_merged,
                   conds = as.character(as.data.frame(colData(object))[,anova.group]),
                   useMC = useMC)
    x.kw <- aldex.kw(x)
    result =x.kw
    result$rank_name = CollTaxaTable[rownames(x.kw),taxaLEVEL]
    result$Level = rep(taxaLEVEL,nrow(result))
    result= result[,c("rank_name","Level","kw.ep","kw.eBH","glm.ep","glm.eBH")]
  }

  if(length(all.vars(object@design)) > 1 ){
    message("multi-factor complex design approaching, contrast shouldn't be provided")
    covariates <- as.data.frame(colData(object))[,all.vars(object@design)]
    mm <- model.matrix(object@design, covariates)

    x <- aldex.clr(count_merged,
                   mm,
                   useMC = useMC)
    glm.test <- aldex.glm(x, mm)
    result = glm.test

    if(DE == "taxa"){
      rank_name = CollTaxaTable[rownames(glm.test),taxaLEVEL]
      Level = rep(taxaLEVEL,nrow(result))
      result = cbind(rank_name,Level,result)
    }
  }

  return(result)
}


setClass("WhamOutput16s",
         contains = "RangedSummarizedExperiment",
         representation = representation(se = "RangedSummarizedExperiment",
                                         design = "ANY",
                                         TaxaInfo = "character",
                                         DE.result = "ANY",
                                         TaxaTable = "data.frame" )
)

#' workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @description R package allowing for exploratory analysis of metagenomics and metatranscriptomic data.  Includes visualization and statistical analysis on the gene family and taxa level. This package is developped based on package ALDEx2.
#' @import ALDEx2 SummarizedExperiment ggplot2 tidyr ggpubr methods ComplexHeatmap
#' @param taxaCount a data frame or a matrix contains the counts generated from the 16s process. The first column of the count table have to be "Taxa" respectively(case insensitive)
#' @param colData metadata table contains experiment design with at least one column. The row number of coldata has to be equal to the columns of count table excluding the first three columns.
#' @param phyloseq_object object from package phyloseq
#' @param design a formula expresses how the counts for each genes depend on the variables.e.g ~ group
#' @param taxa.level required, c("k","p","o","c","f","g","s","otu"). Collapsing the count to provided taxonomic level. options are "k","p","c","o","f","g","s",and "otu. Every taxonimic level has to have "k__","p__","c__","o__","f__","g__",or "s__"
#' @param contrast Specifying what comparison to extract from the WhamInput. a vector indicates the variable, numerator,and denominator. e.g. contrast = c("group","eye","mouth"). The argument contrast, ref.contrast, and All.contrast can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param ref.contrast Specifying a reference element in the variable, All other elements in the variable will be processed one vs one comparison (vs reference variable). A vector indicates the variable, denominator. e.g. contrast = c("cytokine","Buffer"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param All.contrast Specifying a variable. Every two elements in the variable will be processed one vs one comparison. a vector indicates the variable, denominator. e.g. contrast = c("cytokine"). The argument 'contrast', 'ref.contrast', and 'All.contrast' can't be declaimed at the same time. If none of three are declaimed, function will use the first element in the variable as the numerator, and the last element in the variable as the denominator.
#' @param WhamTransformation c("count","WhamTransformation"). Default value is "count". Setting as "count" uses original count. Setting as WhamTransformation converts count to ratio(count divided by the total sum of each sample), and then times one million( x 1000,000)
#' @param useMC use multicore by default(FALSE). Multi core processing will be attempted with the BiocParallel package, then the parallel package. If neither are installed, serial processing will be used.
#' @param aldex.module c("t.test","anova"), Default value is t.test(Welch's t and Wilcoxon rank test). "anova" will perform glm and Kruskal Wallace tests for one-way ANOVA
#' @return a dataframe contains statistic results or a WhamResult object contains statistic results and contrast.
#' @return we.ep - Expected P value of Welch's t test
#' @return we.eBH - Expected Benjamini-Hochberg corrected P value of Welch's t test
#' @return wi.ep - Expected P value of Wilcoxon rank test
#' @return wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test
#' @return rab.all: a vector containing the median clr value for each feature
#' @return rab.win.conditionA: a vector containing the median clr value for each feature in condition A
#' @return rab.win.conditionB: a vector containing the median clr value for each feature in condition B
#' @return diff.btw: a vector containing the per-feature median difference between condition A and B
#' @return diff.win: a vector containing the per-feature maximum median difference between Dirichlet instances within conditions
#' @return effect:  median effect size: diff.btw / max(diff.win) for all instances
#' @return overlap: proportion of effect size that overlaps 0 (i.e. no effect)
#' @return kw.ep a vector containing the expected p-value of the Kruskal-Wallis test for each feature
#' @return kw.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @return glm.ep a vector containing the expected p-value of the glm ANOVA for each feature
#' @return glm.eBH a vector containing the corresponding expected value of the Benjamini-Hochberg corrected p-value for each feature
#' @examples
#' dir <-  system.file("extdata","16s_phyloseq.rds",package = "Wham")
#' phy_object = readRDS(dir)
#' wham16s_object = Wham16s(phyloseq_object =  phy_object,
#'                          design = ~Position,
#'                          taxa.level = "otu"
#'                          )
#'
#' # see the vignette for more details
#'
#' @export

Wham16s <- function(taxaCount =NULL,
                    colData = NULL,
                    design= NULL,
                    phyloseq_object = NULL,
                    taxa.level =  c("k","p","o","c","f","g","s","otu"),
                    contrast,
                    ref.contrast,
                    All.contrast,
                    WhamTransformation = c("count","WhamTransformation"),
                    useMC = F,
                    aldex.module = c("t.test","anova")
){
  output = WhamFrom16s(taxaCount = taxaCount,colData = colData,design = design,phyloseq_object = phyloseq_object)

  DE_result = WhamDEFrom16s(object = output,
                            contrast=contrast,
                            ref.contrast = ref.contrast,
                            All.contrast = All.contrast,
                            WhamTransformation = WhamTransformation,
                            taxa.level = taxa.level,
                            useMC = useMC,
                            aldex.module = aldex.module
  )
  se = SummarizedExperiment(assays <- list(assay(output)),colData = data.frame(colData(output)))
  if (!is(se, "RangedSummarizedExperiment")) {
    if (is(se, "SummarizedExperiment")) {
      se <- as(se, "RangedSummarizedExperiment")
    } else {
      stop("'se' must be a RangedSummarizedExperiment object")
    }
  }

  otu_taxainfo = paste0(output@TaxaInfo,";OTU_",1:length(output@TaxaInfo))
  output@TaxaTable = TaxaToDF_fn(otu_taxainfo,level = "otu")

  result <- new("WhamOutput16s",se,
                design = output@design,
                TaxaInfo = output@TaxaInfo,
                TaxaTable = output@TaxaTable,
                DE.result = DE_result)
  return(result)
}
