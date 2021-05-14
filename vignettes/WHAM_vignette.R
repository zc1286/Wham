## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----echo=FALSE, message=FALSE, warning=FALSE---------------------------------
library(Wham)
dir <-  system.file("extdata","biobakery_sample_input.tsv",package = "Wham")
countdata <-  read.delim(dir)
knitr::kable(countdata[1:5,1:5])

## ----echo=FALSE---------------------------------------------------------------
dir.feature <-  system.file("extdata","EBI_feature_input.tsv",package = "Wham")
featuretable <-  read.delim(dir.feature)
dir.taxa <-  system.file("extdata","EBI_taxa_input.tsv",package = "Wham")
taxatable <-  read.delim(dir.taxa)
knitr::kable(featuretable [1:5,1:5])

## ----echo=FALSE---------------------------------------------------------------
taxatable_short = taxatable[1:5,1:5]
taxatable_short$Taxa = gsub("archaeota|nobacteria|anobacteriales|anobacteriaceae|anobrevibacter|ethanobrevibacter|ethanosphaera","",taxatable_short$Taxa)
knitr::kable(taxatable_short)

## ----fig.height=5, fig.width=5, message=FALSE---------------------------------
dir <-  system.file("extdata","16s_phyloseq.rds",package = "Wham")
phy_object = readRDS(dir)
wham16s_object = WhamFrom16s(phyloseq_object =  phy_object)

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
dir.data <-  system.file("extdata","biobakery_sample_input.tsv",package = "Wham")
countData <-  read.delim(dir.data)
countData <- countData[1:2000,] #we subset the count table to reduce computational cost, 

dir.meta = system.file("extdata","biobakery_sample_metadata.csv",package = "Wham")
metadata <- read.csv(dir.meta,row.names = 1)


 wham_bbk <-  WhamBiobakery(countData = countData,
                         colData = metadata,
                         DE = "taxa",  ##required,choose taxa or feautre(gene family)
                         design = ~ Location,    ##design formula to let function conduct tests on the group of interest,
                         taxa.level = "s",  ##collapse the bacterial rank level when choosing DE as "taxa", default:"s"("speices")."otu" when analyze 16s,
                         contrast = c("Location","Stool","Arm")  ##specify the comparison of interest: Stool(numerator) vs Arm(denominator) in 'Location' variable in metadata.
                         )

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
wham_bbk_anova <-  WhamBiobakery(countData = countData,
                         colData = metadata,
                         DE = "taxa",  ##required,choose taxa or feautre(gene family)
                         design = ~Location,    ##design formula to let function conduct tests on the group of interest,
                         taxa.level = "s",  ##collapse the bacterial rank level when choosing DE as "taxa", default:"s"("speices")."otu" when analyze 16s,
                         aldex.module = "anova"
                         )

## ----message=FALSE, warning=FALSE---------------------------------------------
set.seed(123)
metadata$States = sample(c("LA","NY","CO"),47,replace = T)

 wham_bbk_two_factor <-  WhamBiobakery(countData = countData,
                         colData = metadata,
                         DE = "feature",
                         design = ~ Location + States,    ##complex design 
                         )

## -----------------------------------------------------------------------------
VolcanoPlot(wham_bbk)

## -----------------------------------------------------------------------------
TaxaBarPlot(wham_bbk,
            filter = "all", 
            display.number = 30, ## only visualize top 30 abundant bacteria
            taxa.level = "g", ## visualiazation will be collapsed at "genues level"
)

## -----------------------------------------------------------------------------
TaxaBarPlot(wham_bbk,
            filter = "DE.filter", ## visualize based on statistics analysis
            taxa.level = "g", ## all bacteria whose passed the statistic test will be collapsed at "genus level"
            p.cutoff = 0.05,
            effect.size.range = c(-1,1), ## setting appropriate range cutoff, the default range is c(0,0)
            merge_group = T)

## -----------------------------------------------------------------------------
TaxaHeatmap(wham_bbk,
            taxa.level = "s",
            filter = "DE.filter",
            effect.size.range = c(-0.5,0.5),
            scale = T, ## default setting
            column_split = c(rep("Arm",12),rep("Stool",12)), ## column_split is one of argument in Heatmap in ComplexHeatmap and compitle with other arguments; our demonstration comparison only contains 12 Arm samples and 11 Stool samples.
            row_names_gp = grid::gpar(fontsize = 5)
            )

## -----------------------------------------------------------------------------
TaxaHeatmap(wham_bbk,
            taxa.level = "s",
            filter = "DE.filter",
            effect.size.range = c(-0.5,0.5),
            scale = F, ## Using raw count
            show_row_names = FALSE,
            show_column_names = T)

## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
 wham_bbk_feature <-  WhamBiobakery(countData = countData,
                         colData = metadata, ##required, choose feautre(gene family)
                         DE = "feature",
                         design = ~ Location,    ##design formula to let function conduct tests on the group of interest,
                         ref.contrast = c("Location","Arm")
                         )

FeatureHeatmap(wham_bbk_feature,
               filter = "DE.filter",
               scale = T,
               column_names_gp = grid::gpar(fontsize = 7))

## ----fig.height=6, fig.width=6------------------------------------------------
FeatureCorr(wham_bbk_feature,
            UseDE.result = T, ## we use results from differential analysis
            effect.size.range = c(-1.5,1.5),
            )
            

## ----echo=TRUE, fig.height=6, fig.width=6-------------------------------------
feature_selection = wham_bbk_feature@DE.result$name[1:10] ##select first 10 pathways  
sample_selection = rownames(wham_bbk@colData)[1:25] ##select fist 25 samples
FeatureCorr(wham_bbk_feature,
            feature = feature_selection,
            group = sample_selection
            )

## ----echo=TRUE----------------------------------------------------------------
sessionInfo()

