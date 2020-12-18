###note
#sep_sym = c("(\\.|;|,)")




#######################TaxaToDF_fn####################################
####converting taxa character to data frame. Missing taxa are filled with 'unassiagned'
#' Workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @description converting taxa character to data frame. Missing taxa are filled with 'unassiagned'
#' @param taxainfo a vector store taxonomic information
#' @param level_table decide to use 16s or shotgun sequencing
#' @return return a data frame contain taxonomic information
#'
#' @examples taxainfo = c("Unclassified",
#' @examples              "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales",
#' @examples              "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae")
#' @examples TaxaToDF_fn(taxainfo)
#'
#' @export
TaxaToDF_fn = function(taxainfo,level_table = c("not_otu","otu")){

  if(length(grep("k__|p__|c__|o__|f__|g__|s__|unassigned",taxainfo))==0)
    stop("the taxa info has to contain at least one of following:k__,p__,c__,o__,f__,g__,s__ or unassigned")

  #######identify seperate symbol and unassigned symbol
  sep_sym = c("(\\.|;|,)")
  unassigned_sym = "unassigned"

    ######create table
    TaxaTable = c()

    level_table = match.arg(level_table,choices = c("not_otu","otu"))

    for( i in 1:length(taxainfo)){
      if(level_table == "not_otu"){
        if(length(grep("sk__",as.character(taxainfo[i]))) == 0){ superkingdom_level = unassigned_sym}else{
          superkingdom_level = gsub(paste0("^.*","sk","__"),"",as.character(taxainfo[i]) )
          superkingdom_level = paste0("sk__",superkingdom_level)
          superkingdom_level = gsub( paste0(sep_sym,"k","__",".*$"),"",superkingdom_level,ignore.case = F)
        }

        if(length(grep(paste0("^k__|",sep_sym,"k__"),as.character(taxainfo[i]))) == 0){ kingdom_level = unassigned_sym}else{
          kingdom_level = gsub( paste0("^.*","k","__"),"",as.character(taxainfo[i]) )
          kingdom_level = paste0("k__",kingdom_level)
          kingdom_level = gsub( paste0(sep_sym,"p","__",".*$"),"",kingdom_level)
        }

        if(length(grep("p__",as.character(taxainfo[i]))) == 0){ phylum_level = unassigned_sym}else{
          a = gsub( paste0("^.*","p","__"),"",as.character(taxainfo[i]) )
          a = paste0("p__",a)
          phylum_level = gsub( paste0(sep_sym,"c","__",".*$"),"",a)
        }

        if(length(grep("c__",as.character(taxainfo[i]))) == 0){ class_level = unassigned_sym}else{
          a = gsub( paste0("^.*","c","__"),"",as.character(taxainfo[i]) )
          a = paste0("c__",a)
          class_level = gsub( paste0(sep_sym,"o","__",".*$"),"",a)
        }

        if(length(grep("o__",as.character(taxainfo[i]))) == 0){ order_level = unassigned_sym}else{
          a = gsub( paste0("^.*","o","__"),"",as.character(taxainfo[i]) )
          a = paste0("o__",a)
          order_level = gsub( paste0(sep_sym,"f","__",".*$"),"",a)
        }

        if(length(grep("f__",as.character(taxainfo[i]))) == 0){ family_level = unassigned_sym}else{
          a = gsub( paste0("^.*","f","__"),"",as.character(taxainfo[i]) )
          a = paste0("f__",a)
          family_level = gsub( paste0(sep_sym,"g","__",".*$"),"",a)
        }

        if(length(grep("g__",as.character(taxainfo[i]))) == 0){ genus_level = unassigned_sym}else{
          a = gsub( paste0("^.*","g","__"),"",as.character(taxainfo[i]) )
          a = paste0("g__",a)
          genus_level = gsub( paste0(sep_sym,"s","__",".*$"),"",a)
        }

        if(length(grep("s__",as.character(taxainfo[i]))) == 0){ species_level = unassigned_sym}else{
          species_level = gsub( paste0("^.*","s","__"),"",as.character(taxainfo[i]) )
          species_level = paste0("s__",species_level)
        }

        taxaTable = data.frame(superkingdom = superkingdom_level)
        taxaTable$kingdom = kingdom_level
        taxaTable$phylum = phylum_level
        taxaTable$class = class_level
        taxaTable$order = order_level
        taxaTable$family = family_level
        taxaTable$genus = genus_level
        taxaTable$species = species_level
      }


      if(level_table == "otu"){
        if(length(grep("OTU",as.character(taxainfo[i]))) == 0)
          stop("OTUs are not found")

        otu_level = gsub( paste0("^.*","OTU","_"),"",as.character(taxainfo[i]) )
        otu_level = paste0("OTU_",otu_level)

        taxainfo[i] =  gsub( paste0(sep_sym,"OTU","_",".*$"),"",taxainfo[i])

        if(length(grep("sk__",as.character(taxainfo[i]))) == 0){ superkingdom_level = unassigned_sym}else{
          superkingdom_level = gsub(paste0("^.*","sk","__"),"",as.character(taxainfo[i]) )
          superkingdom_level = paste0("sk__",superkingdom_level)
          superkingdom_level = gsub( paste0(sep_sym,"k","__",".*$"),"",superkingdom_level,ignore.case = F)
        }

        if(length(grep(paste0("^k__|",sep_sym,"k__"),as.character(taxainfo[i]))) == 0){ kingdom_level = unassigned_sym}else{
          kingdom_level = gsub( paste0("^.*","k","__"),"",as.character(taxainfo[i]) )
          kingdom_level = paste0("k__",kingdom_level)
          kingdom_level = gsub( paste0(sep_sym,"p","__",".*$"),"",kingdom_level)
        }

        if(length(grep("p__",as.character(taxainfo[i]))) == 0){ phylum_level = unassigned_sym}else{
          a = gsub( paste0("^.*","p","__"),"",as.character(taxainfo[i]) )
          a = paste0("p__",a)
          phylum_level = gsub( paste0(sep_sym,"c","__",".*$"),"",a)
        }

        if(length(grep("c__",as.character(taxainfo[i]))) == 0){ class_level = unassigned_sym}else{
          a = gsub( paste0("^.*","c","__"),"",as.character(taxainfo[i]) )
          a = paste0("c__",a)
          class_level = gsub( paste0(sep_sym,"o","__",".*$"),"",a)
        }

        if(length(grep("o__",as.character(taxainfo[i]))) == 0){ order_level = unassigned_sym}else{
          a = gsub( paste0("^.*","o","__"),"",as.character(taxainfo[i]) )
          a = paste0("o__",a)
          order_level = gsub( paste0(sep_sym,"f","__",".*$"),"",a)
        }

        if(length(grep("f__",as.character(taxainfo[i]))) == 0){ family_level = unassigned_sym}else{
          a = gsub( paste0("^.*","f","__"),"",as.character(taxainfo[i]) )
          a = paste0("f__",a)
          family_level = gsub( paste0(sep_sym,"g","__",".*$"),"",a)
        }

        if(length(grep("g__",as.character(taxainfo[i]))) == 0){ genus_level = unassigned_sym}else{
          a = gsub( paste0("^.*","g","__"),"",as.character(taxainfo[i]) )
          a = paste0("g__",a)
          genus_level = gsub( paste0(sep_sym,"s","__",".*$"),"",a)
        }

        if(length(grep("s__",as.character(taxainfo[i]))) == 0){ species_level = unassigned_sym}else{
          species_level = gsub( paste0("^.*","s","__"),"",as.character(taxainfo[i]) )
          species_level = paste0("s__",species_level)
        }

        taxaTable = data.frame(superkingdom = superkingdom_level)
        taxaTable$kingdom = kingdom_level
        taxaTable$phylum = phylum_level
        taxaTable$class = class_level
        taxaTable$order = order_level
        taxaTable$family = family_level
        taxaTable$genus = genus_level
        taxaTable$species = species_level
        taxaTable$OTU = otu_level
        }
      TaxaTable = rbind(TaxaTable,taxaTable)
    }

    ########remove superkingdom in case of no superkingdom assignment
    if(all(TaxaTable$superkingdom == unassigned_sym)){
      TaxaTable = TaxaTable[,-1]
    }

  return(TaxaTable)
}

#######################TaxaIden_fn####################################

#' Workflow Hub for Automated Metagenomic Exploration (WHAM!)
#' @description check the format of taxainfo and output taxa.LEVEL
#' @param TaxaInfo a vector store taxonomic information
#' @param taxa.level required, c("k","p","o","c","f","g","s"). Collapsing the count to provided taxonomic level. options are "k","p","c","o","f","g","s".Only need to declaimed when argument DE is setup as "taxa". Every taxonimic level has to have "k__","p__","c__","o__","f__","g__",or "s__"
#' @return a character vector corresponds with taxa.level
#'
#' @examples taxainfo = c("Unclassified",
#' @examples              "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales",
#' @examples              "k__Archaea;p__Euryarchaeota;c__Methanobacteria;o__Methanobacteriales;f__Methanobacteriaceae")
#' @examples TaxaIden_fn(taxainfo,taxa.level = "f")
#'
#' @export

TaxaIden_fn <- function(TaxaInfo,taxa.level =  c("k","p","o","c","f","g","s","otu")){

  taxa.level=match.arg(taxa.level,choices = c("k","p","o","c","f","g","s","otu"))

  sep_sym = c("(\\.|;|,)")

  if(taxa.level == "1"| taxa.level == "k"){
    check_taxa = grep(paste0("^k__|",sep_sym,"k__"),TaxaInfo)
    if(identical(check_taxa, integer(0)))
      stop("there is no kingdom level bacteria, the name has to start with 'k__' or '1__' ")
  }

  if(taxa.level == "2"| taxa.level == "p"){
    check_taxa = grep("p__",TaxaInfo)
    if(identical(check_taxa, integer(0)))
      stop("there is no phylum level bacteria, the name has to start with 'p__' or '2__' ")
  }

  if(taxa.level == "3"| taxa.level == "c"){
    check_taxa = grep("c__",TaxaInfo)
    if(identical(check_taxa, integer(0)))
      stop("there is no class level bacteria, the name has to start with 'c__' or '3__' ")
  }

  if(taxa.level == "4"| taxa.level == "o"){
    check_taxa = grep("o__",TaxaInfo)
    if(identical(check_taxa, integer(0)))
      stop("there is no order level bacteria, the name has to start with 'o__' or '4__' ")
  }

  if(taxa.level == "5"| taxa.level == "f"){
    check_taxa = grep("f__",TaxaInfo)
    if(identical(check_taxa, integer(0)))
      stop("there is no family level bacteria, the name has to start with 'f__' or '5__' ")
  }

  if(taxa.level == "6"| taxa.level == "g"){
    check_taxa = grep("g__",TaxaInfo)
    if(identical(check_taxa, integer(0)))
      stop("there is no genus level bacteria, the name has to start with 'g__' or '6__' ")
  }

  if(taxa.level == "7"| taxa.level == "s"){
    check_taxa = grep("s__",TaxaInfo)
    if(identical(check_taxa, integer(0)))
      stop("there is no species level bacteria, the name has to start with 's__' or '7__' ")
  }

  if( taxa.level == "otu"){
    check_taxa = grep("s__",TaxaInfo)
    if(identical(check_taxa, integer(0)))
      stop("there is no species level bacteria, the name has to start with 's__' or '7__' ")
  }



  if(taxa.level == "1"| taxa.level == "k"){taxaLEVEL = "kingdom"}
  if(taxa.level == "2"| taxa.level == "p"){taxaLEVEL = "phylum"}
  if(taxa.level == "3"| taxa.level == "c"){taxaLEVEL = "class"}
  if(taxa.level == "4"| taxa.level == "o"){taxaLEVEL = "order"}
  if(taxa.level == "5"| taxa.level == "f"){taxaLEVEL = "family"}
  if(taxa.level == "6"| taxa.level == "g"){taxaLEVEL = "genus"}
  if(taxa.level == "7"| taxa.level == "s"){taxaLEVEL = "species"}
  if(taxa.level == "otu"){taxaLEVEL = "OTU"}
  return(taxaLEVEL)
}

