library(RDAVIDWebService)
library(clusterProfiler)
library(AnnotationHub)
library(org.Hs.eg.db)
detach(package:igraph)


##DAVID analysis


david <- DAVIDWebService$new(email="david.stocklausen@ext.uni.lu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

#generate Functional Annotation Cluster Report file for a community
getFunctionalAnnotationClusterReport <- function(community, index, community_set) {
  wd <- getwd()
  setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB/DAVID Result Tables")
  
  if(length(community)!=0) { #skip empty communities to avoid crashing
    setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))
    listname <- paste(community_set, index, sep="_")
    # print(listname)
    result <- addList(david, community, "ENTREZ_GENE_ID", listName = listname, listType = "Gene")
  
    # termCluster <- getClusterReport(david)
    filename <- paste0("termClusterReport", "_", community_set, "_", index, ".tab")
    getClusterReportFile(david, filename)
    filename <- paste0("functionalAnnotationChart", "_", community_set, "_", index, ".tab")
    getFunctionalAnnotationChartFile(david, filename, threshold = 0.1)
    print(paste("Report file for community ", index, " of set ", community_set, " saved successfully"))
  } else { NA }
  setwd(wd)
}


#Cluster reports
#cfg_ids
# mapply(getFunctionalAnnotationClusterReport, community = cfg_ids, index = 1:length(cfg_ids), MoreArgs = list(community_set="cfg_ids"))

#ceb_ids
mapply(getFunctionalAnnotationClusterReport, community = ceb_ids, index = 1:length(ceb_ids), MoreArgs = list(community_set="ceb_ids"))
# mapply(getFunctionalAnnotationClusterReport, community = ceb_ids_genelist, index = 1:length(ceb_ids), MoreArgs = list(community_set="ceb_ids"))

#cml_ids
mapply(getFunctionalAnnotationClusterReport, community = cml_ids, index = 1:length(cml_ids), MoreArgs = list(community_set="cml_ids"))


#Plots
summary(termCluster)
plot2D(termCluster, 3)
davidGOdag <- DAVIDGODag(members(termCluster)[[3]], pvalueCutoff=0.1, "CC")
plotGOTermGraph(g=goDag(davidGOdag), r = davidGOdag, max.nchar = 40, node.shape = "ellipse")




##ClusterProfiler

#enrichGO function (without fold change) (over-representation test)
getEnrichmentCategories <- function(community, index, community_set_name) {
  if(length(community)!=0) { #skip empty communities to avoid crashing
    if(!file.exists(paste("enrichmentCategoriesResult_", community_set_name,"_", index, ".tab", sep=""))){ #skip existing files to save time; delete them in the folder if you want to regenrate them!
    # community = ceb_ids
    # community_set_name = "ceb_ids"
    # index = 1
      index <- as.numeric(index)
      genelist <- get("community")
      # print(index)
      ego <- tryCatch(
        enrichGO(gene     = genelist,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 keyType = "ENTREZID",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 qvalueCutoff  = 0.05
        ),
         error = function(e) {
          error <- paste(paste("Community ", index, ": ", e, sep = ""))
          print(error)
          write(error, file = "Communities where enrichment categories analysis fails.txt", append = T, sep = "")
        }
      )
        # head(ego)
        if(class(ego)[1]=="enrichResult") {
          ego <- setReadable(ego, OrgDb =  org.Hs.eg.db)
          print(paste("enrichmentCategoriesResult_", community_set_name,"_", index, ".tab", sep=""))
          write.table(ego, 
                      file = paste("enrichmentCategoriesResult_", community_set_name,"_", index, ".tab", sep=""), 
                      sep = "\t", 
                      quote = F, 
                      row.names = F)
        }
          # else {print("Enrichment analysis failed for community ", index, sep = "")}
    }
  }
}

#gseGO function (with fold change)
getEnrichmentAnalysis <- function(community, index, community_set_name) {
  if(length(community)!=0) { #skip empty communities to avoid crashing
    if(!file.exists(paste("enrichmentAnalysisResult_", community_set_name,"_", index, ".tab", sep=""))){ #skip existing files to save time; delete them in the folder if you want to regenerate them!
      # community = ceb_ids_genelist
      # community_set_name = "ceb_ids"
      # index = 1
      index <- as.numeric(index)
      # print(index)
      # str(index)
      genelist <- get("community")
      # print(genelist)
      ego <- tryCatch(
        gseGO(gene     = genelist,
                 OrgDb        = org.Hs.eg.db,
                 ont          = "BP",
                 keyType = "ENTREZID",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.05,
                 minGSSize = 5,
                 nPerm = 500
        )
        , error = function(e) {
                    error <- paste(paste("Community ", index, ": ", e, sep = ""))
                    print(error)
                    write(error, file = "Communities where GSEA fails.txt", append = T, sep = "")
                }
       )
      # head(ego)
      if(class(ego)[1]=="gseaResult") {
        ego <- setReadable(ego, OrgDb =  org.Hs.eg.db)
        print(paste("enrichmentAnalysisResult_", community_set_name,"_", index, ".tab", sep=""))
        write.table(ego, 
                    file = paste("enrichmentAnalysisResult_", community_set_name,"_", index, ".tab", sep=""), 
                    sep = "\t", 
                    quote = F,
                    row.names = F)
      }
      # else {print("Enrichment analysis failed for community ", index, sep = "")}
    }
  }
}




perform_ClusterProfiler_analysis <- function(community, community_fc, community_set_name) {
  wd <- getwd()
  setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB/ClusterProfiler Result Tables")
  
  print("Performing Enrichment Categories analysis")
  if(file.exists("Communities where enrichment categories analysis fails.txt")) {
    file.remove("Communities where enrichment categories analysis fails.txt")
  }
  mapply(getEnrichmentCategories, community = community, index = seq_along(community), MoreArgs = list(community_set_name = community_set_name))
  
  # print("Performing GSEA")
  # if(file.exists("Communities where GSEA fails.txt")) {
  #   file.remove("Communities where GSEA fails.txt")
  # }
  # mapply(getEnrichmentAnalysis, community = community_fc, index = seq_along(community_fc), MoreArgs = list(community_set_name = community_set_name))
  setwd(wd)
}

perform_ClusterProfiler_analysis(ceb_ids, ceb_ids_genelist, "ceb_ids")

perform_ClusterProfiler_analysis(cml_ids, cml_ids_genelist, "cml_ids")


# #make a table of the communities and enrichment value 1 for marek (testing)
# create_table <- function(community_set, index, community_set_name) {
#   current_wd <- getwd()
#   setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB/Tables for Marek")
#   # community_set <- ceb_ids[[1]]
#   community_table <- data.frame(unlist(get("community_set")))
#   community_table[["value"]] <- rep("1", length(community_set))
#   # print(community_table)
#   colnames(community_table)[[1]] <- c("entrez gene")
#   # print(community_table)
#   write.table(community_table, file = paste(community_set_name, "_", index, ".tab", sep=""), sep = "\t", row.names = F, col.names = T, quote = F)
#   setwd(current_wd)
# }
# 
# #tables for marek
# mapply(create_table, community_set = ceb_ids, index = seq_along(ceb_ids), MoreArgs = list(community_set_name = "ceb_ids"))



##Merging community data and target data into a nested data structure of communities, their members, drugs etc. that work on them, and enrichment results

get_community_member_data <- function(community_members, index, all_targets){
  gene_info <- get("community_members")[index]
  gene <- list(symbol=AnnotationDbi::select(org.Hs.eg.db, keys = names(gene_info), columns = "SYMBOL", keytype = "ENTREZID")[[2]],
                     entrezID = names(gene_info), 
                     sn_fold_change=unlist(gene_info[[1]]))
  ab_fold_change <- c(unique(ab_overlay_data[ab_overlay_data$Target==gene$symbol, 2]))
  if(length(ab_fold_change)==0) {ab_fold_change <- c(0)}
  drugs <- c(unique(all_targets[(all_targets$Target==gene$symbol & all_targets$Type=="Drug" ), 2]))
  chems <- c(unique(all_targets[(all_targets$Target==gene$symbol & all_targets$Type=="Chemical" ), 2]))
  miRNAs <- c(unique(all_targets[(all_targets$Target==gene$symbol & all_targets$Type=="miRNA" ), 2]))
  gene <- c(gene, list(ab_fold_change=ab_fold_change), list(drugs = drugs), list(chems=chems), list(miRNAs=miRNAs))

  return(gene)
}

#Read David Cluster reports and make them into a list of lists of clusters per community set
#code taken from Marek and made into a function
get_david_clusters <- function(community_set_name) {
  
  setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB/DAVID Result Tables/")
  DAVID_tabs <- shell("cd \"C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB/DAVID Result Tables/\" && dir /b /a-d termClusterReport_" %l%
                        community_set_name %l% "*.tab"
                      ,
                      intern = T)
  
  comm_david <- list()
  for(dt in DAVID_tabs) {
    dtl <- readLines(dt)
    stends <- cbind(grep("^Annotation Cluster", dtl), grep("^$", dtl))
    clusters <- list()
    if(length(stends) > 0) {
      for(i in 1:nrow(stends)) {
        ftbl <- read.table(file = dt, 
                           skip = stends[i,1], nrows = stends[i,2]-stends[i,1]-2,
                           sep = "\t", header = T, stringsAsFactors = F)
        clusters[[length(clusters) + 1]] <- list(name = dtl[stends[i,1]], cluster = ftbl)
      } 
    }
    comm_david[[length(comm_david) + 1]] <- clusters
  }
  names(comm_david) <- sapply(strsplit(DAVID_tabs, split = "\\.|_"), `[[`, 4)
  comm_david <- comm_david[order(as.numeric(names(comm_david)))]
  return(comm_david)
}

get_community_data <- function(community_set, index, all_targets, community_set_name) {
  genelist <- get(paste(community_set_name, "genelist", sep="_"))[index]

  wd <- getwd()
  setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB/ClusterProfiler Result Tables")
  
  if(length(genelist[[1]])!=0) {
    #add the community member data
    set_community_data <- mapply(get_community_member_data, community_members=genelist, index = seq_along(genelist[[1]]), MoreArgs = list(all_targets=all_targets), SIMPLIFY = F)
    names(set_community_data) <- sapply(set_community_data, function(x) x[[1]])
    
    #add the enrichment result tables
    #ClusterProfiler
    setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB/ClusterProfiler Result Tables")
    if(file.exists(paste("enrichmentCategoriesResult_", community_set_name, "_", index, ".tab", sep=""))) {
      CPenrichmentResults <- read.table(paste("enrichmentCategoriesResult_", community_set_name, "_", index, ".tab", sep=""), header = T, sep = "\t")
    } else { 
      CPenrichmentResults <- "ClusterProfiler Enrichment Analysis unavailable for community"
    }
    
    #DAVID; reads the cluster reports list and picks out the one for the specific community
    comm_david <- get(community_set_name %l% "_DAVID")[[as.character(index)]]
    if(length(comm_david)==0) {
      comm_david <- "DAVID Enrichment Analysis unavailable for community"
    }

    
    set_community_data <- c(set_community_data, list(CPenrichmentResults = CPenrichmentResults, DAVIDenrichmentClusters = comm_david))
    # print(set_community_data)
    return(set_community_data)
  } else { return(NULL)}
  setwd(wd)
}

ceb_ids_DAVID <- get_david_clusters("ceb")
cml_ids_DAVID <- get_david_clusters("cml")

ceb_all_community_data <- mapply(get_community_data, community_set = ceb_ids, index = seq_along(ceb_ids), MoreArgs=list(all_targets=all_targets, community_set_name="ceb_ids"), SIMPLIFY = F)

cml_all_community_data <- mapply(get_community_data, community_set = cml_ids, index = seq_along(cml_ids), MoreArgs=list(all_targets=all_targets, community_set_name="cml_ids"), SIMPLIFY = F)


setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB")
save(ceb_all_community_data, file = "ceb_all_community_data.RData")
save(cml_all_community_data, file = "cml_all_community_data.RData")


