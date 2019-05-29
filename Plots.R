library(ggplot2)

# setwd("~/Work/Projects/bioscripts/David Stocklausen")

#The community drug enrichment plot and the necessary functions were made by Marek


# Convenience function for string pasting
"%l%" <- function(a,b) paste(a,b,sep = "")

#enrichment results barplot (with GO categories)
sorted_com_1_df <- ceb_all_community_data[[1]]$CPenrichmentResults[order(ceb_all_community_data[[1]]$CPenrichmentResults$Count, decreasing = T),]
sorted_com_1_df$Description <- factor(sorted_com_1_df$Description, levels = unique(sorted_com_1_df$Description[order(sorted_com_1_df$Count)]))
ggplot(head(sorted_com_1_df, n=20), aes(x=Description, y=Count)) +
  geom_bar(stat="identity", aes(fill=p.adjust)) +
  ggtitle("Community 1 (ceb) top 20 GSEA terms") + 
  coord_flip()


#for cml
ggplot(cml_all_community_data[[1]]$CPenrichmentResults, aes(x=Description, y=Count)) +
  geom_bar(stat="identity", aes(fill=p.adjust)) +
  ggtitle("Community 1 (cml) GSEA results") +
  coord_flip()

#dotplot of # of drugs affecting each community
#get list of drugs acting on each community
get_community_tally <- function(community, type = c("drugs", "chems", "miRNAs") ) {
  type_list <- c()
  if(length(community)!=0) {
    type_list <- sapply(head(community,-1), function(x) x[[type]])
    type_list <- unique(unlist(type_list))
    target_counts <- length(type_list)
  }
  
  else { target_counts <- 0 }
  return(target_counts)
}



# A function to retrieve the drug targets for genes that have a non-zero fold change in Substantia Nigra (for the enrichment - drug plot)
filter_fc_targets <- function(community) {
  bres <- sapply(community[1:(length(community)-2)], function(x) x[["sn_fold_change"]] != 0)
  sapply(community[names(bres)[bres]], function(x) x$drugs)
}


#modified enrichment plot: consider only genes with fold change for the enrichment axis

filter_cp_enrichment <- function(community, cutoff = 1) {
  enr <- community$CPenrichment
  ###geneID
  bres <- sapply(community[1:(length(community)-2)], function(x) x[["sn_fold_change"]] != 0)
  fcgenes <- names(bres)[bres]
  hits <- NULL
  for(gene in fcgenes) {
    begin <- "^" %l% gene %l% "/"
    middle <- "/" %l% gene %l% "/"
    end <- "/" %l% gene %l% "$"
    gpattern <- begin %l% "|" %l% middle %l% "|" %l% end
    hits <- cbind(hits, (rownames(enr) %in% grep(gpattern, enr$geneID)))
  }
  if(length(hits) == 0) { return(data.frame()) }
  enr <- enr[apply(hits, 1, sum) > cutoff,]
  return(enr)
  
}



#Function for printing community data to files (genes with fold change only!) and upload them as overlays to the PD map for judging their quality/get an overview of the spread of members etc.
print_com <- function(communities, prefix, index, token = NULL) {
  base_url <- "https://pdmap.uni.lu/minerva/api/"
  session_token <- login(base_url = base_url, login = "david.stocklausen", password = "ds1234")
  genes <- unique(names(communities[[index]]))
  write.table(genes, quote = F, row.names = F, col.names = F, sep = "\t",
              file = prefix %l% index %l% "_genes.txt")
  write.table(filter_cp_enrichment(communities[[index]]), quote = F, row.names = F, col.names = T, sep = "\t",
              file = prefix %l% index %l% "_enrichment.txt")
  drugs <- table(unlist(filter_fc_targets(communities[[index]])))
  write.table(cbind(drugs = names(drugs), hits = drugs), quote = F, row.names = F, col.names = T, sep = "\t",
              file = prefix %l% index %l% "_drugs.txt")
  ### Create overlay, if token not null
  if(!is.null(token)) {
    content <- "name\tcolor\n"
    for(g in genes) { content <- content %l% g %l% "\t00FF00\n" }
    upl_data <- "content=" %l% content %l% 
      "&description=" %l% prefix %l% index %l% 
      "&filename=" %l% prefix %l% index %l% ".data" %l% 
      "&name=" %l% prefix %l% index %l%
      "&googleLicenseConsent=true"
    resp <- simple_POST("projects/pd_map_winter_19/overlays/", 
                        token, 
                        base_url = "https://pdmap.uni.lu/minerva/api/", 
                        body = upl_data)
    print(resp)
  }
  return(NULL)
}

#function to output the drug_enrichment graph for the wanted community set (unscaled and unfiltered)
community_drug_enrichment_plot <- function(all_community_data, community_set_name) {
  fc_drugs <- list()
  all_community_data <- all_community_data[!sapply(all_community_data, is.null)] #remove empty communities

  for(c in all_community_data) {
    # Create a table with drug frequencies
    fc_drugs[[length(fc_drugs) + 1]] <- table(unlist(filter_fc_targets(c)))
    
  }

  ##Unfiltered plot
  # Create a data frame for plotting
  sumdata <- data.frame(
              drugs = sapply(fc_drugs, length), # Number of drugs
              size = sapply(all_community_data, length), # Size of the community
              enrichment = sapply(all_community_data, function(x) # Number of enriched terms
                ifelse(length(x[["CPenrichmentResults"]]) > 1, nrow(x[["CPenrichmentResults"]]), 0) # ifelse needed, because for no enrichment there's no data frame, only a string
               ))
  
  # Plot the data frame
  ggplot(data = sumdata, aes(x = size, y = enrichment, size = drugs, color = (drugs == 0))) + # color is a simple transformation to differentiate communities with no drug targets
    geom_point() + 
    ggtitle("Community relative size, enrichment and drug data (" %l% community_set_name %l% ")")+
    # scale_size(trans = "log") +
    scale_size(range=c(2,10)) + 
    labs(color="Affected by drugs", size ="N° of drugs") +
    geom_text(mapping = aes(label=rownames(sumdata)), color = "black", size = 3.5, hjust = -0.2, vjust = -0.2)
  
  

  
  # return(plot1, plot2)
}

filtered_drug_enrichment_plot <- function(all_community_data, community_set_name) {
  
  fc_drugs <- list()
  all_community_data <- all_community_data[!sapply(all_community_data, is.null)] #remove empty communities
  
  for(c in all_community_data) {
    # Create a table with drug frequencies
    fc_drugs[[length(fc_drugs) + 1]] <- table(unlist(filter_fc_targets(c)))
    
  }
  
  ##enrichment filtered + scaled plot
  sumdataCPfiltered <- data.frame(
    drugs = sapply(fc_drugs, length), # Number of drugs
    size = sapply(all_community_data, length), # Size of the community
    enrichment = sapply(all_community_data, function(x) # Number of enriched terms
      nrow(filter_cp_enrichment(x)) # ifelse needed, because for no enrichment there's no data frame, only a string
    ))

  ggplot(data = sumdataCPfiltered, aes(x = size, y = enrichment, size = drugs, color = (drugs == 0))) + # color is a simple transformation to differentiate communities with no drug targets
    geom_point() + 
    ggtitle("Community relative size, enrichment and drug data (" %l% community_set_name %l% ") (filtered + scaled)")+
    scale_size(range=c(2,10)) + 
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10') +
    labs(color="Affected by drugs", size ="N° of drugs") +
    geom_text(mapping = aes(label=rownames(sumdataCPfiltered)), color = "black", size = 3, hjust = 0, vjust = 0)

}

#function to output the drug-chem counts graph for the wanted community set
drug_chem_counts_plot <- function(all_community_data, community_set_name) {
  #count number of drugs and chems affecting each community
  drug_tally <- lapply(all_community_data, get_community_tally, type = "drugs")
  chem_tally <- lapply(all_community_data, get_community_tally, type = "chems")
  drug_chem_tally <- do.call(rbind, Map(data.frame, drugs=drug_tally, chems=chem_tally))
  
  #add a row with community sizes to this for dot size
  drug_chem_tally$community_size <- sapply(all_community_data, function(x) length(x))
  
  #plot it
  ggplot(drug_chem_tally, aes(x = drugs, y=chems, color = "red"), show.legend = F) + 
    geom_point(aes(size = community_size), show.legend = T) +
    geom_text(aes(label=row.names(drug_chem_tally)), color = "black", hjust=0.5, vjust=-1, size = 4) +
    scale_size(range=c(0,10)) + 
    labs(color="Affected by drugs", size ="Community size") +
    ggtitle("Drug and chemical counts per community ("%l% community_set_name %l% ")") + guides(color = F)
}



###Plotting the results from the previous scripts

#drug chem count plots
drug_chem_counts_plot(ceb_all_community_data, "ceb")
drug_chem_counts_plot(cml_all_community_data, "cml")

#drug_enrichment plots (by Marek)
community_drug_enrichment_plot(ceb_all_community_data, "ceb")
community_drug_enrichment_plot(cml_all_community_data, "cml")

filtered_drug_enrichment_plot(ceb_all_community_data, "ceb")
filtered_drug_enrichment_plot(cml_all_community_data, "cml")



