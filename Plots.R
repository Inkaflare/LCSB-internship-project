library(ggplot2)

# setwd("~/Work/Projects/bioscripts/David Stocklausen")

#The community drug enrichment plot and the necessary functions were made by Marek


# Convenience function for string pasting
"%l%" <- function(a,b) paste(a,b,sep = "")

#enrichment results barplot (with GO categories)
ggplot(ceb_all_community_data[[1]]$CPenrichmentResults, aes(x=Description, y=Count)) +
  geom_bar(stat="identity", aes(fill=p.adjust)) +
  ggtitle("Community 1 (ceb) GSEA results") + 
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
filter_fc_targets <- function(community, quer = "sn_fold_change") {
  bres <- sapply(community[1:(length(community)-2)], function(x) x[["sn_fold_change"]] != 0)
  sapply(community[names(bres)[bres]], function(x) x$drugs)
}

#function to output the drug_enrichment graph for the wanted community set
community_drug_enrichment_plot <- function(all_community_data, community_set_name) {
  fc_drugs <- list()
  all_community_data <- all_community_data[!sapply(all_community_data, is.null)] #remove empty communities

  for(c in all_community_data) {
    # Create a table with drug frequencies
    fc_drugs[[length(fc_drugs) + 1]] <- table(unlist(filter_fc_targets(c)))
    
  }

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
    geom_text(mapping = aes(label=rownames(sumdata)), color = "black", size = 3.5, hjust = -0.2, vjust = -0.2)
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
    ggtitle("Drug and chemical counts per community ("%l% community_set_name %l% ")") + guides(color = F)
}


###Plotting the results from the previous scripts
#drug_enrichment plots (by Marek)

community_drug_enrichment_plot(ceb_all_community_data, "ceb")
community_drug_enrichment_plot(cml_all_community_data, "cml")


#drug chem count plots
drug_chem_counts_plot(ceb_all_community_data, "ceb")
drug_chem_counts_plot(cml_all_community_data, "cml")
