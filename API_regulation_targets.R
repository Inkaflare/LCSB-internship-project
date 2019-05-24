library(httr)
library(tidyr)
library(dplyr)
library(stringr)

#logging in - obtaining session cookie
login <- "dstocklausen"
password <- "zohajo45"

#login <- "anonymous"  #for users without their own credentials
#password <- ""

#api_url <- "https://curation-sandbox.uni.lu/minerva/api/"
#api_url <- "https://minerva-dev.lcsb.uni.lu/minerva/api/"
api_url <- "http://10.240.6.156/minerva/api/" #use IP in URL to bypass the firewall that makes repeated requests whiff and thus return an incomplete dataset
project_id <- "pd_map_spring_18"

login_request <- POST(paste(api_url, "doLogin/", sep = ""), query = list(login = login, password = password))
content(login_request, "text") #confirm if login was succesful


###Functions

#id request function for adding names to the ids of elements (used in drug/miRNA/chemical search)
id_requesting <- function(id) {
  id_request <- content(GET(paste(api_url, "projects/", project_id, "/models/*/bioEntities/elements/?columns=name&id=", id, sep="")), "parsed")
  return(id_request[[1]][[1]])
}

#Collecting drug targets

drug_search <- function(drug) {
  print(drug)

  try({
  drugsearch_request <- GET(paste(api_url, "projects/", project_id, "/drugs:search", sep =""), query = list(query=drug), timeout(100))
  drug_targets = content(drugsearch_request, "parsed")
  stop_for_status(drugsearch_request)

  #trimming list down to only resource names + ids and converting to dataframe
    if(length(drug_targets[[1]])==8) {
    if(length(drug_targets[[1]][[8]])!=0) {
      drug_targets <- drug_targets[[1]][[8]]
      #remove sublists with empty targetParticipants
      drug_targets_cleaned <- drug_targets[c(sapply(drug_targets, function(x) length(x[[4]])!=0&length(x[[3]])!=0))]
      if(length(drug_targets_cleaned)==0) { return() }
      drug_targets_cleaned <- lapply(drug_targets_cleaned, function(x) x[[3]])
      drug_targets_cleaned <- unlist(drug_targets_cleaned, recursive = FALSE)

      drug_target_names <- lapply(drug_targets_cleaned, function(x) id_requesting(x[["id"]]) )
      drug_targets_cleaned <- lapply(drug_targets_cleaned, function(x) x[[1]])
      drug_targets <- data.frame(cbind(drug_target_names), cbind(drug_targets_cleaned))
      
      colnames(drug_targets) <- c(paste("Drug targets for", drug, sep =" "), "id")
      
      return(drug_targets) #returns the list of targets for the input drug
    }
  }
  })
}

#Collecting miRNA targets

miRNA_search <- function(miRNA) {

  try({
  miRNAsearch_request <- GET(paste(api_url, "projects/", project_id, "/miRnas:search", sep =""), query = list(query=miRNA), timeout(100))
  miRNA_targets <- content(miRNAsearch_request, "parsed")
  stop_for_status(miRNAsearch_request)
  
  if(length(miRNA_targets[[1]])==3) {
    if(length(miRNA_targets[[1]][[3]])!=0) {
      miRNA_targets <- miRNA_targets[[1]][[3]]
      print(miRNA)
      miRNA_targets_cleaned <- miRNA_targets[c(sapply(miRNA_targets, function(x) length(x[[4]])!=0&length(x[[3]])!=0))]
      if(length(miRNA_targets_cleaned)==0) { return() }

      miRNA_targets_cleaned <- lapply(miRNA_targets_cleaned, function(x) x[[3]])
      miRNA_targets_cleaned <- unlist(miRNA_targets_cleaned, recursive = FALSE)
      
      miRNA_target_names <- lapply(miRNA_targets_cleaned, function(x) id_requesting(x[["id"]]) )
      miRNA_targets_cleaned <- lapply(miRNA_targets_cleaned, function(x) x[[1]])
      miRNA_targets <- data.frame(cbind(miRNA_target_names), cbind(miRNA_targets_cleaned))
      
      colnames(miRNA_targets) <- c(paste("Targets for", miRNA, sep =" "), "id")
      
      return(miRNA_targets)
    }
  }
  })
}


#Collecting chemical targets

chem_search <- function(chem) {

  try( {
  chemsearch_request <- GET(paste(api_url, "projects/", project_id, "/chemicals:search", sep =""), query = list(query=chem), timeout(100))
  chem_targets <- content(chemsearch_request, "parsed")
  stop_for_status(chemsearch_request)
  
  if(length(chem_targets)!=0) {
    if(length(chem_targets[[1]])==8) {
      if(length(chem_targets[[1]][[8]])!=0) {
        chem_targets <- chem_targets[[1]][[8]]
        print(chem)
        chem_targets_cleaned <- chem_targets[c(sapply(chem_targets, function(x) length(x[[4]])!=0&length(x[[3]])!=0))]
        if(length(chem_targets_cleaned)==0) {return() }
        
        chem_targets_cleaned <- lapply(chem_targets_cleaned, function(x) x[[3]])
        chem_targets_cleaned <- unlist(chem_targets_cleaned, recursive = FALSE)
        
        chem_target_names <- lapply(chem_targets_cleaned, function(x) id_requesting(x[["id"]]) )
        chem_targets_cleaned <- lapply(chem_targets_cleaned, function(x) x[[1]])
        chem_targets <- data.frame(cbind(chem_target_names), cbind(chem_targets_cleaned))
        
        colnames(chem_targets) <- c(paste("Targets for", chem, sep =" "), "id")
        
        return(chem_targets)
      }
    }
  }

  })
}



#Getting names of everything in the map to search for potential drugs/chemicals/miRNAs that target them

get_element_names <- function() {
  element_request <- GET(paste(api_url, "projects/", project_id, "/models/*/bioEntities/elements/?columns=name,id", sep =""))
  element_names = content(element_request, "parsed")
  stop_for_status(element_request)
  element_names <- data.frame(do.call(rbind.data.frame, element_names))
  element_names <- element_names[match(unique(element_names$name), element_names$name),] #pick first occurence id of each unique name
  return(element_names)
}


#create list of drugs that could act on each target
get_single_query_drug_list <- function(targetID) {
  try( {
  target <- "ALIAS"
  drugsearch_request <- GET(paste(api_url, "projects/", project_id, "/drugs:search?columns=name&target=", target, ":", targetID, sep=""), times = 5, timeout(100))
  drug_list = content(drugsearch_request, "parsed")
  stop_for_status(drugsearch_request)
  drug_list <- lapply(drug_list, toupper)
  return(unique(drug_list))
  } )
} 

get_query_drug_list <- function() {
  element_names <- get_element_names()
  drug_list <- lapply(element_names$id, get_single_query_drug_list)  
  drug_list <- unlist(drug_list)
  drug_list <- lapply(drug_list, function(x) x)
  return(unique(drug_list))

}

#create list of miRNAs that could act on each target
get_single_query_miRNA_list <- function(targetID) {  
  try( {
  target <- "ALIAS"
  miRNAsearch_request <- GET(paste(api_url, "projects/", project_id, "/miRnas:search?columns=name&target=", target, ":", targetID, sep=""), times = 5, timeout(100))
  miRNA_list = content(miRNAsearch_request, "parsed")
  stop_for_status(miRNAsearch_request)
  miRNA_list <- lapply(miRNA_list, tolower)
  return(unique(miRNA_list))
  })
} 

get_query_miRNA_list <- function() {
  element_names <- get_element_names()
  miRNA_list <- lapply(element_names$id, get_single_query_miRNA_list)  
  miRNA_list <- unlist(miRNA_list)
  miRNA_list <- lapply(miRNA_list, function(x) x)
  return(unique(miRNA_list))
}


#create list of chemicals that could act on each target
#for this one the suggestion feature in MINERVA works so the function is simpler
get_query_chem_list <- function() {
  chem_suggestions <- content(GET(paste(api_url, "projects/", project_id, "/chemicals/suggestedQueryList", sep=""), times = 5, timeout(100)), "parsed") 
  return(unique(chem_suggestions))
}


#turn the list of dataframes of targets into a single large dataframe
targets_dataframes_merge <- function(target_list) {
  colnames = c("Target", "id")
  target_list <- lapply(target_list, setNames, colnames)
  target_list <- unique(bind_rows(target_list, .id = 'source'))
  return(target_list)
}




#Collecting overlay information

get_overlay_list <- function() {
  overlay_list <- content(GET(paste(api_url, "projects/", project_id, "/overlays/?publicOverlay=true", sep = "")), "parsed")
  overlay_list <- lapply(overlay_list, function(x) x[c("name", "idObject")])
  overlay_list <- data.frame(do.call(rbind.data.frame, overlay_list))
  return(overlay_list)
}

#Collecting data from a specific overlay (use above list to find the id)

get_overlay_content <- function(overlay_id) {
  overlay_request <- GET(paste(api_url, "projects/", project_id, "/overlays/", overlay_id, "/models/*/bioEntities/", sep =""), times = 5, timeout(100))
  overlay_content = content(overlay_request, "parsed")
  stop_for_status(overlay_request)
  overlay_content <- lapply(overlay_content, function(x) x[[1]][c("idObject", "value")])
  overlay_content <- data.frame(do.call(rbind.data.frame, overlay_content))
  return(overlay_content)
}



#join the lists of drugs/chemicals/miRNAs and bioEntities (up/down regulated proteins) 
# -> determine drugs that act on multiple proteins for example

overlay_merge <- function(overlay_id, targets) {
  # overlay_id <- "7608"
  overlay_content <- get_overlay_content(overlay_id)
  #Joining target and overlay data
  targets[,4] <- as.numeric(targets[,4])
  overlay_content[,2] <- as.numeric(overlay_content[,2])
  overlay_target_data <- merge(targets, overlay_content, by.x = 'id', by.y = 'idObject')
  overlay_target_data <- unique(overlay_target_data[c(4,2,3,5)])
  
  #sort by maximum up/down regulation
  overlay_target_data <- overlay_target_data[order(-abs(overlay_target_data$value)),]
  
  return(overlay_target_data)
}

#Make the merged target dataframe  more compact and readable (spreading)
clean_overlay_target_data <- function(target_data) {
  target_data <- sn_overlay_data
  target_data[,1] <- unlist(target_data[,1])
  target_data <- mutate(target_data, i=row_number())
  
  # failsafe for when one of the columns(drugs/chems/miRNAs) would be empty, and keep this column order
  columns_to_summarize <- c(unique(unlist(target_data[,2])))
  columns_to_summarize <- columns_to_summarize[order(factor(columns_to_summarize, levels = c("Drug", "Chemical", "miRNA")))]
  
  target_data <- spread(target_data, key=Type, value=source) %>% group_by(Target, value) %>% summarise_at(columns_to_summarize, paste, collapse =", ")
  
  remove_NA <- function(df, column) {
    df[[column]] <- gsub("NA, ", "", df[[column]])
    df[[column]] <- gsub(", NA", "", df[[column]])
    df[[column]] <- gsub("NA", "", df[[column]])
    return(df)
  }
  
  for (column in columns_to_summarize) {
    target_data <- remove_NA(target_data, column)
  }
  
  target_data <- apply(target_data, 2, function(x) gsub("^$|^ $", NA, x))  #replace "NA" characters by actual NA data type
  target_data <- as.data.frame(target_data)
  target_data[,2] <- as.numeric(paste(target_data[,2])) #transform up/downreg. values from factor into numeric values
  
  # #failsafe for keeping this column order (needed for the target counts function)
  # col_order <- c("Target", "value")
  # if("Drug" %in% colnames(target_data)) {
  #   col_order <- c(col_order, "Drug")  
  # }
  # if("Chemical" %in% colnames(target_data)) {
  #   col_order <- c(col_order, "Chemical")    }
  # if("miRNA" %in% colnames(target_data)) {
  #   col_order <- c(col_order, "miRNA")  
  # }
  # target_data <- target_data[, col_order]
  
  return(target_data)
}


#tally the total number of targeters per target (if their column exists)
count_targets <- function(target_data) {
  # target_data <- sn_overlay_data
  if("Drug" %in% colnames(target_data)) {
    target_data <- mutate(target_data, total_drugs = str_count(target_data[,"Drug"], ",") + 1)  
  }
  if("Chemical" %in% colnames(target_data)) {
    target_data <- mutate(target_data, total_chems = str_count(target_data[,"Chemical"], ",") + 1)  
  }
  if("miRNA" %in% colnames(target_data)) {
    target_data <- mutate(target_data, total_miRNAs = str_count(target_data[,"miRNA"], ",") + 1)  
  }
  return(target_data)
}


####Function calls


#Drug search
drug_list <- get_query_drug_list()

drug_search_results <- lapply(drug_list, drug_search)
drug_search_results <- drug_search_results[!sapply(drug_search_results, function(x) class(x)!="data.frame")]
drug_search_results <- drug_search_results[unlist(lapply(drug_search_results, length)!=0)]
names(drug_search_results) <- sapply(drug_search_results, function(x) colnames(x)[1]<-paste(unlist(strsplit(colnames(x)[1], split=' '))[-(1:3)], collapse=" "))

#miRNA search
miRNA_list <- get_query_miRNA_list()

miRNA_search_results <- lapply(miRNA_list, miRNA_search)
miRNA_search_results <- miRNA_search_results[!sapply(miRNA_search_results, function(x) class(x)!="data.frame")]
miRNA_search_results2 <- miRNA_search_results[unlist(lapply(miRNA_search_results, length)!=0)]
names(miRNA_search_results) <- sapply(miRNA_search_results, function(x) colnames(x)[1]<-paste(unlist(strsplit(colnames(x)[1], split=' '))[-(1:2)], collapse=" "))

#chemicals
chem_list <- get_query_chem_list()

chem_search_results <- lapply(chem_list, chem_search)
chem_search_results <- chem_search_results[!sapply(chem_search_results, function(x) class(x)!="data.frame")]
chem_search_results <- chem_search_results[unlist(lapply(chem_search_results, length)!=0)]
names(chem_search_results) <- sapply(chem_search_results, function(x) colnames(x)[1]<-paste(unlist(strsplit(colnames(x)[1], split=' '))[-(1:2)], collapse=" "))



 
#Combining into one large dataframe of drug/miRNA/chemical targets
all_targets <- list(drug_search_results, miRNA_search_results, chem_search_results)
names(all_targets) <- c("Drug", "miRNA", "Chemical")
all_targets <- bind_rows(lapply(all_targets, targets_dataframes_merge), .id='Type')




#Getting bioentity data from overlay

overlay_list <- get_overlay_list()

sn_overlay_id <- "7608" #id of the substantia nigra overlay

sn_overlay_data <- overlay_merge(sn_overlay_id, all_targets)
sn_overlay_data <- clean_overlay_target_data(sn_overlay_data)

sn_overlay_data <- count_targets(sn_overlay_data)

#I'm working with the substantia nigra overlay data for the res of the project, but you can also grab the data from the ageing brain overlay

ab_overlay_id <- "7607" #id of the ageing brain overlay
ab_overlay_data <- overlay_merge(ab_overlay_id, all_targets)
ab_overlay_data <- clean_overlay_target_data(ab_overlay_data)

#Saving to file
setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB")
save(sn_overlay_data, file = "sn_overlay_data.Rda")
write.table(sn_overlay_data, file = "sn_overlay_data.tab", sep = "\t", quote = F, row.names = F)

load("sn_overlay_data.Rda")
