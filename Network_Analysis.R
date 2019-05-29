library(igraph)
library(httr)
library(tidyr)
library(rgl)
library(data.table)
detach(package:dplyr) #some conflicts between igraph and dplyr functions

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


#small scale testing can be done with modelid = 15986 (iron metabolism in the pd_map_spring_18)

#Getting all reactions (=edges for the network)
get_reactions <- function() {
  reactions_request <- GET(paste(api_url, "projects/", project_id, "/models/*/bioEntities/reactions/", sep =""), query = list(columns="lines,centerPoint,reactionId,reactants,modifiers,products"))
  reactions <- content(reactions_request, "parsed")
  reaction_names <- lapply(reactions, function(x) x[[6]])
  reactions <- setNames(reactions, reaction_names)
  return(reactions)
}

#Getting all elements (=nodes/vertices for the network)
get_elements <- function() {
  elements_request <- GET(paste(api_url, "projects/", project_id, "/models/*/bioEntities/elements/", sep =""), query = list(columns="name,bounds,id"))
  elements <- content(elements_request, "parsed")
  return(elements)
}

#combine reactions (edges and node ids) with elements (node details) for further processing
add_element_details <- function(reactions, elements) {  
  object_annotating <- function(object) {
      for (element in elements) {
        if(element[['id']]==object[['aliasId']]) {
          object <- append(object, element)
          return(object)
        }
      }
    }
  for (reaction in 1:length(reactions)) {
    for (i in list("modifiers", "products", "reactants")) {
      if (length(reactions[[reaction]][[i]])!=0) {
        for (object in 1:length(reactions[[reaction]][[i]])) {
          reactions[[reaction]][[i]][[object]] <- object_annotating(reactions[[reaction]][[i]][[object]])
        }
      }
    }
  }
  return(reactions)
}


#create bipartite graphs - add center nodes, then add edges between them and all reactants/products
create_reaction_graph <- function(reaction) {
  reaction_graph <- make_empty_graph(n=0)
  reaction_graph <- add_vertices(reaction_graph, 1, name = reaction[["reactionId"]], type = FALSE, color = "black", shape = "square", size = 7)
  colors <- list(modifiers = "blue", products = "green", reactants = "red")
  for (j in list("modifiers", "products", "reactants")) {
    for (i in 1:length(reaction[[j]])) {
      if (length(reaction[[j]])!=0) {
        reaction_graph <- add_vertices(reaction_graph, 1, name = reaction[[j]][[i]][["name"]], id = reaction[[j]][[i]][["id"]],  type = TRUE, color = colors[[j]], shape = "circle", size = 15)
      }#
    }
  }
   for (vertice in 2:length(V(reaction_graph)$id)) {
    if(V(reaction_graph)[[vertice]]$color%in%list("green")) {
      reaction_graph <- add_edges(reaction_graph, c(1, vertice))
    }
    if(V(reaction_graph)[[vertice]]$color%in%list("blue", "red")) {
      reaction_graph <- add_edges(reaction_graph, c(vertice, 1))
    }
  }
  return(reaction_graph)
}


#join community member lists with overlay target data -> obtain drugs working on communities etc. to discover clusters of same drug targets
get_community_target_data <- function(community_object, overlay_target_data) {
  options(stringsAsFactors=FALSE)
  
  # community_object <- ceb
  # overlay_target_data <- sn_overlay_target_data
  community_list <- community_to_list(community_object)
  community_df_list <- lapply(community_list, do.call, what = rbind.data.frame)
  
  
  joined_community_df <- rbindlist(community_df_list, idcol= 'Community')
  colnames(joined_community_df)[2] <- "Target"
  
  joined_community_df <- merge(joined_community_df, overlay_target_data, by = "Target")
  
  return(joined_community_df)
}


#turn community members into lists we can work with (remove reaction nodes)
community_to_list <- function(community_object) {
  # community_object <- ceb
  com_list <- communities(community_object)
  com_list <- lapply(com_list, function(x) lapply(x, function(y) y[!grepl("re[0-9]+", y)]))
  com_list <- lapply(com_list, function(x) x[lapply(x, length)>0])
  # com_list <- com_list[order(-sapply(com_list, length))]
  return(com_list)
}

#get entrez ids for community members #needed for GSEA
get_entrez_ids <- function(symbol) {
  # if(!grepl(" ", symbol)) {
    # symbol = "SLC11 A2"
    url = utils::URLencode(paste(api_url, "projects/", project_id, "/models/*/bioEntities:search?query=", symbol, "&count=1", sep = ""))
    id_request <- content(GET(url), "parsed")
    id <- id_request[[1]][["id"]]
    entrez_request <- content(GET(paste(api_url, "projects/", project_id, "/models/*/bioEntities/elements/?columns=references&id=", id, sep = "")), "parsed")
    entrez_request <- entrez_request[[1]][["references"]]
    if(length(entrez_request)!=0) {
      for(i in 1:length(entrez_request)) {
        if(all(entrez_request[[i]][["type"]]=="ENTREZ")) {
          return(entrez_request[[i]][["resource"]])
        }
      }
    }
  # }
}

#get overlay up/downregulation data from entrez ids (optional parameter for GSEA)
get_expression_change <- function(symbol, overlay_data, base = 0){
    url = utils::URLencode(paste(api_url, "projects/", project_id, "/models/*/bioEntities:search?query=", symbol, "&count=1", sep = ""))
    id_request <- content(GET(url), "parsed")
    id <- id_request[[1]][["id"]]
    valuedf <- overlay_data %>% filter(Target == symbol) %>% select(value)  #get fold change values from the overlay data frame
    value <- as.double(valuedf$value[1])
    if(is.na(value)) { value <- 0}
    if(base==1) {
      value <- value + 1
    }
    # print(value)
    
    entrez_request <- content(GET(paste(api_url, "projects/", project_id, "/models/*/bioEntities/elements/?columns=references&id=", id, sep = "")), "parsed")
    entrez_request <- entrez_request[[1]][["references"]]
    if(length(entrez_request)!=0) {
      for(i in 1:length(entrez_request)) {
        if(all(entrez_request[[i]][["type"]]=="ENTREZ")) {
          return(list("entrez_id" = entrez_request[[i]][["resource"]], "value"=value))
        }
      }
    }
  # }
}

#turn fold change list into proper genelist format as used by clusterProfiler
format_to_genelist <- function(expression_change_list) {
  names(expression_change_list) <- lapply(expression_change_list, function(x) as.character(x[[1]]))
  expression_change_list <- lapply(expression_change_list, function(x) x[[2]])
  return(sort(unlist(expression_change_list, use.names = T), decreasing = T))
}



##Function calls

reactions <- get_reactions()
elements <- get_elements()

reactionsAnnotated <- add_element_details(reactions, elements)

reaction_graph_list <- lapply(reactionsAnnotated, create_reaction_graph)

# plot(reaction_graph_list[[27]])


#Join the graphs into a single network
joined_graph <- do.call(union, reaction_graph_list)

  #adding the attributes back into the joined graph (reaction nodes -> small squares)
joined_attributes <- as.list(get.vertex.attribute(joined_graph))
shape_attributes <- grep("^shape_\\d+$", names(joined_attributes), value = TRUE)
V(joined_graph)$shape <- Reduce(f = function(a, b) ifelse(is.na(a), b, a), x = joined_attributes[shape_attributes])

vertex_sizes  <- numeric(length = length(V(joined_graph)))
for(vertex in 1:length(V(joined_graph))){
  if(V(joined_graph)[[vertex]]$shape=="square"){
    vertex_sizes[[vertex]] <- 0.5 
    # print("0.5")
  } else {
    vertex_sizes[[vertex]] <- 1
    # print("1")
  }
}
  

l = layout_with_graphopt(joined_graph)
plot(joined_graph, vertex.size = 5*vertex_sizes, vertex.label = 0.5, edge.arrow.size = 0.5, layout = l, main = project_id)


#Network analysis
#Density
edge_density(joined_graph)

#Diameter (longest distance in all shortest paths between two nodes)
diameter(joined_graph)

#degree - cross-correlate with regulation targets?
deg <- degree(joined_graph, mode = "all")


hist(deg[!grepl("re[0-9]+", names(deg))], breaks=1:200-1, main="Histogram of node degree")
deg.dist <- degree_distribution(joined_graph, cumulative=T, mode="all")
plot( x=0:max(deg), y=1-deg.dist, pch=19, cex=1.2, col="orange", xlab="Degree", ylab="Cumulative Frequency")

#targets sorted by degree + up/downreg value -> possibly interesting hub nodes + the drugs that apply to them
sn_overlay_data$degree <- deg[sn_overlay_data$Target]
sn_overlay_data <- sn_overlay_data[order(-sn_overlay_data$degree, na.last = NA),]


#closeness - seems fairly pointless
closeness <- closeness(joined_graph, mode="all", weights=NA)
# closeness <- data.frame(key = names(closeness), values = closeness)
# colnames(closeness) <- c("node", "score")
# closeness <- closeness[order(-closeness$score),]

sn_overlay_data$closeness <- closeness[sn_overlay_data$Target]


#betweenness - number of nodes passing through each node
betweenness <- betweenness(joined_graph)

sn_overlay_data$betweenness <- betweenness[sn_overlay_data$Target]
# sn_overlay_data <- sn_overlay_data[order(-sn_overlay_data$betweenness, -abs(sn_overlay_data$value), na.last = NA),]

#scale closeness and betweenness to have more reasonable values there instead of millions and millionths
sn_overlay_data$closeness <- scale(sn_overlay_data$closeness)
sn_overlay_data$betweenness <- scale(sn_overlay_data$betweenness)


#create a table of quantiles for the network metrics I'm interested in
centrality_quantile_table <- data.frame("Degree" = quantile(sn_overlay_data$degree, na.rm = T), 
                                        "Betweenness" = quantile(sn_overlay_data$betweenness, na.rm = T), 
                                        "Closeness" = quantile(sn_overlay_data$closeness, na.rm = T))



#betweenness and degree plot
ggplot(sn_overlay_data, aes(x=degree, y = betweenness)) +
  geom_point() +
  geom_smooth(method = "lm", se =F) +
  ggtitle("Degree and betweenness of network nodes")

# #hubs and authorities
# hubscore <- hub_score(joined_graph)$vector
# sn_overlay_data$hubscore <- hubscore[sn_overlay_data$Target]


#Saving to file after adding network centrality metrics
setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB")
save(sn_overlay_data, file = "sn_overlay_data.Rda")
write.table(sn_overlay_data, file = "sn_overlay_data.tab", sep = "\t", quote = F, row.names = F)


#check whether the network is scale-free
fit_power_law(degree(joined_graph_undir), xmin = 1)

#finding subnetworks
#turn into undirected network for community detection
joined_graph_undir <- as.undirected(joined_graph, mode ="collapse") %>% simplify

l = layout_with_graphopt(joined_graph_undir)
plot(joined_graph_undir, vertex.size = 5*vertex_sizes, vertex.label = NA, arrow.size = 0.01, layout = l, main = project_id)


ceb <- cluster_edge_betweenness(joined_graph_undir)
# plot_dendrogram(ceb, main = "Communities (edge betweenness) - Dendrogram" )
length(ceb)
sizes(ceb)
plot(ceb, joined_graph_undir, vertex.size = 5*vertex_sizes, vertex.label= NA, arrow.size = 0.01, main = "Communities (edge betweenness)", layout = l)
membership(ceb)
modularity(ceb)


clp <- cluster_label_prop(joined_graph_undir)
plot(clp, joined_graph_undir, vertex.size = 5*vertex_sizes, vertex.label=NA, arrow.size = 0.01, main = "Communities (label propagation)", layout = l)


cfg <- cluster_fast_greedy(joined_graph_undir)
plot(cfg, joined_graph_undir, vertex.size = 5*vertex_sizes, vertex.label= NA, arrow.size = 0.01, main = "Communities (fast greedy)", layout = l)
# plot_dendrogram(cfg)

cwt <- cluster_walktrap(joined_graph_undir)
plot(cwt, joined_graph_undir, vertex.size = 5*vertex_sizes, vertex.label= NA, arrow.size = 0.01, main = "Communities (walktrap)", layout = l)

cle <- cluster_leading_eigen(joined_graph_undir)
plot(cle, joined_graph_undir, vertex.size = 5*vertex_sizes, vertex.label= NA, arrow.size = 0.01, main = "Communities (leading eigenvector)", layout = l)

cml <- multilevel.community(joined_graph_undir)
plot(cml, joined_graph_undir, vertex.size = 5*vertex_sizes, vertex.label= NA, arrow.size = 0.01, main = "Communities (multilevel)", layout = l)

ceb_joined_data <- get_community_target_data(ceb, sn_overlay_data)

#for GO analysis:

#edge betweenness
#create community lists, filter out the reaction nodes
ceb_list <- community_to_list(ceb)

#get the entrez ids for analysis with DAVID
ceb_ids <- lapply(ceb_list, function(x) purrr::compact(lapply(x, get_entrez_ids)))

# lapply(ceb_ids[[8]], write, "test.txt", append=TRUE, ncolumns=1000)

library(dplyr)
#get the entrez ids + expression change for analysis with DAVID & ClusterProfiler
ceb_ids_fc <- lapply(ceb_list, function(x) purrr::compact(lapply(x, get_expression_change, overlay_data = sn_overlay_data)))

ceb_ids_genelist <- lapply(ceb_ids_fc, format_to_genelist)


#multilevel
cml_list <- community_to_list(cml)

#get the entrez ids for analysis with DAVID
cml_ids <- lapply(cml_list, function(x) purrr::compact(lapply(x, get_entrez_ids)))

# lapply(ceb_ids[[8]], write, "test.txt", append=TRUE, ncolumns=1000)

library(dplyr)
#get the entrez ids + expression change for analysis with DAVID & ClusterProfiler
cml_ids_fc <- lapply(cml_list, function(x) purrr::compact(lapply(x, get_expression_change, overlay_data = sn_overlay_data)))

cml_ids_genelist <- lapply(cml_ids_fc, format_to_genelist)



#fast greedy
cfg_list <- community_to_list(cfg)

#get the entrez ids for analysis with DAVID
cfg_ids <- lapply(cfg_list, function(x) purrr::compact(lapply(x, get_entrez_ids)))


#fast greedy
cwt_list <- community_to_list(cwt)

#get the entrez ids for analysis with DAVID
cwt_ids <- lapply(cwt_list, function(x) purrr::compact(lapply(x, get_entrez_ids)))





#Saving to file after adding network centrality metrics
setwd("C:/Users/David/Dropbox/Universiteitsdokumenter/Stage LCSB")
save(sn_overlay_data, file = "sn_overlay_data.Rda")
write.table(sn_overlay_data, file = "sn_overlay_data.tab", sep = "\t", quote = F, row.names = F)

load("sn_overlay_data.Rda")
