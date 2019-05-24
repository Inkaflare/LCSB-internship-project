# LCSB-internship-project
Panel search for drugs, chemicals and miRNA affecting genes involved in Parkinson's disease.

This git holds the scripts written as part my project as well as the relevant generated data from them.
The R script files correspond to the major parts of the project (see report). The script files and their corresponding result files are:

1. target search: API_regulation.targets.R  
  * dataframe of all up- or downregulated targets and the targeters affecting them: sn_overlay_data.Rda
  
2. network analysis and community detection: Network_analysis.R  
  Two community sets are used in further steps: ceb (edge betweenness) and cml (multilevel)

3. GSEA: Community_enrichment_analysis.R  
  * enrichment analysis reports in folders:
    * ClusterProfiler Result Tables (for ceb and cml) using CC and BP ontologies
    * DAVID Result Tables (for ceb and cml) for all three ontologies (BP; CC; MF)  
      * functional annotation chart: raw GSEA data  
      * termCluster reports: grouped results  

4. Combining the previous data and plotting the results: Plots.R
  * ceb_all_community_data: combining all target and GSEA data per community into a single dataframe  
  * cml_all_community_data: same as above for cml  
  
