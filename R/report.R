#' Write a report
#'
#' \code{write_report} Write a report
#'
#' @param outputdir Output directory of the report. String. Default is current working directory.
#' @param filename Filename of the report. String. Default is "report".
#' @param author Author of the report. String. Default is "".
#' @param title Title of the report. String. Default is "Power calculation report for microbiome data".
#' @param datatype Data type needed in the report. Vector. Default is c("species", "metagenomics_pathways","metatranscriptomics_pathways").
#' @param N_sample Number of samples in each data set. Vector. Default is c(200,200,200).
#' @param N_feature Number of features in each data set. Vector. Default is c(300,300,300).
#' @param N_repeat Number of repeat measurements for each data set. Vector. Default is c(2,2,2).
#' @param N_covariate Number of covariates in the study. Number. Default is 5.
#' @param contin The type of covariate we are interested. The type includeds "continuous", "discrete" and "inter" (interaction). Strint. Default is c("continuous","continuous","continuous").
#' @param contin_prop The proportion of continuous covariates. Vector. Default is c(0.5,0.5,0.5)
#' @param pos_prop The proportion of features correlated with the variable we are interested in. Vector. Default is c(0.2,0.2,0.2)
#' @param increment A small number to adjust the smoothness of the power curve. The smaller, the more smooth. Default is 0.005.
#' @param q_cutoff The threshold for q value (p value after BH adjustment). A number in (0,1). Default is 0.05.
#' @param target_power The least power the study wants to reach. A number in (0,1). Default is 0.9.
#' @param target_fpr The highest false positive rate of the study. A number in (0,1). Default is 0.05.
#'
#' @return Generate a .Rnw file and a .pdf file.
#' @author Liu Cao
#' @details
#' Generate a .Rnw file and a .pdf file.
#' @import knitr
#' @export
#'

write_report <- function(outputdir,
                         filename,
                         author,
                         title = "Power calculation report for microbiome data",
                         datatype = c("species", "metagenomics_pathways","metatranscriptomics_pathways"),
                         N_sample = c(200,200,200),
                         N_feature = c(300,300,300),
                         N_covariate = c(5,5,5),
                         N_repeat = c(1,1,1),
                         contin = c("continuous","continuous","continuous"),
                         conti_prop = c(0.5,0.5,0.5),
                         pos_prop = c(0.2,0.2,0.2),
                         increment = 0.005,
                         q_cutoff = 0.05,
                         target_power = 0.9,
                         target_fpr = 0.05,
                         mean_distribution_bugs = ref_dat$mean_distribution_bugs,
                         mean_distribution_pathway_dna = ref_dat$mean_distribution_pathway_dna,
                         mean_distribution_pathway_rna = ref_dat$mean_distribution_pathway_rna,
                         sd_quartile_bugs_common = ref_dat$sd_quartile_bugs[1,],
                         sd_quartile_bugs_rare = ref_dat$sd_quartile_bugs[2,],
                         sd_quartile_pathway_dna_common = ref_dat$sd_quartile_pathway_dna[1,],
                         sd_quartile_pathway_dna_rare = ref_dat$sd_quartile_pathway_dna[2,],
                         sd_quartile_pathway_rna_common = ref_dat$sd_quartile_pathway_rna[1,],
                         sd_quartile_pathway_rna_rare = ref_dat$sd_quartile_pathway_rna[2,],
                         source = c("STAR project","STAR project","STAR project")
){
  setwd(outputdir)
  rnwfilename = paste0(outputdir,filename,".Rnw")
  if(sum(datatype %in% c("species", "metagenomics_pathways","metatranscriptomics_pathways")) != length(datatype) ){
    stop("Misspecified datatype.\n datatype should be 'species', 'metagenomics_pathways' or 'metatranscriptomics_pathways'")
    print("datatype error")
  }
  if(sum(contin %in% c("continuous", "discrete", "inter", "time")) != length(contin)){
    stop("Misspecified variable type. \n contin should be 'continous', 'discrete', 'inter' or 'time'")
  }
  write_setup(rnwfilename)
  write_parameters(rnwfilename,
                   datatype,
                   N_sample,
                   N_covariate,
                   N_feature,
                   N_repeat,
                   conti_prop,
                   pos_prop,
                   contin,
                   increment,
                   q_cutoff,
                   target_power,
                   target_fpr,
                   mean_distribution_bugs,
                   mean_distribution_pathway_dna,
                   mean_distribution_pathway_rna,
                   sd_quartile_bugs_common,
                   sd_quartile_bugs_rare,
                   sd_quartile_pathway_dna_common,
                   sd_quartile_pathway_dna_rare,
                   sd_quartile_pathway_rna_common,
                   sd_quartile_pathway_rna_rare,
                   source )
  write_docmeta(rnwfilename,author,title)
  write_intro(rnwfilename,datatype)
  write_results(rnwfilename,datatype)
  write_methods(rnwfilename,datatype)
  write_ending(rnwfilename)
  knit2pdf(rnwfilename,paste0(outputdir,filename,".tex"))
}

# write basic setups for .Rnw files
write_setup <- function(filename){
  header_path = system.file("extdata", "header.txt", package = "bugpower")
  header = scan(header_path,what="character",sep="\n",quiet=TRUE)
  cat(header,file=filename, sep="\n",append=FALSE)
  cat("",file=filename, sep="\n",append=TRUE)
}

# write parameters
write_parameters <- function(filename,
                             datatype = c("species", "metagenomics_pathways","metatranscriptomics_pathways"),
                             N_sample=c(200,200,200),
                             N_covariate=c(8,8,8),
                             N_feature=c(200,200,200),
                             N_repeat=c(1,1,1),
                             conti_prop=c(0.5,0.5,0.5),
                             pos_prop=c(0.2,0.2,0.2),
                             contin=c("continuous","continuous","continuous"),
                             increment=0.005,
                             q_cutoff=0.05,
                             target_power=0.9,
                             target_fpr=0.05,
                             mean_distribution_bugs = ref_dat$mean_distribution_bugs,
                             mean_distribution_pathway_dna = ref_dat$mean_distribution_pathway_dna,
                             mean_distribution_pathway_rna = ref_dat$mean_distribution_pathway_rna,
                             sd_quartile_bugs_common = ref_dat$sd_quartile_bugs[1,],
                             sd_quartile_bugs_rare = ref_dat$sd_quartile_bugs[2,],
                             sd_quartile_pathway_dna_common = ref_dat$sd_quartile_pathway_dna[1,],
                             sd_quartile_pathway_dna_rare = ref_dat$sd_quartile_pathway_dna[2,],
                             sd_quartile_pathway_rna_common = ref_dat$sd_quartile_pathway_rna[1,],
                             sd_quartile_pathway_rna_rare = ref_dat$sd_quartile_pathway_rna[2,],
                             source = c("STAR project","STAR project","STAR project")
){
  cat("<<parameters, echo=FALSE>>=",sep="\n",file=filename,append=TRUE)

  col1 = c("datatype","N_sample","N_covariate","N_feature","N_repeat",
           "conti_prop","pos_prop","contin","increment",
           "q_cutoff","target_power","target_fpr",
           "mean_distribution_bugs","mean_distribution_pathway_dna","mean_distribution_pathway_rna",
           "sd_quartile_bugs_common","sd_quartile_bugs_rare",
           "sd_quartile_pathway_dna_common","sd_quartile_pathway_dna_rare",
           "sd_quartile_pathway_rna_common","sd_quartile_pathway_rna_rare",
           "source")

  datatype = paste0("c(",paste(paste0("\"",datatype,"\""),collapse=","),")")
  N_sample = paste0("c(",paste(N_sample,collapse=","),")")
  N_covariate = paste0("c(",paste(N_covariate,collapse=","),")")
  N_feature = paste0("c(",paste(N_feature,collapse=","),")")
  N_repeat = paste0("c(",paste(N_repeat,collapse=","),")")
  conti_prop = paste0("c(",paste(conti_prop,collapse=","),")")
  pos_prop = paste0("c(",paste(pos_prop,collapse=","),")")
  contin = paste0("c(",paste(paste0("\"",contin,"\""),collapse=","),")")
  mean_distribution_bugs = paste0("c(",paste(mean_distribution_bugs,collapse=","),")")
  mean_distribution_pathway_dna = paste0("c(",paste(mean_distribution_pathway_dna ,collapse=","),")")
  mean_distribution_pathway_rna = paste0("c(",paste(mean_distribution_pathway_rna ,collapse=","),")")
  sd_quartile_bugs_common = paste0("c(",paste(sd_quartile_bugs_common,collapse=","),")")
  sd_quartile_bugs_rare = paste0("c(",paste(sd_quartile_bugs_rare,collapse=","),")")
  sd_quartile_pathway_dna_common = paste0("c(",paste(sd_quartile_pathway_dna_common,collapse=","),")")
  sd_quartile_pathway_dna_rare = paste0("c(",paste(sd_quartile_pathway_dna_rare,collapse=","),")")
  sd_quartile_pathway_rna_common = paste0("c(",paste(sd_quartile_pathway_rna_common,collapse=","),")")
  sd_quartile_pathway_rna_rare = paste0("c(",paste(sd_quartile_pathway_rna_rare,collapse=","),")")
  source = paste0("c(",paste(paste0("\"",source,"\""),collapse=","),")")

  col2 = c(datatype, N_sample, N_covariate, N_feature, N_repeat,
           conti_prop, pos_prop, contin, increment,
           q_cutoff,target_power,target_fpr,
           mean_distribution_bugs,mean_distribution_pathway_dna,mean_distribution_pathway_rna,
           sd_quartile_bugs_common,sd_quartile_bugs_rare,
           sd_quartile_pathway_dna_common,sd_quartile_pathway_dna_rare,
           sd_quartile_pathway_rna_common,sd_quartile_pathway_rna_rare,
           source
  )

  cat("library(bugpower)",file=filename, sep="\n",append=TRUE)
  cat(paste0(col1," = ", col2), file=filename, sep="\n",append=TRUE)
  #cat("print(length(mean_distribution_pathway_dna))",file=filename,sep="\n",append=TRUE)
  #cat("print(length(mean_distribution_pathway_rna))",file=filename,sep="\n",append=TRUE)
  #cat("print(length(mean_distribution_bugs))",file=filename,sep="\n",append=TRUE)

  cat("names(N_sample) = datatype",file=filename, sep="\n",append=TRUE)
  cat("names(N_covariate) = datatype",file=filename, sep="\n",append=TRUE)
  cat("names(N_feature) = datatype",file=filename, sep="\n",append=TRUE)
  cat("names(N_repeat) = datatype",file=filename, sep="\n",append=TRUE)
  cat("names(conti_prop) = datatype",file=filename, sep="\n",append=TRUE)
  cat("names(pos_prop) = datatype",file=filename, sep="\n",append=TRUE)
  cat("names(contin) = datatype",file=filename, sep="\n",append=TRUE)
  cat("names(source) = datatype",file=filename, sep="\n",append=TRUE)

  cat("@",file=filename, sep="\n",append=TRUE)
  cat("\n",file=filename, sep="\n",append=TRUE)

}

# meta data of the report
write_docmeta <- function(filename,author,title){
  cat(paste0("\\title{",title,"}"),file=filename, sep="\n",append=TRUE)
  cat(paste0("\\author{",author,"}"),file=filename, sep="\n",append=TRUE)
  cat("\\maketitle",file=filename, sep="\n",append=TRUE)
  cat("",file=filename, sep="\n",append=TRUE)
}

# write introduction
write_intro <- function(filename, datatype = c("species", "metagenomics_pathways","metatranscriptomics_pathways")){
  cat("\\section*{Introduction}",file=filename, sep="\n",append=TRUE)

  if("species" %in% datatype){
    text = paste0("\\textbf{Taxonomic composition (Metagenome)}. The microbiome whole metagenome shotgun datasets will comprise ",
                  "\\Sexpr{N_sample[\"species\"]}",
                  " total samples. ", "We estimate detectable effect sizes based on the ",
                  "\\Sexpr{source[\"species\"]}",
                  " in which approximately ",
                  "\\Sexpr{length(mean_distribution_bugs)}",
                  " taxa per sample yielded relative abundances with standard deviation quartiles",
                  " [ \\Sexpr{sd_quartile_bugs_rare} ] ",
                  "among rare taxa (present in less than 50\\% of samples) and ",
                  " [ \\Sexpr{sd_quartile_bugs_common} ] ",
                  "among common taxa (present in more than 50\\% of samples).\n")
    cat(text,file=filename, sep="\n",append=TRUE)
    write_sd(filename, datatype="species")
    cat("\n",file=filename, sep="\n",append=TRUE)
  }

  if("metagenomics_pathways" %in% datatype){
    text = paste0("\\textbf{Pathway composition (Metagenome)}. We estimate detectable relative abundance for metagenomics pathways based on the ",
                  "\\Sexpr{source[\"metagenomics_pathways\"]}",
                  " in which the approximately ",
                  " \\Sexpr{length(mean_distribution_pathway_dna)}",
                  " observed pathways estimated per sample yielded relative abundances standard deviation quartiles",
                  " [ \\Sexpr{sd_quartile_pathway_dna_rare} ] ",
                  "among rare pathways ",
                  " [ \\Sexpr{sd_quartile_pathway_dna_common} ] ",
                  "and common pathways (present in less or more than 50\\% of samples).\n")
    cat(text,file=filename, sep="\n",append=TRUE)
    write_sd(filename, datatype="metagenomics_pathways")
    cat("\n",file=filename, sep="\n",append=TRUE)
  }

  if("metatranscriptomics_pathways" %in% datatype){
    text = paste0("\\textbf{Pathway composition (Metatranscriptome)} The microbiome whole metatranscriptome shotgun datasets will comprise ",
                  "\\Sexpr{N_sample[\"metatranscriptomics_pathways\"]}",
                  " total samples.",
                  "We estimate detectable relative abundance for metatranscriptome pathways based on the ",
                  "\\Sexpr{source[\"metatranscriptomics_pathways\"]}",
                  " in which the approximately ",
                  " \\Sexpr{length(mean_distribution_pathway_rna)}",
                  " observed pathways per sample yielded relative abundances with standard deviation quartiles ",
                  " [ \\Sexpr{sd_quartile_pathway_dna_rare} ] ",
                  "among rare pathways ",
                  " [ \\Sexpr{sd_quartile_pathway_dna_common} ] ",
                  "and common pathways (present in less or more than 50\\% of samples).\n")
    cat(text,file=filename, sep="\n",append=TRUE)
    write_sd(filename, datatype="metatranscriptomics_pathways")
    cat("\n",file=filename, sep="\n",append=TRUE)
  }

  cat("\n",file=filename, sep="\n",append=TRUE)
}

# write result
write_results <- function(filename,
                          datatype = c("species", "metagenomics_pathways","metatranscriptomics_pathways")
){
  cat("\\section*{Results}",file=filename, sep="\n",append=TRUE)
  if("species" %in% datatype){
    write_result_plot(filename,datatype="species")
    write_result_table(filename,datatype="species")
  }
  if("metagenomics_pathways" %in% datatype){
    write_result_plot(filename,datatype="metagenomics_pathways")
    write_result_table(filename,datatype="metagenomics_pathways")
  }
  if("metatranscriptomics_pathways" %in% datatype){
    write_result_plot(filename,datatype="metatranscriptomics_pathways")
    write_result_table(filename,datatype="metatranscriptomics_pathways")
  }

}



# write methods
write_methods <- function(filename,
                          datatype = c("species", "metagenomics_pathways","metatranscriptomics_pathways")
){
  cat("\n",file=filename, sep="\n",append=TRUE)
  cat("\\section*{Method}",file=filename, sep="\n",append=TRUE)

  if("species" %in% datatype){
    species = paste0("For taxa profile, there are ",
                     "\\Sexpr{N_sample[\"species\"]}",
                     " samples and ",
                     "\\Sexpr{N_feature[\"species\"]}",
                     " features. ",
                     "Each sample is measured for ",
                     "\\Sexpr{N_repeat[\"species\"]}",
                     " times. ",
                     "\\Sexpr{pos_prop[\"species\"]}",
                     " of the features are correlated with the covariate we are interested in. ",
                     "\\Sexpr{N_covariate[\"species\"]}",
                     " covariates are included, ",
                     "\\Sexpr{conti_prop[\"species\"]}",
                     " of which are continuous variables following standard normal distribution ($N(0,1)$) and ",
                     "\\Sexpr{1-conti_prop[\"species\"]}",
                     " of which are binary variables following Bernolli distribution ($b(1,0.5)$). The effect size for all the covariates are set to constant 0.5 except for the one we are interested in.\n\n")
    cat(species,file=filename, sep="\n",append=TRUE)
  }
  if("metagenomics_pathways" %in% datatype){
    dna_pathway = paste0("For metagenomics pathways profile, there are ",
                         "\\Sexpr{N_sample[\"metagenomics_pathways\"]}",
                         " samples and " ,
                         "\\Sexpr{N_feature[\"metagenomics_pathways\"]}",
                         " features. ",
                         "Each sample is measured for ",
                         "\\Sexpr{N_repeat[\"metagenomics_pathways\"]}",
                         " times. ",
                         "\\Sexpr{pos_prop[\"metagenomics_pathways\"]}",
                         " of the pathways are correlated with the covariate we are interested in. ",
                         "\\Sexpr{N_covariate[\"metagenomics_pathways\"]}",
                         " covariates are included, ",
                         "\\Sexpr{conti_prop[\"metagenomics_pathways\"]}",
                         " of which are continuous variables following standard normal distribution ($N(0,1)$) and ",
                         "\\Sexpr{1-conti_prop[\"metagenomics_pathways\"]}", "of which are binary variables following Bernolli distribution ($b(1,0.5)$). The effect size for all the covariates are set to constant 0.5 except for the one we are interested in.\n\n")
    cat(dna_pathway,file=filename, sep="\n",append=TRUE)
  }
  if("metatranscriptomics_pathways" %in% datatype){
    rna_pathway = paste0("For metatranscriptomics pathways profile, there are ",
                         "\\Sexpr{N_sample[\"metatranscriptomics_pathways\"]}",
                         " samples and " ,
                         "\\Sexpr{N_feature[\"metatranscriptomics_pathways\"]}",
                         " features. ",
                         "Each sample is measured for ",
                         "\\Sexpr{N_repeat[\"metatranscriptomics_pathways\"]}",
                         " times. ",
                         "\\Sexpr{pos_prop[\"metatranscriptomics_pathways\"]}",
                         " of the pathways are correlated with the covariate we are interested in. ",
                         "\\Sexpr{N_covariate[\"metatranscriptomics_pathways\"]}",
                         " covariates are included, ",
                         "\\Sexpr{conti_prop[\"metatranscriptomics_pathways\"]}",
                         " of which are continuous variables following standard normal distribution ($N(0,1)$) and ",
                         "\\Sexpr{1-conti_prop[\"metatranscriptomics_pathways\"]}", "of which are binary variables following Bernolli distribution ($b(1,0.5)$). The effect size for all the covariates are set to constant 0.5 except for the one we are interested in.\n\n")
    cat(rna_pathway,file=filename, sep="\n",append=TRUE)
  }

  second = "We assume that after arcsin square root transformation, the transformed relative abundance for each feature will follow normal distribution. The mean of the transformed relative abundance are set to sum of the linear combination of the covariates and a tranformed relative abundance sampled from reference data. The standard deviation for each simulated feature are set to the transformed standard deviation quartiles of the reference data according to the sampled relative abundance of each simulated feature.\n\n"

  third = paste0("All the estimations incoprate Benjamini Hochberg (BH) algorithm to control false discovery rate (FDR) to ",
                 "\\Sexpr{q_cutoff}",
                 ". The power and false positive rate (FPR) for each type of feature (defined by common or rare, standard deviation) is calculated as the mean power and FPR within each category.\n\n")
  cat(second,file=filename, sep="\n",append=TRUE)
  cat(third,file=filename, sep="\n",append=TRUE)
}

# write ending
write_ending <- function(filename){
  cat("\\end{document}",file=filename, sep="\n",append=TRUE)
}

# write table
write_sd <- function(filename,datatype){
  if(datatype == "species"){
    text = paste("\\begin{table}[h]",
                 "\\centering",
                 "\\begin{tabular}{c|c|c|c}",
                 "\\hline",
                 "& 25\\% & 50\\% & 75\\% \\\\",
                 "\\hline",
                 "\\Sexpr{paste0(\"Common & \", paste(sd_quartile_bugs_common,collapse=\" & \"))} \\\\",
                 "\\Sexpr{paste0(\"Rare & \", paste(sd_quartile_bugs_rare,collapse=\" & \"))} \\\\",
                 "\\hline",
                 "\\end{tabular}",
                 "\\caption{The quartiles of standard deviations of transformed relative abundance \n for the common and rare taxa.}",
                 "\\end{table}",
                 sep="\n")
    cat(text,file=filename, sep="\n",append=TRUE)
  }
  if(datatype == "metagenomics_pathways"){
    text = paste("\\begin{table}[h]",
                 "\\centering",
                 "\\begin{tabular}{c|c|c|c}",
                 "\\hline",
                 "& 25\\% & 50\\% & 75\\% \\\\",
                 "\\hline",
                 "\\Sexpr{paste0(\"Common & \", paste(sd_quartile_pathway_dna_common,collapse=\" & \"))} \\\\",
                 "\\Sexpr{paste0(\"Rare & \", paste(sd_quartile_pathway_dna_rare,collapse=\" & \"))}  \\\\",
                 "\\hline",
                 "\\end{tabular}",
                 "\\caption{The quartiles of standard deviations of transformed relative abundance \n for the common and rare metagenomics pathways.}",
                 "\\end{table}",
                 sep="\n")
    cat(text,file=filename, sep="\n",append=TRUE)
  }
  if(datatype == "metatranscriptomics_pathways"){
    text = paste("\\begin{table}[h]",
                 "\\centering",
                 "\\begin{tabular}{c|c|c|c}",
                 "\\hline",
                 "& 25\\% & 50\\% & 75\\% \\\\",
                 "\\hline",
                 "\\Sexpr{paste0(\"Common & \", paste(sd_quartile_pathway_rna_common,collapse=\" & \"))}  \\\\",
                 "\\Sexpr{paste0(\"Rare & \", paste(sd_quartile_pathway_rna_rare,collapse=\" & \"))}  \\\\",
                 "\\hline",
                 "\\end{tabular}",
                 "\\caption{The quartiles of standard deviations of transformed relative abundance \n for the common and rare metatranscriptomics pathways.}",
                 "\\end{table}",
                 sep="\n")
    cat(text,file=filename, sep="\n",append=TRUE)
  }
}

# write result
write_result_plot <- function(filename,datatype){
  cat(paste0("<<plot_results_", datatype, ", out.width='.45\\\\textwidth',fig.pos='H',fig.keep='all',message=FALSE,warning=FALSE,error=FALSE,echo=FALSE>>=\n"),file=filename, sep="\n",append=TRUE)
  cat("library(bugpower)",file=filename,sep="\n",append=TRUE)
  cat(paste0("a = calc_power_common_rare(N_sample=N_sample[\"",datatype,"\"],",
             "N_covariate=N_covariate[\"",datatype,"\"],",
             "N_feature=N_feature[\"",datatype,"\"],",
             "N_repeat=N_repeat[\"",datatype,"\"],",
             "conti_prop=conti_prop[\"",datatype,"\"],",
             "pos_prop=pos_prop[\"",datatype,"\"],",
             "contin=contin[\"",datatype,"\"],",
             "increment=increment,
             q_cutoff = q_cutoff,
             quartile_common = sd_quartile_bugs_common,
             quartile_rare = sd_quartile_bugs_rare,
             v_mean_transform = mean_distribution_bugs)\n"),file=filename, sep="\n",append=TRUE)
  #cat("print(a)\n",file=filename, sep="\n",append=TRUE)
  cat("plot_power(a$out_common, title = paste0(datatype,' common'))\n",file=filename, sep="\n",append=TRUE)
  cat("plot_power(a$out_rare, title = paste0(datatype,' rare'))\n",file=filename, sep="\n",append=TRUE)

  cat("a_common = result_power(a$out_common)\n",file=filename, sep="\n",append=TRUE)
  cat("a_rare = result_power(a$out_rare)\n",file=filename, sep="\n",append=TRUE)

  cat("@\n",file=filename, sep="\n",append=TRUE)
  cat("\n",file=filename,append=TRUE)
}

write_result_table <- function(filename,datatype){
  if(datatype == "metagenomics_pathways"){
    datatype = "metagenomics pathways"
  }
  if(datatype == "metatranscriptomics_pathways"){
    datatype = "metatranscriptomics pathways"
  }
  text = paste("\\begin{table}[h]",
               "\\centering",
               "\\begin{tabular}{c|c|c|c}",
               "\\hline",
               "& 25\\% & 50\\% & 75\\% \\\\",
               "\\hline",
               "minimum difference & \\Sexpr{a_common[1,1]} & \\Sexpr{a_common[1,2]} & \\Sexpr{a_common[1,3]} \\\\",
               "fpr & \\Sexpr{a_common[2,1]} & \\Sexpr{a_common[2,2]} & \\Sexpr{a_common[2,3]} \\\\",
               "\\hline",
               "\\end{tabular}",
               paste0("\\caption{Summary for common ", datatype, " at power ","\\Sexpr{target_power}", ".}"),
               "\\end{table}",
               sep="\n")
  cat(text,file=filename, sep="\n",append=TRUE)
  cat("\n",file=filename,append=TRUE)
  text = paste("\\begin{table}[h]",
               "\\centering",
               "\\begin{tabular}{c|c|c|c}",
               "\\hline",
               "& 25\\% & 50\\% & 75\\% \\\\",
               "\\hline",
               "minimum difference & \\Sexpr{a_rare[1,1]} & \\Sexpr{a_rare[1,2]} & \\Sexpr{a_rare[1,3]} \\\\",
               "fpr & \\Sexpr{a_rare[2,1]} & \\Sexpr{a_rare[2,2]} & \\Sexpr{a_rare[2,3]} \\\\",
               "\\hline",
               "\\end{tabular}",
               paste0("\\caption{Summary for rare ", datatype, " at power ","\\Sexpr{target_power}", ".}"),
               "\\end{table}",
               sep="\n")
  cat(text,file=filename, sep="\n",append=TRUE)
  cat("\n",file=filename, sep="\n",append=TRUE)
}







