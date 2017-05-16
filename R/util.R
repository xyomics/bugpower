#' @import multtest
#' @export
adjp <- function(x,proc="BH",...){
  tmp <- mt.rawp2adjp(x,proc=proc,...)
  tmp <- tmp$adjp[order(tmp$index),]
  return(tmp)
}

#'@export
target_value <- function(v_power,v_fpr,v_beta,v_diff,target_power=0.9){
  index2 = which(v_power > target_power)[1]
  index1 = tail(which(v_power<= target_power),n=1)
  lambda = (target_power - v_power[index1]) / (v_power[index2] - v_power[index1])
  target_beta = v_beta[index1] + (v_beta[index2] - v_beta[index1])*lambda
  target_diff = v_diff[index1] + (v_diff[index2] - v_diff[index1])*lambda
  target_fpr = v_fpr[index1] + (v_fpr[index2] - v_fpr[index1])*lambda
  return(matrix(c(target_diff, target_fpr),nrow=2))
}

#'@export
is_common <- function(otu_df){
  num_row = dim(otu_df)[1]
  num_col = dim(otu_df)[2]
  otu_barcode = otu_df > 0.0
  common_index = apply(otu_barcode,1,mean) > 0.5
  return(common_index)
}

#'Generate mean distribution and quartiles of standard deviation
#'
#' \code{calc_ref_distribution} calculate reference distribution
#'
#'@param file The full path of the reference data. String.
#'@param common_rare An indicator whether to distinguish the features to be common and rare features. Bool. Default is TRUE.
#'@param sep The seperator of the reference data. String. Default is \\t.
#'
#'@return A list of two or three element: mean_distribution and (sd_quartile or sd_quartile_common, sd_quartile_rare).
#'@author Liu Cao
#'@details
#'Given a reference dataset, perform arcsin square root transformation, and calculate its mean distribution, standard deviation of common and rare.
#'The reference dataset should contain relative abundance.
#'@export
#'

calc_ref_distribution <- function(file,common_rare=TRUE,sep="\t"){
  dat = read.csv(file=file,header=TRUE,row.names=1,sep=sep)
  dat = asin(sqrt(dat))
  mean_distribution = unname(apply(dat,1,mean))
  if(common_rare){
    dat_common = dat[is_common(dat),]
    dat_rare = dat[!is_common(dat),]
    sd_quartile_common = summary(apply(dat_common,1,sd))[c(2,3,5)]
    sd_quartile_rare = summary(apply(dat_rare,1,sd))[c(2,3,5)]
    return(list(mean_distribution = mean_distribution,
                sd_quartile_common=sd_quartile_common,
                sd_quartile_rare=sd_quartile_rare))
  }
  else{
    sd_quartile = summary(apply(dat,1,sd))[c(2,3,5)]
    return(list(mean_distribution = mean_distribution,
                sd_quartile = sd_quartile))
  }
}

#'@export
get_p <- function(lm.obj){
  return(lm.obj$coefficients[2,4])
}
