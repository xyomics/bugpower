#' Power calculation for a given study design
#'
#' \code{calc_power} calculates the power species and pathways.
#'
#' @param N_sample Number of subjects in the study. Integer.
#' @param N_covariate Number of covariates in the study. Integer.
#' @param N_feature Number of expected species or pathways. Integer.
#' @param N_repeat Number of repeat measurements for each data set. Integer.
#' @param conti_prop The proportion of continuous covariates. A number in [0,1]. Default is 0.5.
#' @param pos_prop The proportion of species or pathways that are truely correlated with the covariate we are interested in. A number in [0,1]. Default is 0.2.
#' @param contin The type of covariate we are interested. The type includeds "continuous", "discrete", "inter" (interaction), "time". Strint. Default is "continuous".
#' @param increment A small number to adjust the smoothness of the power curve. The smaller, the more smooth. Default is 0.005.
#' @param q_cutoff The threshold for q value (p value after BH adjustment). A number in (0,1). Default is 0.05.
#' @param quartile_sd Standard deviation quartiles
#' @param v_mean_transform Transformed mean distribution
#'
#' @return A data frame with power and beta.
#' @author Liu Cao
#' @details
#' Returning a data frame containing power, FPR and Beta
#' @export

calc_power <- function(N_sample=200,
                       N_covariate=8,
                       N_feature=300,
                       N_repeat = 1,
                       conti_prop=0.5,
                       pos_prop=0.2,
                       contin="continuous",
                       increment=0.005,
                       q_cutoff = 0.05,
                       quartile_sd,
                       v_mean_transform){
  # return power and corresponding effect size
  #print("return power and corresponding effect size")
  n_continuous = ceiling(conti_prop*N_covariate)
  n_discrete = N_covariate - n_continuous

  n_feature_pos = ceiling(N_feature*pos_prop)
  n_feature_neg = N_feature - n_feature_pos

  power_current_min = 0
  power_current = c()
  fpr_current = c()
  beta_interest = 0
  power_out = c()
  fpr_out = c()
  beta_out = c()
  diff_out = c()
  num_rep = 100

  #print("iteration")
  while(power_current_min < 0.99){
    qvalue_matrix = matrix(rep(0,num_rep*N_feature),nrow=num_rep)
    start.time <- Sys.time()
    #print(power_current_min)
    for( rep.i in 1:num_rep ){
      ## initialize p-value vector for each replicate
      #print("initialize p-value vector for each replicate")
      pvalue_v = c()
      ## simulate covariates
      #print("simulate covariates")
      dat_continuous = matrix(rnorm(n_continuous*N_sample,mean=0,sd=1),nrow=N_sample)
      dat_discrete = matrix(rbinom(n_discrete*N_sample,size=1,prob=0.5),nrow=N_sample)
      if(N_repeat > 1){
        dat_time = rep(1,N_sample)
        dat_continuous_single = dat_continuous
        dat_discrete_single = dat_discrete
        for(i in 2:N_repeat){
          dat_time = c(dat_time,rep(i,N_sample))
          dat_continuous = rbind(dat_continuous,dat_continuous_single)
          dat_discrete = rbind(dat_discrete,dat_discrete_single)
        }
      }

      ## assign coefficient for positive features and negative features separately
      #print("assign coefficient for positive features and negative features separately")
      beta_continuous_neg = rep(0.5,n_continuous)
      beta_discrete_neg = rep(0.5,n_discrete)
      beta_continuous_pos = rep(0.5,n_continuous)
      beta_discrete_pos = rep(0.5,n_discrete)

      ## assign coefficient for covariate of interest
      #print("assign coefficient for covariate of interest")
      if(contin == "continuous"){
        beta_continuous_pos[1] = beta_interest
        beta_continuous_neg[1] = 0
        if(N_repeat == 1){
          beta_all_pos = c(beta_continuous_pos,beta_discrete_pos)
          beta_all_neg = c(beta_continuous_neg,beta_discrete_neg)
          dat_all = cbind(dat_continuous,dat_discrete)
        }
        else{
          beta_all_pos = c(beta_continuous_pos,beta_discrete_pos,0.5)
          beta_all_neg = c(beta_continuous_neg,beta_discrete_neg,0.5)
          dat_all = cbind(dat_continuous,dat_discrete,dat_time)
        }
      }
      if(contin == "discrete"){
        beta_discrete_pos[1] = beta_interest
        beta_discrete_neg[1] = 0
        if(N_repeat == 1){
          beta_all_pos = c(beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(beta_discrete_neg,beta_continuous_neg)
          dat_all = cbind(dat_discrete,dat_continuous)
        }
        else{
          beta_all_pos = c(beta_discrete_pos,beta_continuous_pos,0.5)
          beta_all_neg = c(beta_discrete_neg,beta_continuous_neg,0.5)
          dat_all = cbind(dat_discrete,dat_continuous,dat_time)
        }
      }
      if(contin == "inter"){
        if(N_repeat == 1){
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg)
          inter = dat_discrete[,1]*dat_discrete[,2]
          dat_all = cbind(inter,dat_discrete,dat_continuous)
        }
        else{
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos,0.5)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg,0.5)
          inter = dat_discrete[,1]*dat_discrete[,2]
          dat_all = cbind(inter,dat_discrete,dat_continuous,dat_time)
        }
      }
      if(contin == "time"){
        if(N_repeat == 1){
          stop("No repeat measure for the samples.\n")
        }
        else{
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg)
          dat_all = cbind(dat_time,dat_discrete,dat_continuous)
        }
      }

      ## sample mean
      #print("sample mean")
      v_mean = sample(v_mean_transform,N_feature,replace=TRUE)
      index_pos = sample(1:length(v_mean),n_feature_pos,replace=FALSE)
      v_mean_pos = v_mean[1:length(v_mean) %in% index_pos]
      v_mean_neg = v_mean[!(1:length(v_mean) %in% index_pos)]

      ## calculate p-value for positive features
      #print("calculate p-value for positive features")
      dat_all_df = data.frame(dat_all)
      fix_pos = c(dat_all %*% beta_all_pos)
      v_sd_pos = c(rep(quartile_sd[1],floor(length(v_mean_pos)/3)),
                   rep(quartile_sd[2],floor(length(v_mean_pos)/3)),
                   rep(quartile_sd[3],length(v_mean_pos)-floor(length(v_mean_pos)/3)*2))
      Y_pos = matrix(rep(v_mean_pos,length(fix_pos)),nrow=length(fix_pos),byrow=TRUE) +
        matrix(rnorm(length(fix_pos)*length(v_sd_pos),mean=0,sd=rep(v_sd_pos,length(fix_pos))),nrow=length(fix_pos),byrow=TRUE)
      Y_pos = Y_pos + fix_pos
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_pos~.,data=dat_all_df)),get_p)))

      ## calculate p-value for negative features
      #print("calculate p-value for negative features")
      fix_neg = c(dat_all %*% beta_all_neg)
      v_sd_neg = c(rep(quartile_sd[1],floor(length(v_mean_neg)/3)),
                   rep(quartile_sd[2],floor(length(v_mean_neg)/3)),
                   rep(quartile_sd[3],length(v_mean_neg) - floor(length(v_mean_neg)/3)*2))
      Y_neg = matrix(rep(v_mean_neg,length(fix_neg)),nrow=length(fix_neg),byrow=TRUE) +
        matrix(rnorm(length(fix_neg)*length(v_sd_neg),mean=0,sd=rep(v_sd_neg,length(fix_neg))),nrow=length(fix_neg),byrow=TRUE)
      Y_neg = Y_neg + fix_neg
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_neg~.,data=dat_all_df)),get_p)))


      qvalue_matrix[rep.i,] = adjp(pvalue_v)[,2]
    }

    is.sig = qvalue_matrix < q_cutoff
    feature_power = apply(is.sig,2,mean)
    power_current_1 = mean(feature_power[1:(floor(length(v_mean_pos)/3))])
    power_current_2 = mean(feature_power[(floor(length(v_mean_pos)/3)+1):(floor(length(v_mean_pos)/3)*2)])
    power_current_3 = mean(feature_power[((floor(length(v_mean_pos)/3)*2)+1):(length(v_mean_pos))])
    fpr_current_1 = mean(feature_power[(n_feature_pos+1):(n_feature_pos+floor(length(v_mean_neg)/3))])
    fpr_current_2 = mean(feature_power[(n_feature_pos+floor(length(v_mean_neg)/3)+1):(n_feature_pos+floor(length(v_mean_neg)/3)*2)])
    fpr_current_3 = mean(feature_power[(n_feature_pos+floor(length(v_mean_neg)/3)*2+1):(n_feature_pos+n_feature_neg)])

    power_current = c(power_current_1,power_current_2,power_current_3)
    fpr_current = c(fpr_current_1,fpr_current_2,fpr_current_3)
    power_current_min = min(power_current)

    power_out = rbind(power_out,power_current)
    fpr_out = rbind(fpr_out,fpr_current)
    beta_out = c(beta_out,beta_interest)

    beta_interest = beta_interest + increment
    #print(c(beta_interest,power_current,fpr_current))
    # calculate time for each beta
    #end.time <- Sys.time()
    #print(end.time - start.time)
  }

  ## calculate min difference
  diff_out_common = c()

  for(j in 1:100){
    sample(v_mean_transform,N_feature,replace=TRUE)

    # common
    Y_diff_common = rnorm(dim(dat_all)[1],mean=sample(v_mean_transform,dim(dat_all)[1],replace = TRUE),sd=sample(quartile_sd,dim(dat_all)[1],replace = TRUE))

    dat_all_neg = dat_all
    dat_all_neg[,1] = dat_all_neg[,1] - 1

    diff_out = c()
    for(i in out_common$beta){
      beta_all_pos[1] = i
      fix_neg = dat_all_neg[,1]*i
      fix_pos = dat_all[,1]*i
      Y_diff_pos = Y_diff_common + fix_pos
      Y_diff_neg = Y_diff_common + fix_neg
      diff_out = c(diff_out,mean(sin(Y_diff_pos^2*sign(Y_diff_pos))-mean(sin(Y_diff_neg^2*sign(Y_diff_neg)))))
    }
    diff_out_common = rbind(diff_out_common,diff_out)

  }

  out = cbind(power_out,fpr_out,beta_out,diff_out)
  colnames(out) = c("power_1","power_2","power_3","fpr_1","fpr_2","fpr_3","beta","diff_out")
  rownames(out) = NULL
  out = as.data.frame(out)
  return(out)
}


#' Plot power and fpr curve
#'
#' \code{plot_power} plots the power and fpr curve for common and rare species or pathways.
#'
#' @param power.obj A data frame generated by \code{calc_power} or \code{calc_power_common_rare}
#' @param target_power The least power the study wants to reach. A number in (0,1). Default is 0.9.
#' @param target_fpr The highest fpr of the study. A number in (0,1). Default is 0.05.
#' @param title title for the plot. Default is "".
#' @return None.
#' @author Liu Cao
#' @details
#' ehehehehheheheheheheheeh
#' @import ggplot2 reshape2 gridExtra
#' @export

plot_power <-function(power.obj,target_power=0.9,target_fpr=0.05,title=""){
  #data preparation
  df_power = power.obj[,c("power_1","power_2","power_3","diff_out")]
  colnames(df_power) = c("25%","50%","75%","diff_out")
  df_power = melt(df_power,id="diff_out")
  df_fpr = power.obj[,c("fpr_1","fpr_2","fpr_3","diff_out")]
  colnames(df_fpr) = c("25%","50%","75%","diff_out")
  df_fpr = melt(df_fpr,id="diff_out")

  #plot
  p_power = ggplot(df_power,aes(x=diff_out,y=value)) + geom_line(aes(group=variable,colour=variable)) + geom_hline(yintercept=target_power, linetype="dashed", color = "red") +
    theme_bw() + labs(title=paste(title,"Power Curve"),y="power",x="minimum difference") + scale_y_continuous(breaks=c(0,0.25,0.5,0.75,target_power,1),limits = c(0,1))

  print(p_power)
  p_fpr = ggplot(df_fpr,aes(x=diff_out,y=value)) + geom_line(aes(group=variable,colour=variable)) + geom_hline(yintercept=target_fpr, linetype="dashed", color = "red") +
    theme_bw() + labs(title=paste(title,"FPR Curve"),y="fpr",x="minimum difference") + scale_y_continuous(breaks=c(0,0.02,0.04,0.06,0.07,target_fpr),limits = c(0,0.07))

  print(p_fpr)
}

#' Minimal difference for target power
#'
#' \code{result_power} Calculate minimal difference for common and rare species or pathways at target power.
#'
#' @param power.obj A data frame generated by \code{calc_power}
#' @param target_power The least power the study wants to reach. A number in (0,1). Default is 0.9.
#' @param target_fpr The highest fpr of the study. A number in (0,1). Default is 0.05.
#' @param title title for the plot. Default is "".
#' @return None.
#' @author Liu Cao
#' @details
#' ehehehehheheheheheheheeh
#' @export

result_power <-function(power.obj,target_power=0.9,target_fpr=0.05){
  #data preparation
  out_1 = target_value(power.obj$power_1,power.obj$fpr_1,power.obj$beta,power.obj$diff,target_power = target_power)
  out_2 = target_value(power.obj$power_2,power.obj$fpr_2,power.obj$beta,power.obj$diff,target_power = target_power)
  out_3 = target_value(power.obj$power_3,power.obj$fpr_3,power.obj$beta,power.obj$diff,target_power = target_power)

  out_diff = cbind(out_1,out_2,out_3)
  rownames(out_diff) = c("diff","FPR")
  colnames(out_diff) = c("25%","50%","75%")
  return(out_diff)
}


#' Power calculation for a given study design
#'
#' \code{calc_power_common_rare} calculates the power for common and rare species and pathways. The difference between this function and \code{calc_power} is that \code{calc_power} does not consider common and rare features.
#'
#' @param N_sample Number of subjects in the study. Integer.
#' @param N_covariate Number of covariates in the study. Integer.
#' @param N_feature Number of expected species or pathways. Integer.
#' @param N_repeat Number of repeat measurements for each data set. Integer.
#' @param conti_prop The proportion of continuous covariates. A number in [0,1]. Default is 0.5.
#' @param pos_prop The proportion of species or pathways that are truely correlated with the covariate we are interested in. A number in [0,1]. Default is 0.2.
#' @param contin The type of covariate we are interested. The type includeds "continous", "discrete" and "inter" (interaction). Strint. Default is "continuous".
#' @param increment A small number to adjust the smoothness of the power curve. The smaller, the more smooth. Default is 0.005.
#' @param q_cutoff The threshold for q value (p value after BH adjustment). A number in (0,1). Default is 0.05.
#' @param quartile_common Standard deviation quartile common
#' @param quartile_rare Standard deviation quartile rare
#' @param v_mean_transform Transformed mean distribution
#'
#' @return A data frame with power and beta.
#' @author Liu Cao
#' @details
#' Returning a data frame containing power, FPR, Beta and mininum difference
#' @export

calc_power_common_rare <- function(N_sample=200,
                                   N_covariate=8,
                                   N_feature=300,
                                   N_repeat = 1,
                                   conti_prop=0.5,
                                   pos_prop=0.2,
                                   contin="continuous",
                                   increment=0.005,
                                   q_cutoff = 0.05,
                                   quartile_common,
                                   quartile_rare,
                                   v_mean_transform){
  # return power and corresponding effect size
  #print("return power and corresponding effect size")
  n_continuous = ceiling(conti_prop*N_covariate)
  n_discrete = N_covariate - n_continuous

  n_feature_pos = ceiling(N_feature*pos_prop)
  n_feature_neg = N_feature - n_feature_pos

  ########## common #############
  power_current_min = 0
  power_current = c()
  fpr_current = c()
  beta_interest = 0
  power_out = c()
  fpr_out = c()
  beta_out = c()
  diff_out = c()
  num_rep = 100

  #print("iteration")
  while(power_current_min < 0.99){
    qvalue_matrix = matrix(rep(0,num_rep*N_feature), nrow=num_rep)
    start.time <- Sys.time()
    for( rep.i in 1:num_rep ){
      ## initialize p-value vector for each replicate
      pvalue_v = c()
      ## simulate covariates
      dat_continuous = matrix(rnorm(n_continuous*N_sample,mean=0,sd=1),nrow=N_sample)
      dat_discrete = matrix(rbinom(n_discrete*N_sample,size=1,prob=0.5),nrow=N_sample)
      if(N_repeat > 1){
        dat_time = rep(1,N_sample)
        dat_continuous_single = dat_continuous
        dat_discrete_single = dat_discrete
        for(i in 2:N_repeat){
          dat_time = c(dat_time,rep(i,N_sample))
          dat_continuous = rbind(dat_continuous,dat_continuous_single)
          dat_discrete = rbind(dat_discrete,dat_discrete_single)
        }
      }

      ## assign coefficient for positive features and negative features separately
      beta_continuous_neg = rep(0.5,n_continuous)
      beta_discrete_neg = rep(0.5,n_discrete)
      beta_continuous_pos = rep(0.5,n_continuous)
      beta_discrete_pos = rep(0.5,n_discrete)

      ## assign coefficient for covariate of interest
      if(contin == "continuous"){
        beta_continuous_pos[1] = beta_interest
        beta_continuous_neg[1] = 0
        if(N_repeat == 1){
          beta_all_pos = c(beta_continuous_pos,beta_discrete_pos)
          beta_all_neg = c(beta_continuous_neg,beta_discrete_neg)
          dat_all = cbind(dat_continuous,dat_discrete)
        }
        else{
          beta_all_pos = c(beta_continuous_pos,beta_discrete_pos,0.5)
          beta_all_neg = c(beta_continuous_neg,beta_discrete_neg,0.5)
          dat_all = cbind(dat_continuous,dat_discrete,dat_time)
        }
      }
      if(contin == "discrete"){
        beta_discrete_pos[1] = beta_interest
        beta_discrete_neg[1] = 0
        if(N_repeat == 1){
          beta_all_pos = c(beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(beta_discrete_neg,beta_continuous_neg)
          dat_all = cbind(dat_discrete,dat_continuous)
        }
        else{
          beta_all_pos = c(beta_discrete_pos,beta_continuous_pos,0.5)
          beta_all_neg = c(beta_discrete_neg,beta_continuous_neg,0.5)
          dat_all = cbind(dat_discrete,dat_continuous,dat_time)
        }
      }
      if(contin == "inter"){
        if(N_repeat == 1){
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg)
          inter = dat_discrete[,1]*dat_discrete[,2]
          dat_all = cbind(inter,dat_discrete,dat_continuous)
        }
        else{
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos,0.5)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg,0.5)
          inter = dat_discrete[,1]*dat_discrete[,2]
          dat_all = cbind(inter,dat_discrete,dat_continuous,dat_time)
        }
      }
      if(contin == "time"){
        if(N_repeat == 1){
          stop("No repeat measure for the samples.\n")
        }
        else{
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg)
          dat_all = cbind(dat_time,dat_discrete,dat_continuous)
        }
      }

      ## sample mean

      v_mean = sample(v_mean_transform,N_feature,replace=TRUE)
      v_mean_common = sort(v_mean,decreasing=TRUE)[1:(N_feature/2)]
      index_common_pos = sample(1:length(v_mean_common),floor(pos_prop*length(v_mean_common)),replace=FALSE)
      v_mean_common_pos = v_mean_common[1:length(v_mean_common) %in% index_common_pos]
      v_mean_common_neg = v_mean_common[!(1:length(v_mean_common) %in% index_common_pos)]

      v_mean_rare = sort(v_mean,decreasing=TRUE)[ceiling(N_feature/2+1):N_feature]
      index_rare_pos = sample(1:length(v_mean_rare),floor(pos_prop*length(v_mean_rare)),replace=FALSE)
      v_mean_rare_pos = v_mean_rare[1:length(v_mean_rare) %in% index_rare_pos]
      v_mean_rare_neg = v_mean_rare[!(1:length(v_mean_rare) %in% index_rare_pos)]

      ## calculate p-value for positive features
      dat_all_df = data.frame(dat_all)
      fix_pos = c(dat_all %*% beta_all_pos)
      v_sd_common_pos = c(rep(quartile_common[1],floor(length(v_mean_common_pos)/3)),rep(quartile_common[2],floor(length(v_mean_common_pos)/3)),rep(quartile_common[3],length(v_mean_common_pos)-floor(length(v_mean_common_pos)/3)*2))
      v_sd_rare_pos = c(rep(quartile_rare[1],floor(length(v_mean_rare_pos)/3)),rep(quartile_rare[2],floor(length(v_mean_rare_pos)/3)),rep(quartile_rare[3],length(v_mean_rare_pos)-floor(length(v_mean_rare_pos)/3)*2))

      Y_pos = matrix(rep(v_mean_common_pos,length(fix_pos)),nrow=length(fix_pos),byrow=TRUE) +
        matrix(rnorm(length(fix_pos)*length(v_sd_common_pos),mean=0,sd=rep(v_sd_common_pos,length(fix_pos))),nrow=length(fix_pos),byrow=TRUE)
      Y_pos = Y_pos + fix_pos
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_pos~.,data=dat_all_df)),get_p)))

      Y_pos = matrix(rep(v_mean_rare_pos,length(fix_pos)),nrow=length(fix_pos),byrow=TRUE) +
        matrix(rnorm(length(fix_pos)*length(v_sd_rare_pos),mean=0,sd=rep(v_sd_rare_pos,length(fix_pos))),nrow=length(fix_pos),byrow=TRUE)
      Y_pos = Y_pos + fix_pos
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_pos~.,data=dat_all_df)),get_p)))

      ## calculate p-value for negative features
      fix_neg = c(dat_all %*% beta_all_neg)
      v_sd_common_neg = c(rep(quartile_common[1],floor(length(v_mean_common_neg)/3)),rep(quartile_common[2],floor(length(v_mean_common_neg)/3)),rep(quartile_common[3],length(v_mean_common_neg)-floor(length(v_mean_common_neg)/3)*2))
      v_sd_rare_neg = c(rep(quartile_rare[1],floor(length(v_mean_rare_neg)/3)),rep(quartile_rare[2],floor(length(v_mean_rare_neg)/3)),rep(quartile_rare[3],length(v_mean_rare_neg)-floor(length(v_mean_rare_neg)/3)*2))

      Y_neg = matrix(rep(v_mean_common_neg,length(fix_neg)),nrow=length(fix_neg),byrow=TRUE) +
        matrix(rnorm(length(fix_neg)*length(v_sd_common_neg),mean=0,sd=rep(v_sd_common_neg,length(fix_neg))),nrow=length(fix_neg),byrow=TRUE)
      Y_neg = Y_neg + fix_neg
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_neg~.,data=dat_all_df)),get_p)))

      Y_neg = matrix(rep(v_mean_rare_neg,length(fix_neg)),nrow=length(fix_neg),byrow=TRUE) +
        matrix(rnorm(length(fix_neg)*length(v_sd_rare_neg),mean=0,sd=rep(v_sd_rare_neg,length(fix_neg))),nrow=length(fix_neg),byrow=TRUE)
      Y_neg = Y_neg + fix_neg
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_neg~.,data=dat_all_df)),get_p)))

      qvalue_matrix[rep.i,] = adjp(pvalue_v)[,2]
    }

    is.sig = qvalue_matrix < q_cutoff
    feature_power = apply(is.sig,2,mean)
    power_current_common1 = mean(feature_power[1:(floor(length(v_mean_common_pos)/3))])
    power_current_common2 = mean(feature_power[(floor(length(v_mean_common_pos)/3)+1):(floor(length(v_mean_common_pos)/3)*2)])
    power_current_common3 = mean(feature_power[((floor(length(v_mean_common_pos)/3)*2)+1):(length(v_mean_common_pos))])
    fpr_current_common1 = mean(feature_power[(n_feature_pos+1):(n_feature_pos+floor(length(v_mean_common_neg)/3))])
    fpr_current_common2 = mean(feature_power[(n_feature_pos+floor(length(v_mean_common_neg)/3)+1):(n_feature_pos+floor(length(v_mean_common_neg)/3)*2)])
    fpr_current_common3 = mean(feature_power[(n_feature_pos+floor(length(v_mean_common_neg)/3)*2+1):(n_feature_pos+length(v_mean_common_neg))])

    power_current = c(power_current_common1,power_current_common2,power_current_common3)
    fpr_current = c(fpr_current_common1,fpr_current_common2,fpr_current_common3)
    power_current_min = min(power_current)

    power_out = rbind(power_out,power_current)
    fpr_out = rbind(fpr_out,fpr_current)
    beta_out = c(beta_out,beta_interest)

    beta_interest = beta_interest + increment
    #print(c(beta_interest,power_current,fpr_current))
    # calculate time for each beta
    #end.time <- Sys.time()
    #print(end.time - start.time)
  }

  out_common = cbind(power_out,fpr_out,beta_out)
  colnames(out_common) = c("power_1","power_2","power_3",
                           "fpr_1","fpr_2","fpr_3","beta")
  rownames(out_common) = NULL
  out_common = as.data.frame(out_common)




  ############ rare ##############
  power_current_min = 0
  power_current = c()
  fpr_current = c()
  beta_interest = 0
  power_out = c()
  fpr_out = c()
  beta_out = c()
  diff_out = c()
  num_rep = 100

  #print("iteration")
  while(power_current_min < 0.99){
    qvalue_matrix = matrix(rep(0,num_rep*N_feature),nrow=num_rep)
    #start.time <- Sys.time()
    for( rep.i in 1:num_rep ){
      ## initialize p-value vector for each replicate
      pvalue_v = c()
      ## simulate covariates
      dat_continuous = matrix(rnorm(n_continuous*N_sample,mean=0,sd=1),nrow=N_sample)
      dat_discrete = matrix(rbinom(n_discrete*N_sample,size=1,prob=0.5),nrow=N_sample)
      if(N_repeat > 1){
        dat_time = rep(1,N_sample)
        dat_continuous_single = dat_continuous
        dat_discrete_single = dat_discrete
        for(i in 2:N_repeat){
          dat_time = c(dat_time,rep(i,N_sample))
          dat_continuous = rbind(dat_continuous,dat_continuous_single)
          dat_discrete = rbind(dat_discrete,dat_discrete_single)
        }
      }

      ## assign coefficient for positive features and negative features separately
      beta_continuous_neg = rep(0.5,n_continuous)
      beta_discrete_neg = rep(0.5,n_discrete)
      beta_continuous_pos = rep(0.5,n_continuous)
      beta_discrete_pos = rep(0.5,n_discrete)

      ## assign coefficient for covariate of interest
      if(contin == "continuous"){
        beta_continuous_pos[1] = beta_interest
        beta_continuous_neg[1] = 0
        if(N_repeat == 1){
          beta_all_pos = c(beta_continuous_pos,beta_discrete_pos)
          beta_all_neg = c(beta_continuous_neg,beta_discrete_neg)
          dat_all = cbind(dat_continuous,dat_discrete)
        }
        else{
          beta_all_pos = c(beta_continuous_pos,beta_discrete_pos,0.5)
          beta_all_neg = c(beta_continuous_neg,beta_discrete_neg,0.5)
          dat_all = cbind(dat_continuous,dat_discrete,dat_time)
        }
      }
      if(contin == "discrete"){
        beta_discrete_pos[1] = beta_interest
        beta_discrete_neg[1] = 0
        if(N_repeat == 1){
          beta_all_pos = c(beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(beta_discrete_neg,beta_continuous_neg)
          dat_all = cbind(dat_discrete,dat_continuous)
        }
        else{
          beta_all_pos = c(beta_discrete_pos,beta_continuous_pos,0.5)
          beta_all_neg = c(beta_discrete_neg,beta_continuous_neg,0.5)
          dat_all = cbind(dat_discrete,dat_continuous,dat_time)
        }
      }
      if(contin == "inter"){
        if(N_repeat == 1){
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg)
          inter = dat_discrete[,1]*dat_discrete[,2]
          dat_all = cbind(inter,dat_discrete,dat_continuous)
        }
        else{
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos,0.5)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg,0.5)
          inter = dat_discrete[,1]*dat_discrete[,2]
          dat_all = cbind(inter,dat_discrete,dat_continuous,dat_time)
        }
      }
      if(contin == "time"){
        if(N_repeat == 1){
          stop("No repeat measure for the samples.\n")
        }
        else{
          beta_all_pos = c(beta_interest,beta_discrete_pos,beta_continuous_pos)
          beta_all_neg = c(0,beta_discrete_neg,beta_continuous_neg)
          dat_all = cbind(dat_time,dat_discrete,dat_continuous)
        }
      }

      ## sample mean

      v_mean = sample(v_mean_transform,N_feature,replace=TRUE)
      v_mean_common = sort(v_mean)[1:(N_feature/2)]
      index_common_pos = sample(1:length(v_mean_common),floor(pos_prop*length(v_mean_common)),replace=FALSE)
      v_mean_common_pos = v_mean_common[1:length(v_mean_common) %in% index_common_pos]
      v_mean_common_neg = v_mean_common[!(1:length(v_mean_common) %in% index_common_pos)]

      v_mean_rare = sort(v_mean)[ceiling(N_feature/2+1):N_feature]
      index_rare_pos = sample(1:length(v_mean_rare),floor(pos_prop*length(v_mean_rare)),replace=FALSE)
      v_mean_rare_pos = v_mean_rare[1:length(v_mean_rare) %in% index_rare_pos]
      v_mean_rare_neg = v_mean_rare[!(1:length(v_mean_rare) %in% index_rare_pos)]

      ## calculate p-value for positive features
      dat_all_df = data.frame(dat_all)
      fix_pos = c(dat_all %*% beta_all_pos)
      v_sd_common_pos = c(rep(quartile_common[1],floor(length(v_mean_common_pos)/3)),rep(quartile_common[2],floor(length(v_mean_common_pos)/3)),rep(quartile_common[3],length(v_mean_common_pos)-floor(length(v_mean_common_pos)/3)*2))
      v_sd_rare_pos = c(rep(quartile_rare[1],floor(length(v_mean_rare_pos)/3)),rep(quartile_rare[2],floor(length(v_mean_rare_pos)/3)),rep(quartile_rare[3],length(v_mean_rare_pos)-floor(length(v_mean_rare_pos)/3)*2))

      Y_pos = matrix(rep(v_mean_common_pos,length(fix_pos)),nrow=length(fix_pos),byrow=TRUE) +
        matrix(rnorm(length(fix_pos)*length(v_sd_common_pos),mean=0,sd=rep(v_sd_common_pos,length(fix_pos))),nrow=length(fix_pos),byrow=TRUE)
      Y_pos = Y_pos + fix_pos
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_pos~.,data=dat_all_df)),get_p)))

      Y_pos = matrix(rep(v_mean_rare_pos,length(fix_pos)),nrow=length(fix_pos),byrow=TRUE) +
        matrix(rnorm(length(fix_pos)*length(v_sd_rare_pos),mean=0,sd=rep(v_sd_rare_pos,length(fix_pos))),nrow=length(fix_pos),byrow=TRUE)
      Y_pos = Y_pos + fix_pos
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_pos~.,data=dat_all_df)),get_p)))

      ## calculate p-value for negative features
      fix_neg = c(dat_all %*% beta_all_neg)
      v_sd_common_neg = c(rep(quartile_common[1],floor(length(v_mean_common_neg)/3)),rep(quartile_common[2],floor(length(v_mean_common_neg)/3)),rep(quartile_common[3],length(v_mean_common_neg)-floor(length(v_mean_common_neg)/3)*2))
      v_sd_rare_neg = c(rep(quartile_rare[1],floor(length(v_mean_rare_neg)/3)),rep(quartile_rare[2],floor(length(v_mean_rare_neg)/3)),rep(quartile_rare[3],length(v_mean_rare_neg)-floor(length(v_mean_rare_neg)/3)*2))

      Y_neg = matrix(rep(v_mean_common_neg,length(fix_neg)),nrow=length(fix_neg),byrow=TRUE) +
        matrix(rnorm(length(fix_neg)*length(v_sd_common_neg),mean=0,sd=rep(v_sd_common_neg,length(fix_neg))),nrow=length(fix_neg),byrow=TRUE)
      Y_neg = Y_neg + fix_neg
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_neg~.,data=dat_all_df)),get_p)))

      Y_neg = matrix(rep(v_mean_rare_neg,length(fix_neg)),nrow=length(fix_neg),byrow=TRUE) +
        matrix(rnorm(length(fix_neg)*length(v_sd_rare_neg),mean=0,sd=rep(v_sd_rare_neg,length(fix_neg))),nrow=length(fix_neg),byrow=TRUE)
      Y_neg = Y_neg + fix_neg
      pvalue_v = c(pvalue_v,unlist(lapply(summary(lm(Y_neg~.,data=dat_all_df)),get_p)))

      qvalue_matrix[rep.i,] = adjp(pvalue_v)[,2]
    }

    is.sig = qvalue_matrix < q_cutoff
    feature_power = apply(is.sig,2,mean)
    power_current_rare1 = mean(feature_power[(length(v_mean_common_pos)+1):(length(v_mean_common_pos)+floor(length(v_mean_rare_pos)/3))])
    power_current_rare2 = mean(feature_power[((length(v_mean_common_pos)+floor(length(v_mean_rare_pos)/3))+1):((length(v_mean_common_pos)+floor(length(v_mean_rare_pos)/3)*2))])
    power_current_rare3 = mean(feature_power[((length(v_mean_common_pos)+floor(length(v_mean_rare_pos)/3)*2)+1):(n_feature_pos)])
    fpr_current_rare1 = mean(feature_power[(n_feature_pos+length(v_mean_common_neg)+1):(n_feature_pos+length(v_mean_common_neg)+floor(length(v_mean_rare_neg)/3) )])
    fpr_current_rare2 = mean(feature_power[(n_feature_pos+length(v_mean_common_neg)+floor(length(v_mean_rare_neg)/3)+1):(n_feature_pos+length(v_mean_common_neg)+floor(length(v_mean_rare_neg)/3)*2)])
    fpr_current_rare3 = mean(feature_power[(n_feature_pos+length(v_mean_common_neg)+floor(length(v_mean_rare_neg)/3)*2+1):(n_feature_pos+n_feature_neg)])

    power_current = c(power_current_rare1,power_current_rare2,power_current_rare3)
    fpr_current = c(fpr_current_rare1,fpr_current_rare2,fpr_current_rare3)
    power_current_min = min(power_current)

    power_out = rbind(power_out,power_current)
    fpr_out = rbind(fpr_out,fpr_current)
    beta_out = c(beta_out,beta_interest)

    beta_interest = beta_interest + increment/10

    #print(c(beta_interest,power_current,fpr_current))
    # calculate time for each beta
    #end.time <- Sys.time()
    #print(end.time - start.time)
  }

  out_rare = cbind(power_out,fpr_out,beta_out)
  colnames(out_rare) = c("power_1","power_2","power_3",
                         "fpr_1","fpr_2","fpr_3","beta")
  rownames(out_rare) = NULL
  out_rare = as.data.frame(out_rare)


  ####### calculate min difference ######
  diff_out_common = c()
  diff_out_rare = c()

  for(j in 1:100){
  sample(v_mean_transform,N_feature,replace=TRUE)
  v_mean_common = sort(v_mean_transform,decreasing=TRUE)[1:(length(v_mean_transform)/2)]
  v_mean_rare = sort(v_mean_transform,decreasing = FALSE)[1:(length(v_mean_transform)/2)]

  # common
  Y_diff_common = rnorm(dim(dat_all)[1],mean=sample(v_mean_common,dim(dat_all)[1],replace = TRUE),sd=sample(quartile_common,dim(dat_all)[1],replace = TRUE))

  dat_all_neg = dat_all
  dat_all_neg[,1] = dat_all_neg[,1] - 1

  diff_out = c()
  for(i in out_common$beta){
    beta_all_pos[1] = i
    #fix_neg = c(dat_all_neg %*% beta_all_pos)
    fix_neg = dat_all_neg[,1]*i
    #fix_pos = c(dat_all %*% beta_all_pos)
    fix_pos = dat_all[,1]*i
    Y_diff_pos = Y_diff_common + fix_pos
    Y_diff_neg = Y_diff_common + fix_neg
    diff_out = c(diff_out,mean(sin(Y_diff_pos^2*sign(Y_diff_pos))-mean(sin(Y_diff_neg^2*sign(Y_diff_neg)))))
  }
  diff_out_common = rbind(diff_out_common,diff_out)


  # rare
  Y_diff_rare = rnorm(dim(dat_all)[1],mean=sample(v_mean_rare,dim(dat_all)[1],replace = TRUE),sd=sample(quartile_rare,dim(dat_all)[1],replace = TRUE))
  diff_out = c()
  for(i in out_rare$beta){
    beta_all_pos[1] = i
    #fix_neg = c(dat_all_neg %*% beta_all_pos)
    #fix_pos = c(dat_all %*% beta_all_pos)
    fix_neg = dat_all_neg[,1]*i
    fix_pos = dat_all[,1]*i
    Y_diff_pos = Y_diff_rare + fix_pos
    Y_diff_neg = Y_diff_rare + fix_neg
    diff_out = c(diff_out,mean(sin(Y_diff_pos^2*sign(Y_diff_pos)))-mean(sin(Y_diff_neg^2*sign(Y_diff_neg))))
  }
  diff_out_rare = rbind(diff_out_rare,diff_out)

  }

  out_common$diff_out = apply(diff_out_common,2,mean)
  out_rare$diff_out = apply(diff_out_rare,2,mean)
  return(list(out_common=out_common,
              out_rare=out_rare))
}



