#read in packages
sapply(c("forestplot", #display packages
         "survival","SemiCompRisksBin",#survival analysis packages
         "tidyverse","lubridate", #tidy packages
         "rstan","loo","RColorBrewer"),
       function(x){if(!require(package=x,character.only = TRUE))
       {install.packages(pkgs=x);require(package=x,character.only = TRUE)}})

#code will create (possibly large) files here to store log-likelihood contributions
#if estimating ELPD psis-loo is desired (can be set to null if not. DIC and ELPD is-loo still estimated anyways)

# logll_path <- "[path here]"
logll_path <- NULL

####PLOTTING AND FORMATTING HELPERS####

#color blind friendly colors from here: https://davidmathlogic.com/colorblind/#%23648FFF-%23785EF0-%23DC267F-%23FE6100-%23FFB000
cb_blue <- "#648FFF"; cb_red <- "#DC267F"; cb_purple <- "#785EF0"; cb_orange <- "#FE6100"; cb_grey <- "#CACACA"

#colors used in plots
two_color_cb <- c(cb_blue,cb_red)

three_color <- c("dodgerblue","firebrick3","purple3")
three_color_cb <- c(cb_blue,cb_red,cb_purple)
four_color <- c("lightgray","firebrick3","purple3","dodgerblue")
four_color_cb <- c(cb_grey,cb_red,cb_purple,cb_blue)
four_color_forest <- c("dodgerblue","firebrick3","purple3","magenta")
four_color_forest_cb <- c(cb_blue,cb_red,cb_purple,cb_orange)
five_color <- c("lightgray","firebrick3","magenta","purple3","dodgerblue")
five_color_cb <- c(cb_grey,cb_red,cb_orange,cb_purple,cb_blue)

# RColorBrewer::display.brewer.all(n=4,colorblindFriendly = TRUE)
# color-blind friendly categorical colors
three_color_qual <- RColorBrewer::brewer.pal(n=3,name="Set2")
four_color_qual <- RColorBrewer::brewer.pal(n=4,name="Dark2")
five_color_qual <- RColorBrewer::brewer.pal(n=5,name="Dark2")


#extract matrix with estimate information
get_bayes_ests <- function(x,warmup=0,digits_summary=5){
  invisible(
    capture.output(
      out_mat <- as.matrix(
        rstan::monitor(x$samples,
                       warmup = warmup,digits_summary = digits_summary))[,c("50%","2.5%","97.5%","mean","sd","Rhat","Bulk_ESS")]
    )
  )
  out_mat
}

####function for generating labeled forest plot from ID model fit ("standard")
gen_forest_plot <- function(complete_coefs,type="logit",
                            bounds=c(0.2, 120),
                            ticks=log(c(0.2,0.5,1,1.5,3,5,10,20,50)),
                            title=""){
  tabletext <- rownames(complete_coefs)
  forestplot(tabletext,
             mean = cbind(complete_coefs[,1],complete_coefs[,4],complete_coefs[,7],
                          if(type=="logit") complete_coefs[,10]),
             lower = cbind(complete_coefs[,2],complete_coefs[,5],complete_coefs[,8],
                           if(type=="logit") complete_coefs[,11]),
             upper = cbind(complete_coefs[,3],complete_coefs[,6],complete_coefs[,9],
                           if(type=="logit") complete_coefs[,12]),
             title = title, xlog=TRUE,
             clip =bounds,
             fn.ci_norm = c(function(...) fpDrawPointCI(pch=15,...),
                            function(...) fpDrawPointCI(pch=16,...),
                            function(...) fpDrawPointCI(pch=17,...),
                            if(type=="logit") function(...) fpDrawPointCI(pch=18,...)),
             boxsize = .2, xticks = ticks,
             zero=1,
             xlab=if(type=="logit") "Hazard/Odds Ratio (Estimate, 95% CrI)"
               else "Hazard Ratio (Estimate, 95% CrI)",
             legend=if(type=="logit") c("Preeclampsia\nDiagnosis",
                                        "Delivery without\nPreeclampsia",
                                        "Delivery after\nPreeclampsia",
                                        "Immediate\nDelivery")
               else c("PE Onset", "Delivery without PE", "Delivery after PE"),
             col= if(type=="logit") fpColors(box=four_color_forest_cb,
                                             lines=four_color_forest_cb)
               else fpColors(box=three_color_cb,
                             lines=three_color_cb),
             lty.ci = c(1, 1,1, if(type=="logit") 1),
             lwd.ci = c(1.75,1.75,1.75,
                        if(type=="logit") 1.75),
             txt_gp = fpTxtGp(xlab=gpar(cex=1.2),
                              ticks=gpar(cex=1.2)),
             hrzl_lines=rep(list(gpar(lwd=2, lineend="butt", col="lightgray")),
                            length(tabletext)+1)
  )
}

get_heterogeneity_plot <- function(pred_panel,t_cutoff_vec, title=NULL,
                                   subtitle = NULL, headsub = NA, nrows = 1,
                                   type="logit"){
  # browser()
  plot_frame_panel <- NULL

  nt_rows <- if(type=="logit") 1:3 else 1:2
  #overall risk of nonterminal event
  nt_rank <- rank(rowSums(t(pred_panel[paste0("t",max(t_cutoff_vec)),nt_rows,])), ties.method = "first")
  for(time_ind in 1:length(t_cutoff_vec)){
    time <- t_cutoff_vec[time_ind]
    plot_frame <- as.data.frame(t(pred_panel[paste0("t",t_cutoff_vec[time_ind]),,])) %>%
      mutate(Time=time, subject_id = row_number(), subject_ordered = nt_rank) %>%
      pivot_longer(cols=starts_with("p_"),names_to = "Outcome", values_to = "Probability")
    if(type=="logit"){
      plot_frame_factor <- plot_frame %>%
        mutate(Outcome = recode_factor(Outcome,
                                       "p_neither"="Pregnant without\nPreeclampsia",
                                       "p_tonly"="Delivered without\nPreeclampsia",
                                       "p_both_inst"="Delivered Immediately\nwith Preeclampsia",
                                       "p_both_noinst"="Delivered Non-Immediately\nwith Preeclampsia",
                                       "p_ntonly"="Pregnant with\nPreeclampsia"))
    } else{
      plot_frame_factor <- plot_frame %>%
      mutate(Outcome = recode_factor(Outcome,
                                     "p_neither"="Pregnant without\nPreeclampsia",
                                     "p_tonly"="Delivered without\nPreeclampsia",
                                     "p_both"="Delivered with\nPreeclampsia",
                                     "p_ntonly"="Pregnant with\nPreeclampsia"))
    }
    plot_frame_panel <- bind_rows(plot_frame_panel, plot_frame_factor)
  }

  color_temp <- if(type=="logit") five_color_cb else four_color_cb
  gg <- ggplot(plot_frame_panel %>%
                 filter(subject_ordered > length(nt_rank) - headsub),
               aes(x=subject_ordered, y=Probability)) +
    geom_bar(aes(colour=Outcome, fill=Outcome),#width=1,
             stat="identity") +
    scale_fill_manual(values = color_temp) + scale_color_manual(values = color_temp) +
    theme_classic() + facet_wrap(~as.factor(paste0(Time," Weeks")),nrow = nrows) +
    theme(axis.text.x = element_blank(),axis.ticks.x = element_blank(),
          legend.position = "bottom",legend.title = element_blank()) +
    labs(title = title,
         subtitle = subtitle,
         x=paste0("Individuals Ordered by Overall Preeclampsia Risk by Week ",
                 max(t_cutoff_vec)*10 + 20," (Top ",headsub,")"),
         y = "Probability of Having Experienced Outcome") #  +
  # theme(axis.text=element_text(size=12),axis.title=element_text(size=16),legend.title=element_text(size=14),legend.text=element_text(size=12))
  gg
}

get_sample_plot <- function(pred_selected, t_cutoff_vec, i, type="logit",
                            plot_type = "stacked"){
  # browser()
  #create long-format data to plot
  plot_frame <- cbind(Time=t_cutoff_vec,
                      as.data.frame(pred_selected)) %>%
    pivot_longer(cols=starts_with("p_"), names_to = "Outcome", values_to = "Probability")

  if(type=="logit"){
    plot_frame_factor <- plot_frame %>%
      mutate(Outcome = recode_factor(Outcome,
                                     "p_neither"="Pregnant without\nPreeclampsia",
                                     "p_tonly"="Delivered without\nPreeclampsia",
                                     "p_both_inst"="Delivered Immediately\nwith Preeclampsia",
                                     "p_both_noinst"="Delivered Non-Immediately\nwith Preeclampsia",
                                     "p_ntonly"="Pregnant with\nPreeclampsia"))
  } else{
    plot_frame_factor <- plot_frame %>%
      mutate(Outcome = recode_factor(Outcome,
                                     "p_neither"="Pregnant without\nPreeclampsia",
                                     "p_tonly"="Delivered without\nPreeclampsia",
                                     "p_both"="Delivered with\nPreeclampsia",
                                     "p_ntonly"="Pregnant with\nPreeclampsia"))
  }

  color_temp <- if(type=="logit") five_color_cb else four_color_cb

  gg <- ggplot(plot_frame_factor, aes(x=Time*10 + 20, y=Probability))
  if(plot_type == "line"){
    gg <- gg + geom_line(aes(colour=Outcome)) #+
#      guides(colour = guide_legend(reverse = T))
  } else{
   gg <- gg + geom_area(aes(colour=Outcome, fill=Outcome)) +
     scale_fill_manual(values = color_temp) #+
#     guides(colour = guide_legend(reverse = T),fill = guide_legend(reverse = T))
  }
  gg <- gg +
    scale_color_manual(values = color_temp) +
    scale_x_continuous(breaks = seq(from=25,to=40,by=3),limits = c(25,40)) +
    # scale_y_continuous(labels = scales::percent) +
    labs(tag = LETTERS[i]) +
    theme_classic() + theme(text=element_text(size=14))
}

#https://stackoverflow.com/questions/39011020/ggplot-align-plots-together-and-add-common-labels-and-legend
#extract legend
get_legend <- function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

##********************##
####DATA PREPARATION####
##********************##

N <- N_test <- 3000

#cutoff for number of days considered "immediate"
true_cutoff <- 2/7/10


true_base <- c(kappa1=exp(-5.71 + 1.62),alpha1=exp(1.59),
               kappa2=exp(-9.54 + 0.47),alpha2=exp(2.62),
               kappa3=exp(1.52 - 0.23 + 0.75),
               alpha3=exp(0.087),
               betaD1 = -0.49 - 0.33)

true_theta <- 0.7
true_beta1 <- c(0.5,0.25)
true_beta2 <- c(0.1,0)
true_beta3 <- c(0.2,0.1)
true_betaD <- c(-0.6,0.1)

true_tv_knots <- c(1.2,1.4,1.7)
true_beta3tv <- c(0.75,1.5,3)
true_betaDtv <- c(0.1,1.5,3)


set.seed(123)
xtemp <- cbind(x1=rbinom(n=N,size=1,prob=0.25),
               x2=rnorm(n=N,mean=0,sd=1))
simData <- cbind(simID_logit2(
  x1 = xtemp, x2 = xtemp, x3 = xtemp, xD = cbind(1,xtemp),
  beta1.true = true_beta1, beta2.true = true_beta2,
  beta3.true = true_beta3,

  # beta3tv.true = NULL,
  h3tv_degree = 0, h3tv_knots = c(0,true_tv_knots,Inf), beta3tv.true = true_beta3tv,

  alpha1.true = true_base[2], alpha2.true = true_base[4],
  alpha3.true = true_base[6], kappa1.true = true_base[1],
  kappa2.true = true_base[3], kappa3.true = true_base[5],
  theta.true = true_theta, beta2frail.true = 1, beta3frail.true = 1,
  anyD = TRUE,
  betaD.true = c(true_base[7],true_betaD), betaDfrail.true = 1,

  # betaDtv.true = NULL,
  betaDtv.true = true_betaDtv, Dtv_degree = 0, Dtv_knots = c(0,true_tv_knots,Inf),

  D_eps = true_cutoff, D_dist = "uniform",
  frailty_type ="lognormal", cens = c(0,0)),xtemp)

xtemp_h3 <- cbind(xtemp,as.matrix(simData[,grep(pattern="h3tv",x = colnames(simData))]))





set.seed(1234)
xtemp_test <- cbind(x1=rbinom(n=N_test,size=1,prob=0.25),
                    x2=rnorm(n=N_test,mean=0,sd=1))
simData_test <- cbind(simID_logit2(
  x1 = xtemp_test, x2 = xtemp_test, x3 = xtemp_test, xD = cbind(1,xtemp_test),
  beta1.true = true_beta1, beta2.true = true_beta2,
  beta3.true = true_beta3,

  # beta3tv.true = NULL,
  h3tv_degree = 0, h3tv_knots = c(0,true_tv_knots,Inf), beta3tv.true = true_beta3tv,

  alpha1.true = true_base[2], alpha2.true = true_base[4],
  alpha3.true = true_base[6], kappa1.true = true_base[1],
  kappa2.true = true_base[3], kappa3.true = true_base[5],
  theta.true = true_theta, beta2frail.true = 1, beta3frail.true = 1,
  anyD = TRUE,
  betaD.true = c(true_base[7],true_betaD), betaDfrail.true = 1,

  # betaDtv.true = NULL,
  betaDtv.true = true_betaDtv, Dtv_degree = 0, Dtv_knots = c(0,true_tv_knots,Inf),

  D_eps = true_cutoff, D_dist = "uniform",
  frailty_type ="lognormal", cens = c(0,0)),xtemp_test)

xtemp_h3_test <- cbind(xtemp_test,as.matrix(simData_test[,grep(pattern="h3tv",x = colnames(simData_test))]))



#in practice, numReps_scr should be set to something larger, like 1000000
#because with only 100000 random scans the chains will not have converged
nchain <- 3
numReps_scr = 100000; burninPerc = 0.3; thin = 10


#model setup is all the same too
## Hyperparameters ##
theta.ab <- c(0.7, 0.7) # prior parameters for 1/theta
WB.ab <- WB.ab1 <- WB.ab2 <- WB.ab3 <- c(0.5, 0.01) # prior parameters for alpha
WB.cd <- WB.cd1 <- WB.cd2 <- WB.cd3 <- c(0.5, 0.05) # prior parameters for kappa
m0 <- rep(0, ncol(xtemp_h3)+1) #prior mean for logit parameters
P0diag <- rep(1, ncol(xtemp_h3)+1) #diagonal of prior precision for logit parameters (aka flat)
hyperParams <- list(theta=theta.ab,
                    WB=list(WB.ab=WB.ab,WB.ab=WB.cd,
                            WB.ab1=WB.ab1, WB.ab2=WB.ab2, WB.ab3=WB.ab3,
                            WB.cd1=WB.cd1, WB.cd2=WB.cd2, WB.cd3=WB.cd3),
                    logit=list(logit.m0=m0,
                               logit.P0diag=P0diag))
nlogLHi_save <- N #number of loglikelihood contributions to save (should be all of them)
nGam_save <- 0 #number of individuals' frailties to save (should be all of them)
## Tuning parameters for specific updates
mhProp_theta_var  <- 0.1
mhProp_alpha_var <- 0.001
mhProp_alphag_var <- c(0.001, 0.001, 0.001)
mcmc.WB_scr  <- list(run=list(numReps=numReps_scr,
                              thin=thin, burninPerc=burninPerc),
                     storage=list(nGam_save=nGam_save,
                                  nlogLHi_save=nlogLHi_save),
                     tuning=list(mhProp_theta_var=mhProp_theta_var,
                                 mhProp_alpha_var=mhProp_alpha_var,
                                 mhProp_alphag_var=mhProp_alphag_var))
myModel_scr <- c("semi-Markov", "Weibull")

#form "true" vector
true_vec <- c(true_base[1:6],true_theta,
              true_beta1,true_beta2,true_beta3,true_beta3tv,
              true_base[7],true_betaD,true_betaDtv)

set.seed(12345)
tempfit <-
  Bayes_SCR_logit_eps(data = simData,
                      y1 = simData$y1,
                      y_sm = simData$y_sm,
                      delta1 = simData$delta1,
                      delta1D = simData$deltaD,
                      delta2 = simData$delta2,
                      eps = true_cutoff,
                      vars1 = c("x1","x2"), vars2 =  c("x1","x2"),
                      vars3 = c("x1","x2",colnames(simData)[grep(pattern="h3tv",x = colnames(simData))]),
                      varsD = c("x1","x2",colnames(simData)[grep(pattern="h3tv",x = colnames(simData))]),
                      hyperParams = hyperParams, mcmcParams = mcmc.WB_scr,
                      n_chains = nchain, frailty= TRUE,
                      frail_path = NULL, logLHi_path = logll_path)

#outputs
get_bayes_ests(tempfit)

temp_para <- apply(tempfit$samples,MARGIN = 3,FUN = median)
temp_pred <- predict(tempfit)

####MODEL COMPARISON####

if(!is.null(logll_path)){
  #computing ELPD psis-loo may require slight further thinning if memory constrained
  thin_loo <- 3
  thin_ind <- (1:(numReps_scr * (1-burninPerc) / thin))[-seq(from=1,
                                                             to=numReps_scr * (1-burninPerc) / thin,
                                                             by=thin_loo)]

  loglik_array <- array(dim = c(numReps_scr * (1-burninPerc) / thin * (thin_loo-1)/thin_loo,
                                nchain,
                                N))
  for(i in 1:nchain){
    # print(i); print(Sys.time())
    loglik_array[,i,] <- t(readRDS(paste0(logll_path,"/logLHi_chain",i,".RDS")))[thin_ind,]
    # file.remove(paste0(logll_path,"/logLHi_chain",i,".RDS"))
  }

  rel_n_eff <- relative_eff(exp(loglik_array), cores=1)
  temp_loo <- loo(loglik_array, r_eff = rel_n_eff, cores = 1)

  loglik_array <- array(dim = c(numReps_scr * (1-burninPerc) / thin * (thin_loo-1)/thin_loo,
                                nchain,
                                N))
  for(i in 1:nchain){
    # print(i); print(Sys.time())
    loglik_array[,i,] <- t(readRDS(paste0(logll_path,"/logLHi_marg_chain",i,".RDS")))[thin_ind,]
    # file.remove(paste0(logll_path,"/logLHi_marg_chain",i,".RDS"))
  }

  rel_n_eff <- relative_eff(exp(loglik_array), cores=1)
  temp_loo_marg <- loo(loglik_array, r_eff = rel_n_eff, cores = 1)
}

#ELPD metrics (larger is better)
c(if(!is.null(logll_path)) c(ELPDpsisloo=temp_loo$estimates[1],
                             ELPDpsisloo_marg=temp_loo_marg$estimates[1]),
  ELPDdic=tempfit$diagnostics$DIC/-2,ELPDdic_marg=tempfit$diagnostics$DIC_marg/-2,
  ELPDisloo=tempfit$diagnostics$LPML,ELPDisloo_marg=tempfit$diagnostics$LPML_marg)

####FOREST PLOT####

bounds_temp <- c(0.25, 120)
ticks_temp <- log(c(0.25,0.5,1,1.5,3,5,7.5,20,50))
  temp_coef_mat <- exp(summary(tempfit)$coef)
gen_forest_plot(complete_coefs = temp_coef_mat,
                type="logit",
                bounds = bounds_temp,
                ticks = ticks_temp,
                title = "")

#### BASELINE PLOTS ####
t_seq <- seq(from=0.01, to=2.2, by=0.01)
t_seq_h3 <- c(0.00005,0.0005,seq(from=0.005, to=0.5, by=0.005))

#first make matrix with one row per baseline stratified by piecewise t1 effect in h3
tv_varnames <- colnames(simData)[grep(pattern="h3tv",x = colnames(simData))]
t1cat_pred_h3_mat <-
  cbind(
    matrix(data = 0, ncol = ncol(xtemp), nrow = length(tv_varnames) + 1),
    rbind(0,diag(length(tv_varnames)))
  )
colnames(t1cat_pred_h3_mat) <- colnames(xtemp_h3)

#now make
temp_fac <- as.factor(cut(t_seq,
                          breaks=c(0,true_tv_knots,Inf),
                          right = FALSE,include.lowest = TRUE,labels=FALSE))
t1cat_pred_logit_mat <- cbind(
  matrix(data = 0, ncol = ncol(xtemp), nrow = length(t_seq)),
  model.matrix(~ 0 + temp_fac)[,-1]
)

par(mfrow=c(2,4))
for(ty in c("h","S")){
  for(i in 1:3){
    #skip plotting h3 for t1 categorical, bc we plot it below
    if(i %in% 1:2){
      plot_mat <- predict(tempfit,
                          tseq=if(i==3) t_seq_h3 else t_seq)
      xlab_temp <- switch(i,"Gestational Age at Admission\nwith Preeclampsia (Weeks)",
                          "Gestational Age at Delivery\nwithout Preeclampsia (Weeks)",
                          "Weeks to Delivery from\nAdmission with Preeclampsia")
      matplot(x=t_seq,
              y=plot_mat[[paste0(ty,".",i)]][,-1],
              type="l", lty=c(1,3,3), lwd=c(2,1,1),
              col = three_color_cb[i],
              ylim = if(ty=="S") c(0,1) else NULL,
              ylab=paste0(
                if(i==3) "Conditional " else "Cause-Specific " ,
                if(ty=="h") "Hazard Function" else "Survivor Function"),
              xlab=xlab_temp,
              add = FALSE )
    } else{
      ##excluding the 95% credible intervals
      plot_mat <- sapply(1:NROW(t1cat_pred_h3_mat), function(j){
        pred <- predict(tempfit,
                        tseq=t_seq_h3, x3new = t1cat_pred_h3_mat[j,])
        pred[[paste0(ty,".3")]][,paste0(ty,".3")] })
      #now, plot the curves for every category
      matplot(x=t_seq_h3,
              y=plot_mat, type="l", col = five_color_qual,
              lty = 1, lwd=2,
              ylim = if(ty=="S") c(0,1) else NULL,
              ylab=paste0(
                if(i==3) "Conditional " else "Cause-Specific " ,
                if(ty=="h") "Hazard Function" else "Survivor Function"),
              xlab="Weeks to Delivery from\nTwo Days after Admission with PE", add = FALSE)
    }
    if(ty=="S" & i==3){
      legend(x="topright", cex=0.8, fill = c(four_color_qual),
             legend = c("ref.",tv_varnames))
    }
  }

  #put a blank plot in the first row of hazards, to leave a buffer
  if(ty=="h"){ plot.new(); next}

  aug_instprob_t1cat_bayes_plotmat <- t(sapply(1:length(t_seq), function(j){
    predict_logit(tempfit,
                  xnew = t1cat_pred_logit_mat[j,])
  }))
  matplot(t_seq, aug_instprob_t1cat_bayes_plotmat,
          type="l",lty=c(1,3,3),col=four_color_forest_cb[4],
          ylab="Probability of Delivery Within Two Days",
          xlab="Gestational Age at Admission\nwith Preeclampsia (Weeks)")
}
par(mfrow=c(1,1))


#### PREDICTION HETEROGENEITY PLOTS ####
#note, for now we're simplifying by
#doing predictions using point estimates based on posterior medians

#same set of covariates for every transition (except categorical), so this same matrix used for all inputs
t_cutoff_pred_vec <- c(1.4,1.7)

#first, generate marginal predictions for every subject across time.
    temp_para <- apply(tempfit$samples,MARGIN = 3,FUN = median)
    pred_m_panel <- calc_risk_logitbayes(para = temp_para, logit_ind = TRUE,
                                         Xmat1 = xtemp, Xmat2 = xtemp,
                                         Xmat3 = xtemp, XmatD = xtemp,
                                         D_eps = true_cutoff,
                                         frailty = TRUE,  #say frailty=TRUE even if it's not doesn't matter in this case
                                         model = "semi-markov",
                                         type = "marginal",
                                         gamma = 1,
                                         hazard = "weibull",
                                         t_cutoff = t_cutoff_pred_vec,
                                         tol = 1e-3,
                                         h3_tv = "piecewise",
                                         h3tv_knots = true_tv_knots,
                                         logit_tv = "piecewise",
                                         logit_tv_knots = true_tv_knots)
get_heterogeneity_plot(pred_panel = pred_m_panel,
                       t_cutoff_vec = t_cutoff_pred_vec,
                       subtitle = NULL, type="logit",
                       headsub = 250, nrows = 2)

#compute brier score in test set
brier_list <- brier_list4 <- kl_list <- kl_list4 <- list(conditional=rep(NA,length(t_cutoff_pred_vec)),
                                                         marginal=rep(NA,length(t_cutoff_pred_vec)))
for(j in 1:length(t_cutoff_pred_vec)){
  for(pred_type in c("conditional","marginal")){
    pred_m_panel_test <- calc_risk_logitbayes(para = temp_para, logit_ind = TRUE,
                                         Xmat1 = xtemp_test, Xmat2 = xtemp_test,
                                         Xmat3 = xtemp_test, XmatD = xtemp_test,
                                         D_eps = true_cutoff,
                                         frailty = TRUE,  #say frailty=TRUE even if it's not doesn't matter in this case
                                         model = "semi-markov",
                                         type = pred_type,
                                         gamma = 1,
                                         hazard = "weibull",
                                         t_cutoff = t_cutoff_pred_vec[j], tol = 1e-3,
                                         h3_tv = "piecewise",
                                         h3tv_knots = true_tv_knots,
                                         logit_tv = "piecewise",
                                         logit_tv_knots = true_tv_knots)

    outcome_mat_test <- get_outcome_mat_logit(y1=simData_test$y1,y2=simData_test$y2,
                                              delta1 = simData_test$delta1,delta2 = simData_test$delta2,
                                              deltaD = simData_test$deltaD,t_cutoff = t_cutoff_pred_vec[j])
    brier_list[[pred_type]][j] <- SemiCompRisksBin:::compute_score(outcome_mat=outcome_mat_test,pred_mat = pred_m_panel_test,
                                                                   ipcw_mat = matrix(nrow = NROW(simData_test),ncol = 1,data = 1),
                                                                   score = "brier")
    kl_list[[pred_type]][j] <- SemiCompRisksBin:::compute_score(outcome_mat=outcome_mat_test,pred_mat = pred_m_panel_test,
                                                                ipcw_mat = matrix(nrow = NROW(simData_test),ncol = 1,data = 1),
                                                                score = "kl")

    #just for fun, look at 'performance' reducing to just 4 categories
    outcome_mat_test4 <- cbind(ntonly=outcome_mat_test[,"ntonly"],
                               both=outcome_mat_test[,"both_noinst"] + outcome_mat_test[,"both_inst"],
                               tonly=outcome_mat_test[,"tonly"],
                               neither=outcome_mat_test[,"neither"])
    pred_m_panel_test4 <- cbind(pred_m_panel_test[,"p_ntonly"],
                                pred_m_panel_test[,"p_both_noinst"]+pred_m_panel_test[,"p_both_inst"],
                                pred_m_panel_test[,"p_tonly"],
                                pred_m_panel_test[,"p_neither"])
    brier_list4[[pred_type]][j] <- SemiCompRisksBin:::compute_score(outcome_mat=outcome_mat_test4,
                                                                    pred_mat = pred_m_panel_test4,
                                                                    ipcw_mat = matrix(nrow = NROW(simData_test),
                                                                                      ncol = 1,
                                                                                      data = 1),
                                                                    score = "brier")
    kl_list4[[pred_type]][j] <- SemiCompRisksBin:::compute_score(outcome_mat=outcome_mat_test4,
                                                                 pred_mat = pred_m_panel_test4,
                                                                 ipcw_mat = matrix(nrow = NROW(simData_test),
                                                                                   ncol = 1,
                                                                                   data = 1),
                                                                 score = "kl")
  }
}

brier_list;brier_list4
kl_list;kl_list4

#### SAMPLE PLOTS ####

t_cutoff_plot_vec <- seq(from=0.5,to=2.0,by=0.02)

#now, to choose four sample patients
cov_mat <- cbind("x1"=c(0,1,0,1),
                 "x2"=c(0,0,2,2))
rownames(cov_mat) <- LETTERS[1:4]

###AUGMENTED MODEL###
#first, generate predictions for every sample subject across time, and corresponding plots
    temp_para <- apply(tempfit$samples,MARGIN = 3,FUN = median)
    pred_c_selected <- calc_risk_logitbayes(para = temp_para, logit_ind = TRUE,
                                            Xmat1 = cov_mat, Xmat2 = cov_mat,
                                            Xmat3 = cov_mat, XmatD = cov_mat,
                                            D_eps = true_cutoff,
                                            frailty = TRUE,  #say frailty=TRUE even if it's not doesn't matter in this case
                                            model = "semi-markov",
                                            type = "conditional", gamma = 1, hazard = "weibull",
                                            t_cutoff = t_cutoff_plot_vec, tol = 1e-3,
                                            h3_tv = "piecewise",
                                            h3tv_knots = true_tv_knots,
                                            logit_tv = "piecewise",
                                            logit_tv_knots = true_cutoff)

      pred_m_selected <- calc_risk_logitbayes(para = temp_para, logit_ind = TRUE,
                                              Xmat1 = cov_mat, Xmat2 = cov_mat,
                                              Xmat3 = cov_mat, XmatD = cov_mat,
                                              D_eps = true_cutoff,
                                              frailty = TRUE,  #say frailty=TRUE even if it's not doesn't matter in this case
                                              model = "semi-markov", n_quad=15,
                                              type = "marginal", hazard = "weibull",
                                              t_cutoff = t_cutoff_plot_vec, tol = 1e-3,
                                              h3_tv = "piecewise",
                                              h3tv_knots = true_tv_knots,
                                              logit_tv = "piecewise",
                                              logit_tv_knots = true_tv_knots)

temp_list <- list()
#generate four plots (shown here for marginal)
for(i in 1:4){
  temp_list[[i]] <-
    get_sample_plot(pred_selected = pred_m_selected[,,paste0("i",i)],
                    t_cutoff_vec = t_cutoff_plot_vec, i=i, type = "logit")
}
#get one horizontal legend for everyone and then strip them of legends and labels
gg_legend_horizontal_aug <- get_legend(temp_list[[1]] + theme(legend.position = "bottom",legend.title = element_blank()))
for(i in 1:4){
  temp_list[[i]] <- temp_list[[i]] + theme(legend.position="none",axis.title=element_blank())
}

temp_grid <- cowplot::plot_grid(temp_list[[1]], temp_list[[2]],
                                temp_list[[3]], temp_list[[4]],ncol =2, align="v")
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(temp_grid,
                         bottom=grid::textGrob(label= "Gestational Age (Weeks)",
                                               gp= gpar(fontsize=14,col="black")),
                         left=grid::textGrob(label="Probability of Having Experienced Outcome",
                                             rot=90, gp= gpar(fontsize=14,col="black"))),
  gg_legend_horizontal_aug,heights=c(6,1))



#### FRAILTY-VARYING PLOTS ####

gamma_vec2 <- exp(c(-0.75,-0.25,0.25,0.75))
temp_para <- apply(tempfit$samples,MARGIN = 3,FUN = median)
temp_list <- list()
for(i in 1:length(gamma_vec2)){
  pred_c_strat_gamma <- calc_risk_logitbayes(para = temp_para, logit_ind = TRUE, D_eps=true_cutoff,
                                             Xmat1 = cov_mat[4,,drop=FALSE], Xmat2 = cov_mat[4,,drop=FALSE],
                                             Xmat3 = cov_mat[4,,drop=FALSE], XmatD = cov_mat[4,,drop=FALSE],
                                             frailty = TRUE,  #say frailty=TRUE even if it's not doesn't matter in this case
                                             model = "semi-markov",
                                             type = "conditional", gamma = gamma_vec2[i], hazard = "weibull",
                                             t_cutoff = t_cutoff_plot_vec, tol = 1e-3,
                                             h3_tv = "piecewise",
                                             h3tv_knots = true_tv_knots,
                                             logit_tv = "piecewise",
                                             logit_tv_knots = true_tv_knots)

  temp_list[[i]] <-
    get_sample_plot(pred_selected = pred_c_strat_gamma,t_cutoff_vec = t_cutoff_plot_vec,
                    i = NULL,plot_type = "stacked",type = "logit") +
    labs(subtitle = paste0("v=",log(gamma_vec2[i])))
}

#get one horizontal legend for everyone and then strip them of legends and labels
gg_legend_horizontal_aug <- get_legend(temp_list[[1]] + theme(legend.position = "bottom",legend.title = element_blank()))
for(i in 1:4){
  temp_list[[i]] <- temp_list[[i]] + theme(legend.position="none",axis.title=element_blank())
}

temp_grid <- cowplot::plot_grid(temp_list[[1]], temp_list[[2]],
                   temp_list[[3]], temp_list[[4]],ncol =2, align="v")
#finally, plot the final sample plots
gridExtra::grid.arrange(
  gridExtra::arrangeGrob(temp_grid,
                         bottom=grid::textGrob(label= "Gestational Age (Weeks)",
                                               gp= gpar(fontsize=14,col="black")),
                         left=grid::textGrob(label="Probability of Having Experienced Outcome",
                                             rot=90, gp= gpar(fontsize=14,col="black"))),
  gg_legend_horizontal_aug,heights=c(6,1))


#### PE-CONDITIONAL SAMPLE PREDICTIONS ####
#I'll predict from time 1.3, which sets h3tv1=1
t_cutoff_vec_term_imm <- seq(from=true_cutoff,to=1,length.out = 50)
cov_mat_13 <- cbind(cov_mat,
                    matrix(data=c(1,numeric(length(tv_varnames)-1)),byrow=TRUE,
                           nrow=4,ncol=length(tv_varnames),
                           dimnames = list(NULL,tv_varnames)))

###AUGMENTED MODELS###
#make a matrix of plots
temp_mat <- matrix(nrow=length(t_cutoff_vec_term_imm),
                   ncol=4,dimnames = list(NULL,LETTERS[1:4]))
for(i in 1:4){
  temp_mat[,i] <- predict_term(object = tempfit,
                               x3new = cov_mat_13[i,],
                               xDnew = cov_mat_13[i,],
                               tseq = t_cutoff_vec_term_imm,
                               D_eps=true_cutoff,
                               type = "marginal")$CIF
}

matplot(x=t_cutoff_vec_term_imm,
        y=temp_mat,
        type="l", lty=1, lwd=2,
        col = four_color_qual,
        ylim = c(0,1),
        ylab="Cumulative Probability of Delivery after Preeclampsia at Week 28",
        xlab="Gestational Age (Weeks)", add = FALSE )
legend(x="bottomright",legend = LETTERS[1:4],fill = four_color_qual,
       title = "Individual")

