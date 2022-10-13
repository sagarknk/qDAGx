if(!require("ggplot2",quietly = TRUE)){install.packages("ggplot2")}
if(!require("gridExtra", quietly = TRUE)){install.packages("gridExtra")}
if(!require("dplyr",quietly = TRUE)){install.packages("dplyr")}
if(!require("scales", quietly = TRUE)){install.packages("scales")}

num_of_p = 5
num_of_q = 2
sample_size = 10
KT_val = 0.5

if(num_of_q == 2)
{
  true_th = "0p5"
}else if(num_of_q == 5)
{
  true_th = "1"
}

if(KT_val == 0.5)
{
  kendall_value = "0p5"
}else if(KT_val == 0.25)
{
  kendall_value = "0p25"
}

  file_name = paste0("HS_All_in_one_result_mean_with_th_",true_th,"_",num_of_p,"_q_",num_of_q,"_n_",sample_size,".csv")
  results_HS = read.csv(file_name, header = FALSE)
  results_HS = as.matrix(results_HS)
  
  file_name = paste0("Union_DAG_All_in_one_result_mean_with_th_",true_th,"_",num_of_p,"_q_",num_of_q,"_n_",sample_size,".csv")
  results_ud = read.csv(file_name, header = FALSE)
  results_ud = as.matrix(results_ud)

  file_name = paste0("Misspecified_HS_All_in_one_result_with_KT_",kendall_value,"_mean_th_",true_th,"_",num_of_p,"_q_",num_of_q,"_n_",sample_size,".csv")
  results_misspec_hs = read.csv(file_name, header = FALSE)
  results_misspec_hs = as.matrix(results_misspec_hs)
  
  file_name = paste0("HS_All_in_one_result_std_with_th_",true_th,"_",num_of_p,"_q_",num_of_q,"_n_",sample_size,".csv")
  sd_HS = read.csv(file_name, header = FALSE)
  sd_HS = as.matrix(sd_HS)
  
  file_name =  paste0("Union_DAG_All_in_one_result_std_with_th_",true_th,"_",num_of_p,"_q_",num_of_q,"_n_",sample_size,".csv")
  sd_ud = read.csv(file_name, header = FALSE)
  sd_ud = as.matrix(sd_ud)
  
  file_name = paste0("Misspecified_HS_All_in_one_result_with_KT_",kendall_value,"_std_th_",true_th,"_",num_of_p,"_q_",num_of_q,"_n_",sample_size,".csv")
  sd_misspec_hs = read.csv(file_name, header = FALSE)
  sd_misspec_hs = as.matrix(sd_misspec_hs)
  
  df = data.frame(Method =  factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  TPR_P = c(results_HS[4,], results_ud[1,], results_misspec_hs[4,]), 
                  sd = c(sd_HS[4,],sd_ud[1,], sd_misspec_hs[4,]))
  
  g1 = ggplot(df, aes(x=Tau, y=TPR_P, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = TPR_P - sd, ymax = TPR_P + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression("TPR"[Y]^tau))
  
  df = data.frame(Method =  factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  FPR_P = c(results_HS[5,], results_ud[2,],results_misspec_hs[5,]), 
                  sd = c(sd_HS[5,], sd_ud[2,], sd_misspec_hs[5,]))
  
  g2 = ggplot(df, aes(x=Tau, y=FPR_P, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = FPR_P - sd, ymax = FPR_P + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression("FPR"[Y]^tau))
  
  df = data.frame(Method =  factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  AUC_P = c(results_HS[6,], results_ud[3,],results_misspec_hs[6,]), 
                  sd = c(sd_HS[6,], sd_ud[3,],sd_misspec_hs[6,]))
  
  g3 = ggplot(df, aes(x=Tau, y=AUC_P, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = AUC_P - sd, ymax = AUC_P + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression("AUC"[Y]^tau))
  
  df = data.frame(Method =  factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  TPR_M = c(results_HS[7,], results_ud[4,],results_misspec_hs[7,]), 
                  sd = c(sd_HS[7,], sd_ud[4,],sd_misspec_hs[7,]))
  
  g4 = ggplot(df, aes(x=Tau, y=TPR_M, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = TPR_M - sd, ymax = TPR_M + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression("TPR"[X]^tau))
  
  df = data.frame(Method =  factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  FPR_M = c(results_HS[8,], results_ud[5,],results_misspec_hs[8,]), 
                  sd = c(sd_HS[8,], sd_ud[5,],sd_misspec_hs[8,]))
  
  g5 = ggplot(df, aes(x=Tau, y=FPR_M, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = FPR_M - sd, ymax = FPR_M + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression("FPR"[X]^tau))
  
  df = data.frame(Method =  factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  AUC_M = c(results_HS[9,], results_ud[6,],results_misspec_hs[9,]), 
                  sd = c(sd_HS[9,], sd_ud[6,],sd_misspec_hs[9,]))
  
  g6 = ggplot(df, aes(x=Tau, y=AUC_M, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = AUC_M - sd, ymax = AUC_M + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression("AUC"[X]^tau))
  
  df = data.frame(Method =  factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  MSE = c(results_HS[10,], results_ud[7,], results_misspec_hs[10,]), 
                  sd = c(sd_HS[10,], sd_ud[7,],sd_misspec_hs[10,]))
  
  g7 = ggplot(df, aes(x=Tau, y=MSE, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = MSE - sd, ymax = MSE + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression("MSE"^tau))
 
   df = data.frame(Method =  factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  Theta_diff_norm = c(results_HS[1,],results_ud[9,], results_misspec_hs[1,]), 
                  sd = c(sd_HS[1,], sd_ud[9,],sd_misspec_hs[1,]))
  
  g8 = ggplot(df, aes(x=Tau, y=Theta_diff_norm, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = Theta_diff_norm - sd, ymax = Theta_diff_norm + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression(Delta[F]*theta^tau))
  
  ############# Attempt to get subscripts in legend #############

  
  df = data.frame(Method = factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  Beta_diff_norm = c(results_HS[2,], results_ud[10,], results_misspec_hs[2,]), 
                  sd = c(sd_HS[2,], sd_ud[10,], sd_misspec_hs[2,]))
  
  df = df %>% 
    mutate(Method = recode_factor(Method, `a` = "qDAGx[0]", `b` = "qDAGx",`c`="qDAGx[m]"))
  
  g9 = ggplot(df, aes(x=Tau, y=Beta_diff_norm, group=Method)) +
    geom_line(aes(color=Method), size = 2)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.direction = "horizontal", 
          legend.position = "bottom",
          legend.box = "horizontal",
          legend.text = element_text(size=8, face="bold"),
          legend.title = element_text(size=8))+
    geom_errorbar(aes(ymin = Beta_diff_norm - sd, ymax = Beta_diff_norm + sd), width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(col = expression('Method:'~'')) + 
    scale_colour_discrete(labels = parse_format())
  
  #########################################
  get_legend <- function(p) {
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
  }
  leg = get_legend(g9)
  ############################################
  df = data.frame(Method = factor(rep(c("a","b","c"),each = 9)),
                  Tau = c(rep(seq(0.1,0.9,0.1),3)),
                  Beta_diff_norm = c(results_HS[2,], results_ud[10,], results_misspec_hs[2,]), 
                  sd = c(sd_HS[2,], sd_ud[10,], sd_misspec_hs[2,]))
  
  g9 = ggplot(df, aes(x=Tau, y=Beta_diff_norm, group=Method)) +
    geom_line(aes(color=Method), size = 0.25)+
    geom_point(aes(color=Method), size = 0.125, stroke = 0.2, shape = 19, color="black")+
    theme(legend.position="none", axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold") )+
    geom_errorbar(aes(ymin = Beta_diff_norm - sd, ymax = Beta_diff_norm + sd), 
                  width=0.05, size=0.2, alpha=0.5,
                  position=position_dodge(0.05))+
    labs(x = expression(tau), y = expression(Delta[F]*beta^tau))
  
  
  grid_plot = grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8,g9,leg, layout_matrix = rbind(c(1,2,3),c(4,5,6),c(7,8,9),c(10,10,10)),
               heights = c(1,1,1,0.25))