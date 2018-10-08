#set working dir
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#sets working directory to source file location

source("ms1_linear_models.R")


run_combn_lm <- function(df,  alpha_val=1, yvar = 'W'){
  alpha <- alpha_val #1 is LASSO and 0 is ridge
  dat.reg <- model.matrix(~species*context*P21L*A26T*L28R, data = df)
  dat.reg <- dat.reg[,-which(colnames(dat.reg) == "(Intercept)")]
  
  all_err <- list()
  best_mod <- NULL
  
  for(i in 2:ncol(dat.reg)){
    comb_i <- combn(x = 1:ncol(dat.reg), m = i)
    pb <- txtProgressBar(1, ncol(comb_i), style=3)
    for(j in 1:ncol(comb_i)){
      cv.ij <- cv.glmnet(x = dat.reg[,comb_i[,j]], y = with(df, get(yvar)), nfolds = nrow(dat.reg), family='gaussian', alpha=alpha, standardize = FALSE, keep = TRUE, grouped = FALSE)
      err.ij <- abs(cv.ij$fit.preval[,which.min(cv.ij$cvm)] - with(df, get(yvar)))
      all_err[[paste0(i,"_",j)]] <- err.ij
      
      if(length(best_mod) == 0){
        best_mod <- paste0(i,"_",j)
        best_mods <- best_mod
        best_error <- err.ij
      }else{
        test_stat_ij <- wilcox.test(err.ij, best_error, conf.int = TRUE)
        if(test_stat_ij$p.value > 0.05){
          best_mods <- c(best_mods, paste0(i,"_",j))
        }
        if(test_stat_ij$p.value < 0.05 & mean(err.ij) < mean(best_error)){
          best_mod <- paste0(i,"_",j)
          best_mods <- best_mod
          best_error <- err.ij
        }
      }
      setTxtProgressBar(pb, j)
    }
  }
  return(list(best_mods, best_error, all_err))
}

#####
## Rafael's scrap / playground
#####

#library(olsrr)
#model <- lm(W~species*context*P21L*A26T*L28R, data  = lic)
#ols_step_best_subset(model)

#library(ggeffects)
#licpred <- ggpredict(licglm, c("P21L","A26T","L28R"))
#plot(licpred)
#ggplot(tmp, aes(x=coef, y=glmcoef))+geom_point(aes(color=glmPval <0.01))+ facet_wrap(~degree)

get_glm_pval <- function(longlmDF, dataDF){
  dataDF <- dataDF %>% mutate(P21L = as.factor(P21L), A26T = as.factor(A26T), L28R = as.factor(L28R))
  tmpglm <- glm(logW~species*context*P21L*A26T*L28R, data = dataDF, family='gaussian')
  tmptbl <-as.tibble(summary(tmpglm)$coefficients, rownames = 'term' )
  
  longlmDF %>% 
    mutate(glmcoef = tmptbl$Estimate[match(term,tmptbl$term)], glmPval = tmptbl$`Pr(>|t|)`[match(term,tmptbl$term)])
}

glm_lic <- get_glm_pval(longlm_lic, lic)
tmpDF <- lic %>% mutate(P21L = as.factor(P21L), A26T = as.factor(A26T), L28R = as.factor(L28R))
tmpglm <- glm(logW~species*context*P21L*A26T*L28R, data = tmpDF, family='gaussian')
tmptbl <-as.tibble(summary(tmpglm)$coefficients, rownames = 'term' )


#####
## Sam's scrap / playground
#####

run_aic_bic_lm <- function(df, direction = "both", method = "bic"){
  dat_no_na <- na.omit(df)
  n <- nrow(dat_no_na)
  
  mfull <- lm(y ~ species*context*P21L*A26T*L28R, data = dat_no_na)
  
  if(method == "bic"){
    k <- log(n)
  }else{
    if(method == "aic"){
      k <- 2
    }else{
      if(is.numeric(method) == TRUE){
        k <- method
      }else{
        stop("Please use 'aic' or 'bic' or a number of method")
      }
    }
  }
  step_mod <- step(object = mfull, direction = direction, k = k) #log(n) is bic, 2 is aic
  
  return(summary(step_mod)$coefficients)
}

mod_sel_abu <- run_aic_bic_lm(df = abu, direction = "both", method = "bic")
mod_sel_lic <- run_aic_bic_lm(df = lic, direction = "both", method = "bic")
mod_sel_gro <- run_aic_bic_lm(df = gro, direction = "both", method = "bic")

#figure
abu_orders <- unlist(lapply(strsplit(rownames(mod_sel_abu), ""), FUN = function(x) length(which(x == ":"))))+1

lic_orders <- unlist(lapply(strsplit(rownames(mod_sel_lic), ""), FUN = function(x) length(which(x == ":"))))+1

gro_orders <- unlist(lapply(strsplit(rownames(mod_sel_gro), ""), FUN = function(x) length(which(x == ":"))))+1

orders <- as.factor(c(abu_orders, lic_orders, gro_orders))
names <- factor(c(rep("Abundance", length(abu_orders)), rep("IC50", length(lic_orders)), rep("Drugless Growth", length(gro_orders))), levels = c("Drugless Growth", "Abundance", "IC50"))

dat.plot <- data.frame(orders, names)
p <- ggplot(dat.plot, aes(x = names, fill = orders))

quartz(width = 7, height = 7)
p + geom_bar(position = "dodge") + xlab("Phenotype") + ylab("Count") + scale_fill_manual(values = c("#e0e0e0", "#4d4d4d", "#b2182b"),  guide_legend(title = "Order")) + theme(legend.position = c(0.1,0.8), legend.key = element_rect(fill = "#f0f0f0"), legend.background = element_rect(fill = "#ffffffaa", colour = "black"), panel.background = element_rect(fill = "white", colour = "black"), axis.text.y = element_text(colour = "black", size = 14), axis.text.x = element_text(colour = "black", size = 14), axis.title = element_text(colour = "black", size = 15), panel.grid.minor = element_line(colour = "#00000050",linetype = 3), panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + scale_y_continuous(expand = c(0.01,0.01))


pal <- c("#7570b3","#666666","#66a61e")

cslm_lic$term <- factor(cslm_lic$term, levels = c("P21L", "A26T","L28R","P21L:A26T", "P21L:L28R","A26T:L28R", "P21L:A26T:L28R"))
cslm_lic$context <- factor(cslm_lic$context, levels = c("WT","GroEL","LON"))
quartz(width = 10, height = 5)
p1 <- ggplot(cslm_lic, aes(y = abs(coef), x = term, group = species, fill = species))

p1 + geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~context) + scale_fill_manual(values = pal, guide_legend(title = "Species")) + 
  xlab("Model term") + 
  ylab("Regression coefficient (abs.)") + 
  theme(legend.position = c(0.90, 0.80), 
        legend.key = element_rect(fill = "#f0f0f0"), 
        legend.background = element_rect(fill = "#ffffffaa", colour = "black"), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 14), 
        axis.text.x = element_text(colour = "black", size = 12, angle=45, hjust=1), 
        axis.title = element_text(colour = "black", size = 15), 
        panel.grid.minor = element_line(colour = "#00000050",linetype = 3), 
        panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + 
  labs(title = "IC50")

cslm_abu$term <- factor(cslm_abu$term, levels = c("P21L", "A26T","L28R","P21L:A26T", "P21L:L28R","A26T:L28R", "P21L:A26T:L28R"))
cslm_abu$context <- factor(cslm_abu$context, levels = c("WT","GroEL","LON"))

quartz(width = 10, height = 5)
p2 <- ggplot(cslm_abu, aes(y = abs(coef), x = term, group = species, fill = species))
p2 + geom_bar(stat = "identity", position = "dodge") + 
  facet_wrap(~context) + scale_fill_manual(values = pal, guide_legend(title = "Species")) + 
  xlab("Model term") + 
  ylab("Regression coefficient (abs.)") + 
  theme(legend.position = c(0.90, 0.80), 
        legend.key = element_rect(fill = "#f0f0f0"), 
        legend.background = element_rect(fill = "#ffffffaa", colour = "black"), 
        panel.background = element_rect(fill = "white", colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 14), 
        axis.text.x = element_text(colour = "black", size = 12, angle=45, hjust=1), 
        axis.title = element_text(colour = "black", size = 15), 
        panel.grid.minor = element_line(colour = "#00000050",linetype = 3), 
        panel.grid.major = element_line(colour = "#00000060", linetype = 3)) + 
  labs(title = "Abundance")
