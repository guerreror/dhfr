source("ms1_load_data.R")
library(glmnet)

prep_for_lm <- function(df, phenotype, with_log=TRUE){
 
  if(with_log){
    outdf <- df %>% dplyr::rename(y = logW, yrel = logw) %>% filter(pheno == phenotype)
  } else{
    outdf <- df %>% dplyr::rename(y = W, yrel = w) %>% filter(pheno == phenotype)  
    }
  mean1 <- mean(outdf$y)
  sd1 <- sd(outdf$y)

  outdf %>% 
    separate(binary, c("P21L", "A26T","L28R"), sep = c(1,2))%>%
    transmute(species, context, P21L, A26T, L28R, y, yrel)%>%
    mutate(species = factor(species, levels = c("E. coli", "C. muridarum", "L. grayi")),
           context = factor(context, levels = c("WT", "GroEL", "LON")))%>%
    mutate(normy = (y - mean1)/sd1)
}

# LOG10 IC50 data with replicates
lic <- prep_for_lm(allphenos,'IC50')

# LOG10 Abundance data with replicates
abu <- prep_for_lm(allphenos,'Abundance')

# Growth data with replicates
gro <- prep_for_lm(allphenos,'Growth', with_log = FALSE)

#####
# Total Lasso (glmnet) models (aka. "long") 
#####
run_long_lm <- function(df, alpha_val=1){
  alpha <- alpha_val #1 is LASSO and 0 is ridge
  dat.reg <- model.matrix(~species*context*P21L*A26T*L28R, data = df)
  dat.reg <- dat.reg[,-which(colnames(dat.reg) == "(Intercept)")]
  
  mycv <- cv.glmnet(x = dat.reg, 
                    y = df$normy, 
                    nfolds = nrow(dat.reg), 
                    family='gaussian', 
                    alpha=alpha, 
                    standardize = TRUE,
                    keep = TRUE, 
                    grouped = FALSE,
                    intercept = FALSE)

  coefs <- coef(mycv, s = "lambda.1se")
  ## https://stats.stackexchange.com/questions/70249/feature-selection-model-with-glmnet-on-methylation-data-pn
  ## This is a question about parsimony. The lambda.min option refers to value of λ at the lowest CV error. 
  ## The error at this value of λ is the average of the errors over the k folds and hence this estimate of the error is uncertain. 
  ## The lambda.1se represents the value of λ in the search that was simpler than the best model (lambda.min), 
  ## but which has error within 1 standard error of the best model. In other words, using the value of lambda.1se as the selected
  ## value for λ results in a model that is slightly simpler than the best model but which cannot be distinguished from the best model
  ## in terms of error given the uncertainty in the k-fold CV estimate of the error of the best model.
  #mod <- glmnet(x = dat.reg, y = df$logW, family='gaussian', alpha=alpha, standardize= FALSE, lambda = mycv$lambda.min)
  
  tibble(term = rownames(coefs), coef = unname(coefs[,1]) )%>% arrange(desc(abs(coef)))
}

myalpha <- 1
longlm_lic <- run_long_lm(lic, myalpha)
longlm_abu <- run_long_lm(abu, myalpha)
longlm_gro <- run_long_lm(gro, myalpha)

#####
# Lasso within Context+Species (aka "cs models", "Rafael's models", "profiles")
#####
one_cs_lm <-function(df, alpha_val=1, yvar = 'y') {
  
  mean1 <- mean(with(df, get(yvar)))
  sd1 <- sd(with(df, get(yvar)))
  local_y <- (with(df, get(yvar)) - mean1)/sd1
  
  alpha <- alpha_val 
  dat.reg <- model.matrix(~P21L*A26T*L28R, data = df)
  dat.reg <- dat.reg[,-which(colnames(dat.reg) == "(Intercept)")]
  
  mycv <- cv.glmnet(x = dat.reg, 
                    y = local_y, 
                    nfolds = nrow(dat.reg), 
                    family='gaussian', 
                    alpha=alpha, 
                    standardize = FALSE, 
                    keep = TRUE, 
                    grouped = FALSE,
                    intercept = FALSE)
  
  coefs <- coef(mycv, s = "lambda.1se")
  #https://stats.stackexchange.com/questions/70249/feature-selection-model-with-glmnet-on-methylation-data-pn  
  #mod <- glmnet(x = dat.reg, y = with(df, get(yvar)), family='gaussian', alpha=alpha, standardize= FALSE, lambda = mycv$lambda.min)
  
  tibble(term = rownames(coefs), coef = unname(coefs[,1]) )%>% arrange(desc(abs(coef)))
}

run_cs_lm <- function(df, alpha=1){
  df %>% 
    group_by(species, context)%>%
    nest(.key = vals)%>%
    mutate(model = map(vals, one_cs_lm, alpha_val=alpha))%>%
    unnest(model)%>%
    filter(term !="(Intercept)")%>%
    mutate(term = str_replace_all(term, "1:",":"),
           term = str_remove_all(term, "1$"))
}

cslm_lic <- run_cs_lm(lic, myalpha)
cslm_abu <- run_cs_lm(abu, myalpha)
cslm_gro <- run_cs_lm(gro, myalpha)

