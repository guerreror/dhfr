source('ms1_linear_models.R')
#####
# SupData export
#####

alpha05_lic <- run_long_lm(lic, 0.5) %>% fix_lm_terms()  %>% transmute(Term =term, `Estimate (Alpha = 0.5)`=coef)
alpha05_abu <- run_long_lm(abu, 0.5) %>% fix_lm_terms()  %>% transmute(Term =term, `Estimate (Alpha = 0.5)`=coef)
alpha05_gro <- run_long_lm(gro, 0.5) %>% fix_lm_terms()  %>% transmute(Term =term, `Estimate (Alpha = 0.5)`=coef)

tableS1_lic <- prettyLM_lic %>% transmute(Term =term, `Estimate (Alpha = 1)`=coef) %>% left_join(y=alpha05_lic, by='Term')
tableS1_abu <- prettyLM_abu %>% transmute(Term =term, `Estimate (Alpha = 1)`=coef) %>% left_join(y=alpha05_abu, by='Term')
tableS1_gro <- prettyLM_gro %>% transmute(Term =term, `Estimate (Alpha = 1)`=coef) %>% left_join(y=alpha05_gro, by='Term')
openxlsx::write.xlsx(tableS1_lic, 'table_S1_ic50.xlsx')
openxlsx::write.xlsx(tableS1_abu, 'table_S2_abundance.xlsx')
openxlsx::write.xlsx(tableS1_gro, 'table_S3_growth.xlsx')

openxlsx::write.xlsx(rename(cslm_lic, Term =term, `Estimate`=coef), 'table_S4_ic50_PQC_models.xlsx')
openxlsx::write.xlsx(rename(cslm_abu, Term =term, `Estimate`=coef), 'table_S5_abundance_PQC_models.xlsx')


#####
# Power-Transform linearization (Figs S2, S3)
#####

#get predictions from an additive model (phenotype ~ species + context + genotype)
run_additive_lm <- function(df){
  df$y <- df$y + abs(min(0, min(df$y)))
  lin <- glm(y~species+context+P21L+A26T+L28R, data = df)
  padd <- predict(lin, type='response')
  df%>%mutate(p_add = padd)
}

#get geometric mean
gm_mean <- function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

# This function is the power-transform described in Sailer & Harms 2017
powertrans <- function(padd, p){
  a = p[1]
  b = p[2]
  l = p[3]
  gm <- gm_mean(padd+a)
  up <- (padd + a)^l - 1
  lo <- l * gm^(l-1)
  out <- up/lo + b
  out[is.nan(out)] <- NA
  out
}
# Sums of squared to be minimized
sums_sq <- function(pobs, padd){
  function(p){
    p_trans <- powertrans(padd, p)
    sum((p_trans - pobs)^2)
  }
} 
# This function is Eq. 2 in Sailer & Harms 2017
backtrans <- function(pobs, padd, p){
  a = p[1]
  b = p[2]
  l = p[3]
  gm <- gm_mean(padd+a)
  (l * gm^(l-1) * (pobs-b) + 1)^(1/l) - a
}

plot_supp_fig_1 <- function(df, yvar){
  ggplot(df, aes(x=p_add, y= y))+
    geom_abline(slope=1, intercept = 0, lty=2)+
    geom_point()+
    geom_smooth(method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
    theme_minimal()+
    labs(y=yvar, x="Predicted by additive model")
}

abu_linear <- run_additive_lm(abu)
abu_fig_s1 <- plot_supp_fig_1(abu_linear, "Observed log10(abundance)")

lic_linear <- run_additive_lm(lic)
lic_fig_s1 <- plot_supp_fig_1(lic_linear, "Observed log10(IC50)")

abu_powfn <- sums_sq(abu_linear$y, abu_linear$p_add)
abu_est <- nlm(abu_powfn, c(0,0,1.1), iterlim=10000)$estimate
abu_bc_transf <- abu_linear %>% 
  mutate(nu = map_dbl(y, backtrans, padd=p_add, p = abu_est))

fig_s1 <- cowplot::plot_grid(abu_fig_s1, lic_fig_s1, labels = c("A", "B"), ncol = 1, scale=0.96)
ggsave('s2_powertrans.pdf', fig_s1, width = 5, height = 8, useDingbats=FALSE)

# ggplot(abu_bc_transf, aes(x=p_add, y=nu))+
#   geom_point()+
#   geom_smooth(color='green',method = lm, formula = y ~ splines::bs(x, 3), se = FALSE)+
#   geom_abline(slope=1, intercept = 0, lty=2)+
#   lims(x = c(2.75,5), y=c(2,6))

longlm_abu_nu <- run_long_lm(abu_bc_transf, myalpha,  yname = 'nu')%>%rename(nucoef=coef)

ctrast_abu <- full_join(longlm_abu_nu, longlm_abu, by='term')
cstrast_plot <- ggplot(ctrast_abu)+
  geom_abline(slope=1, intercept = 0, lty=2)+
  geom_point(aes(x=coef, y=nucoef))+
  labs(x='log10(Abundance) LASSO Coefficients', y='Linearized Abundance LASSO Coefficients')+
  theme_minimal()
ggsave('s3_abu_log_v_nu.pdf', cstrast_plot, width = 6, height = 5, useDingbats=FALSE)


# prettyLM_abu_nu <- fix_lm_terms(longlm_abu_nu)
# longlm_abu_nu_betaplot <- long_lm_coef_plot(prettyLM_abu_nu, "Effect on Abundance","")
# abu_nu_square <- square_plot(prettyLM_abu_nu, 100, F)
# fig3_nu <- cowplot::ggdraw()+ 
#   cowplot::draw_plot(longlm_abu_nu_betaplot, x = 0, y = 0, width = 1) +
#   cowplot::draw_plot(abu_nu_square, x = 0.8, y = 0.15, width = 0.2, height=0.25)
# ggsave('abu_nu_abundance_longLM.pdf', fig3_nu, width = 7, height = 5, useDingbats=FALSE)
# 

# cslm_abu_nu <- run_cs_lm(abu_bc_transf, myalpha, yvar = 'nu')
# abu_nu_cs_plot <- cs_plot(cslm_abu_nu, "Effect on Abundance")
# abu_nu_bars <- cs_barplot(cslm_abu_nu) 
# fig6_nu <- cowplot::ggdraw()+ 
#   cowplot::draw_plot(abu_nu_cs_plot, x = 0, y = 0.3, width = 1, height=0.7) +
#   cowplot::draw_plot(abu_nu_bars, x = 0.12, y = 0, width = 0.88, height=0.3)
# ggsave('abu_nu_withinPQC.pdf', fig6_nu, width = 6, height = 3, useDingbats=FALSE)

