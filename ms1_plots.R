source('ms1_linear_models.R')
library(treemapify)
library(patchwork)

# 1) Broad picture of whole data set, across all contexts, for: (a) IC50 (b) abundance, (c) growth
# 2) IC50 glmnet coefficients: a) treemap, b) barplot of effects
# 3) Abundance glmnet coefficients: a) treemap, b) barplot of effects
# 4) Growth glmnet coefficients: a) treemap, b) barplot of effects
# 5) Coefficient scatterplots: a) Abundance vs. IC50, b) Growth vs. IC50 
# 6) Species-Context coefficients (glmnet by background)


pantone <- list('purple'="#66648B",
                'red'="#D34E13",
                'pink'="#EE5C6C",
                'blue'="#3E4F5C",
                'gray'="#D5D5D8",
                'green'="#88B04B",
                'yellow'="#F9A828", 
                "iceberg"="#71a6d2",
                "light_blue"="#8cbed6",
                "light_purple"="#d68bbd")

#####
# 1. Broad view
#####
broad_df <- allphenos %>%
  mutate(context = factor(context, levels = c("WT", "GroEL", "LON")),
         species = factor(species, levels = c("E. coli", "C. muridarum", "L. grayi")),
         binary = factor(binary, levels = c("000","001", "010", "100", "011", "101", "110", "111")),
         pheno = factor(pheno, levels = c("IC50", "Abundance", "Growth")))

broad_plot <- function(df, phenos, use_log=T){
  outdf <- df %>%mutate(W = if_else(pheno=='Growth', 10^W, W))
  if(!is.null(phenos)){
    outdf <- df %>%
    filter(pheno %in% phenos)
  }
   plo <- ggplot(outdf,aes(x=context, y=W))+
     geom_point(aes(color=species,fill=species, shape=species), alpha=0.5, position = position_jitter(0.01))+
     facet_grid(pheno~binary, scales='free_y', switch = "y")+
     stat_summary(fun.y = "mean", geom ="line", aes(group=species, color=species), lty=1)+
     theme(panel.grid = element_blank(), 
           strip.background = element_blank(),
           axis.text.x = element_text(angle = 45, hjust = 1), 
           strip.text.y = element_text(angle = 180),
           strip.placement = "outside",
           legend.position = "none", 
           legend.direction = "horizontal",
           legend.text = element_text(face = "italic"), 
           legend.key = element_blank(), 
           legend.title = element_blank(),
           legend.box.background = element_rect(fill='white', color="white"))+
     scale_color_manual(values=c(pantone$green, pantone$light_blue, pantone$pink))+
     scale_fill_manual(values=c(pantone$green, pantone$light_blue, pantone$pink))+
     scale_shape_manual(values=c(21,23,22))+
     labs(x="", y="")+
     scale_x_discrete(labels=c("WT"="WT","GroEL"="GroEL", "LON"=expression(Delta*italic("lon"))))
   if(use_log){ return(plo + 
                         theme(strip.text.x = element_text(face="bold"))+
                         scale_y_log10(breaks = scales::log_breaks(n = 10), labels = scales::comma) )
     }
   else{return(plo + 
                 theme(strip.text.x = element_blank())+
                 scale_y_continuous(labels = scales::comma) )
     }
}

broad_all <- broad_plot(broad_df, NULL)
ggsave('raw_1_broadview.pdf', broad_all, width = 9, height = 7, useDingbats=FALSE)

broad_gro <- broad_plot(broad_df, "Growth", use_log = F)
ggsave('raw_broadview_gro.pdf', broad_gro, width = 9, height = 2.5, useDingbats=FALSE)


#####
# 2A,3A,4A Total (glmnet) Regression coeffients
#####

fix_lm_terms <- function(df){
  df %>% mutate(spp = str_detect(term, "species"),
                  ctx = str_detect(term, "context"),
                  geno = str_detect(term, "[0-9]+"),
                  degree = str_count(term, ":")+1L)%>%
    mutate(term = str_remove_all(term, 'species'),
           term = str_replace_all(term, '. ', '.'),
           term = str_remove_all(term, 'context'),
           term = str_replace_all(term, "1:",":"),
           term = str_remove_all(term, "1$"))
}

prettyLM_lic <- fix_lm_terms(longlm_lic)
prettyLM_abu <- fix_lm_terms(longlm_abu)
prettyLM_gro <- fix_lm_terms(longlm_gro)

longlm_pal <- with(pantone, c(light_blue, yellow, green, purple, pink))

long_lm_coef_plot <- function(indf, xaxis, main, pal= longlm_pal, ncoef=20){
  df <- indf %>% filter(abs(coef)>0)%>%
    arrange(desc(abs(coef)))%>%
    filter(row_number() < ncoef+1)%>%
    arrange(degree, desc(abs(coef)))
  
  ggplot(df) + geom_col(aes(x=factor(term, levels = term), y=coef, fill=factor(degree) )) +
    labs(y=xaxis, x='',title=main) +
    theme_classic()+
    theme(plot.margin=unit(c(0,0,0,0),"cm"),
          axis.title.y=element_text(angle=0, vjust=0.5),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")+
    scale_fill_manual(values = pal)
}

longlm_lic_betaplot <- long_lm_coef_plot(prettyLM_lic, "Effect on\nIC50", "")
longlm_abu_betaplot <- long_lm_coef_plot(prettyLM_abu, "Effect on\nAbundance","")
longlm_gro_betaplot <- long_lm_coef_plot(prettyLM_gro, "Effect on\nGrowth","")

inset_longlm_plot <- function(indf, pal=longlm_pal){
  df <- indf %>% group_by(degree) %>% summarize(effect = sum(abs(coef)))
  
  ggplot(df)+ geom_col(aes(x=factor(degree), y=effect, fill=factor(degree)))+
    theme_classic()+
    geom_text(aes(x=factor(degree), y=0, label=degree), vjust=0, size=9, color='white', fontface='italic')+
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank())+
    scale_fill_manual(values = pal)+
    labs(x="", y="Total effect")
}
longlm_inset_lic <- inset_longlm_plot(prettyLM_lic)
longlm_inset_abu <- inset_longlm_plot(prettyLM_abu)
longlm_inset_gro <- inset_longlm_plot(prettyLM_gro)

#####
# 2B,3B,4B Treeplots
#####
square_plot <- function(indf, ncoef, termlabels = T){
  pal <- with(pantone, c(light_blue, yellow, green, purple, pink))
  
  df <- indf %>% filter(row_number() < ncoef+1)
  
  plo <- ggplot(df, aes(area=abs(coef), fill= factor(degree), subgroup= degree)) + 
    geom_treemap(size=1.5, color='white')+
    geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.3, color='white', fontface='italic', min.size = 0)+
    geom_treemap_subgroup_border(size=6, color='white')+
    theme(legend.position = "none")+
    scale_fill_manual(values = pal)
  
  if(termlabels){
    return(plo +
             geom_treemap_text(aes(label=term),colour = "white", place = "centre", reflow = T, size=14, padding.x=unit(4, units = "pt"), padding.y=unit(4, units = "pt"))  )
  } else{
    return(plo) 
    }
}

lic_square <- square_plot(prettyLM_lic, 100, F)
abu_square <- square_plot(prettyLM_abu, 100, F)
gro_square <- square_plot(prettyLM_gro, 100, F)

#lic_square_labeld <- square_plot(prettyLM_lic, 20, T)
#abu_square_labeld <- square_plot(prettyLM_abu, 20, T)
#ggsave('ic50_treemap_labelled.pdf', lic_square_labeld, width = 8, height = 8, useDingbats=FALSE)
#ggsave('abundance_treemap_labelled.pdf', abu_square_labeld, width = 8, height = 8, useDingbats=FALSE)

#####
# 2,3,4 Composition
#####
fig2 <- cowplot::ggdraw()+ 
  cowplot::draw_plot(longlm_lic_betaplot, x = 0.02, y = 0, width = 0.6) +
  cowplot::draw_plot(lic_square, x = 0.65, y = 0.28, width = 0.33, height=0.68)

ggsave('2_ic50_longLM.pdf', fig2, width = 10, height = 4, useDingbats=FALSE)

fig3 <- cowplot::ggdraw()+ 
  cowplot::draw_plot(longlm_abu_betaplot, x = 0.02, y = 0, width = 0.6) +
  cowplot::draw_plot(abu_square, x = 0.65, y = 0.28, width = 0.33, height=0.68)

ggsave('3_abundance_longLM.pdf', fig3, width = 10, height = 4, useDingbats=FALSE)

#####
# 5. Pheno vs Pheno scatterplots
#####
versus_plot <- function(df1, df2, x1lab, x2lab, topn=5, pal){
  vsDF <- df1 %>% mutate(x2 = df2$coef[match(term, df2$term)]) %>%
    mutate(spp = str_detect(term, "species"),
           ctx = str_detect(term, "context"),
           geno = str_detect(term, "[0-9]+"),
           degree = str_count(term, ":")+1L)%>%
    mutate(term = str_remove_all(term, 'species'),
           term = str_replace_all(term, '. ', '.'),
           term = str_remove_all(term, 'context'),
           term = str_replace_all(term, "1:",":"),
           term = str_replace_all(term, ":",": "),
           term = str_remove_all(term, "1$"))%>%
    mutate(dif = abs(coef - x2) )%>% arrange(desc(dif))%>%
    mutate(inc = row_number()<topn +1)

ggplot(vsDF, aes(x=coef, y=x2))+ 
  geom_abline(intercept = 0, slope=1, lty=2, alpha=0.5)+
  geom_hline(yintercept = 0)+ geom_vline(xintercept = 0)+
  geom_point(data= filter(vsDF, !inc), color="grey80", alpha=0.7)+
  geom_point(data= filter(vsDF, inc), aes(color=factor(degree)))+
  geom_text(data=filter(vsDF, inc), aes(label=term), hjust = 0, nudge_x = 0.05, size=3)+
  labs(x=paste(x1lab,"Coefficients"), y=paste(x2lab,"Coefficients"))+
  lims(x=c(-2.04, 2), y=c(-2.04, 2))+    
  theme_classic()+
  theme(legend.position = "none")+
  scale_color_manual(values = pal)
}

avl_plot <- versus_plot(longlm_abu, longlm_lic, 'Abundance', 'IC50',  5, longlm_pal)

lvg_plot <- versus_plot(longlm_lic, longlm_gro, 'IC50', 'Growth', 5, longlm_pal)
gva_plot <- versus_plot(longlm_abu, longlm_gro, 'Abundance', 'Growth', 5, longlm_pal)

fig5 <- avl_plot + lvg_plot + plot_layout(ncol=2)

ggsave('5_coefs_scatterplots.pdf', avl_plot, width = 6, height = 4, useDingbats=FALSE)

#modelr::rsquare(glm(abu~0+coef, data=lic_vs_abu),data=lic_vs_abu)
#cor(lic_vs_abu$abu, lic_vs_abu$coef)^2

#####
# 6. linear models within Context*Species
#####

cs_plot <- function(df, label_y){
  ggplot(df, aes(x=factor(term, levels=c("P21L","A26T","L28R","P21L:A26T","A26T:L28R", "P21L:L28R","P21L:A26T:L28R")) ))+ 
    geom_hline(aes(yintercept=0), color='grey99')+
    geom_line(aes(y=coef, group=species, color=species), lty=3)+
    geom_point(aes(y=coef, color=species, fill=species, shape=species))+
    facet_wrap(~factor(context, levels = c("WT", "GroEL", "LON"), labels =c("WT", "GroEL", "Delta*italic(lon)")), labeller=label_parsed)+
    labs(y=label_y, x="")+
    theme(axis.title.y = element_text(angle = 0, vjust=0.5),
          legend.position = "none", 
          panel.grid = element_blank(), 
          panel.spacing.y=unit(2, "lines"),
          axis.text.x = element_text(angle = 45, hjust = 1),
          strip.background = element_blank(),
          strip.placement = "outside")+
    scale_fill_manual(values=c(pantone$green, pantone$light_blue, pantone$pink))+
    scale_shape_manual(values=c(21,23,22))+
    scale_color_manual(values=c(pantone$green, pantone$light_blue, pantone$pink))
}

lic_cs_plot <- cs_plot(cslm_lic, "Effect on\nIC50")+ 
  theme( strip.text.x = element_blank(),
         axis.text.x = element_blank(),
         axis.ticks.x = element_blank())


abu_cs_plot <- cs_plot(cslm_abu, "Effect on\nAbundance")+
  theme(legend.justification=c(1,1),
        legend.position = c(0,0),
        legend.title = element_blank(),
        legend.box.background = element_blank(),
        legend.background = element_blank(),
        legend.key = element_rect(color=NA, fill=NA),
        legend.text = element_text(face="italic"))

fig6 <- lic_cs_plot + abu_cs_plot + patchwork::plot_layout(ncol = 1)

ggsave('6_csLM.pdf', fig6, width = 8, height = 6, useDingbats=FALSE)



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

