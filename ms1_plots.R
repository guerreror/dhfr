source('ms1_linear_models.R')
library(treemapify)

pantone <- list('purple'="#66648B",
                'red'="#D34E13",
                'orange'="#FFA33E",
                'pink'="#EE5C6C",
                'blue'="#3E4F5C",
                'gray'="#D5D5D8",
                'green'="#88B04B",
                'yellow'="#F9A828", 
                "iceberg"="#71a6d2",
                "light_blue"="#8cbed6",
                "light_purple"="#d68bbd")


longlm_pal <- with(pantone, c(light_blue, yellow, green, purple, pink))

white_to_red_palette <- colorRampPalette(colors=c("white","darkred"))
white_to_red_5 <- white_to_red_palette(5)
order_colors <- c("grey90", white_to_red_5[2:5])

#####
# Fig 1. Broad view, hyper-cube -ish
#####
library(igraph, pos='package:base', quietly = T)
library(ggraph, pos='package:base', quietly = T)
library(ggimage)
paths <- function(x){
  if(is.character(x)){
    out <-  switch(x, 
                   "000" = c("001", "010", "100"),
                   "001" = c("011", "101"),
                   "010" = c("011", "110"),
                   "100" = c("101", "110"),
                   "011" = c("111"),
                   "101" = c("111"),
                   "110" = c("111"),
                   character(0)  )
  }else{
    out<- switch(as.character(x), 
                 "0" = c(1, 10, 100),
                 "1" = c(11, 101),
                 "10" = c(11, 110),
                 "100" = c(101, 110),
                 "11" = c(111),
                 "101" = c(111),
                 "110" = c(111),
                 numeric(0)  )
  }
  return(out)
}
bin_allele_count <- function(bin){
  nchar(bin) - nchar(gsub('1', '', bin))
}
bin_int_to_char <- function(x){
  out<- switch(as.character(x), 
               "0" = "000",
               "10" = "010",
               "1" = "001",
               "11" = "011",
               as.character(x) )
  return(out)
}


allbar <- allphenos %>% group_by(pheno, context, species, binary, geno) %>% summarize(W= mean(W), logW=mean(logW)) %>% ungroup()%>%
  mutate(context = factor(context, levels = c("WT", "GroEL", "LON")),
         species = factor(species, levels = c("E. coli", "C. muridarum", "L. grayi")))

mutant_plot <- function(df,xvar='mutcount', yvar='W'){
  
  nodelist <- data_frame(binary = sort(map_chr(c(0, 1, 10, 100, 11, 101, 110, 111), bin_int_to_char))) %>%
    left_join(., df, by='binary')%>%
    mutate(img = paste0("img/", binary, ".png"))
  
  xpos <- with(nodelist, get(xvar))
  ypos <- with(nodelist, get(yvar))
  
  templateNet <- data_frame(from = nodelist$binary)%>%
    mutate(to = map(from, paths) )%>%unnest() 
  
  graf <- igraph::graph_from_data_frame(templateNet, vertices = nodelist$binary)
  lay <- ggraph::create_layout(graf, 
                               layout ='manual',
                               node.positions = data.frame(x = xpos, y = ypos))
  out <- ggraph::ggraph(lay) +
    ggraph::geom_edge_link0(edge_width=0.1, color='grey30')+ 
    geom_image(data=nodelist, aes(x=mutcount, y=W, image = img), size=0.05, by='height')+
    scale_y_log10(labels = scales::trans_format("log10", scales::math_format("10"^.x)))+
    coord_cartesian(clip='off')+
    theme(legend.position = 'none',
          plot.background = element_rect(fill='grey95'),
          plot.margin=unit(c(0,0,0,0), "null"),
          #axis.ticks.x = element_blank(), 
          panel.grid.minor = element_blank(), 
          #axis.text.x = element_blank(), 
          panel.grid.major = element_blank())+
    expand_limits(x=c(-0.1, 3.1))+
    labs(x="", y="")
  out
}

licbar <- allbar %>% 
  filter(pheno=='IC50')%>% 
  arrange(species, context)%>%
  mutate(mutcount = map_int(binary,  bin_allele_count))%>%
  group_by(species, context)%>% nest()%>%
  mutate(mutplot = map(data, mutant_plot))

lic_mutplot <-(licbar$mutplot[[1]] / licbar$mutplot[[4]] / licbar$mutplot[[7]] * expand_limits(y=c(1, 1100))) | 
  (licbar$mutplot[[2]]/ licbar$mutplot[[5]] / licbar$mutplot[[8]] * expand_limits(y=c(1, 1100)) *theme(axis.text.y = element_blank(), axis.ticks = element_blank() )) |
  (licbar$mutplot[[3]]/ licbar$mutplot[[6]] / licbar$mutplot[[9]] * expand_limits(y=c(1, 1100)) *theme(axis.text.y = element_blank(), axis.ticks = element_blank()))
ggsave(filename = '1a_lic_origami.pdf', lic_mutplot, width = 6, height = 6, useDingbats=FALSE)

abubar <- allbar %>% 
  filter(pheno=='Abundance')%>% 
  arrange(species, context)%>%
  mutate(mutcount = map_int(binary, bin_allele_count))%>%
  group_by(species, context)%>% nest()%>%
  mutate(mutplot = map(data, mutant_plot))

abu_mutplot <-(abubar$mutplot[[1]] / abubar$mutplot[[4]] / abubar$mutplot[[7]] * expand_limits(y=c(300, 100000))) | 
  (abubar$mutplot[[2]]/ abubar$mutplot[[5]] / abubar$mutplot[[8]] * expand_limits(y=c(300, 100000)) *theme(axis.text.y = element_blank(), axis.ticks = element_blank())) |
  (abubar$mutplot[[3]]/ abubar$mutplot[[6]] / abubar$mutplot[[9]] * expand_limits(y=c(300, 100000)) *theme(axis.text.y = element_blank(), axis.ticks = element_blank()))
ggsave(filename = '1b_abu_origami.pdf', abu_mutplot, width = 6, height = 6, useDingbats=FALSE)



#####
# Fig 2A Total (glmnet) Regression coeffients
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

long_lm_coef_plot <- function(indf, xaxis="", main="", pal= order_colors, ncoef=20){
  df <- indf %>% filter(abs(coef)>0)%>%
    arrange(desc(abs(coef)))%>%
    filter(row_number() < ncoef+1)
  
  ggplot(df) +
    geom_col(aes(x=factor(term, levels = rev(term)), y=coef, fill=factor(degree))) +
    geom_hline(yintercept = 0, color="grey50", size = 0.1)+
    labs(y=xaxis, x='',title=main) +
    theme_classic()+
    scale_y_continuous(breaks = c(-1, 1))+
    theme(axis.title.y=element_text(angle=0, vjust=0.5),
          panel.grid.major = element_line(color = 'grey95', size = 0.3),
          axis.ticks.y = element_line(color="grey95", size = 0.3),
          axis.line.y = element_blank(),
          axis.ticks.x = element_line(color="grey50", size = 0.1),
          axis.line.x = element_line(color="grey50", size = 0.1),
          legend.position = "none")+
    scale_fill_manual(values = pal)+
    coord_flip()+
    NULL
}

longlm_lic_betaplot <- long_lm_coef_plot(prettyLM_lic, "Effect on IC50", "")
longlm_abu_betaplot <- long_lm_coef_plot(prettyLM_abu, "Effect on Abundance","")
longlm_gro_betaplot <- long_lm_coef_plot(prettyLM_gro, "Effect on\nGrowth","")

inset_longlm_plot <- function(indf, pal=order_colors){
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
# Fig 2B Treeplots
#####
square_plot <- function(indf, ncoef, termlabels = T, pal=order_colors){

  df <- indf %>% filter(row_number() < ncoef+1) %>%
    group_by(degree)%>% summarize(magnitude = sum(abs(coef)))
  
  plo <- ggplot(df, aes(area= magnitude, fill= factor(degree))) + 
    treemapify::geom_treemap(size=2, color='white')+
    treemapify::geom_treemap_text(aes(label=degree), place = "centre", grow = T, alpha = 0.5, color='white', fontface='italic', min.size = 0)+
    theme_void()+
    theme(legend.position = "none")+
    scale_fill_manual(values = pal)+
    coord_fixed()+
    NULL
  
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
# Fig 2AB Composition
#####
fig2a <- cowplot::ggdraw()+ 
  cowplot::draw_plot(longlm_lic_betaplot, x = 0, y = 0, width = 1) +
  cowplot::draw_plot(lic_square, x = 0.8, y = 0.15, width = 0.2, height=0.25)

#ggsave('2a_ic50_longLM.pdf', fig2a, width = 7, height = 5, useDingbats=FALSE)

fig2b <- cowplot::ggdraw()+ 
  cowplot::draw_plot(longlm_abu_betaplot, x = 0, y = 0, width = 1) +
  cowplot::draw_plot(abu_square, x = 0.8, y = 0.15, width = 0.2, height=0.25)

#ggsave('2b_abundance_longLM.pdf', fig2b, width = 7, height = 5, useDingbats=FALSE)

fig2comp<- cowplot::ggdraw()+ 
  cowplot::draw_plot(fig2a, x = 0, y = 0, width = 0.48) +
  cowplot::draw_plot(fig2b, x = 0.52, y = 0, width = 0.48)

ggsave('2_longLM.pdf', fig2comp, width = 9, height = 5, useDingbats=FALSE)
#####
# Fig 3. Pheno vs Pheno scatterplots
#####
versus_plot <- function(df1, df2, x1lab, x2lab, topn=5){
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
  geom_point(data= filter(vsDF, inc), color="darkred")+
  geom_text(data=filter(vsDF, inc), aes(label=term), hjust = 0, nudge_x = 0.05, size=3)+
  labs(x=paste(x1lab,"Coefficients"), y=paste(x2lab,"Coefficients"))+
  lims(x=c(-2.04, 2), y=c(-2.04, 2))+    
  theme_classic()+
  theme(legend.position = "none")
}

avl_plot <- versus_plot(longlm_abu, longlm_lic, 'Abundance', 'IC50',  5)
ggsave('3_var_correlation.pdf', avl_plot, width = 6, height = 4, useDingbats=FALSE)

lvg_plot <- versus_plot(longlm_lic, longlm_gro, 'IC50', 'Growth', 5, order_colors)
gva_plot <- versus_plot(longlm_abu, longlm_gro, 'Abundance', 'Growth', 5, order_colors)


#modelr::rsquare(glm(abu~0+coef, data=lic_vs_abu),data=lic_vs_abu)
#cor(lic_vs_abu$abu, lic_vs_abu$coef)^2


#####
# Fig 4A PQC model coefficients
#####
cs_plot <- function(indf, label_y=""){
  df <- indf %>% mutate(
    term = factor(term, levels=rev(c("P21L","A26T","L28R","P21L:A26T","A26T:L28R", "P21L:L28R","P21L:A26T:L28R"))) ,
    context = factor(context, levels = c("WT", "GroEL", "LON"), labels =c("WT", "GroEL", "Delta*italic(lon)")))
  
  ggplot(df, aes(x= term ))+ 
    geom_line(aes(y=coef, group=context, color=context), lty=3, size=0.4)+
    geom_hline(yintercept = 0, color='grey50', size=0.1)+
    geom_point(aes(y=coef, color=context, fill=context, shape=context))+
    facet_wrap(~species)+
    scale_y_continuous(breaks=c(-1, 0, 1))+
    labs(y=label_y, x="")+
    theme(axis.title.y = element_text(angle = 0, vjust=0.5),
          panel.grid.major = element_line(color = 'grey95', size = 0.3),
          panel.background = element_rect(fill = 'white'),
          panel.spacing.x = unit(2, 'lines'),
          axis.ticks.y = element_blank(),
          axis.line.x = element_line(color='grey50', size=0.1),
          legend.position = "none", 
          panel.grid = element_blank(), 
          axis.text.x = element_text(color = "grey50"),
          axis.ticks.x = element_line(color="grey50", size = 0.1),
          axis.text.y = element_text(color = "grey50"),
          strip.text = element_text(face='bold.italic', size = 10),
          strip.background = element_blank())+
    scale_fill_manual(values=c(pantone$green, pantone$orange, pantone$light_blue))+
    scale_shape_manual(values=c(21,22,24))+
    scale_color_manual(values=c(pantone$green, pantone$orange, pantone$light_blue))+
    coord_flip()+
    NULL
}

lic_cs_plot <- cs_plot(cslm_lic,"Effect on IC50")
abu_cs_plot <- cs_plot(cslm_abu, "Effect on Abundance")

#####
# Fig 4B PQC Barplots
#####

cs_barplot <- function(indf, label_y="", pal=order_colors[1:3]){
  df <- indf %>% mutate(degree = factor(str_count(term, ":")+1L), dir =  ifelse(coef >0, "(+)","(-)"))%>%
    mutate(
      term = factor(term, levels=c("P21L","A26T","L28R","P21L:A26T","A26T:L28R", "P21L:L28R","P21L:A26T:L28R")),
      context = factor(context, levels = rev(c("WT", "GroEL", "LON")), labels =rev(c("WT", "GroEL", "Delta*italic(lon)"))))%>%
    group_by(degree, context, species)%>%
    summarise(effect = sum(abs(coef)))%>%
    mutate(Order = factor(degree, levels=c("3","2","1")), labels=c("3rd", "2nd", "1st"))

  ggplot(df, aes(x= context ))+
    geom_bar(aes(y=effect, fill= Order), stat="identity", position='fill')+
    geom_text(data=group_by(df, species, context)%>%filter(row_number()==1), aes(y=0.02, label=context), hjust=0, size=2.5, color = "grey50", parse=TRUE)+
    scale_fill_manual(values = rev(pal))+
    scale_y_continuous(breaks=c(0,0.5,1), labels = c("0", "0.5", "1"))+
    facet_wrap(~species, nrow=1)+
    labs(y=label_y, x="")+
    theme(axis.text.x = element_text(color = 'grey50'),
          axis.text.y = element_blank(),
          panel.background = element_blank(),
          panel.spacing.x = unit(2, 'lines'),
          panel.grid = element_blank(),
          strip.text = element_blank(),
          axis.ticks.x = element_line(color = "grey50", size=0.1),
          axis.ticks.y = element_blank(),
          legend.position = 'none')+
    coord_flip()+
    NULL
}

licbars <- cs_barplot(cslm_lic)
abubars <- cs_barplot(cslm_abu) 

#####
# Fig 4AB Composition
#####
fig4a <-cowplot::ggdraw()+ 
  cowplot::draw_plot(lic_cs_plot, x = 0, y = 0.3, width = 1, height=0.7) +
  cowplot::draw_plot(licbars, x = 0.12, y = 0, width = 0.88, height=0.3)


ggsave('4A_ic50_withinPQC.pdf', fig4a, width = 6, height = 3, useDingbats=FALSE)

fig4b <- cowplot::ggdraw()+ 
  cowplot::draw_plot(abu_cs_plot, x = 0, y = 0.3, width = 1, height=0.7) +
  cowplot::draw_plot(abubars, x = 0.12, y = 0, width = 0.88, height=0.3)

ggsave('4B_abu_withinPQC.pdf', fig4b, width = 6, height = 3, useDingbats=FALSE)

