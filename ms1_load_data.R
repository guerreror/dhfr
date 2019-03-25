library(tidyverse)

######
## load data
######
read_single_sheet <- function(x, f){
  readxl::read_xlsx(f, sheet = x) 
}

# read all phenotypic data, includes IC50
datafile <- "GGG Whole Raw Set Clean.xlsx"

grodata <- read_csv('zero_growth_rate_reps.csv', col_types = '-ccccc-cd') %>%
  unite(binary, contains('locus'), sep="") %>% 
  transmute(pheno='Growth', context, species, geno, binary, W = growth)%>%na.omit()

raw_data <- tibble( shnames= readxl::excel_sheets(datafile)) %>% 
  mutate(content = map(shnames, read_single_sheet, f = datafile))%>%
  separate(shnames, c("pheno", "context"),sep = " ")%>%
  unnest() %>% 
  dplyr::rename(species = Bacteria, geno = Genotype, binary = Binary)%>%
  mutate(species =  str_replace_all(species, "miridarum", "muridarum"))

# Makes relative fitness within context+species
get_refbar <- function(df, ph, sp, ctx, bin, x = "logW"){
  if(!is.null(ph)) {  out<- df%>%filter(species==sp & binary == bin & context == ctx & pheno == ph)}
    else{out<- df%>%filter(species==sp & binary == bin & context == ctx)}
  mean(with(out, get(x)))
} 

allphenos <- raw_data %>%
  gather(key = 'rep', value='W', starts_with("Rep"))%>%
  mutate(rep = NULL)%>%
  bind_rows(., grodata)%>%
  mutate(logW = log10(W))%>%
  na.omit()
  
allphenos <- allphenos %>% 
  mutate(logrefbar = pmap_dbl(list(ph=pheno, sp=species, ctx = context), get_refbar, df = allphenos, bin="000"))%>%
  mutate(refbar = pmap_dbl(list(ph=pheno, sp=species, ctx = context), get_refbar, df = allphenos, bin="000", x="W"))%>%
  group_by(pheno)%>%
  mutate(logw = logW/logrefbar,
         w = W/refbar,
         refbar=NULL, logrefbar=NULL)%>%
  ungroup()

## after checking normality, we're going with log(IC50) and log(abundance)
#qqnorm(allphenos%>%filter(pheno=="IC50")%>%.$W)
#qqnorm(allphenos%>%filter(pheno=="IC50")%>%.$logW)
#qqnorm(allphenos%>%filter(pheno=="Abundance")%>%.$W)
#qqnorm(allphenos%>%filter(pheno=="Abundance")%>%.$logW)
