library(reshape2)
library(reshape)
library(dplyr)

# upload molten gall network data
gall_net_melt <- read.csv("~/Documents/Genotype_Networks/data/gall_network_data.csv")
gall_net_melt <- tbl_df(gall_net_melt)

##### Does the proportion of parasitized vLG galls vary among genotypes?

### contingency table analyses for each gall species separately. 

# Iteomyia salicisverruca
vLG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("vLG")) %>%
  mutate(vLG_ptoid_attack = vLG_Eulo.fem + vLG_Eulo.mal + vLG_Mesopol + vLG_Mymarid + vLG_Platy + vLG_Ptero.2 + vLG_Tory.fem + vLG_Tory.mal, vLG_prop_ptoid_attack = vLG_ptoid_attack/(vLG_vLG.pupa + vLG_ptoid_attack)) %>% # omitting exit hole because of potential overlap with other ectoparasitoids. Need to probably resolve this by not using "unk" gall.id
  select(Genotype, vLG_ptoid_attack, vLG_vLG.pupa, vLG_prop_ptoid_attack)

vLG_for_chi <- vLG_geno_parasitism %>%
  filter(vLG_ptoid_attack > 0 & vLG_vLG.pupa > 0) %>% # remove all observations with zeros in either column
  select(vLG_ptoid_attack, vLG_vLG.pupa)

chisq.test(vLG_for_chi)
chisq.test(vLG_for_chi, simulate.p.value = TRUE, B = 10000) # same results as above

# Rabdophaga salicisbrassicoides
rG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("rG")) %>%
  mutate(rG_ptoid_attack = rG_diff.or.larv + rG_Eulo.fem + rG_Eulo.mal + rG_Lestodip + rG_Mesopol + rG_Platy + rG_Tory.fem + rG_Tory.mal + rG_unk.ptoid, rG_prop_ptoid_attack = rG_ptoid_attack/(rG_rG.larv + rG_ptoid_attack)) %>% # omitting exit hole because of potential overlap with other ectoparasitoids. Need to probably resolve this by not using "unk" gall.id
  select(Genotype, rG_ptoid_attack, rG_rG.larv, rG_prop_ptoid_attack)

rG_for_chi <- rG_geno_parasitism %>%
  filter(rG_ptoid_attack > 0 & rG_rG.larv > 0) %>% # remove all observations with zeros in either column
  select(rG_ptoid_attack, rG_rG.larv)

chisq.test(rG_for_chi) # qualitatively the same even if I only remove observations with > 0 in either ptoid attack or rG.larva
chisq.test(rG_for_chi, simulate.p.value = TRUE, B = 10000)

# Rabdophaga salicisbattatus
SG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("SG")) %>%
  mutate(SG_prop_ptoid_attack = SG_Platy/(SG_Platy + SG_SG.larv)) %>% # omitting exit hole because of potential overlap with other ectoparasitoids. Need to probably resolve this by not using "unk" gall.id. SG_Platy was the only species reared from SG
  select(Genotype, SG_Platy, SG_SG.larv, SG_prop_ptoid_attack)

SG_for_chi <- SG_geno_parasitism %>%
  filter(SG_Platy > 0 & SG_SG.larv > 0) %>% # unable to run the anlaysis unless I permit there to be zeros in at least one column
  select(SG_Platy, SG_SG.larv)

chisq.test(SG_for_chi) # I don't know how much I trust this...appears to be driven by one data point.
chisq.test(SG_for_chi, simulate.p.value = TRUE, B = 10000)

# Cecidomyiidae sp. A (aSG)
aSG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("aSG")) %>%
  mutate(aSG_ptoid_attack = aSG_Tory.fem + aSG_Tory.mal, aSG_prop_ptoid_attack = aSG_ptoid_attack/(aSG_ptoid_attack + aSG_aSG.larv)) %>% 
  select(Genotype, aSG_ptoid_attack, aSG_aSG.larv, aSG_prop_ptoid_attack)

aSG_for_chi <- aSG_geno_parasitism %>%
  filter(aSG_ptoid_attack > 0 & aSG_aSG.larv > 0) %>% # 
  select(aSG_ptoid_attack, aSG_aSG.larv)

chisq.test(aSG_for_chi) # qualitatively the same results no matter how I decide to retain the data.
chisq.test(aSG_for_chi, simulate.p.value = TRUE, B = 10000)

# Pontania californica
rsLG_geno_parasitism <- gall_net_melt %>%
  dcast(Genotype ~ gall.sp + gall_contents, sum) %>%
  select(Genotype, starts_with("rsLG")) %>%
  mutate(rsLG_surv = rsLG_Pont.ad + rsLG_Pont.prep, rsLG_ptoid_attack = rsLG_Eury.fem + rsLG_Eury.mal + rsLG_Lathro.fem + rsLG_Lathro.mal, rsLG_prop_ptoid_attack = rsLG_ptoid_attack/(rsLG_surv + rsLG_ptoid_attack)) %>%
  select(rsLG_ptoid_attack, rsLG_surv, rsLG_prop_ptoid_attack)

rsLG_for_chi <- rsLG_geno_parasitism %>%
  filter(rsLG_ptoid_attack > 0 & rsLG_surv > 0) %>% # 
  select(rsLG_ptoid_attack, rsLG_surv)

chisq.test(rsLG_for_chi) # qualitatively the same results no matter how I decide to retain the data.
chisq.test(rsLG_for_chi, simulate.p.value = TRUE, B = 10000)

