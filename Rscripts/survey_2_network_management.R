##############################################################################
# Description: This script organizes and manages data collected in the fall of 
#              2012 in the willow garden at HBNWR in Eureka, California.
#
# Manager: Matt Barbour
# Email: barbour@zoology.ubc.ca
##############################################################################

# load libraries
library(reshape)
library(reshape2)
library(dplyr)

#### Data management

# Shoot count estimates as well as plant positions where zero galls were collected
shoot.countEst <- read.csv("~/Documents/Genotype_Networks/data/survey_2_stem_diams_shootEsts.csv")
shoot.count_df <- tbl_df(shoot.countEst)
shoot.count_df <- shoot.count_df %>%
  select(plant.position = Plant.Position, Galls.found, shootEst.all, shootEst.no18)

# upload plant position data with details on plant identity to eventually merge with gall network data
plant.position.info <- read.csv("~/Documents/Genotype_Willow_Community/datasets_&_Rscripts/Willow Garden Positions.csv")
plant.position.info <- tbl_df(plant.position.info)
plant.position.info <- plant.position.info %>%
  select(Genotype, Gender, Row = Row.., plant.position = Plant.Position)

# merge plant position info with shoot estimate data
plant.info_df <- left_join(shoot.count_df, plant.position.info) %>%
  select(Gender, Genotype, Row, plant.position:shootEst.no18)

# upload gall network data (survey #2)
Gall_Network_Data_Survey_2_2012 <- read.csv("~/Documents/Genotype_Networks/data/Gall_Network_Data_Survey_2_raw.csv", skip=1, stringsAsFactors = FALSE, strip.white = TRUE)
Gall_Network_Data_Survey_2_2012 <- tbl_df(Gall_Network_Data_Survey_2_2012)
Gall_Network_Data_Survey_2_2012[which(Gall_Network_Data_Survey_2_2012$plant.position == 49.8),"plant.position"] <- rep(167, 7) # 49.8 label came from the stem diameter for this plant, but we also knew it was genotype "O". I also noticed I was missing plant 167 from the dataset and plant was Genotype "O" and had a stem diameter of 49.8, therefore I input 167 as the plant position.

# replaces all NAs in gall content matrix with zeros.  These are biologically meaningful zeros and do not represent missing data.
g.contents_matrix <- as.matrix(select(Gall_Network_Data_Survey_2_2012, vLG.pupa:exit.hole))
head(g.contents_matrix)
g.contents_matrix[is.na(g.contents_matrix)] <- 0 

# merge the updated gall content matrix back with the original survey data
gall_net <- cbind.data.frame(Gall_Network_Data_Survey_2_2012[ ,c("survey","plant.position","gall.id","g.id.unk","gall.sp","g.sp.unsure","g.height","g.width","g.3meas","point.count", "exit.size")], g.contents_matrix) # excludes  "data.page", "found.in.bag", "need.to.review" and "notes" from original data
gall_net <- tbl_df(gall_net)

# merge the gall network data and plant position information
#gall_net <- left_join(gall_net, plant.position.info)
#gall_net <- select(gall_net, survey, Gender, Genotype, Row, plant.position:exit.hole) # reorder

#### Examine data for errors and make appropriate corrections. Since I never found any "twG" larva, I'm ommitting this as a possible gall species (although it may related in some way to aSG).

gall_net[which(gall_net$g.id.unk == "unk"), ] # represents all galls in which we could not directly link the parasitoid to the exact gall for g.height, g.width, g.3meas, point.count, and exit hole. Therefore, it should not count toward gall abundance, but the parasitoid information is useful.


# gall.sp identification data errors and adjustment
table(gall_net$gall.sp) # examine approximate number of galls collected for each species. Not entirely accurate because there may be multiple larva within SG and vLG galls and it doesn't adjust for g.id.unk.

gall_net[which(gall_net$gall.sp == "VLG"),c("gall.sp")] <- c("vLG","vLG") # double checked with original data

gall_net[which(gall_net$gall.sp == ""),c("gall.sp")] <- "vLG" # notes indicate that the parasitoids were found in the vial, and only vLG were in the vial. However, these parasitoids could not be tied to a particular gall, which is why there were no height measurements.

gall_net[which(gall_net$gall.sp == "twG" & gall_net$plant.position == "257" ), "gall.sp"] <- "aSG" # after looking at the picture, I have confirmed that was an aSG.

gall_net[which(gall_net$g.sp.unsure == 1),]  # find all galls with "unsure" labels
gall_net[which(gall_net$gall.sp == "aSG/twG"),"gall.sp"] <- "aSG" # may need to double check this. Should be okay though, since I no longer am considering twG as a possible gall sp.
gall_net[which(gall_net$gall.sp == "aSG" & gall_net$g.sp.unsure == 1), "g.sp.unsure"] <- 0 #see above note

gall_net[which(gall_net$gall.sp == "SG/PG"),c("gall.sp")] <- "SG" # checked noted and picture on my phone and confirmed that this was a Stem Gall.
gall_net[which(gall_net$gall.sp == "SG" & gall_net$g.sp.unsure == 1), "g.sp.unsure"] <- 0 # see above note

gall_net[which(gall_net$gall.sp == "SG2"),"gall.sp"] <- rep("SG",7) # I belive SG2 is simply a deformed SG so I've decided to lump these together. This conclusion is further supported by the fact that SG and SG2 larva looked identical.

gall_net <- gall_net[-which(gall_net$gall.sp == "PG/bud gall"), ] # removed these from the dataset (were plant position #63), because I'm not confident about their identification.
gall_net <- gall_net[-which(gall_net$gall.sp == "PG"), ] # removed this sample from the dataset, because I'm not confident about its identification
gall_net <- gall_net[-which(gall_net$gall.sp == "twG"),] # removed from dataset, because I'm not confident that this was actually a gall species.
gall_net <- gall_net[-which(gall_net$gall.sp == "rsLG" & gall_net$g.sp.unsure == 1), ] # unusual case and I'm not sure about rsLG identification, so I removed it from the data set

gall_net[which(gall_net$g.sp.unsure == 1),] # no longer unsure gall.sp

table(gall_net$gall.sp) # only 5 gall.sp present


# gall-ptoid interaction data errors and adjustments
gall_net <- gall_net[-which(gall_net$Pont.ad == 1 & gall_net$gall.sp == "rG"),] # removed this datapoint, because it was a sawfly adult found in the vial, and not explicitly linked as emerging from an rG (g.id.unk == unk).
gall_net[which(gall_net$gall.sp == "rG" & gall_net$Pont.prep == 1),"Pont.prep"] <- 0 # in notes, larva was "unknown", so I'm unsure why it was included as "Pont.prep" here.
gall_net <- gall_net[-which(gall_net$gall.sp == "rsLG" & gall_net$Tory.mal == 1),] # I'm removing this data point and classifying this as a data error for the following reasons. This parasitoid was found within a vial (i.e. g.id.unk = unk) that only had rsLG, but both of these galls were "Solid rsLG with no larva". Therefore, the parasitoid could not have been from either of these galls. That, and we have no other recordings of a Torymus coming from rsLG, plus no evidence of rsLG being associated with Torymus in the literature (Caltagirone 1964)
gall_net[which(gall_net$gall.sp == "rsLG" & gall_net$Eulo.fem == 1),"gall.sp"] <- "rG" # although paper copy of data suggests this is "rsLG", this parasitoid was put in the "rG" collection bag for Eulophid females. This piece of evidence as well as the fact that this is the only occassion where a Eulophid was associated with rsLG, and no evidence from the literature of Eulophids attacking rsLG (Caltagirone 1964), compels me to switch the gall.sp id to "rG".

gall_net[which(gall_net$gall.sp == "rsLG" & gall_net$Platy == 1),"rsLG.solid"] <- 1 # although paper copy of data suggests this was sp. 12 ("Platy"), the note says "no signs of life", plus the sp.12 category was right next to the "nothing" category. Additionally, this gall was a very small "rsLG" which are typically solid galls. Furthermore, this is the only "rsLG-Platy" association in the entire data. This evidence compels me to change this to "rsLG.solid = 1" for this data entry. (PERHAPS EVEN SOLID...)
gall_net[which(gall_net$gall.sp == "rsLG" & gall_net$Platy == 1),"Platy"] <- 0 # see note above.

gall_net[which(gall_net$gall.sp == "vLG" & gall_net$Eury.mal == 1),"Tory.mal"] <- 1 # data entry error. Original data specify this is a "Tory.mal"
gall_net[which(gall_net$gall.sp == "vLG" & gall_net$Eury.mal == 1),"Eury.mal"] <- 0 #see note above.

gall_net[which(gall_net$gall.sp == "rG" & gall_net$Mesopol == 1),] # this datapoint seems legit, seeing as how a new bag was explicitly created for this association. However, it is the only sample of this interaction.

gall_net[which(gall_net$gall.sp == "rG" & gall_net$plant.position == 626 & gall_net$gall.id == 13),"exit.hole" ] <- 0 # after referencing notes, I don't know if exit hole was the accurate label, because the notes said the "gall was already opened when first observed). I decided to label this simply as "unsure
gall_net[which(gall_net$gall.sp == "rG" & gall_net$plant.position == 626 & gall_net$gall.id == 13), "unsure"] <- 1 # see note above

gall_net[which(gall_net$gall.sp == "aSG" & gall_net$plant.position == 374 & gall_net$gall.id %in% c(5.1,5.2)), "exit.hole"] <- c(0,0) # I'm unsure about these labels, primarily because one of the exit hole sizes in the hard copy of data seemed rather large (big than a typical ectoparasitoid exit hole) and I have personally never seen an exit hole on an aSG. Therefore, I've decided to change these to "unsure" labels.
gall_net[which(gall_net$gall.sp == "aSG" & gall_net$plant.position == 374 & gall_net$gall.id %in% c(5.1,5.2)), "unsure"] <- c(1,1)

gall_net[which(gall_net$gall.sp == "SG" & gall_net$plant.position == 63 & gall_net$gall.id == 4), ] # I saw what looks like Josh's picture of this gall with an exit hole and it seems to be pretty clearly an exit hole from an unknown parasitoid, therefore I've decided to keep it.

# removing galls collected from ground or survey 1 or non-target branches. I'm also removing all galls on plant position 656, because all of these galls were missed from survey 1.
gall_net <- gall_net[-which(gall_net$plant.position == 656), ] # all of these galls were missed on survey 1, and therefore shouldn't count as part of survey 2.
gall_net <- gall_net[-which(gall_net$plant.position == 495 & gall_net$gall.id == 1.2), ] # collected from non-target branch
gall_net <- gall_net[-which(gall_net$plant.position == 495 & gall_net$gall.id == 1.1), ] # collected from non-target branch
gall_net <- gall_net[-which(gall_net$plant.position == 495 & gall_net$g.id.unk == "unk"), ] # associated with galls collected from non-target branch

gall_net <- gall_net[-which(gall_net$plant.position == 494), ] # notes indicate that both galls were collected from a non-target branch so I should not consider them in these density estimates. Therefore, galls found should be changed to "zero" for the shoot estimate data.
gall_net <- gall_net[-which(gall_net$plant.position == 73), ] # removing from this data set because it was collected from the ground.
gall_net <- gall_net[-which(gall_net$plant.position == 69),] # removing from this data set because it was collected from the ground
gall_net <- gall_net[-which(gall_net$plant.position == 71 & gall_net$gall.id %in% c(1.1,1.2)),] # removing from this data set because it was collected from the ground
gall_net <- gall_net[-which(gall_net$plant.position == 75 & gall_net$g.height == 6.31),]# removing from this data set because it was collected from the ground
gall_net <- gall_net[-which(gall_net$plant.position == 65 & gall_net$gall.id %in% c("B,1")),] # only found one gall in dataset, although field notes suggested two galls were collected from the ground under plant #65.
gall_net <- gall_net[-which(gall_net$plant.position == 402),] # this gall was actually collected from a survey #1, so I removed it from this dataset.


# more gall-ptoid interaction data errors and adjustments
gall_net[which(gall_net$unk.larv == 1 & gall_net$plant.position == 365), "diff.or.larv"] <- 1 # confirmed that this should be classified as "diff.or.larv"
gall_net[which(gall_net$unk.larv == 1 & gall_net$plant.position == 365), "unk.larv"] <- 0 # switched to diff.or.larv (note above)

gall_net[which(gall_net$unk.larv == 1 & gall_net$plant.position == 509), "rG.wh.larv"] <- 3 # added unknown larv to this designation, after I double checked it under the microscope
gall_net[which(gall_net$unk.larv == 1 & gall_net$plant.position == 509), "unk.larv"] <- 0 # see note above

gall_net[which(gall_net$unk.larv == 1 & gall_net$plant.position == 270), "vLG.pupa"] <- 1 # although it wasn't found in a white papery cocoon, it appeared to be a fully grown Iteomyia larva.
gall_net[which(gall_net$unk.larv == 1 & gall_net$plant.position == 270), "unk.larv"] <- 0 # see note above

gall_net[which(gall_net$unsure == 1 & gall_net$plant.position == 59), "vLG.pupa"] <- 3 # data entry error. Also, moved "unsure" to vLG pupa based off note that I found a cocoon of the vLG, but didn't find the larva inside...
gall_net[which(gall_net$unsure == 1 & gall_net$plant.position == 59), "unsure"] <- 0 

gall_net[which(gall_net$unsure == 1 & gall_net$plant.position == 370), ] # unsure label is accurate. Forgot to check for exit hole, so I don't know if the missing larva was a result of inquiline damage, parasitism, or something else.
gall_net[which(gall_net$unsure == 1 & gall_net$plant.position == 63), ] # "unsure" label is accurate. I was not able to find a picture or the contents to double check, so I'm maintaining the "unsure" label.
gall_net[which(gall_net$nothing == 1 & gall_net$exit.hole == 1), "nothing"] <- 0 # changed it to zero, so I would only have exit hole marked down (because if there is an exit hole, then there is nothing inside).
gall_net[which(gall_net$exit.hole == 1 & gall_net$Eulo.fem == 1), "Eulo.fem" ] <- 0 # data entry error
gall_net[which(gall_net$exit.hole == 1 & gall_net$Eulo.mal == 1 & gall_net$plant.position == 62), "Eulo.mal"] <- 0 # data entry error
gall_net[which(gall_net$exit.hole == 1 & gall_net$vLG.pupa == 1),] # this scenario is okay because there were multiple larva in the gall, so its possible to have both a vLG.pupa and an exit hole.

gall_net[which(gall_net$plant.position == 461 & gall_net$gall.id == 17.1),"Eulo.mal" ] <- 0 # moved to gall.id 17 since we couldn't actually assign this parasitoid to an exact gall.
gall_net[which(gall_net$plant.position == 461 & gall_net$gall.id == 17), "Eulo.mal"] <- 1 # see above note

gall_net[which(gall_net$unk.larv == 1 & gall_net$plant.position == 535), "unk.larv"] <- 0 # hard copy of data page missing...But, I looked at the vial and it appeared to be a sawfly larva, suggesting that it was associated with the outside of the 'rG' gall. Therefore I changed this to a zero.

gall_net[which(gall_net$moth.ad1 == 1 & gall_net$exit.hole == 1), "exit.hole" ] <- 0 # confirmed with hard copy of data that this "exit hole" was actually likely a result of moth adult 1 and would have been better categorized as inquiline damage. However, there is no need to also put inquiline damage because this would be redundant.

gall_net[which(gall_net$moth.larv == 1 & gall_net$plant.position == 471 & gall_net$gall.id == 10 & gall_net$exit.hole == 1), ] # note that this overlap is okay. Looking at the notes, the exit hole appears to have been due to parasitoid emergence and not confusion between inquiline damage and an exit hole.

gall_net[which(gall_net$plant.position == 325 & gall_net$gall.id == 4), "nothing"] <- 0 # may have supposed to have been exit hole, but it is really unclear to me why this was labeled as such, therefore I simply changed it to "unsure"
gall_net[which(gall_net$plant.position == 325 & gall_net$gall.id == 4), "unsure"] <- 1 # see note above

gall_net[which(gall_net$plant.position == 370 & gall_net$gall.id == 1), "unsure" ] <- 0 # have inq.dam was sufficient

gall_net[which(gall_net$plant.position == 75 & gall_net$gall.id == 7), "nothing"] <- 0 # data entry error

# mystery...
gall_net[which(gall_net$plant.position == 77),] #77 is a bit of a mystery (field notes suggest no galls, but there is a gall collection and there is no note explicitly stating that these were from survey #1). Right now, I've decided to keep it in the dataset. Notably though, I did sample a big branch and it does have the galls typically associated with Genotype T...



# Notes below are for more abiguous categories. I believe I should initially include them to see if they vary systematically at all among the genotypes. If they do not, then I think I am justified in removing them from the analysis.
gall_net[which(gall_net$rsLG.solid == 1), ] # looks good
gall_net[which(gall_net$moth.larv == 1), ] # moth.larv does not appear to always be associated with parasitism so I think I should remove it. Note that "unsure" column may need to be filled in for instances where no other mortality/survival can be attributed. Especially, if this is "rG", it should maybe be nothing CHECK NOTES
gall_net[which(gall_net$unk.ptoid == 1), ] # may be able to determine whether this was a Torymus (but unknown male/female)
gall_net[which(gall_net$Pont.larv == 1), ] # may be due to being paralyzed by an ectoparasitoid or not permitted to fully develop.
gall_net[which(gall_net$exit.hole == 1), ] # this is a difficult category, because this information may only be useful for identifying parasitism or attack from vLG. For example, rsLG exit holes may have been due to prepupa chewing a large hole. Exit holes can't be accurately identified on rG galls. Also, I think Josh may have confused inquiline damage with "exit holes", so even this may be difficult. I may be able to sort some of these out, but I would need to match these up with the appropriate exit hole size and I won't be able to do this for all of them.  
gall_net[which(gall_net$diff.or.larv == 1), ] # need to decide what to do here...
gall_net[which(gall_net$nothing == 1), ] # need to sort out the meaning behind these. Its possible for 'rG' that some of the galls collected were "old galls" and therefore would have already been exited? No, I am pretty sure I restricted my sampling to fresh galls...
# 
gall_net[which(gall_net$moth.ad1 == 1), ] # it appears that "moth.ad1" was not consistently a mortality factor since I reared a Platy from this gall that only had 1 vLG larva (point.count = 1). Therefore, I'm no longer going to consider 'moth.ad1' in the dataset
gall_net[which(gall_net$incidental.sp.239 == 2), ] # although sp.239 is a parasitoid, it appears to just be using the habitat created by rG to pupate. Therefore, it should not be considered as a column of useful data for further analyses.
gall_net[which(gall_net$inq.dam == 1), ] # inq.dam may also associated with "unsure" for vLG and "nothing" for rsLG. 
gall_net[which(gall_net$moth.ad3 == 1),] # datapoint #906 suggests that moth.ad3 wasn't always a mortality agent since Pont.ad emerged. Consider removing this data.

# looking at some other aspects of the data 
colSums(select(filter(gall_net, exit.hole == 1), vLG.pupa:exit.hole)) # check to see where there is overlap between exit holes and other category designations. Some overlap with moth damage and exit holes
table(gall_net$gall.sp, gall_net$exit.hole)

colSums(select(gall_net, vLG.pupa:exit.hole)) # quick glance of the abundance of different categories

colSums(gall_net[which(gall_net$moth.ad1 > 0 | gall_net$moth.larv > 0 | gall_net$moth.ad3 > 0 | gall_net$inq.dam > 0), 15:48]) # it doesn't appear that inquline attack was always associated with mortality, therefore, I feel justified in removing these columns from the dataset to focus on gall-parasitoid associations.


# add a new gall.id variable so each gall will be accurately nested within a particular plant position. In retrospect, I could have probably done this easier to make it easier to check the data for errors.
gall_net <- gall_net %>%
  mutate(gall.id.nest = 1:dim(gall_net)[1], rG.larv = rG.wh.larv + rG.or.larv) %>% # combined rG.wh.larv and rG.or.larv because they appear to be the same species
  select(survey:gall.id, gall.id.nest, g.id.unk:vLG.pupa, rG.larv, SG.larv:exit.hole, -c(moth.ad1, moth.ad3, inq.dam, moth.larv, incidental.sp.239, unk.larv, twG.larv)) # removed moth damage for several reasons. First, they were not always a source of mortality for the galler when present. Second, they appear to make up a small portion of the dataset. Finally, this enables me to focus the analysis on gall-parasitoid associations. Also, after my data checks, there are no more "unk.larv", plus incidental.sp.239 was associated with the outside of the gall and not a mortality agent. twG.larv is no longer a category either as I converted these to all aSG (and never found a larva inside a twG for sure)

colSums(select(gall_net, vLG.pupa:exit.hole)) # quick glance of the abundance of different categories again

dim(shoot.count_df)[1] - length(table(gall_net$plant.position)) - dim(filter(shoot.count_df, Galls.found == 0))[1] # couldn't find galls associated with 7 of the trees from the original survey...
shoot.count.df
filter(shoot.count_df, Galls.found == 0)

# melt the data frame for easier management
gall_net_melt <- gall_net %>%
  select(plant.position, gall.id:gall.sp, g.height:exit.hole) %>%
  melt(id.vars = c("plant.position", "gall.id", "gall.id.nest", "g.id.unk", "g.height", "g.width", "g.3meas", "point.count", "exit.size", "gall.sp"), variable_name = "gall_contents") #"Gender","Genotype","Row",
gall_net_melt <- tbl_df(gall_net_melt)
gall_net_melt <- filter(gall_net_melt, value > 0 & plant.position < 9999 & plant.position != 310) # remove unecessary data. Note that this has more rows than the gall.net data which is likely because there were some galls with multiple types of contents. Also dropped the unidentified plant position (#9999). Also dropped plant position 310 because I found out that it was actually a duplicate of 510 (which was already in the data set), plus this tree was never surveyed...

cbind(names(table(gall_net_melt$plant.position)), names(table(gall_net$plant.position))) # looks like plant.position 521 was the only plant that had no "values" associated with it after removing selecting the subset of columns. Used to have an rsLG-moth.ad3 connection, but I removed these from the dataset. However, I'm going to make it so this one had no galls associated with it to be consistent with the treatment of other interactions.


# identify duplicates for gall id. Note that many of these are okay, since they are often associated with vLG which may contain multiple larva. Some of the moth related ones may end up being thrown out too. I need to make sure all of the rG ones are okay.
id.duplicates <- gall_net_melt %>%
  filter(duplicated(gall_net_melt$gall.id.nest) == TRUE) %>%
  select(plant.position, gall.id, gall.id.nest, point.count, g.id.unk, gall.sp, gall_contents, value)

duplicate_df = gall_net_melt[gall_net_melt$gall.id.nest %in% id.duplicates$gall.id.nest, ]
write.csv(sort_df(duplicate_df, "gall.id.nest"), "~/Documents/Genotype_Networks/data/gall_network_duplicates_to_check.csv") # all duplicates seem okay, because there are often cases where multiple specimens may come from the same galls. I have double checked all of these.

### Add plant position where no galls were collected, even though the tree was surveyed for galls.
plant.positions.no.galls <- shoot.count_df %>%
  filter(Galls.found == 0 & plant.position != 77) %>% # removed 77, because we do have galls linked to this, although it is a bit of a mystery...
  select(plant.position)
no_gall_df <- data.frame(plant.position = c(plant.positions.no.galls$plant.position,521),
                         gall.id = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         gall.id.nest = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         g.id.unk = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         g.height = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         g.width = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         g.3meas = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         point.count = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         exit.size = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         gall.sp = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         gall_contents = rep("NA", length(plant.positions.no.galls$plant.position)+1),
                         value = rep("NA", length(plant.positions.no.galls$plant.position)+1))

gall_net_melt_with_no_gall_data <- rbind(gall_net_melt, no_gall_df)


### Merge in data with estimated shoot counts and add Genotype, Gender, and row information
gall_net_melt_plant_info <- left_join(gall_net_melt_with_no_gall_data, plant.info_df) %>%
  select(Gender:shootEst.no18, plant.position:value)

#### Data has been checked and appears to be error free.

write.csv(gall_net_melt_plant_info,"~/Documents/Genotype_Networks/data/gall_network_data.csv")

