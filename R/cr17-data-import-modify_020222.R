####################
## Importing data ##
####################
bac_ko <- bind_rows(read.csv("base_data/bac_ko1.csv"), 
                    read.csv("base_data/bac_ko2.csv"))
abun <- read.csv("base_data/bin_abundance_table.csv")
metapath <- read.csv("base_data/bin_metapath.csv")
geochem <- read.csv("base_data/SubductCR_bac_sample_table2.csv")
bin_quality <- read.csv("base_data/bin_quality.csv")
Fe_genie_data <- read_csv("base_data/FeGenie-heatmap-data.csv")
auto_gene_abun <- read_csv("base_data/autotroph_gene_abundance.csv") %>% rename("gene_id"= "Reference Contig")
MAGQuality <- read_csv("base_data/MAGQualityAssessment.csv") %>%
  rename('bins' = "organism" )
#######################################
## Initial Modifications to the data ##
#######################################
#Create metapath
metapath <- unite(metapath, bin.taxa.name, c(1:8), sep=", ", remove=FALSE)

#Now I need to get rid of the extra commas
metapath$bin.taxa.name <- gsub(",", " ", metapath$bin.taxa.name)
metapath$bin.taxa.name <- trimws(metapath$bin.taxa.name, which = c("both"))
metapath$bin.taxa.name <- gsub("  ", ", ", metapath$bin.taxa.name)

metapath <- inner_join(bin_quality,
                       metapath, by = "bins") # Adding the bin quality to metapath data

bin_with_taxa <- subset(metapath,
                        select = c(bins, bin.taxa.name, kingdom, phylum, class, 
                                   order, family, genus, species)) # Subseting taxa names from metapath
abun <- inner_join(bin_with_taxa, 
                   abun, by = "bins") # Adding taxa names from metapath to abun


# Changing from site intials to site identifiers
abun <- abun %>% rename(
  "MTF" = MTF_, "BRF1" = BRF1_, "RSF" = RS_S9_, "SIS" = SIS_, "QHF2" = QH2_S7_,
  "CYF" = CY_S3_, "SLF" = SL_S11_, "BRS2" = BRS2_, "CYS" = CYS_S13_, "RVF" = RVF_, 
  "RSS" = RSS_, "TCF" = TC_S12_, "QNF" = QN_S8_, "EPS" = EPS_S14_, "EPF" = EP_S4_, 
  "SLS" = SLS_, "QNS" = QNS_, "BRF2" = BRF2_, "TCS" = TCS_, "FAS" = FAS_S15_, 
  "BQS" = BQS_, "ETS" = ETS_, "QHS2" = QHS2_, "SIF" = SI_S10_,"BQF" = BQF_, 
  "BRS1" = BRS1_, "ESF9" = ESF9_
)

# Removing sites that we are not analysizing: ALL_READS, PF, PBS, PFS, PGS, PG, and ARS
abun <- abun[setdiff(colnames(abun), c("BQ_S1_", "BR1_S2_", "PGS_", "PBS_", "ARS_", "ALL_READS", 
                                       "PG_S6_", "PFS_", "PF_S5_", "PLS_", "QHS1_S13_"))]
abun <- inner_join(bin_quality, abun, by = "bins") #Adding bin quality to abun data

##### Site Names for d


###############################################################################
## Figure 1A: Modifying data to constructed Dendrogram of sample sites/types ##
###############################################################################
abun_py <- abun
#Find the the highest abundance of the lowest abundant MAG and use anything less than that as a cut off
abun_py[,15:41][abun_py[,15:41] < .03] <- 0 #This value was made based on the low abun (bin_296 < .03)

cuttoff_site_hc <- as.dendrogram(
  hclust(d = dist(cor(abun_py[,-1:-14], method = "spearman")))) #Spearman rank coorelation gave a better clustering pattern than Person

#######################################################################
## Figure 1B: Modifying data to constructed Dendrogram of MAGs(bins) ##
#######################################################################
cuttoff_site_order <- order.dendrogram(cuttoff_site_hc) #Look at this output and arrange the columns based on it
cuttoff_sites <- abun_py[, c(1:14, (cuttoff_site_order + length(1:14)))] #Based on cuttoff_site_order, this puts the sites in order of Spearman Rank Correlation

#Finding spearman correlations between bins
cut_off_correlation <- cuttoff_sites # To preserve 'cuttoff_sites' object, I renamed it here for further analysis
row.names(cut_off_correlation) <- cut_off_correlation$bin.taxa.name # This designates the 'bin.taxa.name' column as the row names
cut_off_correlation <- cut_off_correlation[,-1:-14] # This removes taxonomic and bin assignments as these are not required for correlation analysis. We have already assigned 'bin.taxa.names' as row names
cut_off_correlation2 <- as.data.frame(t(cut_off_correlation)) # Here I am transposing the 'cut_off_correlation' so I can run a Spearman correlation between the MAGs
cut_off_hclust <- as.dendrogram(hclust(d = dist(cor(cut_off_correlation2, method = "spearman")))) # Constructing the dendrogram based on the above Spearman correlation

#######################################################################################
## Figure 1C: Modifying data to constructed Heat Map Figure based on the dendrograms ##
#######################################################################################
#Naming the spearman bin correlation order
cut_off_correlation.order <- order.dendrogram(cut_off_hclust)

cut_off_long <- gather(cuttoff_sites, site, abundance, 
                       "QHF2":"TCS", 
                       factor_key = TRUE) # Converting to long data


cut_off_long$bin.taxa.name  <- factor(x = cut_off_long$bin.taxa.name ,
                                      levels = cuttoff_sites$bin.taxa.name[cut_off_correlation.order], 
                                      ordered = TRUE) # This preserves the MAG orders for the figure

###############################################################################################
## Figure 1D-1E Modifing the data to identify MAGs with and without Carbon Fixation Pathways ##
###############################################################################################



  ##### Autotrophs ####
#Making a dataframe with autotrophs only
autotrophs <- metapath %>% subset(CBB.Cycle >= .6 | 
                                             (rTCA.Cycle >= 1 & TCA.Cycle >= .5) |
                                             Wood.Ljungdahl >= .6) %>% 
  mutate(.,"CBB.present" = case_when(CBB.Cycle >= .6 ~ 0, TRUE ~ 0),# CBB.Cycle is set to zero for below filtration
         "rTCA.present" = case_when(rTCA.Cycle == 1 ~ 1, TRUE ~ 0),
         "Wood.Ljungdahl.present" = case_when(Wood.Ljungdahl >= .6 ~ 0, TRUE ~ 0))# WL is set to zero for below filtration

# I need to identify MAGs that have both phosphoribulokinase and RuBisCo
#KO for CBB
cbb_ko <- autotrophs %>% subset(., CBB.present == 0) %>% 
  subset(., select = c(bins)) %>% 
  inner_join(., bac_ko, by = "bins")
cbb_phospho <- cbb_ko[cbb_ko$ko == "K00855",] # K00855 is the KO identifier for phosphoribulokinase
cbb_phospho_bins <- data.frame(unique(cbb_phospho$bins, ' ')) %>% 
  rename(., bins = unique.cbb_phospho.bins......)
cbb_phospho_ko <- cbb_phospho_bins %>% inner_join(., bac_ko, by = "bins") #Get all KO numbers of those with phosphoribulokinase to see if we can find Rubisco

cbb_rubisco_01_ko <- cbb_phospho_ko[cbb_phospho_ko$ko == 'K01601',] # Here I used the bins with phosphoribulokinase and looked to see if they had the large subunit, 'K01601'
cbb_rubisco_01_ko_unique <- data.frame(unique(cbb_rubisco_01_ko$bins, ' ')) # All of them have the large subunit.
cbb_rubisco_01_ko_unique <- cbb_rubisco_01_ko_unique %>% rename(., bins = unique.cbb_rubisco_01_ko.bins......)

cbb_rubisco_02_ko <- cbb_phospho_ko[cbb_phospho_ko$ko == 'K01602',] # Here I used the bins with phosphoribulokinase and looked to see if they had the small subunit, 'K01602'
cbb_rubisco_02_ko_unique <- data.frame(unique(cbb_rubisco_02_ko$bins, ' ')) #31 bins with RubisCo Small subunit
cbb_rubisco_02_ko_unique <- cbb_rubisco_02_ko_unique %>% rename(., bins = unique.cbb_rubisco_02_ko.bins......)

autotrophs$CBB.present[match(cbb_phospho_bins$bins, autotrophs$bins)]<-1 #Using the ko analysis above, I was able to find the bins with CBB



# I need to identify MAGs that have Key enzymes for WL
#KO for WL
wl_ko <- autotrophs %>% subset(., Wood.Ljungdahl.present == 0) %>%
  subset(., select = c(bins)) %>% 
  inner_join(., bac_ko, by = "bins")
wl_co_methly_acetyl_CoA_syn <- wl_ko[wl_ko$ko == "K00192" | wl_ko$ko == "K14138",] # These are the 2 catalytic subunits for acetyl-CoA synthase

wl_co_methly_acetyl_CoA_syn_bins <- data.frame(unique(wl_co_methly_acetyl_CoA_syn$bins, ' ')) #78 bins with CO-methylating acetyl-CoA Synthase 
wl_co_methly_acetyl_CoA_syn_bins <- wl_co_methly_acetyl_CoA_syn_bins %>% 
  rename(., bins = unique.wl_co_methly_acetyl_CoA_syn.bins......)#Start right here

wl_co_methly_acetyl_CoA_syn_ko <- wl_co_methly_acetyl_CoA_syn_bins %>%
  inner_join(., bac_ko, by = "bins")

wl_CO_dehydrogenase_anaero_ko <- wl_co_methly_acetyl_CoA_syn_ko[wl_co_methly_acetyl_CoA_syn_ko$ko == 'K00198',]
wl_CO_dehydrogenase_anaero_ko_unique <- data.frame(unique(wl_CO_dehydrogenase_anaero_ko$bins, ' ')) #67 bins with anaerobic CO dehydrogenase Large subunit

wl_CO_dehydrogenase_aero_ko <- wl_co_methly_acetyl_CoA_syn_ko[wl_co_methly_acetyl_CoA_syn_ko$ko == 'K03520',]
wl_CO_dehydrogenase_aero_ko_unique <- data.frame(unique(wl_CO_dehydrogenase_aero_ko$bins, ' ')) #7 bins with aerobic carbon-monoxide dehydrogenase large subunit
wl_CO_dehydrogenase_aero_ko_unique <- wl_CO_dehydrogenase_aero_ko_unique %>% 
  rename(., bins = unique.wl_CO_dehydrogenase_aero_ko.bins......)

#Combined wl_anaero and aero
wl_CO_dehydrogenase_anaeroandaero_ko <- wl_co_methly_acetyl_CoA_syn_ko[wl_co_methly_acetyl_CoA_syn_ko$ko == 'K03520' | wl_co_methly_acetyl_CoA_syn_ko$ko == 'K00198',]

wl_CO_dehydrogenase_anaeroandaero_ko_unique <- data.frame(unique(wl_CO_dehydrogenase_anaeroandaero_ko$bins, ' '))
wl_CO_dehydrogenase_anaeroandaero_ko_unique <- wl_CO_dehydrogenase_anaeroandaero_ko_unique %>% 
  rename(., bins = unique.wl_CO_dehydrogenase_anaeroandaero_ko.bins......)

autotrophs$Wood.Ljungdahl.present[match(wl_CO_dehydrogenase_anaeroandaero_ko_unique$bins, autotrophs$bins)]<-1 #Using the ko analysis above, I was able to find the bins with WL
#write_csv(autotrophs, "output_data/data/autotrophs_010821.csv")



  ##### Heterotrophs ####
#Creating the start of heterotrophs
heterotrophs <- subset(metapath, CBB.Cycle < .6 & (rTCA.Cycle < 1 | TCA.Cycle < .5) & Wood.Ljungdahl < .6)
heterotrophs <- heterotrophs %>% 
  mutate(.,"CBB.present" = case_when(CBB.Cycle < .6 ~ 0, TRUE ~ 0),
         "rTCA.present" = case_when(rTCA.Cycle < 1 ~ 0, TRUE ~ 0),
         "Wood.Ljungdahl.present" = case_when(Wood.Ljungdahl < .6 ~ 0, TRUE ~ 0))

  ##### combining the pre-trophy dataframes #####
cfp <- bind_rows(autotrophs, heterotrophs)
true_autotrophs <- cfp %>% 
  subset(CBB.present ==1 | rTCA.present == 1 | Wood.Ljungdahl.present == 1)
true_autotrophs <- cbind(Autotrophy = 1, true_autotrophs)
true_autotrophs <- true_autotrophs %>%
  select(-Autotrophy,Autotrophy)


true_heterotrophs <- cfp %>% 
  subset(CBB.present ==0 & rTCA.present == 0 & Wood.Ljungdahl.present == 0)
true_heterotrophs <- cbind(Heterotrophy = 1, true_heterotrophs)
true_heterotrophs <- true_heterotrophs %>% 
  select(-Heterotrophy ,Heterotrophy)

true_cfp <- bind_rows(true_autotrophs, true_heterotrophs)
true_cfp[is.na(true_cfp)] = 0

# Add Fe genie predictions to dataframe
Fe_genie_data <- column_to_rownames(Fe_genie_data, var = "X")
Fe_genie_data_flipped <- as.data.frame(t(Fe_genie_data))
Fe_genie_data_flipped  <- rownames_to_column(Fe_genie_data_flipped, var = "bins ")
Fe_genie_data_flipped <- as.data.frame(Fe_genie_data_flipped)
Fe_genie_data_flipped <- Fe_genie_data_flipped %>% rename("bins" = `bins `)
cfp_genie <- left_join(true_cfp, Fe_genie_data_flipped, by = "bins")
cfp_genie[,"iron_oxidation"]=sapply(cfp_genie[,"iron_oxidation"],function(x) replace(x,x>0,1))
cfp_genie[,"iron_reduction"]=sapply(cfp_genie[,"iron_reduction"],function(x) replace(x,x>0,1))
cfp_genie$bin.taxa.name <- factor(x = cfp_genie$bin.taxa.name ,
                                  levels = cuttoff_sites$bin.taxa.name[cut_off_correlation.order], 
                                  ordered = TRUE)
cfp_long <- gather(cfp_genie, trophic, presence, "Autotrophy":"Heterotrophy", factor_key = TRUE)

##################################################
## Figure 1E: Carbon Fixation Prediction Figure. #
##################################################
# Carbon Fixation prediction Map
#This is the metabolic map based on the above bin organization. This will likely be the final

#Add in CBB groups based on RuBisCo Form to cfpgenie
cbb_FI <- c("bin_3", "bin_252", "bin_224", "bin_71", "bin_145", #Bin 3 also has two unknown forms of rubisco, but we know that form I is involved in CBB, so we will focus on that one.
            "bin_132", "bin_265", "bin_322", "bin_280",
            "bin_194", "bin_381", "bin_365", "bin_207", 
            "bin_379","bin_387", "bin_20", "bin_130", 
            "bin_16", "bin_12", "bin_169", "bin_365",  "bin_264",
            "bin_75", "bin_198", "bin_143", "bin_301")

cbb_FII <- c("bin_261", "bin_5", "bin_100", "bin_136",
             "bin_14", "bin_161", "bin_308", "bin_87", 
             "bin_110")

cbb_FI_FII <- c("bin_153", "bin_254", "bin_49")

cbb_unknown <- c("bin_386", "bin_263", "bin_251", "bin_55")

cfp_genie <- cfp_genie %>% 
  mutate(.,"CBB_FormI" = case_when(bins %in% cbb_FI ~ 1, TRUE ~ 0),
         "CBB_FormII" = case_when(bins %in% cbb_FII ~ 1, TRUE ~ 0),
         "CBB_FormI_and_II" = case_when(bins %in% cbb_FI_FII ~ 1, TRUE ~ 0),
         "CBB_Form_Unknown" = case_when(bins %in% cbb_unknown ~ 1, TRUE ~ 0))



cfpathways_long <- gather(cfp_genie, carbon.fixation.pathway, carbon.pathway.present, 
                          "CBB.present":"Wood.Ljungdahl.present", factor_key = TRUE)
cfpathways_long$cfp.color <- ifelse(cfpathways_long$carbon.pathway.present < 1, 0, 1)
##################################################################
## Figure 1F-1G: Electron Donor and Acceptor Predictions Figure ##
##################################################################
#####'Hydrogen Oxidation'as donor ####
clean_metapath <- cfp_genie %>% mutate(.,"NiFe.hydrogenase_present" = 
                                         case_when(NiFe.hydrogenase > 0  ~ 1, TRUE ~ 0), 
                                       "NiFe.hydrogenase.Hyd.1_present" = 
                                         case_when(NiFe.hydrogenase.Hyd.1 > 0 ~ 0, TRUE ~ 0),
                                       "hydrogen.quinone.oxidoreductase_present" = 
                                         case_when(hydrogen.quinone.oxidoreductase > 0 ~ 1, TRUE ~ 0),
                                       "NAD.reducing.hydrogenase_present" = 
                                         case_when(NAD.reducing.hydrogenase > 0 ~ 0, TRUE ~ 0),
                                       "NADP.reducing.hydrogenase_present" = 
                                         case_when(NAD.reducing.hydrogenase > 0 ~ 0, TRUE ~ 0))
#Finding those with Hyd 1 Hydrogenase
NiFe.hydrogenase.Hyd.1_target <- c("K06282", "K06281")
NiFe.hydrogenase.Hyd.1_present_ko <- filter(bac_ko, ko %in% NiFe.hydrogenase.Hyd.1_target)
NiFe.hydrogenase.Hyd.1_present_ko.1 <- NiFe.hydrogenase.Hyd.1_present_ko %>% mutate(yesno = 1) %>% 
  distinct %>% spread(ko, yesno, fill = 0) 
NiFe.hydrogenase.Hyd.1_present_ko.2 <- NiFe.hydrogenase.Hyd.1_present_ko.1%>% 
  subset(., K06282 > 0) %>% subset(., select = c(bins, K06282))
NiFe.hydrogenase.Hyd.1_present_ko.3 <- NiFe.hydrogenase.Hyd.1_present_ko.1%>%
  subset(.,K06281 > 0) %>% subset(., select = c(bins, K06281))
NiFe.hydrogenase.Hyd.1_present_ko.4 <- NiFe.hydrogenase.Hyd.1_present_ko.3 %>% 
  inner_join(NiFe.hydrogenase.Hyd.1_present_ko.2, ., by = "bins")


NiFe.hydrogenase.Hyd.1_present_ko_bins <- data.frame(unique(NiFe.hydrogenase.Hyd.1_present_ko.4$bins, ' '))
NiFe.hydrogenase.Hyd.1_present_ko_bins <- NiFe.hydrogenase.Hyd.1_present_ko_bins %>% 
  rename(., bins = unique.NiFe.hydrogenase.Hyd.1_present_ko.4.bins......)
clean_metapath$NiFe.hydrogenase.Hyd.1_present[match(NiFe.hydrogenase.Hyd.1_present_ko_bins$bins, clean_metapath$bins)]<-1


#Fingind those with NAD reducing hydrogenase
NAD.reducing.hydrogenase_target <- c("K00436", "K18007")
NAD.reducing.hydrogenase_present_ko <- filter(bac_ko, ko %in% NAD.reducing.hydrogenase_target)
NAD.reducing.hydrogenase_present_ko.1 <- NAD.reducing.hydrogenase_present_ko %>% 
  mutate(yesno = 1) %>% 
  distinct %>% spread(ko, yesno, fill = 0)

NAD.reducing.hydrogenase_present_ko.2 <- NAD.reducing.hydrogenase_present_ko.1%>% 
  subset(., K00436 > 0) %>% subset(., select = c(bins, K00436))
NAD.reducing.hydrogenase_present_ko.3 <- NAD.reducing.hydrogenase_present_ko.1%>% 
  subset(.,K18007 > 0) %>% subset(., select = c(bins, K18007))
NAD.reducing.hydrogenase_present_ko.4 <- NAD.reducing.hydrogenase_present_ko.3 %>% 
  inner_join(NAD.reducing.hydrogenase_present_ko.2, ., by = "bins")


NAD.reducing.hydrogenase_present_ko_bins <- data.frame(unique(NAD.reducing.hydrogenase_present_ko.4$bins, ' '))
NAD.reducing.hydrogenase_present_ko_bins <- NAD.reducing.hydrogenase_present_ko_bins %>% 
  rename(., bins = unique.NAD.reducing.hydrogenase_present_ko.4.bins......)
clean_metapath$NAD.reducing.hydrogenase_present[match(NAD.reducing.hydrogenase_present_ko_bins$bins, clean_metapath$bins)]<-1

#Finding those with NADP reducing hydrogenase. This hydrogenase is a FeFe hydrogenase. Usually FeFe hydrogenases are associated with 'Hydrogen production and not oxidation. However, this is an exception. NADP reducing hydrogenase is capable of both produciton and oxidation of H2
NADP.reducing.hydrogenase_target <- c("K18332")
NADP.reducing.hydrogenase_present_ko <- filter(bac_ko, ko %in% NADP.reducing.hydrogenase_target)
NADP.reducing.hydrogenase_present_ko_bins <- data.frame(unique(NADP.reducing.hydrogenase_present_ko$bins, ' '))
NADP.reducing.hydrogenase_present_ko_bins <- NADP.reducing.hydrogenase_present_ko_bins %>% 
  rename(., bins = unique.NADP.reducing.hydrogenase_present_ko.bins......)
clean_metapath$NADP.reducing.hydrogenase_present[match(NADP.reducing.hydrogenase_present_ko_bins$bins, clean_metapath$bins)]<-1

#Creating long plot for 'Hydrogen Oxidation' and sulfur cycle components. 

H2as_donor_long_test <- clean_metapath %>% 
  mutate(.,"Hydrogen Oxidation" = 
           case_when(NiFe.hydrogenase_present > 0 | 
                       NiFe.hydrogenase.Hyd.1_present > 0 |
                       hydrogen.quinone.oxidoreductase_present > 0 | 
                       NAD.reducing.hydrogenase_present > 0 |
                       NADP.reducing.hydrogenase_present > 0 ~ 1 , TRUE ~ 0)) %>%
  add_column("rdsrAB" = 0)

  #####rdsrABas Donor ####
#rdsr added to rdsrABas a donor. dsrAB can be a donor and accepter. I looked for dsr genes and then constructed a gene tree with oxidative and reductive representatives. From that I was able to identify the MAGs that oxidize and those that reduce SO4. 
#As of Dec 29, there is`Nitric Oxide Reduction`need to adjust this further
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_145"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_264"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_110"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_12"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_255"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_265"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_280"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_5"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_136"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_14"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_169"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_38"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_16"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_322"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_143"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_198"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_301"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_87"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_153"] <- 1
H2as_donor_long_test$rdsrAB[H2as_donor_long_test$bins == "bin_254"] <- 1

#Adding in parts of the nitrogen cycle. Also seperating out Micro and Areobic cytochromes

clean_metapath2 <- H2as_donor_long_test %>% 
  mutate(., "Nitrogen Oxidation" = case_when(ammonia.oxidation..amo.pmmo. > .33 ~1, TRUE ~0)) %>%
  mutate(., "Nitrogen Reduction*" = case_when(nitric.oxide.reduction > .33 | 
                                        nitrous.oxide.reduction > .33 | 
                                        nitrite.reduction > .33 | 
                                        dissim.nitrate.reduction > .33 ~ 1, TRUE ~ 0)) %>%
  rename(., "Iron Oxidation" = iron_oxidation)  %>% 
  rename(., "Iron Reduction" = iron_reduction) %>% 
  mutate(., Microaerobic = case_when(Cytochrome.c.oxidase..cbb3.type > .25 | 
                                       Cytochrome.bd.complex > .25~ 0, TRUE ~ 0),
         Aerobic = case_when(Cytochrome.c.oxidase > .25 | 
                               Cytochrome.o.ubiquinol.oxidase > .25 ~ 0,TRUE ~ 0)) %>% #Cytochromes are set to zero to so I can look for correct subunits
  mutate(., "Sulfur Oxidation*" = case_when(thiosulfate.oxidation >= .32 |
                                alt.thiosulfate.oxidation.tsdA >= .32 |
                                alt.thiosulfate.oxidation.doxAD >= .32 |
                                sulfide.oxidation > 0 |
                                rdsrAB > 0 ~ 1, TRUE ~ 0)) %>%
  mutate(., "Sulfur Reduction" = case_when(dissimilatory.sulfite.....sulfide > 0 ~ 1, TRUE ~ 0))


  ##### Identifying Cytochromes (O2 Reduction) ####
#Cytochrome.c.oxidase..cbb3.type: catalytic subunit = ccoN
cbb_cataylic_subunit_target <- c("K00404") #This is correct
cbb_cataylic_subunit_ko <- filter(bac_ko, ko %in% cbb_cataylic_subunit_target)
cbb_cataylic_subunit_ko_bins <- data.frame(unique(cbb_cataylic_subunit_ko$bins, ' '))
cbb_cataylic_subunit_ko_bins <- cbb_cataylic_subunit_ko_bins %>% rename(., bins = unique.cbb_cataylic_subunit_ko.bins......)
clean_metapath2$Microaerobic[match(cbb_cataylic_subunit_ko_bins$bins, clean_metapath2$bins)]<-1

#Cytochrome.bd.complex: catalytic subunit = cydAB
bd_cataylic_subunit_target <- c("K00425", "K00426") #This is correct
bd_cataylic_subunit_ko <- filter(bac_ko, ko %in% bd_cataylic_subunit_target)
bd_cataylic_subunit_ko.1 <- bd_cataylic_subunit_ko %>% mutate(yesno = 1) %>% 
  distinct %>% spread(ko, yesno, fill = 0) 

bd_cataylic_subunit_ko.2 <- bd_cataylic_subunit_ko.1 %>% subset(., K00425 > 0) %>% subset(., select = c(bins, K00425))
bd_cataylic_subunit_ko.3 <- bd_cataylic_subunit_ko.1 %>% subset(., K00426 > 0) %>% subset(., select = c(bins, K00426))
bd_cataylic_subunit_ko.4 <- bd_cataylic_subunit_ko.3 %>% inner_join(bd_cataylic_subunit_ko.2, ., var = "bins")

bd_cataylic_subunit_ko_bins <- data.frame(unique(bd_cataylic_subunit_ko.4$bins, ' '))
bd_cataylic_subunit_ko_bins <- bd_cataylic_subunit_ko_bins %>% rename(., bins = unique.bd_cataylic_subunit_ko.4.bins......)
clean_metapath2$Microaerobic[match(bd_cataylic_subunit_ko_bins$bins, clean_metapath2$bins)]<-1

#Cytochrome.c.oxidase: catalytic subunit = coxA 
Cytochrome_c_cataylic_subunit_target <- c("K02274") #This is correct
Cytochrome_c_subunit_ko <- filter(bac_ko, ko %in% Cytochrome_c_cataylic_subunit_target)
Cytochrome_c_subunit_ko_bins <- data.frame(unique(Cytochrome_c_subunit_ko$bins, ' '))
Cytochrome_c_subunit_ko_bins <- Cytochrome_c_subunit_ko_bins %>% 
  rename(., bins = unique.Cytochrome_c_subunit_ko.bins......)
clean_metapath2$Aerobic[match(Cytochrome_c_subunit_ko_bins$bins, clean_metapath2$bins)]<-1

#Cytochrome.o.ubiquinol.oxidase: catalytic subunit = cyoB
Cytochrome_o_cataylic_subunit_target <- c("K02298") #This is correct
Cytochrome_o_subunit_ko <- filter(bac_ko, ko %in% Cytochrome_o_cataylic_subunit_target)
Cytochrome_o_subunit_ko_bins <- data.frame(unique(Cytochrome_o_subunit_ko$bins, ' '))
Cytochrome_o_subunit_ko_bins <- Cytochrome_o_subunit_ko_bins %>% 
  rename(., bins = unique.Cytochrome_o_subunit_ko.bins......)
clean_metapath2$Aerobic[match(Cytochrome_o_subunit_ko_bins$bins, clean_metapath2$bins)]<-1

  ##### Identifing dsrAB Reduction (Removing MAGs with dsrA for HS2 oxidation) ####
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_145"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_264"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_110"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_12"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_255"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_265"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_280"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_5"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_136"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_14"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_169"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_38"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_16"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_322"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_143"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_198"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_301"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_87"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_153"] <- 0
clean_metapath2$`Sulfur Reduction`[clean_metapath2$bins == "bin_254"] <- 0




##### Constructing final df for Donors and Acceptors ####
clean_metapath3 <- clean_metapath2 %>%  
  subset(., select = c("bins", "bin.taxa.name", "kingdom",
                       "phylum", "class", "order", "family", 
                       "genus", "Iron Oxidation", "Nitrogen Oxidation", 
                       "Hydrogen Oxidation", "Sulfur Oxidation*", 
                       "Sulfur Reduction",  "Nitrogen Reduction*",
                       "Iron Reduction", "Microaerobic", "Aerobic",
                       "CBB.present",  "CBB_FormI", "CBB_FormII", 
                       "CBB_FormI_and_II", "CBB_Form_Unknown", 
                       "rTCA.present", "Wood.Ljungdahl.present",
                       "Autotrophy", "Heterotrophy"))

########### Cytochrome c oxidase aa3 as possible O2 detox
# Cytochrome c oxidase aa3 requires subunits I,II,and III or I+III(fuse) to function properly. Below I will remove those that do not contain one or more of these
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_144"] <- 0 #missing sub2
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_149"] <- 0 #missing sub2
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_157"] <- 0 #missing sub2
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_178"] <- 0 #missing sub2 and sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_190"] <- 0 #missing sub2 and sub3: WL
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_49"] <- 0  #missing sub2 and sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_75"] <- 0  #missing sub2 and sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_84"] <- 0  #missing sub2 and sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_328"] <- 0 #missing sub2 and sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_365"] <- 0 #missing sub2 and sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_373"] <- 0 #missing sub2
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_113"] <- 0 #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_118"] <- 0 #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_120"] <- 0 #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_124"] <- 0 #missing sub3: WL
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_157"] <- 0 #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_168"] <- 0 #missing sub3: WL
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_175"] <- 0 #missing sub3: WL
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_176"] <- 0 #missing sub3: WL
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_3"] <- 0   #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_40"] <- 0  #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_47"] <- 0  #missing sub3: WL
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_56"] <- 0  #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_224"] <- 0 #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_252"] <- 0 #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_257"] <- 0 #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_369"] <- 0 #missing sub3
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_404"] <- 0 #missing sub3: WL

#all those that contain the fusion of I+III also contain I and III separate
## Now I need to remove the WL guys that have an anerobic organization 
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_106"] <- 0 #anerobic organization of subunits: Lamrabet 2011
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_109"] <- 0 #anerobic organization of subunits: kublanov 2017 ex geneome Caldithrix abyssi
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_364"] <- 0 #anerobic organization of subunits: kublanov 2017 ex geneome Caldithrix abyssi

##Now to remove all Anaerolineae with aa3 cytochrome as these are known stricked anaerobes and likely use this as detoxification
clean_metapath3$Aerobic[clean_metapath3$class == "Anaerolineae"] <- 0 #31 Anaerolineae MAGs
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_240"] <- 0 #Chloroflexota with WL 
clean_metapath3$Aerobic[clean_metapath3$bins == "bin_402"] <- 0 #Chloroflexota with WL

##Now to remove predictions for all Bipolaricaulota with aa3. 
clean_metapath3$Aerobic[clean_metapath3$phylum == "Bipolaricaulota"] <- 0 

#In Bipolaricaulota, this could be used for detoxification. Youssef et al 2019 suggested that this could indicate the ablity to use O2, but this was not tested and is unusual for an organism with WL.
##removing predictions for all FirmicutesH with aa3
clean_metapath3$Aerobic[clean_metapath3$phylum == "FirmicutesH"] <- 0 

#All isolated memebers of the order Thermodesulfovibrionales are obligate anaerobes Umezawa., et al., 2021 and are known to have c-type cytochromes
clean_metapath3$Aerobic[clean_metapath3$order == "Thermodesulfovibrionales"] <- 0 

#Nearly all evaluations of Desulfobacterota consider them to be anaerobic
clean_metapath3$Aerobic[clean_metapath3$phylum == "Desulfobacterota"] <- 0 

# cytochrome bd-type is known to function as a O2 detox
clean_metapath3$Microaerobic[clean_metapath3$class == "Anaerolineae"] <- 0 #this is known to be a o2 detox in Anaerolineae
clean_metapath3$Microaerobic[clean_metapath3$phylum == "Desulfobacterota"] <- 0 
clean_metapath3$Microaerobic[clean_metapath3$phylum == "Bipolaricaulota"] <- 0 

#Set all Microaerobic with WL to 0
WL_micro <- clean_metapath3 %>% subset(., Microaerobic == 1 & Wood.Ljungdahl.present == 1)
clean_metapath3$Microaerobic[clean_metapath3$bins %in% WL_micro$bins] <- 0

clean_metapath4 <- clean_metapath3 %>% 
  subset(., select = c("bins", "bin.taxa.name", "kingdom",
                       "phylum", "class", "order", "family", 
                       "genus", "Iron Oxidation", "Nitrogen Oxidation",
                       "Hydrogen Oxidation", "Sulfur Oxidation*", 
                       "Sulfur Reduction",  "Nitrogen Reduction*",
                       "Iron Reduction", "Microaerobic", "Aerobic",
                       "CBB.present",  "CBB_FormI", "CBB_FormII", 
                       "CBB_FormI_and_II", "CBB_Form_Unknown", 
                       "rTCA.present", "Wood.Ljungdahl.present",
                       "Autotrophy", "Heterotrophy")) %>%
  gather(., "electron.donors", "percent.of.electron.donor", 
         "Iron Oxidation", "Nitrogen Oxidation", 
         "Hydrogen Oxidation", "Sulfur Oxidation*") %>%
  gather(., "electron.acceptor", "percent.of.electron.acceptor", 
         "Sulfur Reduction", "Nitrogen Reduction*", "Iron Reduction", 
         "Microaerobic", "Aerobic")

##################################
# Modifying DF's for Figures 2-3 #
##################################
#### Septerating subsurface from shallow community ####
cutoff_clusters <- data.frame(cut_off_correlation, cluster = cutree(cut_off_hclust, h = 10))# Original set to ten

bin_with_taxa2 <- abun %>% 
  subset(., select = c(bin.taxa.name, bins, kingdom, phylum, 
                       class, order, family, genus, species))
cutoff_clusters2 <- setDT(cutoff_clusters, keep.rownames = TRUE)[]
cutoff_clusters2 <- cutoff_clusters2 %>% rename(bin.taxa.name = rn)
cutoff_clusters2 <- bin_with_taxa2 %>% inner_join(cutoff_clusters2, bin_with_taxa2,  by = "bin.taxa.name")


new_cluster1 <- cutoff_clusters2 %>% subset(cluster == 1) #Arc/Forearc
new_cluster2 <- cutoff_clusters2 %>% subset(cluster == 2) #Prevalent
new_cluster3 <- cutoff_clusters2 %>% subset(cluster == 5) #Outer Forearc
new_cluster4 <- cutoff_clusters2 %>% subset(cluster == 4) #Prevalent
new_cluster5 <- cutoff_clusters2 %>% subset(cluster == 3) #Prevalent

prevalent <- bind_rows(new_cluster2, new_cluster4)
prevalent2 <- bind_rows(prevalent, new_cluster5)

ProvinceSpec <- bind_rows(new_cluster1, new_cluster3)
#write.csv(ProvinceSpec, "Dropbox (Lloyd Lab)/Costa Rica Project/tj_analysis/Costa-Rica-2017/output_data/ProvinceSpec.csv", row.names=FALSE)
#write.csv(prevalent2, "Dropbox (Lloyd Lab)/Costa Rica Project/tj_analysis/Costa-Rica-2017/output_data/prevalent.csv", row.names=FALSE)

#### Outerfore arc analysis (cluster 3) ####
unique(new_cluster3$phylum)
count_phylum <- as.data.frame(table(new_cluster3$phylum)) # Creates a dataframe that tells me the number of MAGs in a phylum

percent_abun_cluster3 <- new_cluster3 %>% subset(., select=c(-1:-3,-5:-9, -37)) %>% 
  group_by(., phylum) %>% summarise_all(sum) 
percent_abun_cluster3$Total <- rowSums(percent_abun_cluster3[-1])
percent_abun_cluster3$Percent <- percent_abun_cluster3$Total/colSums(percent_abun_cluster3['Total'])*100

cluster3Metab <- new_cluster3 %>% inner_join(., clean_metapath3)
cluster3Metab[is.na(cluster3Metab)] <- 0
#write_csv(cluster3Metab, "output_data/cluster3Metab.csv")

metabol_list <- c("Iron Oxidation", 
                     #"Nitrogen Oxidation",
                     "Hydrogen Oxidation", 
                     "Sulfur Oxidation*", 
                     "Sulfur Reduction", 
                     "Nitrogen Reduction*",
                     "Iron Reduction", "Microaerobic", "Aerobic", "CBB.present", 
                     'rTCA.present', 'Wood.Ljungdahl.present')   

#"Ammonium Oxidation", "Denitrification", and "Comammox" are not present
for (i in metabol_list){
  x <- nrow(subset(cluster3Metab, cluster3Metab[,i] >0))#/nrow(cluster3Metab)*100
  print(paste(i, " is in", x, "% of the MAGs"))
}
 #FIX THIS
for (i in metabol_list){
  x <- sum(subset(cluster3Metab, cluster3Metab[,i] >0)[,10:36])/sum(cluster3Metab[,10:36])*100
  print(paste(i, " is in", x, "% of the Reads"))
}



#### Forearc and Arc analysis (cluster 1) ####
ncol(as.data.frame(unique(new_cluster1$phylum)))
cluster1count_phylum <- as.data.frame(table(new_cluster1$phylum)) # Creates a dataframe that tells me the number of MAGs in a phylum

percent_abun_cluster1 <- new_cluster1 %>% subset(., select=c(-1:-3,-5:-9, -37)) %>% 
  group_by(., phylum) %>% summarise_all(sum) 
percent_abun_cluster1$Total <- rowSums(percent_abun_cluster1[-1])
percent_abun_cluster1$Percent <- percent_abun_cluster1$Total/colSums(percent_abun_cluster1['Total'])*100

cluster1Metab <- new_cluster1 %>% inner_join(., clean_metapath3)
cluster1Metab[is.na(cluster1Metab)] <- 0
#write_csv(cluster1Metab, "output_data/cluster1Metab.csv")

metabol_list <- c("Iron Oxidation", "Nitrogen Oxidation", "Hydrogen Oxidation", 
                     "Sulfur Oxidation*", "Sulfur Reduction",  "Nitrogen Reduction*",
                     "Iron Reduction", "Microaerobic", "Aerobic", "CBB.present", 
                     'rTCA.present', 'Wood.Ljungdahl.present')   

for (i in metabol_list){
  x <- nrow(subset(cluster1Metab, cluster1Metab[,i] >0))
  print(paste(i, " is in", x, " MAGs"))
}

for (i in metabol_list){
  x <- sum(subset(cluster1Metab, cluster1Metab[,i] >0)[,10:36])/sum(cluster1Metab[,10:36])*100
  print(paste(i, " is in", x, "% of the Reads"))
}

#### Prevalent cluster 3 ####
nrow(as.data.frame(unique(prevalent2$phylum)))
percent_abun_new_prevalent2 <- prevalent2 %>% subset(., select=c(-1:-3,-5:-9, -37)) %>% 
  group_by(., phylum) %>% summarise_all(sum) 
percent_abun_new_prevalent2$Total <- rowSums(percent_abun_new_prevalent2[-1])
percent_abun_new_prevalent2$Percent <- percent_abun_new_prevalent2$Total/colSums(percent_abun_new_prevalent2['Total'])*100




new_cluster3count_phylum <- as.data.frame(table(prevalent2$phylum)) # Creates a dataframe that tells me the number of MAGs in a phylum

PrevalentMetab <- prevalent2 %>% inner_join(., clean_metapath3)
PrevalentMetab[is.na(PrevalentMetab)] <- 0
#write_csv(PrevalentMetab, "output_data/PrevalentMetab.csv")

metabol_list <- c("Iron Oxidation", "Nitrogen Oxidation", "Hydrogen Oxidation", 
                     "Sulfur Oxidation*", "Sulfur Reduction",  "Nitrogen Reduction*",
                     "Iron Reduction", "Microaerobic", "Aerobic", "CBB.present", 
                     'rTCA.present', 'Wood.Ljungdahl.present') 

for (i in metabol_list){
  x <- nrow(subset(PrevalentMetab, PrevalentMetab[,i] >0))/nrow(PrevalentMetab)*100
  print(paste(i, " is in", x, "% of the MAGs"))
}
nrow(PrevalentMetab %>% subset(CBB_FormI ==1 | CBB_FormI_and_II ==1))
for (i in metabol_list){
  x <- sum(subset(PrevalentMetab, PrevalentMetab[,i] >0)[,10:36])/sum(PrevalentMetab[,10:36])*100
  print(paste(i, " is in", x, "% of the Reads"))
}




#### Dividing up clusters based on prev/provence ####
prevalent <- prevalent2
outerforearc_cutoff_guys <- new_cluster3
forearc_cutoff_guys <- new_cluster1

#Preping above for data output
prev_cutoff_contam_genie <- clean_metapath3[(clean_metapath3$bins %in% prevalent$bins),]

new_cluster3_taxa <- prevalent %>% subset(., select = bins)
cutoff_geniedeep <- clean_metapath3[!(clean_metapath3$bins %in% new_cluster3_taxa$bins),]

ssub <- prev_cutoff_contam_genie
dsub <- cutoff_geniedeep

#Seperate abundance data based on province
outer_sites_abun <- abun_py %>% subset(., select = c('bins', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'QHF2', 'QHS2',
                                                     'EPS','EPF', 'ESF9', 'RSS', 'RSF', 'SIS', 'SIF')) #These site are in the outer forearc

#outer_sites_abun[336,10] = .000000001 #This is only added because the for loop to create the box plot does not work without a value included when making the outer fore arc plot. I need to remove this value from the final data frames that make the plot

fore_sites_abun <- abun_py %>% subset(., select = c('bins','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','CYF', 'CYS',
                                                    'FAS','SLF', 'MTF', 'SLS'))

arc_sites_abun <- abun_py %>% subset(., select = c('bins','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','BQS', 'BQF',
                                                   'BRF1', 'BRS2', 'BRF2', 'BRS1', 'ETS', 'QNF', 'QNS', 'RVF', 'TCF', 'TCS'))#These site are in the arc

shallow_abun <- abun_py %>% subset(., select = c('bins','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','QHF2', 'QHS2',
                                                 'EPS','EPF', 'ESF9', 'RSS', 'RSF', 'SIS', 'SIF', 'CYF', 'CYS', 'FAS', 'SLF', 'SLS',
                                                 'BQS', 'BQF', 'BRF1', 'BRS2', 'BRF2', 'BRS1', 'ETS', 'MTF', 'QNF', 'QNS', 'RVF', 'TCF',
                                                 'TCS'))#These site are in the surface/shallow subsurface
all_site_abun <- abun_py %>% subset(., select = c('bins','kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species','QHF2', 'QHS2',
                                                  'EPS','EPF', 'ESF9', 'RSS', 'RSF', 'SIS', 'SIF', 'CYF', 'CYS', 'FAS', 'SLF', 'SLS',
                                                  'BQS', 'BQF', 'BRF1', 'BRS2', 'BRF2', 'BRS1', 'ETS', 'MTF', 'QNF', 'QNS', 'RVF', 'TCF',
                                                  'TCS'))
# Seperate deep data based on predicted genes/metabolic pathways

# Donors
# Iron Oxidation Donor
deep_Fe2_bins <- dsub %>% subset(.,`Iron Oxidation`> 0 ) %>%
  subset(., select = c(bins)) %>% as.data.frame(.)
# `Hydrogen Oxidation` Donor
deep_H2_bins <- dsub %>% subset(.,`Hydrogen Oxidation`> 0 ) %>%
  subset(., select = c(bins)) %>% as.data.frame(.)
# `Sulfur Oxidation*` Donor
deep_SulfurOxi_bins <- dsub %>% subset(., `Sulfur Oxidation*` > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
#`Nitrogen Oxidation`Donor
deep_NitrogenOxi_bins <- dsub %>% subset(.,`Nitrogen Oxidation`> 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)



#Acceptors
#Sulfur Reduction
deep_SulfurReduc_bins <- dsub %>% subset(., `Sulfur Reduction` > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
# Nitrogen Reduction
deep_NitrogenReduc_bins <- dsub %>% subset(., `Nitrogen Reduction*` > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
#Fe acceptor
deep_IronReduc_bins <- dsub %>% subset(., `Iron Reduction` > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
#O2-cbb3 acceptor
deep_Aerobic_bins <- dsub %>% subset(., Aerobic > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
#O2-bd acceptor
deep_Microaerobic_bins <- dsub %>% subset(., Microaerobic > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)




#Info on metabolics of Carbon fixers only
auto_clean_metapath <- clean_metapath3 %>%
  subset(., CBB.present > 0 | CBB_FormI > 0 | CBB_FormII > 0 | CBB_FormI_and_II > 0 | CBB_Form_Unknown > 0 | rTCA.present > 0 | Wood.Ljungdahl.present > 0) %>%
  subset(., select = c(bins, `Iron Oxidation`:Wood.Ljungdahl.present))

#Carbon fixation
#Calvin Cycle
deep_CBB_bins <- dsub %>% subset(., CBB.present > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
deep_CBB_FI_bins <- dsub %>% subset(., CBB_FormI > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.) #There are 5 bins with form I rubisco cbb in deep
#deep_CBB_FII_bins <- dsub %>% subset(., CBB_FormII > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.) #There is no form II in the deep
#deep_CBB_FI_FII_bins <- dsub %>% subset(., CBB_FormI_and_II > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.) #There is no bin with both form I and II in the deep
deep_CBB_bins_unknown <- dsub %>% subset(., CBB_Form_Unknown > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.) #There are 2 bins with unknown rubisco cbb in deep

#WL Cycle
deep_WL_bins <- dsub %>% subset(., Wood.Ljungdahl.present > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)

#r-TCA 
deep_rtca_bins <- dsub %>% subset(., rTCA.present > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)



# Seperate shallow data based on predicted genes/metabolic pathways
# Donors
# Iron Oxidation Donor
shallow_Fe2_bins <- ssub %>% subset(.,`Iron Oxidation`> 0 ) %>%
  subset(., select = c(bins)) %>% as.data.frame(.)
# `Hydrogen Oxidation` Donor
shallow_H2_bins <- ssub %>% subset(.,`Hydrogen Oxidation`> 0 ) %>%
  subset(., select = c(bins)) %>% as.data.frame(.)
# `Sulfur Oxidation*` Donor
shallow_SulfurOxi_bins <- ssub %>% subset(., `Sulfur Oxidation*` > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
#`Nitrogen Oxidation`Donor
shallow_NitrogenOxi_bins <- ssub %>% subset(.,`Nitrogen Oxidation`> 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)



#Acceptors
#Sulfur Reduction
shallow_SulfurReduc_bins <- ssub %>% subset(., `Sulfur Reduction` > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
# Nitrogen Reduction
shallow_NitrogenReduc_bins <- ssub %>% subset(., `Nitrogen Reduction*` > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
#Fe acceptor
shallow_IronReduc_bins <- ssub %>% subset(., `Iron Reduction` > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
#O2-cbb3 acceptor
shallow_Aerobic_bins <- ssub %>% subset(., Aerobic > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)
#O2-bd acceptor
shallow_Microaerobic_bins <- ssub %>% subset(., Microaerobic > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)



#Carbon fixation
#Calvin Cycle
shallow_CBB_bins <- ssub %>% subset(., CBB.present > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.) # Bin_38's RuBisCo Could not be identified for the gene tree, so we axed it
shallow_CBB_FI_bins <- ssub %>% subset(., CBB_FormI > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.) 
shallow_CBB_FII_bins <- ssub %>% subset(., CBB_FormII > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.) 
shallow_CBB_FI_FII_bins <- ssub %>% subset(., CBB_FormI_and_II > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.) 
shallow_CBB_bins_unknown <- ssub %>% subset(., CBB_Form_Unknown > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)  

#WL Cycle
shallow_WL_bins <- ssub %>% subset(., Wood.Ljungdahl.present > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)

#r-TCA 
shallow_rtca_bins <- ssub %>% subset(., rTCA.present > 0) %>% subset(., select = c(bins)) %>% as.data.frame(.)

#All CBB
AllCBB <- auto_clean_metapath[!(auto_clean_metapath$CBB.present == 0), ]

#Form I
cbb_formI <- auto_clean_metapath[!(auto_clean_metapath$CBB_FormI == 0), ]
cbbfI_meta_and_site_totals <- as.data.frame(colSums(cbb_formI[, -1]))
cbbfI_meta_and_site_totals <- setDT(cbbfI_meta_and_site_totals, keep.rownames = TRUE)[]
cbbfI_meta_and_site_totals <- cbbfI_meta_and_site_totals %>% rename(., "sample_and_metabolism" = rn, "cbbI" = `colSums(cbb_formI[, -1])`)

#Form II
cbb_formII <- auto_clean_metapath[!(auto_clean_metapath$CBB_FormII == 0), ] 
cbbfII_meta_and_site_totals <- as.data.frame(colSums(cbb_formII[, -1]))
cbbfII_meta_and_site_totals <- setDT(cbbfII_meta_and_site_totals, keep.rownames = TRUE)[]
cbbfII_meta_and_site_totals <- cbbfII_meta_and_site_totals %>% rename(., "sample_and_metabolism" = rn, "cbbII" = `colSums(cbb_formII[, -1])`)

#Form I and II
cbb_formI_and_II <- auto_clean_metapath[!(auto_clean_metapath$CBB_FormI_and_II == 0), ] 
cbb_formI_and_II_meta_and_site_totals <- as.data.frame(colSums(cbb_formI_and_II[, -1]))
cbb_formI_and_II_meta_and_site_totals <- setDT(cbb_formI_and_II_meta_and_site_totals, keep.rownames = TRUE)[]
cbb_formI_and_II_meta_and_site_totals <- cbb_formI_and_II_meta_and_site_totals %>% rename(., "sample_and_metabolism" = rn, "cbbI.and.II" = `colSums(cbb_formI_and_II[, -1])`)

#Form Unknown
cbb_form_unknown <- auto_clean_metapath[!(auto_clean_metapath$CBB_Form_Unknown == 0), ]
cbb_form_unknown_meta_and_site_totals <- as.data.frame(colSums(cbb_form_unknown[, -1]))
cbb_form_unknown_meta_and_site_totals <- setDT(cbb_form_unknown_meta_and_site_totals, keep.rownames = TRUE)[]
cbb_form_unknown_meta_and_site_totals <- cbb_form_unknown_meta_and_site_totals %>% rename(., "sample_and_metabolism" = rn, "cbb.unknown" = `colSums(cbb_form_unknown[, -1])`)

#r-TCA
rtca_only <- auto_clean_metapath[!(auto_clean_metapath$rTCA.present == 0), ]
rtca_only_meta_and_site_totals <- as.data.frame(colSums(rtca_only[, -1]))
rtca_only_meta_and_site_totals <- setDT(rtca_only_meta_and_site_totals, keep.rownames = TRUE)[]
rtca_only_meta_and_site_totals <- rtca_only_meta_and_site_totals %>% rename(., "sample_and_metabolism" = rn, "rtca" = `colSums(rtca_only[, -1])`)


wl_only <- auto_clean_metapath[!(auto_clean_metapath$Wood.Ljungdahl.present == 0), ]
wl_only_meta_and_site_totals <- as.data.frame(colSums(wl_only[, -1]))
wl_only_meta_and_site_totals <- setDT(wl_only_meta_and_site_totals, keep.rownames = TRUE)[]
wl_only_meta_and_site_totals <- wl_only_meta_and_site_totals %>% rename(., "sample_and_metabolism" = rn, "wl" = `colSums(wl_only[, -1])`)

combined_cfp_totals <- left_join(cbbfI_meta_and_site_totals, cbbfII_meta_and_site_totals, by = "sample_and_metabolism")
combined_cfp_totals <- left_join(combined_cfp_totals, cbb_formI_and_II_meta_and_site_totals, by = "sample_and_metabolism")
combined_cfp_totals <- left_join(combined_cfp_totals, cbb_form_unknown_meta_and_site_totals, by = "sample_and_metabolism")
combined_cfp_totals <- left_join(combined_cfp_totals, rtca_only_meta_and_site_totals, by = "sample_and_metabolism")
combined_cfp_totals <- left_join(combined_cfp_totals, wl_only_meta_and_site_totals, by = "sample_and_metabolism")

combined_cfp_totals <- combined_cfp_totals %>% column_to_rownames(., var  = "sample_and_metabolism")
combined_cfp_totals2 <-   as.data.frame(t(combined_cfp_totals)) %>% subset(., select = `Hydrogen Oxidation`:Microaerobic) 
combined_cfp_totals2 <- as.data.frame(t(combined_cfp_totals2)) 
combined_cfp_totals2 <-  setDT(combined_cfp_totals2, keep.rownames = "metabolism")[] 
combined_cfp_totals2 <- gather(combined_cfp_totals2, cfp, number, "cbbI":"wl")



# List for abundance function
#Bin Lists


deep_bins_list <- list(deep_H2_bins, deep_SulfurOxi_bins, deep_Fe2_bins, deep_NitrogenOxi_bins,
                       deep_SulfurReduc_bins, deep_NitrogenReduc_bins, 
                       deep_IronReduc_bins, deep_Microaerobic_bins, deep_Aerobic_bins,
                       deep_CBB_bins, deep_rtca_bins, deep_WL_bins)#This is a list with all the deep bin list included. 




shallow_bin_list <- list(shallow_H2_bins, shallow_SulfurOxi_bins, shallow_Fe2_bins, shallow_NitrogenOxi_bins, 
                         shallow_SulfurReduc_bins, shallow_NitrogenReduc_bins, 
                         shallow_IronReduc_bins, shallow_Microaerobic_bins, shallow_Aerobic_bins, 
                         shallow_CBB_bins, shallow_rtca_bins, shallow_WL_bins) #This is a list with all the shallow bin list included.


#These are the names of the metabolic predictions
meta_list <- list("Hydrogen Oxidation",
                  "Sulfur Oxidation*",
                  "Iron Oxidation", 
                  "Nitrogen Oxidation", 
                  "Sulfur Reduction",  
                  "Nitrogen Reduction*",
                  "Iron Reduction", 
                  "Microaerobic", 
                  "Aerobic", "CBB.present", 
                  'rTCA.present', 'Wood.Ljungdahl.present')



my_colors <- c("#7FC97F", "#BEAED4", "#FDC086", 
               "#FFFF99", "#386CB0", "#F0027F", 
               "#BF5B17", "#666666", "#1B9E77", "#D95F02",
               "#7570B3", "#E7298A")
#and are in the same order as the above two lists





outer_meta <- list()
fore_meta <- list()
arc_meta <- list()
shallow_outer_meta <- list()
shallow_fore_meta <- list()
shallow_arc_meta <- list()
shallow_all_meta <- list()

#Outer
for (i in 1:length(deep_bins_list)) {
  outer_meta[[i]] <- meta_abun(deep_bins_list[[i]], outer_sites_abun[,c(1, 9:ncol(outer_sites_abun))], meta_list[[i]])
  
}
outer_meta2 <- bind_rows(outer_meta)

#Insuring that the order of the plot is the correct order
outer_meta2$Metabolism <- as.character(outer_meta2$Metabolism)
outer_meta2$Metabolism <- factor(outer_meta2$Metabolism, levels=unique(outer_meta2$Metabolism))



#Forearc data set for plot
for (i in 1:length(deep_bins_list)) {
  fore_meta[[i]] <- meta_abun(deep_bins_list[[i]], fore_sites_abun[,c(1, 9:ncol(fore_sites_abun))], meta_list[[i]])
  
}
fore_meta2 <- bind_rows(fore_meta)

#Insuring that the order of the plot is the correct order
fore_meta2$Metabolism <- as.character(fore_meta2$Metabolism)
fore_meta2$Metabolism <- factor(fore_meta2$Metabolism, levels=unique(fore_meta2$Metabolism))

#Arc data set for plot
for (i in 1:length(deep_bins_list)) {
  arc_meta[[i]] <- meta_abun(deep_bins_list[[i]], arc_sites_abun[,c(1, 9:ncol(arc_sites_abun))], meta_list[[i]])
  
}
arc_meta2 <- bind_rows(arc_meta)

#Insuring that the order of the plot is the correct order
arc_meta2$Metabolism <- as.character(arc_meta2$Metabolism)
arc_meta2$Metabolism <- factor(arc_meta2$Metabolism, levels=unique(arc_meta2$Metabolism))

############################
#Shallow data set for plot #
############################
# All sites
for (i in 1:length(shallow_bin_list)) {
  shallow_all_meta[[i]] <- meta_abun(shallow_bin_list[[i]], all_site_abun[,c(1, 9:ncol(all_site_abun))], meta_list[[i]])
  
}
shallow_all_meta2 <- bind_rows(shallow_all_meta)

#Insuring that the order of the plot is the correct order
shallow_all_meta2$Metabolism <- as.character(shallow_all_meta2$Metabolism)
shallow_all_meta2$Metabolism <- factor(shallow_all_meta2$Metabolism, levels=unique(shallow_all_meta2$Metabolism))

# Shallow outer
for (i in 1:length(shallow_bin_list)) {
  shallow_outer_meta[[i]] <- meta_abun(shallow_bin_list[[i]], outer_sites_abun[,c(1, 9:ncol(outer_sites_abun))], meta_list[[i]])
  
}
shallow_outer_meta2 <- bind_rows(shallow_outer_meta)

#Insuring that the order of the plot is the correct order
shallow_outer_meta2$Metabolism <- as.character(shallow_outer_meta2$Metabolism)
shallow_outer_meta2$Metabolism <- factor(shallow_outer_meta2$Metabolism, levels=unique(shallow_outer_meta2$Metabolism))

# Shallow Forearc
for (i in 1:length(shallow_bin_list)) {
  shallow_fore_meta[[i]] <- meta_abun(shallow_bin_list[[i]], fore_sites_abun[,c(1, 9:ncol(fore_sites_abun))], meta_list[[i]])
  
}
shallow_fore_meta2 <- bind_rows(shallow_fore_meta)

#Insuring that the order of the plot is the correct order
shallow_fore_meta2$Metabolism <- as.character(shallow_fore_meta2$Metabolism)
shallow_fore_meta2$Metabolism <- factor(shallow_fore_meta2$Metabolism, levels=unique(shallow_fore_meta2$Metabolism))

#Shallow Arc
for (i in 1:length(shallow_bin_list)) {
  shallow_arc_meta[[i]] <- meta_abun(shallow_bin_list[[i]], arc_sites_abun[,c(1, 9:ncol(arc_sites_abun))], meta_list[[i]])
  
}
shallow_arc_meta2 <- bind_rows(shallow_arc_meta)

#Insuring that the order of the plot is the correct order
shallow_arc_meta2$Metabolism <- as.character(shallow_arc_meta2$Metabolism)
shallow_arc_meta2$Metabolism <- factor(shallow_arc_meta2$Metabolism, levels=unique(shallow_arc_meta2$Metabolism))



##########################################
# Looking for genes in mags across sites #
##########################################
## Data needed
  # Metabolic gene dataframe with KO numbers
ko_meta <- read.csv("base_data/ko_for_cfp_mags.csv")
  # bac_ko modified with only MAGs with CFP's
cfp_ko_meta <- bac_ko %>% inner_join(., ko_meta, by = "ko") %>% 
  subset(., select = c(bins, Metabolism, gene_id, gene, ko))
cfp_ko_meta2 <-   cfp_ko_meta[(cfp_ko_meta$bins %in% autotrophs$bins),]
#write_csv(cfp_ko_meta2, "base_data/cfp_ko_meta.csv")

cfp_mag_id <- as.data.frame(autotrophs$bins)
#write_csv(cfp_mag_id, "base_data/cfp_mag_id.csv")
  

###################################################################
# Looking at what cytochromes are found within MAGs containing WL #
###################################################################
WLMAG <- metapath %>% subset(., Wood.Ljungdahl > .6)

######################################################################
# data frames for stacked bargraph of metabolisms at the order level #
######################################################################

outer_meta_order <- list()
fore_meta_order <- list()
arc_meta_order <- list()
shallow_outer_meta_order <- list()
shallow_fore_meta_order <- list()
shallow_arc_meta_order <- list()
shallow_all_meta_order <- list()



#Outer
for (i in 1:length(deep_bins_list)) {
  outer_meta_order[[i]] <- meta_of_orders(deep_bins_list[[i]], outer_sites_abun, meta_list[[i]])
  
}
outer_meta_order2 <- bind_rows(outer_meta_order)

#Insuring that the order of the plot is the correct order
outer_meta_order2$Metabolism <- as.character(outer_meta_order2$Metabolism)
outer_meta_order2$Metabolism <- factor(outer_meta_order2$Metabolism, levels=unique(outer_meta_order2$Metabolism))
outer_meta_order2$class <- reorder(outer_meta_order2$class, outer_meta_order2$Average_Abundance, sum)


#Forearc data set for plot
for (i in 1:length(deep_bins_list)) {
  fore_meta_order[[i]] <- meta_of_orders(deep_bins_list[[i]], fore_sites_abun, meta_list[[i]])
  
}
fore_meta_order2 <- bind_rows(fore_meta_order)

#Insuring that the order of the plot is the correct order
fore_meta_order2$Metabolism <- as.character(fore_meta_order2$Metabolism)
fore_meta_order2$Metabolism <- factor(fore_meta_order2$Metabolism, levels=unique(fore_meta_order2$Metabolism))
fore_meta_order2$class <-  reorder(fore_meta_order2$class, fore_meta_order2$Average_Abundance, sum)


#Arc data set for plot
for (i in 1:length(deep_bins_list)) {
  arc_meta_order[[i]] <- meta_of_orders(deep_bins_list[[i]], arc_sites_abun, meta_list[[i]])
  
}
arc_meta_order2 <- bind_rows(arc_meta_order)

#Insuring that the order of the plot is the correct order
arc_meta_order2$Metabolism <- as.character(arc_meta_order2$Metabolism)
arc_meta_order2$Metabolism <- factor(arc_meta_order2$Metabolism, levels=unique(arc_meta_order2$Metabolism))
arc_meta_order2$class <-  reorder(arc_meta_order2$class, arc_meta_order2$Average_Abundance, sum)

############################
#Shallow data set for plot #
############################
# All

for (i in 1:length(shallow_bin_list)) {
  shallow_all_meta_order[[i]] <- meta_of_orders(shallow_bin_list[[i]], all_site_abun, meta_list[[i]])
  
}
shallow_all_meta_order2 <- bind_rows(shallow_all_meta_order)

#Insuring that the order of the plot is the correct order
shallow_all_meta_order2$Metabolism <- as.character(shallow_all_meta_order2$Metabolism)
shallow_all_meta_order2$Metabolism <- factor(shallow_all_meta_order2$Metabolism, levels=unique(shallow_all_meta_order2$Metabolism))
shallow_all_meta_order2$class <- reorder(shallow_all_meta_order2$class, shallow_all_meta_order2$Average_Abundance, sum)

# Shallow outer
for (i in 1:length(shallow_bin_list)) {
  shallow_outer_meta_order[[i]] <- meta_of_orders(shallow_bin_list[[i]], outer_sites_abun, meta_list[[i]])
  
}
shallow_outer_meta_order2 <- bind_rows(shallow_outer_meta_order)

#Insuring that the order of the plot is the correct order
shallow_outer_meta_order2$Metabolism <- as.character(shallow_outer_meta_order2$Metabolism)
shallow_outer_meta_order2$Metabolism <- factor(shallow_outer_meta_order2$Metabolism, levels=unique(shallow_outer_meta_order2$Metabolism))
shallow_outer_meta_order2$class <- reorder(shallow_outer_meta_order2$class, shallow_outer_meta_order2$Average_Abundance, sum)


for (i in 1:length(shallow_bin_list)) {
  shallow_fore_meta_order[[i]] <- meta_of_orders(shallow_bin_list[[i]], fore_sites_abun, meta_list[[i]])
  
}
shallow_fore_meta_order2 <- bind_rows(shallow_fore_meta_order)

#Insuring that the order of the plot is the correct order
shallow_fore_meta_order2$Metabolism <- as.character(shallow_fore_meta_order2$Metabolism)
shallow_fore_meta_order2$Metabolism <- factor(shallow_fore_meta_order2$Metabolism, levels=unique(shallow_fore_meta_order2$Metabolism))
shallow_fore_meta_order2$class <-  reorder(shallow_fore_meta_order2$class, shallow_fore_meta_order2$Average_Abundance, sum)



for (i in 1:length(shallow_bin_list)) {
  shallow_arc_meta_order[[i]] <- meta_of_orders(shallow_bin_list[[i]], arc_sites_abun, meta_list[[i]])
  
}
shallow_arc_meta_order2 <- bind_rows(shallow_arc_meta_order)

#Insuring that the order of the plot is the correct order
shallow_arc_meta_order2$Metabolism <- as.character(shallow_arc_meta_order2$Metabolism)
shallow_arc_meta_order2$Metabolism <- factor(shallow_arc_meta_order2$Metabolism, levels=unique(shallow_arc_meta_order2$Metabolism))
shallow_arc_meta_order2$class <-  reorder(shallow_arc_meta_order2$class, shallow_arc_meta_order2$Average_Abundance, sum)


########################################
# Assigning colors to specific classes #
########################################
all_orders <- bind_rows(outer_meta_order2, fore_meta_order2, arc_meta_order2, shallow_outer_meta_order2, shallow_fore_meta_order2, shallow_arc_meta_order)

#### This is something to remember!!! 
all_unique_orders <- as.vector(unique(all_orders$class)) #This finds unique orders
## Generate ramdom colors ###
# set.seed(1234) #This sets the seed for the random color generator
# color <- c(distinctColorPalette(length(all_unique_orders))) #This randomly picks colors

## Make those colors permenate
# cat(paste(shQuote(color, type="cmd"), collapse=", "))

#This is a list of mags that were missed due to not having metabolic predictions
added_names <- c("WOR~3, WOR~3", "Micrarchaeota, Micrarchaeia", 
                 "Crenarchaeota, Bathyarchaeia", "Micrarchaeota, Iainarchaeia",    
                 "Nanoarchaeota, Woesearchaeia", "Nanoarchaeota, Aenigmarchaeia", 
                 "unknown phylum, unknown class", "Patescibacteria, Microgenomatia")

#Add MAGs that had`Nitric Oxide Reduction`metabolic predictions
all_unique_orders <- sort(append(all_unique_orders, added_names, after = 78))

colors <- c("#CCF6D4", "#D5AB90", "#8AA2EB", "#58EDAB", "#E4C343",
            "#C4DAC4", "#4BED41", "#C0A9E2", "#EDB972", "#ED4A54",
            "#8A5E55", "#95BEBB", "#C29F55", "#95EAE1", "#AFB280",
            "#EEDFC9", "#E0E67F", "#F0E6E8", "#C038EF", "#F0B7E4",
            "#B461E0", "#E2962F", "#DC9EEB", "#DD9286", "#9DD2EE",
            "#4FE2B7", "#DF7790", "#EBEAA3", "#E3EEBA", "#E37DDD",
            "#D2E55C", "#E44484", "#E5E637", "#5FACBB", "#3E036B",
            "#59ABE1", "#E3CBE5", "#53BF5B", "#9DEE41", "#AE9BA4",
            "#C8D3E6", "#648182", "#8EDE9E", "#BCF397", "#9F7ADC",
            "#E744AE", "#45AB7F", "#4FA190", "#8CC630", "#708EB4", 
            "#C9EFED", "#DA6A3A", "#EBD0A2", "#517DD1", "#EBB2B7", 
            "#BDDE97", "#5968E8", "#EA8CC3", "#7040CB", "#ABE66A", 
            "#92E9C0", "#663B88", "#6F2537", "#BE4850", "#AB7341", 
            "#52CFC4", "#6339EF", "#95BA99", "#84F193", "#56D1ED", 
            "#55EDF1", "#93B04D", "#4AEE7D", "#EF986C", "#9D6BA0", 
            "#B1BFF0", "#A73CA0", "#5DF2DC", "#A73CA0", "#046063", 
            "#809b9c", "#08c986", "#4e8773", "#6a7052", "#4e025c", 
            "#592f16")


names(colors) <- all_unique_orders # This assigns randomly picked colors to unique orders
colorsdf1 <- as.data.frame(stack(colors))

names(colorsdf1)[1] <- 'ClassColors'
names(colorsdf1)[2] <- 'class'

redox_list <- c("Iron Oxidation", "Nitrogen Oxidation", "Hydrogen Oxidation", 
                  "Sulfur Oxidation*", "Sulfur Reduction",  "Nitrogen Reduction*",
                  "Iron Reduction", "Microaerobic", "Aerobic")   



###############
# MAG Quality #
###############
MAGQuality <- MAGQuality %>% 
  inner_join(., bin_with_taxa2, by = "bins")

MAGQuality[is.na(MAGQuality)] <- 0 #Convert NA's to 0's

MAGQuality.1 <- MAGQuality %>% 
  mutate_if(is.numeric, ~1 * (. > 0)) #Convert all numbers >0 to 1

MAGQuality.1$trnaTotal <- rowSums(MAGQuality.1[, c(5:24)]) #Find total number of individual tRNA's

MAGQuality.1TRNA <- MAGQuality.1 %>% 
  subset(., trnaTotal >= 18) #Find those MAGs with the required >= 18 tRNA's

HiqhQualMAGs <- MAGQuality.1TRNA %>% 
  subset(., `16s` == 1 & `23s` == 1 & `5s` == 1) #Of those with >= 18 tRNA's, find those with all 3 required rRNA

# Making MAG quality table for suplemental document
MAGQuality.1 <- MAGQuality.1 %>% 
  inner_join(., bin_quality[,c("bins", "size", "completeness", "contamination")])
MAGQualTable <- MAGQuality.1 %>% 
  subset(., select = c(bins, bin.taxa.name, size, completeness,
                       contamination, `23s`, `16s`, `5s`, trnaTotal))
MAGQualTable$bin.taxa.name <- trimws(MAGQualTable$bin.taxa.name, whitespace = "^[^,]+,\\s*")
MAGQualTable <- MAGQualTable %>% 
  rename("Bin #" = bins, "Taxa ID" = bin.taxa.name,
         "Completeness (%)" = completeness, "Contamination (%)" = contamination,
         "23S rRNA" = `23s`, "16S rRNA" = `16s`, "5S rRNA" = `5s`, "# of individual tRNA" = trnaTotal)
#write_csv(MAGQualTable, "ManuscriptFiles/Supplementals/MAGQualTable.csv")







