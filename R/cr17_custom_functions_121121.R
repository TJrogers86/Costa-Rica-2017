all_site <- c("DCO_LLO_Bv4v5..MTF_MT170219", "DCO_LLO_Bv4v5..BRF1_BR170218_1",
              "DCO_LLO_Bv4v5..RSF_RS170216", "DCO_LLO_Bv4v5..SIS_SI170217",
              "DCO_LLO_Bv4v5..QHF2_QH170212_2", "DCO_LLO_Bv4v5..CYF_CY170214",
              "DCO_LLO_Bv4v5..SLF_SL170214", "DCO_LLO_Bv4v5..BRS2_BR170218_2",
              "DCO_LLO_Bv4v5..CYS_CY170214", "DCO_LLO_Bv4v5..RVF_RV170221",
              "DCO_LLO_Bv4v5..RSS_RS170216", "DCO_LLO_Bv4v5..TCF_TC170221",
              "DCO_LLO_Bv4v5..QNF_QN170220", "DCO_LLO_Bv4v5..EPS_EP170215",
              "DCO_LLO_Bv4v5..EPF_EP170215", "DCO_LLO_Bv4v5..SLS_SL170214_S",
              "DCO_LLO_Bv4v5..QNS_QN170220", "DCO_LLO_Bv4v5..BRF2_BR170218_2",
              "DCO_LLO_Bv4v5..TCS_TC170221_2", "DCO_LLO_Bv4v5..FAS_FA170219",
              "DCO_LLO_Bv4v5..BQS_BQ170220", "DCO_LLO_Bv4v5..ETS_ET170220",
              "DCO_LLO_Bv4v5..QHS2_QH170213_2", "DCO_LLO_Bv4v5..SIF_SI170217",
              "DCO_LLO_Bv4v5..BQF_BQ170220", "DCO_LLO_Bv4v5..BRS1_BR170218_1",
              "DCO_LLO_Bv4v5..ESF9_ES_F9")





#############################
# Function to make total KO #
#############################
TotalKO <- function(file1 = "", file2 = ""){
  ko1 <- read_csv(file1)
  ko2 <- read_csv(file2)
  bac_ko <- bind_rows(ko1, ko2)
  return(bac_ko)
}

############################
# Calculate Standard Error #
############################
#Function to calculate standard error
standard_error <- function(x) sd(x) / sqrt(length(x))

#################################################################
# Constructed DF of Order level metabolic prediction abundances #
#################################################################
# This function is used to construct the data frame for the metabolisms within provinces of the order level
meta_of_orders <- function(bin_df, abun_df, name_of_matabolism, custom_colors){
  one = abun_df %>% subset(., abun_df$bins %in% bin_df$bins) %>% 
    remove_rownames %>% column_to_rownames(var="bins") #This finds the abundance for bins specific to that province
  two = one[as.logical(rowSums(one[,8:ncol(one)] != 0)), ] #This fliters out any row that sums to zero
  three = setDT(two)[class=="", class:= "unknown class"]
  four = unite(three, class, c(phylum, class), sep=", ", remove=FALSE)
  five = setDT(four)[order=="", order:= "unknown order"]
  six = unite(five, order, c(class, order), sep=", ", remove=FALSE)
  six$Number = 1
  seven = six %>% subset(., select = c(class, 8:ncol(six))) %>% group_by(class) %>%
    summarise_each(funs(sum))
  seven$Average_Abundance = rowMeans(seven[2:(length(seven)-1)])
  eight = seven %>% subset(., select = c(class,Average_Abundance,Number)) %>%
    add_column(Metabolism = name_of_matabolism, .before = 1) 
  return(eight)
}

######################################################################
# Constructed DF of the abundances of metabolic prediction (Boxplot) #
######################################################################
# This function is used to construct the data frame for the Boxplots
meta_abun <- function(bin_df, abun_df, name_of_matabolism, custom_colors){
  one = abun_df %>% subset(., abun_df$bins %in% bin_df$bins) %>% 
    remove_rownames %>% column_to_rownames(var="bins")
  two = one[as.logical(rowSums(one != 0)), ]
  site_abun <- as.data.frame(t(colSums(two))) %>%
    add_column(Metabolism = name_of_matabolism, .before = 1) 
  long_abun = site_abun %>% gather(., Site, Abundance, -Metabolism)
  site_num = as.data.frame(t(colSums(two != 0))) %>%
    add_column(Metabolism = name_of_matabolism, .before = 1)
  long_num = site_num %>% gather(., Site, Number, -Metabolism)
  final_df = inner_join(long_abun, long_num, var = Metabolism)#working here
  return(final_df)
}

##################################
# Functions for CFP REDOX graphs #
##################################
# This function constructs the data frame for the deep community barplots for carbon fixers
deep <- function(cfp_accep_don_df, abun_df){
  one = cfp_accep_don_df %>% inner_join(., abun_df)#This joins the don/acc df with site abun info of the bins # might delete
  two = one[!(one$bins %in% new_cluster3_taxa$bins),]
  meta_abund = two %>% mutate_at(vars(one_of(redox_list)),
                                  funs(case_when(. > 0 ~ rowMeans(two[,all_site]), TRUE ~ 0))) %>% # This looks for the prediction of a don/acc and replaces it with the total abun of bin in province
    mutate(Standard_error = select(., all_site) %>% apply(.,1,standard_error)) %>%
    subset(., select = c(bins:Aerobic, Standard_error))# We dont need the data frame with site abun values, so we select only predection abun columns
  indx <- sapply(meta_abund, is.factor)
  meta_abund[indx] <- lapply(meta_abund[indx], function(x) as.numeric(as.character(x)))
  test <- gather(meta_abund, Metabolism, meta.ave.abundance, `Iron Oxidation`:Aerobic)
  test <- test %>% mutate(., Occurance = case_when(meta.ave.abundance > 0 ~ 1 , TRUE ~ 0))
 return(test)
}



# This function constructs the data frame for the deep community barplots for carbon fixers
shallow <- function(cfp_accep_don_df, abun_df){
  one = cfp_accep_don_df %>% inner_join(., abun_df)#This joins the don/acc df with site abun info of the bins # might delete
  two = one[one$bins %in% new_cluster3_taxa$bins,]
  meta_abund = two %>% mutate_at(vars(one_of(redox_list)),
                                 funs(case_when(. > 0 ~ rowMeans(two[,all_site]), TRUE ~ 0))) %>% # This looks for the prediction of a don/acc and replaces it with the total abun of bin in province
    mutate(Standard_error = select(., all_site) %>% apply(.,1,standard_error)) %>%
    subset(., select = c(bins:Aerobic, Standard_error))# We dont need the data frame with site abun values, so we select only predection abun columns
  indx <- sapply(meta_abund, is.factor)
  meta_abund[indx] <- lapply(meta_abund[indx], function(x) as.numeric(as.character(x)))
  test <- gather(meta_abund, Metabolism, meta.ave.abundance, `Iron Oxidation`:Aerobic)
  test <- test %>% mutate(., Occurance = case_when(meta.ave.abundance > 0 ~ 1 , TRUE ~ 0))
  return(test)
}

# This function constructs the data frame for the province/prevalent community barplots for carbon fixers
CFP.Redox.provence <- function(sites, CFP, group){
  cfp <- sites %>% inner_join(CFP, ., by = "bins") %>% semi_join(., group, by= "bins")
  setDT(cfp)[class=="", class:= "unknown class"]
  cfp2 <- unite(cfp, class, c(phylum, class), sep=", ", remove=FALSE)
  cfp3 = cfp2 %>% mutate_at(vars(one_of("Iron Oxidation",         "Nitrogen Oxidation",     "Hydrogen Oxidation",     "Sulfur Oxidation*",
                                        "Sulfur Reduction",       "Nitrogen Reduction*",   "Iron Reduction",         "Microaerobic",           "Aerobic",
                                        "CBB.present",            "CBB_FormI",             "CBB_FormII",             "CBB_FormI_and_II",       "CBB_Form_Unknown",
                                        "rTCA.present",           "Wood.Ljungdahl.present")),funs(case_when(. > 0 ~ rowMeans(cfp2[,25:ncol(cfp2)]), TRUE ~ 0)))
  cfp4 <- gather(cfp3, Metabolism, meta.ave.abundance,"Iron Oxidation":"Aerobic")
  cfp5 <- cfp4 %>% mutate(., Occurance = case_when(meta.ave.abundance > 0 ~ 1 , TRUE ~ 0))
  cfp5$Metabolism <- factor(cfp5$Metabolism, levels= redox_list)
  return(cfp5)
}







# Order for the above function
CFPStackOrder <- function(sites, CFP, group){
  cfp <- sites %>% inner_join(CFP, ., by = "bins") %>% semi_join(., group, by= "bins")
  setDT(cfp)[class=="", class:= "unknown class"]
  cfp2 <- unite(cfp, class, c(phylum, class), sep=", ", remove=FALSE)
  orderdf <- cfp2[,-c(1:18,20:24)]
  orderdf2 <- orderdf %>%
    mutate(Total = rowSums(.[,2:ncol(.)]), .after = 1) %>%
    mutate(Ave = rowMeans(.[,2:ncol(.)]), .after = 2) %>%
    subset(., select = c(1:3)) %>% group_by(class) %>%
    summarize_all(sum) %>% arrange(desc(Ave)) %>%
    pull(class)
  return(orderdf2)
}





# Function needed for abundance and number of metabolisms associated with cfps
MAG_bar_plots <- function(df, # Data Frame to be ploted
                          name = "", # Name of over all plot for x-axis label 
                          y_height, # This is the wanted height value of the plot
                          type = "Abundance of MAGs",
                          legpos = "none"
                          ){
  plt = if(type =="Abundance of MAGs" ){
    ggplot(df, aes(x=Metabolism, y=meta.ave.abundance, fill = class)) +
      geom_bar(stat='identity')+ scale_fill_manual(values = colors) +
      theme(axis.text.x=element_text(angle=45,hjust=1))+
      xlab(name)+
      ylab(type) + theme(text = element_text(size=30),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"))+
      theme(legend.position = legpos)+ 
      theme(axis.title.y = element_blank(),
            axis.title.x.top = element_text(size = 30),
            axis.ticks.y=element_blank()
      ) + coord_cartesian(ylim = c(0, y_height))
    
  }else{ggplot(df, aes(x=Metabolism, y = Occurance, fill = class)) +
      geom_bar(stat = 'identity') + scale_fill_manual(values = colors) +
      theme(axis.text.x=element_text(angle=45,hjust=1))+
      xlab(name)+
      ylab(type) + theme(text = element_text(size=30),
                         panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank(),
                         panel.background = element_blank(), 
                         axis.line = element_line(colour = "black"))+
      theme(legend.position = legpos)+ 
      theme(axis.title.y = element_blank(),
                                             axis.title.x.top = element_text(size = 12),
                                             axis.ticks.y=element_blank()
      ) + coord_cartesian(ylim = c(0, y_height))
  }
  return(plt)
}

#########################################
# Assigning colors to top10 class level #
#########################################
TopOther <- function(cluster, coldf, filter){
  Top10 <-  filter %>%
    semi_join(.,cluster) %>%
    mutate(., "Total Abundance" = rowSums((.[,9:ncol(.)]), na.rm = TRUE)) #Finds MAG totals in the province area
  
  Top10two <-  setDT(Top10)[class=="", class:= "unknown class"]
  Top10three <- unite(Top10two, class, c(phylum, class), sep=", ", remove=FALSE) %>% subset(., select = c("class", "Total Abundance")) %>%
    group_by(class) %>% summarise_all(sum) #This summerizes all class level totals
  
  Top10four <- Top10three %>% column_to_rownames(., "class") %>% arrange(desc(`Total Abundance`)) 
  names <- Top10four %>% rownames()
  Top10four2 <- Top10four %>% slice(1:10) %>% rownames() #This finds the top10
  top10cols <- colors %>% subset(., names(.) %in% Top10four2) #This adds color to the top10
  
  #Other colors
  other <- Top10three %>% column_to_rownames(., "class") %>% arrange(desc(`Total Abundance`)) %>% slice(11:n()) %>% rownames()
  othercols <- replicate(length(other), "#dbdbdb")
  names(othercols) <- other
  allcols <- c(top10cols, othercols)
  return(list("allcols" = allcols, "names" = names))
}

###########################################
# Calling in and modifying labels of tree #
###########################################

GetSingleGenusSpecies <- function(x) {
  return(paste(strsplit(x, " |_")[[1]][1:2], collapse=" "))
}

GetAllGenusSpecies <- function(x) {
  y <- sapply(x, GetSingleGenusSpecies)
  gsub("\\..*","",y)
}

GetTreeWithNameProcessing <- function(treefile) {
  phy <- read.tree(treefile)
  phy$tip.label <- unname(GetAllGenusSpecies(phy$tip.label))
  phy$tip.label <- make.unique(phy$tip.label)
  phy <- multi2di(phy)
  phy$edge.length[phy$edge.length==0]<-1e-16
  return(phy)
}




















