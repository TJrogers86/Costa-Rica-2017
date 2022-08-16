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
  meta_abund = two %>% mutate_at(vars(one_of("SO4", "Fe3", "N2O", "NO", "NO3",
                                              "O2_aerobic", "O2_micro",
                                              "H2", "H2S", "S2O3", "Fe2", "NH3")),
                                  funs(case_when(. > 0 ~ rowMeans(two[,21:ncol(two)]), TRUE ~ 0))) %>% # This looks for the prediction of a don/acc and replaces it with the total abun of bin in province
    mutate(Standard_error = select(., DCO_LLO_Bv4v5..MTF_MT170219:DCO_LLO_Bv4v5..ESF9_ES_F9) %>% apply(.,1,standard_error)) %>%
    subset(., select = c(bins:NH3, Standard_error))# We dont need the data frame with site abun values, so we select only predection abun columns
  indx <- sapply(meta_abund, is.factor)
  meta_abund[indx] <- lapply(meta_abund[indx], function(x) as.numeric(as.character(x)))
  test <- gather(meta_abund, Metabolism, meta.ave.abundance, SO4:NH3)
  test <- test %>% mutate(., Occurance = case_when(meta.ave.abundance > 0 ~ 1 , TRUE ~ 0))
 return(test)
}

# This function constructs the data frame for the deep community barplots for carbon fixers
shallow <- function(cfp_accep_don_df, abun_df){
  one = cfp_accep_don_df %>% inner_join(., abun_df)#This joins the don/acc df with site abun info of the bins # might delete
  two = one[one$bins %in% new_cluster3_taxa$bins,]
  meta_abund = two %>% mutate_at(vars(one_of("SO4", "Fe3", "N2O", "NO", "NO3",
                                             "O2_aerobic", "O2_micro",
                                             "H2", "H2S", "S2O3", "Fe2", "NH3")),
                                 funs(case_when(. > 0 ~ rowMeans(two[,21:ncol(two)]), TRUE ~ 0))) %>% # This looks for the prediction of a don/acc and replaces it with the total abun of bin in province
    mutate(Standard_error = select(., DCO_LLO_Bv4v5..MTF_MT170219:DCO_LLO_Bv4v5..ESF9_ES_F9) %>% apply(.,1,standard_error)) %>%
    subset(., select = c(bins:NH3, Standard_error))# We dont need the data frame with site abun values, so we select only predection abun columns
  indx <- sapply(meta_abund, is.factor)
  meta_abund[indx] <- lapply(meta_abund[indx], function(x) as.numeric(as.character(x)))
  test <- gather(meta_abund, Metabolism, meta.ave.abundance, SO4:NH3)
  test <- test %>% mutate(., Occurance = case_when(meta.ave.abundance > 0 ~ 1 , TRUE ~ 0))
  return(test)
}

# This function constructs the data frame for the province/prevalent community barplots for carbon fixers
CFP.Redox.provence <- function(sites, CFP, group){
  cfp <- sites %>% inner_join(CFP, ., by = "bins") %>% semi_join(., group, by= "bins")
  setDT(cfp)[class=="", class:= "unknown class"]
  cfp2 <- unite(cfp, class, c(phylum, class), sep=", ", remove=FALSE)
  cfp3 = cfp2 %>% mutate_at(vars(one_of('Fe2', 'NH3', 'H2', 'H2S', 'S', 'S2O3',
                                        'SO3', 'SO4', 'N2O', 'NO',  'NO2_reduction', 'NO3',
                                        'Fe3', 'O2_micro', 'O2_aerobic')),
                            funs(case_when(. > 0 ~ rowMeans(cfp2[,31:ncol(cfp2)]), TRUE ~ 0)))
  cfp4 <- gather(cfp3, Metabolism, meta.ave.abundance, SO4:NH3) 
  cfp5 <- cfp4 %>% mutate(., Occurance = case_when(meta.ave.abundance > 0 ~ 1 , TRUE ~ 0))
  cfp5$Metabolism <- factor(cfp5$Metabolism, levels= c('Fe2', 'NH3', 'H2', 'H2S', 'S', 'S2O3',
                                                       'SO3', 'SO4', 'N2O', 'NO',  'NO2_reduction', 'NO3',
                                                       'Fe3', 'O2_micro', 'O2_aerobic'))
  return(cfp5)
}

# Order for the above function
CFPStackOrder <- function(sites, CFP, group){
  cfp <- sites %>% inner_join(CFP, ., by = "bins") %>% semi_join(., group, by= "bins")
  setDT(cfp)[class=="", class:= "unknown class"]
  cfp2 <- unite(cfp, class, c(phylum, class), sep=", ", remove=FALSE)
  orderdf <- cfp2[,-c(1:24,26:30)]
  orderdf2 <- orderdf %>%
    mutate(Total = rowSums(.[,2:ncol(.)]), .after = 1) %>%
    mutate(Ave = rowMeans(.[,3:ncol(.)]), .after = 2) %>%
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



##################################
# Functions for the CCA analysis #
##################################

# I want to create specific colors for each phylum
Unique_phylum <- as.vector(unique(abun_py$phylum)) #This finds unique orders

## Generate unique colors for each phylum ###
Phycols <- c("#BC76E5", "#84B969", "#D75F48", "#CEBA7C", "#E9BFE8", "#E6D0B5", "#6D4AA5", "#68DC68", "#62ACCD", "#7BEF42", "#B49EDE", "#59F2AD", "#6B86EB", "#ACB6AD", "#98585E", "#A6F1D0", "#8A5AE4",
             "#E6E547", "#9FC4A1", "#BDED8E", "#A0D5EC", "#CFEBE6", "#DF5DAE", "#52ECDB", "#54936C", "#556E72", "#59D3B0", "#6A73A9", "#EAA4B1", "#8AE8E4", "#68A9EB", "#DE9A7B", "#D839E9", "#A8E04D",
             "#55D7E9", "#B2999D", "#7527E0", "#6EB8B4", "#9EE9AB", "#E692E9", "#DB89B5", "#E9D8E1", "#E4E284", "#E0EFBD", "#D89F3D", "#ACBAE0", "#E84170", "#DE51D3")

names(Phycols) <- Unique_phylum
Phycolsdf1 <- as.data.frame(stack(Phycols))
names(Phycolsdf1)[1] <- 'PhyColors'
names(Phycolsdf1)[2] <- 'phylum'


#Function to assign species correlations and color to the species names, and create class item
#Assign Phycolors
ModifySpecDF <- function(moddf,
                         colordf,#colors for unique classes
                         phycolordf, #colors for unique phylum
                         abuns, #Full abundf for the group
                         geo
){
  #Assigning Class color Scaling with species
  classes.ScaleSpecies <- as.data.frame(scores(moddf, scaling = "species")$species) %>% 
    rownames_to_column(., var = "bins") %>%
    inner_join(., abuns[,c(2:5)],
               by = 'bins')
  classes.ScaleSpecies$class <- classes.ScaleSpecies$class %>% replace_na("unknown class")
  classes.ScaleSpecies <- classes.ScaleSpecies %>% inner_join(., phycolordf, by = "phylum")
  classes.ScaleSpecies <- unite(classes.ScaleSpecies, class, c(phylum, class), sep=", ", remove=FALSE)
  classes.ScaleSpecies <- classes.ScaleSpecies %>% inner_join(., colordf, by = "class")
  
  
  #MAG Size as a sum of abundance across volcanic arc
  abuns$BinSum <- rowMeans(abuns[,c(10:36)])
  abuns$ScaledAver <- scales::rescale(abuns$BinSum, to = c(1, 20))
  BinSums <- abuns[,c(2,39)]
  classes.ScaleSpecies <- classes.ScaleSpecies %>% inner_join(., BinSums, by = "bins")
  
  #Assigning Class color Scaling with species
  classes.ScaleSites <- as.data.frame(scores(moddf, scaling = "sites")$species) %>% 
    rownames_to_column(., var = "bins") %>%
    inner_join(., abuns[,c(2:5)],
               by = 'bins')
  classes.ScaleSites$class <- classes.ScaleSites$class %>% replace_na("unknown class")  
  classes.ScaleSites <- classes.ScaleSites %>% inner_join(., phycolordf, by = "phylum")
  classes.ScaleSites <- unite(classes.ScaleSites, class, c(phylum, class), sep=", ", remove=FALSE)
  classes.ScaleSites <- classes.ScaleSites %>% inner_join(., colordf, by = "class")
  
  #MAG Size as a sum of abundance across volcanic arc
  abuns$BinSum <- rowMeans(abuns[,c(10:36)])
  abuns$ScaledAver <- scales::rescale(abuns$BinSum, to = c(1, 20))
  BinSums <- abuns[,c(2,39)]
  classes.ScaleSites <- classes.ScaleSites %>% inner_join(., BinSums, by = "bins")
  
  #Assigning Province Colors Species scaling
  Sites.ScalingSpecies <- as.data.frame(scores(moddf, scaling = "species")$sites) %>% rownames_to_column(., var = "sites")
  Location <- geo[,c(1,2,5)] %>% rownames_to_column(., var = "sites")
  Sites.ScalingSpecies2 <- Sites.ScalingSpecies %>% inner_join(., Location, by = "sites") %>% column_to_rownames(., var = "sites")
  Sites.ScalingSpecies3 <- Sites.ScalingSpecies2 %>% unite(., abbrev, c(code, type), sep = "", remove=TRUE)
  Sites.ScalingSpeciescols <- with(Sites.ScalingSpecies3,
                                   data.frame(province = levels(t(province)),
                                              color = c(arc = "orange", forearc = "darkblue", outer_forearc = "green")
                                   ))
  Sites.ScalingSpecies3$Cols <- Sites.ScalingSpecies3 %>% inner_join(., Sites.ScalingSpeciescols, by = "province")
  
  #Assigning Province Colors Species scaling
  Sites.ScalingSites <- as.data.frame(scores(moddf, scaling = "sites")$sites) %>% rownames_to_column(., var = "sites")
  Location <- geo[,c(1,2,5)] %>% rownames_to_column(., var = "sites")
  Sites.ScalingSites2 <- Sites.ScalingSites %>% inner_join(., Location, by = "sites") %>% column_to_rownames(., var = "sites")
  Sites.ScalingSites3 <- Sites.ScalingSites2 %>% unite(., abbrev, c(code, type), sep = "", remove=TRUE)
  Sites.ScalingSitescols <- with(Sites.ScalingSites3,
                                 data.frame(province = levels(t(province)),
                                            color = c(arc = "orange", forearc = "darkblue", outer_forearc = "green")
                                 ))
  
  Sites.ScalingSites3$Cols <- Sites.ScalingSites3 %>% inner_join(., Sites.ScalingSitescols, by = "province")
  return(list("classes.ScaleSpecies" = classes.ScaleSpecies, 
              "classes.ScaleSites" = classes.ScaleSites, 
              "Sites.ScalingSites3" = Sites.ScalingSites3,
              "Sites.ScalingSpecies3" = Sites.ScalingSpecies3,
              "modd" = moddf))
}




















