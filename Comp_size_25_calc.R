# Pre-processing of data for paper "Microbial community analyses of composts are influenced by
# particle size fraction and DNA extraction method"


# Author: Anja Logo
# Started: 21.01.25
# Last changes: 02.12.25

# README

# qPCR can be independently run from the other parts

# Packages and R env specifications
source(file = "Comp_size_25_setup.R")

# Functions----
# Function to summarize ASV tables
smry_asv_seq <-function(data){
  print(paste("Number of sequences", sum(data)))
  print(paste(" Number of ASVs", nrow(data)))
  print(paste("Number of sequences per sample"))
  print(data %>% colSums() %>% sort() %>% as.data.frame() %>% summary()) # sequences per samples
  print(paste("Number of ASVs per sample"))
  print(ifelse(data >0, 1,0) %>% colSums() %>% sort() %>% as.data.frame() %>% summary())
  
}

# Remove percentages
remove_percentage <- function(column) {
  sub("\\(.*", "", column)
}

# Average beta diversity
avg.vegdist = function(x, method ="bray") {
  data = merge(design, x, by= 0) # merge with design file
  rownames(data) = data[,1] # change rownames
  data[,1] = NULL
  data2 = as.data.frame(aggregate(data[,(ncol(design)+1):ncol(data)], list(data$treatment), mean))
  rownames(data2) <-data2[,1]
  data2[,1] <- NULL
  vegdist(data2, method = "bray")
}

# ANOVA summary
ANOVA_summary <- function(data = design.alpha, x1 = "compost" , x2= "size", x3= "kit", y = "sobs") {
  
  # Fit ANOVA model
  ANOVA <- aov(as.formula(paste(y, "~", x1,"*", x2,"*", x3)), data = data)
  
  # Summarize ANOVA
  ANOVA_summary <- summary(ANOVA)
  
  # Tukey HSD test
  TUKEY <- TukeyHSD(ANOVA, conf.level = 0.95)
  
  # Compact letters for Tukey HSD
  compact_letters <- multcompView::multcompLetters4(ANOVA, TUKEY)
  
  # Convert compact letters to data frame
  compact_letters_df <- as.data.frame.list(compact_letters[[paste0(x1,":",x2, ":", x3)]])
  
  # Summarize data
  data_sum <- data %>%  
    group_by(.data[[x1]], .data[[x2]], .data[[x3]]) %>%
    summarise(
      mean = mean(.data[[y]], na.rm = TRUE),
      max = max(.data[[y]], na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    arrange(desc(mean)) %>%    
    add_column(letters = compact_letters_df$Letters)
  
  # Add position for text labels
  data_sum$pos.letter <- data_sum$max+ 0.03*data_sum$max  # Adjust as needed
  
  # Return results
  list(
    ANOVA_summary = ANOVA_summary,
    TUKEY = TUKEY,
    compact_letters = compact_letters,
    data_sum = data_sum
  )
}

# ASV table-----
  ## FUNGI
  asv.F.all <- read.table(file = "sequencing_data/9.all.ASV_map_F.txt", header =TRUE)
  # Replace "." by "-" (I do not know why this is changed when reading the data)
  colnames(asv.F.all) <- gsub(x = colnames(asv.F.all), pattern = "\\.", replacement = "-")
  # Get rid of the appendix in the name
  colnames(asv.F.all) <- gsub(x = colnames(asv.F.all), pattern = "\\_S.*", replacement = "")
  # Get rid of the F
  colnames(asv.F.all) <- gsub(x = colnames(asv.F.all), pattern = "F_", replacement = "")
  
  # Adapt to names for bacteria
  colnames(asv.F.all) <- gsub(x = colnames(asv.F.all), pattern = "_", replacement = "")
  colnames(asv.F.all) <- gsub(x = colnames(asv.F.all), pattern = "-", replacement = "_")
  rownames(asv.F.all) <- asv.F.all$ASV; asv.F.all$ASV <-NULL
  
  ## BACTERIA
  # ASV table fungi
  asv.B.all <-read.table(file="sequencing_data/9.all.ASV_map_B.txt", header=TRUE)
  rownames(asv.B.all) <- asv.B.all$ASV; asv.B.all$ASV <-NULL
  # Get rid of the appendix in the name
  colnames(asv.B.all) <- gsub(x = colnames(asv.B.all), pattern = "\\_S.*", replacement = "")
  
  # colnames(asv.B.all) ==  colnames(asv.F.all) # Check if the colnames are the same
  
# tax table-----

  ## FUNGI
  column_tax <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
  tax.F.all <- read.table(file = "sequencing_data/8.all.ASV_ITSx.ITS2.UNITEv10_sh_dynamic_all.wang.taxonomy", header=F, sep='\t',row.names = 1)
  rownames <-rownames(tax.F.all)
  tax.F.all<-as.data.frame(str_split_fixed(tax.F.all[,1],';',8))
  tax.F.all$V8<-NULL
  colnames(tax.F.all)<- column_tax
  rownames(tax.F.all) <- rownames; rm(rownames)
  
  # unclassified to unclassifed
  
  expected_starts <- c(Kingdom = "k_", Phyla = "p_", Class = "c_", Order = "o_", Family = "f_", Genus = "g_", Species = "s_")
  
  for (col_name in names(expected_starts)) {
    # Check if the column exists in your dataframe to avoid errors
    if (col_name %in% colnames(tax.F.all)) {
      # Get the expected start letter for the current column
      expected_letter <- expected_starts[col_name]
      
      # Loop through each row in the current column
      for (row_index in 1:nrow(tax.F.all)) {
        # Extract the first letter of the current cell
        current_letter <- substr(tax.F.all[[col_name]][row_index], 1, 2)
        
        # Check if the first letter matches the expected letter
        if (current_letter != expected_letter) {
          # If not, replace the cell value with "unclassified"
          tax.F.all[[col_name]][row_index] <- "unclassified"
        }}}}
  
  # Lowest classfication
  tax.F.all <-tax.F.all %>%
    dplyr::mutate(lowest = case_when(
      Species != "unclassified" ~ Species,
      Genus != "unclassified" ~ Genus,
      Family != "unclassified" ~ Family,
      Order != "unclassified" ~ Order,
      Class != "unclassified" ~ Class,
      Phyla != "unclassified" ~ Phyla,
      Kingdom != "unclassfied" ~ Kingdom,
      TRUE ~ "Unclassified"
    ))

  # Remove percentages
  tax.F.all[c(column_tax, "lowest")] <- lapply(tax.F.all[c(column_tax, "lowest")], remove_percentage)

  # Remove "x__" prior to the name for all columns except for lowest

  tax.F.all[column_tax] <-lapply(tax.F.all[column_tax], function(column) {sub("*.__", "", column)})
  
  # Remove non-fungal sequences
  tax.F<-tax.F.all[grep('Fungi',tax.F.all$Kingdom),]
  asv.F <-asv.F.all[rownames(tax.F),]
  
  ## BACTERIA
  tax.B.all <- read.table(file = "sequencing_data/8.all.ASV_metaxa.silva_v138_2_var34.wang.taxonomy", header=F, sep='\t',row.names = 1)
  rownames <-rownames(tax.B.all)
  tax.B.all<-as.data.frame(str_split_fixed(tax.B.all[,1],';',8))
  tax.B.all$V8<-NULL; tax.B.all$V7<-NULL
  colnames(tax.B.all)<- column_tax[1:6]
  rownames(tax.B.all) <- rownames; rm(rownames)
  
  
  # Add prefixes conditionally for each column
  for (col in column_tax[1:6]) {
    tax.B.all[[col]] <- ifelse(
      !grepl("unclassified|unknown", tax.B.all[[col]], ignore.case = TRUE),
      paste0(tolower(substr(col, 1, 1)), "__", tax.B.all[[col]]), # Lowercase prefix
      tax.B.all[[col]]
    )
  }

  # unclassified to unclassifed
  
  expected_starts <- c(Kingdom = "k_", Phyla = "p_", Class = "c_", Order = "o_", Family = "f_", Genus = "g_")
  
  for (col_name in names(expected_starts)) {
    # Check if the column exists in your dataframe to avoid errors
    if (col_name %in% colnames(tax.B.all)) {
      # Get the expected start letter for the current column
      expected_letter <- expected_starts[col_name]
      
      # Loop through each row in the current column
      for (row_index in 1:nrow(tax.B.all)) {
        # Extract the first letter of the current cell
        current_letter <- substr(tax.B.all[[col_name]][row_index], 1, 2)
        
        # Check if the first letter matches the expected letter
        if (current_letter != expected_letter) {
          # If not, replace the cell value with "unclassified"
          tax.B.all[[col_name]][row_index] <- "unclassified"
        }}}}
  
  # Lowest classfication
  tax.B.all <-tax.B.all %>%
    dplyr::mutate(lowest = case_when(
      Genus != "unclassified" ~ Genus,
      Family != "unclassified" ~ Family,
      Order != "unclassified" ~ Order,
      Class != "unclassified" ~ Class,
      Phyla != "unclassified" ~ Phyla,
      Kingdom != "unclassfied" ~ Kingdom,
      TRUE ~ "Unclassified"
    ))
  

  # Remove percentages
  tax.B.all[c(column_tax[1:6], "lowest")] <- lapply( tax.B.all[c(column_tax[1:6], "lowest")], remove_percentage)
  
  # Remove "x__" prior to the name for all columns except for lowest
  
  tax.B.all[column_tax[1:6]] <-lapply(tax.B.all[column_tax[1:6]], function(column) {sub("*.__", "", column)})
  
  # Remove non-bacterial sequences
  tax.B <- tax.B.all[tax.B.all$Kingdom =="Bacteria",]
  tax.n.B <- tax.B.all[tax.B.all$Kingdom !="Bacteria",]
  tax.n.C <- tax.B[tax.B$Order == "Chloroplast",]
  tax.B <- tax.B[tax.B$Order != "Chloroplast",]
  tax.n.M = tax.B[tax.B$Family == "Mitochondria",] 
  tax.B <- tax.B[tax.B$Family != "Mitochondria",] 
  
  #sink(file="output/Number_of_seq_non_bacteria.txt")
  paste("Number of non-bacterial sequences",asv.B.all[rownames(tax.n.B),] %>% sum())
  paste("Number of non-bacterial ASVs", tax.n.B %>% nrow())
  paste("Number of Chloroplast sequences",asv.B.all[rownames(tax.n.C),] %>% sum())
  paste("Number of Chloroplast ASVs", tax.n.C %>% nrow())
  paste("Number of Mitrochondria sequences",asv.B.all[rownames(tax.n.M),] %>% sum())
  paste("Number of Mitrochondria ASVs", tax.n.M %>% nrow())
  #sink()
  rm(tax.n.B, tax.n.C, tax.n.M)
  
  # Only bacterial ASVs
  asv.B <-asv.B.all[rownames(tax.B),]
  
  rm(row_index, expected_letter, expected_starts, column_tax, col, col_name,
     current_letter)
  
# Normalization----
  design <- read.csv(file ="data/Particle_size_design.csv", sep =";")
  design <- design[order(design$ID, asv.F %>% colnames()),]
  design$ID == (asv.F %>% colnames())
  rownames(design) <- design$ID
  
  ## FUNGI-----
  
  # Rarefy abundance table 100 times
  ISS.iters.F <- mclapply(as.list(1:100), function(x) rrarefy(t(asv.F), min(colSums(asv.F))),
                        mc.cores = 1)
  
  # For alpha diversity calculations & beta diversity based on Ackermann methods
  ISS.array <- laply(ISS.iters.F, as.matrix)
  ISS <- apply(ISS.array, 2:3, median) # median of all iterations
  rm(ISS.array)
  
  # Rarefied and robustly detected ASVs only composts
  ISS <- ISS %>% t() %>% as.data.frame()
  ISS01 <- ifelse(ISS > 0, 1, 0) %>% t() # presence/absence
  ISS01.meta <-merge(design["treatment"], ISS01, by = 0)  # merge with design file
  row.names(ISS01.meta)<-ISS01.meta$Row.names;ISS01.meta$Row.names<-NULL
  ISS_01_rob <-as.data.frame(aggregate(ISS01.meta[,2:ncol(ISS01.meta)], list(ISS01.meta$treatment), median)) # calculate median

  rownames(ISS_01_rob) =ISS_01_rob$Group.1 ; ISS_01_rob$Group.1  <- NULL
  select <-colSums(ISS_01_rob) >0 # Only select ASVs that are at least robustly detected in one compost
  ISS_01_rob_red <-ISS_01_rob[,select] # Select only robustly detected ASVs (unrarefied data)
  sum(colSums(ISS_01_rob_red)== 0) # check if really all the 0 are gone
  ISS.rob <- ISS[colnames(ISS_01_rob_red),] %>% t() %>% as.data.frame()
  rownames(ISS.rob) == rownames(design) # Check order
  ISS.rob.avg <- ISS.rob %>% aggregate(list(design$treatment), mean) # Calculate the average of the four samples
  rownames(ISS.rob.avg) <-ISS.rob.avg[,1]; ISS.rob.avg[,1] <- NULL
  ISS.F <- ISS[rowSums(ISS) != 0,]
  ISS.rob.F <- ISS.rob[rowSums(ISS.rob) != 0,]
  ISS.rob.F <- ISS.rob.F %>% t() %>% as.data.frame()
  ISS.rob.avg.F <- ISS.rob.avg[rowSums(ISS.rob.avg) != 0,]
  ISS.rob.avg.F <- ISS.rob.avg.F %>% t() %>% as.data.frame()
  
  rm(ISS, ISS01, ISS01.meta, ISS_01_rob, select, ISS_01_rob_red, ISS.rob, ISS.rob.avg)
  
  # All sequences
  #write.table(asv.F.all, file ="data/asv.all.F.txt")
  # Only fungi
  #write.table(asv.F, file ="data/asv.F.txt")
  # Fungi after rarefying
  #write.table(ISS.F, file ="data/ISS.F.txt")
  # Only robustly detected
  #write.table(ISS.rob.F, file = "data/ISS.rob.F.txt")
  # Average of replicates
  #write.table(ISS.rob.avg.F, file = "data/ISS.rob.avg.F.txt") 
  # tax file only with fungi
  #write.table(tax.F, file = "data/tax.f.txt")
  
  
  #sink("output/Number_of_seq.F.txt")
  paste("Fungi")
  paste("All")
  smry_asv_seq(asv.F.all)
  paste("Only fungi")
  smry_asv_seq(asv.F)
  paste("after rarefying")
  smry_asv_seq(ISS.F)
  paste("only robustly detected 2/3 replicates")
  smry_asv_seq(ISS.rob.F)
  #sink()
  
  ## BACTERIA-----
    # Rarefy abundance table 100 times
  ISS.iters.B <- mclapply(as.list(1:100), function(x) rrarefy(t(asv.B), min(colSums(asv.B))),
                        mc.cores = 1)
  
  # For alpha diversity calculations & beta diversity based on Ackermann methods
  ISS.array <- laply(ISS.iters.B, as.matrix)
  ISS <- apply(ISS.array, 2:3, median) # median of all iterations
  rm(ISS.array)
  
  # Rarefied and robustly detected ASVs only composts
  ISS <- ISS %>% t() %>% as.data.frame()
  ISS01 <- ifelse(ISS > 0, 1, 0) %>% t() # presence/absence
  ISS01.meta <-merge(design["treatment"], ISS01, by = 0)  # merge with design file
  row.names(ISS01.meta)<-ISS01.meta$Row.names;ISS01.meta$Row.names<-NULL
  ISS_01_rob <-as.data.frame(aggregate(ISS01.meta[,2:ncol(ISS01.meta)], list(ISS01.meta$treatment), median)) # calculate median

  rownames(ISS_01_rob) <-ISS_01_rob$Group.1 ; ISS_01_rob$Group.1  <- NULL
  select <-colSums(ISS_01_rob) >0 # Only select ASVs that are at least robustly detected in one compost
  ISS_01_rob_red <-ISS_01_rob[,select] # Select only robustly detected ASVs (unrarefied data)
  sum(colSums(ISS_01_rob_red)== 0) # check if really all the 0 are gone
  ISS.rob <- ISS[colnames(ISS_01_rob_red),] %>% t() %>% as.data.frame()
  rownames(ISS.rob) == rownames(design) # Check order
  ISS.rob.avg <- ISS.rob %>% aggregate(list(design$treatment), mean) # Calculate the average of the four samples
  rownames(ISS.rob.avg) <-ISS.rob.avg[,1]; ISS.rob.avg[,1] <- NULL
  ISS.B <- ISS[rowSums(ISS) != 0,]
  ISS.rob.B <- ISS.rob[rowSums(ISS.rob) != 0,]
  ISS.rob.B <- ISS.rob.B %>% t() %>% as.data.frame()
  ISS.rob.avg.B <- ISS.rob.avg[rowSums(ISS.rob.avg) != 0,]
  ISS.rob.avg.B <- ISS.rob.avg.B %>% t() %>% as.data.frame()
  
  rm(ISS, ISS01, ISS01.meta, ISS_01_rob, select, ISS_01_rob_red, ISS.rob, ISS.rob.avg)
  
  # All sequences
  #write.table(asv.B.all, file ="data/asv.all.B.txt")
  # Only fungi
  #write.table(asv.B, file ="data/asv.B.txt")
  # Bacteria after rarefying
  #write.table(ISS.B, file ="data/ISS.B.txt")
  # Only robustly detected
  #write.table(ISS.rob.B, file = "data/ISS.rob.B.txt")
  # Average of replicates
  #write.table(ISS.rob.avg.B, file = "data/ISS.rob.avg.B.txt") 
  # tax file only with fungi
  #write.table(tax.B, file = "data/tax.B.txt")
  
  
  #sink("output/Number_of_seq.B.txt")
  paste("Bacteria")
  paste("All")
  smry_asv_seq(asv.B.all)
  paste("Only bacteria")
  smry_asv_seq(asv.B)
  paste("after rarefying")
  smry_asv_seq(ISS.B)
  paste("only robustly detected 2/3 replicates")
  smry_asv_seq(ISS.rob.B)
  #sink()
  
  
# Rarefaction curves----
  
  ## FUNGI
  #pdf(file="figures/F_rarefaction_fract.pdf", width = 6, height = 4)
  rarecurve(t(asv.F), step = 100, col = "blue", cex = 0.6) 
  #dev.off()
  
  ## BACTERIA
  #pdf(file="figures/B_rarefaction_fract.pdf", width = 6, height = 4)
  rarecurve(t(asv.B), step = 100, col = "blue", cex = 0.6) 
  #dev.off()
  
# Alpha diversity----
  
  ## FUNGI
  # Observed richness
  sobs <- mclapply(ISS.iters.F, function(x) specnumber(x), mc.cores = 1)
  sobs <- apply(laply(sobs, as.matrix), 2, mean)
  
  # Shannon diversity
  shannon <- mclapply(ISS.iters.F, function(x) diversity(x, index = "shannon"))
  shannon <- apply(laply(shannon, as.matrix), 2, mean)
  
  # Inverse Simpson diversity
  invsimpson <- mclapply(ISS.iters.F, function(x) diversity(x, index = "invsimpson"))
  invsimpson <- apply(laply(invsimpson, as.matrix), 2, mean)
  
  # Pielou's evenness
  evenness <- shannon/log(sobs)
  
  # Bind to data.frame
  ISS.alpha.F <- data.frame(cbind(sobs, shannon, invsimpson, evenness))
  rm(sobs, shannon, invsimpson, evenness)
  colnames(ISS.alpha.F) <-  c("sobs.F", "shannon.F", "invsimpson.F", "evenness.F")
  
  ## BACTERIA
  # Observed richness
  sobs <- mclapply(ISS.iters.B, function(x) specnumber(x), mc.cores = 1)
  sobs <- apply(laply(sobs, as.matrix), 2, mean)
  
  # Shannon diversity
  shannon <- mclapply(ISS.iters.B, function(x) diversity(x, index = "shannon"))
  shannon <- apply(laply(shannon, as.matrix), 2, mean)
  
  # Inverse Simpson diversity
  invsimpson <- mclapply(ISS.iters.B, function(x) diversity(x, index = "invsimpson"))
  invsimpson <- apply(laply(invsimpson, as.matrix), 2, mean)
  
  # Pielou's evenness
  evenness <- shannon/log(sobs)
  
  # Bind to data.frame
  ISS.alpha.B <- data.frame(cbind(sobs, shannon, invsimpson, evenness))
  rm(sobs, shannon, invsimpson, evenness)
  colnames(ISS.alpha.B) <-  c("sobs.B", "shannon.B", "invsimpson.B", "evenness.B")
  
  ISS.alpha <-merge(ISS.alpha.F, ISS.alpha.B, by=0)
  rownames(ISS.alpha) <- rownames(ISS.alpha.B); ISS.alpha$Row.names <-NULL
  write.csv(ISS.alpha, file ="data/ISS.alpha.F.B.csv")
  
# Beta diversity-----
  #data  <- read.table("ISS.bray.txt") %>% as.matrix() %>% as.dist() # How to read in the data again
  
  ## FUNGI
  rownames(ISS.iters.F[[1]]) == rownames(design) # Check if the rownames are the same

  # Bray-curtis no transformation
  ISS.iters.bray <- mclapply(ISS.iters.F, function(x) vegdist(x, method = "bray"), mc.cores = 1)
  ISS.array.bray <- laply(ISS.iters.bray, as.matrix)
  ISS.bray.F <- as.dist(apply(ISS.array.bray, 2:3, mean))
  rm(ISS.iters.bray, ISS.array.bray)
  #write.table(ISS.bray.F %>% as.matrix() %>% as.data.frame(), file ="data/ISS.bray.F.txt")
  
  # Bray-curtis with sqrt
  ISS.iters.bray <- mclapply(ISS.iters.F, function(x) vegdist(sqrt(x), method = "bray"), mc.cores = 1)
  ISS.array.bray <- laply(ISS.iters.bray, as.matrix)
  ISS.bray.sqrt.F <- as.dist(apply(ISS.array.bray, 2:3, mean))
  rm(ISS.iters.bray, ISS.array.bray)
  #write.table(ISS.bray.sqrt.F %>% as.matrix() %>% as.data.frame(), file ="data/ISS.bray.sqrt.F.txt")
  
  # Jaccard
  ISS.iters.jac <- mclapply(ISS.iters.F, function(x) vegdist((x > 0) * 1, method = "jaccard"),
                            mc.cores = 1)
  ISS.array.jac <- laply(ISS.iters.jac, as.matrix)
  ISS.jac.F <- as.dist(apply(ISS.array.jac, 2:3, mean))
  rm(ISS.iters.jac, ISS.array.jac)
  #write.table(ISS.jac.F %>% as.matrix() %>% as.data.frame(), file ="data/ISS.jac.F.txt")
  
  # Euclidean distance without transformation
  ISS.iters.euc <- mclapply(ISS.iters.F, function(x) vegdist(x, method = "euclidean"), mc.cores = 1)
  ISS.array.euc <- laply(ISS.iters.euc, as.matrix)
  ISS.euc.F <- as.dist(apply(ISS.array.euc, 2:3, mean))
  rm(ISS.iters.euc, ISS.array.euc)
  #write.table(ISS.euc.F %>% as.matrix() %>% as.data.frame(), file ="data/ISS.euc.F.txt")
  
  # Euclidean distance with sqrt transformation
  ISS.iters.euc <- mclapply(ISS.iters.F, function(x) vegdist(sqrt(x), method = "euclidean"), mc.cores = 1)
  ISS.array.euc <- laply(ISS.iters.euc, as.matrix)
  ISS.euc.sqrt.F <- as.dist(apply(ISS.array.euc, 2:3, mean))
  rm(ISS.iters.euc, ISS.array.euc)
  #write.table(ISS.euc.F %>% as.matrix() %>% as.data.frame(), file ="data/ISS.euc.sqrt.F.txt")
  

  # Mantel test among different community metrics
  
  bray.F <-vegdist(asv.F %>% t(), method="bray")
  jac.F <-vegdist(asv.F %>% t(), method="jaccard")
  euc.F <-vegdist(asv.F %>% t(), method="euclidean")

  #sink(paste0("output/Mantel_normalization_F.txt"))
  paste("Comparison without and with data normalization")
  mantel(bray.F, ISS.bray.F)
  mantel(jac.F, ISS.jac.F)
  mantel(euc.F, ISS.euc.F)
  
  paste("comparison without and with sqrt transformation")
  mantel(ISS.bray.F, ISS.bray.sqrt.F)
  mantel(ISS.euc.F, ISS.euc.sqrt.F)
  
  paste("Comparison between different methods")
  mantel(ISS.bray.F, ISS.jac.F)
  mantel(ISS.bray.F, ISS.euc.F)
  mantel(ISS.jac.F, ISS.euc.F)
  #sink()

  ## BACTERIA
  rownames(ISS.iters.B[[1]]) == rownames(design) # Check if the rownames are the same
  
  # Bray-curtis no transformation
  ISS.iters.bray <- mclapply(ISS.iters.B, function(x) vegdist(x, method = "bray"), mc.cores = 1)
  ISS.array.bray <- laply(ISS.iters.bray, as.matrix)
  ISS.bray.B <- as.dist(apply(ISS.array.bray, 2:3, mean))
  rm(ISS.iters.bray, ISS.array.bray)
  #write.table(ISS.bray.B %>% as.matrix() %>% as.data.frame(), file ="data/ISS.bray.B.txt")
  
  # Bray-curtis with sqrt
  ISS.iters.bray <- mclapply(ISS.iters.B, function(x) vegdist(sqrt(x), method = "bray"), mc.cores = 1)
  ISS.array.bray <- laply(ISS.iters.bray, as.matrix)
  ISS.bray.sqrt.B <- as.dist(apply(ISS.array.bray, 2:3, mean))
  rm(ISS.iters.bray, ISS.array.bray)
  write.table(ISS.bray.sqrt.B %>% as.matrix() %>% as.data.frame(), file ="data/ISS.bray.sqrt.B.txt")
  
  # Jaccard
  ISS.iters.jac <- mclapply(ISS.iters.B, function(x) vegdist((x > 0) * 1, method = "jaccard"),
                            mc.cores = 1)
  ISS.array.jac <- laply(ISS.iters.jac, as.matrix)
  ISS.jac.B <- as.dist(apply(ISS.array.jac, 2:3, mean))
  rm(ISS.iters.jac, ISS.array.jac)
  #write.table(ISS.jac.B %>% as.matrix() %>% as.data.frame(), file ="data/ISS.jac.B.txt")
  
  # Euclidean distance without transformation
  ISS.iters.euc <- mclapply(ISS.iters.B, function(x) vegdist(x, method = "euclidean"), mc.cores = 1)
  ISS.array.euc <- laply(ISS.iters.euc, as.matrix)
  ISS.euc.B <- as.dist(apply(ISS.array.euc, 2:3, mean))
  rm(ISS.iters.euc, ISS.array.euc)
  write.table(ISS.euc.B %>% as.matrix() %>% as.data.frame(), file ="data/ISS.euc.B.txt")
  
  # Euclidean distance with sqrt transformation
  ISS.iters.euc <- mclapply(ISS.iters.B, function(x) vegdist(sqrt(x), method = "euclidean"), mc.cores = 1)
  ISS.array.euc <- laply(ISS.iters.euc, as.matrix)
  ISS.euc.sqrt.B <- as.dist(apply(ISS.array.euc, 2:3, mean))
  rm(ISS.iters.euc, ISS.array.euc)
  #write.table(ISS.euc.B %>% as.matrix() %>% as.data.frame(), file ="data/ISS.euc.sqrt.B.txt")
  
  # Mantel test among different community metrcis
  
  bray.B <-vegdist(asv.B %>% t(), method="bray")
  jac.B <-vegdist(asv.B %>% t(), method="jaccard")
  euc.B <-vegdist(asv.B %>% t(), method="euclidean")
  
  #sink(paste0("output/Mantel_normalization_B.txt"))
  paste("Comparison without and with data normalization")
  mantel(bray.B, ISS.bray.B)
  mantel(jac.B, ISS.jac.B)
  mantel(euc.B, ISS.euc.B) # Jere has a big effect why?
  
  paste("comparison without and with sqrt transformation")
  mantel(ISS.bray.B, ISS.bray.sqrt.B)
  mantel(ISS.euc.B, ISS.euc.sqrt.B)
  
  paste("Comparison between different methods")
  mantel(ISS.bray.B, ISS.jac.B)
  mantel(ISS.bray.B, ISS.euc.B)
  mantel(ISS.jac.B, ISS.euc.B)
  #sink()
  
# qPCR normalization-----
  
# Load data
  data_0525 <- read_excel("data/250521_qPCR_particle_fraction_16S_ITS.xlsx", sheet ="data_qPCR")
  data_0306 <- read_excel("data/250521_qPCR_particle_fraction_16S_ITS.xlsx", sheet ="data_qPCR_rep16_0306")
  design_05_25 <- read_excel("data/250521_qPCR_particle_fraction_16S_ITS.xlsx", sheet ="qPCR_design")
  
  data_0525_meta <- merge(data_0525, design_05_25, by ="Well")
  data_0525_meta_ITS <- data_0525_meta %>% filter(gene =="ITS")
  
  data_0306_meta <- merge(data_0306, design_05_25, by ="Well")
  data_0306_meta_16S <- data_0306_meta %>% filter(gene =="16S")
  
  # Tecan 1:20 diluted raw DNA, 16S also 1:20, ITS2 1:10 diluted
  
  # Design file for analysis 
  design <- read.csv(file="data/Particle_size_design.csv", sep=";")  
  rownames(design) <- design$ID
  design <-dplyr::mutate(design, compost = case_when(
    compost == "K20" ~ "A",
    compost == "K21" ~ "B",
    compost == "K22" ~ "C"
  ))
  design$compost <- as.factor(design$compost)
  design$kit <- as.factor(design$kit)
  design$size <- as.factor(design$size)
  design$dna23 <- design$DNA_conc; design$DNA_conc <-NULL
  
  # Sequencing data
  asv.F <-read.table(file="data/asv.F.txt") # unrarefied only fungi
  asv.B <-read.table(file="data/asv.B.txt") # unrarefied only bacteria
  
  # Tecan 1:20 diluted raw DNA
  data_DNA <- read_excel("data/250521_qPCR_particle_fraction_16S_ITS.xlsx", sheet ="data_tecan")
  
  ## DNA concentration Tecan reader 21.05.25-----
  # Calculate std curve
  std_curve <-data_DNA %>% filter(name_F =="Std")
  data_DNA <- data_DNA %>% filter(name_F !="Std")

  model <-list(lm(measure ~conc, data = std_curve))

  text <- paste("R-square:", round(summary(model[[1]])$r.squared,3))
  ggplot(data= std_curve, aes(x = conc, y = measure))+ geom_smooth(method ="lm")+geom_point() +
    xlab("Concentrations [ng DNA/mL]")+ ylab("Relative fluorescence units")+ theme_classic()+
    annotate("text", x= 0.5, y=0, hjust=0, label=text)
  
  #ggsave(filename ="figures/qPCR/250521_DNA_measurments_std_curve.pdf", height =5, width = 5)
  
  model <- lm(measure ~conc, data = std_curve)
  data_DNA$dna_conc <- ((1/ model$coefficients[2])*data_DNA$measure - ( model$coefficients[1]/ model$coefficients[2]))*20
  
  # Measure the mean of the two replicates
  
  data_DNA_mean <- data_DNA %>%group_by(name_F) %>% dplyr::summarize(dna_25 = mean(dna_conc))
  
  design2 <- left_join(design, data_DNA_mean, by =c("F.number"="name_F"))
  
  ggplot(data = design2, aes(x = dna23, y = dna_25)) + geom_point() +geom_smooth(method ="lm")+
    stat_cor(method = "pearson", label.x = 20, label.y = 50)
 
  #ggsave(filename = "figures/qPCR/DNA_measumrents_23_25.pdf", height =5, width = 5)
  
  # Compare DNA concentration across treatments
  
  
  results <-ANOVA_summary(data = design2, x1 = "compost", x2 = "size", x3 = "kit", y= "dna_25")
  
  ggplot(data = design2, aes(x = size, y = dna_25, colour = compost, fill=compost)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = interaction(compost, size)), position = position_dodge(width = 0.5)) +
    facet_wrap(~kit, labeller = labeller(kit = c("M" = "NucleoMag", "S" = "NucleoSpin"))) +
    geom_text(
      data = results$data_sum, 
      aes(x = size, y = pos.letter, label = letters, fill=compost),  # Use `colour = compost` to match group
      position = position_dodge(width = 0.5),  # Apply consistent dodging
      vjust = -0.5,  # Place text slightly above points,
      color="darkred",
      size = 5, 
      inherit.aes = FALSE
    ) +
    ylab("DNA concentration [ng/qL]") +
    scale_color_manual(values = color.compost) +
    scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
    bg_theme +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())
  
  #ggsave(filename = "figures/DNA_conc/Comp_particle_DNA_conc25_metabarcoding.pdf", height = 5, width = 10)
  
  ANOVA <-aov(dna_25  ~ compost*size*kit, data = design2)
  #pdf(file="figures/DNA_conc/ANOVA_DNA_conc_residual_plots_25.pdf", height =6, width = 6)
  par(mfrow=c(2,2))
  plot(ANOVA)
  #dev.off()
  
  #sink(file="output/ANOVA_DNA_conc_25.txt")
  summary(ANOVA)
  #sink()
  
  ## ITS 21.05.25 -----
  data_F_r <-data_0525_meta_ITS %>% filter(grepl("B", sample))
  data_F_r$Cq <- as.numeric(data_F_r$Cq)
  data_F_r_std <- data_0525_meta_ITS[!data_0525_meta_ITS$conc =="NA",]
  data_F_r_std <- data_F_r_std %>% filter(conc != "NTC" & conc !="1" & conc !="2") # We had some contamination in the standard curve made with apotheken water
  data_F_r_std$Cq <- data_F_r_std$Cq %>% as.numeric()
  data_F_r_std$conc <- data_F_r_std$conc %>% as.numeric()
  
  efficiency_data <- data_F_r_std %>%
    summarise(
      model = list(lm(Cq ~ conc)),  # Store model
      slope = coef(model[[1]])[2],  # Extract slope
      intercept = coef(model[[1]])[1],
      efficiency = (10^(-1/slope) - 1) * 100,  # Calculate efficiency
      r_squared = summary(model[[1]])$r.squared,  # Extract R²
      .groups = "drop"
    ) 
  
  results_lm1 <- paste0("Slope: ", round(efficiency_data$slope,2),", Intercept: ", round(efficiency_data$intercept,1))
  results_lm2 <- paste0("Efficiency: ", round(efficiency_data$efficiency,1),"%, R-square: ", round(efficiency_data$r_squared, 3))
  
  
  g1 <-ggplot(data= data_F_r_std, aes(x = conc, y = Cq))+  geom_smooth(method ="lm") +geom_point()+
    ylab("Ct")+ xlab("Log10(gBlock ITS2 copies/μL)")+ theme_classic()+
    annotate("text", x = 3, y = 11, label =results_lm1, hjust =0)+
    annotate("text", x = 3, y = 10, label =results_lm2, hjust =0)
  
  #ggsave(g1, filename = "figures/qPCR/210525_ITS_r_gradient_linear.pdf", height =5, width = 5)
  # Efficiency of 95% and R2 of 0.997 is good
  
  # Formula 10^((Cq-b)/m), still DNA copies still log10 transformed
  
  data_F_r$logcopy <- ((data_F_r$Cq-efficiency_data$intercept)/efficiency_data$slope)
  data_F_r$F.number <-gsub("B", "F", data_F_r$sample)
  
  
  data_F_r_meta <-merge(data_F_r %>% group_by(F.number) %>% dplyr::summarize(logcopy = mean(logcopy)),
                        design2 %>% select(ID, F.number, compost, rep, kit, size, dna_25), by = "F.number")
  
  # Normalize to 1 ng/mL
  data_F_r_meta$copies_per_ngDNA <-((10^data_F_r_meta$logcopy)*10)/data_F_r_meta$dna_25
 
  # ANOVA to test differences
  
  # Without adapting DNA concentration
  ggplot(data = data_F_r_meta, aes(x = size, y = copies_per_ngDNA, colour = compost, fill=compost)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    bg_theme+
    geom_line(aes(group = interaction(compost, size)), position = position_dodge(width = 0.5)) +
    facet_wrap(~kit, labeller = labeller(kit = c("M" = "NucleoMag", "S" = "NucleoSpin"))) +
    scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())+
    scale_color_manual(values = color.compost)+
    ylab("Copies ITS2/ ng DNA")
  

  results <-ANOVA_summary(data = data_F_r_meta, x1 = "compost", x2 = "size", x3 = "kit", y= "copies_per_ngDNA")
  
  ggplot(data = data_F_r_meta, aes(x = size, y = copies_per_ngDNA, colour = compost, fill=compost)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = interaction(compost, size)), position = position_dodge(width = 0.5)) +
    facet_wrap(~kit, labeller = labeller(kit = c("M" = "NucleoMag", "S" = "NucleoSpin"))) +
    geom_text(
      data = results$data_sum, 
      aes(x = size, y = pos.letter, label = letters, fill=compost),  # Use `colour = compost` to match group
      position = position_dodge(width = 0.5),  # Apply consistent dodging
      vjust = -0.5,  # Place text slightly above points,
      color="darkred",
      size = 5, 
      inherit.aes = FALSE
    ) +
    ylab("Copies ITS2/ ng DNA") +
    scale_color_manual(values = color.compost) +
    scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
    bg_theme +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())
  
  #ggsave(filename = "figures/qPCR/210525_ITS_CQ_norm_1ng_qL.pdf", height = 5, width = 10)
  
  ANOVA <-aov(copies_per_ngDNA  ~ compost*size*kit, data = data_F_r_meta)
  #pdf(file="figures/qPCR/ANOVA_copies_per_ng_DNA_plots_F.pdf", height =6, width = 6)
  par(mfrow=c(2,2))
  plot(ANOVA)
  #dev.off()
  
  #sink(file="output/qPCR/ANOVA_copy_per_ng_DNA_F.txt")
  summary(ANOVA)
  #sink()

  ## 16S 03.06.25------
  
  data_B_r <-data_0306_meta_16S %>% filter(grepl("F", sample))
  data_B_r$Cq <- as.numeric(data_B_r$Cq)
  data_B_r_std <- data_0306_meta_16S[!data_0306_meta_16S$conc =="NA" ,]
  
  data_B_r_std <- data_B_r_std %>% filter(gene =="16S")
  data_B_r_std$Cq <- data_B_r_std$Cq %>% as.numeric()
  data_B_r_std$conc <- data_B_r_std$conc %>% as.numeric()
  # Exclude higher dilutions 
  data_B_r_std <- data_B_r_std %>% filter(conc >3) 
 
  efficiency_data <- data_B_r_std %>%
    summarise(
      model = list(lm(Cq ~ conc)),  # Store model
      slope = coef(model[[1]])[2],  # Extract slope
      intercept = coef(model[[1]])[1],
      efficiency = (10^(-1/slope) - 1) * 100,  # Calculate efficiency
      r_squared = summary(model[[1]])$r.squared,  # Extract R²
      .groups = "drop"
    ) 
  
  results_lm1 <- paste0("Slope: ", round(efficiency_data$slope,2),", Intercept: ", round(efficiency_data$intercept,1))
  results_lm2 <- paste0("Efficiency: ", round(efficiency_data$efficiency,1),"%, R-square: ", round(efficiency_data$r_squared, 3))
  
  
  g1 <-ggplot(data= data_B_r_std, aes(x = conc, y = Cq))+  geom_smooth(method ="lm") +geom_point()+
    ylab("Ct")+ xlab("Log10(gBlock 16S copies/μL)")+ theme_classic()+
    annotate("text", x = 4, y = 11, label =results_lm1, hjust =0)+
    annotate("text", x = 4, y = 10, label =results_lm2, hjust =0)
  
  #ggsave(g1, filename = "figures/qPCR/030625_16S_gradient_linear_red.pdf", height =5, width = 5)
  #ggsave(g1, filename = "figures/qPCR/030625_16S_gradient_linear_all_conc.pdf", height =5, width = 5) # Wihtou excluding high dilutions
  
  # Formula 10^((Cq-b)/m), still DNA copies still log10 transformed
  
  data_B_r$logcopy <- ((data_B_r$Cq-efficiency_data$intercept)/efficiency_data$slope)
  data_B_r$F.number <- data_B_r$sample
  
  
  data_B_r_meta <-merge(data_B_r %>% group_by(F.number) %>% dplyr::summarize(logcopy = mean(logcopy)),
                        design2 %>% select(F.number, ID, compost, rep, kit, size, dna_25), by = "F.number")
  
  # Normalize to 1 ng/mL
  data_B_r_meta$copies_per_ngDNA <-((10^data_B_r_meta$logcopy)*10)/data_B_r_meta$dna_25
  
  
  # Without adapting DNA concentration
  ggplot(data = data_B_r_meta, aes(x = size, y = copies_per_ngDNA, colour = compost, fill=compost)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    bg_theme+
    geom_line(aes(group = interaction(compost, size)), position = position_dodge(width = 0.5)) +
    facet_wrap(~kit, labeller = labeller(kit = c("M" = "NucleoMag", "S" = "NucleoSpin"))) +
    scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())+
    scale_color_manual(values = color.compost)+
    ylab("Copies 16S/ ng DNA")
  
  
  results <-ANOVA_summary(data = data_B_r_meta, x1 = "compost", x2 = "size", x3 = "kit", y= "copies_per_ngDNA")
  
  ggplot(data = data_B_r_meta, aes(x = size, y = copies_per_ngDNA, colour = compost, fill=compost)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = interaction(compost, size)), position = position_dodge(width = 0.5)) +
    facet_wrap(~kit, labeller = labeller(kit = c("M" = "NucleoMag", "S" = "NucleoSpin"))) +
    geom_text(
      data = results$data_sum, 
      aes(x = size, y = pos.letter, label = letters, fill=compost),  # Use `colour = compost` to match group
      position = position_dodge(width = 0.5),  # Apply consistent dodging
      vjust = -0.5,  # Place text slightly above points,
      color="darkred",
      size = 5, 
      inherit.aes = FALSE
    ) +
    ylab("Copies 16S/ ng DNA") +
    scale_color_manual(values = color.compost) +
    scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
    bg_theme +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())
  
  #ggsave(filename = "figures/qPCR/210525_16S_CQ_norm_1ng_qL.pdf", height = 5, width = 10)
  ANOVA <-aov(copies_per_ngDNA ~ compost*size*kit, data = data_B_r_meta)
  #pdf(file="figures/qPCR/ANOVA_copies_per_ng_DNA_plots_B.pdf", height =6, width = 6)
  par(mfrow=c(2,2))
  plot(ANOVA)
  #dev.off()
  
  #sink(file="output/qPCR/ANOVA_copy_per_ng_DNA_B.txt")
  summary(ANOVA)
  #sink()
  
  
  # Comparing bacteria and fungi
  par(mfrow =c(1,1))

  #data_B_r_meta$F.number == data_F_r_meta$F.number # Check the order
  
  plot(data_B_r_meta$copies_per_ngDNA, data_F_r_meta$copies_per_ngDNA)
  cor.test(data_B_r_meta$copies_per_ngDNA, data_F_r_meta$copies_per_ngDNA)

  # Correlation DNA conc and ITS/16S abundance
  
  #sink("output/DNA_concentration/Correlation_ITS_16S_with_dna_conc.txt")
  cor.test(data_B_r_meta$copies_per_ngDNA, data_B_r_meta$dna_25)
  cor.test(data_F_r_meta$copies_per_ngDNA, data_F_r_meta$dna_25)
  #sink()
  
  rm(g1, results, model, std_curve, efficiency_data, text, results_lm1, results_lm2,
     ANOVA)
  
  ## Merge files and save
  
  data_B_select <-data_B_r_meta %>% select(F.number, logcopy, copies_per_ngDNA)
  colnames(data_B_select) <- c("F.number","logcopy16s", "copies_per_ngDNA_16S")
  
  data_F_select <- data_F_r_meta
  colnames(data_F_select) <- c("F.number","logcopyITS2", "ID", "compost", "rep","kit", "size","dna_25","copies_per_ngDNA_ITS2")
  data_qPCR_meta <-left_join(data_F_select, data_B_select, by ="F.number")
  
  # Calcuate B & F ratio
  data_qPCR_meta$ratioF_B <- data_qPCR_meta$copies_per_ngDNA_ITS2/data_qPCR_meta$copies_per_ngDNA_16S
  
  results <-ANOVA_summary(data = data_qPCR_meta, x1 = "compost", x2 = "size", x3 = "kit", y= "ratioF_B")
  
  ggplot(data = data_qPCR_meta, aes(x = size, y = ratioF_B, colour = compost, fill=compost)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = interaction(compost, size)), position = position_dodge(width = 0.5)) +
    facet_wrap(~kit, labeller = labeller(kit = c("M" = "NucleoMag", "S" = "NucleoSpin"))) +
    geom_text(
      data = results$data_sum, 
      aes(x = size, y = pos.letter, label = letters, fill=compost),  # Use `colour = compost` to match group
      position = position_dodge(width = 0.5),  # Apply consistent dodging
      vjust = -0.5,  # Place text slightly above points,
      color="darkred",
      size = 5, 
      inherit.aes = FALSE
    ) +
    ylab("Copies ITS2/ Copies 16S")+
    scale_color_manual(values = color.compost) +
    scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
    bg_theme +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())
  
  #ggsave(filename = "figures/qPCR/Ratio_ITS_16S_copy.pdf", height = 5, width = 10)
  ANOVA <-aov(ratioF_B ~ compost*size*kit, data = data_qPCR_meta)
  #pdf(file="figures/qPCR/ANOVA_ratio_ITS_16S_copy.pdf", height =6, width = 6)
  par(mfrow=c(2,2))
  plot(ANOVA)
  #dev.off()
  
  #sink(file="output/qPCR/ANOVA_ratio_ITS_16S_copy.txt")
  summary(ANOVA)
  #sink()
  
  # Save file to data
  #write.csv(data_qPCR_meta, file="data/080725_qPCR_data_calc.csv", row.names = F )
  
  ## Normalize unrarfied sequence count----
  # Bacteria

  rownames(data_qPCR_meta) <- data_qPCR_meta$ID
  all(rownames(data_qPCR_meta) %in% colnames(asv.B)) # Check if the rownames are the same
  
  # Bacteria
  asv.B.rel <- sweep(asv.B, 2, colSums(asv.B), FUN = "/")
  copies_16S <- data_qPCR_meta[colnames(asv.B), "copies_per_ngDNA_16S"]
  asv.B.absolut <- sweep(asv.B.rel, 2, copies_16S, FUN = "*")
  #write.table(asv.B.absolut, file ="data/asv.B.absolut.txt")
  
  # Fungi
  asv.F.rel <- sweep(asv.F, 2, colSums(asv.F), FUN = "/")
  copies_ITS2<- data_qPCR_meta[colnames(asv.F), "copies_per_ngDNA_ITS2"]
  asv.F.absolut <- sweep(asv.F.rel, 2, copies_ITS2, FUN = "*")
  #write.table(asv.F.absolut, file ="data/asv.F.absolut.txt")

  
  
  