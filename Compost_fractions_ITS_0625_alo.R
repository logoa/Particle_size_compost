# Statistical analysis and visualizations for paper "Microbial community analyses of composts are influenced by
# particle size fraction and DNA extraction method"

# Analysis of bacterial and fungal communities
# Started 25.06.25
# Last change 02.12.25
# Author: Anja Logo
# Agroscope, Switzerland

# Summary of approach:

# For alpha and beta diversity rarefied data is used (lowest read number per sample)
# For taxonomic composition and Differential abundance analysis unrarefied data + qPCR normalization is used
# Differential abundance analysis: only taxa which have been found in 2/3 replicates, relative abundance cut-off 0.01\%

# Load data, packages, functions-------
source(file = "Comp_size_25_setup.R")

#sink(file="Session_Info_021225.txt")
#sessionInfo()
#sink()

  ## Functions----
  
    # ANOVA with plot
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
  
  # Function to draw beta diversity plot
    PcoA_beta <- function(distance_matrix, output_file, save =TRUE) {
    # Perform PCoA
    pco_result <- cmdscale(distance_matrix, eig = TRUE, add = TRUE)
    
    # Calculate PCO1 and PCO2 percentages
    pco1 <- (pco_result$eig / sum(pco_result$eig))[1]
    pco2 <- (pco_result$eig / sum(pco_result$eig))[2]
    pco1_label <- paste("PCO1 [", format(round(pco1 * 100, 1), nsmall = 1), "%]", sep = "")
    pco2_label <- paste("PCO2 [", format(round(pco2 * 100, 1), nsmall = 1), "%]", sep = "")
    
    # Extract PCoA points
    pco_points <- pco_result$points
    colnames(pco_points) <- c("pco1", "pco2")
    
    # Merge with design metadata
    design_pcoa <- merge(design.alpha, pco_points, by = 0)
    rownames(design_pcoa) <- design_pcoa$Row.names
    design_pcoa$Row.names <- NULL
    
    # Create hulls for polygons
    find_hull <- function(df) df[chull(df$pco1, df$pco2), ]
    hulls <- ddply(design_pcoa, "compost", find_hull)
    
    # Plot
    p <- ggplot(data = design_pcoa, aes(x = pco1, y = pco2, colour = size, shape = kit)) +
      geom_polygon(data = hulls, alpha = 0.8, aes(x = pco1, y = pco2, fill = compost), inherit.aes = FALSE) +
      geom_point(size = 2) +
      scale_colour_manual(values = c("blue", "orange"), name = "Particle fraction", labels = c("fine", "coarse")) +
      theme_classic(base_size = 14) +
      coord_fixed() +
      ylab(pco2_label) + xlab(pco1_label) +
      scale_fill_manual(values = color.compost, name = "Compost") +
      scale_shape_manual(values = c(15,17), name = "DNA extraction", labels = c("Bead-based", "Column-based")) +
      theme(legend.position = "none")
    
    if (save ==TRUE) {
      ggsave(output_file, plot = p, height = 6, width = 8)
    } 
    return(p)
    
    # Save plot
  
  }
  
  # Random color generator
    generate_random_color <- function() {
    rgb(runif(1), runif(1), runif(1))
  } # Function to generate a random RGB color
  
  # Generate distinct colors that can be differentiate better than random colors
    generate_distinct_colors <- function(n) {
      colors <- vector("character", n)  # Initialize an empty vector for colors
      
      # Generate n distinct colors
      for (i in 1:n) {
        # Generate a random color
        new_color <- rgb(runif(1), runif(1), runif(1) * 0.8)  # RGB with some constraints
        # Check for uniqueness
        while (any(colors == new_color)) {
          new_color <- rgb(runif(1), runif(1), runif(1) * 0.8)
        }
        colors[i] <- new_color
      }
  
      return(colors)
    }
  
  # Draw venn diagram function
    draw.venn.diagram = function(data, venn.factor, title, custom_colors) {
    data.meta = merge(design, data, by = 0)
    rownames(data.meta) = data.meta$Row.names
    data.meta$Row.names <- NULL 
    factor = levels(as.factor(data.meta[,venn.factor]))
    asv_tables<- lapply(factor, function(method) {
      subset_data <- data.meta[data.meta[, venn.factor] == method, ]
      asv_table_subset <- subset_data %>% select(contains("ASV"))  # Remove metadata columns
      return(asv_table_subset)
    })
    names(asv_tables) <- factor
    
    discard_missing_asvs <- function(asv_table) {
      present_columns <- apply(asv_table, 2, function(col) any(col > 0))
      return(asv_table[, present_columns])
    }
    
    # Apply the function to each site's ASV table
    asv_tables_filtered <- lapply(asv_tables, discard_missing_asvs)
    
    # Extract ASV names from each site
    asv_names <- lapply(asv_tables_filtered, colnames)
    # Create a list of sets for the Venn diagram
    
    
    sets <- lapply(asv_names, function(x) as.character(x))
    venn.plot <- VennDiagram::venn.diagram(
      x = sets,
      category.names = names(asv_tables_filtered),
      filename = NULL,
      output = TRUE,
      fill = custom_colors, 
      main = title
      # Adjust the position of category names
    )
    # Display the Venn diagram without colors
    grid::grid.newpage()
    grid::grid.draw(grid::grobTree(grobs = venn.plot))
    
    return(sets)
    
  }
  
  # Function to search for indicators
  
    ica_levels <- function(rar_abund, tax_file, tax_level, factor, compost_list, microbe, min_rel= 0.1){
    
    if (tax_level != "ASV") {
      formula <- reformulate(tax_level, response = "as.matrix(rar_abund)")
      rar_abund_abb <- aggregate(formula , data = tax_file[rownames(rar_abund),], FUN = sum)  # aggregate by phylum
      rownames(rar_abund_abb) <-rar_abund_abb[,tax_level]; rar_abund_abb[,tax_level] <-NULL
    } else{rar_abund_abb <- rar_abund}
    design.rar_abund <-merge(t(rar_abund_abb), design[,c(factor, "compost")], by =0)
    rownames(design.rar_abund) <-design.rar_abund$Row.names; design.rar_abund$Row.names <-NULL
    if (compost_list == "all") {
      data <- design.rar_abund
    }else{data <- design.rar_abund %>% dplyr::filter(compost == compost_list)}
    
    groups <-as.character(data[,factor])
    data[,c(factor, "compost")] <-NULL
    
    # Filter rare ASVs
    min.read <-data %>% rowSums() %>% mean()/(100/min_rel) # Number of sequences for a minimum read count of 0.1% in this samples
    data <-data[,data %>% colSums() > min.read]
    print(paste0("Analaysis: ", microbe, "; ", factor, "; ", tax_level, "; Compost ", compost_list,"; Number of ", tax_level, " ",ncol(data)))
    
    # Point biserial correltion
    indicators <- multipatt(data, groups, func = "r.g", control = how(nperm = 999))
    output <- indicators$sign
    # BH correction for multiple testing
    output$p.adj<- p.adjust(output$p.value, method = "BH")
    output$significance <- ifelse(output$p.adj < 0.05, "Significant", "Not Significant")
    
    # Group
    group <-indicators$cluster %>% unique()
    output$group <- ifelse(output$index ==1,  group[1],  group[2]) # Classify the results
    
    if (tax_level == "ASV") {
      output_tax <-merge(output,tax_file[,c("Phyla", "Family", "Genus")], by = 0)
      rownames(output_tax) <- output_tax$Row.names; output_tax$Row.names <- NULL
      
      # Calculate the abundance in the composts
      output_tax$count <-rar_abund[rownames(output_tax),rownames(data)] %>% rowSums()
      output_tax$count_rel <- output_tax$count / (rar_abund[rownames(output_tax),rownames(data)] %>% colSums() %>% mean())
      output_tax$count_rel_log <- output_tax$count_rel %>% log10()
    } else {output_tax <- output}
    
    write.csv(output_tax, file = paste0("output/Indicator_analysis/", factor, "/", microbe, tax_level, "/",
                                        microbe, tax_level, "_",compost_list, "_", factor, "_PBC.csv"))
  }
  
  # Differential abundance function
    aldex_levels <-function(asv_table, tax_file, tax_level, factor, compost_list, microbe, min_rel = 0.1){
    
    # aggregate by taxonomic classification
    if (tax_level != "ASV") {
      formula <- reformulate(tax_level, response = "as.matrix(asv_table)")
      asv.F.phyla <- aggregate(formula , data = tax_file[rownames(asv_table),], FUN = sum)  
      
      rownames(asv.F.phyla) <-asv.F.phyla[, tax_level]; asv.F.phyla[, tax_level] <-NULL # HERE
    } else {asv.F.phyla <- asv_table}
    
    # subset only one compost
    if (compost_list == "all") {
      asv.F.phyla.sub <- asv.F.phyla
    } else{ 
      asv.F.phyla.sub <- asv.F.phyla[, grepl(compost_list, colnames(asv.F.phyla))]
    }
    # Define the levels
    if (factor == "size") {
      conds <- ifelse(grepl("10", colnames(asv.F.phyla.sub)), "10", "2")
    } else {
      conds <- ifelse(grepl("M", colnames(asv.F.phyla.sub)), "M", "S") }
    
    # Filter rare taxa
    min.read <-asv.F.phyla.sub %>% colSums() %>% mean()/(100/min_rel) # Number of sequences for a minimum read count of 0.1% 
    asv.F.phyla.sub <-asv.F.phyla.sub[asv.F.phyla.sub %>% rowSums() > min.read,]
    print(paste0("Analysis: ", microbe, "; ", factor, "; ", tax_level, "; Compost ", compost_list,"; Number of ", tax_level, " ", nrow(asv.F.phyla.sub)))
    
    # Performe the difference abundance analysis
    x <- aldex.clr(reads = asv.F.phyla.sub, conds = conds, mc.samples = 128, denom = "all", verbose = FALSE)
    aldex_res <- aldex.ttest(x, paired.test = FALSE, verbose =FALSE) # t.tests
    aldex_effect <-aldex.effect(x, CI =TRUE, verbose =FALSE) # calculate the effect size
    aldex_out <- data.frame(aldex_res, aldex_effect)
    
    groups <- unique(x@conds) # Extract the unique levels of the conds vector
    group1 <- groups[1]  # Reference group (negative effect size)
    group2 <- groups[2]  # Comparison group (positive effect size)
    
    aldex_out$group <- ifelse(aldex_out$effect > 0, group2, group1) # Classify the results
    aldex_out$significance <- ifelse(aldex_out$wi.eBH < 0.05, "Significant", "Not Significant")
    
    if (tax_level =="ASV") {
      aldex_tax <-merge(aldex_out,tax_file[,c("Phyla", "Family", "Genus")], by = 0)
      rownames(aldex_tax) <- aldex_tax$Row.names; aldex_tax$Row.names <- NULL
    } else{aldex_tax <- aldex_out}
    
    aldex_tax$count <-asv_table[rownames(aldex_out),colnames(asv.F.phyla.sub)] %>% rowSums()
    aldex_tax$count_rel <- aldex_tax$count / (asv.B[rownames(aldex_out),colnames(asv.F.phyla.sub)] %>% colSums() %>% mean())
    aldex_tax$count_rel_log <- aldex_tax$count_rel %>% log10()
    
    write.csv(aldex_tax, file = paste0("output/Indicator_analysis/", factor, "/", microbe, tax_level, "/",
                                       microbe, tax_level, "_",compost_list, "_", factor, "_ALDEX.csv"))
  
  }
  
  # Differential abundance analysis paired
    
    #asv_table <- asv.B.absolut
    #tax_file  <- tax.B
    #tax_level <- "ASV"
    #factor <- "size"
    #compost_list <-"K22"
    #microbe <-"B"
    #sub_factor <- "M"
    
    aldex_levels_paired <-function(asv_table, tax_file, tax_level, factor, compost_list, microbe, min_rel = 0.0001, sub_factor){
    
    # aggregate by taxonomic classification
    if (tax_level != "ASV") {
      formula <- reformulate(tax_level, response = "as.matrix(asv_table)")
      asv.F.phyla <- aggregate(formula , data = tax_file[rownames(asv_table),], FUN = sum)  
      
      rownames(asv.F.phyla) <-asv.F.phyla[, tax_level]; asv.F.phyla[, tax_level] <-NULL # HERE
    } else {asv.F.phyla <- asv_table}
    
    # subset only one compost
    if (compost_list == "all") {
      asv.F.phyla.sub <- asv.F.phyla
    } else{ 
      asv.F.phyla.sub <- asv.F.phyla[, grepl(compost_list, colnames(asv.F.phyla))]
    }
    
    # kit or size & condition levels
    if (factor == "size" & sub_factor =="M") {
      asv.F.phyla.sub <- asv.F.phyla.sub[, grep("M", colnames(asv.F.phyla.sub))]
      conds <- ifelse(grepl("10", colnames(asv.F.phyla.sub)), "10", "2")
    } 
    
    if(factor == "size" & sub_factor =="S"){
      asv.F.phyla.sub <- asv.F.phyla.sub[, grep("S", colnames(asv.F.phyla.sub))]
      conds <- ifelse(grepl("10", colnames(asv.F.phyla.sub)), "10", "2")
    }
    
    if (factor =="kit" & sub_factor =="2") {
      asv.F.phyla.sub <- asv.F.phyla.sub[, !grepl("10", colnames(asv.F.phyla.sub))]
      conds <- ifelse(grepl("S", colnames(asv.F.phyla.sub)), "S", "M")
    }
    
    if (factor == "kit" & sub_factor =="10") {
      asv.F.phyla.sub <- asv.F.phyla.sub[, !grepl("10", colnames(asv.F.phyla.sub))]
      conds <- ifelse(grepl("S", colnames(asv.F.phyla.sub)), "S", "M")
     }
    
    print(table(conds))
    try(if(exists("conds") ==FALSE) stop("Factor and sub_factor are not valid"))
    
    # Create replicate ID by stripping treatment info from column names
    replicate_ids <- substr(colnames(asv.F.phyla.sub), 1,5) 
  
    # Filter raretaxa which have are not more abundant than 0.01% in at least one of the samples included in the analysis!
    rel_abund <- sweep(asv.F.phyla.sub, 2, colSums(asv.F.phyla.sub), "/")  # relative abundances
    keep_taxa <- apply(rel_abund, 1, function(x) any(x >= 0.0001))  # ≥0.01% in at least one sample
    asv.F.phyla.sub <- asv.F.phyla.sub[keep_taxa,]
  
    print(paste0("Analysis: ", microbe, "; ", factor, "; ", sub_factor,";", tax_level, "; Compost ", compost_list,"; Number of ", tax_level, ": ", nrow(asv.F.phyla.sub)))
    
      # Perform ALDEx2 analysis with paired test
    x <- aldex.clr(reads = asv.F.phyla.sub, conds = conds, mc.samples = 128, denom = "all", verbose = FALSE)
    aldex_res <- aldex.ttest(x, paired.test = TRUE, verbose = FALSE)
    aldex_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)
    aldex_out <- data.frame(aldex_res, aldex_effect)
    
    groups <- unique(x@conds) # Extract the unique levels of the conds vector
    group1 <- groups[1]  # Reference group (positive effect size)
    group2 <- groups[2]  # Comparison group (negative effect size)
    
    print(groups)
    
    aldex_out$group <- ifelse(aldex_out$effect > 0, group1, group2) # Classify the results
    aldex_out$significance <- ifelse(aldex_out$wi.eBH < 0.05 & abs(aldex_out$effect) > 1, "Significant", "Not Significant")
    
    hist(aldex_out$effect, breaks = 50, main = "Effect size distribution", xlab = "Effect size")
    
    if (tax_level =="ASV") {
      aldex_tax <-merge(aldex_out,tax_file[,c("Phyla", "Family", "Genus")], by = 0)
      rownames(aldex_tax) <- aldex_tax$Row.names; aldex_tax$Row.names <- NULL
      aldex_tax$count <-asv_table[rownames(aldex_out),colnames(asv.F.phyla.sub)] %>% rowSums()
      aldex_tax$count_rel <- aldex_tax$count / (asv.B[rownames(aldex_out),colnames(asv.F.phyla.sub)] %>% colSums() %>% mean())
      aldex_tax$count_rel_log <- aldex_tax$count_rel %>% log10()
    } else{aldex_tax <- aldex_out}
    
     write.csv(aldex_tax, file = paste0("output/Indicator_analysis_ab_ab/", factor, "/ALDEX_",microbe, "_", tax_level, "_",compost_list, "_", factor, "_", sub_factor, ".csv"))
    return(aldex_tax)
    }
    
    # Function to run Wilcoxon signed-rank test for each ASV
    consistency_levels_paired <- function(asv_table, tax_file, tax_level, factor, compost_list, microbe, min_rel = 0.0001, sub_factor) {
      
      # Aggregate by taxonomic classification if needed
      if (tax_level != "ASV") {
        formula <- reformulate(tax_level, response = "as.matrix(asv_table)")
        asv.F.phyla <- aggregate(formula, data = tax_file[rownames(asv_table), ], FUN = sum)
        rownames(asv.F.phyla) <- asv.F.phyla[, tax_level]; asv.F.phyla[, tax_level] <- NULL
      } else {
        asv.F.phyla <- asv_table
      }
        asv.F.phyla.sub <- asv.F.phyla[, grepl(compost_list, colnames(asv.F.phyla))]
  
      
      # Define conditions based on factor and sub_factor
      if (factor == "size" & sub_factor == "M") {
        asv.F.phyla.sub <- asv.F.phyla.sub[, grepl("M", colnames(asv.F.phyla.sub))]
        conds <- ifelse(grepl("10", colnames(asv.F.phyla.sub)), "10", "2")
      }
      
      if (factor == "size" & sub_factor == "S") {
        asv.F.phyla.sub <- asv.F.phyla.sub[, grepl("S", colnames(asv.F.phyla.sub))]
        conds <- ifelse(grepl("10", colnames(asv.F.phyla.sub)), "10", "2")
      }
      
      if (factor == "kit" & sub_factor == "2") {
        asv.F.phyla.sub <- asv.F.phyla.sub[, !grepl("10", colnames(asv.F.phyla.sub))]
        conds <- ifelse(grepl("S", colnames(asv.F.phyla.sub)), "S", "M")
      }
      
      if (factor == "kit" & sub_factor == "10") {
        asv.F.phyla.sub <- asv.F.phyla.sub[, grepl("10", colnames(asv.F.phyla.sub))]
        conds <- ifelse(grepl("S", colnames(asv.F.phyla.sub)), "S", "M")
      }
      
      if (!exists("conds")) stop("Factor and sub_factor are not valid")
      
      print(conds)
      # Replicate IDs (first 5 characters of sample names assumed to indicate replicate)
      replicate_ids <- substr(colnames(asv.F.phyla.sub), 1, 5)
      
      # Filter rare taxa
      rel_abund <- sweep(asv.F.phyla.sub, 2, colSums(asv.F.phyla.sub), "/")
      keep_taxa <- apply(rel_abund, 1, function(x) any(x >= min_rel))
      asv.F.phyla.sub <- asv.F.phyla.sub[keep_taxa, ]
      
      print(paste0("Analysis: ", microbe, "; ", factor, "; ", sub_factor, "; ", tax_level, "; Compost ", compost_list, "; Number of ", tax_level, ": ", nrow(asv.F.phyla.sub)))
      
      print(colnames(asv.F.phyla.sub))
      
      # Run Wilcoxon signed-rank test
      sample_order <- order(replicate_ids, conds)
      asv.F.phyla.sub <- asv.F.phyla.sub[, sample_order]
      conds <- conds[sample_order]
      replicate_ids <- replicate_ids[sample_order]
      
      group_levels <- unique(conds)
      group1_idx <- which(conds == group_levels[1])
      group2_idx <- which(conds == group_levels[2])
      
      if (length(group1_idx) != length(group2_idx)) {
        stop("Unequal sample sizes in paired design.")
      }
      
      consistency_results <- apply(asv.F.phyla.sub, 1, function(asv_row) {
        diffs <- asv_row[group2_idx] - asv_row[group1_idx]
        
        all_positive <- all(diffs > 0)
        all_negative <- all(diffs < 0)
        
        consistency_flag <- if (all_positive) {
          "Consistently Increased"
        } else if (all_negative) {
          "Consistently Decreased"
        } else {
          "Inconsistent Change"
        }
        
        # Median of absolute differences
        median_diff <- median(diffs)
        
        # Median log2 fold-change (with pseudo-count of 1 to avoid division by zero)
        log2fc_samples <- log2((asv_row[group2_idx] + 1) / (asv_row[group1_idx] + 1))
        log2fc <- median(log2fc_samples)
        
        c(consistency_flag = consistency_flag, 
          median_diff = median_diff, 
          log2FC = log2fc)
      })
      
      consistency_results <- as.data.frame(t(consistency_results))
      consistency_results$median_diff <- as.numeric(as.character(consistency_results$median_diff))
      consistency_results$log2FC <- as.numeric(as.character(consistency_results$log2FC))
      
      if (tax_level == "ASV") {
        consistency_results <- merge(consistency_results, tax_file[, c("Phyla", "Family", "Genus")], by = 0)
        rownames(consistency_results) <- consistency_results$Row.names
        consistency_results$Row.names <- NULL
        consistency_results$count <- rowMeans(asv_table[rownames(consistency_results), colnames(asv.F.phyla.sub)])
      }
      
      consistency_results <- consistency_results %>% filter(consistency_flag != "Inconsistent Change")
      #write.csv(consistency_results, file = paste0("output/Indicator_analysis_ab_ab/", factor, "/Consistency_", microbe, "_", tax_level, "_", compost_list, "_", factor, "_", sub_factor, ".csv"))
      
      return(consistency_results)
      
    }
    
    # Compoare across composts which ASVs are enriched
    differ_analysis_plot <-function(res1, res2){
      
      d1 <-results[[res1]] %>% filter(abs(median_diff) >=10 & abs(log2FC) >=1)
      d2 <-results[[res2]] %>% filter(abs(median_diff) >=10 & abs(log2FC) >=1)
      
      factor1 <- intersect(rownames(d1[d1$log2FC >0,]), rownames(d2[d2$log2FC >0,]))
      tax.B[factor1,] %>% pull(Phyla) %>% table() %>% sort()
      factor2 <- intersect(rownames(d1[d1$log2FC <0,]), rownames(d2[d2$log2FC <0,]))
      tax.B[factor2,] %>% pull(Phyla) %>% table() %>% sort()
      
      plot_data1 <-d1[c(factor1, factor2),]
      plot_data2 <-d2[c(factor1, factor2),c("median_diff", "log2FC", "count")]
      
      colnames(plot_data2) <-c("median_diff_2", "log2FC2", "count2")
      
      plot_data  <-plot_data1 %>% merge(plot_data2, by =0)
      rownames(plot_data) <- plot_data$Row.names; plot_data$Row.names <-NULL
      plot_data$log2FCmean <- log2((2^plot_data$log2FC + 2^plot_data$log2FC2)/2)
      plot_data$countmean <- (plot_data$count + plot_data$count2)/2
  
      return(plot_data)
    }
    
    # for mixed linear models
    fit_and_report <- function(data, response, fixed_formula, random_effect) {
      
      # Define formulas
      formula_mixed <- as.formula(paste(response, "~", fixed_formula, "+ (1 |", random_effect, ")"))
      formula_fixed <- as.formula(paste(response, "~", fixed_formula))
      
      # Fit models
      model_mixed <- lmerTest::lmer(formula_mixed, data = data, REML = FALSE)
      singular <- isSingular(model_mixed)
      
      if (!singular) {
        message("Random effect meaningful: using lmer")
        final_model <- model_mixed
      } else {
        message("Random effect singular: using lm")
        final_model <- lm(formula_fixed, data = data)
      }
      
      return(final_model)
    }
    
  ## Data------
    ### General-----
    bg_theme <- theme_classic() +
      theme(
        legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.3, vjust = 0, size = 14, color ="black"),
        axis.text.y = element_text(size = 14, color ="black"),
        text = element_text(size = 14, color ="black"),
        axis.title = element_text(size=14),
        legend.text = element_text(size=14)
      )
    color.compost <- c("grey80", "grey40", "grey10")
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
    
    design.all <- read.csv(file ="data/design_all_seq_alpha_quanty.csv")
    rownames(design.all) <- design.all$X; design.all$X <-NULL
    
    color.genus <- c("lightblue", "#8D8F29", "#FFB7B4", "#453020", "#916195", "#FDC817",  "#498A57", "#1BA29F", "indianred3",
                            "#403368", "red4" ,"#1E6640", "darkorange", "lightsalmon",  "blue4" ,"#45FF8A" ,"#544C40" ,"#A2217A" ,"#76AD24" ,"moccasin")
    ### Sequencing data------
  
  ## FUNGI
    asv.F.all <-read.table(file= "data/asv.all.F.txt") # raw
    asv.F <-read.table(file="data/asv.F.txt") # unrarefied only fungi
    ISS.F <-read.table(file ="data/ISS.F.txt") # rarefied
    ISS.rob.F <- read.table(file= "data/ISS.rob.F.txt") # Only robustly deteceted
    tax.F <-read.table(file="data/tax.F.txt") # Unrarefied only fungi
    
    # BACTERIA
    asv.B.all <-read.table(file= "data/asv.all.B.txt") # raw
    asv.B <-read.table(file="data/asv.B.txt") # unrarefied only bacteria
    ISS.B <-read.table(file ="data/ISS.B.txt") # Unrarefied only bacteria
    ISS.rob.B <- read.table(file= "data/ISS.rob.B.txt") # only robustly detected
    tax.B <-read.table(file="data/tax.B.txt") # unrarefied only bacteria
    
    # Alpha diversity
    alpha <- read.csv(file="data/ISS.alpha.F.B.csv", row.names = 1)
    design.alpha <- merge(design, alpha, by=0)
    rownames(design.alpha)<- design.alpha$Row.names; design.alpha$Row.names<-NULL
    design.alpha$size <- as.factor(design.alpha$size)
    design.alpha$compost <- as.factor(design.alpha$compost)
    design.alpha$kit <- as.factor(design.alpha$kit)
    rm(alpha)
    
    # Beta diversity
    ISS.bray.F <- read.table(file="data/ISS.bray.F.txt")  %>% as.matrix() %>% as.dist()
    ISS.bray.sqrt.F <- read.table(file= "data/ISS.bray.sqrt.F.txt")  %>% as.matrix() %>% as.dist()
    ISS.jac.F <- read.table(file="data/ISS.jac.F.txt")  %>% as.matrix() %>% as.dist()
    ISS.euc.F <- read.table(file="data/ISS.euc.F.txt") %>% as.matrix() %>% as.dist()
    
    ISS.bray.B <- read.table(file= "data/ISS.bray.B.txt") %>% as.matrix() %>% as.dist()
    ISS.bray.sqrt.B <- read.table(file= "data/ISS.bray.sqrt.B.txt") %>% as.matrix() %>% as.dist()
    ISS.jac.B <- read.table(file= "data/ISS.jac.B.txt")  %>% as.matrix() %>% as.dist()
    ISS.euc.B <- read.table(file="data/ISS.euc.B.txt")  %>% as.matrix() %>% as.dist()
    
    # Summarizing the most common phyla
    
    # FUNGI
    ISS.prop <- prop.table(as.matrix(ISS.F), margin = 2) * 100  # get proportions
    ISS.p <- aggregate(ISS.prop ~ Phyla, data = tax.F[rownames(ISS.F),], FUN = sum)  # aggregate by phylum
    rownames(ISS.p) <- ISS.p$Phyla; ISS.p$Phyla <- NULL
    phyla.sum <-(rowSums(ISS.p)/ncol(ISS.p)) %>% sort() %>% as.data.frame()
    colnames(phyla.sum) <-c("summ_ab")
    phyla.ab <-phyla.sum %>% filter(summ_ab > 0.5)
    
    tax.F.r <- tax.F[, c("Phyla", "Family", "Genus")]
    
    for (i in 1: nrow(tax.F.r)) {
      if(tax.F[i,]$Phyla %in% rownames(phyla.ab)){
        tax.F.r[i,]$Phyla <-tax.F[i,]$Phyla
      } else{tax.F.r[i,]$Phyla <- "other"}
    }
    
    # BACTERIA
    ISS.prop <- prop.table(as.matrix(ISS.B), margin = 2) * 100  # get proportions
    ISS.p <- aggregate(ISS.prop ~ Phyla, data = tax.B[rownames(ISS.B),], FUN = sum)  # aggregate by phylum
    rownames(ISS.p) <- ISS.p$Phyla; ISS.p$Phyla <- NULL
    phyla.sum <-(rowSums(ISS.p)/ncol(ISS.p)) %>% sort() %>% as.data.frame()
    colnames(phyla.sum) <-c("summ_ab")
    phyla.ab <-phyla.sum %>% filter(summ_ab > 0.5)
    tax.B.r <- tax.B[, c("Phyla", "Family", "Genus")]
    
    for (i in 1: nrow(tax.B.r)) {
      if(tax.B[i,]$Phyla %in% rownames(phyla.ab)){
        tax.B.r[i,]$Phyla <-tax.B[i,]$Phyla
      } else{tax.B.r[i,]$Phyla <- "other"}
    }
    rm(ISS.prop, ISS.p, phyla.sum, phyla.ab,i)
    
    # qPCR ITS & 16S
    
      asv.B.absolut <- read.table(file ="data/asv.B.absolut.txt")
      asv.F.absolut <- read.table(file ="data/asv.F.absolut.txt")
  
      data_qpcr_meta <-read.csv(file ="data/080725_qPCR_data_calc.csv")
      rownames(data_qpcr_meta) <- data_qpcr_meta$ID
  
# Number of sequences per sample (Figure S2)----

  num_seq <-asv.F %>% t() %>% rowSums() %>% as.data.frame()
  colnames(num_seq) <- "num_seq.F"
  design.nseq <- merge(design, num_seq , by =0)
  rownames(design.nseq) <- design.nseq$Row.names
  design.nseq$Row.names <- NULL
  num_seq <-asv.B %>% t() %>% rowSums() %>% as.data.frame()
  colnames(num_seq) <- "num_seq.B"
  design.nseq <- merge(design.nseq, num_seq , by =0)
  rownames(design.nseq) <- design.nseq$Row.names
  design.nseq$Row.names <- NULL
  
  design.nseq$size <- as.factor(design.nseq$size)
  design.nseq$compost <- as.factor(design.nseq$compost)
  design.nseq$kit <- as.factor(design.nseq$kit)
  rm(num_seq)

  ## FUNGI------ 
  ANOVA_num <- aov(num_seq.F  ~size*compost*kit, data = design.nseq)

  #sink(file="output/Number_of_seq_sample/ANOVA_number_seq_F.txt")
  summary(ANOVA_num)
  #sink()

  #pdf(file="figures/ANOVA_number_seq_residual_plots_F.pdf", height =6, width = 6)
  par(mfrow=c(2,2))
  plot(ANOVA_num)
  #dev.off()
  
  results <-ANOVA_summary(data = design.nseq, x1 = "compost", x2 = "size", x3 = "kit", y= "num_seq.F")
  
  ggplot(data = design.nseq, aes(x = size, y = num_seq.F, colour = compost, fill=compost)) +
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
    ylab("Number of sequences")+
    scale_color_manual(values = color.compost) +
    scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
    bg_theme +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())
  
  
  #ggsave(filename = "figures/Number_of_seq_sample/Comp_particle_number_seq_F.pdf", height = 5, width = 10)
  rm(ANOVA_num)
  
  ## BACTERIA------
  ANOVA_num <- aov(num_seq.B  ~size*compost*kit, data = design.nseq)
  
  #sink(file="output/Number_of_seq_sample/ANOVA_number_seq_B.txt")
  summary(ANOVA_num)
  #sink()
  
  #pdf(file="figures/Number_of_seq_sample/ANOVA_number_seq_residual_plots_B.pdf", height =6, width = 6)
  par(mfrow=c(2,2))
  plot(ANOVA_num)
  #dev.off()
  
  results <-ANOVA_summary(data = design.nseq, x1 = "compost", x2 = "size", x3 = "kit", y= "num_seq.B")
  
  ggplot(data = design.nseq, aes(x = size, y = num_seq.B, shape = kit)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = interaction(size, kit)), position = position_dodge(width = 0.5)) +
    facet_wrap(~compost) +
    geom_text(
      data = results$data_sum, 
      aes(x = size, y = pos.letter, label = letters, fill=kit),  # Use `colour = compost` to match group
      position = position_dodge(width = 0.5),  # Apply consistent dodging
      vjust = -0.5,  # Place text slightly above points,
      color="darkred",
      size = 5, 
      inherit.aes = FALSE
    ) +
    ylab("Number of sequences")+
    scale_color_manual(values = color.compost) +
    scale_y_continuous(labels = scales::label_scientific(digits = 1))+
    scale_shape_discrete(labels = c("NucleoMag", "NucleoSpin"), name ="DNA extraction")+
    scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
    bg_theme +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank())
  
  
  #ggsave(filename = "figures/Number_of_seq_sample/Comp_particle_number_seq_B.pdf", height = 5, width = 10)

# Sequencing overview (Table S3)------
  
  # Create an overview of number of sequences, ASV with ASV filtering
   
  seq_overview <-matrix(data =NA, ncol =7, nrow =12) %>% as.data.frame()
  colnames(seq_overview) <- c("Group", "Step", "Type", "Total", "Min", "Max", "Mean")
  seq_overview$Group <- c(rep("Bacteria",6), rep("Fungi",6))
  seq_overview$Step <-c(rep(c("High-quality", "Only bacteria", "Rarefied"),2),
                        rep(c("High-quality", "Only fungi", "Rarefied"),2))
  seq_overview$Type <- c(rep(c(rep("Sequences",3), rep("ASVs",3)),2))
  
  seq_overview[1, "Total"] <-asv.B.all %>% colSums() %>% sum()
  seq_overview[1, "Min"] <-asv.B.all %>% colSums() %>% min()
  seq_overview[1, "Max"] <-asv.B.all %>% colSums() %>% max()
  seq_overview[1, "Mean"] <-asv.B.all %>% colSums() %>% mean() %>% round()
  
  seq_overview[2, "Total"] <-asv.B %>% colSums() %>% sum()
  seq_overview[2, "Min"] <-asv.B %>% colSums() %>% min()
  seq_overview[2, "Max"] <-asv.B %>% colSums() %>% max()
  seq_overview[2, "Mean"] <-asv.B %>% colSums() %>% mean() %>% round()
  
  seq_overview[3, "Total"] <-ISS.B %>% colSums() %>% sum()
  seq_overview[3, "Min"] <-ISS.B %>% colSums() %>% min()
  seq_overview[3, "Max"] <-ISS.B %>% colSums() %>% max()
  seq_overview[3, "Mean"] <-ISS.B %>% colSums() %>% mean() %>% round()
  
  seq_overview[4, "Total"] <-asv.B.all %>% nrow()
  seq_overview[4, "Min"] <-ifelse(asv.B.all >0,1,0) %>% colSums() %>% min()
  seq_overview[4, "Max"] <-ifelse(asv.B.all >0,1,0) %>% colSums() %>% max()
  seq_overview[4, "Mean"] <-ifelse(asv.B.all >0,1,0)  %>% colSums() %>% mean() %>% round()
  
  seq_overview[5, "Total"] <-asv.B %>% nrow()
  seq_overview[5, "Min"] <-ifelse(asv.B>0,1,0) %>% colSums() %>% min()
  seq_overview[5, "Max"] <-ifelse(asv.B>0,1,0) %>% colSums() %>% max()
  seq_overview[5, "Mean"] <-ifelse(asv.B >0,1,0)  %>% colSums() %>% mean() %>% round()
  
  seq_overview[6, "Total"] <-ISS.B[rowSums(ISS.B)!=0,] %>% nrow()
  seq_overview[6, "Min"] <-ifelse(ISS.B>0,1,0) %>% colSums() %>% min()
  seq_overview[6, "Max"] <-ifelse(ISS.B>0,1,0) %>% colSums() %>% max()
  seq_overview[6, "Mean"] <-ifelse(ISS.B >0,1,0)  %>% colSums() %>% mean() %>% round()
  
  seq_overview[7, "Total"] <-asv.F.all %>% colSums() %>% sum()
  seq_overview[7, "Min"] <-asv.F.all %>% colSums() %>% min()
  seq_overview[7, "Max"] <-asv.F.all %>% colSums() %>% max()
  seq_overview[7, "Mean"] <-asv.F.all %>% colSums() %>% mean() %>% round()
  
  seq_overview[8, "Total"] <-asv.F %>% colSums() %>% sum()
  seq_overview[8, "Min"] <-asv.F %>% colSums() %>% min()
  seq_overview[8, "Max"] <-asv.F %>% colSums() %>% max()
  seq_overview[8, "Mean"] <-asv.F %>% colSums() %>% mean() %>% round()
  
  seq_overview[9, "Total"] <-ISS.F %>% colSums() %>% sum()
  seq_overview[9, "Min"] <-ISS.F %>% colSums() %>% min()
  seq_overview[9, "Max"] <-ISS.F %>% colSums() %>% max()
  seq_overview[9, "Mean"] <-ISS.F %>% colSums() %>% mean() %>% round()
  
  seq_overview[10, "Total"] <-asv.F.all %>% nrow()
  seq_overview[10, "Min"] <-ifelse(asv.F.all >0,1,0) %>% colSums() %>% min()
  seq_overview[10, "Max"] <-ifelse(asv.F.all >0,1,0) %>% colSums() %>% max()
  seq_overview[10, "Mean"] <-ifelse(asv.F.all >0,1,0)  %>% colSums() %>% mean() %>% round()
  
  seq_overview[11, "Total"] <-asv.F %>% nrow()
  seq_overview[11, "Min"] <-ifelse(asv.F>0,1,0) %>% colSums() %>% min()
  seq_overview[11, "Max"] <-ifelse(asv.F>0,1,0) %>% colSums() %>% max()
  seq_overview[11, "Mean"] <-ifelse(asv.F >0,1,0)  %>% colSums() %>% mean() %>% round()
  
  seq_overview[12, "Total"] <-ISS.F[rowSums(ISS.F)!=0,] %>% nrow()
  seq_overview[12, "Min"] <-ifelse(ISS.F>0,1,0) %>% colSums() %>% min()
  seq_overview[12, "Max"] <-ifelse(ISS.F>0,1,0) %>% colSums() %>% max()
  seq_overview[12, "Mean"] <-ifelse(ISS.F >0,1,0)  %>% colSums() %>% mean() %>% round()
  
  seq_overview[,c(4:7)] <-round(seq_overview[,c(4:7)])
  
  #write.csv(seq_overview, file = "output/Number_of_seq_sample/Number_of_seq_overview.csv", row.names = F)
  
## Goods coverage
  
  # Percentage of samples that are singletons in a sample. Above 0.99 in generally good.
  Goods <- data.frame(sample =colnames(asv.F),
                      gc.F = rep(NA, ncol(asv.F)),
                      gc.B = rep(NA, ncol(asv.F)))
  
  for (i in 1:ncol(asv.F)) {
    Goods$gc.F[i] = 1 - (sum(asv.F[,i] == 1) / sum(asv.F[,] > 0))
    Goods$gc.B[i] = 1 - (sum(asv.B[,i] == 1) / sum(asv.B[,] > 0))
  }
  
  #write.csv(Goods, file = "output/Goods_coverage_F_B.csv")
# DNA concentration metabarcoding & qPCR (S)------
  # See plot calc file, correlates with 
  
# qPCR results 16S & ITS2 (Table S6)-----
  # See calc file

# Correlation analysis sample quality (S)-----
  
  # Make a large meta_file of three design files
  
  d1 <-design.alpha %>% select(sobs.F, shannon.F, evenness.F, sobs.B, shannon.B, evenness.B)
  d2 <-design.nseq %>% select(num_seq.F, num_seq.B, DNA_conc)
  colnames(d2) <- c("num_seq.F", "num_seq.B", "dna_23")
  d3 <-data_qpcr_meta %>% select(compost, rep, kit, size, dna_25, logcopyITS2, copies_per_ngDNA_ITS2, logcopy16s, copies_per_ngDNA_16S, ratioF_B)

  f1 <-merge(d3, d2, by =0)
  rownames(f1) <- f1$Row.names; f1$Row.names <- NULL
  design.all <- merge(f1, d1, by =0)
  rownames(design.all) <- design.all$Row.names; design.all$Row.names <- NULL
  rm(f1, d1, d2, d3)
  
  testRes <- cor.mtest(design.all[, c(5:19)], conf.level = 0.95, method ="spearman")
  
  #pdf(file = "figures/Correlation_plot_diversity_quanity.pdf", height = 10, width = 10)
  corrplot(cor(design.all[, c(5:19)], method ="pearson"), p.mat =testRes$p,
           sig.level = 0.05,order="hcl", insig ="blank", addCoef.col ="black",
           number.cex = 0.8, diag =FALSE, method = "color")
  
  #dev.off()
  testRes$p
  
  # check specific correlations
  ggplot(data =design.all, aes(sobs.B, copies_per_ngDNA_16S, color=compost))+ geom_point()
  ggplot(data =design.all, aes(sobs.F, copies_per_ngDNA_ITS2, color=compost))+ geom_point()
  
  ggplot(data =design.all %>% filter(compost =="C"), aes(sobs.B, copies_per_ngDNA_16S, color=compost))+ geom_point()
  
  ggplot(data =design.all %>% filter(compost =="C"), aes(sobs.F, copies_per_ngDNA_ITS2, color=compost))+ geom_point()
  
# Main analysis-------
  ## Alpha diversity-----
    ### Mean values & perc change-----
  # Alpha diversity
   mean.alpha.comp <-design.alpha %>% group_by(compost) %>% dplyr::summarize(
    mean.sobs.B = mean(sobs.B),
    mean.even.B = mean(evenness.B),
    mean.sobs.F = mean(sobs.F),
    mean.even.F = mean(evenness.F))
  
  (1/mean.alpha.comp[1,"mean.sobs.B"])*mean.alpha.comp[3,"mean.sobs.B"] # +12.3%
  (1/mean.alpha.comp[1,"mean.sobs.F"])*mean.alpha.comp[3,"mean.sobs.F"] # + 12.1%
  (1/mean.alpha.comp[3,"mean.even.B"])*mean.alpha.comp[1,"mean.even.B"] # -2.4%
  (1/mean.alpha.comp[1,"mean.even.F"])*mean.alpha.comp[3,"mean.even.F"] # 6
  
  mean.alpha.kit <-design.alpha %>% group_by(kit) %>% dplyr::summarize(
    mean.sobs.B = mean(sobs.B),
    mean.even.B = mean(evenness.B),
    mean.sobs.F = mean(sobs.F),
    mean.even.F = mean(evenness.F))
  
  
  (1/mean.alpha.kit[2,"mean.sobs.B"])*mean.alpha.kit[1,"mean.sobs.B"]
  (1/mean.alpha.kit[2,"mean.even.B"])*mean.alpha.kit[1,"mean.even.B"]
  
  mean.alpha.size <-design.alpha %>% group_by(size) %>% dplyr::summarize(
    mean.sobs.B = mean(sobs.B),
    mean.even.B = mean(evenness.B),
    mean.sobs.F = mean(sobs.F),
    mean.even.F = mean(evenness.F))
  
  (1/mean.alpha.size[1,"mean.sobs.B"])*mean.alpha.size[2,"mean.sobs.B"] # 4.6%
  (1/mean.alpha.size[1,"mean.even.B"])*mean.alpha.size[2,"mean.even.B"] # 3.4%
  
  # Quantity
  
  mean.quant.comp <- data_qpcr_meta %>% group_by(compost) %>% dplyr::summarize(
    mean.quant.B = mean(copies_per_ngDNA_16S),
    mean.quant.F = mean(copies_per_ngDNA_ITS2)
  )
  
  (1/mean.quant.comp[1, "mean.quant.B"])*mean.quant.comp[3, "mean.quant.B"] # + 66.7%
  (1/mean.quant.comp[3, "mean.quant.F"])*mean.quant.comp[1, "mean.quant.F"] # -4.2%
  
  mean.quant.kit <- data_qpcr_meta %>% group_by(kit) %>% dplyr::summarize(
    mean.quant.B = mean(copies_per_ngDNA_16S),
    mean.quant.F = mean(copies_per_ngDNA_ITS2)
  )
  
  (1/mean.quant.kit[2, "mean.quant.B"])*mean.quant.kit[1, "mean.quant.B"] # + 33.0%
  (1/mean.quant.kit[2, "mean.quant.F"])*mean.quant.kit[1, "mean.quant.F"] # + 14.1%
  
  mean.quant.size <- data_qpcr_meta %>% group_by(size) %>% dplyr::summarize(
    mean.quant.B = mean(copies_per_ngDNA_16S),
    mean.quant.F = mean(copies_per_ngDNA_ITS2)
  )
  
  (1/mean.quant.size[2, "mean.quant.B"])*mean.quant.size[1, "mean.quant.B"] # 35.1%
  (1/mean.quant.size[2, "mean.quant.F"])*mean.quant.size[1, "mean.quant.F"] # 43.8%
  
  
    ### FUNGI (Table 1, Table S6/S7)-----
  #sink(file="output/ANOVA_alpha_F.txt")
  ANOVA_sob =aov(sobs.F  ~ compost*size*kit, data = design.alpha)
  paste("Observed richness")
  summary(ANOVA_sob)
  ANOVA_shan =aov(shannon.F  ~ compost*size*kit, data = design.alpha)
  paste("Shannon diversity")
  summary(ANOVA_shan)
  ANOVA_ev <-aov(sqrt(evenness.F) ~ compost*size*kit , data = design.alpha)
  paste("Evenness")
  summary(ANOVA_ev)
  #sink()
  
  #pdf(file="figures/alpha_diversity/ANOVA_alpha_residual_plots_F.pdf", height =6, width = 6)
  par(mfrow=c(2,2))
  plot(ANOVA_sob)
  plot(ANOVA_shan)
  plot(ANOVA_ev)
  #dev.off()
  
  ylab <- c("Observed richness", "Shannon diversity", "Shannon evenness")
  factor <- c("sobs.F", "shannon.F", "evenness.F")
  for (i in 1:3) {
    design.alpha$y <- design.alpha[, factor[i]]
    results <-ANOVA_summary(data = design.alpha, x1 = "compost", x2 = "size", x3 = "kit", y= factor[i])
    
    p <-ggplot(data = design.alpha, aes(x = size, y = y, colour = compost, fill=compost)) +
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
      ylab(ylab[i]) +
      scale_color_manual(values = color.compost) +
      scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
      bg_theme +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            legend.position = "bottom", axis.title.x = element_blank())
    print(p)
    #ggsave(filename = paste0("figures/alpha_diversity/Comp_alpha_", factor[i], ".pdf"), height = 5, width = 10)
  }
  
  
    ### BACTERIA (Table 1, Table S6/S7)-----
  #sink(file="output/ANOVA_alpha_B.txt")
  ANOVA =aov(sobs.B  ~ compost*size*kit, data = design.alpha)
  paste("Observed richness")
  summary(ANOVA)
  ANOVA =aov(shannon.B  ~ compost*size*kit, data = design.alpha)
  paste("shannon diversity")
  summary(ANOVA)
  ANOVA =aov(invsimpson.B  ~ compost*size*kit, data = design.alpha)
  paste("Inversed Simpson")
  summary(ANOVA)
  ANOVA =aov(evenness.B  ~ compost*size*kit, data = design.alpha)
  paste("Evenness")
  summary(ANOVA)
  #sink()
  
  #pdf(file="figures/alpha_diversity/ANOVA_alpha_residual_plots_B.pdf", height =6, width = 6)
  par(mfrow=c(2,2))
  plot(ANOVA_sob)
  plot(ANOVA_shan)
  plot(ANOVA_ev)
  #dev.off()
  
  ylab <- c("Observed richnes", "Shannon diversity", "Shannon evenness")
  factor <- c("sobs.B", "shannon.B", "evenness.B")
  for (i in 1:3) {
    design.alpha$y <- design.alpha[, factor[i]]
    results <-ANOVA_summary(data = design.alpha, x1 = "compost", x2 = "size", x3 = "kit", y= factor[i])
    
    p <-ggplot(data = design.alpha, aes(x = size, y = y, colour = compost, fill=compost)) +
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
      ylab(ylab[i]) +
      scale_color_manual(values = color.compost) +
      scale_x_discrete(labels = c("0-2mm", "2-10mm")) +
      bg_theme +
      theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
            legend.position = "bottom", axis.title.x = element_blank())
    print(p)
    #ggsave(filename = paste0("figures/alpha_diversity/Comp_alpha_", factor[i], ".pdf"), height = 5, width = 10)
  }
  
  
    ### (FIGURE 2) Combined plot -----
  
  design.all$replicate_ID <-substr(rownames(design.all), 1,5) 
  design.all$kit <- as.factor(design.all$kit)
  design.all$size <- as.factor(design.all$size)
  design.all$compost <- as.factor(design.all$compost)
  
  #write.csv(design.all, file ="data/design_all_seq_alpha_quanty.csv", row.names = TRUE)
  
  factor.list <-c("copies_per_ngDNA_16S", "copies_per_ngDNA_ITS2", "sobs.B", "sobs.F",
    "evenness.B", "evenness.F")

  anova_list <- list()
  results <- list()
  
  # Run loop with all factors
  
  for (i in 1:length(factor.list)) {
    
  # Fit model
  model <- fit_and_report(data = design.all, factor.list[i], "compost * kit * size", "replicate_ID")

  # ANOVA table
  anova_table <- if (inherits(model, "lmerMod")) {
    Anova(model, type = 3)  # For mixed models
  } else {
    anova(model)            # For lm
  }
  anova_list[[i]] <-print(anova_table)
  
  # Test model assumptions
  
  #pdf(file = paste0("figures/alpha_diversity/LMM_residual_plots_", factor.list[i], ".pdf"))
  if (inherits(model, "lmerMod")) {
    par(mfrow =c(2,2))
    
    resid <- residuals(model)
    hist(resid)
    qqnorm(resid); qqline(resid)
    
    fitted_vals <- fitted(model)
    plot(fitted_vals, resid, 
         xlab = "Fitted values", ylab = "Residuals")
    abline(h = 0, col = "red")
    
    re <- ranef(model)$replicate_ID[[1]]  # or appropriate random effect
    qqnorm(re, main = "QQ plot of Random Effects")
    qqline(re)
  } else{
    par(mfrow =c(2,2))
    plot(model)
  }
  #dev.off()
  par(mfrow = c(1,1))
  
  emm <- emmeans(model, ~ compost * size * kit)
  cld_letters <- cld(emm, alpha = 0.05, Letters = letters, adjust = "tukey")
  cld_df <- cld_letters %>%
    as.data.frame() %>%
    select(compost, size, kit, .group) %>%
    dplyr::rename(letter = .group)
  
  
  summary_data <- design.all %>%
    group_by(compost, size, kit) %>%
    summarise(
      mean = mean(.data[[factor.list[i]]], na.rm = TRUE),
      max = max(.data[[factor.list[i]]], na.rm = TRUE),
      sd = sd(.data[[factor.list[i]]], na.rm = TRUE),
      se = sd / sqrt(n()),
      .groups = 'drop'
    ) %>%
    left_join(cld_df, by = c("compost", "size", "kit"))
  
  summary_data$pos.letter <- summary_data$max+ 0.01*summary_data$max  # Adjust as needed
  summary_data$size <- as.factor(summary_data$size)
  
  summary_data$letter <- gsub(" ", "", summary_data$letter)
  
  results[[i]] <- summary_data
  }

  
  # Extract ANOVA type 3 statistics
  
  anova_df <- list()
  for (i in 1:6) {
    data <- anova_list[[i]] %>% as.data.frame()
    data$factor <- rownames(data); rownames(data) <-NULL
    data$response <- factor.list[i]
    anova_df[[i]] <- data
  }
  
  data_lm <-rbind(anova_df[[1]], anova_df[[5]], anova_df[[6]]) %>% select("F value", "Pr(>F)", "factor", "response")
  colnames(data_lm) <- c("F", "p", "factor", "response")
  
  data_lmer <-rbind(anova_df[[2]], anova_df[[3]], anova_df[[4]]) %>% select("Chisq", "Pr(>Chisq)", "factor", "response")
  colnames(data_lmer) <- c("F", "p", "factor", "response")
  
  # Save output
  
  for (i in 1:length(factor.list)) {
    write.table(anova_df[[i]], file = paste0("output/alpha_diversity/ANOVA_alpha_B_mixed_", factor.list[i], ".csv"),  sep = ";", row.names = FALSE)
  }
  
  
  # Beta diversity
  beta_B <-adonis2(ISS.bray.B ~ compost*size*kit,
          data = design.all,
          permutations = 999,
          method = "bray",
          strata = design.all$replicate_uid) %>% as.data.frame()
  
  beta_F <-adonis2(ISS.bray.F ~ compost*size*kit,
          data = design.all,
          permutations = 999,
          method = "bray",
          strata = design.all$replicate_uid) %>% as.data.frame()
  
  
  beta_B$factor <-rownames(beta_B); rownames(beta_B) <-NULL
  beta_B$response <- "bstructure"
  
  beta_F$factor <-rownames(beta_F); rownames(beta_F) <-NULL
  beta_F$response <- "Fstructure"
  
  data_perma <-rbind(beta_B, beta_F) %>% select("F", "Pr(>F)", "factor", "response")
  colnames(data_perma) <-  c("F", "p", "factor", "response")
  
  data_F <-rbind(data_lm, data_lmer, data_perma) %>% filter(!factor %in% c("Residuals", "(Intercept)", "Residual", "Total")) %>% select("F", "factor", "response")
  #data_F <-data_F %>% pivot_wider(names_from = factor, values_from = F)
  
  data_p <-rbind(data_lm, data_lmer, data_perma) %>% filter(!factor %in% c("Residuals", "(Intercept)", "Residual", "Total")) %>% select("p", "factor", "response")
  #data_p <-data_p %>% pivot_wider(names_from = factor, values_from = p)
  
  data_plot <- merge(data_F, data_p, by = c("factor", "response"))
  
  data_sig <- subset(data_plot, p < 0.05)
  
  
  response_labels <- c("Bacterial quantity",  "Bacterial richness", "Bacterial evenness", "Bacterial structure",
                       "Fungal quantity", "Fungal richness",  "Fungal evenness","Fungal structure")
  
  data_sig$response<- factor(data_sig$response, 
                              levels = c("copies_per_ngDNA_16S", "sobs.B", "evenness.B","bstructure",
                                         "copies_per_ngDNA_ITS2","sobs.F", "evenness.F",  "Fstructure"),
                             labels = response_labels)
  
  data_sig$factor<- factor(data_sig$factor, 
                             levels = c("compost", "size", "kit", "compost:size",
                                        "compost:kit", "kit:size", "compost:kit:size"))
  

  #ggplot(data_sig, aes(x = response, y = factor, fill = p)) +
  #  geom_tile(color = "white") +
  #  geom_text(aes(label = sprintf("%.2f", F)), color = "black", size = 4) +
  #  scale_fill_gradient(low = "darkgreen", high = "lightgreen", name = "p-value") +
  #  theme_minimal() +
  #  theme(axis.text.x = element_text(angle = 45, hjust = 1),
   #       panel.grid.major = element_blank())+
   # scale_y_discrete(limits = rev(levels(data_sig$factor))) +xlab("") + ylab("")
  
  #ggsave(filename ="figures/Figure_4_Statistics.pdf", height =4, width=6 )
  
  data <- design.all
  
  ylabel <- c("Copies 16S/ ng DNA", "Copies ITS2/ ng DNA", "Bacterial richness", "Fungal richness", "Bacterial evenness",
              "Fungal evenness")
  
  
  plot_list <-list()
  ylower <- c(0,0,2000,100, 0.6, 0)
  yupper <- c(700000, 6000,3000, 300,0.75, 0.6)
  
  mixed_label <- function(threshold = 10000, digits = 3) {
    function(x) {
      sapply(x, function(val) {
        if (is.na(val)) {
          return(NA)
        } else if (abs(val) < threshold) {
          format(val, scientific = FALSE, big.mark = ",")
        } else {
          format(val, scientific = TRUE, digits = digits)
        }
      })
    }
  }

  for (i in 1:length(factor.list)) {
  
  data$y <-design.all[, factor.list[i]]
  plot_list[[i]] <-ggplot(data = data, aes(x = size, y = y, shape = kit, color = size)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_line(aes(group = interaction(size, kit)), position = position_dodge(width = 0.5)) +
    facet_wrap(~compost) +
    geom_text(
      data = results[[i]], 
      aes(x = size, y = pos.letter, label = letter, fill =kit),  
      position = position_dodge(width = 0.5),  # Apply consistent dodging
      vjust = -0.5,  # Place text slightly above points,
      color="darkred",
      size = 5, 
      inherit.aes = FALSE
    ) +
    ylab(ylabel[i])+
    scale_color_manual(values = c("blue", "orange"), name ="Particle fraction",
                       labels = c("fine", "coarse") ) +
    scale_y_continuous(labels = mixed_label(), limits = c(ylower[i], yupper[i]))+
    scale_shape_manual(values=c(15,17), labels = c("Bead-based", "Column-based"), name ="DNA extraction")+
    scale_x_discrete(labels = c("", "")) +
    bg_theme +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
          legend.position = "bottom", axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  }
  
  
  ggarrange(plotlist = plot_list, common.legend = TRUE, nrow=3, ncol =2, labels = c("A", "B", "C", "D", "E", "F"))
  
   
   
   combined_plot <- wrap_plots(plot_list, nrow = 3, ncol = 2, byrow = TRUE) +
     plot_annotation(tag_levels = 'A') +
     plot_layout(guides = "collect") & 
     theme(
       axis.title.y = element_text(angle = 90),
       axis.text.y = element_text(hjust = 1),
       plot.margin = margin(5, 5, 5, 5),
       legend.position = "bottom"
     )
   
   print(combined_plot)
   ggsave(filename = "figures/Figure_2_quantity_diversity.pdf", height =14, width = 11)
  
  ## Beta diversity, (FIGURE 3, Table 1, Table S8)-----
  
  ### FUNGI
  ISS.bray.F <- vegdist(ISS.F %>% t(), method="bray")
  asv.bray.ab.F <-vegdist(asv.F.absolut %>% t(), method ="bray") # Absolute numbers
  
  # Check the order!!
  (rownames(as.matrix(ISS.bray.F))==rownames(design.alpha)) %>% sum()
  (rownames(as.matrix(asv.bray.ab.F))==rownames(design.alpha)) %>% sum()
  
  #sink(file="output/Beta_diversity/PERMANOVA_comp_size_kit_F.txt")
  adonis2(ISS.bray.F~compost*size*kit, design.alpha)
  adonis2(asv.bray.ab.F~compost*size*kit, design.alpha)
  #sink()
  
  ## BACTERIA
  
  ISS.bray.B <- vegdist(ISS.B %>% t(), method="bray")
  asv.bray.ab.B <- vegdist(asv.B.absolut %>% t(), method="bray")
  
  # Check the order!!
  (rownames(as.matrix(ISS.bray.B))==rownames(design.alpha)) %>% sum()
  (rownames(as.matrix(asv.bray.ab.B))==rownames(design.alpha)) %>% sum()
  
  #sink(file="output/Beta_diversity/PERMANOVA_comp_size_kit_B.txt")
  adonis2(ISS.bray.B~compost*size*kit, design.alpha)
  adonis2(asv.bray.ab.B~compost*size*kit, design.alpha)
  #sink()  
  
  # Calculate PCoA
  
  p1 <-PcoA_beta(ISS.bray.B, output_file ="figures/beta_diversity/ISS.bray.B.PcoA.pdf", save = FALSE)+ labs(tag ="A")
  PcoA_beta(asv.bray.ab.B, output_file ="figures/beta_diversity/asv.bray.ab.B.PcoA.pdf", save = TRUE)
  
  p2 <- PcoA_beta(ISS.bray.F, output_file ="figures/beta_diversity/ISS.bray.F.PcoA.pdf", save = FALSE)+ labs(tag ="B")
  PcoA_beta(asv.bray.ab.F, output_file ="figures/beta_diversity/asv.bray.ab.F.PcoA.pdf", save = TRUE) 

  
  plot_combined <- p1 / p2 + plot_layout(guides = "collect") & theme(legend.position = "right")
 
  
  #, FGggsave(plot_combined, file ="figures/Figure_2_ISS_bray_B_F_combined.pdf", height =8, width = 8)
  

  ## Taxonomic composition (Figure S4)-----
    ### FUNGI-----
  #ISS.prop <- prop.table(as.matrix(ISS.F), margin = 2) * 100  # get proportions 
  ISS.prop <- asv.F.absolut %>% as.matrix()  # get proportions 
  
  design.order <- design.alpha[with(design.alpha, order(compost,size, kit)),]
  # design.order <- design.alpha[with(design.alpha, order(compost, size, rep, kit)),]
  
  ## Phylum
  
  ISS.p <- aggregate(ISS.prop ~ Phyla, data = tax.F[rownames(asv.F.absolut),], FUN = sum)  # aggregate by phylum
  ISS.p <- as.matrix(data.frame(ISS.p[, -1], row.names = ISS.p[, 1], check.names = FALSE))  # move first col to row.names
  ISS.p <- ISS.p[order(rowSums(ISS.p), decreasing = TRUE), ]  # order by rowSums
  ISS.p.nc <- as.data.frame(ISS.p[rownames(ISS.p) != "unclassified", ]) 
  ISS.p.nc <- as.data.frame(ISS.p[!grepl("unclassified", rownames(ISS.p)), ])
  ISS.p.nc["unclassified", ] <- ISS.p[rownames(ISS.p) == "unclassified", drop = F]  # add unclassified at the end
  ISS.p.plot <- ISS.p.nc
  ISS.p.plot <- as.matrix(ISS.p.plot)
  ISS.p.plot <- ISS.p.plot[, rownames(design.order)]
  
  # plot
  #pdf(file="figures/beta_diversity/taxonomy_phylum_F_absolute.pdf", height=5, width=8)
  par(omi = c(0, 0, 0, 0), mai = c(1.5, 0.9, 0.3, 1.8))
  color <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F",
             "#FF7F00", "#CAB2D6", "#000000")
  barplot(ISS.p.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
          width = 1, cex.axis = 0.8, cex.names = 0.5, axisnames = TRUE)
  par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1.4, 0.9, 0.3, 0.1), new = TRUE)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("topright", legend = rownames(ISS.p.plot), pch = 22, cex = 0.8, bty = "n",
         ncol = 1, pt.bg = color, col = NA)
  #dev.off()
  rm(color, ISS.p, ISS.p.nc, ISS.p.plot) 

  # Aggregate by family level
  ISS.f <- aggregate(ISS.prop ~ Family, data = tax.F[rownames(asv.F.absolut),], FUN = sum)  # aggregate by family
  ISS.f <- as.matrix(data.frame(ISS.f[, -1], row.names = ISS.f[, 1], check.names = FALSE))  # move first col to row.names
  ISS.f <- ISS.f[order(rowSums(ISS.f), decreasing = TRUE), ]  # order by rowSums
  ISS.f.nc <- as.data.frame(ISS.f[rownames(ISS.f) != "unclassified", ])  # remove unclassified group
  ISS.f.nc["others", ] <- colSums(ISS.f.nc[min(nrow(ISS.f.nc), 21):nrow(ISS.f.nc), ])  # merge rare groups to others
  ISS.f.nc["unclassified", ] <- ISS.f[rownames(ISS.f) == "unclassified", drop = F]  # add unclassified at the end
  ISS.f.plot <- ISS.f.nc[-(min(nrow(ISS.f.nc), 21):(nrow(ISS.f.nc) - 2)), ]  # remove rare groups
  ISS.f.plot <- as.matrix(ISS.f.plot)
  ISS.f.plot = ISS.f.plot[, rownames(design.order)]
  
  # plot
  pdf(file="figures/beta_diversity/taxonomy_family_F_absolute.pdf", height=5, width=12)
  par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 2.5))
  color <-  c(replicate(20, generate_random_color()), "#999999", "#000000")
  barplot(ISS.f.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
          width = 1, cex.axis = 0.8, cex.names = 0.5, axisnames = TRUE)
  par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1.4, 0.9, 0.3, 0.1), new = TRUE)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("topright", legend = rownames(ISS.f.plot), pch = 22, cex = 0.8, bty = "n",
         ncol = 1, pt.bg = color, col = NA)
  dev.off()
  rm(color, ISS.f, ISS.f.nc, ISS.f.plot) 
  
  # Aggregate by genus level
  ISS.g <- aggregate(ISS.prop ~ Genus, data = tax.F[rownames(asv.F.absolut),], FUN = sum)  # aggregate by family
  ISS.g <- as.matrix(data.frame(ISS.g[, -1], row.names = ISS.g[, 1], check.names = FALSE))  # move first col to row.names
  ISS.g <- ISS.g[order(rowSums(ISS.g), decreasing = TRUE), ]  # order by rowSums
  ISS.g.nc <- as.data.frame(ISS.g[rownames(ISS.g) != "unclassified", ])  # remove unclassified group
  ISS.g.nc["others", ] <- colSums(ISS.g.nc[min(nrow(ISS.g.nc), 21):nrow(ISS.g.nc), ])  # merge rare groups to others
  ISS.g.nc["unclassified", ] <- ISS.g[rownames(ISS.g) == "unclassified", drop = F]  # add unclassified at the end
  ISS.g.plot <- ISS.g.nc[-(min(nrow(ISS.g.nc), 21):(nrow(ISS.g.nc) - 2)), ]  # remove rare groups
  ISS.g.plot <- as.matrix(ISS.g.plot)
  ISS.g.plot = ISS.g.plot[, rownames(design.order)]
  
  # plot
  pdf(file="figures/beta_diversity/taxonomy_genus_F_absolute.pdf", height=6, width=11)
  par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 1.5))
  color <-  c(color.genus, "#999999", "#000000")
  barplot(ISS.g.plot, col = color, ylab = "Absolute abundances", border = NA, las = 2,
          width = 1, cex.axis = 0.9, cex.names = 0.8, axisnames = TRUE)
  par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1.4, 0.9, 0.3, 0.1), new = TRUE)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("bottomright", legend = rev(rownames(ISS.g.plot)), pch = 22, cex = 0.9, bty = "n",
         ncol = 1, pt.bg = rev(color), col = NA)
  dev.off()
  rm(color, ISS.g, ISS.g.nc, ISS.g.plot) 
  
    ### BACTERIA----
  
  ISS.prop <- asv.B.absolut %>% as.matrix()
  design.order <- design.alpha[with(design.alpha, order(compost,size, kit)),]
  ## Phylum
  
  ISS.p <- aggregate(ISS.prop ~ Phyla, data = tax.B[rownames(asv.B.absolut),], FUN = sum)  # aggregate by phylum
  ISS.p <- as.matrix(data.frame(ISS.p[, -1], row.names = ISS.p[, 1], check.names = FALSE))  # move first col to row.names
  ISS.p <- ISS.p[order(rowSums(ISS.p), decreasing = TRUE), ]  # order by rowSums
  ISS.p.nc <- as.data.frame(ISS.p[rownames(ISS.p) != "unclassified", ]) 
  ISS.p.nc["others", ] <- colSums(ISS.p.nc[min(nrow(ISS.p.nc), 21):nrow(ISS.p.nc), ])  # merge rare groups to others
  ISS.p.nc <- as.data.frame(ISS.p[!grepl("unclassified", rownames(ISS.p)), ])
  ISS.p.nc["unclassified", ] <- ISS.p[rownames(ISS.p) == "unclassified", drop = F]  # add unclassified at the end
  ISS.p.plot <- ISS.p.nc[-(min(nrow(ISS.p.nc), 21):(nrow(ISS.p.nc) - 2)), ]  # remove rare groups
  ISS.p.plot <- as.matrix(ISS.p.plot)
  ISS.p.plot = ISS.p.plot[, rownames(design.order)]
  
  # plot
  pdf(file="figures/beta_diversity/taxonomy_phylum_B_absolute.pdf", height=5, width=8)
  par(omi = c(0, 0, 0, 0), mai = c(1.5, 0.9, 0.3, 2))
  color <-  c(replicate(20, generate_random_color()), "#999999", "#000000")
  barplot(ISS.p.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
          width = 1, cex.axis = 0.8, cex.names = 0.5, axisnames = TRUE)
  par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(1.4, 0.9, 0.3, 0.1), new = TRUE)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("topright", legend = rownames(ISS.p.plot), pch = 22, cex = 0.8, bty = "n",
         ncol = 1, pt.bg = color, col = NA)
  dev.off()
  rm(color, ISS.p, ISS.p.nc, ISS.p.plot) 
  
  
  # Aggregate by family level
  ISS.f <- aggregate(ISS.prop ~ Family, data = tax.B[rownames(asv.B.absolut),], FUN = sum)  # aggregate by family
  ISS.f <- as.matrix(data.frame(ISS.f[, -1], row.names = ISS.f[, 1], check.names = FALSE))  # move first col to row.names
  ISS.f <- ISS.f[order(rowSums(ISS.f), decreasing = TRUE), ]  # order by rowSums
  ISS.f.nc <- as.data.frame(ISS.f[rownames(ISS.f) != "unclassified", ])  # remove unclassified group
  ISS.f.nc["others", ] <- colSums(ISS.f.nc[min(nrow(ISS.f.nc), 21):nrow(ISS.f.nc), ])  # merge rare groups to others
  ISS.f.nc["unclassified", ] <- ISS.f[rownames(ISS.f) == "unclassified", drop = F]  # add unclassified at the end
  ISS.f.plot <- ISS.f.nc[-(min(nrow(ISS.f.nc), 21):(nrow(ISS.f.nc) - 2)), ]  # remove rare groups
  ISS.f.plot <- as.matrix(ISS.f.plot)
  ISS.f.plot = ISS.f.plot[, rownames(design.order)]
  
  # plot
  pdf(file="figures/beta_diversity/taxonomy_family_B_absolute.pdf", height=5, width=13)
  par(omi = c(0, 0, 0, 0), mai = c(1, 0.9, 0.3, 3.2))
  color <-  c(replicate(20, generate_random_color()), "#999999", "#000000")
  barplot(ISS.f.plot, col = color, ylab = "Relative abundance [%]", border = NA, las = 2,
          width = 1, cex.axis = 0.8, cex.names = 0.5, axisnames = TRUE)
  par(fig = c(0, 1, 0, 1), omi = c(0, 0, 0, 0), mai = c(0.5, 0.9, 0.3, 0.1), new = TRUE)
  plot(0, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  legend("topright", legend = rownames(ISS.f.plot), pch = 22, cex = 0.8, bty = "n",
         ncol = 1, pt.bg = color, col = NA)
  dev.off()
  rm(color, ISS.f, ISS.f.nc, ISS.f.plot) 
  
  # Aggregate by genus level
  ISS.g <- aggregate(ISS.prop ~ Genus, data = tax.B[rownames(asv.B.absolut),], FUN = sum)  # aggregate by family
  ISS.g <- as.matrix(data.frame(ISS.g[, -1], row.names = ISS.g[, 1], check.names = FALSE))  # move first col to row.names
  ISS.g <- ISS.g[order(rowSums(ISS.g), decreasing = TRUE), ]  # order by rowSums
  ISS.g.nc <- as.data.frame(ISS.g[rownames(ISS.g) != "unclassified", ])  # remove unclassified group
  ISS.g.nc["others", ] <- colSums(ISS.g.nc[min(nrow(ISS.g.nc), 21):nrow(ISS.g.nc), ])  # merge rare groups to others
  ISS.g.nc["unclassified", ] <- ISS.g[rownames(ISS.g) == "unclassified", drop = F]  # add unclassified at the end
  ISS.g.plot <- ISS.g.nc[-(min(nrow(ISS.g.nc), 21):(nrow(ISS.g.nc) - 2)), ]  # remove rare groups
  ISS.g.plot <- as.matrix(ISS.g.plot)
  ISS.g.plot = ISS.g.plot[, rownames(design.order)]
  
  
inserta 
  

  rm(color, ISS.g, ISS.g.nc, ISS.g.plot) 
  
# Differential abundance analysis-------
  ## BACTERIA------
   ### ASV level (FIGURE 4, Table S5)------
  # Calculate the differences
  results <- list()
  
  #sink(file ="output/Indicator_analysis_ab_ab/Bacteria_overview_number_of_ASV.txt")
  for (main_factor in c("size", "kit")) {
    
  for (compost in c("K20", "K21", "K22")) {
    
    if(main_factor =="size"){
      side_factor_levels <- c("M","S")
    }
      
     if(main_factor =="kit"){
       side_factor_levels <- c("10", "2")
     }
      for (sidefactor in side_factor_levels) {
        results[[paste0(main_factor,compost, sidefactor)]] <- consistency_levels_paired(asv.B.absolut, tax.B, "ASV", main_factor, compost, "B", min_rel = 0.0001, sidefactor)
          }
     }
  }
  #sink()
  

  d1 <-differ_analysis_plot(res1 ="sizeK20M", res2 ="sizeK20S")
  d2 <-differ_analysis_plot(res1 ="sizeK21M", res2 ="sizeK21S")
  d3 <-differ_analysis_plot(res1 ="sizeK22M", res2 ="sizeK22S")
  d4 <-differ_analysis_plot(res1 ="kitK2010", res2 ="kitK202")
  d5 <-differ_analysis_plot(res1 ="kitK2110", res2 ="kitK212")
  d6 <-differ_analysis_plot(res1 ="kitK2210", res2 ="kitK222")
  
  d1$ASV <- rownames(d1); rownames(d1) <-NULL; d1$group <- "sizeK20"
  d2$ASV <- rownames(d2); rownames(d2) <-NULL; d2$group <- "sizeK21"
  d3$ASV <- rownames(d3); rownames(d3) <-NULL; d3$group <- "sizeK22"
  d4$ASV <- rownames(d4); rownames(d4) <-NULL; d4$group <- "kitK20"
  d5$ASV <- rownames(d5); rownames(d5) <-NULL; d5$group <- "kitK21"
  d6$ASV <- rownames(d6); rownames(d6) <-NULL; d6$group <- "kitK22"
  
  d <-rbind(d1, d2, d3, d4, d5, d6)
  
  
  unique_phyla <-d$Phyla %>% unique()
  phylum_colors <- c(paletteMartin[2:15], lighten(paletteMartin[2:15], 0.3))
  phylum_colors <- setNames(phylum_colors, unique_phyla)
  phylum_shape <- c(rep(16,14), rep(18, 14))
  phylum_shape <- setNames(phylum_shape, unique_phyla)
  
  p1 <-ggplot(d %>% filter(group =="sizeK20"), aes(x = log2FCmean, y = log10(countmean + 1), color =Phyla, shape = Phyla)) +
    geom_point(size =2) + 
    labs(x = "Log2 Fold Change", y = "Log10 Mean Abundance") +
    scale_color_manual(values = phylum_colors) +
    scale_shape_manual(values= phylum_shape)+
    theme_classic()+ geom_vline(xintercept = 0)+ theme(legend.position = "none")+
    xlim(-6,6) + ylim(1,4)+
    annotate("text", x = c(-5,5), y = c(4,4),
             label = c("coarse", "fine"), size = 3.5)+
    ggtitle("Compost A") + theme(plot.title = element_text(hjust =0.5))
  
  p2 <- ggplot(d %>% filter(group =="sizeK21"), aes(x = log2FCmean, y = log10(countmean + 1), color =Phyla, shape = Phyla)) +
    geom_point(size =2) + 
    labs(x = "Log2 Fold Change", y = "Log10 Mean Abundance") +
    scale_color_manual(values = phylum_colors) +
    scale_shape_manual(values= phylum_shape)+
    theme_classic()+ geom_vline(xintercept = 0)+ theme(legend.position = "none")+
    xlim(-6,6)+ ylim(1,4)+
    annotate("text", x = c(-5,5), y = c(4,4),
             label = c("coarse", "fine"), size = 3.5)+
    ggtitle("Compost B") + theme(plot.title = element_text(hjust =0.5))
  
  p3 <- ggplot(d %>% filter(group =="sizeK22"), aes(x = log2FCmean, y = log10(countmean + 1), color =Phyla, shape = Phyla)) +
    geom_point(size =2) + 
    labs(x = "Log2 Fold Change", y = "Log10 Mean Abundance") +
    scale_color_manual(values = phylum_colors) +
    scale_shape_manual(values= phylum_shape)+
    theme_classic()+ geom_vline(xintercept = 0)+ theme(legend.position = "none")+
    xlim(-6,6)+ ylim(1,4)+
    annotate("text", x = c(-5,5), y = c(4,4),
             label = c("coarse", "fine"), size = 3.5)+
    ggtitle("Compost C") + theme(plot.title = element_text(hjust =0.5))

  
  p4 <-ggplot(d %>% filter(group =="kitK20"), aes(x = log2FCmean, y = log10(countmean + 1), color =Phyla, shape = Phyla)) +
    geom_point(size =2) + 
    labs(x = "Log2 Fold Change", y = "Log10 Mean Abundance") +
    scale_color_manual(values = phylum_colors) +
    scale_shape_manual(values= phylum_shape)+
    theme_classic()+ geom_vline(xintercept = 0)+ theme(legend.position = "none")+
    xlim(-6,6)+ ylim(1,4)+
    annotate("text", x = c(-4.8,5.2), y = c(4,4),
             label = c("beads", "column"), size = 3.5,
             color = c("black", "black"))+
    ggtitle("Compost A") + theme(plot.title = element_text(hjust =0.5))
  
  p5 <- ggplot(d %>% filter(group =="kitK21"), aes(x = log2FCmean, y = log10(countmean + 1), color =Phyla, shape = Phyla)) +
    geom_point(size =2) + 
    labs(x = "Log2 Fold Change", y = "Log10 Mean Abundance") +
    scale_color_manual(values = phylum_colors) +
    scale_shape_manual(values= phylum_shape)+
    theme_classic()+ geom_vline(xintercept = 0)+ theme(legend.position = "none")+
    xlim(-6,6)+ ylim(1,4)+
    annotate("text", x = c(-4.8,5.2), y = c(4,4),
             label = c("beads", "column"), size = 3.5,
             color = c("black", "black"))+
    ggtitle("Compost B") + theme(plot.title = element_text(hjust =0.5))
  
  p6 <- ggplot(d %>% filter(group =="kitK22"), aes(x = log2FCmean, y = log10(countmean + 1), color =Phyla, shape = Phyla)) +
    geom_point(size =2) + 
    labs(x = "Log2 Fold Change", y = "Log10 Mean Abundance") +
    scale_color_manual(values = phylum_colors) +
    scale_shape_manual(values= phylum_shape)+
    theme_classic()+ geom_vline(xintercept = 0)+ theme(legend.position = "none")+
    xlim(-6,6)+ ylim(1,4)+
    annotate("text", x = c(-4.8,5.2), y = c(4,4),
             label = c("beads", "column"), size = 3.5,
             color = c("black", "black"))+
    ggtitle("Compost C") + theme(plot.title = element_text(hjust =0.5))
  
  p_dummy <- ggplot(d, aes(x = log2FCmean, y = countmean, color = Phyla, shape = Phyla)) +
    geom_point() +
    scale_color_manual(values = phylum_colors) +
    scale_shape_manual(values= phylum_shape) +
    theme_classic()+ theme(legend.position = "bottom",
                           legend.title=element_blank(),
                          legend.key.height = unit(0.4, "lines"))+
    guides(
      color = guide_legend(override.aes = list(size = 2.5), nrow= 5),
      shape = guide_legend(override.aes = list(size = 2.5), nrow=5)
    )
  
  legend <- ggpubr::get_legend(p_dummy)
  
  plots <- ggarrange(p1, p2, p3, p4, p5, p6, ncol = 3, nrow = 2)
  final_plot <- plot_grid(plots, legend, ncol = 1, rel_heights = c(1, 0.2))
  final_plot
  #ggsave(filename = "figures/Figure3_ASV_diff_ab.pdf", height = 8, width = 10)
  
  #write.csv(d, file = "output/Indicator_analysis_ab_ab/All_differential_abundance_ASVs_bacteria.csv")
  
  # Calculate summed relative abundance
  
  d_select <-d %>% select(ASV, group, consistency_flag)
  d_select$compost <- str_sub(d_select$group, -3)
  levels <-d_select$group %>% unique()
  
  asv.B.prop <-prop.table(as.matrix(asv.B), margin = 2) * 100 
  
  
  
  for (i in levels) {
  print(i)
  d_select_d <-d_select %>% filter(group ==i & consistency_flag =="Consistently Increased")
  print("fine or column")
  if (nrow(d_select_d) >1) {
  asv.B.prop[d_select_d$ASV, grepl(d_select_d$compost[1], colnames(asv.B.prop))] %>% colSums() %>% mean() %>% round(2) %>% print()
  } else {
    
    asv.B.prop[d_select_d$ASV, grepl(d_select_d$compost[1], colnames(asv.B.prop))]%>% mean() %>% round(2) %>% print()
  }
    
  d_select_i <-d_select %>% filter(group ==i & consistency_flag =="Consistently Decreased") 
  print("coarse or beads")
  if (nrow(d_select_i) >1) {
    asv.B.prop[d_select_i$ASV, grepl(d_select_i$compost[1], colnames(asv.B.prop))] %>% colSums() %>% mean() %>% round(2) %>% print()
  } else {
    asv.B.prop[d_select_i$ASV, grepl(d_select_i$compost[1], colnames(asv.B.prop))] %>% mean() %>% round(2) %>% print()
  }
  }
  
  rm(d_select, levels, d_select_d, d_select_i)
  # ASVs enriched in all?
  
  ### SIZE
  d1 <-differ_analysis_plot(res1 ="sizeK20M", res2 ="sizeK20S", plot =FALSE)
  d2 <-differ_analysis_plot(res1 ="sizeK21M", res2 ="sizeK21S", plot =FALSE)
  d3 <-differ_analysis_plot(res1 ="sizeK22M", res2 ="sizeK22S", plot =FALSE)
  
  d1_in <-d1 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d2_in <-d2 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d3_in <-d3 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()

  d1_de <-d1 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d2_de <-d2 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d3_de <-d3 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  
  # Create a excle with all the differential abundance ASVs
  d_in <-data.frame(ASV = c(d1_in, d2_in, d3_in), factor = "2mm",
             compost = c(rep("A", length(d1_in)), rep("B", length(d2_in)), rep("C", length(d3_in))) )
  
  d_de <-data.frame(ASV = c(d1_de, d2_de, d3_de), factor = "10mm",
                    compost = c(rep("A", length(d1_de)), rep("B", length(d2_de)), rep("C", length(d3_de))) )
  
  d_size <-rbind(d_in, d_de)
  
  tax.B.adapt <- tax.B
  tax.B.adapt$ASV <- rownames(tax.B.adapt); rownames(tax.B.adapt) <- NULL
  
  d_size <-d_size %>% left_join(tax.B.adapt[, c("Phyla", "Family", "Genus", "ASV")], by ="ASV")
  
  summary_size <-d_size %>%
    group_by(factor, compost, Phyla) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Phyla, values_from = count, values_fill = 0) %>% as.data.frame()
  
  rownames(summary_size) <- paste0(summary_size$factor, "_", summary_size$compost)
  summary_size$factor <- NULL; summary_size$compost <- NULL
  
  #summary_size %>% t() %>% write.csv(file= "output/Indicator_analysis_ab_ab/size/ASV/Bacteria_size_summary_ASV_phyla_all.csv")
  
  #sink(file = "output/Indicator_analysis_ab_ab/size/B_shared_size_among_composts.txt")
  paste("Enriched in 2mm")
  paste("A and B");tax.B[intersect(d1_in, d2_in),c("Phyla", "Family", "lowest")]
  paste("A and C");tax.B[intersect(d1_in, d3_in),c("Phyla", "Family", "lowest")]
  paste("B and C");tax.B[intersect(d2_in, d3_in),c("Phyla", "Family", "lowest")]
  paste("All three");tax.B[intersect(d1_in, intersect(d2_in, d3_in)), c("Phyla", "Family", "lowest")]
  paste("Enriched in 10mm")
  paste("A and B");tax.B[intersect(d1_de, d2_de),c("Phyla", "Family", "lowest")]
  paste("A and C");tax.B[intersect(d1_de, d3_de),c("Phyla", "Family", "lowest")]
  paste("B and C");tax.B[intersect(d2_de, d3_de),c("Phyla", "Family", "lowest")]
  paste("All three");tax.B[intersect(d1_de, intersect(d2_de, d3_de)), c("Phyla", "Family", "lowest")]
  #sink()
  

  # KIT
  d1 <-differ_analysis_plot(res1 ="kitK2010", res2 ="kitK202", plot =FALSE)
  d2 <-differ_analysis_plot(res1 ="kitK2110", res2 ="kitK212", plot =FALSE)
  d3 <-differ_analysis_plot(res1 ="kitK2210", res2 ="kitK222", plot =FALSE)
  
  d1_in <-d1 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d2_in <-d2 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d3_in <-d3 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  
  d1_de <-d1 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d2_de <-d2 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d3_de <-d3 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  
  # Create a excle with all the differential abundance ASVs
  d_in <-data.frame(ASV = c(d1_in, d2_in, d3_in), factor = "NucleoSpin",
                    compost = c(rep("A", length(d1_in)), rep("B", length(d2_in)), rep("C", length(d3_in))) )
  
  d_de <-data.frame(ASV = c(d1_de, d2_de, d3_de), factor = "NucleoMag",
                    compost = c(rep("A", length(d1_de)), rep("B", length(d2_de)), rep("C", length(d3_de))) )
  
  d_kit<-rbind(d_in, d_de)
  
  d_kit <-d_kit %>% left_join(tax.B.adapt[, c("Phyla", "Family", "Genus", "ASV")], by ="ASV")
  
  summary_kit <-d_kit %>%
    group_by(factor, compost, Phyla) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Phyla, values_from = count, values_fill = 0) %>% as.data.frame()
  
  rownames(summary_kit) <- paste0(summary_kit$factor, "_", summary_kit$compost)
  summary_kit$factor <- NULL; summary_kit$compost <- NULL
  
  #summary_kit %>% t() %>% write.csv(file= "output/Indicator_analysis_ab_ab/kit/ASV/Bacteria_kit_summary_ASV_phyla_all.csv")
  
  #sink(file = "output/Indicator_analysis_ab_ab/kit/B_shared_kit_among_composts.txt")
  paste("Enriched in NucleoSpin")
  paste("A and B");tax.B[intersect(d1_in, d2_in),c("Phyla", "Family", "lowest")]
  paste("A and C");tax.B[intersect(d1_in, d3_in),c("Phyla", "Family", "lowest")]
  paste("B and C");tax.B[intersect(d2_in, d3_in),c("Phyla", "Family", "lowest")]
  paste("All three");tax.B[intersect(d1_in, intersect(d2_in, d3_in)), c("Phyla", "Family", "lowest")]
  paste("Enriched in NucleoMag")
  paste("A and B");tax.B[intersect(d1_de, d2_de),c("Phyla", "Family", "lowest")]
  paste("A and C");tax.B[intersect(d1_de, d3_de),c("Phyla", "Family", "lowest")]
  paste("B and C");tax.B[intersect(d2_de, d3_de),c("Phyla", "Family", "lowest")]
  paste("All three");tax.B[intersect(d1_de, intersect(d2_de, d3_de)), c("Phyla", "Family", "lowest")]
  #sink()

   ### On phylum level------
  
  # Calculate the differences
  results <- list()
  
  #sink(file ="output/Indicator_analysis_ab_ab/Bacteria_overview_number_of_phylum.txt")
  for (main_factor in c("size", "kit")) {
    
    for (compost in c("K20", "K21", "K22")) {
      
      if(main_factor =="size"){
        side_factor_levels <- c("M","S")
      }
      
      if(main_factor =="kit"){
        side_factor_levels <- c("10", "2")
      }
      for (sidefactor in side_factor_levels) {
        results[[paste0(main_factor,compost, sidefactor)]] <- consistency_levels_paired(asv.B.absolut, tax.B, "Phyla", main_factor, compost, "B", min_rel = 0.001, sidefactor)
      }
    }
  }
  #sink()

  # SIZE
  
  ### SIZE
  d1 <-differ_analysis_plot(res1 ="sizeK20M", res2 ="sizeK20S", plot =FALSE)
  d2 <-differ_analysis_plot(res1 ="sizeK21M", res2 ="sizeK21S", plot =FALSE)
  d3 <-differ_analysis_plot(res1 ="sizeK22M", res2 ="sizeK22S", plot =FALSE)
  
  d1_in <-d1 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d2_in <-d2 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d3_in <-d3 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  
  d1_de <-d1 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d2_de <-d2 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d3_de <-d3 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  
  # Create a excle with all the differential abundance ASVs
  d_in <-data.frame(ASV = c(d1_in, d2_in, d3_in), factor = "2mm",
                    compost = c(rep("A", length(d1_in)), rep("B", length(d2_in)), rep("C", length(d3_in))) )
  
  d_de <-data.frame(ASV = c(d1_de, d2_de, d3_de), factor = "10mm",
                    compost = c(rep("A", length(d1_de)), rep("B", length(d2_de)), rep("C", length(d3_de))) )
  
  d_size <-rbind(d_in, d_de)
  
  d_size$compost_label <- c("K20", "K20", "K21", "K22", "K22", "K20")
  d_size <-d_size[order(d_size$compost_label),]
  # Check specific phyla and plot them
  
  list <- list()
  for (i in 1:nrow(d_size)) {
    data <-asv.B.absolut[tax.B %>% filter(Phyla ==d_size[i, "ASV"]) %>% rownames(),
                         grepl(d_size[i, "compost_label"], colnames(asv.B.absolut))] %>%
      colSums() %>% as.data.frame()
    colnames(data) <- c("count")
    data$factor <- ifelse(grepl("10", rownames(data)), "10mm", "2mm")
    
    # Wilcoxon test (non-parametric)
    test_result <- wilcox.test(count ~ factor, data = data, exact = FALSE)
    pval <- signif(test_result$p.value, 3)
    
    list[[i]] <- ggplot(data = data, aes(x = factor, y = count)) + 
      geom_boxplot() + 
      theme_classic() +
      ggtitle(paste0("Compost ", d_size[i, "compost"], ": ", d_size[i, "ASV"])) + 
      ylab("Absolute abundance") + 
      xlab("") +
      annotate("text", x = 1.5, y = max(data$count)*1.05, 
               label = paste0("p = ", pval), size = 4)
  }
  ggarrange(plotlist = list, ncol =3, nrow =2)
  #ggsave(filename = "output/Indicator_analysis_ab_ab/size/Phylum/B_Plot_phylum_diff_size.pdf", height =8, width =11)
  
  
  # KIT
  d1 <-differ_analysis_plot(res1 ="kitK2010", res2 ="kitK202", plot =FALSE)
  d2 <-differ_analysis_plot(res1 ="kitK2110", res2 ="kitK212", plot =FALSE)
  d3 <-differ_analysis_plot(res1 ="kitK2210", res2 ="kitK222", plot =FALSE)
  
  d1_in <-d1 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d2_in <-d2 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d3_in <-d3 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  
  d1_de <-d1 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d2_de <-d2 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d3_de <-d3 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  
  # Create a excle with all the differential abundance ASVs
  # No enriched in NucleoSpin
  #d_in <-data.frame(ASV = c(d1_in, d2_in, d3_in), factor = "NucleoSpin",
  #                  compost = c(rep("A", length(d1_in)), rep("B", length(d2_in)), rep("C", length(d3_in))) )
  
  d_de <-data.frame(ASV = c(d1_de, d2_de, d3_de), factor = "NucleoMag",
                    compost = c(rep("A", length(d1_de)), rep("B", length(d2_de)), rep("C", length(d3_de))) )
  
  d_kit<- d_de
  
  d_kit$compost_label <- c("K20", "K20", "K22")
  # Check specific phyla and plot them
  
  list <- list()
  for (i in 1:nrow(d_kit)) {
    data <-asv.B.absolut[tax.B %>% filter(Phyla ==d_kit[i, "ASV"]) %>% rownames(),
                         grepl(d_kit[i, "compost_label"], colnames(asv.B.absolut))] %>%
      colSums() %>% as.data.frame()
    colnames(data) <- c("count")
    data$factor <- ifelse(grepl("10", rownames(data)), "NucleoMag", "NucleoSpin")
    
    # Wilcoxon test (non-parametric)
    test_result <- wilcox.test(count ~ factor, data = data, exact = FALSE)
    pval <- signif(test_result$p.value, 3)
    
    list[[i]] <- ggplot(data = data, aes(x = factor, y = count)) + 
      geom_boxplot() + 
      theme_classic() +
      ggtitle(paste0("Compost ", d_kit[i, "compost"], ": ", d_kit[i, "ASV"])) + 
      ylab("Absolute abundance") + 
      xlab("") +
      annotate("text", x = 1.5, y = max(data$count)*1.05, 
               label = paste0("p = ", pval), size = 4)
  }
  ggarrange(plotlist = list, ncol =3, nrow =1)
  #ggsave(filename = "output/Indicator_analysis_ab_ab/kit/Phylum/B_Plot_phylum_diff_kit.pdf", height =4, width =11)
  

  
  ## FUNGI------
   ### ASV level (Table S5)-----
    # Calculate the differences
  results <- list()
  
  #sink(file ="output/Indicator_analysis_ab_ab/Fungi_overview_number_of_ASV.txt")
  for (main_factor in c("size", "kit")) {
    
    for (compost in c("K20", "K21", "K22")) {
      
      if(main_factor =="size"){
        side_factor_levels <- c("M","S")
      }
      
      if(main_factor =="kit"){
        side_factor_levels <- c("10", "2")
      }
      for (sidefactor in side_factor_levels) {
        results[[paste0(main_factor,compost, sidefactor)]] <- consistency_levels_paired(asv.F.absolut, tax.F, "ASV", main_factor, compost, "F", min_rel = 0.0001, sidefactor)
      }
    }
  }
  #sink()
  
  # Which ASVs are enriched independant of secondary factor
  unique_phyla <- unique(tax.F$Phyla)
  phylum_colors <- c(
    "#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928", "#a6cee3",
    "#b2df8a", "#fb9a99", "#fdbf6f")
  
  phylum_colors <- setNames(phylum_colors, unique_phyla)

  
  # Ploting funciton ist veraltet und muss erneurt werden
  #p1 <-differ_analysis_plot(res1 ="sizeK20M", res2 ="sizeK20S", plot =TRUE)
  #p2 <-differ_analysis_plot(res1 ="sizeK21M", res2 ="sizeK21S", plot =TRUE)
  #p3 <-differ_analysis_plot(res1 ="sizeK22M", res2 ="sizeK22S", plot =TRUE)
  
  #ggarrange(p2,p3, ncol = 3, labels =c("A", "B", "C"))
  
  #ggsave(filename = "output/Indicator_analysis_ab_ab/size/F_plot_size_diff_ab.pdf", height = 5, width = 12)
  
  #p1 <-differ_analysis_plot(res1 ="kitK2010", res2 ="kitK202", plot =TRUE)
  #p2 <-differ_analysis_plot(res1 ="kitK2110", res2 ="kitK212", plot =TRUE)
  #p3 <-differ_analysis_plot(res1 ="kitK2210", res2 ="kitK222", plot =TRUE)
  
  #ggarrange(p1, p2,p3, ncol = 3, common.legend = TRUE, labels =c("A", "B", "C"))
  
  #ggsave(filename = "output/Indicator_analysis_ab_ab/kit/F_plot_kit_diff_ab.pdf", height = 5, width = 12)
  
  
  # ASVs enriched in all?
  
  d2 <-differ_analysis_plot(res1 ="sizeK21M", res2 ="sizeK21S")
  d3 <-differ_analysis_plot(res1 ="sizeK22M", res2 ="sizeK22S")
  
  d2$ASV <- rownames(d2); rownames(d2) <-NULL; d2$group <- "sizeK21"
  d3$ASV <- rownames(d3); rownames(d3) <-NULL; d3$group <- "sizeK22"
  
  
  d <-rbind(d2, d3)
  #write.csv(d, file = "output/Indicator_analysis_ab_ab/All_differential_abundance_ASVs_fungi.csv")
  
  # Calculate summed relative abundance
  
  d_select <-d %>% select(ASV, group, consistency_flag)
  d_select$compost <- str_sub(d_select$group, -3)
  levels <-d_select$group %>% unique()
  
  asv.F.prop <-prop.table(as.matrix(asv.F), margin = 2) * 100 
  
  
  for (i in levels) {
    print(i)
    d_select_d <-d_select %>% filter(group ==i & consistency_flag =="Consistently Increased")
    print("fine or column")
    if (nrow(d_select_d) >1) {
      asv.F.prop[d_select_d$ASV, grepl(d_select_d$compost[1], colnames(asv.F.prop))] %>% colSums() %>% mean() %>% round(2) %>% print()
    } else {
      
      asv.F.prop[d_select_d$ASV, grepl(d_select_d$compost[1], colnames(asv.F.prop))]%>% mean() %>% round(2) %>% print()
    }
    
    d_select_i <-d_select %>% filter(group ==i & consistency_flag =="Consistently Decreased") 
    print("coarse or beads")
    if (nrow(d_select_i) >1) {
      asv.F.prop[d_select_i$ASV, grepl(d_select_i$compost[1], colnames(asv.F.prop))] %>% colSums() %>% mean() %>% round(2) %>% print()
    } else {
      asv.F.prop[d_select_i$ASV, grepl(d_select_i$compost[1], colnames(asv.F.prop))] %>% mean() %>% round(2) %>% print()
    }
  }
  
  rm(d_select, levels, d_select_d, d_select_i)
  
  ### SIZE
  d1 <-differ_analysis_plot(res1 ="sizeK20M", res2 ="sizeK20S")
  d2 <-differ_analysis_plot(res1 ="sizeK21M", res2 ="sizeK21S")
  d3 <-differ_analysis_plot(res1 ="sizeK22M", res2 ="sizeK22S")
  
  rbind(d1,d2, d3)
  
  d1_in <-d1 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d2_in <-d2 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d3_in <-d3 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  
  d1_de <-d1 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d2_de <-d2 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d3_de <-d3 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  
  # Create a excle with all the differential abundance ASVs
  d_in <-data.frame(ASV = c(d1_in, d2_in, d3_in), factor = "2mm",
                    compost = c(rep("A", length(d1_in)), rep("B", length(d2_in)), rep("C", length(d3_in))) )
  
  d_de <-data.frame(ASV = c(d1_de, d2_de, d3_de), factor = "10mm",
                    compost = c(rep("A", length(d1_de)), rep("B", length(d2_de)), rep("C", length(d3_de))) )
  
  d_size <-rbind(d_in, d_de)
  
  tax.F.adapt <- tax.F
  tax.F.adapt$ASV <- rownames(tax.F.adapt); rownames(tax.F.adapt) <- NULL
  
  d_size <-d_size %>% left_join(tax.F.adapt[, c("Phyla", "Family", "Genus", "ASV")], by ="ASV")
  
  summary_size <-d_size %>%
    group_by(factor, compost, Phyla) %>%
    summarise(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = Phyla, values_from = count, values_fill = 0) %>% as.data.frame()
  
  rownames(summary_size) <- paste0(summary_size$factor, "_", summary_size$compost)
  summary_size$factor <- NULL; summary_size$compost <- NULL
  
  #summary_size %>% t() %>% write.csv(file= "output/Indicator_analysis_ab_ab/size/ASV/Fungi_size_summary_ASV_phyla_all.csv")
  
    #sink(file = "output/Indicator_analysis_ab_ab/size/F_Shared_size_among_composts.txt")
  paste("Enriched in 2mm")
  paste("A and B");intersect(d1_in, d2_in)
  paste("A and C"); intersect(d1_in, d3_in)
  paste("B and C");intersect(d2_in, d3_in)
  intersect(d1_in, intersect(d2_in, d3_in))
  paste("Enriched in 10mm")
  paste("A and B"); intersect(d1_de, d2_de)
  paste("A and c");intersect(d1_de, d3_de)
  paste("B and C");intersect(d2_de, d3_de)
  paste("All three");intersect(d1_de, intersect(d2_de, d3_de))
  tax.F[intersect(d1_de, intersect(d2_de, d3_de)),]
  #sink()
  
  # KIT unnecessary, since there was no ASV with different abundance!
  
   ### On phylum level------
  
  # Calculate the differences
  results <- list()
  
  #sink(file ="output/Indicator_analysis_ab_ab/Fungi_overview_number_of_phylum.txt")
  for (main_factor in c("size", "kit")) {
    
    for (compost in c("K20", "K21", "K22")) {
      
      if(main_factor =="size"){
        side_factor_levels <- c("M","S")
      }
      
      if(main_factor =="kit"){
        side_factor_levels <- c("10", "2")
      }
      for (sidefactor in side_factor_levels) {
        results[[paste0(main_factor,compost, sidefactor)]] <- consistency_levels_paired(asv.F.absolut, tax.F, "Phyla", main_factor, compost, "F", min_rel = 0.001, sidefactor)
      }
    }
  }
  #sink()
  
  # SIZE
  
  ### SIZE
  d1 <-differ_analysis_plot(res1 ="sizeK20M", res2 ="sizeK20S", plot =FALSE)
  d2 <-differ_analysis_plot(res1 ="sizeK21M", res2 ="sizeK21S", plot =FALSE)
  d3 <-differ_analysis_plot(res1 ="sizeK22M", res2 ="sizeK22S", plot =FALSE)
  
  d1_in <-d1 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d2_in <-d2 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  d3_in <-d3 %>% filter(consistency_flag =="Consistently Increased") %>% rownames()
  
  d1_de <-d1 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d2_de <-d2 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  d3_de <-d3 %>% filter(consistency_flag =="Consistently Decreased") %>% rownames()
  
  # Create a excle with all the differential abundance ASVs
  d_in <-data.frame(ASV = c(d1_in, d2_in, d3_in), factor = "2mm",
                    compost = c(rep("A", length(d1_in)), rep("B", length(d2_in)), rep("C", length(d3_in))) )
  
  # No decrease
  
  d_size <- d_in
  
  d_size$compost_label <- c("K22", "K22")
  
  # Check specific phyla and plot them
  
  list <- list()
  for (i in 1:nrow(d_size)) {
    data <-asv.F.absolut[tax.F %>% filter(Phyla ==d_size[i, "ASV"]) %>% rownames(),
                         grepl(d_size[i, "compost_label"], colnames(asv.F.absolut))] %>%
      colSums() %>% as.data.frame()
    colnames(data) <- c("count")
    data$factor <- ifelse(grepl("10", rownames(data)), "10mm", "2mm")
    
    # Wilcoxon test (non-parametric)
    test_result <- wilcox.test(count ~ factor, data = data, exact = FALSE)
    pval <- signif(test_result$p.value, 3)
    
    list[[i]] <- ggplot(data = data, aes(x = factor, y = count)) + 
      geom_boxplot() + 
      theme_classic() +
      ggtitle(paste0("Compost ", d_size[i, "compost"], ": ", d_size[i, "ASV"])) + 
      ylab("Absolute abundance") + 
      xlab("") +
      annotate("text", x = 1.5, y = max(data$count)*1.05, 
               label = paste0("p = ", pval), size = 4)
  }
  ggarrange(plotlist = list, ncol =3, nrow =1)
  ggsave(filename = "output/Indicator_analysis_ab_ab/size/Phylum/F_Plot_phylum_diff_size.pdf", height =4, width =11)
  
  
  # KIT unnecessary since no sig. effect

# Factor-specific ASVs-----
  ## BACTERIA (not in paper)-----
  
  smr_fct_sp <- data.frame(matrix(ncol = 5, nrow= 7))
  colnames(smr_fct_sp) <- c("factor", "nASVs", "sum_rel_ab", "max_rel_ab", "median")
  smr_fct_sp$factor <- c("A", "B", "C", "2mm", "10mm", "Mag", "Spin")
  
  title = ""
  # KIT
  kit.list <- list()
  venn.factor = "kit"
  custom_colors <- c("blue", "orange")
  
  data = asv.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_kit_asv_B.pdf", height=3, width =5)
  kit.list[[1]] <-draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_kit_ISS_B.pdf", height=3, width =5)
  kit.list[[2]] <-draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data <- ISS.rob.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_kit_ISS_rob_B.pdf", height=3, width =5)
  kit.list[[3]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  # KOMPOST
  compost.list <- list()
  venn.factor = "compost"
  custom_colors <- color.compost
  
  data = asv.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_compost_asv_B.pdf", height=5, width =5)
  compost.list[[1]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  
  data = ISS.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_compost_ISS_B.pdf", height=5, width =5)
  compost.list[[2]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.rob.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_compost_ISS_rob_B.pdf", height=5, width =5)
  compost.list[[3]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  # SIZE
  size.list <-list()
  venn.factor = "size"
  custom_colors <- c("lightblue", "lightgreen")
  
  data = asv.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_size_asv_B.pdf", height=3, width =5)
  size.list[[1]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_size_ISS_B.pdf", height=3, width =5)
  size.list[[2]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.rob.B %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_size_ISS_rob_B.pdf", height=3, width =5)
  size.list[[3]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  # summary
  
  # relative abundance
  ISS.rob.prop.B <- prop.table(as.matrix(ISS.rob.B), margin = 2) * 100 
  # rank based on ASV abundance in the full data set
  sum_abundance_B <-ISS.rob.prop.B %>% rowSums() %>% as.data.frame()
  colnames(sum_abundance_B) <- c("sum_ab")
  sum_abundance_B$rank <-rank(-sum_abundance_B$sum_ab)
  
  # composts
  #sets <-compost.list[[3]]
  sets <-compost.list[[1]]
  A <-setdiff(sets$A, union(sets$B, sets$C))
  B <-setdiff(sets$B, union(sets$A, sets$C))
  C <-setdiff(sets$C, union(sets$A, sets$B))
  
  smr_fct_sp[1:3,"nASVs"] <- c(length(A), length(B), length(C))
  smr_fct_sp[1:3,"sum_rel_ab"] <-c(ISS.rob.prop.B[A,] %>% colSums() %>% mean() %>% round(3), 
                                   ISS.rob.prop.B[B,] %>% colSums() %>% mean() %>% round(3),
                                   ISS.rob.prop.B[C,] %>% colSums() %>% mean() %>% round(3))
  
  smr_fct_sp[1:3,"max_rel_ab"] <-c(ISS.rob.prop.B[A,] %>% max() %>% mean() %>% round(3), 
                                   ISS.rob.prop.B[B,] %>% max() %>% mean() %>% round(3),
                                   ISS.rob.prop.B[C,] %>% max() %>% mean() %>% round(3))
  
  smr_fct_sp[1:3,"median"] <-c(sum_abundance_B[A,"rank"] %>% median() %>% round(0), 
                               sum_abundance_B[B,"rank"] %>% median() %>% round(0),
                               sum_abundance_B[C,"rank"] %>% median() %>% round(0))
  
  # Small analysis how many are found in all three analysis
  (asv.B[intersect(intersect(sets$A, sets$B), sets$C),] %>% sum())
  
  
  # size
  sets <-size.list[[3]]
  size2 <- setdiff(sets$'2', sets$'10')
  size10 <-setdiff(sets$'10', sets$'2')
  
  smr_fct_sp[4:5,"nASVs"] <- c(length(size2), length(size10))
  smr_fct_sp[4:5,"sum_rel_ab"] <-c(ISS.rob.prop.B[size2,] %>% colSums() %>% mean() %>% round(3), 
                                   ISS.rob.prop.B[size10,] %>% colSums() %>% mean() %>% round(3))
  
  smr_fct_sp[4:5,"max_rel_ab"] <-c(ISS.rob.prop.B[size2,] %>% max() %>% mean() %>% round(3), 
                                   ISS.rob.prop.B[size10,] %>% max() %>% mean() %>% round(3))
  
  smr_fct_sp[4:5,"median"] <-c(sum_abundance_B[size2,"rank"] %>% median() %>% round(0), 
                               sum_abundance_B[size10,"rank"] %>% median() %>% round(0))
  
  
  sets <-kit.list[[3]]
  kitS <-setdiff(sets$S, sets$M)
  kitM <-setdiff(sets$M, sets$S)
  
  smr_fct_sp[6:7,"nASVs"] <- c(length(kitM), length(kitS))
  smr_fct_sp[6:7,"sum_rel_ab"] <-c(ISS.rob.prop.B[kitM,] %>% colSums() %>% mean() %>% round(3), 
                                   ISS.rob.prop.B[kitS,] %>% colSums() %>% mean() %>% round(3))
  
  smr_fct_sp[6:7,"max_rel_ab"] <-c(ISS.rob.prop.B[kitM,] %>% max() %>% round(3), 
                                   ISS.rob.prop.B[kitS,] %>% max()  %>% round(3))
  
  smr_fct_sp[6:7,"median"] <-c(sum_abundance_B[kitM,"rank"] %>% median() %>% round(0), 
                               sum_abundance_B[kitS,"rank"] %>% median() %>% round(0))
  
  #write.csv(smr_fct_sp, file="output/factor_specific_ASVs_B.csv")
  
  factor_specific_B <-list(kitS_B = kitS, kitM_B =kitM, size2_B = size2, size10_B = size10 )
  
  # Taxonomy
  
 # sink(file ="output/Factor_specific_ASVS_tax_B.txt")
  paste("All bacteria")
  tax.B$Phyla %>% table() %>% sort() %>% as.data.frame() %>% dplyr::select(Freq) %>% sum()
  
  paste("2mm and 10 mm")
  tax.B[size2,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  tax.B[size10,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  
  paste("compost A, B and C")
  tax.B[A,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  tax.B[B,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  tax.B[C,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  
  paste("Extractio with Nucleospin or Nucleomag")
  tax.B[kitS,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  tax.B[kitM,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  #sink()
  
  
  #sink(file ="output/Factor_specific_ASVS_tax_B_Genus.txt")
  paste("The 50 most common bacteria classfication on genus level")
  tax.B$Genus %>% table() %>% sort() %>% as.data.frame() %>% tail(50)
  
  paste("2mm and 10 mm")
  tax.B[size2,]$Genus %>% table() %>% sort() %>% as.data.frame()
  tax.B[size10,]$Genus %>% table() %>% sort() %>% as.data.frame()
  
  paste("Extraction with Nucleospin or Nucleomag")
  tax.B[kitS,]$Genus %>% table() %>% sort() %>% as.data.frame()
  tax.B[kitM,]$Genus %>% table() %>% sort() %>% as.data.frame()
  #sink()
  
  ## FUNGI (not in paper)----
  
  kit.list <- list()
  title = ""
  venn.factor = "kit"
  custom_colors <- c("blue", "orange")
  
  data = asv.F %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_kit_asv_F.pdf", height=3, width =5)
  kit.list[[1]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.F %>% t()
  # pdf(file="figures/venn_factor/Venn_diagram_kit_ISS_F.pdf", height=3, width =5)
  kit.list[[2]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.rob.F %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_kit_ISS_rob_F.pdf", height=3, width =5)
  kit.list[[3]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  # KOMPOST
  compost.list <- list()
  venn.factor = "compost"
  custom_colors <- color.compost
  
  data = asv.F %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_compost_asv_F.pdf", height=5, width =5)
  compost.list[[1]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.F %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_compost_ISS_F.pdf", height=5, width =5)
  compost.list[[2]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.rob.F %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_compost_ISS_rob_F.pdf", height=5, width =5)
  compost.list[[3]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  
  # SIZE
  size.list <-list()
  venn.factor = "size"
  custom_colors <- c("lightblue", "lightgreen")
  
  data = asv.F %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_size_asv_F.pdf", height=3, width =5)
  size.list[[1]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.F %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_size_ISS_F.pdf", height=3, width =5)
  size.list[[2]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  data = ISS.rob.F %>% t()
  #pdf(file="figures/venn_factor/Venn_diagram_size_ISS_rob_F.pdf", height=3, width =5)
  size.list[[3]] =draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  # summary
  
  # relative abundance
  ISS.rob.prop.F = prop.table(as.matrix(ISS.rob.F), margin = 2) * 100 
  # rank based on ASV abundance in the full data set
  sum_abundance_F <-ISS.rob.prop.F %>% rowSums() %>% as.data.frame()
  colnames(sum_abundance_F) <- c("sum_ab")
  sum_abundance_F$rank <-rank(-sum_abundance_F$sum_ab)
  
  # composts
  #sets <-compost.list[[3]]
  sets <-compost.list[[1]]
  A <-setdiff(sets$A, union(sets$B, sets$C))
  B <-setdiff(sets$B, union(sets$A, sets$C))
  C <-setdiff(sets$C, union(sets$A, sets$B))
  
  smr_fct_sp[1:3,"nASVs"] <- c(length(A), length(B), length(C))
  smr_fct_sp[1:3,"sum_rel_ab"] <-c(ISS.rob.prop.F[A,] %>% colSums() %>% mean() %>% round(3), 
                                   ISS.rob.prop.F[B,] %>% colSums() %>% mean() %>% round(3),
                                   ISS.rob.prop.F[C,] %>% colSums() %>% mean() %>% round(3))
  
  smr_fct_sp[1:3,"max_rel_ab"] <-c(ISS.rob.prop.F[A,] %>% max() %>% mean() %>% round(3), 
                                   ISS.rob.prop.F[B,] %>% max() %>% mean() %>% round(3),
                                   ISS.rob.prop.F[C,] %>% max() %>% mean() %>% round(3))
  
  smr_fct_sp[1:3,"median"] <-c(sum_abundance_F[A,"rank"] %>% median() %>% round(0), 
                               sum_abundance_F[B,"rank"] %>% median() %>% round(0),
                               sum_abundance_F[C,"rank"] %>% median() %>% round(0))
  
  
  (asv.F[intersect(intersect(sets$A, sets$B), sets$C),] %>% sum())/(asv.F %>% sum())
  # size
  sets <-size.list[[3]]
  size2 <- setdiff(sets$'2', sets$'10')
  size10 <-setdiff(sets$'10', sets$'2')
  
  smr_fct_sp[4:5,"nASVs"] <- c(length(size2), length(size10))
  smr_fct_sp[4:5,"sum_rel_ab"] <-c(ISS.rob.prop.F[size2,] %>% colSums() %>% mean() %>% round(3), 
                                   ISS.rob.prop.F[size10,] %>% colSums() %>% mean() %>% round(3))
  
  smr_fct_sp[4:5,"max_rel_ab"] <-c(ISS.rob.prop.F[size2,] %>% max() %>% mean() %>% round(3), 
                                   ISS.rob.prop.F[size10,] %>% max() %>% mean() %>% round(3))
  
  smr_fct_sp[4:5,"median"] <-c(sum_abundance_F[size2,"rank"] %>% median() %>% round(0), 
                               sum_abundance_F[size10,"rank"] %>% median() %>% round(0))
  
  
  sets <-kit.list[[3]]
  kitS <-setdiff(sets$S, sets$M)
  kitM <-setdiff(sets$M, sets$S)
  
  smr_fct_sp[6:7,"nASVs"] <- c(length(kitM), length(kitS))
  smr_fct_sp[6:7,"sum_rel_ab"] <-c(ISS.rob.prop.F[kitM,] %>% colSums() %>% mean() %>% round(3), 
                                   ISS.rob.prop.F[kitS,] %>% colSums() %>% mean() %>% round(3))
  
  smr_fct_sp[6:7,"max_rel_ab"] <-c(ISS.rob.prop.F[kitM,] %>% max() %>% mean() %>% round(3), 
                                   ISS.rob.prop.F[kitS,] %>% max() %>% mean() %>% round(3))
  
  smr_fct_sp[6:7,"median"] <-c(sum_abundance_F[kitM,"rank"] %>% median() %>% round(0), 
                               sum_abundance_F[kitS,"rank"] %>% median() %>% round(0))
  
  #write.csv(smr_fct_sp, file="output/factor_specific_ASVs_F.csv")
  
  factor_specific_F <-list(kitS_F = kitS, kitM_F =kitM, size2_F = size2, size10_F = size10 )
  
  
  # Taxonomy
  
  # Phylum
  #sink(file ="output/Factor_specific_ASVS_tax_F.txt")
  paste("All Fungi")
  tax.F$Phyla %>% table() %>% sort() %>% as.data.frame()
  
  paste("2mm and 10 mm")
  tax.F[size2,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  tax.F[size10,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  
  paste("compost A, B and C")
  tax.F[A,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  tax.F[B,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  tax.F[C,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  
  paste("Extractio with Nucleospin or Nucleomag")
  tax.F[kitS,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  tax.F[kitM,]$Phyla %>% table() %>% sort() %>% as.data.frame()
  #sink()
  
  # Genus
  
  sink(file ="output/Factor_specific_ASVS_tax_B_Funig.txt")
  paste("The 20 most common fungi classfication on genus level")
  tax.F$Genus %>% table() %>% sort() %>% as.data.frame() %>% tail(20)
  
  paste("2mm and 10 mm")
  tax.F[size2,]$Genus %>% table() %>% sort() %>% as.data.frame()
  tax.F[size10,]$Genus %>% table() %>% sort() %>% as.data.frame()
  
  paste("Extraction with Nucleospin or Nucleomag")
  tax.F[kitS,]$Genus %>% table() %>% sort() %>% as.data.frame()
  tax.F[kitM,]$Genus %>% table() %>% sort() %>% as.data.frame()
  sink()
  

  ## Factor-specific ASVs compost separately (in paper, Table S4)----
  ISS.rob.F.01 <- ifelse(ISS.rob.F>0, 1,0) %>% t() %>% as.data.frame()
  ISS.rob.F.01$sample <- rownames(ISS.rob.F.01)
  ISS.rob.F.01.des <- merge(design %>% select(compost, size), ISS.rob.F.01, by.y = "sample", by.x = "row.names")
  
  rownames(ISS.rob.F.01.des)<-ISS.rob.F.01.des$Row.names; ISS.rob.F.01.des$Row.names <-NULL
  ISS.rob.F.01.des$compostsize <- paste0(ISS.rob.F.01.des$compost, "_", ISS.rob.F.01.des$size)
  ISS.rob.F.01.des$size <- NULL; ISS.rob.F.01.des$compost <-NULL
  
  
  # For each group, calculate presence (1) if ASV is present in at least one sample
  summary_df <- ISS.rob.F.01.des %>%
    group_by(compostsize) %>%
    summarise(across(where(is.numeric), ~ as.integer(any(. == 1))), .groups = "drop") %>%
    as.data.frame()
  rownames(summary_df) <- summary_df$compostsize; summary_df$compostsize <- NULL
  
  summary_df_t <- summary_df %>% t() %>% as.data.frame()
  upset(summary_df_t,
        sets = colnames(summary_df_t),
        keep.order = TRUE,
        nintersects = 20,        # Show only top 10 intersections
        order.by = "freq")  

  
  ISS.rob.B.01 <- ifelse(ISS.rob.B>0, 1,0) %>% t() %>% as.data.frame()
  ISS.rob.B.01$sample <- rownames(ISS.rob.B.01)
  ISS.rob.B.01.des <- merge(ISS.rob.B.01, design %>% select(compost, kit, size), by.x = "sample", by.y = "row.names")
  
  
  # For each group, calculate presence (1) if ASV is present in at least one sample
  summary_df <- ISS.rob.B.01.des %>%
    group_by(compost, size) %>%
    summarise(across(where(is.numeric), ~ as.integer(any(. == 1))), .groups = "drop") %>%
    as.data.frame()
  
  summary_df$compostsize <- paste0(summary_df$compost, "_", summary_df$size)
  summary_df$size <- NULL; summary_df$compost <-NULL
  rownames(summary_df) <- summary_df$compostsize; summary_df$compostsize <- NULL
  
  summary_df_t <- summary_df %>% t() %>% as.data.frame()
  upset(summary_df_t,
        sets = colnames(summary_df_t),
        keep.order = TRUE,
        nintersects = 20,        # Show only top 10 intersections
        order.by = "freq")  
  
  # Compost specific
  
  custom_colors <- c("blue", "orange")
  ASV_names <- list()
  ASV_summary <- list()
for (compost in c("K20", "K21", "K22")) {
 
  # Bacteria size
  venn.factor <- "size"
  title <- paste("Particle size fraction",compost)
  data <- ISS.rob.B %>% select(contains(compost)) %>% t()
  #pdf(file=paste0("figures/venn_factor/Venn_diagram_size_ISS_rob_B_", compost, ".pdf"), height=3, width =5)
  results <-draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  mean_rel_ab <- rbind(c("coarse","fine","both"), c(NA, NA, NA))
  mean_rel_ab[2,1] <- ISS.rob.prop.B[rownames(ISS.rob.prop.B) %in% setdiff(results$'10', results$'2'),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  mean_rel_ab[2,2] <- ISS.rob.prop.B[rownames(ISS.rob.prop.B) %in% setdiff(results$'2', results$'10'),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  mean_rel_ab[2,3] <-ISS.rob.prop.B[rownames(ISS.rob.prop.B) %in% intersect(results$'2', results$'10'),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  ASV_summary[[paste0("B", venn.factor, compost)]] <- mean_rel_ab
  ASV_names[[paste0("B", venn.factor, compost, "2")]] <- setdiff(results$'2', results$'10') 
  ASV_names[[paste0("B", venn.factor, compost, "10")]] <- setdiff(results$'10', results$'2') 
  
  # Fungi size
  venn.factor <- "size"
  title <- paste("Particle size fraction",compost)
  data <- ISS.rob.F %>% select(contains(compost)) %>% t()
  #pdf(file=paste0("figures/venn_factor/Venn_diagram_size_ISS_rob_F_", compost, ".pdf"), height=3, width =5)
  results <-draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  mean_rel_ab <- rbind(c("coarse","fine","both"), c(NA, NA, NA))
  mean_rel_ab[2,1] <- ISS.rob.prop.F[rownames(ISS.rob.prop.F) %in% setdiff(results$'10', results$'2'),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  mean_rel_ab[2,2] <- ISS.rob.prop.F[rownames(ISS.rob.prop.F) %in% setdiff(results$'2', results$'10'),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  mean_rel_ab[2,3] <-ISS.rob.prop.F[rownames(ISS.rob.prop.F) %in% intersect(results$'2', results$'10'),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  ASV_summary[[paste0("F", venn.factor, compost)]] <- mean_rel_ab
  ASV_names[[paste0("F", venn.factor, compost, "2")]] <- setdiff(results$'2', results$'10') 
  ASV_names[[paste0("F", venn.factor, compost, "10")]] <- setdiff(results$'10', results$'2') 
  
  # Bacteria kit
  venn.factor <- "kit"
  title <- paste("Kit fraction",compost)
  data <- ISS.rob.B %>% select(contains(compost)) %>% t() 
  #pdf(file=paste0("figures/venn_factor/Venn_diagram_kit_ISS_rob_B_", compost, ".pdf"), height=3, width =5)
  results <-draw.venn.diagram(data, venn.factor, title, custom_colors)
  #dev.off()
  
  mean_rel_ab <- rbind(c("mag","spin","both"), c(NA, NA, NA))
  mean_rel_ab[2,1] <- ISS.rob.prop.B[rownames(ISS.rob.prop.B) %in% setdiff(results$M, results$S),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  mean_rel_ab[2,2] <- ISS.rob.prop.B[rownames(ISS.rob.prop.B) %in% setdiff(results$S, results$M),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  mean_rel_ab[2,3] <-ISS.rob.prop.B[rownames(ISS.rob.prop.B) %in% intersect(results$M, results$S),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  ASV_summary[[paste0("B", venn.factor, compost)]] <- mean_rel_ab
  ASV_names[[paste0("B", venn.factor, compost, "M")]] <- setdiff(results$M, results$S) 
  ASV_names[[paste0("B", venn.factor, compost, "S")]] <- setdiff(results$S, results$M) 
  
  # Fungi kit
  venn.factor <- "kit"
  title <- paste("Kit fraction",compost)
  data <- ISS.rob.F %>% select(contains(compost)) %>% t() 
  pdf(file=paste0("figures/venn_factor/Venn_diagram_kit_ISS_rob_F_", compost, ".pdf"), height=3, width =5)
  results <-draw.venn.diagram(data, venn.factor, title, custom_colors)
  dev.off()
  
  mean_rel_ab <- rbind(c("mag","spin","both"), c(NA, NA, NA))
  mean_rel_ab[2,1] <- ISS.rob.prop.F[rownames(ISS.rob.prop.F) %in% setdiff(results$M, results$S),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  mean_rel_ab[2,2] <- ISS.rob.prop.F[rownames(ISS.rob.prop.F) %in% setdiff(results$S, results$M),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  mean_rel_ab[2,3] <-ISS.rob.prop.F[rownames(ISS.rob.prop.F) %in% intersect(results$M, results$S),] %>% as.data.frame() %>%
    select(contains(compost)) %>% colSums() %>% mean() %>% round(2)
  
  ASV_summary[[paste0("F", venn.factor, compost)]] <- mean_rel_ab
  ASV_names[[paste0("F", venn.factor, compost, "M")]] <- setdiff(results$M, results$S) 
  ASV_names[[paste0("F", venn.factor, compost, "S")]] <- setdiff(results$S, results$M) 
}
  
  # Compare across composts
  
  #sink(file = "output/Factor_specific_ASVs_comparison_compost.txt")
  for (microbe in c("B", "F")) {
    for (factor in c("size", "kit")) {
      
      ifelse(factor =="size", level_list <-c("10", "2"), level_list <-c("M", "S"))
      for (level in level_list) {
        AB <-intersect(ASV_names[[paste0(microbe, factor, "K20", level)]], ASV_names[[paste0(microbe, factor, "K21", level)]]) %>% length()
        BC <-intersect(ASV_names[[paste0(microbe, factor, "K21", level)]], ASV_names[[paste0(microbe, factor, "K22", level)]]) %>% length()
        AC <-intersect(ASV_names[[paste0(microbe, factor, "K20", level)]], ASV_names[[paste0(microbe, factor, "K22", level)]]) %>% length()
        ABC <-intersect(intersect(ASV_names[[paste0(microbe, factor, "K20", level)]], ASV_names[[paste0(microbe, factor, "K22", level)]]),
                        ASV_names[[paste0(microbe, factor, "K21", level)]]) %>% length()
        print(paste0(microbe, factor, level))
        print(c(AB, BC, AC, ABC))
        }
    }
  }
  #sink()
  
  
  ## Overlap Differential abundance and factor-specific----
    ### BACTERIA-----
  factor_specific_B
  ISS.rob.V.01$sample <- rownames(ISS.rob.F.01)
  ISS.rob.F.01.des <- merge(ISS.rob.F.01, design %>% select(compost, kit, size), by.x = "sample", by.y = "row.names")
  
  
  # For each group, calculate presence (1) if ASV is present in at least one sample
  summary_df <- ISS.rob.F.01.des %>%
    group_by(compost, size) %>%
    summarise(across(where(is.numeric), ~ as.integer(any(. == 1))), .groups = "drop") %>%
    as.data.frame()
  
  summary_df$compostsize <- paste0(summary_df$compost, "_", summary_df$size)
  summary_df$size <- NULL; summary_df$compost <-NULL
  rownames(summary_df) <- summary_df$compostsize; summary_df$compostsize <- NULL
  
  library(UpSetR)
  summary_df_t <- summary_df %>% t() %>% as.data.frame()
  upset(summary_df_t,
        sets = colnames(summary_df_t),
        keep.order = TRUE,
        nintersects = 20,        # Show only top 10 intersections
        order.by = "freq")  
  
  
  ISS.rob.B.01 <- ifelse(ISS.rob.B>0, 1,0) %>% t() %>% as.data.frame()
  
  
  # Sparceity?
  ifelse(ISS.rob.B[ factor_specific_B$kitS_B,]>0, 1,0) %>% rowSums()
  ifelse(ISS.rob.B[ factor_specific_B$kitM_B,]>0, 1,0) %>% rowSums()
  ifelse(ISS.rob.B[ factor_specific_B$size2_B,]>0, 1,0) %>% rowSums()
  ifelse(ISS.rob.B[ factor_specific_B$size10_B,]>0, 1,0) %>% rowSums()
  
  diff_ab_B <- read.csv(file="output/Indicator_analysis_ab_ab/All_differential_abundance_ASVs_bacteria.csv")
  
  # SIZE
  intersect(factor_specific_B[["size2_B"]],diff_ab_B %>% filter(grepl("size", group) & log2FCmean>0) %>% pull(ASV))
  intersect(factor_specific_B[["size10_B"]],diff_ab_B %>% filter(grepl("size", group) & log2FCmean<0) %>% pull(ASV))
  
  asv.B["ASV1817",] # Only found in 10 mm fraction, but only for compost A!
  # KIT
  
  intersect(factor_specific_B[["kitS_B"]], diff_ab_B %>% filter(grepl("kit", group) & log2FCmean>0) %>% pull(ASV))
  intersect(factor_specific_B[["kitM_B"]],diff_ab_B %>% filter(grepl("kit", group) & log2FCmean<0) %>% pull(ASV))
  
  asv.B["ASV981",] # Ralstonia only found with NucleoMag and in all three composts!

  # Count differentiabl abundant ones
  
  print("2mm")
  diff_ab_B %>% filter(consistency_flag =="Consistently Increased" & grepl("size", group)) %>% pull(ASV) %>% unique() %>% length()
  print("10mm")
  diff_ab_B %>% filter(consistency_flag =="Consistently Decreased" & grepl("size", group)) %>% pull(ASV) %>% unique() %>% length()
  print("Column")
  diff_ab_B %>% filter(consistency_flag =="Consistently Increased" & grepl("kit", group)) %>% pull(ASV) %>% unique() %>% length()
  print("Beads")
  diff_ab_B %>% filter(consistency_flag =="Consistently Decreased" & grepl("kit", group)) %>% pull(ASV) %>% unique() %>% length()
  
  
  print("2mm")
  diff_ab_F %>% filter(consistency_flag =="Consistently Increased" & grepl("size", group)) %>% pull(ASV) %>% unique() %>% length()
  print("10mm")
  diff_ab_F %>% filter(consistency_flag =="Consistently Decreased" & grepl("size", group)) %>% pull(ASV) %>% unique() %>% length()
  
  # Compare abundance of differently abundant ASVs
  
  diff_ab_B_select <-diff_ab_B %>% select(ASV, countmean, group) %>% filter(grepl("size", group) & countmean >= 1000)
  
  
  K20_high <-(asv.B.absolut[diff_ab_B_select %>% filter(grepl("K20", group)) %>%pull(ASV), grepl("K20", colnames(asv.B.absolut))] %>% t() / 
    (asv.B.absolut[, grepl("K20", colnames(asv.B.absolut))] %>% colSums() %>% as.vector())) %>% t() %>% rowMeans()

  K21_high <-(asv.B.absolut[diff_ab_B_select %>% filter(grepl("K21", group)) %>%pull(ASV), grepl("K21", colnames(asv.B.absolut))] %>% t() / 
                (asv.B.absolut[, grepl("K21", colnames(asv.B.absolut))] %>% colSums() %>% as.vector())) %>% t() %>% rowMeans()
  
  K22_high <-(asv.B.absolut[diff_ab_B_select %>% filter(grepl("K22", group)) %>%pull(ASV), grepl("K22", colnames(asv.B.absolut))] %>% t() / 
                (asv.B.absolut[, grepl("K22", colnames(asv.B.absolut))] %>% colSums() %>% as.vector())) %>% t() %>% rowMeans()
  
  K20_high[order(K20_high)]*100
  K21_high[order(K21_high)]*100
  K22_high[order(K22_high)]*100
  
  tax.B[intersect(K21_high %>% names(), K22_high %>% names()),]
  
  
  # Which Bacillota are increased?
  
  diff_ab_B_Bac <-diff_ab_B %>% filter(Phyla =="Bacillota" & consistency_flag =="Consistently Decreased")
  diff_ab_B_Bac <-diff_ab_B_Bac %>% filter(group == "kitK20" | group=="kitK22") 
  diff_ab_B_Bac$Genus %>% table() %>% sort() %>% as.data.frame()
  diff_ab_B_Bac$Family %>% table() %>% sort() %>% as.data.frame()
  
  
  
    ### FUNGI------
  
  # Sparceity?
  ifelse(ISS.rob.F[ factor_specific_F$kitS_F,]>0, 1,0) %>% rowSums()
  ifelse(ISS.rob.F[ factor_specific_F$kitM_F,]>0, 1,0) %>% rowSums()
  ifelse(ISS.rob.F[ factor_specific_F$size2_F,]>0, 1,0) %>% rowSums()
  ifelse(ISS.rob.F[ factor_specific_F$size10_F,]>0, 1,0) %>% rowSums()
  
  diff_ab_F <- read.csv(file="output/Indicator_analysis_ab_ab/All_differential_abundance_ASVs_fungi.csv")
  
  # SIZE
  intersect(factor_specific_F[["size2_F"]],diff_ab_F %>% filter(grepl("size", group) & log2FCmean>0) %>% pull(ASV))
  intersect(factor_specific_F[["size10_F"]],diff_ab_F %>% filter(grepl("size", group) & log2FCmean<0) %>% pull(ASV))
  
  # For the fungi we do not see an overlap. Fungi the factor specific ASVs are very rare.
  # However, highly abundant fungi are affected by differential abundance!
  

  ## Check Bacillota perc ASVs in composts----
  
    compostA <-asv.B.absolut %>% select(contains("K20"))
    compostA <-compostA[compostA %>% rowSums() != 0,] 
    tax.B[rownames(compostA),"Phyla"] %>% table() %>% sort() %>% as.data.frame()
    (1/(tax.B[rownames(compostA),"Phyla"] %>% length()) *1565*100) %>% round(2)
    
      
    compostB <-asv.B.absolut %>% select(contains("K21"))
    compostB <-compostB[compostB %>% rowSums() != 0,] 
    tax.B[rownames(compostB),"Phyla"] %>% table() %>% sort() %>% as.data.frame()
    (1/(tax.B[rownames(compostB),"Phyla"] %>% length()) *1504*100) %>% round(2)
    
    compostC <-asv.B.absolut %>% select(contains("K22"))
    compostC <-compostC[compostC %>% rowSums() != 0,]  
    tax.B[rownames(compostC),"Phyla"]  %>% table() %>% sort() %>% as.data.frame()
    (1/(tax.B[rownames(compostC),"Phyla"] %>% length()) *1520*100) %>% round(2)

# ASV tables with taxonomy for Supplementary csv files-----
  
  # Read original tax file with classification percentages
    
    column_tax <- c("Kingdom", "Phyla", "Class", "Order", "Family", "Genus", "Species")
    tax.F.all <- read.table(file = "sequencing_data/8.all.ASV_ITSx.ITS2.UNITEv10_sh_dynamic_all.wang.taxonomy", header=F, sep='\t',row.names = 1)
    rownames <-rownames(tax.F.all)
    tax.F.all<-as.data.frame(str_split_fixed(tax.F.all[,1],';',8))
    tax.F.all$V8<-NULL
    colnames(tax.F.all)<- column_tax
    rownames(tax.F.all) <- rownames; rm(rownames)
    
    tax.B.all <- read.table(file = "sequencing_data/8.all.ASV_metaxa.silva_v138_2_var34.wang.taxonomy", header=F, sep='\t',row.names = 1)
    rownames <-rownames(tax.B.all)
    tax.B.all<-as.data.frame(str_split_fixed(tax.B.all[,1],';',8))
    tax.B.all$V8<-NULL; tax.B.all$V7<-NULL
    colnames(tax.B.all)<- column_tax[1:6]
    rownames(tax.B.all) <- rownames; rm(rownames)
    
    # Merge with ISS.rob.B and ISS.rob.F
    tax.F.all <-tax.F.all[rownames(ISS.F),]
    
    tax.B.all <-tax.B.all[rownames(ISS.B),]
    
    tax.F.all.asv <-merge(tax.F.all, ISS.F, by=0)
    tax.B.all.asv <-merge(tax.B.all, ISS.B, by=0)
    
    write.csv(tax.F.all.asv, file="output/tax.F.all.asv.csv")
    write.csv(tax.B.all.asv, file="output/tax.B.all.asv.csv")
    
    rm(tax.F.all, tax.F.all.asv, tax.B.all, tax.B.all.asv, column_tax)
    