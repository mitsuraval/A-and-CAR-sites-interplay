


# Run the classic ANOVA that assumes normality
anova_result <- aov(hydrogen_bonds ~ car_site * w_position * a_site, data = A_filtered_hydrogen_bonds)
summary(anova_result)
anova_result <- aov(hydrogen_bonds ~ car_site * w_position * a_site, data = CAR_filtered_hydrogen_bonds)
summary(anova_result)
anova_result <- aov(stack_value ~ car_site * w_position * a_site, data = stack_pair_data)
summary(anova_result)

############################################################################################

#ART-ANOVA that can handle non-normal data

library(ARTool)
library(flextable)
library(officer)

# Convert categorical variables to factors
A_filtered_hydrogen_bonds$car_site <- factor(A_filtered_hydrogen_bonds$car_site)
A_filtered_hydrogen_bonds$a_site <- factor(A_filtered_hydrogen_bonds$a_site)
A_filtered_hydrogen_bonds$w_position <- factor(A_filtered_hydrogen_bonds$w_position)
art_model <- art(hydrogen_bonds ~ car_site * a_site * w_position, data = A_filtered_hydrogen_bonds)

CAR_filtered_hydrogen_bonds$car_site <- factor(CAR_filtered_hydrogen_bonds$car_site)
CAR_filtered_hydrogen_bonds$a_site <- factor(CAR_filtered_hydrogen_bonds$a_site)
CAR_filtered_hydrogen_bonds$w_position <- factor(CAR_filtered_hydrogen_bonds$w_position)
art_model <- art(hydrogen_bonds ~ car_site * a_site * w_position, data = CAR_filtered_hydrogen_bonds)


stack_pair_data$car_site <- factor(stack_pair_data$car_site)
stack_pair_data$a_site <- factor(stack_pair_data$a_site)
stack_pair_data$w_position <- factor(stack_pair_data$w_position)
art_model <- art(stack_value ~ car_site * a_site * w_position, data = stack_pair_data)

# Perform ART ANOVA
anova_results <- anova(art_model)
print(anova_results)


#################################################################
#Post hoc comparisons with bonferroni correction


# Run ART-C for post-hoc contrasts
posthoc_car_site <- art.con(art_model, "car_site")
posthoc_a_site <- art.con(art_model, "a_site")
posthoc_w_position <- art.con(art_model, "w_position")
posthoc_car_a <- art.con(art_model, "car_site:a_site")
posthoc_car_w <- art.con(art_model, "car_site:w_position")
posthoc_a_w <- art.con(art_model, "a_site:w_position")

# Function to extract p-values and adjust them
adjust_p_values <- function(art_contrast_result) {
  contrast_table <- as.data.frame(art_contrast_result)  # Convert to dataframe
  contrast_table$p.value <- p.adjust(contrast_table$`p.value`, method = "bonferroni")  # Adjust p-values
  return(contrast_table)
}

# Apply Bonferroni correction
posthoc_car_site_adj <- adjust_p_values(posthoc_car_site)
posthoc_a_site_adj <- adjust_p_values(posthoc_a_site)
posthoc_w_position_adj <- adjust_p_values(posthoc_w_position)
posthoc_car_a_adj <- adjust_p_values(posthoc_car_a)
posthoc_car_w_adj <- adjust_p_values(posthoc_car_w)
posthoc_a_w_adj <- adjust_p_values(posthoc_a_w)

# Create a Word document and add the results
doc <- read_docx()

# Function to add tables to the document
add_table_to_doc <- function(doc, table_data, title) {
  doc <- body_add_par(doc, title, style = "heading 1")
  doc <- body_add_flextable(doc, flextable(table_data))
  return(doc)
}

# Add ANOVA results
anova_table <- as.data.frame(anova_results)
doc <- add_table_to_doc(doc, anova_table, "ART ANOVA Results")

# Add Post-hoc tables (could make into a loop)
doc <- add_table_to_doc(doc, posthoc_car_site_adj, "Post-hoc Comparisons: car_site")
doc <- add_table_to_doc(doc, posthoc_a_site_adj, "Post-hoc Comparisons: a_site")
doc <- add_table_to_doc(doc, posthoc_w_position_adj, "Post-hoc Comparisons: w_position")
doc <- add_table_to_doc(doc, posthoc_car_a_adj, "Post-hoc Comparisons: car_site × a_site")
doc <- add_table_to_doc(doc, posthoc_car_w_adj, "Post-hoc Comparisons: car_site × w_position")
doc <- add_table_to_doc(doc, posthoc_a_w_adj, "Post-hoc Comparisons: a_site × w_position")

# Save the document
output_file <- "ART_ANOVA_Posthoc_Results.docx"
print(paste("Saving results to:", output_file))
doc <- print(doc, target = output_file)

