

# ------
# SETUP
# ------

rm(list=ls())

library(broom)
library(optparse)
library(scales)

options(digits = 3, scipen = -2)

# ----------
# FUNCTIONS
# ----------

heaD <- function(x, ...){
  head(x, ...)
}

round_if <- function(x, round_place = 3){
  x_new <- vector()
  for(i in seq(x)){
    num <- as.numeric(x[i])
    if(is.na(num)){
      x_new[i] <- x[i]
    }else if(num < 0.01){
      x_new[i] <- scales::scientific(num)
    }else{
      x_new[i] <- round(num, round_place)
    }
  }
  x_new
}


# Arguments ----

option_list = list(
  make_option(c("-m", "--metadata_file"), type="character", default=NULL,
              help="input location and name of metadata file", metavar="character"),
  make_option(c("-u", "--mutations_data_file"), type="character", default=NULL,
              help="input location and name of mutations file", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="input location and name of output file", metavar="character")
  );

# make_option(c("-s", "--STUDY_ACCESSION"), type="character", default=NULL,
#             help="input study accession number", metavar="character"),

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

print("ARGUMENTS:")
print(opt)
print("---")
print(str(opt))

# ------
# PATHS
# ------

# metadata_path <- "~/Documents/metadata/"
# mutations_data_path <- "~/Documents/comp_mut/metadata/"

# ------
# FILES
# ------

# metadata_file <- paste0(metadata_path, "tb_data_18_02_2021.csv")
# ahpc_mut_file <- paste0(mutations_data_path, "novel_ahpc_mutations.txt")
# outfile <- paste0(mutations_data_path, "ahpc_regression_results.csv")

metadata_file <- opt$metadata_file
ahpc_mut_file <- opt$mutations_data_file
outfile <- opt$outfile

# -------------
# READ IN DATA
# -------------

# Read in data
metadata <- read.csv(metadata_file, header = T)
ahpc_mut <- read.delim(ahpc_mut_file, sep = "\t", header = T)

# --------
# WRANGLE
# --------

# Subset metadata
metadata <- metadata[, c("wgs_id", "genotypic_drtype", "main_lineage", "sub_lineage", "isoniazid")]
names(metadata) <- c("wgs_id", "genotypic_drtype", "main_lineage", "sub_lineage", "inh_dst")

# Split out mutation results by mutation
ahpc_mut_split <- split(ahpc_mut, ahpc_mut$change)
ahpc_mut_split <- ahpc_mut_split[sapply(ahpc_mut_split, function(x){nrow(x) >= 10})]

# Make a list storing the binary values of whether the samples have the ahpc mutations or not
meta_change_list <- matrix(nrow = nrow(metadata), ncol = length(ahpc_mut_split))
for(i in seq(ahpc_mut_split)){
  # meta_change_list[[i]] <- metadata
  meta_change_list[,i] <- ifelse(metadata$wgs_id %in% ahpc_mut_split[[i]]$wgs_id, 1, 0)
}

# cbind this data to the metadata
metadata <- cbind(metadata, setNames(data.frame(meta_change_list), names(ahpc_mut_split)))
# Tidy up the names otherwise the formulae in the glms can't handle the col names
names(metadata) <- gsub(">", "_to_", names(metadata))
names(metadata) <- gsub("-", "MINUS", names(metadata))

# Do the same for the ahpc_mut_split list
names(ahpc_mut_split) <- gsub(">", "_to_", names(ahpc_mut_split))
names(ahpc_mut_split) <- gsub("-", "MINUS", names(ahpc_mut_split))

# -------
# MODELS
# -------

# Do the models 
mutations <- names(ahpc_mut_split)
model_list <- list()
for(i in seq(mutations)){
  mod <- as.formula(sprintf("inh_dst ~ %s", mutations[i]))
  model_list[[i]] <- glm(formula = mod, data = metadata, family = binomial)
}
names(model_list) <- mutations

# for(i in model_list){
#   print(summary(i))
# }

# Tidy up and store as df
model_list_results <- data.frame()
for(i in seq(model_list)){
  mutation <- names(model_list[i])
  x <- broom::tidy(model_list[[i]])
  model_list_results <- data.frame(rbind(model_list_results, x))
}
model_list_results <- model_list_results[!(model_list_results[, "term"] == "(Intercept)"), ]
model_list_results$p.value <- round_if(model_list_results$p.value)

# ---

# Add main lineage as co-variate
mut_plus_lin_model_list <- list()
for(i in seq(mutations)){
  mod <- as.formula(sprintf("inh_dst ~ %s + main_lineage", mutations[i]))
  mut_plus_lin_model_list[[i]] <- glm(formula = mod, data = metadata, family = binomial)
}
names(mut_plus_lin_model_list) <- mutations

# for(i in mut_plus_lin_model_list){
#   print(summary(i))
# }

# Tidy up and store as df
mut_plus_lin_model_list_results <- data.frame()
for(i in seq(mut_plus_lin_model_list)){
  mutation <- names(mut_plus_lin_model_list[i])
  x <- broom::tidy(mut_plus_lin_model_list[[i]])
  mut_plus_lin_model_list_results <- data.frame(rbind(mut_plus_lin_model_list_results, x))
}
mut_plus_lin_model_list_results <- mut_plus_lin_model_list_results[!(mut_plus_lin_model_list_results[, "term"] == "(Intercept)"), ]
mut_plus_lin_model_list_results$p.value <- round_if(mut_plus_lin_model_list_results$p.value)


mut_plus_lin_model_list_results <- mut_plus_lin_model_list_results[-(grep("lineage", 
                                                                          mut_plus_lin_model_list_results$term)), ]


# Put together
nms <- names(mut_plus_lin_model_list_results)[2:ncol(mut_plus_lin_model_list_results)]
names(mut_plus_lin_model_list_results)[names(mut_plus_lin_model_list_results) %in% nms] <- paste0(nms, "_lin")
models_results <- merge(model_list_results, mut_plus_lin_model_list_results, by = "term")
models_results$p_diff <- as.numeric(models_results$p.value) - as.numeric(models_results$p.value_lin)
models_results$p_diff <- round_if(models_results$p_diff)

# ---

# Check with anova
# https://stats.stackexchange.com/questions/59879/logistic-regression-anova-chi-square-test-vs-significance-of-coefficients-ano
anova_results <- data.frame()
for(i in seq(mut_plus_lin_model_list)){
  x <- broom::tidy(anova(model_list[[i]], mut_plus_lin_model_list[[i]], test='Chisq'))
  anova_results <- data.frame(rbind(anova_results, x))
}
# Tidy up
anova_results <- anova_results[seq(2, nrow(anova_results), by = 2), ]
anova_results$mutation <- mutations
anova_results <- anova_results[, c("mutation", names(anova_results)[2:(ncol(anova_results)-1)])]

# ---

# Filter out the models that are non-sig with the lineage as a covar and that are no different with the lineage 
non_sig_models <- anova_results[anova_results$p.value >= 0.05, "mutation"]
models_results <- models_results[as.numeric(models_results[,"p.value_lin"]) < 0.05 
                                 & !(models_results[,"term"] %in% non_sig_models), ]

# Take out ones with negative estimate
models_results <- models_results[!(models_results[, "estimate_lin"] < 0), ]

# Clean
models_results$term <- gsub("_to_", ">", models_results$term)
models_results$term <- gsub("MINUS", "-", models_results$term)

write.csv(outfile, row.names = F)










