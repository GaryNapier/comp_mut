

heaD <- function(x, ...){
  head(x, ...)
}

metadata_file <- "~/Documents/metadata/tb_data_18_02_2021.csv"
ahpc_mut_file <- "~/Documents/comp_mut/metadata/novel_ahpc_mutations.txt"

# Read in data
metadata <- read.csv(metadata_file, header = T)
ahpc_mut <- read.delim(ahpc_mut_file, sep = "\t", header = T)

# Subset metadata
metadata <- metadata[, c("wgs_id", "genotypic_drtype", "sub_lineage", "isoniazid")]
names(metadata) <- c("wgs_id", "genotypic_drtype", "lineage", "inh_dst")

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

# Do the models 
mutations <- names(ahpc_mut_split)
model_list <- list()
for(i in seq(mutations)){
  mod <- as.formula(sprintf("inh_dst ~ %s", mutations[i]))
  model_list[[i]] <- glm(formula = mod, data = metadata, family = binomial)
}
names(model_list) <- mutations

for(i in model_list){
  print(summary(i))
}











