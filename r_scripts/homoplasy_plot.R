
rm(list = ls())

setwd("~/Documents/comp_mut/")

source("https://raw.githubusercontent.com/GaryNapier/Packages_functions/master/Functions.R")

library(ggplot2)
library(plyr)
library(dplyr)
library(tidyr)
library(stringr)

clean_binary_table <- function(x){
  x <- data.frame(apply(x, 2, function(x){gsub("\\[\\]", NA, as.character(x))}))
  x <- data.frame(apply(x, 2, function(x){ gsub("\\[||\\]", "", as.character(x)) }))
  x <- data.frame(apply(x, 2, function(x){ gsub("\\), \\(", "; ", as.character(x)) }))
  x <- data.frame(apply(x, 2, function(x){ gsub("\\', \\'", "-", as.character(x)) }))
  x <- data.frame(apply(x, 2, function(x){ gsub("\\(||\\)", "", as.character(x)) }))
  x <- data.frame(apply(x, 2, function(x){ gsub("\\'", "", as.character(x)) }))
  x
}

find_pos <- function(x){
  stringr::str_extract(x, "-[0-9]+|[0-9]+")
}

results_path <- "results/"
lineage_path <- "../spolpred/results/"

binary_table_file <- paste0(results_path, "isoniazid_binary_table.csv")
lineage_data_file <- paste0(lineage_path, "lineage_file.csv")

binary_table <- read.csv(binary_table_file)
lineage_data <- read.csv(lineage_data_file)

# Remove all the "[]" present in the binary table caused by python's lists
binary_table <- clean_binary_table(binary_table)

# Clean lineage column 
binary_table$main_lineage <- gsub("lineage", "", binary_table$main_lineage)
binary_table$sublin <- gsub("lineage", "", binary_table$sublin)

# Trim white space
binary_table <- trimws_df(binary_table)

# Make binary cols for CM KRM and PRM
binary_table <- data.frame(cbind(binary_table, 
                   CM_bin = ifelse(is.na(binary_table$CM), "absent", "present"), 
                   KRM_bin = ifelse(is.na(binary_table$KRM), "absent", "present"), 
                   PRM_bin = ifelse(is.na(binary_table$PRM), "absent", "present"), 
                   KRM_katG_bin = ifelse(is.na(binary_table$KRM_katG), "absent", "present")))

# Clean del/ins/dup
binary_table$KRM_katG <- clean_del_ins_dup(binary_table$KRM_katG)
binary_table$PRM <- clean_del_ins_dup(binary_table$PRM)


# Clean lineage data
lineage_data$sublin <- gsub("lineage", "", lineage_data$sublin)

# Merge in sublin
binary_table <- merge(binary_table, select(lineage_data, id, sublin), by.x = "wgs_id", by.y = "id", all.x = T, sort = F)
binary_table <- rename(binary_table, sublin = sublin.y)

# KRM ----
known_data <- select(binary_table, wgs_id, sublin, KRM_katG)
# Split out
known_data <- data.frame(known_data %>% mutate(KRM_katG = strsplit(as.character(KRM_katG), "; ")) %>% unnest(KRM_katG))
# Remove NA
known_data <- subset(known_data, !(is.na(KRM_katG)))

# Get counts
known_n <- reshape2::dcast(known_data, KRM_katG ~ 'n', value.var = "KRM_katG", fun.aggregate = length)
known_sublin_n <- reshape2::dcast(reshape2::dcast(known_data, KRM_katG + sublin ~ 'n', value.var = "KRM_katG", fun.aggregate = length), 
                                  KRM_katG ~ 'n_sublin', value.var = "KRM_katG", fun.aggregate = length)

known_pivot <- merge(known_n, known_sublin_n, by = "KRM_katG", all.x = T, sort = F)

# Get pos 
known_pivot$pos <- find_pos(known_pivot$KRM_katG)
# Clean
known_pivot$pos <- as.numeric(known_pivot$pos)
# known_pivot$pos <- ifelse(known_pivot$pos < 0, -100, known_pivot$pos)
known_pivot$status <- rep("known", nrow(known_pivot))

# PRM ----
PRM_data <- select(binary_table, wgs_id, sublin, PRM)
# Split out
PRM_data <- data.frame(PRM_data %>% mutate(PRM = strsplit(as.character(PRM), "; ")) %>% unnest(PRM))
# Remove NA
PRM_data <- subset(PRM_data, !(is.na(PRM)))

# Get counts
PRM_n <- reshape2::dcast(PRM_data, PRM ~ 'n', value.var = "PRM", fun.aggregate = length)
PRM_sublin_n <- reshape2::dcast(reshape2::dcast(PRM_data, PRM + sublin ~ 'n', value.var = "PRM", fun.aggregate = length), 
                                PRM ~ 'n_sublin', value.var = "PRM", fun.aggregate = length)

PRM_pivot <- merge(PRM_n, PRM_sublin_n, by = "PRM", all.x = T, sort = F)
# Get pos 
PRM_pivot$pos <- find_pos(PRM_pivot$PRM)
# Clean
PRM_pivot$pos <- as.numeric(PRM_pivot$pos)
PRM_pivot$pos <- ifelse(PRM_pivot$pos < 0, -100, PRM_pivot$pos)
PRM_pivot$status <- rep("unknown", nrow(PRM_pivot))

# Put together
data <- setNames(rbind_force(known_pivot, PRM_pivot), c("mutation", "n", "n_sublin", "pos", "status"))

data <- subset(data, pos > 0)

homoplasy_plot <- ggplot()+
  geom_point(data = data, aes(x = pos, y = log(n_sublin), size = log(n), colour = status))+
  geom_text(data=subset(data, mutation == "katG-p.Ser315Thr"),
            aes(x = pos, y = log(n_sublin), label = "Ser315Thr"), 
            hjust = 1.2, size = 3)+
  xlab("position")+
  ylab("log(n sublineages)")+
  theme_bw()

ggsave("results/homoplasy_plot.png", homoplasy_plot, width = 1100/5, height = 700/5, units = "mm")














