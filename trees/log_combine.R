# Aggregate statistics generated from phylogentic/sequence analyses.
analysis_name <- 'PDM-HA-USA-09-21-Mar2022-merged-cds'  # Should coincide with the name in "run-analyses_multi.sh"

# Overall human-to-swine spillovers:
log_name <- 'spillover_stats.csv'
combined <- read.table(paste('results/', analysis_name, '-', 1, '/', log_name, sep=""), header=TRUE, sep=",")
for (i in 2:20) {
  log_i <- paste('results/', analysis_name, '-', i, '/', log_name, sep="")
  data_i <- read.table(log_i, header=TRUE, sep=",")
  combined <- rbind(combined, data_i)
}

# Overall number of human-to-swine spillovers:
range(combined$hu.to.sw)
median(combined$hu.to.sw)

#==========================

# Analyze identified variants
log_name <- 'variants.csv'
combined <- read.table(paste('results/', analysis_name, '-', 1, '/', log_name, sep=""), header=TRUE, sep=",")
for (i in 2:20) {
  log_i <- paste('results/', analysis_name, '-', i, '/', log_name, sep="")
  data_i <- read.table(log_i, header=TRUE, sep=",")
  combined <- rbind(combined, data_i)
}

variants <- as.data.frame(table(paste(combined$hu.strain, combined$spillover_season)))
# variants <- as.data.frame(table(combined$hu.strain))
top_vars <- variants[variants$Freq>=20,]
View(top_vars)
# View(combined[combined$hu.strain %in% top_vars$Var1,])

#==========================

# Analyze long hu-to-sw spillovers + plot overall spillover data
log_name <- 'long_spillovers_t1.csv'
combined <- read.table(paste('results/', analysis_name, '-', 1, '/', log_name, sep=""), header=TRUE, sep=",", check.names = F)
for (i in 2:20) {
  log_i <- paste('results/', analysis_name, '-', i, '/', log_name, sep="")
  data_i <- read.table(log_i, header=TRUE, sep=",", check.names = F)
  combined <- rbind(combined, data_i)
}
seasons <- c('9-10', '10-11', '11-12', '12-13', '13-14', '14-15', '15-16', '16-17',
             '17-18', '18-19', '19-20')
rolling_avg <- 0
# Prints out table 1 + additional statistics.
cat("season median_spills min_spills max_spills percent_long longest_spill avg_spill num_long")
for (i in 0:10) {
  cat(paste(c(seasons[i+1], median(combined[,i * 5 + 1]),
            range(combined[,i * 5 + 1]), median(combined[, i*5+2]),
            median(combined[, i*5+4]), median(combined[, i*5+5]),
            median(combined[,i * 5 + 1]) * median(combined[, i*5+2]) / 100, '\n')))
  hu_to_sw <- c(hu_to_sw, median(combined[,i * 5 + 1]))
  rolling_avg <- rolling_avg + median(combined[,i * 5 + 5]) * median(combined[, i*5+1])
}
rolling_avg / 371  # The average spillover length overall.


seasons <- c('10-11', '11-12', '12-13', '13-14', '14-15', '15-16', '16-17',
             '17-18', '18-19', '19-20')
sw_to_hu <- c(3, 0, 0, 0, 0, 2, 1, 0, 4, 2)  # From the above analysis.
hu_to_sw <- c()
for (i in 1:10) {
  hu_to_sw <- c(hu_to_sw, median(combined[,i * 5 + 1]))
}
spillovers_by_season <- data.frame(season=seasons, 'sw-to-hu'=sw_to_hu, check.names=F)
spillovers_by_season <- cbind(spillovers_by_season, 'hu-to-sw'=hu_to_sw)

library(ggplot2)
library(reshape2)
df <- melt(spillovers_by_season, id="season")

str_pad_custom <- function(labels){
  new_labels <- stringr::str_pad(labels, 15, "right")
  return(new_labels)
}

ggplot(data=df, aes(x=season, y=value, fill=variable)) +
  scale_fill_manual(labels=str_pad_custom, values=c('black', 'darkorange')) +
  geom_bar(stat="identity") +
  labs(x=element_blank(), y="spillovers")+
  theme(legend.position="bottom", legend.title=element_blank(),
        text=element_text(size=15), axis.text.x=element_text(angle=45, size=14, vjust = 0.8, hjust=1))


#==========================

# Analyze 2020-21 contributions from different spillover-seasons
log_name <- 'spillover_contributions.csv'
combined <- read.table(paste('results/', analysis_name, '-', 1, '/', log_name, sep=""), header=TRUE, sep=",", check.names = F)
for (i in 2:20) {
  log_i <- paste('results/', analysis_name, '-', i, '/', log_name, sep="")
  data_i <- read.table(log_i, header=TRUE, sep=",", check.names = F)
  combined <- rbind(combined, data_i)
}

for (i in 1:dim(combined)[2]) {
  cat(paste(names(combined)[i], median(combined[,i]), '\n'))
}
