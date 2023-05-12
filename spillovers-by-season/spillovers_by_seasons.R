library(ggplot2)
library(viridis)
library(reshape2)

# Analyze seasonal contributions from different spillover-seasons.
analysis_name <- 'PDM-HA-USA-09-21-Mar2022-merged-cds'
log_name <- 'spillover_contributions.csv'
combined <- read.table(paste('../trees/results/', analysis_name, '-', 1, '/', log_name, sep=""), header=TRUE, sep=",", check.names = F)
for (i in 2:20) {
  log_i <- paste('../trees/results/', analysis_name, '-', i, '/', log_name, sep="")
  data_i <- read.table(log_i, header=TRUE, sep=",", check.names = F)
  combined <- rbind(combined, data_i)
}

for (i in 1:dim(combined)[2]) {
  cat(paste(names(combined)[i], median(combined[,i]), '\n'))
}

# Aggregate statistics over all 20 phylogenetic replicates.
medians_contribs <- apply(combined, 2, median)

seasons <- c("2009-10", "2010-11", "2011-12", "2012-13", "2013-14", "2014-15",
            "2015-16", "2016-17", "2017-18", "2018-19", "2019-20", "2020-21")


seasonal_df <- data.frame(season=seasons[-1])
for (i in 1:12) {
  seasonal_df = cbind(seasonal_df, rep(0.0, 11))
}
colnames(seasonal_df) <- c("season", seasons)

i = 1
for (season_i in 2:12) {
  season <- seasons[season_i]
  total_seqs = sum(medians_contribs[i:(i + season_i - 1)])
  for (contrib_i in 1:12) {
    if (contrib_i <= season_i) {
      seasonal_df[season_i - 1,contrib_i + 1] <- medians_contribs[i]
      i = i + 1
    } else {
      seasonal_df[season_i - 1,contrib_i + 1] <- 0
    }
  }
}

seasonal_melt <- melt(seasonal_df, 1)
colnames(seasonal_melt) <- c("season", "contrib_season", "contribution")

# drop out 20-21 from the legend as it does not contribute to other seasons.
seasonal_melt <- seasonal_melt[seasonal_melt$contrib_season != "2020-21",]  

ggplot(data=seasonal_melt, aes(x=season, y=contribution, fill=contrib_season))+
  geom_col(color="white") +
  scale_fill_viridis(discrete = T, option="inferno", name="") +
  labs(y = "H1pdm detections by season of spillover", x = element_blank(), color="test") +
  # scale_x_continuous(limits = c(2, 12), expand = c(0, 0.2)) +
  # scale_y_continuous(limits = c(0,1), expand = c(0, 0)) +
  # theme_minimal() +
  theme(legend.position = "bottom", 
        axis.title.y = element_text(hjust=1),
        axis.title = element_text(size=12,face="bold"),
        axis.text = element_text(size=11,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust=0.9),
        legend.title = element_text(size=12,),
        legend.text = element_text(size=12),
  ) +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
