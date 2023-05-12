library(ggplot2)
library(viridis)
library(reshape2)

surveillance <- read.table('swine-surveillance-data_2021_H1.txt', header=TRUE, sep="\t", check.names = F)
surveillance$Date <- as.Date(surveillance$Date)
surveillance <- na.omit(surveillance)

season_bounds <- function(season_year, is_swine) {
  if (is_swine) {
    # swine: November - October
    start_date <- as.Date(paste(20, season_year, "-11-01", sep=""))
    end_date <- as.Date(paste(20, season_year + 1, "-10-30", sep=""))
  } else {
    # human: October - April
    start_date <- as.Date(paste(20, season_year, "-10-01", sep=""))
    end_date <- as.Date(paste(20, season_year + 1, "-4-30", sep=""))
  }
  c(start_date, end_date)
}

count_detections <- function(season_year) {
  bounds <- season_bounds(season_year, T)
  h1s <- sum(bounds[1] <= surveillance$Date & surveillance$Date <= bounds[2])
  h1pdms <- sum(bounds[1] <= surveillance$Date & surveillance$Date <= bounds[2] &
               surveillance$H1 == 'pandemic')
  c(h1pdms / h1s, (h1s - h1pdms) / h1s)
}

seasons <- 10:20
season_names <- c("2010-11", "2011-12", "2012-13", "2013-14", "2014-15",
             "2015-16", "2016-17", "2017-18", "2018-19", "2019-20", "2020-21")

detections_by_season <- as.data.frame(t(sapply(seasons, count_detections)))
colnames(detections_by_season) <- c('sw-H1-pdms', 'sw-H1-other')
detections_by_season$seasons <- season_names

melt_detect <- melt(detections_by_season, id=c('seasons'))
melt_detect$variable <- factor(melt_detect$variable,
                               levels = c("sw-H1-other", "sw-H1-pdms"))
colnames(melt_detect) <- c('seasons', 'H1 clade', 'value')

ggplot(data=melt_detect, aes(x=seasons, y=value, fill=`H1 clade`, order=`H1 clade`))+
  geom_col(position="stack")+
  scale_fill_manual(values=c('black', 'darkorange'))+
  scale_y_continuous(limits=c(0,1), n.breaks = 6) +
  labs(y="fraction of H1pdms among swine H1s", x=element_blank())+
  theme(legend.position = "bottom",
        axis.title = element_text(size=12,face="bold"),
        axis.text = element_text(size=11,face="bold"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust=0.9),
        legend.title = element_blank(),
        legend.text = element_text(size=12),
        legend.spacing.x = unit(0.5, 'cm')
  )


