# Read WHO FluNet human statistics:
flunet_orig <- read.table('FluNet_2010-1_2021-34.csv', sep=',', header=T)
flunet_orig$SDATE <- as.Date(flunet_orig$SDATE, format="%m/%d/%y")
flunet_orig$EDATE <- as.Date(flunet_orig$EDATE, format="%m/%d/%y")

# Read Iowa State University FLUture swine statistics:
fluture_daily <- read.table("FLUture_2022-10-13_daily.csv",
                            header=T, sep=",")
fluture_daily$date <- as.Date(fluture_daily$date, format="%m/%d/%y")

# Aggregation by seasons.
seasons <- c(10:20)  # 2010-11 up to 2020-21

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

aggregate_flunet <- function(start_date, end_date) {
  start_week = -1
  accumulated <- data.frame(start_date=as.Date(character()), end_date=as.Date(character()),
                            pdm=integer(), h3=integer(), ibv=integer(),
                            unsubtyped=integer(), pdm_adj=integer())
  
  for (i in 1:nrow(flunet_orig)) {
    week_report <- flunet_orig[i,]
    week <- week_report$Week
    if (week_report$SDATE < start_date || week_report$EDATE > end_date) {
      next
    }
    
    if (start_week < 0) {
      accumulated <- data.frame(start_date=c(week_report$SDATE),
                                end_date=c(week_report$EDATE),
                                pdm=c(week_report$AH1N12009),
                                h3=c(week_report$AH3),
                                ibv=c(week_report$INF_B),
                                unsubtyped=c(week_report$ANOTSUBTYPED),
                                pdm_adj=c(0))
      start_week <- week
    } else {
      accumulated$end_date <- week_report$EDATE
      accumulated$pdm <- accumulated$pdm + week_report$AH1N12009
      accumulated$h3 <- accumulated$h3 + week_report$AH3
      accumulated$unsubtyped <- accumulated$unsubtyped + week_report$ANOTSUBTYPED
      accumulated$ibv <- accumulated$ibv + week_report$INF_B
    }
  }
  accumulated
}

aggregate_fluture <- function(start_date, end_date) {
  accumulated <- data.frame(sw_pdm=sum(fluture_daily[fluture_daily$date >= start_date&
                                                       fluture_daily$date <= end_date,]$pdmH1))
}

detections <- data.frame()
for (season in seasons) {
  flunet_season <- aggregate_flunet(season_bounds(season, F)[1], season_bounds(season, F)[2])
  fluture_season <- aggregate_fluture(season_bounds(season, T)[1], season_bounds(season, T)[2])
  flu_season <- cbind(flunet_season, fluture_season, labels=paste(season, season+1, sep='-'))
  detections <- rbind(detections, flu_season)
}
# Account for the ubsybtyped human detections:
detections$pdm_adj <- detections$pdm + detections$unsubtyped * (detections$pdm / (detections$pdm + detections$h3))

# Compute correlation between human and swine detections stats:
cor(detections$sw_pdm, detections$pdm_adj, method="pearson")
cor(detections$sw_pdm, detections$pdm_adj, method="spearman")

# Build a regression model:
reg <- lm(sw_pdm~pdm_adj, data=detections[-c(11),]) # Predictor: human. Remove 2020-21 season

# Optional: check whether residuals are auto-correlated
acf(resid(reg))

# Build a regression graph (Figure 3A).
library(ggrepel)
plot_data <- detections[-c(11),]
ggplot(plot_data, aes(pdm_adj, sw_pdm, label=labels)) +
  geom_text_repel(size=ifelse(plot_data$labels=='20-21', 5, 5),  point.padding=0.2) +
  scale_y_continuous(name="Swine detections", limits=c(0,NA))+
  scale_x_continuous(name="Human detections")+
  geom_point(color = ifelse(plot_data$labels=='20-21', "red", "blue"),
             size=ifelse(plot_data$labels=='20-21', 3, 2.5))+
  theme(text = element_text(size=18))+
  geom_abline(slope=coef(reg)[2], intercept=coef(reg)[1])
