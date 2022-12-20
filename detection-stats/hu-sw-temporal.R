# Code for the temporal analysis in supplemental Figure S1.

period = 4 # in weeks
# Temporal bounds (Jan 2018 - November 2020):
global_start_date <- as.Date("2018-01-01")
global_end_date <- as.Date("2020-12-01")

# Aggregate FluNet data into 4-week intervals
flunet_orig <- read.table('FluNet_2010-1_2021-34.csv', sep=',', header=T)
flunet <- data.frame(start_date=as.Date(character()), end_date=as.Date(character()),
                     pdm=integer(), h3=integer(), ibv=integer(),
                     unsubtyped=integer(), pdm_adj=integer())
flunet_orig$SDATE <- as.Date(flunet_orig$SDATE, format="%m/%d/%y")
flunet_orig$EDATE <- as.Date(flunet_orig$EDATE, format="%m/%d/%y")
start_week <- -1
accumulated <- data.frame(start_date=as.Date(character()), end_date=as.Date(character()),
                          pdm=integer(), h3=integer(), ibv=integer(),
                          unsubtyped=integer(), pdm_adj=integer())
for (i in 1:nrow(flunet_orig)) {
  week_report <- flunet_orig[i,]
  week <- week_report$Week
  if (week_report$SDATE < global_start_date || week_report$SDATE > global_end_date) {
    next
  }
  
  if (start_week < 0 ||
      as.integer(week_report$SDATE - start_date) >= period * 7) {
    flunet <- rbind(flunet, accumulated)
    accumulated <- data.frame(start_date=c(week_report$SDATE),
                              end_date=c(week_report$EDATE),
                              pdm=c(week_report$AH1N12009),
                              h3=c(week_report$AH3),
                              ibv=c(week_report$INF_B),
                              unsubtyped=c(week_report$ANOTSUBTYPED),
                              pdm_adj=c(0))
    start_week <- week
    start_date <- week_report$SDATE
  } else {
    accumulated$end_date <- week_report$EDATE
    accumulated$pdm <- accumulated$pdm + week_report$AH1N12009
    accumulated$h3 <- accumulated$h3 + week_report$AH3
    accumulated$unsubtyped <- accumulated$unsubtyped + week_report$ANOTSUBTYPED
    accumulated$ibv <- accumulated$ibv + week_report$INF_B
  }
}
if (accumulated$pdm > 0) {
  accumulated$end_date <- accumulated$start_date + period * 7
  flunet[nrow(flunet) + 1,] <- accumulated[1,]
}
flunet$pdm_adj <- flunet$pdm + flunet$unsubtyped * (flunet$pdm / (flunet$pdm + flunet$h3))

# Accumulate the respective fluture (swine) data
fluture_daily <- read.table("FLUture_2022-10-13_daily.csv",
                            header=T, sep=",")
fluture_daily$date <- as.Date(fluture_daily$date, format="%m/%d/%y")
fluture <- data.frame(sw_pdm=integer())
for (i in 1:nrow(flunet)) {
  start_date <- flunet[i,]$start_date
  end_date <- flunet[i,]$end_date
  fluture[i,] <- sum(fluture_daily[fluture_daily$date >= start_date&
                                     fluture_daily$date <= end_date,]$pdmH1)
}
detections <- cbind(flunet, sw_pdm=fluture)
if (period >= 52) {
  detections <- cbind(detections, labels=paste(format(detections$start_date, '%y'),
                                               format(detections$end_date, '%y'),
                                               sep="-"))
} else {
  detections <- cbind(detections,
                      labels=format(detections$start_date+
                                      as.integer(detections$end_date - detections$start_date)/2,
                                    '%m/%y'))
}

#-------Time-plot-------
library(ggplot2)
total_pdm = sum(detections$pdm)
total_sw_pdm = sum(detections$sw_pdm)
detections <- cbind(detections, pdm_norm=detections$pdm / total_pdm)
detections <- cbind(detections, sw_pdm_norm=detections$sw_pdm / total_sw_pdm)
detections <- cbind(detections, mid_date=detections$start_date + (period/2)*7)

ggplot(detections, aes(x=mid_date, y=pdm)) +
  geom_line(aes(x=mid_date, y=pdm_norm, color='human')) +
  geom_line(aes(x=mid_date, y=sw_pdm_norm, color='swine')) +
  xlab("") +
  ylab("detection density") +
  scale_x_date(date_labels = "%b%y", breaks=detections$mid_date[seq(1, length(detections$mid_date), 2)]) +
  scale_colour_manual("", 
                      breaks = c("human", "swine"),
                      values = c("blue", "red")) +
  theme(text = element_text(size=15))
