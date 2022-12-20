library(ggplot2)
library(reshape2)

HI_table <- read.table("GK_HI_avg_distance_July2022.csv", header=T,
                              sep=',', check.names = F, allowEscapes=T)
HI_table$strain <- factor(HI_table$strain, levels=HI_table$strain)
HI_list <- melt(HI_table, id="strain")
p <- ggplot(HI_list, aes(x=variable, fill=strain, y=value)) +
    # scale_fill_manual(values=palette()[1:6])+
    # scale_colour_manual(values=palette()[1:6])+
    # geom_bar(aes(colour=strain, fill=strain), stat="identity", position = position_dodge(width=0.55), width=0.5)+
    scale_shape_manual(values = c(21:25, 8))+
    geom_point(aes(color=strain, shape=strain), position=position_dodge(width=0.65),
               size=5)+
    # geom_line(aes(group=strain, color=strain), position=position_dodge(width=0.3),
    #           linetype=2)+
    labs(x=NULL, y="log2 titer ratio") +
    theme(legend.position="top", legend.title=element_blank(),
          text=element_text(size=15), axis.text.x=element_text(size=13))+
    geom_vline(xintercept = c(0.5,1.5,2.5,3.5,4.5, 5.5))
p

