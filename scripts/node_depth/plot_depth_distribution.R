#!/usr/bin/env Rscript
library(ggplot2)
library(ggpmisc)

args<-commandArgs(T)
file_i = args[1]

data = read.table(file_i)
data$V1 = as.numeric(data$V1)
data = subset(data, V1<200)

p = ggplot(data, aes(x = V1)) +
geom_density(color = 'black')

pb <- ggplot_build(p)

x = as.data.frame(pb[["data"]][[1]][["x"]])
y = as.data.frame(pb[["data"]][[1]][["density"]])
z = cbind(x,y)
colnames(z) = c('x','y')

q = ggplot(z, aes(x=x, y=y)) +
  geom_line() +
  stat_peaks(span = 10, colour = "red", size=3) +
  stat_peaks(span = 10, geom = "text", colour = "red", vjust = -0.4, label.fmt = "%.1f") +
  stat_valleys(span = 10, colour = "blue", size=3) +
  stat_valleys(span = 10, geom = "text", colour = "blue", vjust = -0.4, label.fmt = "%.1f")

ggsave(q, file=paste0(file_i,".pdf"), width=6, height=5)
