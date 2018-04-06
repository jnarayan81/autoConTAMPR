#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

xdata <- read.table( args[1], sep="\t", header=TRUE)
head(xdata);

#odata <- read.table( args[2], sep="\t", header=TRUE)
#head(xdata);

# Load package ggplot2
library(ggplot2)

#svg(args[2],width = 8, height = 10.64)

pdf(args[2], width = 11, height = 20 )

# Plot the basic frame of the stacked bar chart.
ggplot(data = xdata, aes(x = name, y = cnt, fill = freq)) + 
    geom_bar(stat="identity") + coord_flip()

#ggplot() + geom_bar(aes(y = cnt, x = name, fill = freq), data = xdata,
#                           stat="identity")

#p1 <- ggplot(odata, aes(x = cnt, y = name))

#p1 + geom_point(aes(color = freq)) +
#  geom_line(aes(y = cnt))

dev.off()
