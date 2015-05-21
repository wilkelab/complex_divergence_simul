rm(list=ls())

library(ggplot2)
library(grid)
library(cowplot)

a.mutants <- read.table('A_mutants.dat', head=T)
a.control <- read.table('A_Control.dat', head=T)

c.mutants <- read.table('C_mutants.dat', head=T)
c.control <- read.table('C_Control.dat', head=T)

plot.binding <- function(data) {
  p <- ggplot(data, aes(x=count, y=binding)) + 
    geom_point(alpha=0.25) + 
    geom_smooth() + 
    ylab('Binding Energy') + 
    xlab('Number of Mutations') + 
    scale_y_continuous(limits=c(-15, 10))
}

plot.a.stability <- function(data, variable) {
  data <- data.frame(x = data$count, y = data[[variable]])
  p <- ggplot(data, aes(x=x, y=y)) + 
    geom_point(alpha=0.25) + 
    geom_smooth() + 
    ylab('Stability') + 
    xlab('Number of Mutations') + 
    scale_y_continuous(limits=c(-40, 175))
}

plot.c.stability <- function(data, variable) {
  data <- data.frame(x = data$count, y = data[[variable]])
  p <- ggplot(data, aes(x=x, y=y)) + 
    geom_point(alpha=0.25) + 
    geom_smooth() + 
    ylab('Stability') + 
    xlab('Number of Mutations') + 
    scale_y_continuous(limits=c(-20, 50))
}

pa.1 <- plot.binding(a.mutants)
pa.2 <- plot.binding(a.control)
pa.3 <- plot.a.stability(a.mutants, "stability1")
pa.4 <- plot.a.stability(a.control, "stability1")

pa.manuscript <- ggdraw() + 
  draw_plot(pa.1, 0.0, 0.5, 0.5, 0.5) +
  draw_plot(pa.2, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(pa.3, 0.0, 0.0, 0.5, 0.5) +
  draw_plot(pa.4, 0.5, 0.0, 0.5, 0.5) +
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.5, 0.5), size = 18)

pc.1 <- plot.binding(c.mutants)
pc.2 <- plot.binding(c.control)
pc.3 <- plot.c.stability(c.mutants, "stability2")
pc.4 <- plot.c.stability(c.control, "stability2")

pc.manuscript <- ggdraw() + 
  draw_plot(pc.1, 0.0, 0.5, 0.5, 0.5) +
  draw_plot(pc.2, 0.5, 0.5, 0.5, 0.5) +
  draw_plot(pc.3, 0.0, 0.0, 0.5, 0.5) +
  draw_plot(pc.4, 0.5, 0.0, 0.5, 0.5) +
  draw_plot_label(c("A", "B", "C", "D"), c(0, 0.5, 0, 0.5), c(1, 1, 0.5, 0.5), size = 18)

ggsave(pa.manuscript, filename='humanizing_a.pdf', height=12, width=12, useDingbats=F)
ggsave(pc.manuscript, filename='humanizing_c.pdf', height=12, width=12, useDingbats=F)