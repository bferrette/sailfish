library(tidyverse)
library(patchwork)

ATL <- read.table('ATL_BSP.tsv', header=TRUE, sep = "\t")
IDWP <- read.table('IDWP_BSP.tsv', header=TRUE, sep = "\t")

p01 <- ggplot(ATL, aes(x=time)) +
       geom_line(aes(y=median), linetype='solid') +
       geom_ribbon(aes(ymin=lower, ymax=upper), fill="blue",
                   size=.2, alpha=.25)+
       scale_x_continuous(trans='log10', labels = scales::comma) +
       scale_y_continuous(trans='log10') +
       labs(title = "Atlantic", x = "Mya", y = "Population") +
       theme_bw()
p01

p02 <- ggplot(IDWP, aes(x=time)) +
       geom_line(aes(y=median), linetype='solid') +
       geom_ribbon(aes(ymin=lower, ymax=upper),
                   size=.2, alpha=.25, fill="green") +
       scale_x_continuous(trans='log10', labels = scales::comma) +
       scale_y_continuous(trans='log10') +
       labs(title = "Indo-Pacific", x = "Mya", y = "Population") +
       theme_bw()
p02
# Patchwork
(p01 / p02) +  plot_annotation(tag_levels = 'a')
# save plot in '.svg' format
ggsave("SFA_BSP.png", device = "png", width = 25, height = 15, units = "cm", dpi = 600)
ggsave("SFA_BSP.svg", device = "svg", width = 25, height = 15, units = "cm", dpi = 600)

p1 <- ggplot(ATL, aes(x=Time, y=log(Median)))+
      geom_line()+
      geom_ribbon(aes(x=Time, ymax=log(Upper),ymin=log(Lower)),
                  alpha=0.2, fill="red")+
      theme_bw()
p1

p2 <- ggplot(INWP, aes(x=Time, y=log(Median)))+
      geom_line()+
      geom_ribbon(alpha=0.2,fill="green",
                  aes(x=Time, ymax=log(Upper),
                  ymin=log(Lower)))+
      theme_bw()
p2

p3 <- ggplot(ECP, aes(x=Time, y=log(Median)))+
      geom_line()+
      geom_ribbon(alpha=0.2,fill="purple",
                  aes(x=Time, ymax=log(Upper),
                  ymin=log(Lower)))+
      theme_bw()
p3

p4 <- ggplot(LAT, aes(x=Time, y=log(Median)))+
      geom_line()+
      geom_ribbon(alpha=0.2,fill="orange",
              aes(x=Time, ymax=log(Upper),
                  ymin=log(Lower)))+
      theme_bw()
p4

p5 <- ggplot(LGL, aes(x=Time, y=log(Median)))+
      geom_line()+
      geom_ribbon(alpha=0.2,fill="blue",
                  aes(x=Time, ymax=log(Upper),
                  ymin=log(Lower)))+
      theme_bw()
p5

# Patchwork
(p1 + p2 + p3) / (p4 + p5)+  plot_annotation(tag_levels = 'a')
# save plot in '.svg' format
ggsave("SFA_BSP.svg", device = "svg", width = 50, height = 25, units = "cm", dpi = 600)
