library(ggplot2)

setwd("C:/Users/rapha/Downloads")
df <- data.frame(read.table("site_depth.global.sampled"))

median <- median(df$V1)
mad <- mad(df$V1)

ggplot(df, aes(x = V1)) +
  geom_bar() +
  geom_vline(xintercept = median,
             color = "green",
             lty = 2) +
  geom_vline(xintercept = c(median-(mad*5),
                            median+(mad*5)),
             color = "red",
             lty = 2) +
  xlab("Site depth") +
  ylab("Count") +
  scale_x_continuous(limit = c(0, median+mad*10),
                     breaks = c(seq(0, median+mad*10, 200))) +
  theme(axis.title = element_text(size = 10),
        axis.text = element_text(size = 7))

ggsave("global_site_depth.pdf", width = 174, height = 89, units = "mm")
