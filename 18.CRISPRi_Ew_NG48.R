# Load libraries ---------------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(tibble)
library(extrafont)
loadfonts()

# Load qPCR dataset ------------------------------------------------------------
qpcr.dat <- read.delim("...input.txt", sep="\t", header=T)
dim(qpcr.dat)

# Change the format of the table------------------------------------------------
qpcr.dat = data.frame(t(qpcr.dat))
colnames(qpcr.dat) = paste0(qpcr.dat["Target",], "_", qpcr.dat["Replicate",])
qpcr.dat = qpcr.dat[3:nrow(qpcr.dat),]
qpcr.dat = data.frame(t(qpcr.dat))

qpcr.dat = qpcr.dat %>%
  rownames_to_column(var = "Row.names") %>%
  pivot_longer(.,cols = -Row.names, names_to = "Sample.long.name", values_to = "2-DeltaDeltaCT")


qpcr.dat = qpcr.dat %>%
  mutate(Target = case_when(grepl("COL13A1", Row.names) ~ "COL13A1",
                            grepl("Ew_NG48", Row.names) ~ "Ew_NG48",
                            grepl("EWSR1-FLI1", Row.names) ~ "EWSR1-FLI1",
                            grepl("RPLP0", Row.names) ~ "RPLP0",
                            grepl("Proliferation", Row.names) ~ "Proliferation"),
         Replicate_qPCR = case_when(grepl("Replicate_1", Row.names) ~ "1",
                                    grepl("Replicate_2", Row.names) ~ "2",
                                    grepl("Replicate_3", Row.names) ~ "3"),
         Treatment = case_when(Sample.long.name=="A673.non.transfected" ~ "non transfected",
                               Sample.long.name=="sg_control1" ~ "sg_control1",
                               Sample.long.name=="sg_control2" ~ "sg_control2",
                               Sample.long.name=="enhancers3.6" ~ "Enhancers",
                               Sample.long.name=="Prom1.enhancers3.6" ~ "Promoter_1_enhancers",
                               Sample.long.name=="Prom2.enhancers3.6" ~ "Promoter_2_enhancers",
                               Sample.long.name=="Prom1" ~ "Promoter_1",
                               Sample.long.name=="Prom2" ~ "Promoter_2"))

qpcr.dat$`2-DeltaDeltaCT` = as.numeric(qpcr.dat$`2-DeltaDeltaCT`)

vec.col <- paletteer::paletteer_d("trekcolors::breen2", n=8, direction =-1) 

# Ew_NG48 expression -----------------------------------------------------------
Ew_NG48 <- qpcr.dat %>%
  filter(Target=="Ew_NG48")

Ew_NG48$Treatment = factor(Ew_NG48$Treatment, levels=c("non transfected", 
                                                       "sg_control1", 
                                                       "sg_control2", 
                                                       "Enhancers", 
                                                       "Promoter_1", 
                                                       "Promoter_2",
                                                       "Promoter_1_enhancers",
                                                       "Promoter_2_enhancers"))

p <- ggplot(Ew_NG48, aes(x = Treatment, 
                         y = `2-DeltaDeltaCT`,
                         fill = Treatment)) +
  
  # Barres = moyenne
  stat_summary(fun = mean,
               geom = "bar",
               color = "black",
               width = 0.7) +
  
  # Barres d'erreur = SEM
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.2,
               linewidth = 0.5) +
  
  # Points individuels
  geom_jitter(width = 0.1,
              alpha = 0.7,
              shape = 21,
              color = "black",
              size = 2) +
  
  scale_y_continuous(limits = c(0, 150)) +
  
  scale_fill_manual(values = vec.col) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 10,
                               angle = 90,
                               hjust = 1,
                               family = "Helvetica"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10,
                                family = "Helvetica"),
    axis.line = element_line(colour = "black",
                             linewidth = 0.5),
    axis.ticks.length = unit(0.2, "cm")
  ) +
  
  labs(
    y = expression(
      "Ew_NG48 expression (" * 2^{-Delta * Delta * Ct} * ")"
    ))

p

ggsave(plot=p, device = "pdf", dpi = 320, units = "cm", width =8, height = 8, 
       filename  = ".../qPCR_Ew_NG48_dark_proteome.pdf")

# No stats because 1v1
