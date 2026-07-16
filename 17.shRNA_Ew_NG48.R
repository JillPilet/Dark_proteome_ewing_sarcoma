### 0. Load libraries ----------------------------------------------------------
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridis)
library(tibble)
library(extrafont)
loadfonts()
library(emmeans)
library(ggpubr)
library(rstatix)

# Load qPCR dataset ------------------------------------------------------------
qpcr.dat <- read.delim(".../input_23052024.txt", sep="\t", header=T)
dim(qpcr.dat)

# Reformat the table------------------------------------------------------------
qpcr.dat = data.frame(t(qpcr.dat))
colnames(qpcr.dat) = paste0(qpcr.dat["Target",], "_", qpcr.dat["Replicate",])
qpcr.dat = qpcr.dat[3:nrow(qpcr.dat),]
qpcr.dat = data.frame(t(qpcr.dat))
qpcr.dat$probe = rownames(qpcr.dat)

qpcr.dat = qpcr.dat %>%
  rownames_to_column(var = "Row.names") %>%
  pivot_longer(.,cols = -Row.names, names_to = "Sample.long.name", values_to = "2-DeltaDeltaCT")


qpcr.dat = qpcr.dat %>%
  mutate(Sample = case_when(grepl("17985", Sample.long.name) ~ "sh17985",
                            grepl("17986", Sample.long.name) ~ "sh17986",
                            grepl("17987", Sample.long.name) ~ "sh17987",
                            grepl("CTL", Sample.long.name) ~ "shcontrol",
                            grepl("wt", Sample.long.name) ~ "ASP14 wt"),
         Treatment = case_when(grepl("untreated", Sample.long.name) ~ "untreated",
                               grepl("DOX", Sample.long.name) ~ "+DOX",
                               grepl("0.5mM", Sample.long.name) ~ "+0.5mM IPTG"),
         Target = case_when(grepl("COL13A1", Row.names) ~ "COL13A1",
                            grepl("Neoantisens", Row.names) ~ "Ew_NG48",
                            grepl("EF1", Row.names) ~ "EWSR1-FLI1",
                            grepl("RPLP0", Row.names) ~ "RPLP0",
                            grepl("Proliferation", Row.names) ~ "Proliferation"),
         Time = case_when(grepl("Day_3", Sample.long.name) ~ "72 hours",
                          grepl("Day_6", Sample.long.name) ~ "6 days"),
         Replicate_qPCR = case_when(grepl("Replicate_1", Row.names) ~ "1",
                                    grepl("Replicate_2", Row.names) ~ "2",
                                    grepl("Replicate_3", Row.names) ~ "3",
                                    grepl("No_Replicate", Row.names) ~ "No_Replicate"))

qpcr.dat$`2-DeltaDeltaCT` = as.numeric(qpcr.dat$`2-DeltaDeltaCT`)

qpcr.dat = qpcr.dat %>%
  mutate(BioRep = case_when(
    grepl("_rep1", Sample.long.name) ~ "rep1",
    grepl("_rep2", Sample.long.name) ~ "rep2",
    grepl("_rep3", Sample.long.name) ~ "rep3"
  ))

qpcr.mean = qpcr.dat %>%
  group_by(Target, Time, Sample, Treatment, BioRep) %>%
  summarise(
    mean_expression = mean(`2-DeltaDeltaCT`, na.rm = TRUE),
    .groups = "drop"
  )


# plot Ew_NG48 expression after 72h---------------------------------------------
Ew_NG48.72h = qpcr.mean %>%
  filter(Target == "Ew_NG48", Time == "72 hours")

table(Ew_NG48.72h$Treatment)
Ew_NG48.72h$Treatment = factor(Ew_NG48.72h$Treatment, levels=c("untreated", "+DOX", "+0.5mM IPTG"))
Ew_NG48.72h$Sample = factor(Ew_NG48.72h$Sample, levels=c("ASP14 wt", "shcontrol", "sh17985", "sh17986", "sh17987"))

p <- ggplot(Ew_NG48.72h,
            aes(x = Sample,
                y = mean_expression,
                fill = Sample)) +
  
  # Barplot
  stat_summary(
    fun = mean,
    geom = "bar",
    color = "black",
    width = 0.7
  ) +
  
  # SEM error bars
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.5
  ) +
  
  # Biological replicates only (n = 3)
  geom_jitter(
    width = 0.08,
    alpha = 0.8,
    shape = 21,
    color = "black",
    size = 2
  ) +
  
  # Colors
  scale_fill_manual(values = c(
    "#ADB7C0FF",
    "#94C5CCFF",
    "#F4ADB3FF",
    "#EEBCB1FF",
    "#ECD89DFF"
  )) +
  
  # Facets
  facet_grid(cols = vars(Treatment)) +
  
  # Theme
  theme_classic() +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 10,
      angle = 90,
      hjust = 1,
      family = "Helvetica"
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 10,
      family = "Helvetica",
      margin = margin(5,0,0,0)
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.2, "cm"),
    
    strip.text.x = element_text(
      size = 10,
      family = "Helvetica"
    ),
    
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5,
      linetype = "solid"
    )
  ) +
  
  labs(
    x = "",
    y = expression(
      "Ew_NG48 expression (" * 2^{-Delta * Delta * Ct} * ")"
    )
  )

p

ggsave(plot=p, device = "pdf", dpi = 320, units = "cm", width =10, height = 6, 
       filename  = ".../all_sh_72h_Ew_NG48_papier_dark_proteome.pdf")

# Log2 transform for statistics-------------------------------------------------

Ew_NG48.72h = Ew_NG48.72h %>%
  mutate(log2_expression = log2(mean_expression))

# Stats ------------------------------------------------------------------------
stat.test = Ew_NG48.72h %>%
  group_by(Treatment) %>%
  t_test(
    log2_expression ~ Sample,
    comparisons = list(
      c("shcontrol", "ASP14 wt"),
      c("shcontrol", "sh17985"),
      c("shcontrol", "sh17986"),
      c("shcontrol", "sh17987")
    ),
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj")

stat.test

stat.test %>%
  select(
    Treatment,
    group1,
    group2,
    p,
    p.adj,
    p.adj.signif
  )

# Stats per construct group ----------------------------------------------------
table(Ew_NG48.72h$Treatment)
stat.test = Ew_NG48.72h %>%
  group_by(Sample) %>%
  t_test(
    log2_expression ~ Treatment,
    comparisons = list(
      c("untreated", "+DOX"),
      c("untreated", "+0.5mM IPTG"),
      c("+DOX", "+0.5mM IPTG")
    ),
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj")

stat.test

stat.test %>%
  select(
    Sample,
    group1,
    group2,
    p,
    p.adj,
    p.adj.signif
  )

#  Ew_NG48 6days----------------------------------------------------------------

Ew_NG48.6days = qpcr.mean %>%
  filter(Target == "Ew_NG48", Time == "6 days")

table(Ew_NG48.6days$Treatment)
Ew_NG48.6days$Treatment = factor(Ew_NG48.6days$Treatment, levels=c("untreated", "+DOX", "+0.5mM IPTG"))
Ew_NG48.6days$Sample = factor(Ew_NG48.6days$Sample, levels=c("ASP14 wt", "shcontrol", "sh17985", "sh17986", "sh17987"))

p <- ggplot(Ew_NG48.6days,
            aes(x = Sample,
                y = mean_expression,
                fill = Sample)) +
  
  # Barplot
  stat_summary(
    fun = mean,
    geom = "bar",
    color = "black",
    width = 0.7
  ) +
  
  # SEM error bars
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.5
  ) +
  
  # Biological replicates only (n = 3)
  geom_jitter(
    width = 0.08,
    alpha = 0.8,
    shape = 21,
    color = "black",
    size = 2
  ) +
  
  # Colors
  scale_fill_manual(values = c(
    "#ADB7C0FF",
    "#94C5CCFF",
    "#F4ADB3FF",
    "#EEBCB1FF",
    "#ECD89DFF"
  )) +
  
  # Facets
  facet_grid(cols = vars(Treatment)) +
  
  # Theme
  theme_classic() +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 10,
      angle = 90,
      hjust = 1,
      family = "Helvetica"
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 10,
      family = "Helvetica",
      margin = margin(5,0,0,0)
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.2, "cm"),
    
    strip.text.x = element_text(
      size = 10,
      family = "Helvetica"
    ),
    
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5,
      linetype = "solid"
    )
  ) +
  
  labs(
    x = "",
    y = expression(
      "Ew_NG48 expression (" * 2^{-Delta * Delta * Ct} * ")"
    )
  )

p

ggsave(plot=p, device = "pdf", dpi = 320, units = "cm", width =10, height = 6, 
       filename  = ".../all_sh_6days_Ew_NG48_papier_dark_proteome.pdf")

# Log2 transform for statistics-------------------------------------------------

Ew_NG48.6days = Ew_NG48.6days %>%
  mutate(log2_expression = log2(mean_expression))

# Stats ------------------------------------------------------------------------
stat.test = Ew_NG48.6days %>%
  group_by(Treatment) %>%
  t_test(
    log2_expression ~ Sample,
    comparisons = list(
      c("shcontrol", "ASP14 wt"),
      c("shcontrol", "sh17985"),
      c("shcontrol", "sh17986"),
      c("shcontrol", "sh17987")
    ),
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj")

stat.test

stat.test %>%
  select(
    Treatment,
    group1,
    group2,
    p,
    p.adj,
    p.adj.signif
  )

# Stats per construct group ----------------------------------------------------
table(Ew_NG48.6days$Treatment)
stat.test = Ew_NG48.6days %>%
  group_by(Sample) %>%
  t_test(
    log2_expression ~ Treatment,
    comparisons = list(
      c("untreated", "+DOX"),
      c("untreated", "+0.5mM IPTG"),
      c("+DOX", "+0.5mM IPTG")
    ),
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj")

stat.test

stat.test %>%
  select(
    Sample,
    group1,
    group2,
    p,
    p.adj,
    p.adj.signif
  )

# plot proliferation 72h--------------------------------------------------------
Prolif.72h = qpcr.mean %>%
  filter(Target == "Proliferation", Time == "72 hours")

table(Prolif.72h$Treatment)
Prolif.72h$Treatment = factor(Prolif.72h$Treatment, levels=c("untreated", "+DOX", "+0.5mM IPTG"))
Prolif.72h$Sample = factor(Prolif.72h$Sample, levels=c("ASP14 wt", "shcontrol", "sh17985", "sh17986", "sh17987"))

p <- ggplot(Prolif.72h,
            aes(x = Sample,
                y = mean_expression,
                fill = Sample)) +
  
  # Barplot
  stat_summary(
    fun = mean,
    geom = "bar",
    color = "black",
    width = 0.7
  ) +
  
  # SEM error bars
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.5
  ) +
  
  # Biological replicates only (n = 3)
  geom_jitter(
    width = 0.08,
    alpha = 0.8,
    shape = 21,
    color = "black",
    size = 2
  ) +
  
  # Colors
  scale_fill_manual(values = c(
    "#ADB7C0FF",
    "#94C5CCFF",
    "#F4ADB3FF",
    "#EEBCB1FF",
    "#ECD89DFF"
  )) +
  
  # Facets
  facet_grid(cols = vars(Treatment)) +
  
  # Theme
  theme_classic() +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 10,
      angle = 90,
      hjust = 1,
      family = "Helvetica"
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 10,
      family = "Helvetica",
      margin = margin(5,0,0,0)
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.2, "cm"),
    
    strip.text.x = element_text(
      size = 10,
      family = "Helvetica"
    ),
    
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5,
      linetype = "solid"
    )
  ) +
  
  labs(
    x = "",
    y = expression(
      "Proliferation"
    )
  )

p

ggsave(plot=p, device = "pdf", dpi = 320, units = "cm", width =10, height = 6, 
       filename  = ".../all_sh_72h_Proliferation_papier_dark_proteome.pdf")

# Log2 transform for statistics-------------------------------------------------

Prolif.72h = Prolif.72h %>%
  mutate(log2_expression = log2(mean_expression))

# Stats ------------------------------------------------------------------------
stat.test = Prolif.72h %>%
  group_by(Treatment) %>%
  t_test(
    log2_expression ~ Sample,
    comparisons = list(
      c("shcontrol", "ASP14 wt"),
      c("shcontrol", "sh17985"),
      c("shcontrol", "sh17986"),
      c("shcontrol", "sh17987")
    ),
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj")

stat.test

stat.test %>%
  select(
    Treatment,
    group1,
    group2,
    p,
    p.adj,
    p.adj.signif
  )

# Stats per construct group ----------------------------------------------------
table(Prolif.72h$Treatment)
stat.test = Prolif.72h %>%
  group_by(Sample) %>%
  t_test(
    log2_expression ~ Treatment,
    comparisons = list(
      c("untreated", "+DOX"),
      c("untreated", "+0.5mM IPTG"),
      c("+DOX", "+0.5mM IPTG")
    ),
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj")

stat.test

stat.test %>%
  select(
    Sample,
    group1,
    group2,
    p,
    p.adj,
    p.adj.signif
  )

#  Prolif 6days-----------------------------------------------------------------

Prolif.6d = qpcr.mean %>%
  filter(Target == "Proliferation", Time == "6 days")

table(Prolif.6d$Treatment)
Prolif.6d$Treatment = factor(Prolif.6d$Treatment, levels=c("untreated", "+DOX", "+0.5mM IPTG"))
Prolif.6d$Sample = factor(Prolif.6d$Sample, levels=c("ASP14 wt", "shcontrol", "sh17985", "sh17986", "sh17987"))

p <- ggplot(Prolif.6d,
            aes(x = Sample,
                y = mean_expression,
                fill = Sample)) +
  
  # Barplot
  stat_summary(
    fun = mean,
    geom = "bar",
    color = "black",
    width = 0.7
  ) +
  
  # SEM error bars
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.5
  ) +
  
  # Biological replicates only (n = 3)
  geom_jitter(
    width = 0.08,
    alpha = 0.8,
    shape = 21,
    color = "black",
    size = 2
  ) +
  
  # Colors
  scale_fill_manual(values = c(
    "#ADB7C0FF",
    "#94C5CCFF",
    "#F4ADB3FF",
    "#EEBCB1FF",
    "#ECD89DFF"
  )) +
  
  # Facets
  facet_grid(cols = vars(Treatment)) +
  
  # Theme
  theme_classic() +
  
  theme(
    legend.position = "none",
    
    axis.text.x = element_text(
      size = 10,
      angle = 90,
      hjust = 1,
      family = "Helvetica"
    ),
    
    axis.title.x = element_blank(),
    
    axis.title.y = element_text(
      size = 10,
      family = "Helvetica",
      margin = margin(5,0,0,0)
    ),
    
    axis.line = element_line(
      colour = "black",
      linewidth = 0.5
    ),
    
    axis.ticks.length = unit(0.2, "cm"),
    
    strip.text.x = element_text(
      size = 10,
      family = "Helvetica"
    ),
    
    strip.background = element_rect(
      colour = "black",
      fill = "white",
      linewidth = 0.5,
      linetype = "solid"
    )
  ) +
  
  labs(
    x = "",
    y = expression(
      "Proliferation"
    )
  )

p

ggsave(plot=p, device = "pdf", dpi = 320, units = "cm", width =10, height = 6, 
       filename  = ".../all_sh_6days_proliferation_papier_dark_proteome.pdf")

# Log2 transform for statistics-------------------------------------------------

Prolif.6d = Prolif.6d %>%
  mutate(log2_expression = log2(mean_expression))

# Stats per treatment group ----------------------------------------------------
stat.test = Prolif.6d %>%
  group_by(Treatment) %>%
  t_test(
    log2_expression ~ Sample,
    comparisons = list(
      c("shcontrol", "ASP14 wt"),
      c("shcontrol", "sh17985"),
      c("shcontrol", "sh17986"),
      c("shcontrol", "sh17987")
    ),
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj")

stat.test

stat.test %>%
  select(
    Treatment,
    group1,
    group2,
    p,
    p.adj,
    p.adj.signif
  )

# Stats per construct group ----------------------------------------------------
table(Prolif.6d$Treatment)
stat.test = Prolif.6d %>%
  group_by(Sample) %>%
  t_test(
    log2_expression ~ Treatment,
    comparisons = list(
      c("untreated", "+DOX"),
      c("untreated", "+0.5mM IPTG"),
      c("+DOX", "+0.5mM IPTG")
    ),
    p.adjust.method = "BH"
  ) %>%
  add_significance("p.adj")

stat.test

stat.test %>%
  select(
    Sample,
    group1,
    group2,
    p,
    p.adj,
    p.adj.signif
  )
