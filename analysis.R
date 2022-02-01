# Analyse MANGA results (BETTINA module)
# https://github.com/jbathmann/pyMANGA

# Required packages
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("ggpubr")) install.packages("ggpubr")  # ggarrange

# Read simulation results (csv-files) and merge them into one dataframe
path = "./simulation_results"
files = list.files(path, pattern = "SAL_")

trees = data.frame()
for (file in files){
   d = read.csv(file = paste(path, file, "/TreeOutput/Population.csv", sep = "/"),
                sep = "\t")
   
   d$salinity = 10 * as.numeric(tail(strsplit(file, "_")[[1]], 1))
   trees = rbind(trees, d)
}

# Add attributes to dataframe
# Year: Simulation time in years
# tree_height: Tree height
# psi_height: Tree height related water potential reduction
# water_l_d: Individual tree water use in liter per day
# water_pot_l_d: Potential individual tree water use in liter per day
delta_t = 864000 # 10 days
psi_leaf = -7860000  # Pa
trees = trees %>% 
   mutate(year = time / 365.25 / 24 / 3600,
          tree_height = h_stem + 2*r_crown,
          psi_height = -9810 * (tree_height),
          water_l_d = bg_resources / delta_t * 24 * 3600 * 10^3,
          water_pot_l_d = -(psi_leaf - psi_height) / 
             ((xylem_resistance + root_surface_resistance) * pi) * 3600 * 24 * 10^3)

#### Figure 2 ####

# Create dataframe with hypothetical water use of tree with allometry
# adopted to 40 psu
trees_40 = trees %>% 
   filter(salinity == 40) %>% 
   filter(complete.cases(.))
for (sal in unique(trees$salinity)){
   print(sal)
   var_name = paste("water_l_d_", sal, sep = "")
   
   # Calculate individual tree water use with allometry of 40 psu tree
   trees_40 = trees_40 %>% rowwise() %>% 
      mutate(!!var_name := -(psi_leaf - psi_height - (-85000*sal)) / 
                ((xylem_resistance + root_surface_resistance) * pi) * 3600 * 24 * 10^3)
   print(max(trees_40[, var_name]))
}

# Transform dataframe to long-format
trees_40 = trees_40 %>% 
   gather(., Key, water_ld_40, water_l_d_0:water_l_d_80) %>% 
   rowwise() %>% 
   mutate(salinity_new = as.numeric(tail(strsplit(Key, "_")[[1]], 1)))


# Create Figure 2
x.axis.label = expression(Tree~italic(dbh)~~(cm))
y.axis.label = expression(Individual~tree~water~use~(L~H[2]~O~day^{-1}))
col.label = expression(Porewater~salinity~(psu))
font.size = 12

x11()

fig.A = trees %>% 
         ggplot(., aes(x = r_stem*2*100, y = water_l_d,
                       col = salinity, group = salinity)) +
         geom_line(size = 1) +
         geom_line(trees_40, mapping=aes(x = r_stem*2*100, y = water_ld_40,
                                         col = salinity_new, group = salinity_new), 
                   linetype = "dashed") +
         scale_color_gradient2(low = "royalblue3", mid = "deeppink4",
                               high = "tomato2", midpoint = 35) +
         geom_segment(aes(x = 0, xend = 75, y = 0, yend = 1.6*75),
                      col = "#222222", size = 1.2) +
         annotate('text', size = 3,
                  x = 65, y = 75, label = expression(1.6~L~day^-1~cm^-1),
                  col = "#222222") +
         labs(x = x.axis.label,
              y = y.axis.label,
              col = col.label) +
         theme_classic() +
         theme(text = element_text(size = font.size))

fig.B = trees %>%
         filter(time == max(time)) %>% 
         mutate(rel_0_salinity = water_l_d / max(water_l_d) * 100,
                act.pot = water_l_d / water_pot_l_d) %>% 
         ggplot(., aes(x = salinity, y = water_l_d, fill = salinity)) +
         geom_smooth(method = "lm", se=T, color="darkgrey", formula = y ~ x) +
         geom_point(size = 3, shape = 21, alpha = 0.8) +
         stat_regline_equation(label.y = 120, label.x.npc = 0.17, size = 3) +
         stat_regline_equation(label.y = 130, aes(label=..adj.rr.label..),
                               label.x.npc = 0.17, size = 3) +
         scale_fill_gradient2(low = "royalblue3", mid = "deeppink4",
                              high = "tomato2", midpoint = 35) +
         labs(x = col.label,
              y = y.axis.label,
              col = col.label) +
         theme_classic() +
         theme(text = element_text(size = font.size))


ggarrange(fig.A, fig.B, ncol = 2, labels = c("(a)", "(b)"),
          align = "hv",
          legend = "bottom", common.legend = T)

# Save figure
ggsave(filename = "Fig_2.jpg",
       width = 10, height = 5, dpi = 900)


# Alternative figure 2: second axis shows water use of each tree relative to
# water use of tree growing at 0 psu

f = trees %>% filter(time == max(time)) %>% 
   distinct(m = max(water_l_d))
fig.B.alt = fig.B +
   scale_y_continuous(sec.axis = sec_axis(~ . / f$m*100,
                                          name = "Relative water use (%)")) 

ggarrange(fig.A, fig.B.alt, ncol = 2, labels = c("(a)", "(b)"),
          align = "hv",
          legend = "bottom", common.legend = T)

ggsave(filename = "Fig_2_alt.jpg",
       width = 10, height = 5, dpi = 900)


#### Supplementary Figure 2 ####

x11()
trees %>% 
   filter(time == max(time)) %>% 
   mutate(act.pot = water_l_d / water_pot_l_d) %>% 
   ggplot(., aes(x = salinity, y = act.pot)) +
   geom_smooth(method = "lm", se=T, color="darkgrey", formula = y ~ x) +
   stat_regline_equation(label.y = 0.7) +
   stat_regline_equation(label.y = 0.8, 
                         aes(label=..adj.rr.label..)) +
   geom_point(aes(size = tree_height, col = 2*r_stem*100)) +
   scale_color_continuous(low = "lightgrey", high = "black") +
   labs(x = col.label,
        y = "Actual:Potential transpiration",
        size = "Tree height (m)",
        col = x.axis.label) +
   theme_classic() +
   theme(text = element_text(size = font.size))


ggsave(filename = "Fig_S2.jpeg",
       width = 7, height = 5, dpi = 900)
