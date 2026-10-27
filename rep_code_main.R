#### libraries and functions ####

# Panel data & regression
library(plm)
library(sandwich)
library(lmtest)
library(fixest)

# Data wrangling
library(data.table)
library(dplyr)
library(plyr)
library(tidyr)
library(reshape)
library(readxl)
library(dataPreparation)

# Visualization
library(ggplot2)
library(ggforce)
library(cowplot)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(grid)

# Output / tables
library(openxlsx)
library(stargazer)
library(xtable)
library(corrplot)

# functions needed
source("functions_rep.R")

#### load data ####

data <- read.csv("em_ine_nexus_data.csv")
data <- data[-c(1)]

#### descriptive statistics (Table 1 and Figure 2) ####

dscr<-data[c("co2_pc_wb", "gini_wid", "palma", "r20", "p80p100", "p90p100", "gini_wid_disp", "palma_d", "r20_d", "p80p100d", "p90p100d","serv_managr", "gdp_pc", "pop", "urban_pop", "services", "prim_en", "over_65", "final_hous_cons", 
             "trade_gdp", "fdi_i", "man_empl", "top10_em", "cons_cf_p0p50", "cons_cf_p90p100", "cons_cf_p99p100", "inv_cf_p90p100", "inv_cf_p99p100", "str_index")]

stargazer(dscr, type="text", title="Descriptive statistics", 
          header=FALSE, # to get rid of r package output text
          iqr = T,
          no.space = FALSE, # to remove the spaces after each line of coefficients
          column.sep.width = "2pt", # to reduce column width
          font.size = "normalsize",
          covariate.labels = c("CO2 emissions per capita","Gini coefficient (pre-tax)", "Palma ratio (pre-tax)", "20:20 ratio (pre-tax)", "Income share of Top 20% (pre-tax)","Income share of Top 10% (pre-tax)",
                               "Gini coefficient (post-tax)", "Palma ratio (post-tax)", "20:20 ratio (post-tax)", "Income share of Top 20% (post-tax)","Income share of Top 10% (post-tax)",
                               "Tertiarization","GDP per capita","Population","Urban population","Services, value added",
                               "Primary enrollment","Population ages 65 and above","Final households consumption",
                               "Trade","FDI inflows","Employment in manufacturing","Carbon footprint of Top 10%",
                               "Consumption carbon footprint of Bottom 50%","Consumption carbon footprint of Top 10%","Consumption carbon footprint of Top 1%",
                               "Investment carbon footprint of Top 10%","Investment carbon footprint of Top 1%","Environmental Policy Stringency Index")          
)


corr <- cor(data[c("l_co2_pc_wb", "l_gini_wid", "l_palma", "l_r20", "l_p80p100", "l_p90p100", "l_gini_wid_disp", "l_palma_d", "l_r20_d", "l_p80p100d", "l_p90p100d","l_serv_managr", "l_gdp_pc", "l_pop", "l_urban_pop", 
                   "l_services", "l_prim_en", "l_over_65", "l_final_hous_cons", "l_trade_gdp", "l_fdi_i")], use = "complete.obs")
colnames(corr) <- c("CO2 emissions per capita (ln)","Gini coefficient  (ln) (pre-tax)", "Palma ratio (pre-tax) (ln)", "20:20 ratio (pre-tax) (ln)", "Income share of Top 20% (ln) (pre-tax)",
                    "Income share of Top 10% (ln) (pre-tax)","Gini coefficient (ln) (post-tax)", "Palma ratio (ln) (post-tax)", "20:20 ratio (ln) (post-tax)", "Income share of Top 20% (ln) (post-tax)",
                    "Income share of Top 10% (post-tax)","Tertiarization (ln)","GDP per capita (ln)","Population (ln)", 
                    "Urban population (ln)","Services, value added (ln)","Primary enrollment (ln)","Population ages 65 and above (ln)",
                    "Final households consumption (ln)","Trade (ln)","FDI inflows (ln)")
rownames(corr) <- c("CO2 emissions per capita (ln)","Gini coefficient  (ln) (pre-tax)", "Palma ratio (pre-tax) (ln)", "20:20 ratio (pre-tax) (ln)", "Income share of Top 20% (ln) (pre-tax)",
                    "Income share of Top 10% (ln) (pre-tax)","Gini coefficient (ln) (post-tax)", "Palma ratio (ln) (post-tax)", "20:20 ratio (ln) (post-tax)", "Income share of Top 20% (ln) (post-tax)",
                    "Income share of Top 10% (post-tax)","Tertiarization (ln)","GDP per capita (ln)","Population (ln)", 
                    "Urban population (ln)","Services, value added (ln)","Primary enrollment (ln)","Population ages 65 and above (ln)",
                    "Final households consumption (ln)","Trade (ln)","FDI inflows (ln)")
cor_log<-corrplot(corr, type = "upper",
                  #method = "number",
                  #col = gray.colors(100),
                  tl.col = "black", tl.srt = 45, tl.cex = 1.5
)

#### FIGURE 1 ####

# serie storiche per gruppi di paesi
paesi <- c("Italy","India","Argentina","Australia","China","Germany","Nigeria","United Kingdom","United States")
data_paesi <- data %>% filter(Country.Name %in% paesi)

data_paesi %>% dplyr::group_by(Country.Code) %>% dplyr::summarise(
  `mean` = mean(serv_managr, na.rm = TRUE),
  .groups = "drop"
) 

serie <- ggplot(data_paesi, aes(Time, serv_managr, colour = Country.Name))
serie + geom_line(size = 1.06) + labs(title = "", x = "Year", y = "Tertiarization") + theme_classic(base_size = 20) +
  theme(legend.position = "bottom", legend.justification = "center", legend.box = "horizontal") + guides(colour=guide_legend("Country"))
#pdf -> 10 x 25


#### main regression results ####

# baseline specifications 
t2 <- co2 ~ inequality * l_serv_managr + l_gdp_pc + l_pop + l_urban_pop + l_services
t3 <- co2 ~ inequality * l_serv_managr + l_gdp_pc + l_pop + l_urban_pop + l_services + l_prim_en + l_over_65
t4 <- co2 ~ inequality * l_serv_managr + l_gdp_pc + l_pop + l_urban_pop + l_services + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i

# prints tables A3,A4, 2
ine_var <-  c("l_gini_wid", "l_palma", "l_r20", "l_p80p100", "l_p90p100") #pre tax
pre_tax_baseline <- run_regressions(data, "l_co2_pc_wb", ine_var, formulas = c(t2,t3,t4), weight = F, start_year = 1980, min_obs = 1, oecd_only = F, output_type = "text")
# mean of turning points 
mean(pre_tax_baseline$all_thresholds) #2.99


# prints tables A5,A6,3
ine_var <-  c("l_gini_wid_disp", "l_palma_d", "l_r20_d", "l_p80p100d", "l_p90p100d") #post tax
post_tax_baseline <- run_regressions(data, "l_co2_pc_wb", ine_var, formulas = c(t2,t3,t4), weight = F, start_year = 1980, min_obs = 1, oecd_only = F, output_type = "text")
# mean of turning points 
mean(post_tax_baseline$all_thresholds) #2.96


#### other result ####

# NOTE
# Not exhaustive, playing with the parameters you can get all the results in the paper

# baseline specifications 
t2 <- co2 ~ inequality * l_serv_managr + l_gdp_pc + l_pop + l_urban_pop + l_services
t3 <- co2 ~ inequality * l_serv_managr + l_gdp_pc + l_pop + l_urban_pop + l_services + l_prim_en + l_over_65
t4 <- co2 ~ inequality * l_serv_managr + l_gdp_pc + l_pop + l_urban_pop + l_services + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i

# other variables used in some specifactions
# l_empl_serv_managr l_serv_empl

#table 4 and 5
t5b <- co2 ~ inequality * l_serv_managr  + l_gdp_pc + l_pop + l_urban_pop  + l_manufacturing + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i
t5c <- co2 ~ inequality * l_serv_managr  + l_gdp_pc + l_pop + l_urban_pop + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i
t5d <- co2 ~ inequality * l_serv_managr  + l_gdp_pc + l_pop + l_urban_pop  + l_services + l_manufacturing + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i
t5e <- co2 ~ inequality * l_serv_managr + l_pop + l_urban_pop + l_services + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i
t5i <- co2 ~ inequality * l_serv_managr  + l_gdp_pc + l_pop + l_urban_pop  + l_services + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i + str_index
# Excluding not High-Income countries above 2.99 
# dajeee<-frequancy_country_above_tp(input_data = data, tp = 2.99, tp_variable = "serv_managr") 
# dajeee<-dajeee%>%filter(ccc>0.2)
# data<-data %>% filter(!(income_group != "High income" & Country.Code %in% dajeee$pa))

# table A2
tA2a <- co2 ~ inequality + l_serv_managr + l_gdp_pc + l_pop + l_urban_pop + l_services
tA2b <- co2 ~ inequality + l_serv_managr + l_gdp_pc + l_pop + l_urban_pop + l_services + l_prim_en + l_over_65
tA2c <- co2 ~ inequality + l_serv_managr + l_gdp_pc + l_pop + l_services + l_urban_pop + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i

#ONLY INE LAGGED - table A.13, A.14
l3 <- co2 ~ lag(inequality,1) * l_serv_managr + l_gdp_pc + l_pop + l_urban_pop  + l_services + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i

# --- Temporal structure of the inequality-emissions relation, tables A.15, A.16 ---
# Short-run persistence: inequality at t and t-1, each moderated by tertiarization at matching time.
f_short_lagged <- co2 ~ inequality * l_serv_managr + lag(inequality,1) * lag(l_serv_managr,1) + l_gdp_pc + l_pop + l_urban_pop + l_services + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i
# Extended persistence: inequality at t, t-1, t-2, each moderated by tertiarization at matching time.
f_extended_lagged <- co2 ~ inequality * l_serv_managr + lag(inequality,1) * lag(l_serv_managr,1) + lag(inequality,2) * lag(l_serv_managr,2) + l_gdp_pc + l_pop + l_urban_pop + l_services + l_prim_en + l_over_65 + l_final_hous_cons + l_trade_gdp + l_fdi_i


# definition of variables
emission_var <- c("l_co2_pc_wb") # l_co2_pc_wb - l_eora_cons - l_eora_prod

ine_var <- c("l_gini_wid_disp", "l_palma_d", "l_r20_d", "l_p80p100d", "l_p90p100d")
# c("l_gini_wid", "l_palma", "l_r20", "l_p80p100", "l_p90p100")
# c("l_gini_wid_disp", "l_palma_d", "l_r20_d", "l_p80p100d", "l_p90p100d")

# run
results <- run_regressions(data, emission_var, ine_var, formulas = c(l3), weight = F, start_year = 1980, min_obs = 1, oecd_only = F, output_type = "text")
all_thresholds <- results$all_thresholds
all_c_a_tps <- results$all_c_a_tps


#### FIGURE 3 and 4 ####

ine_var2 <- c("l_gini_wid","l_gini_wid_disp") #l_p90p100, l_p90p100d, l_gini_wid, l_gini_wid_disp
emission_var <- c("l_co2_pc_wb")

results2 <- run_regressions2(data, emission_var, ine_var2, formulas = c(t4), weight = F, start_year = 1980, min_obs = 1, oecd_only = F, output_type = "text")

selected_models <- subset(
  results2$model_info,
  inequality_variable %in% c(
    "l_gini_wid",
    "l_gini_wid_disp"
  )
)

# Calculate marginal effects for the selected pre-tax and post-tax models

marginal_effects_list <- lapply(
  selected_models$model_id,
  function(current_model_id) {
    
    # Retrieve the row containing information about the current model.
    current_model_info <- selected_models[
      selected_models$model_id == current_model_id,
      ,
      drop = FALSE
    ]
    
    
    # Retrieve the corresponding fitted plm model.
    current_model <- results2$models[[current_model_id]]
    
    
    # Retrieve the corresponding clustered covariance matrix.
    current_vcov <- results2$vcov_matrices[[current_model_id]]
    
    
    # Calculate the marginal effects over the full support.
    calculate_marginal_effects(
      model = current_model,
      V_cluster = current_vcov,
      model_information = current_model_info,
      moderator = "l_serv_managr",
      grid_points = 200,
      confidence_level = 0.9
    )
  }
)

# combine
marginal_effects_data <- do.call(
  rbind,
  marginal_effects_list
)

rownames(marginal_effects_data) <- NULL

# Check that the marginal effect is zero at the turning point
turning_point_check <- unique(
  marginal_effects_data[
    ,
    c(
      "model_id",
      "panel_label",
      "beta_inequality",
      "beta_interaction",
      "turning_point_log",
      "turning_point_level"
    )
  ]
)

turning_point_check$marginal_effect_at_turning_point <-
  turning_point_check$beta_inequality +
  turning_point_check$beta_interaction *
  turning_point_check$turning_point_log

turning_point_check

unique(
  marginal_effects_data[
    ,
    c(
      "panel_label",
      "beta_inequality",
      "beta_interaction",
      "turning_point_level"
    )
  ]
)

# prepare a dataset with one turning-point observation per panel:
turning_point_data <- unique(
  marginal_effects_data[
    ,
    c(
      "model_id",
      "panel_label",
      "turning_point_level",
      "support_min_level",
      "support_max_level"
    )
  ]
)

# only retain a turning point when it lies within the observed support:
turning_point_data <- subset(
  turning_point_data,
  is.finite(turning_point_level) &
    turning_point_level >= support_min_level &
    turning_point_level <= support_max_level
)

turning_point_data$turning_point_label <- paste0(
  "Turning point = ",
  round(turning_point_data$turning_point_level, 2)
)

#figure 
figure_3_new <- ggplot(
  data = marginal_effects_data,
  mapping = aes(
    x = tertiarization_level,
    y = marginal_effect
  )
) +
  
  # Pointwise 90% confidence interval.
  geom_ribbon(
    mapping = aes(
      ymin = conf_low,
      ymax = conf_high
    ),
    fill = "#3399FF",
    alpha = 0.22,
    colour = NA
  ) +
  
  # Estimated marginal effect.
  geom_line(
    colour = "#002BFF",
    linewidth = 1.85
  ) +
  
  # Zero marginal effect.
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "black",
    linewidth = 0.9
  ) +
  
  # Point-estimate turning point.
  #
  # This identifies the point where the estimated marginal effect crosses
  # zero. It does not necessarily identify the point where the effect becomes
  # statistically significant.
  geom_vline(
    data = turning_point_data,
    mapping = aes(
      xintercept = turning_point_level
    ),
    inherit.aes = FALSE,
    linetype = "dotted",
    colour = "#D55E00",
    linewidth = 0.9
  ) +
  
  # Two panels: pre-tax and post-tax Gini.
  facet_wrap(
    facets = vars(panel_label),
    nrow = 1,
    scales = "fixed"
  ) +
  # facet_wrap(
  #   vars(panel_label),
  #   labeller = as_labeller(
  #     c(
  #       "l_p90p100" = "A. Top 10% (pre-tax)",
  #       "l_p90p100d" = "B. Top 10% (post-tax)"
  #     )
  #   )
  # )+
  
  
  labs(
    x = "Tertiarization",
    y = expression(
      "Marginal effect (%) of the Gini index \n      (ln) on CO2 (ln) emissions" #Gini, income share of \nthe Top 10% (ln)
    )
  ) +
  
  geom_text(
    data = turning_point_data,
    aes(
      x = turning_point_level,
      y = Inf,
      label = "Turning point"#turning_point_label
    ),
    inherit.aes = FALSE,
    angle = 90,
    vjust = 1.4,
    hjust = 1.1,
    colour = "#D55E00",
    size = 9
  ) +
  
  theme_classic(
    base_size = 25
  ) +
  
  scale_x_continuous(
    breaks = c(0, 5, 10, 15, 20, 25)
  ) +
  
  theme(
    strip.background = element_blank(),
    
    strip.text = element_text(
      size = 28,
      face = "bold",
      margin = margin(
        b = 10
      )
    ),
    
    axis.title = element_text(
      size = 29
    ),
    
    axis.title.x = element_text(
      margin = margin(
        t = 16
      )
    ),
    
    axis.title.y = element_text(
      margin = margin(
        r = 16
      )
    ),
    
    axis.text = element_text(
      size = 24
    ),
    
    panel.spacing = grid::unit(
      1.5,
      "cm"
    ),
    
    plot.margin = margin(
      t = 15,
      r = 20,
      b = 15,
      l = 35
    )
  )

figure_3_new
# ggsave("C:/Users/perniol/Desktop/code+data DEF/marginal_effects.pdf", width = 20, height = 10, units = "in")

#### FIGURE 5 ####

graph <- countries_above_tp_per_year(input_data = data, tp_variable = "serv_managr", tp = 2.99)

fig5b <- ggplot(aes(x = an, y = pa), data = graph) + geom_line(linewidth = 1.5) + labs (x = "Year", y = "Number of countries over 2.96") + theme_classic(base_size = 15)
#ggsave("Desktop/Progetto/Sprint finale/Grafici DEF/quanti_sp.pdf", width = 15, height = 10, units = "in")


graph2 <- countries_above_tp_per_year(input_data = data, tp_variable = "serv_managr", tp = 2.96)
graph2 <- merge(graph,graph2,by=c("an") )

graph2_long <- graph2 %>%
  pivot_longer(cols = c(pa.x, pa.y), names_to = "group", values_to = "pa")

fig5b <- ggplot(graph2_long, aes(x = an, y = pa, color = group)) +
  geom_line(linewidth = 1.5) + # default è 0.5
  labs(x = "Year", y = "Number of tertiarized countries") +
  theme_classic(base_size = 16) +
  scale_color_manual(name = "", values = c("pa.x" = "blue", "pa.y" = "red"),
                     labels = c("pa.x" = "Pre-tax", "pa.y" = "Post-tax"))


geo <- frequancy_country_above_tp(input_data = data, tp = 2.99, tp_variable = "serv_managr")  
geo <- geo[geo$ccc > 0.5,] # geo$ccc > 0.5   geo$ccc != 0
geo$ccc <- geo$ccc*100

#mappa 
world <- ne_countries(scale = "medium", returnclass = "sf")
world_data <- left_join(world, geo, by=c("iso_a3"="pa"))
fig5a <- world_data%>% 
  ggplot() +
  theme_light() +
  geom_sf(aes(fill = ccc)) +  
  scale_fill_distiller(palette = "Oranges", direction = 1, na.value = "white", guide = guide_colourbar(direction = "horizontal")) + #PuBuGn PRGn Oranges YlOrRd
  ggtitle("% of years in which the country is \n above tertiarization threshold") +
  #guides(fill = guide_legend(nrow = 1)) +
  theme(
    plot.title = element_text(hjust = 0.5, size=18),
    legend.title = element_blank(),
    legend.position = c(0.15,0.25),   #c(0.15,0.25) / "none"
    legend.key.size = unit(0.7, "cm"),
    legend.background = element_rect(fill="transparent"),
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  ) 

# Figure 5 toghether 
cowplot::plot_grid(
  fig5a,
  fig5b,
  nrow = 1, labels = c("A","B"), rel_widths = c(0.5,0.5))
#pdf -> 4.5 x 12

#### mechanism regressions (Table 6, A20,A21,A22) ####

list_paesi_over <-contries_splitted(over_in_at_least = 0.70, input_data = data, tp_variable = "serv_managr", tp = 2.99, output = "list")
data$over <- ifelse(data$Country.Code %in% list_paesi_over,1,0)

# c("gini_wid", "p90p100")
# c("gini_wid_disp", "p90p100d")
ineq_vars_list <- c("gini_wid", "p90p100")

y_list <- c("man_empl","inv_cf_p99p100")

table6 <-run_mechanism_regressions(data = data, y_list, ineq_vars_list, interaction_var = "over" , output_type = "text", vcov = F)

#### Figure 6,7 (default),A.2,A.3 ####

#divide the dataset based on countries below/above the tp in at least over_in_at_least % of the years
paesi_over <- contries_splitted(over_in_at_least = 0.7, input_data = data, tp_variable = "serv_managr", tp = 2.96)$paesi_over
paesi_under <- contries_splitted(over_in_at_least = 0.7, input_data = data, tp_variable = "serv_managr", tp = 2.96)$paesi_under

#remove the outliers (gini >)85 
paesi_over <- paesi_over %>% filter(gini_wid_disp < 85)

#scatterplot man_empl
empl_l <- ggplot(paesi_under, aes(gini_wid_disp, man_empl))  # gini_wid - gini_wid_disp
empl_l <- empl_l + geom_point(alpha = 1, size = 0.8) + geom_smooth(method = lm,se = T) + labs (x = "Gini (post-tax)", y = "Employment in industry") + theme_classic(base_size = 40)  

empl_m <- ggplot(paesi_over, aes(gini_wid_disp, man_empl)) #gini_wid - gini_wid_disp
empl_m <- empl_m + geom_point(alpha = 1, size = 0.8) + geom_smooth(method = lm,se = T) + labs (x = "Gini (post-tax)", y = "Employment in industry") + theme_classic(base_size = 40)


cowplot::plot_grid(
  empl_l,
  empl_m,
  nrow = 1, labels = c("A","B"), rel_widths = c(0.5,0.5), label_size = 40,  label_x = 0.9, label_y = 0.95)
#pdf -> 10 x 30

#### FIGURE 8 ####
Elasticity <- c(0.492, 0.492 + 1.673, 0.339, 0.339 + 1.014,
                0.642, 0.642 + 3.230, 0.456, 0.456 + 1.483)

group <- c("Below tertiarization threshold", "Above tertiarization threshold  ",
           "Below tertiarization threshold", "Above tertiarization threshold  ",
           "Below tertiarization threshold", "Above tertiarization threshold  ",
           "Below tertiarization threshold", "Above tertiarization threshold  ")

category <- c("Consumption \n(pre-tax)", "Consumption \n(pre-tax)",
              "Consumption \n(post-tax)", "Consumption \n(post-tax)",
              "Investment \n(pre-tax)", "Investment \n(pre-tax)",
              "Investment \n(post-tax)", "Investment \n(post-tax)")

confidence_interval <- c(0.104, 0.4377774, 0.092, 0.2841816,
                         0.132, 1.060512, 0.143, 0.8940524)

df <- c(2572, 2572, 2572, 2572, 2558, 2558, 2558, 2558)

data_ef <- data.frame(group, category, Elasticity, confidence_interval, df)

data_ef$lower <- data_ef$Elasticity - qt(0.95, data_ef$df) * data_ef$confidence_interval
data_ef$upper <- data_ef$Elasticity + qt(0.95, data_ef$df) * data_ef$confidence_interval

# Keep only the 4 investment entries and reset x-axis levels
data_ef_plot <- data_ef %>%
  dplyr::filter(!grepl("Consumption", category)) %>%
  dplyr::mutate(
    category = factor(
      category,
      levels = c("Investment \n(pre-tax)", "Investment \n(post-tax)")
    )
  )

data_ef_plot %>%
  ggplot(aes(x = category, y = Elasticity, color = as.factor(group))) +
  geom_point(size = 6) +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
  theme_light() +
  theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    axis.title.y = element_text(size = 22),
    legend.text = element_text(size = 22),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box = "horizontal"
  ) +
  scale_color_manual(values = c("blue", "red")) +
  scale_x_discrete(drop = TRUE)
