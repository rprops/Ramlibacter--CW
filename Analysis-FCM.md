---
title: "Analysis-FCM"
author: "Ruben Props"
date: "31 January, 2018"
output:
  html_document:
    code_folding: show
    highlight: haddock
    keep_md: yes
    theme: united
    toc: yes
    toc_float:
      collapsed: no
      smooth_scroll: yes
      toc_depth: 2
    css: report_styles.css
    df_print: paged
editor_options: 
  chunk_output_type: console
---





# Cell density 




```r
p_counts <- ggplot(counts, aes(x = ExactTime, y = Total.cells, fill = NutrientCondition))+
  geom_line(aes(color = NutrientCondition))+
  geom_point(shape = 21, size = 4)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.direction = "horizontal",legend.position = "bottom")+
  ylab("Cell density (cells/µL)")+
  xlab("Time (h)")+
  labs(title="Total population")+
  guides(color = FALSE)

p_counts_log <- ggplot(counts, aes(x = ExactTime, y = Total.cells, fill = NutrientCondition))+
  geom_line(aes(color = NutrientCondition))+
  geom_point(shape = 21, size = 4)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.direction = "horizontal",legend.position = "bottom")+
  ylab(expression("Cell density - cells µL"^"-1"))+
  xlab("Time (h)")+
  labs(title="Total population")+
  guides(color = FALSE)+
  scale_y_continuous(trans='log2', breaks = c(10, 100,1000,5000,1e4, 1.5e4), limits = c(1,1.6e4))

p_HNA <- ggplot(counts, aes(x = ExactTime, y = HNA.cells, fill = NutrientCondition))+
  geom_line(aes(color = NutrientCondition))+
  geom_point(shape = 21, size = 4)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.direction = "horizontal",legend.position = "bottom")+
  ylab(expression("Cell density - cells µL"^"-1"))+
  xlab("Time (h)")+
  labs(title="HNA population")+
  guides(color = FALSE)


p_LNA <- ggplot(counts, aes(x = ExactTime, y = LNA.cells, fill = NutrientCondition))+
    geom_line(aes(color = NutrientCondition))+
  geom_point(shape = 21, size = 4)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.direction = "horizontal",legend.position = "bottom")+
  ylab(expression("Cell density - cells µL"^"-1"))+
  xlab("Time (h)")+
  labs(title="LNA population")+
  guides(color = FALSE)

grid_arrange_shared_legend(p_counts, p_HNA, p_LNA, ncol = 3)
```

<img src="Figures-FCM/cached/fcm-density-data-2-1.png" style="display: block; margin: auto;" />

```r
p_HNA_pct <- ggplot(counts, aes(x = ExactTime, y = 100*pct_HNA.cells, fill = NutrientCondition))+
  geom_line(aes(color = NutrientCondition))+
  geom_point(shape = 21, size = 4)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=14),
        legend.direction = "horizontal",legend.position = "bottom")+
  ylab("%HNA cells")+
  xlab("Time (h)")+
  labs(title="HNA population")+
  guides(color = FALSE)
```


```r
plot_grid(p_counts_log, p_HNA_pct,nrow = 2, align = 'v')
```

<img src="Figures-FCM/cached/fcm-density-data-3-1.png" style="display: block; margin: auto;" />

```r
plot_grid(p_HNA, p_HNA_pct, nrow = 2, align = 'v')
```

<img src="Figures-FCM/cached/fcm-density-data-3-2.png" style="display: block; margin: auto;" />

# Diversity dynamics


```r
# Resample to atleast 10,000 cells for diversity assay
# flowData_transformed_subs <- FCS_resample(flowData_transformed, sample = 10000, replace = TRUE)

# Calculate phenotypic alpha diversity
diversity_fcm <- Diversity_rf(flowData_transformed, param = param, cleanFCS = FALSE, parallel = TRUE,
                              ncores = 10)
```

```
## -------------------------------------------------------------------------------------------------
## Wed Jan 31 15:59:10 2018 --- Normalizing your FCS data based on maximum FL1-H value
## --- Maximum FL1-H before normalizing: 14.95
## --- Maximum FL3-H before normalizing: 13.27
## --- Maximum SSC-H before normalizing: 17.33
## --- Maximum FSC-H before normalizing: 16.96
## -------------------------------------------------------------------------------------------------
## --- Maximum FL1-H after normalizing: 1
## --- Maximum FL3-H after normalizing: 0.89
## --- Maximum SSC-H after normalizing: 1.16
## --- Maximum FSC-H after normalizing: 1.13
## -------------------------------------------------------------------------------------------------
##  
## Wed Jan 31 15:59:49 2018 --- Using 10 cores for calculations
## Wed Jan 31 16:24:51 2018 --- Closing workers
## Wed Jan 31 16:24:51 2018 --- Alpha diversity metrics (D0,D1,D2) have been computed after 100 bootstraps
## -----------------------------------------------------------------------------------------------------
## 
```



```r
# Add metadata to phenotypic diversity estimate
diversity_fcm <- dplyr::left_join(diversity_fcm, counts, by = c("Sample_names" =  "Samples"))
diversity_fcm <- diversity_fcm %>% dplyr::filter(Timepoint > 5) # Only consider after first 5 samples due to bleaching of tubing

# Plot results
p_div <- ggplot(diversity_fcm, aes(x = ExactTime, y = D2, fill = NutrientCondition))+
    geom_line(aes(color = NutrientCondition))+
    geom_point(shape = 21, size = 4)+
    theme_bw()+
    scale_fill_brewer("Nutrient condition", palette = "Accent")+
    scale_color_brewer("Nutrient condition", palette = "Accent")+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom",
        strip.text = element_text(size = 16))+
    ylab(expression("Phenotypic diversity - D"[2]))+
    facet_grid(~NutrientCondition)+
    xlab("Time (h)")+
    guides(fill = FALSE, color = FALSE)+
  geom_ribbon(aes(ymin = D2 - sd.D2, ymax = D2 + sd.D2), alpha = 0.3)

print(p_div)
```

<img src="Figures-FCM/cached/diversity-fcm-2-1.png" style="display: block; margin: auto;" />


```r
# Reshape combined data
diversity_fcm_long <- tidyr::gather(diversity_fcm, population, Density, Total.cells, LNA.cells, HNA.cells, D2)
diversity_fcm_long$population <- plyr::revalue(diversity_fcm_long$population,  c("Total.cells"="Whole population", "HNA.cells"="HNA population", "LNA.cells"="LNA population", "D2" = "Phenotypic diversity"))

# Combine diversity and count plot
p_count2 <- ggplot(counts, aes(x = ExactTime, y = Total.cells, fill = NutrientCondition))+
  geom_line(aes(color = NutrientCondition))+
  geom_point(shape = 21, size = 4)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom",
        strip.text = element_text(size = 16))+
  ylab("Cell density (cells/µL)")+
  xlab("Time (h)")+
  facet_grid(~NutrientCondition)+
  guides(color = FALSE, fill = FALSE)

p_HNA2 <- ggplot(counts, aes(x = ExactTime, y = HNA.cells, fill = NutrientCondition))+
  geom_line(aes(color = NutrientCondition))+
  geom_point(shape = 21, size = 4)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom",
        strip.text = element_text(size = 16))+
  ylab("Cell density (cells/µL)")+
  xlab("Time (h)")+
  facet_grid(~NutrientCondition)+
  guides(color = FALSE, fill = FALSE)

p_LNA2 <- ggplot(counts, aes(x = ExactTime, y = LNA.cells, fill = NutrientCondition))+
  geom_line(aes(color = NutrientCondition))+
  geom_point(shape = 21, size = 4)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  scale_color_brewer(palette = "Accent")+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=20),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom",
        strip.text = element_text(size = 16))+
  ylab("Cell density (cells/µL)")+
  xlab("Time (h)")+
  facet_grid(~NutrientCondition)+
  guides(color = FALSE, fill = FALSE)

grid.arrange(p_count2, p_div, nrow = 2)
```

<img src="Figures-FCM/cached/count-diversity-fcm-1.png" style="display: block; margin: auto;" />

```r
grid.arrange(p_HNA2, p_div, nrow = 2)
```

<img src="Figures-FCM/cached/count-diversity-fcm-2.png" style="display: block; margin: auto;" />

```r
grid.arrange(p_LNA2, p_div, nrow = 2)
```

<img src="Figures-FCM/cached/count-diversity-fcm-3.png" style="display: block; margin: auto;" />

# Growth rate estimation

## Grofit


```r
# Format data for grofit
counts <- counts %>% group_by(NutrientCondition) %>% 
  dplyr::mutate(ratio_TotalCells = log(Total.cells/min(Total.cells)))

## Samples as rows and timepoints as columns
counts_tmp <- counts %>% dplyr::filter(Timepoint > 5) # Only consider after first 5 samples due to bleaching of tubing
count_wide <- spread(counts_tmp[, c(6,2,8)], ExactTime, Total.cells)
count_wide2 <- spread(counts_tmp[, c(6,9,8)], ExactTime, ratio_TotalCells)
time_wide <- as.matrix(spread(counts_tmp[, c(6,8)], ExactTime, ExactTime)[,-1])

## Add two additional columns (mandatory for grofit to run properly)
count_wide <- data.frame(ReactorName = c("R1", "R2","R3"), 
                         NutrientLevel = c("Low", "Medium", "High"), count_wide)
count_wide2 <- data.frame(ReactorName = c("R1", "R2","R3"), 
                         NutrientLevel = c("Low", "Medium", "High"), count_wide2)

## Fit grofit models
Fit <- grofit::gcFit(time_wide, count_wide, grofit.control(nboot.gc = 1000, smooth.gc=0.5, interactive = FALSE))
```

```
## 
## 
## = 1. growth curve =================================
## ----------------------------------------------------
## --> Try to fit model logistic--> Try to fit model richards--> Try to fit model gompertz--> Try to fit model gompertz.exp
## 
## 
## = 2. growth curve =================================
## ----------------------------------------------------
## --> Try to fit model logistic--> Try to fit model richards--> Try to fit model gompertz--> Try to fit model gompertz.exp
## 
## 
## = 3. growth curve =================================
## ----------------------------------------------------
## --> Try to fit model logistic--> Try to fit model richards--> Try to fit model gompertz--> Try to fit model gompertz.exp
```

```r
Fit2 <- grofit::gcFit(time_wide, count_wide2, grofit.control(nboot.gc = 1000, smooth.gc=0.5, interactive = FALSE))
```

```
## 
## 
## = 1. growth curve =================================
## ----------------------------------------------------
## --> Try to fit model logistic--> Try to fit model richards--> Try to fit model gompertz--> Try to fit model gompertz.exp
## 
## 
## = 2. growth curve =================================
## ----------------------------------------------------
## --> Try to fit model logistic--> Try to fit model richards--> Try to fit model gompertz--> Try to fit model gompertz.exp
## 
## 
## = 3. growth curve =================================
## ----------------------------------------------------
## --> Try to fit model logistic--> Try to fit model richards--> Try to fit model gompertz--> Try to fit model gompertz.exp
```

```r
sum <- summary(Fit)
par(mfrow=c(3,1))
plot(Fit, opt = "s", slope = TRUE, colSpline = 4, cex = 2, title = "Spline fit") # show spline fits
par(mfrow=c(1,1))

# We continue to work with the spline fits as they seem to be the best
growth_results <- summary(Fit)[, c("concentration","mu.bt", "lambda.bt", "A.bt",
                          "ci95.mu.bt.up", "ci95.mu.bt.lo",
                          "ci95.lambda.bt.up", "ci95.lambda.bt.lo",
                          "ci95.A.bt.up", "ci95.A.bt.lo")]
```


```r
# Plot lag-phase length, maximum growth rate, etc etc.
growth_results <- growth_results %>% dplyr::mutate(concentration_value = as.numeric(gsub("mg/L R2A", "", concentration)))

# Lag phase
p_lambda <- ggplot(growth_results, aes(x = concentration_value, y = lambda.bt, fill = concentration))+
  geom_point(shape = 21, size = 5)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom",
        strip.text = element_text(size = 16),
        axis.text.x = element_text(size = 16))+
  ylab(expression(lambda ~ "(h)"))+
  xlab(expression("Medium concentration - mg L"^"-1"))+
  guides(color = FALSE, fill = FALSE)+
  geom_errorbar(aes(ymax = ci95.lambda.bt.up, ymin = ci95.lambda.bt.lo), width = 0.025)+
  scale_y_continuous(trans='log10', breaks = seq(0, 80, 20), minor_breaks = seq(0, 80, 10), limits = c(NA,80))+
  scale_x_continuous(trans='log10', breaks = c(1, 10, 100), minor_breaks = c(1,seq(10, 100, 10)), limits = c(NA, 100))+
  ggtitle("C.")

# print(p_lambda)

# maximum growth rate
p_mu <- ggplot(growth_results, aes(x = concentration_value, y = mu.bt, fill = concentration))+
  geom_point(shape = 21, size = 5)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=16),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom",
        strip.text = element_text(size = 16),
        axis.text.x = element_text(size = 16))+
  ylab(expression(mu["max"] ~ "- cells" ~ "mL"^"-1"~"h"^"-1"))+
  xlab(expression("Medium concentration - mg L"^"-1"))+
  guides(color = FALSE, fill = FALSE)+
  geom_errorbar(aes(ymax = ci95.mu.bt.up, ymin = ci95.mu.bt.lo), width = 0.025)+
  scale_y_continuous(trans='log10', breaks = seq(0,1800, 250), minor_breaks = seq(0,1800, 125), limits = c(250,1800))+
  scale_x_continuous(trans='log10', breaks = c(1, 10, 100), minor_breaks = c(1,seq(10, 100, 10)), limits = c(NA, 100))+
  ggtitle("A.")

# print(p_mu)

# Carrying capacity
p_A <- ggplot(growth_results, aes(x = concentration_value, y = A.bt, fill = concentration))+
  geom_point(shape = 21, size = 5)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
  theme(axis.text=element_text(size=16), axis.title=element_text(size=16),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom",
        strip.text = element_text(size = 16),
        axis.text.x = element_text(size = 16))+
  ylab(expression("Carrying capacity" ~ "- cells µL"^"-1"))+
  xlab(expression("Medium concentration - mg L"^"-1"))+
  guides(color = FALSE, fill = FALSE)+
  geom_errorbar(aes(ymax = ci95.A.bt.up, ymin = ci95.A.bt.lo), width = 0.025)+
  scale_y_continuous(trans='log10', breaks = seq(0,16000, 2000),minor_breaks = seq(0,16000, 1000), limits = c(NA, 16000)) +
  scale_x_continuous(trans='log10', breaks = c(1, 10, 100),minor_breaks = c(1,seq(10, 100, 10)), limits = c(NA, 100))+
  ggtitle("D.")

# print(p_A)
```

## Growthcurver


```r
# Create wide format
count_gc_wide <- spread(counts[, c(6,2,8)], NutrientCondition, Total.cells)

# Run growthcurver
gc_fit <- list()
gc_fit[[1]] <- SummarizeGrowth(count_gc_wide$ExactTime[!is.na(count_gc_wide$`1 mg/L R2A`)], count_gc_wide$`1 mg/L R2A`[!is.na(count_gc_wide$`1 mg/L R2A`)])
gc_fit[[2]] <- SummarizeGrowth(count_gc_wide$ExactTime[!is.na(count_gc_wide$`10 mg/L R2A`)], count_gc_wide$`10 mg/L R2A`[!is.na(count_gc_wide$`10 mg/L R2A`)])
gc_fit[[3]] <- SummarizeGrowth(count_gc_wide$ExactTime[!is.na(count_gc_wide$`100 mg/L R2A`)], count_gc_wide$`100 mg/L R2A`[!is.na(count_gc_wide$`100 mg/L R2A`)])

gc_results <- data.frame(concentration = colnames(count_gc_wide)[-1], 
                         rbind(gc_fit[[1]]$vals, gc_fit[[2]]$vals, gc_fit[[3]]$vals)) 

# Upper error on tgen
round(log(2)/(do.call(rbind, gc_results$r)+do.call(rbind, gc_results$r_se)),2)
```

```
##      [,1]
## [1,] 1.73
## [2,] 2.46
## [3,] 1.88
```

```r
# Lower error on tgen
round(log(2)/(do.call(rbind, gc_results$r)-do.call(rbind, gc_results$r_se)),2)
```

```
##      [,1]
## [1,] 1.85
## [2,] 2.62
## [3,] 2.13
```

```r
# Store in grofit dataframe
growth_results$mu_spec <- do.call(rbind, gc_results$r)
growth_results$mu_se <- do.call(rbind, gc_results$r_se)

# specific growth rate
p_mu_spec <- ggplot(growth_results, aes(x = concentration_value, y = mu_spec, fill = concentration))+
  geom_point(shape = 21, size = 5)+
  theme_bw()+
  scale_fill_brewer("Nutrient condition", palette = "Accent")+
    theme(axis.text=element_text(size=16), axis.title=element_text(size=16),
        title=element_text(size=20), legend.text=element_text(size=16),
        legend.direction = "horizontal",legend.position = "bottom",
        strip.text = element_text(size = 16),
        axis.text.x = element_text(size = 16))+
  ylab(expression(mu["spec"] ~ "- h"^"-1"))+
  xlab(expression("Medium concentration - mg L"^"-1"))+
  guides(color = FALSE, fill = FALSE)+
  geom_errorbar(aes(ymax = mu_spec+mu_se, ymin = mu_spec-mu_se), width = 0.025)+
  scale_y_continuous(trans='log10', breaks = seq(0, 0.45, 0.1), minor_breaks = seq(0, 0.45, 0.05), limits = c(0.2,0.45))+
  scale_x_continuous(trans='log10', breaks = c(1, 10, 100), minor_breaks = c(1,seq(10, 100, 10)), limits = c(NA, 100))+
  ggtitle("B.")

# Plot all 4 parameters
plot_grid(p_mu, p_mu_spec, 
             p_lambda, p_A, ncol = 2, align = 'hv')
```

<img src="Figures-FCM/cached/growthcurver-analysis-1.png" style="display: block; margin: auto;" />


