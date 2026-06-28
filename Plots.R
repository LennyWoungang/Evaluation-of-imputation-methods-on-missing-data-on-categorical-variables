library(ggplot2)
library(dplyr)
library(tidyr)
library(ggridges)
library(viridis)
library(scales)


#1. Bias + error bars 
table5_eval  <- table5_eval  %>% mutate(missingness = "5%",  method = "CCA")
table15_eval <- table15_eval %>% mutate(missingness = "15%", method = "CCA")
table30_eval <- table30_eval %>% mutate(missingness = "30%", method = "CCA")

eval_df <- bind_rows(table5_eval, table15_eval, table30_eval)
eval_df$missingness <- factor(eval_df$missingness,
                              levels = c("5%", "15%", "30%"))


ggplot(eval_df, aes(x = parameter, y = bias)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  
  # point for mean bias
  geom_point(aes(color = "Mean Bias"), size = 2) +
  
  # error bars for empirical SD variability
  geom_errorbar(aes(ymin = bias - emp_sd, 
                    ymax = bias + emp_sd, 
                    color = "Bias ± Empirical SD"),
                width = 0.2) +
  
  facet_grid(method ~ missingness) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.title = element_blank()
  ) +
  labs(
    title = "Bias by Method and Missingness",
    y = "Bias",
    x = "Parameter"
  ) +
  
  # choose colors for legend entries
  scale_color_manual(values = c(
    "Mean Bias" = "black",
    "Bias ± Empirical SD" = "blue"
  ))

### How close is the average estimate to the truth, and how much does it fluctuate
### Bias is literally the distance from the true parameter.
### Plotting bias on the y-axis makes the truth the reference point.
### Together, bias + SD tells you if the method consistently returns values close to the truth.



#3. Maybe compare them by percentage and then also by individual method, so two sets of plots


#4. Coverage heat map

ggplot(eval_df, aes(x = missingness, y = parameter, fill = coverage)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(limits = c(0,1)) +
  theme_bw() +
  labs(
    title = "Coverage Heatmap (CCA)",
    x = "Missingness Level",
    y = "Parameter",
    fill = "Coverage"
  )

#5a.Cloud graph around 
library(tidyverse)

make_long <- function(res, missing_label, true_params) {
  as.data.frame(res$betas) %>%
    mutate(replication = row_number()) %>%
    pivot_longer(cols = -replication, names_to = "parameter", values_to = "estimate") %>%
    mutate(
      true_value = true_params[parameter],
      missingness = missing_label
    )
}

df5  <- make_long(res5,  "5%",  true_params)
df15 <- make_long(res15, "15%", true_params)
df30 <- make_long(res30, "30%", true_params)

plot_df <- bind_rows(df5, df15, df30)

ggplot(plot_df, aes(x = replication, y = estimate)) +
  geom_point(alpha = 0.4, size = 1, color = "steelblue") +
  geom_hline(aes(yintercept = true_value), linetype = "dashed", color = "red") +
  facet_grid(missingness ~ parameter, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 8),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Estimates Across Replications by Parameter and Missingness",
    x = "Replication",
    y = "Estimate"
  )


#5b.Cloud graph around  (individual)
geom_true_line <- function() {
  list(
    geom_hline(aes(yintercept = true_value), 
               linetype = "dashed", color = "red"),
    geom_text(
      aes(x = -Inf, y = true_value, 
          label = paste0("True = ", round(true_value, 3))),
      color = "red", hjust = -0.1, vjust = -0.5, size = 3
    )
  )
}

p5 <- ggplot(df5, aes(x = replication, y = estimate)) +
  geom_point(alpha = 0.4, size = 1, color = "steelblue") +
  geom_true_line() +
  facet_wrap(~ parameter, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Estimates Across Replications – 5% Missingness",
    x = "Replication",
    y = "Estimate"
  )
p5


p15 <- ggplot(df15, aes(x = replication, y = estimate)) +
  geom_point(alpha = 0.4, size = 1, color = "steelblue") +
  geom_true_line() +
  facet_wrap(~ parameter, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Estimates Across Replications – 15% Missingness",
    x = "Replication",
    y = "Estimate"
  )
p15


p30 <- ggplot(df30, aes(x = replication, y = estimate)) +
  geom_point(alpha = 0.4, size = 1, color = "steelblue") +
  geom_true_line() +
  facet_wrap(~ parameter, scales = "free_y") +
  theme_bw() +
  theme(
    strip.text = element_text(size = 9),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Estimates Across Replications – 30% Missingness",
    x = "Replication",
    y = "Estimate"
  )
p30








