# ============================================================
# Human FSI extended PNN analysis using WFA and IMARIS
# Author: Florian Wildner
# Date: 10.03.2026
# ============================================================


# =========================
# --- Load libraries ---
# =========================
library(tidyverse)
library(dplyr)
library(janitor)
library(ggplot2)

library(rstatix)
library(viridis)
library(scico)
library(gridExtra)

library(lme4)
library(lmerTest)
library(emmeans)

library(ggsignif)
library(ggpubr)
library(tibble)
library(here)


# =========================
# --- Output directory ---
# =========================
# Save all figures to: output/IMARIS/<current_date>/
date_str <- format(Sys.Date(), "%y%m%d")
out_dir  <- file.path("output", "IMARIS", date_str)

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

options(fig_out_dir = out_dir)

# Helper function to save figures
save_fig <- function(plot,
                     name,
                     w = 3/2.54,
                     h = 3/2.54,
                     device = "pdf") {
  
  out_dir <- getOption("fig_out_dir")
  
  ggsave(
    filename = file.path(out_dir, paste0(name, ".", device)),
    plot     = plot,
    width    = w,
    height   = h,
    units    = "in",
    device   = device
  )
}


# =========================
# --- Load IMARIS data ---
# =========================
# Each CSV corresponds to one patient (pseudonymized)

df_pat1 <- read_csv(file.path("data", "WFA_HUMAN_QUANT_IMARIS_PatA.csv")) %>%
  select(where(~ any(!is.na(.)))) %>% 
  mutate(pat = "pat_a")

df_pat2 <- read_csv(file.path("data", "WFA_HUMAN_QUANT_IMARIS_PatB.csv")) %>%
  select(where(~ any(!is.na(.)))) %>% 
  mutate(pat = "pat_b")

df_pat3 <- read_csv(file.path("data", "WFA_HUMAN_QUANT_IMARIS_PatC.csv")) %>%
  select(where(~ any(!is.na(.)))) %>% 
  mutate(pat = "pat_c")


# Combine all patients
df_imaris <- bind_rows(df_pat1, df_pat2, df_pat3)


# =========================
# --- Column name mapping ---
# =========================
# Create lookup table for original vs cleaned names (for plotting)

name_dictionary <- tibble(
  original_name = colnames(df_imaris),
  clean_name    = make_clean_names(colnames(df_imaris))
) %>%
  mutate(
    clean_name = ifelse(clean_name == "net_volume_mm3",
                        "net_volume_um3",
                        clean_name),
    original_name = ifelse(clean_name == "pat",
                           "Patient",
                           original_name)
  )

label_lookup <- setNames(
  name_dictionary$original_name,
  name_dictionary$clean_name
)


# =========================
# --- Data preprocessing ---
# =========================
df_imaris <- df_imaris %>%
  clean_names() %>%
  mutate(
    layer = as.factor(layer),
    pat   = as.factor(pat),
    pv    = as.factor(pv)
  ) %>%
  rename(net_volume_um3 = net_volume_mm3) %>%
  select(pat, layer, everything())


# ============================================================
# --- Analysis + plotting function (core pipeline) ---
# ============================================================
analyze_parameter <- function(data,
                              response_var,
                              normal = TRUE,
                              save_plot = TRUE,
                              w = 3/2.54,
                              h = 3/2.54) {
  
  # -------------------------
  # Prepare data
  # -------------------------
  data <- data %>%
    select(pat, layer, all_of(response_var)) %>%
    filter(!is.na(layer))
  
  # -------------------------
  # Normality testing
  # -------------------------
  normality <- data %>%
    group_by(layer) %>%
    summarise(
      shapiro_p = shapiro.test(.data[[response_var]])$p.value,
      .groups = "drop"
    )
  
  cat("\n===== Shapiro-Wilk Normality Test =====\n")
  print(normality)
  
  log_transform <- !normal && any(normality$shapiro_p < 0.05)
  
  if (log_transform) {
    cat("\nLog-transforming data.\n")
    data[[response_var]] <- log10(data[[response_var]] + 1)
  }
  
  # -------------------------
  # Mixed-effects model
  # -------------------------
  formula <- as.formula(paste(response_var, "~ layer + (1 | pat)"))
  model   <- lmer(formula, data = data)
  
  cat("\n===== ANOVA =====\n")
  print(anova(model, type = 3))
  
  cat("\n===== Fixed Effects =====\n")
  print(summary(model)$coefficients)
  
  # -------------------------
  # Random effects + ICC
  # -------------------------
  re_var <- as.data.frame(VarCorr(model))
  
  var_pat <- re_var$vcov[re_var$grp == "pat"]
  var_res <- re_var$vcov[re_var$grp == "Residual"]
  icc     <- var_pat / (var_pat + var_res)
  
  cat("\nICC:", round(icc, 4), "\n")
  
  # -------------------------
  # Post hoc tests
  # -------------------------
  emm      <- emmeans(model, ~ layer)
  pairwise <- pairs(emm, adjust = "tukey")
  
  # -------------------------
  # Summary stats for plotting
  # -------------------------
  summary_data <- data %>%
    group_by(layer) %>%
    summarise(
      mean = if (log_transform) median(.data[[response_var]], na.rm = TRUE)
      else mean(.data[[response_var]], na.rm = TRUE),
      sem  = sd(.data[[response_var]], na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(x_num = as.numeric(layer))
  
  # -------------------------
  # Significance annotation
  # -------------------------
  y_max   <- max(data[[response_var]], na.rm = TRUE)
  y_range <- diff(range(data[[response_var]], na.rm = TRUE))
  
  sig_labels <- as.data.frame(pairwise) %>%
    mutate(
      layer1 = sub(" - .*", "", contrast),
      layer2 = sub(".* - ", "", contrast),
      label = case_when(
        p.value < 0.0001 ~ "***",
        p.value < 0.001  ~ "**",
        p.value < 0.05   ~ "*",
        TRUE             ~ NA_character_
      )
    ) %>%
    filter(!is.na(label)) %>%
    mutate(
      group1 = layer1,
      group2 = layer2,
      y.position = y_max + row_number() * 0.07 * y_range
    )
  
  # -------------------------
  # Plot
  # -------------------------
  p <- ggplot() +
    geom_jitter(
      data = data,
      aes_string(x = "layer", y = response_var),
      width = 0.15, size = 0.2, alpha = 0.3
    ) +
    geom_errorbar(
      data = summary_data,
      aes(x = layer, ymin = mean - sem, ymax = mean + sem),
      width = 0.3, linewidth = 0.15
    ) +
    geom_segment(
      data = summary_data,
      aes(x = x_num - 0.25, xend = x_num + 0.25, y = mean, yend = mean),
      linewidth = 0.3
    ) +
    stat_pvalue_manual(sig_labels, label = "label", hide.ns = TRUE) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 6) +
    labs(
      x = label_lookup["layer"],
      y = label_lookup[response_var]
    )
  
  print(p)
  
  if (save_plot) {
    save_fig(p, response_var, w, h)
  }
  
  return(list(
    model = model,
    ICC   = icc,
    emmeans = emm,
    pairwise = pairwise,
    plot = p
  ))
}


# =========================
# --- Run analyses ---
# =========================
model_mean_intensity   <- analyze_parameter(df_imaris, "mean_wfa_intensity")
model_net_sphericity   <- analyze_parameter(df_imaris, "net_sphericity")
model_integr_intensity <- analyze_parameter(df_imaris, "integrated_wfa_intensity", normal = FALSE, w = 3.65/2.54)
model_intensity_cv     <- analyze_parameter(df_imaris, "intensity_cv")
model_volume           <- analyze_parameter(df_imaris, "net_volume_um3", normal = FALSE)


# ============================================================
# --- PV effect on WFA intensity ---
# ============================================================
df_imaris_clean <- df_imaris %>%
  filter(!is.na(pv), !is.na(layer))

model_wfa <- lmer(
  mean_wfa_intensity ~ layer * pv + (1 | pat),
  data = df_imaris_clean
)

summary(model_wfa)
anova(model_wfa)

emm_layer_pv <- emmeans(model_wfa, ~ pv | layer)


# =========================
# --- Plot EMMs ---
# =========================
emm_df <- as.data.frame(emm_layer_pv)

plot_emm <- ggplot() +
  geom_jitter(
    data = df_imaris_clean,
    aes(x = pv, y = mean_wfa_intensity),
    width = 0.15, alpha = 0.4, size = 0.3
  ) +
  geom_segment(
    data = emm_df,
    aes(
      x = as.numeric(factor(pv)) - 0.3,
      xend = as.numeric(factor(pv)) + 0.3,
      y = emmean,
      yend = emmean
    ),
    color = "red",
    linewidth = 0.3
  ) +
  geom_errorbar(
    data = emm_df,
    aes(x = pv, y = emmean, ymin = lower.CL, ymax = upper.CL),
    width = 0.3, color = "red"
  ) +
  facet_wrap(~layer, strip.position = "bottom") +
  scale_x_discrete(labels = c("no" = "PV-", "yes" = "PV+")) +
  theme_classic(base_size = 6)

print(plot_emm)
save_fig(plot_emm, "Fig_4_G_WFA_Intensity", w = 4.5/2.54, h = 3.6/2.54)
