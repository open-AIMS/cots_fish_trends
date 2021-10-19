##################
# HELPER FUNCTIONS
##################
us_tolower <- function(x) {
  tolower(gsub(" ", "_", x))
}

read_file <- function(filename, ...) {
  read.csv(filename, header = TRUE, stringsAsFactors = FALSE, ...)
}

rounded <- function(value, precision = 1) {
  sprintf(paste0("%.", precision, "f"), round(value, precision))
}

#####################
# QDAF FISHERIES DATA
# RELATED TO FIGURE 1
#####################
split_groups <- function(x) {
  strsplit(x, ", ", fixed = TRUE)[[1]]
}

split_name <- function(vec) {
  vec %>%
    sapply(function(x)strsplit(x, "\\.")[[1]][2], USE.NAMES = FALSE)
}

sum_data_repetitions <- function(data) {
  plyr::ddply(data, .(year, grid_site, target), function(df_) {
    df_$kg <- sum(df_$kg, na.rm = TRUE)
    df_$num <- sum(df_$num, na.rm = TRUE)
    df_[1, ]
  }) %>%
    dplyr::select(-sector, -fishery, -species_common_name, -species_caab_code)
}

lag_match_data <- function(data, t_lag) {
  check_d <- expand.grid(y1 = data$year,
                         y2 = data$year) %>%
    dplyr::mutate(diff = y2 - y1) %>%
    dplyr::filter(diff == t_lag)
  if (nrow(check_d) > 0) {
    fish_density <- data %>%
      dplyr::filter(year %in% check_d$y1)
    cots_density <- data %>%
      dplyr::filter(year %in% check_d$y2)
    type <- ifelse(unique(data$type) == "weight", "kg", "num")
    out <- cots_density %>%
      dplyr::select(year, mean_cots, mean_live) %>%
      cbind(fish_density %>%
              dplyr::select(tidyselect::all_of(type))) %>%
      dplyr::mutate(fish_type = tidyselect::all_of(type),
                    t_lag = t_lag)
    names(out)[names(out) == type] <- "mean_fish"
    out
  }
}

compile_qdaf_model <- function() {
  data <- data.frame(mean_cots = rgamma(10, 2, 0.1),
                     ln_x = rnorm(10),
                     mean_live = plogis(rlogis(10)),
                     grid_site = rep(letters[1:5], each = 2),
                     stringsAsFactors = FALSE)
  priors <- prior(normal(0, 5), class = "b") +
    prior(normal(0, 5), class = "b", dpar = "hu") +
    prior(logistic(0, 1), class = "b", coef = "Intercept",
          dpar = "hu") +
    prior(gamma(2, 0.5), class = "sd") +
    prior(gamma(2, 0.5), class = "shape")
  brms::brm(bf(mean_cots ~ 0 + Intercept + ln_x + mean_live +
                 (1 | grid_site),
               hu ~ 0 + Intercept + mean_live), data = data,
            family = hurdle_gamma, prior = priors, chains = 0,
            sample_prior = "yes")
}

run_r_t_lag <- function(data, model) {
  tag <- unique(data$target)
  lg <- unique(data$t_lag)
  message("Model for group ", tag, "; time lag = ", lg)
  update(model, newdata = data, chains = 4, iter = 5e3,
         control = list(max_treedepth = 20,
                        adapt_delta = 0.99))
}

stack_qdaf_models <- function(df_, base_qdaf_model) {
  model_list <- plyr::dlply(df_, .(target, t_lag), run_r_t_lag,
                            base_qdaf_model)
  all_groups <- sort(unique(df_$target))
  all_stacked <- vector(mode = "list", length = length(all_groups))
  names(all_stacked) <- all_groups
  for (i in seq_along(all_stacked)) {
    to_find <- paste0("^", all_groups[i], "\\.")
    to_loop <- grep(to_find, names(model_list), value = TRUE)
    all_stacked[[i]] <- list()
    n <- 0
    for (j in seq_along(to_loop)) {
      if (!is.null(model_list[[to_loop[j]]])) {
        n <- n + 1
        all_stacked[[i]][[n]] <- model_list[[to_loop[j]]]
        names(all_stacked[[i]])[n] <- split_name(to_loop[j])
      }
    }
  }
  all_stacked
}

run_qdaf_power_analysis <- function(mod_list_) {
  purrr::map_dfr(mod_list_, function(mod_) {
    pred_df <- brms::posterior_predict(mod_, nsamples = 100)
    out_list <- vector(mode = "list", length = nrow(pred_df))
    for (i in seq_along(out_list)) {
      message("iteration ", i)
      tmp_mod <- mod_$data %>%
        dplyr::mutate(mean_cots = pred_df[i, ]) %>%
        update(mod_, newdata = .)
      to_keep <- c("b_Intercept", "b_hu_Intercept", "b_ln_x", "b_mean_live",
                   "b_hu_mean_live", "sd_grid_site__Intercept", "shape")
      out_list[[i]] <- brms::posterior_samples(tmp_mod, pars = to_keep,
                                               fixed = TRUE) %>%
        dplyr::sample_n(1e3) %>%
        dplyr::mutate(iter = i, type = "Simulated")
    }
    do.call("rbind.data.frame", out_list)
  }, .id = "t_lag")
}

named_vec <- function(x) {
  nms <- x %>%
    as.character %>%
    Hmisc::capitalize() %>%
    gsub("_", " ", x = .)
  names(x) <- nms
  x
}

is_signif <- function(x, y, tol = .Machine$double.eps^0.5) {
  x * y > tol
}

make_out_data <- function(par, test, test_name, name, one = FALSE) {
  if (!is.null(test)) {
    if (one) {
      prob_ <- (sum(test(par, 1)) / length(par) * 100) %>%
        rounded(2)
    } else {
      prob_ <- (sum(test(par, 0)) / length(par) * 100) %>%
        rounded(2)
    }
  } else {
    prob_ <- "-"
  }
  qts <- quantile(par, probs = c(0.025, 0.975))
  data.frame(Parameter = name, Estimate = rounded(mean(par), 2),
             Q2.5 = rounded(qts[[1]], 2), Q97.5 = rounded(qts[[2]], 2),
             `Test name` = test_name, `Test probability` = prob_,
             stringsAsFactors = FALSE, check.names = FALSE)
}

wrangle_qdaf_data <- function(fgroups_path, qdaf_path) {
  qdafg <- fgroups_path %>%
    file.path("data", .) %>%
    read_file(na.strings = "", check.names = FALSE) %>%
    dplyr::select(-NOTES) %>%
    dplyr::mutate(group = seq_len(dplyr::n()))
  qdafg_gen <- qdafg %>%
    dplyr::filter(!is.na(`SPECIES COMMON NAME`))
  qdafg_fam <- qdafg %>%
    dplyr::filter(is.na(`SPECIES COMMON NAME`))
  qdaf_path <- qdaf_path %>%
    file.path("data", .)
  qdafc <- qdaf_path %>%
    readxl::read_excel(sheet = "Long data_LTMP") %>%
    dplyr::rename_with(us_tolower) %>%
    dplyr::mutate(dplyr::across(where(is.character), us_tolower))
  qdaff <- qdaf_path %>%
    readxl::read_excel(sheet = "Long data_Fisheries") %>%
    dplyr::mutate(group_gen = 0, group_fam = 0, type = "")
  qdaff$FAMILY[is.na(qdaff$FAMILY)] <- "None"
  for (j in seq_len(nrow(qdafg_gen))) {
    all_sectors <- split_groups(qdafg_gen$SECTOR[j])
    all_fisheries <- split_groups(qdafg_gen$FISHERY[j])
    if (is.na(qdafg_gen$FAMILY[j])) {
      t_ <- qdaff$`SPECIES COMMON NAME` == qdafg_gen$`SPECIES COMMON NAME`[j] &
              qdaff$SECTOR %in% all_sectors & qdaff$FISHERY %in% all_fisheries
      qdaff[t_, "group_gen"] <- qdafg_gen$group[j]
      qdaff[t_, "type"] <- qdafg_gen$TYPE[j]
    } else {
      all_spp <- split_groups(qdafg_gen$`SPECIES COMMON NAME`[j])
      if (any(grepl("Coral trout", all_spp))) {
        all_spp <- c("Coral trout", "Cod - coronation trout",
                     "Trout - blue spot", "Trout - island",
                     "Trout - passionfruit")
      }
      t_ <- qdaff$`SPECIES COMMON NAME` %in% all_spp &
              qdaff$FAMILY == qdafg_gen$FAMILY[j] &
                qdaff$SECTOR %in% all_sectors &
                  qdaff$FISHERY %in% all_fisheries
      qdaff[t_, "group_gen"] <- qdafg_gen$group[j]
      qdaff[t_, "type"] <- qdafg_gen$TYPE[j]
    }
  }
  for (j in seq_len(nrow(qdafg_fam))) {
    all_sectors <- split_groups(qdafg_fam$SECTOR[j])
    all_fisheries <- split_groups(qdafg_fam$FISHERY[j])
    t_ <- qdaff$FAMILY == qdafg_fam$FAMILY[j] & qdaff$SECTOR %in% all_sectors &
            qdaff$FISHERY %in% all_fisheries
    qdaff[t_, "group_fam"] <- qdafg_fam$group[j]
    qdaff[t_, "type"] <- qdafg_fam$TYPE[j]
  }
  qdaff <- qdaff %>%
    dplyr::rename_with(us_tolower) %>%
    dplyr::mutate(dplyr::across(where(is.character), us_tolower)) %>%
    dplyr::select(-`retained_weight_(t)`) %>%
    dplyr::rename(num = `retained_number`, kg = `retained_weight_(kg)`,
                  discarded = `discarded_number`) %>%
    dplyr::filter(!(is.na(kg) & type == "weight"),
                  !(is.na(num) & type == "abundance"))
  # check
  # > sum(is.na(qdaff$kg) & qdaff$type == "weight")
  # > sum(is.na(qdaff$num) & qdaff$type == "abundance")
  qdaf_gen <- qdaff %>%
    dplyr::filter(group_gen != 0) %>%
    dplyr::select(-group_fam, target = group_gen) %>%
    sum_data_repetitions %>%
    dplyr::left_join(qdafc, by = c("year", "grid_site")) %>%
    dplyr::filter(!is.na(openorclosed))
  qdaf_fam <- qdaff %>%
    dplyr::filter(group_fam != 0) %>%
    dplyr::select(-group_gen, target = group_fam) %>%
    sum_data_repetitions %>%
    dplyr::left_join(qdafc, by = c("year", "grid_site")) %>%
    dplyr::filter(!is.na(openorclosed))
  qdaf <- rbind(qdaf_gen, qdaf_fam) %>%
    dplyr::filter(openorclosed == "o")
  all_t_lags <- 1:6
  names(all_t_lags) <- paste0("x = ", all_t_lags)
  purrr::map_dfr(all_t_lags, function(x, df_) {
    # mean_cots needs to be standardised to individuals / minute
    # original data sent by Frederieke is in individuals / 2 minutes
    # remove abundance-based fisheries data, really sparse and unreliable
    df_ %>%
      plyr::ddply(.(grid_site, target), lag_match_data,
                  t_lag = x) %>%
      dplyr::mutate(mean_cots = mean_cots / 2, ln_y = log(mean_cots + 1),
                    ln_x = log(mean_fish + 1), mean_live = mean_live / 1e2)
  }, df_ = qdaf, .id = "lag_lab") %>%
    dplyr::filter(target %in% c(3, 4, 5, 6, 7, 9))
}

#####################
# FISH LTMP DATA
# RELATED TO FIGURE 2
#####################
run_zone_model <- function(data, model) {
  tag <- unique(data$Group)
  message("Model for group: ", tag)
  update(model, newdata = data, chains = 4, iter = 5e3,
         control = list(max_treedepth = 20, adapt_delta = 0.99))
}

run_all_fish_models <- function(data) {
  kw <- c("wrasses", "emperors", "lutjanidae", "plectropomus", "serranidae")
  to_keep <- kw %>%
    sapply(grepl, data$Group, ignore.case = TRUE) %>%
    apply(1, any)
  df_ <- data %>%
    dplyr::filter(to_keep)
  # Biomass
  df_bm_ <- df_ %>%
    dplyr::filter(grepl("biomass", Group, ignore.case = TRUE)) %>%
    dplyr::mutate(Value = Value / 1000)
  priors <- prior(normal(0, 5), class = "b") +
    prior(gamma(1, 0.1), class = "sd") +
    prior(gamma(1, 0.1), class = "shape") +
    prior(normal(0, 5), class = "b", dpar = "hu") +
    prior(logistic(0, 1), class = "b", coef = "Intercept", dpar = "hu")
  mod_bm_ <- brms::bf(Value ~ 0 + Intercept + Zone + (1 | cYear) + (1 | Pair) +
                        (1 | Reef) + (1 | Site), hu ~ 0 + Intercept + Zone)
  mod_b_ <- brms::brm(mod_bm_, data = df_bm_, family = hurdle_gamma,
                      prior = priors, sample_prior = "yes", chains = 0)
  models_biomass <- df_bm_ %>%
    split(f = df_bm_$Group) %>%
    purrr::map(run_zone_model, model = mod_b_)
  # Length
  df_lg_ <- df_ %>%
    dplyr::filter(grepl("length", Group, ignore.case = TRUE),
                  !is.na(Value), Value > 0) %>%
    dplyr::mutate(Value = Value / 100)
  priors <- prior(normal(0, 5), class = "b") +
    prior(gamma(1, 0.1), class = "sd") +
    prior(gamma(1, 0.1), class = "shape")
  mod_l_ <- brms::brm(Value ~ 0 + Intercept + Zone + (1 | cYear) + (1 | Pair) +
                        (1 | Reef) + (1 | Site),
                      data = df_lg_, family = Gamma(link = "log"),
                      prior = priors, sample_prior = "yes", chains = 0)
  models_length <- df_lg_ %>%
    split(f = df_lg_$Group) %>%
    purrr::map(run_zone_model, model = mod_l_)
  # Abundance (i.e. density per 1,000 m2)
  df_ab_ <- df_ %>%
    dplyr::filter(!grepl("length", Group, ignore.case = TRUE),
                  !grepl("biomass", Group, ignore.case = TRUE)) %>%
    dplyr::mutate(Value = Value / 100)
  priors <- prior(normal(0, 5), class = "b") +
    prior(gamma(1, 0.1), class = "sd") +
    prior(gamma(1, 0.1), class = "shape") +
    prior(normal(0, 5), class = "b", dpar = "hu") +
    prior(logistic(0, 1), class = "b", coef = "Intercept", dpar = "hu")
  mod_ab_ <- brms::bf(Value ~ 0 + Intercept + Zone + (1 | cYear) + (1 | Pair) +
                        (1 | Reef) + (1 | Site), hu ~ 0 + Intercept + Zone)
  mod_a_ <- brms::brm(mod_ab_, data = df_ab_, family = hurdle_gamma,
                      prior = priors, sample_prior = "yes", chains = 0)
  models_abun <- df_ab_ %>%
    split(f = df_ab_$Group) %>%
    purrr::map(run_zone_model, model = mod_a_)
  list(models = list(biomass = models_biomass, length = models_length,
                     abun = models_abun))
}

#####################
# COTS LTMP DATA
# RELATED TO FIGURE 3
#####################
predict_cots_counts <- function(df_) {
  priors <- prior(normal(0, 5), class = "b") +
    prior(gamma(1, 0.1), class = "sd") +
    prior(gamma(1, 0.1), class = "shape")
  nb_form <- brms::bf(cot_count ~ 0 + Intercept + status + live_coral_prop +
                        (1 | report_year:fullreef_id))
  brms::brm(nb_form, data = df_, family = negbinomial(), prior = priors,
            iter = 1e4, control = list(max_treedepth = 15, adapt_delta = 0.9),
            sample_prior = "yes")
}

make_brms_predictions <- function(df_, model) {
  seqs <- round(seq(1, nrow(df_), length.out = 22))
  seqs[length(seqs)] <- seqs[length(seqs)] + 1
  out <- list()
  for (i in seq_along(seqs)[-length(seqs)]) {
    count <- seqs[i]:(seqs[i + 1] - 1)
    out[[i]] <- brms::posterior_predict(model, subset = 1:500,
                                        newdata = df_[count, ]) %>%
      t()
  }
  do.call(rbind.data.frame, out)
}

make_dharma <- function(preds, df_) {
  fitted_median <- apply(preds, 1, median)
  DHARMa::createDHARMa(simulatedResponse = as.matrix(preds),
                       observedResponse = df_$cot_count,
                       fittedPredictedResponse = fitted_median,
                       integerResponse = TRUE)
}

#####################
# GENERAL PAPER STATS
#####################
fold_change_all <- function(fish_models) {
  fold_change <- function(x) {
    brms::posterior_samples(x, pars = "b_ZoneOpen") %>%
      data.frame %>%
      dplyr::pull(b_ZoneOpen) %>% {
        1 / exp(.)
      } %>%
      ggdist::mean_hdi()
  }
  purrr::map_dfr(fish_models$models, function(x) {
    purrr::map_dfr(x, fold_change, .id = "fish") %>%
      dplyr::filter(ymin >= 1)
  }, .id = "Response")
}

show_n_signif_models <- function(fisheries_path, ltmpfish_path) {
  calc_n_signif <- function(df_, par_, prob_) {
    df_ %>%
      dplyr::filter(Parameter == par_, `Test probability` > prob_) %>%
      nrow
  }
  fisheries <- read_file(fisheries_path, check.names = FALSE,
                         na.strings = c("", "-", NA))
  ltmpfish <- read_file(ltmpfish_path, check.names = FALSE,
                        na.strings = c("", "-", NA))
  list(
    n_80_fisheries = calc_n_signif(fisheries, "Fish density slope (O)", 0.8),
    n_95_fisheries = calc_n_signif(fisheries, "Fish density slope (O)", 0.95),
    n_80_ltmpfish = calc_n_signif(ltmpfish, "Delta Response (F - U)", 0.8),
    n_95_ltmpfish = calc_n_signif(ltmpfish, "Delta Response (F - U)", 0.95)
  )
}

show_cots_stats <- function(cots_model) {
  pars_ <- c("b_statuso", "b_live_coral_prop",
             "sd_report_year:fullreef_id__Intercept")
  brms::posterior_samples(cots_model, pars = pars_, fixed = TRUE) %>%
    data.frame %>%
    dplyr::mutate(b_statuso = exp(b_statuso)) %>%
    apply(2, ggdist::median_hdci)
}
