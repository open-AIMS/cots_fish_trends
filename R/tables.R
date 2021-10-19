make_ed_table_2 <- function(stacked_models, csv_out_folder, csv_out_path) {
  ed_tab_2 <- purrr::map_dfr(stacked_models, function(y) {
    df <- purrr::map_dfr(y, function(x) {
      posts <- brms::posterior_samples(x)
      all_pars <- list(
        int_open = posts$b_Intercept,
        cc = posts$b_mean_live,
        fish_open = posts$b_ln_x,
        hu_int_open = (posts$b_hu_Intercept + 0.5 *
          posts$b_hu_mean_live) %>% plogis,
        hu_cc = posts$b_hu_mean_live
      )
      test_names <- list(
        int_open = "-", cc = "< 0", fish_open = "> 0", hu_int_open = "-",
        hu_cc = "> 0"
      )
      all_tests <- list(
        int_open = NULL, cc = `<`, fish_open = `>`, hu_int_open = NULL,
        hu_cc = `>`
      )
      all_names <- list(
        int_open = "CoTS density (O)", cc = "CC slope",
        fish_open = "Fish density slope (O)",
        hu_int_open = "Prob. 0 (O)", hu_cc = "Hurdle CC slope"
      )
      if (!all(names(all_pars) == names(all_tests)) |
            !all(names(all_pars) == names(all_names))) {
        stop("parameter names do not match among lists")
      }
      out_list <- vector(mode = "list", length = length(all_pars))
      for (i in seq_along(all_pars)) {
        one <- ifelse(names(all_pars)[i] == "diff_hu_int_oc", TRUE, FALSE)
        out_list[[i]] <- make_out_edtab3_data(all_pars[[i]], all_tests[[i]],
                                              test_names[[i]], all_names[[i]],
                                              one = one)
      }
      do.call("rbind.data.frame", out_list)
    }, .id = "Time lag") %>%
      dplyr::mutate(`Time lag` = ifelse(duplicated(`Time lag`), "",
                                        `Time lag`))
  }, .id = "Analysis") %>%
    dplyr::mutate(Analysis = as.character(Analysis),
                  Analysis = ifelse(duplicated(Analysis), "", Analysis))
  col_names <- c(
    `3` = "Plectropomus and Variola spp. (Coral trout)",
    `4` = "Serranidae (Rockcods)",
    `5` = paste0("Lethrinus miniatus and L. nebulosus ",
                 "(Redthroat and Spangled emperors)"),
    `6` = "Lethrinidae (Emperors)", `7` = "Lutjanidae (Tropical snappers)",
    `9` = "Labridae (Wrasses)", ""
  )
  ed_tab_2$Analysis <- col_names[ed_tab_2$Analysis]
  ed_tab_2$Analysis[is.na(ed_tab_2$Analysis)] <- ""
  write.csv(ed_tab_2, file.path("output/tables", csv_out_path),
            row.names = FALSE)
  ed_tab_2
}

get_model <- function(tar_, models_list) {
  models_list[[grep(tar_, names(models_list), fixed = TRUE)]]
}

make_out_edtab3_data <- function(par, test, test_name, name, one = FALSE) {
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

make_ed_table_3 <- function(fish_models, csv_out_folder, csv_out_path) {
  models_abun <- fish_models$models$abun
  models_biomass <- fish_models$models$biomass
  models_length <- fish_models$models$length
  names_ <- names(models_abun)
  all_stacked <- vector(mode = "list", length = length(names_))
  names(all_stacked) <- names_
  for (i in seq_along(names_)) {
    all_stacked[[i]] <- list(Density = get_model(names_[i], models_abun),
                             Biomass = get_model(names_[i], models_biomass),
                             Length = get_model(names_[i], models_length))
  }
  ed_tab_3 <- purrr::map_dfr(all_stacked, function(y) {
    df <- purrr::map_dfr(y, function(x) {
      posts <- brms::posterior_samples(x)
      all_pars <- list(
        int_closed = posts$b_Intercept,
        int_open = posts$b_Intercept + posts$b_ZoneOpen,
        diff_int_oc = posts$b_ZoneOpen
      )
      test_names <- list(
        int_closed = "-", int_open = "-", diff_int_oc = "< 0"
      )
      all_tests <- list(
        int_closed = NULL, int_open = NULL, diff_int_oc = `<`
      )
      all_names <- list(
        int_closed = "Response (U)", int_open = "Response (F)",
        diff_int_oc = "Delta Response (F - U)"
      )
      if (x$family$family == "hurdle_gamma") {
        all_pars <- c(all_pars, list(
          hu_int_closed = plogis(posts$b_hu_Intercept),
          hu_int_open = plogis(posts$b_hu_Intercept + posts$b_hu_ZoneOpen),
          diff_hu_int_oc = plogis(posts$b_hu_Intercept + posts$b_hu_ZoneOpen) /
            plogis(posts$b_hu_Intercept)
        ))
        test_names <- c(test_names, list(
          hu_int_closed = "-", hu_int_open = "-", diff_hu_int_oc = "> 1"
        ))
        all_tests <- c(all_tests, list(
          hu_int_closed = NULL, hu_int_open = NULL, diff_hu_int_oc = `>`
        ))
        all_names <- c(all_names, list(
          hu_int_closed = "Prob. 0 (U)", hu_int_open = "Prob. 0 (F)",
          diff_hu_int_oc = "Odds ratio prob. 0 (F / U)"
        ))
      }
      if (!all(names(all_pars) == names(all_tests)) |
            !all(names(all_pars) == names(all_names))) {
        stop("parameter names do not match among lists")
      }
      out_list <- vector(mode = "list", length = length(all_pars))
      for (i in seq_along(all_pars)) {
        one <- ifelse(names(all_pars)[i] == "diff_hu_int_oc", TRUE, FALSE)
        out_list[[i]] <- make_out_edtab3_data(all_pars[[i]], all_tests[[i]],
                                              test_names[[i]], all_names[[i]],
                                              one = one)
      }
      do.call("rbind.data.frame", out_list)
    }, .id = "Response") %>%
      dplyr::mutate(`Response` = ifelse(duplicated(`Response`), "",
                                        `Response`))
  }, .id = "Analysis") %>%
    dplyr::mutate(Analysis = as.character(Analysis),
                  Analysis = ifelse(duplicated(Analysis), "", Analysis))
  write.csv(ed_tab_3, file.path("output/tables", csv_out_path),
            row.names = FALSE)
  ed_tab_3
}
