######################
######################
######################
###### WRAPPERS ######
######################
######################
######################
fig_wrap <- function(folder, ...) {
  ggsave(..., units = "in", dpi = 300, dev = "png")
}

######################
######################
######################
#### MAIN FIGURES ####
######################
######################
######################
annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf,
                               ymin = -Inf, ymax = Inf, data) {
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = ggplot2:::GeomCustomAnn, inherit.aes = TRUE,
        params = list(grob = grob, xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax))
}

make_ranges <- function(fish_pic, df_, px, py, plot_a_r, x_r = 0.2) {
  usr <- unlist(df_[1, ])
  x_p <- usr[1] + px * (usr[2] - usr[1])
  x_p <- c(x_p, x_p + x_r * (usr[2] - usr[1]))
  y_p <- usr[3] + py * (usr[4] - usr[3])
  pic_pin <- dim(fish_pic)
  a_r <- pic_pin[1] / pic_pin[2]
  y_p <- c(y_p, y_p + (usr[4] - usr[3]) * x_r * a_r / plot_a_r)
  data.frame(x_min = x_p[1], x_max = x_p[2], y_min = y_p[1], y_max = y_p[2])
}

make_fig_1 <- function(model_list) {
  id_names <- list(
    `3` = substitute(atop(italic("Plectropomus") * " + " * italic("Variola") *
                            " spp.", "(Coral trout)")),
    `4` = "Serranidae (Rockcods)",
    `5` = substitute(atop(italic("Lethrinus miniatus") * " + " *
                            italic("nebulosus"),
                          "(Redthroat and Spangled emperors)")),
    `6` = "Lethrinidae (Emperors)", `7` = "Lutjanidae (Tropical snappers)",
    `9` = "Labridae (Wrasses)"
  )

  dp_id <- purrr::map_chr(id_names, function(x) {
    paste0(deparse(x), collapse = "")
  })
  dp_id <- factor(dp_id, levels = dp_id[c(6, 4, 3, 5, 2, 1)])
  posts <- purrr::map_dfr(model_list, function(y) {
    purrr::map_dfr(y, function(x) {
      brms::posterior_samples(x) %>%
        dplyr::mutate(hu_int_open = (b_hu_Intercept + 0.5 * b_hu_mean_live) %>%
                        plogis,
                      int_open = exp(b_Intercept)) %>%
        dplyr::select(int_open = int_open, cc = b_mean_live,
                      fish_open = b_ln_x, hu_int_open, hu_cc = b_hu_mean_live)
    }, .id = "Time lag")
  }, .id = "Fish group") %>%
    plyr::mutate(id_names = dp_id[`Fish group`])
  out <- ggplot(data = posts) +
    geom_density(mapping = aes(x = fish_open, fill = `Time lag`,
                 colour = `Time lag`), adjust = 2, alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey30") +
    facet_wrap(~id_names, nrow = 4, scales = "free",
               labeller = label_parsed) +
    labs(x = "Fish-biomass-removal scaling exponent", y = "Density") +
    scale_colour_brewer(name  = "Time lag (yrs)",
                        palette = "GnBu", direction = -1) +
    scale_fill_brewer(name  = "Time lag (yrs)",
                        palette = "GnBu", direction = -1) +
    theme_classic() +
    theme(axis.title = element_text(size = 15),
          strip.background = element_blank(),
          aspect.ratio = 0.5) +
    coord_cartesian(clip = "off")
  # should be in alphabetical orders of labels
  limits <- ggplot_build(out)$layout$panel_params %>%
    purrr::map_dfr(function(x) {
      data.frame(x_min = x$x.range[1], x_max = x$x.range[2],
                 y_min = x$y.range[1], y_max = x$y.range[2])
    })
  fish_groups <- levels(posts$id_names)
  fish_pics <- vector(mode = "list", length = length(fish_groups))
  plot_a_r <- ggplot_build(out)$plot$theme$aspect.ratio
  for (i in seq_along(fish_groups)) {
    df_ <- grep(fish_groups[i], posts$id_names, fixed = TRUE)[1]
    fish_pic <- dir("fish_pics", pattern = posts$`Fish group`[df_],
                    full.names = TRUE) %>%
      png::readPNG()
    ranges <- make_ranges(fish_pic, limits[i, ], px = 0.02, py = 0.75, plot_a_r,
                          x_r = 0.3)
    fish_pics[[i]] <- fish_pic %>%
      grid::rasterGrob(interpolate = TRUE) %>%
      annotation_custom2(xmin = ranges$x_min, xmax = ranges$x_max,
                         ymin = ranges$y_min, ymax = ranges$y_max,
                         data = posts[df_, ])
  }
  out + fish_pics[[1]] + fish_pics[[2]] + fish_pics[[3]] +
    fish_pics[[4]] + fish_pics[[5]] + fish_pics[[6]]
}

annotation_custom2 <- function(grob, xmin = -Inf, xmax = Inf,
                               ymin = -Inf, ymax = Inf, data) {
  layer(data = data, stat = StatIdentity, position = PositionIdentity,
        geom = ggplot2:::GeomCustomAnn, inherit.aes = TRUE,
        params = list(grob = grob, xmin = xmin, xmax = xmax,
                      ymin = ymin, ymax = ymax))
}

post_change <- function(model) {
  data.frame(Zone = c("Open", "Closed")) %>%
    brms::posterior_epred(model, newdata = ., re_formula = NA) %>%
    data.frame %>%
    dplyr::rename(Open = X1, Closed = X2) %>%
    dplyr::mutate(change = (Closed - Open) / Open * 100)
}

change_quants <- function(df_) {
  probs <- c(0.025, 0.2, 0.5, 0.8, 0.975)
  qts <- df_$change %>%
    quantile(prob = probs) %>%
    data.frame(quantiles = ., Group = unique(df_$Group))
  dens <- density(df_$change, adjust = 1.5) %>%
    {
      data.frame(x = .$x, y = .$y / max(.$y))
    } %>%
    dplyr::mutate(Group = unique(df_$Group),
                  quant = factor(findInterval(x, qts$quantiles)))
  y_dash <- dens %>%
    dplyr::group_by(quant) %>%
    dplyr::summarise(y = y[which.max(x)])
  qts$y <- y_dash$y[1:5]
  list(dens_data = dens, quantile = qts)
}

re_code_groups <- function(df_, add_ = FALSE, ...) {
  new_levels <- c(
    "Labridae (Wrasses)", "Lutjanidae (Tropical Snappers)", 
    "Lethrinidae (Emperors)", 
    "Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)", 
    "Serranidae (Rockcods)", "Plectropomus and Variola spp (Coral trout)"
  )
  if (add_) {
    new_levels <- paste0(new_levels, " ", ...)
  }
  df_$Group <- factor(df_$Group, levels = new_levels)
  df_
}

png_placement <- function(..., x_min, x_max, df_) {
  amp_ <- x_max - x_min
  x_max <- x_min - (amp_ * 0.02)
  x_min <- x_min - (amp_ * 0.25)
  png::readPNG(...) %>%
    grid::rasterGrob(interpolate = TRUE) %>%
    annotation_custom2(xmin = x_min, xmax = x_max, ymin = 0.05,
                       ymax = 0.95, data = df_)
}

filter_data <- function(df_, tag_) {
  df_ %>%
    dplyr::filter(grepl(tag_, Group))
}

make_facet_panel <- function(model, subtit_, tit_, return_data = FALSE, ...) {
  posts_list_ <- model %>%
    purrr::map_dfr(post_change, .id = "Group") %>%
    split(f = .$Group) %>%
    purrr::map(change_quants)
  posts_resp_ <- posts_list_ %>%
    purrr::map_dfr("dens_data") %>%
    re_code_groups(...)
  posts_qts_ <- posts_list_ %>%
    purrr::map_dfr("quantile") %>%
    re_code_groups(...)
  y_max_ <- posts_resp_ %>%
    split(f = .$Group) %>%
    purrr::map_dfr(function(df_) {
      data.frame(ymax = max(df_$y) + 0.2 * max(df_$y))
    }, .id = "Group") %>%
    re_code_groups(...)
  blues <- RColorBrewer::brewer.pal(9, "Blues")
  blues <- blues[c(2, 5, 9, 9, 5, 2)]
  out <- ggplot() +
    geom_segment(data = posts_qts_, color = "grey30", size = 0.1,
                 mapping = aes(y = 0, yend = y, x = quantiles,
                               xend = quantiles)) +
    geom_ribbon(data = posts_resp_,
                mapping = aes(x = x, ymin = 0, ymax = y, fill = quant)) +
    scale_fill_manual(values = blues, guide = "none") +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey30") +
    scale_y_continuous(breaks = 0.5) +
    facet_wrap(~Group, ncol = 1, scales = "free_y") +
    labs(x = "", y = "", title = tit_,
         subtitle = subtit_) +
    theme_classic() +
    theme(axis.text.x = element_text(size = 12),
          axis.text.y = element_blank(),
          axis.title = element_text(size = 15),
          axis.title.x = element_text(vjust = -0.5),
          axis.title.y = element_text(vjust = 1.5,
                                      margin = margin(t = 0, r = 1, b = 0,
                                                      l = 0, unit = "in")),
          axis.line = element_line(size = 0.5, colour = "black"),
          strip.text = element_blank(),
          plot.title = element_text(face = "bold")) +
    coord_cartesian(clip = "off")
  if (return_data) {
    return(ggplot2::ggplot_build(out)$data[[2]])
  }
  limits <- ggplot_build(out)$layout$panel_params %>%
    purrr::map_dfr(function(x) {
      data.frame(x_min = x$x.range[1], x_max = x$x.range[2],
                 y_min = x$y.range[1], y_max = x$y.range[2])
    }) %>%
    dplyr::slice(1)
  lab <- png_placement("fish_pics/9_Choerodon_anchorago_bw.png",
                       x_min = limits$x_min, x_max = limits$x_max,
                       df_ = filter_data(posts_resp_, "Labridae"))
  lut <- png_placement("fish_pics/7_Lutjanidae_Lutjanus.gibbus.png",
                       x_min = limits$x_min, x_max = limits$x_max,
                       df_ = filter_data(posts_resp_, "Lutjanidae"))
  let <- png_placement("fish_pics/6_Lethrinidae.png",
                       x_min = limits$x_min, x_max = limits$x_max,
                       df_ = filter_data(posts_resp_, "Lethrinidae"))
  let_ <- png_placement("fish_pics/5_Lethrinus_miniatus_vec.png",
                        x_min = limits$x_min, x_max = limits$x_max,
                        df_ = filter_data(posts_resp_, "Lethrinus"))
  ser <- png_placement("fish_pics/4_Serranidae.png",
                       x_min = limits$x_min, x_max = limits$x_max,
                       df_ = filter_data(posts_resp_, "Serranidae"))
  ctr <- png_placement("fish_pics/3_plectropomus_leopardus.png",
                       x_min = limits$x_min, x_max = limits$x_max,
                       df_ = filter_data(posts_resp_, "Plectropomus"))
  out + lab + lut + let + let_ + ser + ctr
}

make_fig_2 <- function(fish_models, return_data = FALSE) {
  models_biomass <- fish_models$models$biomass
  models_abun <- fish_models$models$abun
  models_length <- fish_models$models$length
  a <- make_facet_panel(models_biomass, subtit_ = "Biomass", tit_ = "a",
                        return_data = return_data, add_ = TRUE, "Biomass")
  b <- make_facet_panel(models_abun, subtit_ = "Density", tit_ = "b",
                        return_data = return_data)
  c <- make_facet_panel(models_length, subtit_ = "Length", tit_ = "c",
                        return_data = return_data, add_ = TRUE, "Length")
  if (return_data) {
    plots_ <- list(Biomass = a, Density = b, Length = c)
    return(purrr::map_dfr(plots_, identity, .id = "Response variable"))
  }
  plots_ <- list(a, b, c)
  tg <- grid::textGrob
  gg <- grid::gpar
  vp <- grid::viewport
  gu <- grid::unit
  my_margin <- theme(plot.margin = gu(rep(0.05, 4), "in"))
  x_labb <- tg("Percent change (U - F / F) (%)",
               gp = gg(fontsize = 20, font = 8))
  y_labb <- tg("Posterior density", rot = 90, gp = gg(fontsize = 20, font = 8))
  ggsave("output/figures/fig_2.png",
         gridExtra::grid.arrange(grobs = lapply(plots_, "+", my_margin),
                          bottom = x_labb, left = y_labb,
                          ncol = 1, nrow = length(plots_),
                          padding = gu(0.5, "line"),
                          vp = vp(width = 0.95, height = 0.95)),
         width = 6.5, height = 12, dpi = 300, dev = "png")
}

make_fig_3 <- function(model, return_data = FALSE) {
  fixefs <- brms::fixef(model)
  cc_fix <- fixefs["live_coral_prop", c("Estimate", "Q2.5", "Q97.5")] %>%
    sapply(rounded, 2)
  plots <- plot(brms::conditional_effects(model),
                plot = FALSE, line_args = list(colour = "black",
                                               linetype = 2))
  if (return_data) {
    return(
      purrr::map_dfr(plots, function (x) {
        x$data %>%
          dplyr::select(status, cot_count, live_coral_prop, report_year,
                        fullreef_id, cond__, effect1__, estimate__, se__,
                        lower__, upper__) %>%
          dplyr::mutate(effect1__ = as.character(effect1__))
      }, .id = "plot")
    )
  }
  ypos <- plots[[1]]$data$estimate__
  a <- plots[[1]] +
    labs(x = "Reef zoning status", y = "", title = "a") +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          plot.title = element_text(face = "bold")) +
    scale_x_discrete(labels = c("Unfished", "Fished")) +
    annotate("text", x = 1.02, y = ypos[1] + 0.00002,
             label = expression(paste("" == italic("e")^{italic(beta[0]) - 0.45 %.% 0.22})),
             hjust = 0, vjust = 0.5, size = 4) +
    annotate("text", x = 1.98, y = ypos[2] + 0.00002,
             label = expression(paste(italic("e")^{italic(beta[0]) + italic(beta[1]) - 0.45 %.% 0.22} == "")),
             hjust = 1, vjust = 0.5, size = 4)
  int <- round(plots[[2]]$data$estimate__[1], 3)
  b <- plots[[2]] +
    labs(x = "Manta tow coral cover (prop)", y = "", title = "b") +
    theme_bw() +
    ylim(c(0.0004, 0.0018)) +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          plot.title = element_text(face = "bold")) +
    annotate("text", x = 0.32, y = 0.0017,
             label = deparse(substitute(italic(y) == a %.% italic(x)^beta[2],
                                        list(a = int))),
             hjust = 0, vjust = 0, size = 4, parse = TRUE) +
    annotate("text", x = 0.32, y = 0.00155,
             label = deparse(substitute(beta[2] == a,
                                        list(a = cc_fix[["Estimate"]]))),
             hjust = 0, vjust = 0, size = 4, parse = TRUE) +
    annotate("text", x = 0.32, y = 0.0014,
             label = deparse(substitute("95% C.I." == a - b,
                                        list(a = cc_fix[["Q2.5"]],
                                             b = cc_fix[["Q97.5"]]))),
             hjust = 0, vjust = 0, size = 4, parse = TRUE)
  tg <- grid::textGrob
  ga <- gridExtra::grid.arrange
  gg <- grid::gpar
  y_lab <- tg(substitute("Predicted CoTS counts (2 mins"^-1 * ")"),
              rot = 90, gp = gg(fontsize = 15, font = 8))
  ga(a, b, ncol = 2, left = y_lab)
}

#######################
#######################
#######################
# SUPPLEMENTARY FIGURES
#######################
#######################
#######################
make_fig_s2 <- function(model_list) {
  id_names <- list(`3` = substitute(italic("Plectropomus") * " and " *
                                      italic("Variola") * " spp. " *
                                        "(Coral trout)"),
                   `4` = "Serranidae (Rockcods)",
                   `5` = substitute(italic("Lethrinus miniatus") * " and " *
                                      italic("L. nebulosus") *
                                        " (Redthroat and Spangled emperors)"),
                   `6` = "Lethrinidae (Emperors)",
                   `7` = "Lutjanidae (Tropical snappers)",
                   `9` = "Labridae (Wrasses)")
  plot_data <- model_list %>%
    plyr::ldply(function(x) {
      plyr::ldply(x, function(z) {
        a <- try(make_ggdata(z))
        if (!inherits(a, "try-error")) {
          a
        }
      }, .id = "t_lag")
    }, .id = "fish_group") %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character)) %>%
    dplyr::mutate(t_lag = paste0("x = ", t_lag))
  # Fish order:
  # Labridae (wrasses)[9]
  # Lethrinidae (emperors)[6]
  # Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)[5]
  # Lutjanidae (tropical snappers)[7] 
  # Serranidae (rockcods)[4]
  # Plectropomus and Variola spp (Coral trout)[3]
  fish_groups <- c("9", "6", "5", "7", "4", "3")
  plots_ <- list()
  for (i in seq_along(fish_groups)) {
    message("Plots for fish group ", fish_groups[i], "\n")
    plots_[[i]] <- plot_data %>%
      dplyr::filter(fish_group == fish_groups[i]) %>%
      group_main_effs_fig_s2(id_names = id_names)
  }
  tg <- grid::textGrob
  gg <- grid::gpar
  vp <- grid::viewport
  ac <- ggplot2::annotation_custom
  gt <- grid::grid.text
  gs <- grid::grid.points
  gl <- grid::grid.lines
  gu <- grid::unit
  my_margin <- theme(plot.margin = gu(rep(0.05, 4), "in"))
  x_labb <- tg(substitute("ln[Fish removal + 1] @ yr t"),
               gp = gg(fontsize = 20, font = 8))
  y_labb <- tg(substitute("ln[CoTS density] (min"^-1 * ") @ yr t + x"),
               rot = 90, gp = gg(fontsize = 20, font = 8))
  gridExtra::grid.arrange(grobs = lapply(plots_, "+", my_margin),
                          bottom = x_labb, left = y_labb,
                          ncol = 1, nrow = length(model_list),
                          padding = gu(0.5, "line"),
                          vp = vp(width = 0.95, height = 0.95))
}

group_main_effs_fig_s2 <- function(data, id_names) {
  all_ys <- data %>%
    dplyr::mutate(ln_cots = log(mean_cots)) %>%
    dplyr::select(ln_cots, y_pred, y_pred_l, y_pred_u, y_pols) %>%
    unlist %>% {
      range(.[!is.na(.) & . > -Inf], na.rm = TRUE)
    }
  y_max <- max(all_ys) + (max(all_ys) - min(all_ys)) * 0.1
  ggti <- id_names[[unique(data$fish_group)]]
  g_ <- ggplot() +
    geom_point(data = data %>% dplyr::filter(mean_cots > 0, !is.na(mean_cots)),
               mapping = aes(x = ln_x, y = log(mean_cots)),
               shape = 22, fill = "white", size = 1.5, colour = "grey30",
               alpha = 0.8) +
    geom_polygon(data = data %>% dplyr::filter(!is.na(y_pols)),
                 mapping = aes(x = x_pols, y = y_pols),
                 fill = "white", alpha = 0.5) +
    geom_line(data = data %>% dplyr::filter(!is.na(y_pred)),
              mapping = aes(x = x_pred, y = y_pred), linetype = 2,
              size = 0.5) +
    geom_line(data = data %>% dplyr::filter(!is.na(y_pred_l)),
              mapping = aes(x = x_pred, y = y_pred_l),
              linetype = 1, size = 0.1) +
    geom_line(data = data %>% dplyr::filter(!is.na(y_pred_u)),
              mapping = aes(x = x_pred, y = y_pred_u),
              linetype = 1, size = 0.1) +
    geom_text(data = data %>% dplyr::filter(!is.na(r2)),
              mapping = aes(label = t_lag),
              x = Inf, y = Inf, hjust = 6, vjust = 2, size = 3) +
    geom_text(data = data %>% dplyr::filter(!is.na(r2)),
              mapping = aes(label = r2_fct(r2)), parse = TRUE,
              x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3)
  if (any(!is.na(data$y_mark))) {
    g_ <- g_ +
      geom_line(data = data %>% dplyr::filter(is_mark, mark_type == "x"),
                mapping = aes(x = x_mark, y = y_mark), linetype = 2,
                colour = "tomato", lwd = 0.5) +
      geom_line(data = data %>% dplyr::filter(is_mark, mark_type == "y"),
                mapping = aes(x = x_mark, y = y_mark), linetype = 2,
                colour = "tomato", lwd = 0.5)
  }
  g_ +
    ylim(c(min(all_ys), y_max)) +
    labs(x = "", y = "", title = ggti, parse = TRUE) +
    facet_wrap(~t_lag, ncol = 6) +
    theme_bw() +
    theme(axis.text = element_text(size = 9),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.position = "none",
          panel.spacing = unit(1.5, "lines"),
          strip.text = element_blank()) +
    coord_cartesian(clip = "off")
}

r2_fct <- function(x) {
  lapply(x, function(z) deparse(substitute(italic(R^2) == a, list(a = z))))
}

make_ggdata <- function(model, hurdle = FALSE) {
  r2 <- brms::bayes_R2(model) %>%
    data.frame() %>%
    sapply(rounded, 2)
  df <- model$data
  obrk_threshold <- calc_outbrk_cutoff(model)
  if (!is.na(obrk_threshold)) {
    mark <- data.frame(type = rep(c("x", "y"), each = 2),
                       x = c(rep(obrk_threshold, 2),
                             -Inf, obrk_threshold),
                       y = c(-Inf, log(0.05), rep(log(0.05), 2)),
                       is_mark = TRUE, stringsAsFactors = FALSE)
  } else {
    mark <- data.frame(type = rep(c("x", "y"), each = 2),
                       x = rep(NA, 4), y = rep(NA, 4),
                       is_mark = TRUE, stringsAsFactors = FALSE)
  }
  if (hurdle) {
    df$mean_cots <- ifelse(df$mean_cots == 0, 0, 1)
  } else {
    df <- df %>%
      add_empty_cols("r2") %>%
      add_empty_rows(1) %>%
      assign_new_values(col = "r2", new_values = r2[["Estimate"]])
    mean_cc <- mean(df$mean_live, na.rm = TRUE)
  }
  if (hurdle) {
    df_i <- df %>%
      dplyr::filter(!is.na(mean_cots))
    nd <- data.frame(mean_live = seq(min(df_i$mean_live),
                                     max(df_i$mean_live),
                                     length.out = 50),
                     ln_x = mean(df_i$ln_x))
    fits <- fitted(model, dpar = "hu", newdata = nd,
                   robust = TRUE, re_formula = NA) %>%
      cbind(nd) %>%
      dplyr::mutate(x = mean_live)
  } else {
    df_i <- df %>%
      dplyr::filter(mean_cots > 0, !is.na(mean_cots))
    nd <- data.frame(ln_x = seq(min(df_i$ln_x, na.rm = TRUE),
                                max(df_i$ln_x, na.rm = TRUE),
                                length.out = 50),
                     mean_live = mean_cc)
    fits <- fitted(model, newdata = nd, dpar = "mu", scale = "linear",
                   robust = TRUE, re_formula = NA) %>%
      cbind(nd) %>%
      dplyr::mutate(x = ln_x)
  }
  pols <- fits %>% {
    data.frame(x = c(.$x, rev(.$x)),
               y = c(.$Q2.5, rev(.$Q97.5)))
  }
  df_ <- add_empty_rows(df, 1, bind = FALSE) %>%
    add_empty_cols("y_pred") %>%
    add_empty_cols("x_pred") %>%
    add_empty_cols("y_pred_l") %>%
    add_empty_cols("y_pred_u") %>%
    add_empty_rows(nrow(fits)) %>%
    assign_new_values(col = "y_pred", new_values = fits$Estimate) %>%
    assign_new_values(col = "x_pred", new_values = fits$x) %>%
    assign_new_values(col = "y_pred_l", new_values = fits$Q2.5) %>%
    assign_new_values(col = "y_pred_u", new_values = fits$Q97.5) %>%
    add_empty_cols("y_pols") %>%
    add_empty_cols("x_pols") %>%
    add_empty_rows(nrow(pols)) %>%
    assign_new_values(col = "y_pols", new_values = pols$y) %>%
    assign_new_values(col = "x_pols", new_values = pols$x) %>%
    add_empty_cols("y_mark") %>%
    add_empty_cols("x_mark") %>%
    add_empty_cols("is_mark") %>%
    add_empty_cols("mark_type") %>%
    add_empty_rows(nrow(mark)) %>%
    assign_new_values(col = "y_mark", new_values = mark$y) %>%
    assign_new_values(col = "x_mark", new_values = mark$x) %>%
    assign_new_values(col = "is_mark", new_values = mark$is_mark) %>%
    assign_new_values(col = "mark_type", new_values = mark$type) %>%
    dplyr::slice(-1)
  df %>%
    add_empty_cols("y_pred") %>%
    add_empty_cols("x_pred") %>%
    add_empty_cols("y_pred_l") %>%
    add_empty_cols("y_pred_u") %>%
    add_empty_cols("y_pols") %>%
    add_empty_cols("x_pols") %>%
    add_empty_cols("y_mark") %>%
    add_empty_cols("x_mark") %>%
    add_empty_cols("is_mark") %>%
    add_empty_cols("mark_type") %>%
    rbind(df_)
}

calc_outbrk_cutoff <- function(model, cut_off = log(0.05)) {
  # outbreak levels start at 0.1 CoTS per tow (2 min), but data is at the 1
  # min scale, so it needs to be 0.05 cots / tow, on the log scale
  # that is why cut_off = log(0.05) in function `calc_outbrk_cutoff`
  df <- model$data %>%
    dplyr::filter(mean_cots > 0, !is.na(mean_cots))
  nd <- data.frame(ln_x = seq(min(df$ln_x, na.rm = TRUE),
                              max(df$ln_x, na.rm = TRUE), length.out = 1e3),
                   mean_live = mean(df$mean_live))
  fits <- fitted(model, newdata = nd, dpar = "mu", scale = "linear",
                 robust = TRUE, re_formula = NA) %>%
    cbind(nd)
  tmp_ <- abs(fits[, "Estimate"] - cut_off)
  # use 0.01 diff as tolerance
  if (tmp_[which.min(tmp_)] < 0.01) {
    nd$ln_x[which.min(tmp_)]
  } else {
    NA
  }
}

assign_new_values <- function(data, pos, col, new_values) {
  if (missing(pos)) {
    start <- nrow(data) - length(new_values)
    pos <- (start + 1):(start + length(new_values))
  }
  data[pos, col] <- new_values
  data
}

add_empty_cols <- function(data, colname) {
  data[, colname] <- NA
  data
}

add_empty_rows <- function(data, nrows, bind = TRUE) {
  x <- data[seq_len(nrows), ]
  x[ ] <- NA
  if (bind) {
    rbind(data, x)
  } else {
    x
  }
}

make_fig_s3 <- function(model_list) {
  id_names <- list(`3` = substitute(italic("Plectropomus") * " and " *
                                      italic("Variola") * " spp. " *
                                        "(Coral trout)"),
                   `4` = "Serranidae (Rockcods)",
                   `5` = substitute(italic("Lethrinus miniatus") * " and " *
                                      italic("L. nebulosus") *
                                        " (Redthroat and Spangled emperors)"),
                   `6` = "Lethrinidae (Emperors)",
                   `7` = "Lutjanidae (Tropical snappers)",
                   `9` = "Labridae (Wrasses)")
  plot_data <- model_list %>%
    plyr::ldply(function(x) {
      plyr::ldply(x, function(z) {
        a <- try(make_ggdata(z, hurdle = TRUE))
        if (!inherits(a, "try-error")) {
          a
        }
      }, .id = "t_lag")
    }, .id = "fish_group") %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character)) %>%
    dplyr::mutate(t_lag = paste0("x = ", t_lag))
  # Fish order:
  # Labridae (wrasses)[9]
  # Lethrinidae (emperors)[6]
  # Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)[5]
  # Lutjanidae (tropical snappers)[7] 
  # Serranidae (rockcods)[4]
  # Plectropomus and Variola spp (Coral trout)[3]
  fish_groups <- c("9", "6", "5", "7", "4", "3")
  plots_ <- list()
  for (i in seq_along(fish_groups)) {
    plots_[[i]] <- plot_data %>%
      dplyr::filter(fish_group == fish_groups[i]) %>%
      group_hurdle_effs_fig_s3(id_names = id_names)
  }
  tg <- grid::textGrob
  gg <- grid::gpar
  vp <- grid::viewport
  ac <- ggplot2::annotation_custom
  gt <- grid::grid.text
  gs <- grid::grid.points
  gl <- grid::grid.lines
  gu <- grid::unit
  my_margin <- theme(plot.margin = gu(rep(0.05, 4), "in"))
  x_labb <- tg("Mean coral cover (proportion)",
               gp = gg(fontsize = 20, font = 8))
  y_labb <- tg("Probability of CoTS density being zero",
               rot = 90, gp = gg(fontsize = 20, font = 8))
  gridExtra::grid.arrange(grobs = lapply(plots_, "+", my_margin),
                          bottom = x_labb, left = y_labb,
                          ncol = 1, nrow = length(model_list),
                          padding = gu(0.5, "line"),
                          vp = vp(width = 0.9, height = 0.9))
}

group_hurdle_effs_fig_s3 <- function(data, id_names) {
  ggti <- id_names[[unique(data$fish_group)]]
  ggplot() +
    geom_point(data = data %>% dplyr::filter(!is.na(mean_cots)),
               mapping = aes(x = mean_live, y = mean_cots), shape = 22,
               fill = "white", size = 1.5, colour = "grey30", alpha = 0.8) +
    geom_polygon(data = data %>% dplyr::filter(!is.na(y_pols)),
                 mapping = aes(x = x_pols, y = y_pols),
                 fill = "white", alpha = 0.5) +
    geom_line(data = data %>% dplyr::filter(!is.na(y_pred)),
              mapping = aes(x = x_pred, y = y_pred), linetype = 2,
              size = 0.5) +
    geom_line(data = data %>% dplyr::filter(!is.na(y_pred_l)),
              mapping = aes(x = x_pred, y = y_pred_l),
              linetype = 1, size = 0.1) +
    geom_line(data = data %>% dplyr::filter(!is.na(y_pred_u)),
              mapping = aes(x = x_pred, y = y_pred_u),
              linetype = 1, size = 0.1) +
    geom_text(data = unique(data[, "t_lag", drop = FALSE]),
              mapping = aes(label = t_lag),
              x = Inf, y = Inf, hjust = 5.5, vjust = 2, size = 3) +
    ylim(c(0, 1.2)) +
    labs(x = "", y = "", title = ggti, parse = TRUE) +
    facet_wrap(~t_lag, ncol = 6) +
    theme_bw() +
    theme(axis.text = element_text(size = 9),
          axis.title = element_text(size = 15),
          legend.title = element_text(size = 15),
          legend.text = element_text(size = 12),
          legend.position = "none",
          panel.spacing = unit(1.5, "lines"),
          strip.text = element_blank()) +
    coord_cartesian(clip = "off")
}

make_fig_s4 <- function(model_list) {
  id_names <- list(`3` = substitute(italic("Plectropomus") * " and " *
                                      italic("Variola") * " spp. " *
                                        "(Coral trout)"),
                   `4` = "Serranidae (Rockcods)",
                   `5` = substitute(italic("Lethrinus miniatus") * " and " *
                                      italic("L. nebulosus") *
                                        " (Redthroat and Spangled emperors)"),
                   `6` = "Lethrinidae (Emperors)",
                   `7` = "Lutjanidae (Tropical snappers)",
                   `9` = "Labridae (Wrasses)")
  plot_data <- model_list %>%
    plyr::ldply(function(x) {
      plyr::ldply(x, function(z) {
        out <- brms::pp_check(z, type = "scatter_avg")
        out$data
      }, .id = "t_lag")
    }, .id = "fish_group") %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character)) %>%
    dplyr::mutate(t_lag = paste0("x = ", t_lag))
  # Fish order:
  # Labridae (wrasses)[9]
  # Lethrinidae (emperors)[6]
  # Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)[5]
  # Lutjanidae (tropical snappers)[7] 
  # Serranidae (rockcods)[4]
  # Plectropomus and Variola spp (Coral trout)[3]
  fish_groups <- c("9", "6", "5", "7", "4", "3")
  plots_ <- list()
  for (i in seq_along(fish_groups)) {
    plots_[[i]] <- plot_data %>%
      dplyr::filter(fish_group == fish_groups[i]) %>%
      group_pps_fig_s4(id_names = id_names)
  }
  tg <- grid::textGrob
  gg <- grid::gpar
  vp <- grid::viewport
  ac <- ggplot2::annotation_custom
  gt <- grid::grid.text
  gs <- grid::grid.points
  gl <- grid::grid.lines
  gu <- grid::unit
  my_margin <- theme(plot.margin = gu(rep(0.05, 4), "in"))
  y_labb <- tg("Observed", rot = 90, gp = gg(fontsize = 20, font = 7))
  x_labb <- tg("Predicted", gp = gg(fontsize = 20, font = 7))
  gridExtra::grid.arrange(grobs = lapply(plots_, "+", my_margin),
                          left = y_labb, ncol = 1, nrow = length(model_list),
                          bottom = x_labb, padding = gu(0.5, "line"),
                          vp = vp(width = 0.95, height = 0.95))
}

group_pps_fig_s4 <- function(data, id_names) {
  ggti <- id_names[[unique(data$fish_group)]]
  ggplot() +
    geom_abline(slope = 1, linetype = 2) +
    geom_point(data = data,
               mapping = aes(x = avg_y_rep, y = y), fill = "grey60",
               colour = "black", alpha = 0.5, shape = 21) +
    labs(x = "", y = "", title = ggti, parse = TRUE) +
    facet_wrap(~t_lag, ncol = 6, scales = "free") +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          panel.spacing = unit(1.5, "lines"),
          strip.text = element_blank())
}

make_figs_s5_s40 <- function(qdaf_models, qdaf_panalysis) {
  id_names <- list(
    `3` = "Plectropomus and Variola spp. (Coral trout)",
    `4` = "Serranidae (Rockcods)",  
    `5` = "Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)",
    `6` = "Lethrinidae (Emperors)",
    `7` = "Lutjanidae (Tropical Snappers)",
    `9` = "Labridae (Wrasses)"
  )
  n_ <- 4
  for (i in seq_along(qdaf_models)) {
    for (j in seq_along(qdaf_models[[i]])) {
      n_ <- n_ + 1
      name_i_ <- names(qdaf_models)[i]
      name_i_j_ <- names(qdaf_models[[i]])[j]
      tmp_mod_benchmark_ <- qdaf_models[[i]][[j]]
      tmp_sim_df_ <- qdaf_panalysis %>%
        dplyr::filter(target == name_i_, t_lag == name_i_j_)
      on <- paste0("output/figures/fig_s", n_, ".png")
      p_title <- id_names[[i]]
      p_subtitle <- paste0("Time lag: ", j, " year", ifelse(j == 1, "", "s"))
      fig_wrap(
        fig_out_folder, on,
        make_ind_fig_s5_s40(
          mod_benchmark_ = tmp_mod_benchmark_, sim_df_ = tmp_sim_df_,
          title = p_title, subtitle = p_subtitle
        ),
        width = 10.7, height = 12.6
      )
    }
  }
  invisible()
}

make_ind_fig_s5_s40 <- function(mod_benchmark_, sim_df_, ...) {
  pars_ <- c(`CoTS Intercept` = "b_Intercept",
             `CoTS ~ Biomass slope` = "b_ln_x",
             `CoTS ~ Coral slope` = "b_mean_live",
             `CoTS Intercept (hurdle)` = "b_hu_Intercept",
             `CoTS ~ Coral slope (hurdle)` = "b_hu_mean_live",
             `Grid site S.D.` = "sd_grid_site__Intercept",
             `Shape` = "shape")
  all_panel_df <- purrr::map_dfr(pars_, concat_sim_posts, df_ = sim_df_,
                                 mod_ = mod_benchmark_, .id = "parameter") %>%
    dplyr::mutate(parameter = forcats::fct_relevel(parameter, names(pars_)))
  he <- ggplot(data = all_panel_df, mapping = aes(x = value, fill = type)) +
    geom_density(alpha = 0.5, adjust = 1.5, trim = TRUE) +
    labs(x = "", y = "", title = "a)", fill = "Source: ") +
    scale_fill_manual(values = c("#E64B35FF", "#4DBBD5FF")) +
    facet_wrap(~parameter, ncol = 1, scales = "free") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom")
  hy <- post_hyp_plots(mod_benchmark_, pars_, "b)")
  de <- post_dens_plots(mod_benchmark_, pars_, "c)")
  tr <- post_chains_plots(mod_benchmark_, pars_, "d)")
  p_lay <- "
  ABCD
  "
  he + hy + de + tr +
    patchwork::plot_layout(design = p_lay) +
    patchwork::plot_annotation(...)
}

concat_sim_posts <- function(par_, df_, mod_) {
  tmp <- df_ %>%
    dplyr::select(tidyselect::all_of(par_), iter, type)
  brms::posterior_samples(mod_, pars = par_, fixed = TRUE) %>%
    dplyr::mutate(iter = NA, type = "Original") %>%
    rbind(tmp) %>%
    dplyr::rename("value" = tidyselect::all_of(par_))
}

post_hyp_plots <- function(mod_benchmark_, pars_, index_) {
  all_plots_df <- purrr::map_dfr(pars_, hyp_wrapper, model = mod_benchmark_,
                                 .id = "parameter") %>%
    dplyr::mutate(parameter = forcats::fct_relevel(parameter, names(pars_)))
  ggplot(data = all_plots_df, mapping = aes(x = values, fill = Type)) +
    geom_density(alpha = 0.5, adjust = 1.5) +
    labs(x = "", y = "", title = index_, fill = "Type: ") +
    scale_fill_manual(values = c("grey30", "white")) +
    facet_wrap(~parameter, ncol = 1, scales = "free") +
    theme_classic() +
    theme(axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position = "bottom")
}

hyp_wrapper <- function(x, model) {
  hyp_df <- hypothesis(model, paste0(x, " = 0"), class = "") %>%
    plot(plot = FALSE) %>%
    `[[`(1) %>%
    `[[`("data")
  # constrain x axis in plot for SD and Shape
  is_sd <- grepl("^sd_", x)
  if (is_sd | x == "shape") {
    hyp_df <- hyp_df %>%
      filter(values <= 10)
  }
  hyp_df
}

post_dens_plots <- function(mod_benchmark_, pars_, index_) {
  de_df_ <- brms::posterior_samples(mod_benchmark_, pars = pars_, fixed = TRUE) %>%
    tidyr::pivot_longer(cols = tidyselect::everything(), names_to = "parameter",
                        values_to = "value") %>%
    dplyr::mutate(parameter = names(pars_)[match(parameter, pars_)],
                  parameter = forcats::fct_relevel(parameter, names(pars_)))

  ggplot(data = de_df_, mapping = aes(x = value)) +
    geom_density(alpha = 0.5, adjust = 1.5, fill = "darkgrey") +
    labs(x = "", y = "", title = index_) +
    facet_wrap(~parameter, ncol = 1, scales = "free") +
    theme_classic()
}

post_chains_plots <- function(mod_benchmark_, pars_, index_) {
  tr_rs_ <- rstan::stan_trace(mod_benchmark_$fit, pars = pars_, ncol = 1)
  my_cols <- ggplot_build(tr_rs_)$data[[1]]["colour"] %>%
    dplyr::distinct() %>%
    dplyr::pull()
  tr_df_ <- tr_rs_ %>%
    `[[`("data") %>%
    dplyr::mutate(parameter = names(pars_)[match(parameter, pars_)],
                  parameter = forcats::fct_relevel(parameter, names(pars_)))
  ggplot(data = tr_df_, mapping = aes(x = iteration, y = value)) +
    geom_line(mapping = aes(colour = chain)) +
    scale_colour_manual(values = my_cols) +
    labs(x = "", y = "", colour = "Chain: ", title = index_) +
    facet_wrap(~parameter, ncol = 1, scales = "free") +
    theme_classic() +
    theme(legend.position = "bottom")
}

make_fig_s41 <- function(fish_models) {
  models_biomass <- fish_models$models$biomass
  models_abun <- fish_models$models$abun
  models_length <- fish_models$models$length
  names_ <- names(models_abun)
  all_stacked <- vector(mode = "list", length = length(names_))
  names(all_stacked) <- names_
  for (i in seq_along(names_)) {
    all_stacked[[i]] <- list(Density = get_model(names_[i], models_abun),
                             Biomass = get_model(names_[i], models_biomass),
                             Length = get_model(names_[i], models_length))
  }
  id_names <- list(
    `Labridae (Wrasses)` = "Labridae (Wrasses)",
    `Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)` =
      substitute(italic("Lethrinus miniatus") * " and " *
        italic("L. nebulosus") * " (Redthroat and Spangled emperors)"),
    `Lutjanidae (Tropical Snappers)` = "Lutjanidae (Tropical snappers)",
    `Plectropomus and Variola spp (Coral trout)` =
      substitute(italic("Plectropomus") * " and " * italic("Variola") *
                   " spp. " * "(Coral trout)"),
    `Serranidae (Rockcods)` = "Serranidae (Rockcods)",
    `Lethrinidae (Emperors)` = "Lethrinidae (Emperors)"
  )
  main_effs_fig_s41(all_stacked, id_names)
}

main_effs_fig_s41 <- function(model_list, ...) {
  mod_data <- model_list %>%
    purrr::map_dfr(function(x) {
      purrr::map_dfr(x, "data", .id = "var")
    }, .id = "fish_group") %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character),
                  Zone = recode(Zone, Open = "Fished", Closed = "Unfished"))
  pred_data <- model_list %>%
    purrr::map_dfr(function(x) {
      purrr::map_dfr(x, get_mean_preds, .id = "var")
    }, .id = "fish_group") %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character),
                  Zone = recode(Zone, Open = "Fished", Closed = "Unfished"))
  r2_data <- model_list %>%
    purrr::map_dfr(function(x) {
      purrr::map_dfr(x, r2_to_df, .id = "var")
    }, .id = "fish_group") %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character),
                  lab_ = sapply(Estimate, deparsed_r2))
  plots_ <- list()
  fish_groups <- c(
    "Labridae (Wrasses)", "Lutjanidae (Tropical Snappers)",
    "Lethrinidae (Emperors)",
    "Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)",
    "Serranidae (Rockcods)", "Plectropomus and Variola spp (Coral trout)"
  )
  for (i in seq_along(fish_groups)) {
    plots_[[i]] <- group_main_effs_fig_s41(
      filter_group(pred_data, fish_groups[i]),
      filter_group(mod_data, fish_groups[i]),
      filter_group(r2_data, fish_groups[i]),
      ...
    )
  }
  tg <- grid::textGrob
  gg <- grid::gpar
  vp <- grid::viewport
  gu <- grid::unit
  my_margin <- theme(plot.margin = gu(rep(0.05, 4), "in"))
  x_labb <- tg("Zoning status",
               gp = gg(fontsize = 20, font = 8))
  y_labb <- tg("Value",
               rot = 90, gp = gg(fontsize = 20, font = 8))
  gridExtra::grid.arrange(grobs = lapply(plots_, "+", my_margin),
                          bottom = x_labb, left = y_labb,
                          ncol = 1, nrow = length(model_list),
                          padding = gu(0.5, "line"),
                          vp = vp(width = 0.95, height = 0.95))
}

get_model <- function(tar_, models_list) {
  models_list[[grep(tar_, names(models_list), fixed = TRUE)]]
}

get_mean_preds <- function(model) {
  plot(brms::conditional_effects(model, effect = "Zone"),
       plot = FALSE)[[1]]$data
}

filter_group <- function(df_, tar_) {
  df_ %>%
    dplyr::filter(fish_group == tar_)
}

r2_to_df <- function(model) {
  brms::bayes_R2(model) %>%
    data.frame %>%
    sapply(rounded, 2)
}

deparsed_r2 <- function(x) {
  deparse(substitute(italic(R^2) == a, list(a = x)))
}

group_main_effs_fig_s41 <- function(pred_data, mod_data, r2_data, id_names) {
  ggti <- id_names[[unique(mod_data$fish_group)]]
  plot_ <- ggplot() +
    geom_jitter(data = mod_data,
                mapping = aes(x = Zone, y = Value, shape = Zone),
                size = 2, fill = "grey60", colour = "grey60",
                alpha = 0.8, width = 0.2, show.legend = FALSE) +
    geom_errorbar(data = pred_data,
                  mapping = aes(x = Zone, ymin = lower__, ymax = upper__),
                  colour = "black", lwd = 1.5, width = 0.2) +
    geom_point(data = pred_data, show.legend = FALSE,
               mapping = aes(x = Zone, y = estimate__, shape = Zone),
               colour = "black", fill = "tomato", size = 4) +
    scale_shape_manual(values = 21:22) +
    geom_text(data = r2_data, mapping = aes(label = lab_),
              x = -Inf, y = Inf, parse = TRUE,
              hjust = -0.3, vjust = 2) +
    labs(x = "", y = "", title = ggti, parse = TRUE) +
    facet_wrap(~var, ncol = 3, scales = "free") +
    scale_y_continuous(trans = "log1p") +
    theme_classic() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15)) +
    coord_cartesian(clip = "off")
}

make_fig_s42 <- function(fish_models) {
  models_biomass <- fish_models$models$biomass
  models_abun <- fish_models$models$abun
  models_length <- fish_models$models$length
  names_ <- names(models_abun)
  all_stacked <- vector(mode = "list", length = length(names_))
  names(all_stacked) <- names_
  for (i in seq_along(names_)) {
    all_stacked[[i]] <- list(Density = get_model(names_[i], models_abun),
                             Biomass = get_model(names_[i], models_biomass),
                             Length = get_model(names_[i], models_length))
  }
  id_names <- list(
    `Labridae (Wrasses)` = "Labridae (Wrasses)",
    `Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)` =
      substitute(italic("Lethrinus miniatus") * " and " *
        italic("L. nebulosus") * " (Redthroat and Spangled emperors)"),
    `Lutjanidae (Tropical Snappers)` = "Lutjanidae (Tropical snappers)",
    `Plectropomus and Variola spp (Coral trout)` =
      substitute(italic("Plectropomus") * " and " * italic("Variola") *
                   " spp. " * "(Coral trout)"),
    `Serranidae (Rockcods)` = "Serranidae (Rockcods)",
    `Lethrinidae (Emperors)` = "Lethrinidae (Emperors)"
  )
  plot_data <- purrr::map_dfr(all_stacked, function(x) {
      purrr::map_dfr(x, function(z) {
        out <- brms::pp_check(z, type = "scatter_avg")
        out$data
      }, .id = "var")
    }, .id = "fish_group") %>%
    dplyr::mutate(dplyr::across(where(is.factor), as.character))
  plots_ <- list()
  fish_groups <- c(
    "Labridae (Wrasses)", "Lutjanidae (Tropical Snappers)",
    "Lethrinidae (Emperors)",
    "Lethrinus miniatus and L. nebulosus (Redthroat and Spangled emperors)",
    "Serranidae (Rockcods)", "Plectropomus and Variola spp (Coral trout)"
  )
  for (i in seq_along(fish_groups)) {
    plots_[[i]] <- plot_data %>%
      dplyr::filter(fish_group == fish_groups[i]) %>%
      group_pps_fig_s42(id_names = id_names)
  }
  tg <- grid::textGrob
  gg <- grid::gpar
  vp <- grid::viewport
  ac <- ggplot2::annotation_custom
  gt <- grid::grid.text
  gs <- grid::grid.points
  gl <- grid::grid.lines
  gu <- grid::unit
  my_margin <- theme(plot.margin = gu(rep(0.05, 4), "in"))
  y_labb <- tg("Observed", rot = 90, gp = gg(fontsize = 20, font = 7))
  x_labb <- tg("Predicted", gp = gg(fontsize = 20, font = 7))
  gridExtra::grid.arrange(grobs = lapply(plots_, "+", my_margin),
                          left = y_labb, ncol = 1, nrow = length(all_stacked),
                          bottom = x_labb, padding = gu(0.5, "line"),
                          vp = vp(width = 0.95, height = 0.95))
}

group_pps_fig_s42 <- function(data, id_names) {
  ggti <- id_names[[unique(data$fish_group)]]
  ggplot() +
    geom_abline(slope = 1, linetype = 2) +
    geom_point(data = data,
               mapping = aes(x = avg_y_rep, y = y), fill = "grey60",
               colour = "black", alpha = 0.5, shape = 21) +
    labs(x = "", y = "", title = ggti, parse = TRUE) +
    geom_text(data = data %>% dplyr::distinct(var),
              mapping = aes(label = var),
              x = Inf, y = Inf, hjust = 1.3, vjust = 2, size = 3) +
    facet_wrap(~var, ncol = 3, scales = "free") +
    theme_bw() +
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 15),
          panel.spacing = unit(1.5, "lines"),
          strip.text = element_blank())
}

make_figs_s43_s60 <- function(fish_models) {
  fish_models <- fish_models$models
  n_ <- 42
  for (i in seq_along(fish_models)) {
    for (j in seq_along(fish_models[[i]])) {
      n_ <- n_ + 1
      name_i_ <- Hmisc::capitalize(names(fish_models)[i])
      name_i_j_ <- names(fish_models[[i]])[j]
      tmp_mod_benchmark_ <- fish_models[[i]][[j]]
      on <- paste0("output/figures/fig_s", n_, ".png")
      fig_wrap(
        fig_out_folder, on,
        make_ind_fig_s43_s60(
          mod_benchmark_ = tmp_mod_benchmark_, type = name_i_,
          title = name_i_j_, subtitle = name_i_
        ),
        width = 8.025, height = ifelse(name_i_ == "Length", 9, 12.6)
      )
    }
  }
  invisible()
}

make_ind_fig_s43_s60 <- function(mod_benchmark_, type, ...) {
  pars_ <- c(
    `CoTS Intercept` = "b_Intercept", `CoTS Intercept F` = "b_ZoneOpen",
     `CoTS Intercept (hurdle)` = "b_hu_Intercept",
     `CoTS Intercept F (hurdle)` = "b_hu_ZoneOpen",
     `Year S.D.` = "sd_cYear__Intercept", `Pair S.D.` = "sd_Pair__Intercept",
     `Reef S.D.` = "sd_Reef__Intercept", `Site S.D.` = "sd_Site__Intercept",
     `Shape` = "shape"
  )
  if (type == "Length") {
    pars_ <- pars_[-grep("hu", pars_)]
  }
  hy <- post_hyp_plots(mod_benchmark_, pars_, "a)")
  de <- post_dens_plots(mod_benchmark_, pars_, "b)")
  tr <- post_chains_plots(mod_benchmark_, pars_, "c)")
  hy + de + tr +
    patchwork::plot_layout(design = "ABC") +
    patchwork::plot_annotation(...)
}

make_fig_s61 <- function(df_, predictions, model_r2) {
  df_$pred_count_mean <- apply(predictions, 1, mean)
  modr2 <- data.frame(model_r2) %>%
    sapply(rounded, 2)
  ggplot(data = df_) +
    geom_point(mapping = aes(x = cot_count,
                             y = pred_count_mean),
               show.legend = FALSE) +
    scale_y_continuous(trans = "log1p") +
    scale_x_continuous(trans = "log1p") +
    labs(x = substitute("Observed CoTS counts (2 mins"^-1 * ")"),
         y = substitute("Predicted CoTS counts (2 mins"^-1 * ")")) +
    geom_abline() +
    theme_bw() +
    theme(axis.text = element_text(size = 10),
          axis.title = element_text(size = 15)) +
    annotate("text", x = 100, y = 0.1,
             label = deparse(substitute(italic(R^2) == a,
                             list(a = modr2[["Estimate"]]))),
             hjust = 0, vjust = 0, size = 5, parse = TRUE)
}

make_fig_s62 <- function(mod_benchmark_) {
  pars_ <- c(`CoTS Intercept` = "b_Intercept",
             `CoTS Intercept F` = "b_statuso",
             `CoTS ~ Coral slope` = "b_live_coral_prop",
             `Reef-year S.D.` = "sd_report_year:fullreef_id__Intercept",
             `Shape` = "shape")
  hy <- post_hyp_plots(mod_benchmark_, pars_, "a)")
  de <- post_dens_plots(mod_benchmark_, pars_, "b)")
  tr <- post_chains_plots(mod_benchmark_, pars_, "c)")
  hy + de + tr +
    patchwork::plot_layout(design = "ABC")
}

make_fig_s63 <- function (sim_list, ...) {
  count_zeros <- function(x) sum(x == 0)
  test_generic(
    sim_list = sim_list, summary = count_zeros,
    method_name = "DHARMa zero-inflation test via comparison to expected zeros with simulation under H0 = fitted model", ...
  )
}

test_generic <- function(sim_list, summary,
                         alternative = c("two.sided", "greater", "less"),
                         plot = TRUE,
                         method_name = "DHARMa generic simulation test") {
  ensure_dharma <- DHARMa:::ensureDHARMa
  get_p_val <- DHARMa:::getP
  out <- list()
  out$data.name <- deparse(substitute(sim_list))
  sim_list <- ensure_dharma(sim_list, convert = "Model")
  alternative <- match.arg(alternative)
  observed <- summary(sim_list$observedResponse)
  simulated <- apply(sim_list$simulatedResponse, 2, summary)
  p <- get_p_val(simulated = simulated, observed = observed,
                 alternative = alternative)
  out$statistic <- c(ratioObsSim = observed / mean(simulated))
  out$method <- method_name
  out$alternative <- alternative
  out$p.value <- p
  class(out) <- "htest"
  if (plot) {
    x_lab_ <- paste("Simulated values, red line = fitted model. p-value (",
                    out$alternative, ") = ", out$p.value, sep = "")
    plot_title <- gsub("(.{1,50})(\\s|$)", "\\1\n", method_name)
    x_range <- range(simulated, observed, na.rm = TRUE)
    data.frame(simulated = simulated) %>%
      ggplot(data = .) +
        geom_histogram(mapping = aes(x = simulated), colour = "grey60",
                       bins = max(round(sim_list$nSim / 5), 20)) +
        geom_vline(xintercept = observed, colour = "tomato") +
        labs(x = x_lab_, y = "", title = plot_title,
             subtitle = "Fitted vs. Simulated") +
        xlim(x_range) +
        theme_bw() +
        theme(plot.title = element_text(size = 9),
              plot.subtitle = element_text(size = 9))
  } else {
    out
  }
}
