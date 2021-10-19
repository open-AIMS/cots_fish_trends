plan <- drake::drake_plan(
  ##################
  ###### Data ######
  ##################
  ltmp_fish = read_file("data/ltmp_fish_zoning.csv") %>%
    dplyr::mutate(across(c(Sector, Pair, Reef,
                           Site, Zone, cYear), as.factor)),
  ltmp_cots = read_file("data/ltmp_cots.csv"),
  qdaf_data = wrangle_qdaf_data(
    fgroups_path = "CoTS_Fisheries_Combined_12Jan21_groups.csv",
    qdaf_path = "CoTS_Fisheries_Combined_12Jan21.xlsx"
  ),

  ##################
  #### Analyses ####
  ##################
  base_qdaf_model = compile_qdaf_model(),
  qdaf_models = stack_qdaf_models(qdaf_data, base_qdaf_model),
  
  fish_models = run_all_fish_models(ltmp_fish),

  cots_model = predict_cots_counts(ltmp_cots),
  cots_predictions = make_brms_predictions(ltmp_cots, cots_model),
  cots_model_r2 = brms::bayes_R2(cots_model, subset = 1:1000),
  cots_model_residuals = make_dharma(cots_predictions, ltmp_cots),

  qdaf_panalysis = purrr::map_dfr(qdaf_models, run_qdaf_power_analysis,
                                  .id = "target"),

  ##################
  ##### Figures ####
  ##################
  fig_out_folder = dir.create("output/figures", recursive = TRUE,
                              showWarnings = FALSE),
  fig_1 = fig_wrap(fig_out_folder, "output/figures/fig_1.png",
                   make_fig_1(qdaf_models), width = 8.9, height = 8.4),
  fig_2 = {
    fig_out_folder
    make_fig_2(fish_models)
  },
  fig_3 = fig_wrap(fig_out_folder, "output/figures/fig_3.png",
                   make_fig_3(cots_model), width = 9, height = 4),
  fig_s2 = fig_wrap(fig_out_folder, "output/figures/fig_s2.png",
                    make_fig_s2(qdaf_models), width = 12.7, height = 14.4),
  fig_s3 = fig_wrap(fig_out_folder, "output/figures/fig_s3.png",
                    make_fig_s3(qdaf_models), width = 12.7, height = 14.4),
  fig_s4 = fig_wrap(fig_out_folder, "output/figures/fig_s4.png",
                    make_fig_s4(qdaf_models),
                    width = 12.7, height = 14.4),
  figs_s5_s40 = make_figs_s5_s40(qdaf_models, qdaf_panalysis),
  fig_s41 = fig_wrap(fig_out_folder, "output/figures/fig_s41.png",
                     make_fig_s41(fish_models), width = 10, height = 20),
  fig_s42 = fig_wrap(fig_out_folder, "output/figures/fig_s42.png",
                     make_fig_s42(fish_models), width = 10, height = 20),
  figs_s43_s60 = make_figs_s43_s60(fish_models),
  fig_s61 = fig_wrap(fig_out_folder, "output/figures/fig_s61.png",
                     make_fig_s61(ltmp_cots, cots_predictions, cots_model_r2),
                     width = 6, height = 5.5),
  fig_s62 = fig_wrap(fig_out_folder, "output/figures/fig_s62.png",
                     make_fig_s62(cots_model), width = 8.025, height = 9),
  fig_s63 = fig_wrap(fig_out_folder, "output/figures/fig_s63.png",
                     make_fig_s63(cots_model_residuals), width = 6, height = 6),
  
  ##################
  ##### Tables #####
  ##################
  csv_out_folder = dir.create("output/tables", recursive = TRUE,
                              showWarnings = FALSE),
  ed_table_2 = make_ed_table_2(qdaf_models, csv_out_folder, "ed_table_2.csv"),
  ed_table_3 = make_ed_table_3(fish_models, csv_out_folder, "ed_table_3.csv"),

  ############################
  ##### Manuscript stats #####
  ############################
  fish_model_fold_change = fold_change_all(fish_models),
  n_signif_models = {
    ed_table_2
    ed_table_3
    show_n_signif_models("output/tables/ed_table_2.csv",
                         "output/tables/ed_table_3.csv")
  },
  cots_stats = show_cots_stats(cots_model)
)
