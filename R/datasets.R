# This code was originally written to download the LTMP dataset from
# internal servers via SQL query. This is maintained for the sake of
# transparency only. It will not run on any other machine outside of AIMS
# servers; The datasets were downloaded on the 2nd of September 2020.
check_and_load_file <- function(type, file_tag, fish_str, manta_str) {
  sql_out <- file.path("output", "sql", paste0(file_tag, ".sql"))
  csv_out <- file.path("data", paste0(file_tag, ".csv"))
  if (!file.exists(csv_out)) {
    sql_call(sql_out, csv_out,
             ifelse(type == "fish",
                    fish_str,
                    manta_str))
  }
  csv_out %>%
    read.csv(strip.white = TRUE, stringsAsFactors = FALSE) %>%
    dplyr::rename_with(us_tolower) %>%
    dplyr::mutate(dplyr::across(where(is.character), us_tolower))
}

sql_call <- function(sql_out, csv_out, the_str) {
  dir.create(file.path("output", "sql"),
             recursive = TRUE,
             showWarnings = FALSE)
  writeLines(the_str, sql_out)
  system(paste("java -jar data/dbExport/dbExport.jar",
               sql_out,
               csv_out,
               "REEF REEFMON"))
}

read_ltmp <- function(type, ...) {
  miss_r <- c("11049s" = "c", "11162s" = "c", "mcsweeney_reef" = "o",
              "monsoon_reef" = "o", "mantis_reef" = "o",
              "turtle_group_reef_(no_2)" = "c", "15022s" = "o",
              "long_reef(15019)" = "o",
              "pasco_reef" = "o", "rocky_island_reef(15054)" = "c",
              "15023s" = "o", "15037s" = "c", "harrier_reef" = "o",
              "13093a" = "o", "corbett_reef" = "c", "hedge_reef" = "o",
              "curlew_reef" = "o", "knight_reef" = "o",
              "pine_peak_reef" = "o", "prudhoe_island_reef" = "o",
              "temple_island_reef" = "o", "temple_shoal" = "o")
  data <- check_and_load_file(type, ...) %>%
    dplyr::mutate(openorclosed = ifelse(reef_name %in% names(miss_r),
                                        miss_r[reef_name],
                                        openorclosed),
                  status = ifelse(report_year > 2004,
                                  openorclosed_after2004,
                                  openorclosed))
  if (type == "fish") {
    data %>%
      dplyr::filter(site_no != "", fish_code != "") %>%
      dplyr::mutate(phase1 = ifelse(phase == "t", "t", ""))
  } else if (type == "cots") {
    coral_props <- c("0" = 0, "1" = 0.05, "1l" = 0.025,
                     "1u" = 0.075, "2" = 0.2, "2l" = 0.15,
                     "2u" = 0.25, "3" = 0.4, "3l" = 0.35,
                     "3u" = 0.45, "4" = 0.625, "4l" = 0.5625,
                     "4u" = 0.6875, "5" = 0.875, "5l" = 0.8125,
                     "5u" = 0.9375)
    data %>%
      dplyr::filter(p_code != "ts",
                    shelf != "") %>%
      dplyr::mutate(scar = ifelse(scar == "", "a", scar),
                    live_coral = ifelse(live_coral == "", "0",
                                        live_coral),
                    live_coral_prop = coral_props[live_coral])
  }
}

dump_fish_str <- function() {
  "select p_code,S.Sample_Id,a_sector,shelf,
     reef_name,fullreef_id,reef_lat,reef_long,
     openorclosed,openorclosed_after2004,
     report_year,sample_date,site_no,transect_no,
     fish_code,genus,species,abundance,length,phase
   from rm_fish05 f, v_rm_sample s
   where f.sample_id=s.sample_id
     and (fish_code Not Like 'AET_ROGA' And
          fish_code Not Like 'ANY_LEUC' And
          fish_code Not Like 'APR_VIRE' And
          fish_code Not Like 'ACA_SP' And
          fish_code Not Like 'CHA_OXYC' And
          fish_code Not Like 'CHA_SEME' And
          fish_code Not Like 'CHL_LABI' And
          fish_code Not Like 'CHS_FRON' And
          fish_code Not Like 'CRO_ALTI' And
          fish_code Not Like 'CTE_BINO' And
          fish_code Not Like 'ddc8df90-601f-11e5-8242-ecf4bb65a2a2' And
          fish_code Not Like 'DIP_BIFA' And
          fish_code Not Like 'GNA_AURO' And
          fish_code Not Like 'GRA_ALBI' And
          fish_code Not Like 'GYN_SPP' And
          fish_code Not Like 'IST_DECO' And
          fish_code Not Like 'LET_ERUS' And
          fish_code Not Like 'NAS_ANBR' And
          fish_code Not Like 'PMS_OLIG' And
          fish_code Not Like 'POM_PAVO' And
          fish_code Not Like 'POM_UN' And
          fish_code Not Like 'PSU_TUKA' And
          fish_code Not Like 'SAR_RUBR' And
          fish_code Not Like 'SCA_SP')
     and visit_no is not null
     and p_code in ('RM','RAP','RMRAP')
     and report_year>1994
   order by a_sector,shelf,reef_name,report_year,
   site_no,transect_no,fish_code,length,phase"
}

dump_manta_str <- function() {
  "select p_code,S.Sample_Id,a_sector,shelf,
     reef_name,fullreef_id,reef_lat,reef_long,
     openorclosed,openorclosed_after2004,
     report_year,sample_date,visit_no,tow_seq_no,
     cot_count,live_coral,scar
   from v_rm_sample s inner join
     rm_manta m on s.sample_id = m.sample_id
   where sample_class in ('K', 'C', 'G') or
     sample_class is null
     and report_year>1984
   order by a_sector,shelf,reef_name,
            report_year,tow_seq_no"
}

fish_str <- dump_fish_str()
manta_str <- dump_manta_str()
ltmp_cots <- read_ltmp("cots", "ltmp_cots_3250", fish_str, manta_str)
write.csv(ltmp_cots, "data/ltmp_cots.csv", row.names = FALSE)

load("data/biomass_data.RData")
write.csv(data, "data/ltmp_fish_zoning.csv", row.names = FALSE)
