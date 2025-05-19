################
# Preparations #
################
# this section creates the formulae for the analysis
# and the parameter table which we'll then loop over

create_frm <- function(covars) {
  as.formula(paste("Surv(t_start, t_end, OutcomeRow)", paste(covars, collapse=" + "), sep=" ~ "))
}

covars <- c("Age", "Gender", "Charlson" ,"bad_functional_state", "ImmunosupressionDiagnosis", "PPI")

sep_ab_bin <- c("(Aminoglycosides >= 1)"
                ,"(Augmentin >= 1)"
                ,"(CS1_2 >= 1)"
                ,"(CS3_4 >= 1)"
                ,"(Carbapenem >= 1)"
                ,"(Clindamycin >= 1)"
                ,"(Penicellins >= 1)"
                ,"(Pip_Tazo >= 1)"
                ,"(Quinolone >= 1)"
                ,"(Vanco_IV >= 1)"
                ,"(Other >= 1)"
)

sep_ab_cum <- c("Aminoglycosides"
                ,"Augmentin"
                ,"CS1_2"
                ,"CS3_4"
                ,"Carbapenem"
                ,"Clindamycin"
                ,"Penicellins"
                ,"Pip_Tazo"
                ,"Quinolone"
                ,"Vanco_IV"
                ,"Other"
)

screen_crude <- create_frm("CD_ScreenResult")
screen_adj <- create_frm(c(covars, "CD_ScreenResult"))

dot_crude <- create_frm("DoT")
dot_adj <- create_frm(c(covars, "DoT"))

any_ab_bin_crude <- create_frm("any_Ab")
any_ab_bin_adj <- create_frm(c(covars, "any_Ab"))

any_ab_cum_crude <- create_frm("num_Ab_cumulative")
any_ab_cum_adj <- create_frm(c(covars, "num_Ab_cumulative"))

sep_ab_bin_crude <- create_frm(sep_ab_bin)
sep_ab_bin_adj <- create_frm(c(covars, sep_ab_bin))

sep_ab_cum_crude <- create_frm(sep_ab_cum)
sep_ab_cum_adj <- create_frm(c(covars, sep_ab_cum))

# param table
pop <- c("everyone", "pos_screen", "neg_screen")
adjustment <- c("crude", "adjusted")
exposure <- c("screen", "DoT", "any_ab_bin", "any_ab_cum", "sep_ab_bin", "sep_ab_cum")

param_table <- 
  expand.grid(pop = pop, exposure = exposure, adjustment = adjustment, stringsAsFactors = F) %>% 
  filter(!(exposure == "screen" & pop != "everyone"))


#################
# Analysis Code #
#################
# this function performs the actual analysis per each row in the parameter table
do_one_analysis <- function(i, this_data) {
  
  pop <- param_table[i,1]
  exposure <- param_table[i,2]
  adjustment <- param_table[i,3]
  
  # get pop and outcome
  if(pop == "pos_screen") {this_data <- this_data %>% filter(CD_ScreenResult == 1)}
  if(pop == "neg_screen") {this_data <- this_data %>% filter(CD_ScreenResult == 0)}
  
  # get formula
  if (exposure == "screen" & adjustment == "crude") {this_frm <- screen_crude}
  if (exposure == "screen" & adjustment == "adjusted") {this_frm <- screen_adj}
  if (exposure == "DoT" & adjustment == "crude") {this_frm <- dot_crude}
  if (exposure == "DoT" & adjustment == "adjusted") {this_frm <- dot_adj}
  if (exposure == "any_ab_bin" & adjustment == "crude") {this_frm <- any_ab_bin_crude}
  if (exposure == "any_ab_bin" & adjustment == "adjusted") {this_frm <- any_ab_bin_adj}
  if (exposure == "any_ab_cum" & adjustment == "crude") {this_frm <- any_ab_cum_crude}
  if (exposure == "any_ab_cum" & adjustment == "adjusted") {this_frm <- any_ab_cum_adj}
  if (exposure == "sep_ab_bin" & adjustment == "crude") {this_frm <- sep_ab_bin_crude}
  if (exposure == "sep_ab_bin" & adjustment == "adjusted") {this_frm <- sep_ab_bin_adj}
  if (exposure == "sep_ab_cum" & adjustment == "crude") {this_frm <- sep_ab_cum_crude}
  if (exposure == "sep_ab_cum" & adjustment == "adjusted") {this_frm <- sep_ab_cum_adj}
  
  # do analysis
  this_fit <- coxph(formula = this_frm
                    ,data = this_data, x = T, y = T, model = T, id = ID, cluster = ID)
  
  print(this_fit)
  
  res <- 
    tidy(this_fit, exponentiate = T, conf.int = T) %>% 
    transmute(
      pop = pop
      ,exposure = exposure
      ,adjustment = adjustment 
      ,term
      ,est = glue("{round(estimate, 2)} ({round(conf.low, 2)}-{round(conf.high, 2)})")
    )
  
  if(exposure %in% c("screen", "DoT", "any_ab_bin", "any_ab_cum") & adjustment == "crude") {res <- res %>% slice(1)}
  if(exposure %in% c("screen", "DoT", "any_ab_bin", "any_ab_cum") & adjustment == "adjusted") {res <- res %>% slice(7)}
  if(exposure %in% c("sep_ab_bin", "sep_ab_cum") & adjustment == "crude") {res <- res %>% slice(1:11)}
  if(exposure %in% c("sep_ab_bin", "sep_ab_cum") & adjustment == "adjusted") {res <- res %>% slice(7:17)}
  
  res
  
}


################
# Run Analysis #
################
# in this section we call the analysis function for each row in the parameter table
# we do this 3 times:

# 1. For the full analysis
res <- bind_rows(lapply(1:nrow(param_table), FUN = do_one_analysis, this_data = data))

# 2. for the toxin outcome
res <- bind_rows(lapply(1:nrow(param_table), FUN = do_one_analysis, this_data = data_toxin))

# 3. and only for the first hospitalization of each patient
res <- bind_rows(lapply(1:nrow(param_table), FUN = do_one_analysis, this_data = data_one))

# results are exported after each run