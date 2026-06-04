####### Analysis shell script #####
# calls the four analyses
library(dplyr)
library(readr)
library(tidyr)
library(metafor) 
library(ggplot2) 
library(janitor)
library(forestplot)
library(clubSandwich) 

##### Preparatory work #######################
bhl_meta_full <- read_csv("data/bhl_meta.csv")

# Restrict to emotional eating outcomes:
bhl_meta_tmp <- bhl_meta_full %>% filter(StE_grp1 == 0)
bhl_meta_remainder <- bhl_meta_full %>% filter(StE_grp1 != 0)

# add a column with new id for each effect size
bhl_meta_tmp$effect_id <- 1:nrow(bhl_meta_tmp)

# add ids.
bhl_meta_tmp$intervention_id <- paste(bhl_meta_tmp$covidence_num,
                                      bhl_meta_tmp$intervention_number,
                                      sep = ".")
bhl_meta_tmp$outcome_id <- paste(bhl_meta_tmp$covidence_num,
                                 bhl_meta_tmp$intervention_number,
                                 bhl_meta_tmp$outcome_number,
                                 sep = ".")
bhl_meta_tmp$timepoint_id <- paste(bhl_meta_tmp$covidence_num,
                                   bhl_meta_tmp$intervention_number,
                                   bhl_meta_tmp$outcome_number,
                                   bhl_meta_tmp$timepoint_number,
                                   sep = ".")
bhl_meta_tmp$full_id <- paste(bhl_meta_tmp$outcome_id, 
                              bhl_meta_tmp$timepoint_number, sep=".")

## Create a flag 1365.1.1.1 - this is the only case where pre and post didn't
## match the sign.
bhl_meta_tmp <- bhl_meta_tmp %>% 
  mutate(mismatch_flag = case_when(full_id=="#1365.1.1.1" ~ T,
    T ~ F))
check_flag <- bhl_meta_tmp %>% filter(mismatch_flag==1) %>%
  dplyr::select(full_id, outcome, mean_post, mean_pre, final_cohens_d,
                final_hedges_g) 

bhl_meta <- bhl_meta_tmp %>% 
  arrange(covidence_num, intervention_number, 
          outcome_number, timepoint_number)  

# Count number of studies.
length(unique(bhl_meta$covidence_num))
# Count number of effects.
nrow(bhl_meta)

# What number of these have a .5 pre_post? To use in sensitivity analysis.
check_pp <- bhl_meta %>% filter(r_pre_post==.5)
length(unique(check_pp$covidence_num))

# Summarize study nesting:
bhl_summary <- bhl_meta %>% 
  distinct(covidence_num,Dep_MGrps, Dep_MTPs, Dep_MOut) %>%
  group_by(Dep_MGrps, Dep_MOut, Dep_MTPs) %>%
  summarize(study_count = n())



# Code follow-up category
bhl_meta <- bhl_meta %>% 
  mutate(follow_grp = case_when(FollowUpTP==0 ~ 0,
                                FollowUpTP %in% c(1, 2, 3, 4, 5)~ 1,
                                FollowUpTP %in% c(6, 7) ~ 2,
                                FollowUpTP %in% c(8, 9) ~3))

bhl_meta <- bhl_meta %>% 
  mutate(follow_grp2 = case_when(FollowUpTP==0 ~ 0,
                                 FollowUpTP==1 ~ 1,
                                FollowUpTP %in% c(2, 3, 4, 5)~ 2,
                                FollowUpTP %in% c(6, 7) ~ 3,
                                FollowUpTP %in% c(8, 9) ~ 4))

bhl_meta <- bhl_meta %>% 
  mutate(follow_grp3 = case_when(FollowUpTP==0 ~ 0, 
                                 FollowUpTP %in% c(1, 2, 3, 4)~ 1,
                                 FollowUpTP %in% c(5, 6, 7) ~ 2,
                                 FollowUpTP %in% c(8, 9) ~ 3))


# Compare across designs
table(bhl_meta$FollowUpTP, bhl_meta$design)
table(bhl_meta$follow_grp, bhl_meta$design)
table(bhl_meta$follow_grp2, bhl_meta$design)
table(bhl_meta$follow_grp3, bhl_meta$design)
 
#### Run overall analysis: pooled timepoints, pooled + moderators
# The leave one out analyses make this take a few seconds to run.
source("code/analysis_overall.R")

# Run each individual time point group analysis and sensitivity
source("code/analysis_tp0.R")
source("code/analysis_tp1.R")
source("code/analysis_tp2.R")
source("code/analysis_tp3.R")
source("code/timepoint_regroup1/analysis_tp1.1.R")
source("code/timepoint_regroup1/analysis_tp2.1.R")
source("code/timepoint_regroup2/analysis_tp1.2.R")
source("code/timepoint_regroup2/analysis_tp2.2.R")

source("code/analysis_sensitivity.R")
