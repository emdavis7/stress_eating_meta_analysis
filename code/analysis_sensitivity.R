####### Sensitivity Analyses #####

#### With the flagged study #####
overall_effect_flag <- rma.mv(final_hedges_g, var_hedges_g, 
                              random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                            ~ TP_weeks |outcome_id),
                              struct = "CAR",
                              tdist = T,
                              data = bhl_meta,
                              control=list(iter.max = 1000, rel.tol = 10^-10) )
robust(overall_effect_flag, cluster = covidence_num, clubSandwich = T) 

# get I-sq
W <- diag(1/overall_effect_flag$vi)
X <- model.matrix(overall_effect_flag)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
overall_flag_I2 <- 100*sum(overall_effect_flag$sigma2) / 
  (sum(overall_effect_flag$sigma2) + (overall_effect_flag$k - overall_effect$p)/sum(diag(P)))
overall_flag_I2_comps <- 100*overall_effect_flag$sigma2 / 
  (sum(overall_effect_flag$sigma2) + (overall_effect_flag$k - overall_effect_flag$p)/sum(diag(P)))

overall_flag_I2

###### Remove mid-intervention ########
overall_nomid <- rma.mv(final_hedges_g, var_hedges_g, 
                        random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                      ~ TP_weeks |outcome_id),
                        struct = "CAR",
                        tdist = T,
                        data = bhl_meta,
                        control=list(iter.max = 1000, rel.tol = 10^-10), 
                        slab = full_id,
                        subset =(mismatch_flag==F & follow_grp != 0))
robust(overall_nomid, cluster = covidence_num, clubSandwich = T) 


# get I-sq
W <- diag(1/overall_nomid$vi)
X <- model.matrix(overall_nomid)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
overall_nomid_I2 <- 100*sum(overall_nomid$sigma2) / 
  (sum(overall_nomid$sigma2) + (overall_nomid$k - overall_nomid$p)/sum(diag(P)))
overall_nomid_I2_comps <- 100*overall_nomid$sigma2 / 
  (sum(overall_nomid$sigma2) + (overall_nomid$k - overall_nomid$p)/sum(diag(P)))
# Overall I^2 is 81% : percent of total variance due to heterogeneity.
# 67.6% of this is attributable to variation between studies.
overall_nomid_I2


funnel(overall_nomid, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot (Excluding Mid-Intervention Effects)")
ranktest(overall_nomid)

# Egger - run the same analysis, but include the standard error as
# a moderator.
egger_multi_nomid <- rma.mv(final_hedges_g, var_hedges_g, 
                            random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                          ~ TP_weeks |outcome_id),
                            struct = "CAR",
                            mods = ~se_hedges_g,
                            tdist = T,
                            data = bhl_meta,
                            control=list(iter.max = 1000, rel.tol = 10^-10), 
                            subset =(mismatch_flag==F &  follow_grp != 0))
coef_test(egger_multi_nomid, vcov = "CR2")

# including moderators.
overall_nomid_m <- rma.mv(final_hedges_g, var_hedges_g, 
                          random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                        ~ TP_weeks |outcome_id),
                          mods = ~design,
                          struct = "CAR",
                          tdist = T,
                          data = bhl_meta,
                          control=list(iter.max = 1000, rel.tol = 10^-10), 
                          subset =(mismatch_flag==F & follow_grp != 0))
robust(overall_nomid_m, cluster = covidence_num, clubSandwich = T) 
predict(overall_nomid_m, newmods = 1)


##### Try removing those with .5 as a sensitivity check #####
overall_no5 <- rma.mv(final_hedges_g, var_hedges_g, 
                      random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                    ~ TP_weeks |outcome_id),
                      struct = "CAR",
                      tdist = T,
                      data = bhl_meta,
                      control=list(iter.max = 1000, rel.tol = 10^-10), 
                      subset =(mismatch_flag==F & 
                                 (r_pre_post != .5|is.na(r_pre_post)==T)))
robust(overall_no5, cluster = covidence_num, clubSandwich = T) 

# get I-sq
W <- diag(1/overall_no5$vi)
X <- model.matrix(overall_no5)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
overall_no5_I2 <- 100*sum(overall_no5$sigma2) / 
  (sum(overall_no5$sigma2) + (overall_no5$k - overall_effect$p)/sum(diag(P)))
overall_no5_I2_comps <- 100*overall_no5$sigma2 / 
  (sum(overall_no5$sigma2) + (overall_no5$k - overall_no5$p)/sum(diag(P)))
overall_no5_I2

# including moderators.

overall_no5_m <- rma.mv(final_hedges_g, var_hedges_g, 
                        random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                      ~ TP_weeks |outcome_id),
                        mods = ~design,
                        struct = "CAR",
                        tdist = T,
                        data = bhl_meta,
                        control=list(iter.max = 1000, rel.tol = 10^-10), 
                        subset =(mismatch_flag==F & 
                                   (r_pre_post != .5|is.na(r_pre_post)==T)))
robust(overall_no5_m, cluster = covidence_num, clubSandwich = T) 
predict(overall_no5_m, newmods = 1)

# get I-sq
W <- diag(1/overall_no5$vi)
X <- model.matrix(overall_no5)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
overall_no5_I2 <- 100*sum(overall_no5$sigma2) / 
  (sum(overall_no5$sigma2) + (overall_no5$k - overall_effect$p)/sum(diag(P)))
overall_no5_I2_comps <- 100*overall_no5$sigma2 / 
  (sum(overall_no5$sigma2) + (overall_no5$k - overall_no5$p)/sum(diag(P)))

