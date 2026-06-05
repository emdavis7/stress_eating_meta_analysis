##### Analysis at timepoint 1 ######
bhl_meta %>% filter((mismatch_flag==F & follow_grp3 == 1)) %>% nrow()

first <- rma.mv(final_hedges_g, var_hedges_g, 
                random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                              ~ TP_weeks |outcome_id),
                struct = "CAR",
                tdist = T,
                data = bhl_meta,
                control=list(iter.max = 1000, rel.tol = 10^-10), 
                subset =(mismatch_flag==F & follow_grp3 == 1),
                slab = full_id)
first_rob <- robust(first, cluster = covidence_num, clubSandwich = T)
first_rob

# I2, when restricted to the first timepoint:
W <- diag(1/first$vi)
X <- model.matrix(first)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
first_I2 <- 100*sum(first$sigma2) / 
  (sum(first$sigma2) + (first$k - first$p)/sum(diag(P)))
first_I2_comp <- 100*first$sigma2 / 
  (sum(first$sigma2) + (first$k - first$p)/sum(diag(P)))

first_I2

#### Pub bias check #####
png(filename="images/funnel_sens_tp2_1.png", width=1200, height=1200, units = "px",
    family="sans",
    res = 200)
funnel(first, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot of Short-Term Effects
(Up to 1 Month Post-Intervention)")
dev.off()

ranktest(first)
egger_first <- rma.mv(final_hedges_g, var_hedges_g, 
                      random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                    ~ TP_weeks |outcome_id),
                      struct = "CAR",
                      mods = ~se_hedges_g,
                      tdist = T,
                      data = bhl_meta,
                      control=list(iter.max = 1000, rel.tol = 10^-10), 
                      subset =(mismatch_flag==F & follow_grp==1))
coef_test(egger_first, vcov = "CR2")
