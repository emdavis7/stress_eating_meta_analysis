##### Analysis at timepoint 2 ######
bhl_meta %>% filter((mismatch_flag==F & follow_grp2 == 2)) %>% nrow()

second <- rma.mv(final_hedges_g, var_hedges_g, 
                random = list(~ 1 | covidence_num/intervention_id/outcome_id, 
                              ~ TP_weeks |outcome_id),
                struct = "CAR",
                tdist = T,
                data = bhl_meta,
                control=list(iter.max = 2000, rel.tol = 10^-10), 
                subset =(mismatch_flag==F & follow_grp2 == 2),
                slab = full_id)
second_rob <- robust(second, cluster = covidence_num, clubSandwich = T)
second_rob
# I2, when restricted to the second timepoint:
W <- diag(1/second$vi)
X <- model.matrix(second)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
second_I2 <- 100*sum(second$sigma2) / 
  (sum(second$sigma2) + (second$k - second$p)/sum(diag(P)))
second_I2_comp <- 100*second$sigma2 / 
  (sum(second$sigma2) + (second$k - second$p)/sum(diag(P)))
second_I2

#### Pub bias check ####
png(filename="images/funnel_sens_tp2_1.png", width=1200, height=1200, units = "px",
    family="sans",
    res = 200)
funnel(second, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot of Short-Term Effects 
(1 week to 1 month post-intervention)")
dev.off()


ranktest(second)
egger_second <- rma.mv(final_hedges_g, var_hedges_g, 
                      random = list(~ 1 | covidence_num/outcome_id),
                      mods = ~se_hedges_g,
                      tdist = T,
                      data = bhl_meta,
                      control=list(iter.max = 1000, rel.tol = 10^-10), 
                      subset =(mismatch_flag==F & follow_grp2 ==2))
coef_test(egger_second, vcov = "CR2")

####### Forest plot ###########################################################
png(filename="images/forest_plot_shortterm.png", width=1500, height=2000, units = "px",
    family="sans",
    res = 200)
forest(second_rob,
       order = "fit",
       shade = T)
title(main = "Forest Plot of Short-Term Effects (1 week to 1 month post-intervention)")
dev.off()
