##### Analysis at mid-int ######
bhl_meta %>% filter((mismatch_flag==F & follow_grp== 0)) %>% nrow() 


zero <- rma.mv(final_hedges_g, var_hedges_g, 
                random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                              ~ TP_weeks |outcome_id),
                struct = "CAR",
                tdist = T,
                data = bhl_meta,
                control=list(iter.max = 1000, rel.tol = 10^-10), 
                subset =(mismatch_flag==F & follow_grp== 0),
                slab = full_id)
zero_rob <- robust(zero, cluster = covidence_num, clubSandwich = T)
zero_rob

# I2, when restricted to the zero timepoint:
W <- diag(1/zero$vi)
X <- model.matrix(zero)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
zero_I2 <- 100*sum(zero$sigma2) / 
  (sum(zero$sigma2) + (zero$k - zero$p)/sum(diag(P)))
zero_I2_comp <- 100*zero$sigma2 / 
  (sum(zero$sigma2) + (zero$k - zero$p)/sum(diag(P)))
zero_I2
zero_I2_comp

#### Pub bias check #####
png(filename="images/funnel_midint.png", width=1200, height=1200, units = "px",
    family="sans",
    res = 200)
funnel(zero, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot of Effects Mid-Intervention")
dev.off()

ranktest(zero)

egger_zero <- rma.mv(final_hedges_g, var_hedges_g, 
                      random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                    ~ TP_weeks |outcome_id),
                      struct = "CAR",
                      mods = ~se_hedges_g,
                      tdist = T,
                      data = bhl_meta,
                      control=list(iter.max = 1000, rel.tol = 10^-10), 
                      subset =(mismatch_flag==F & follow_grp==0))
coef_test(egger_zero, vcov = "CR2")



####### Forest plot ###########################################################
png(filename="images/forest_plot_zero.png", width=1500, height=2000, units = "px",
    family="sans",
    res = 200)
forest(zero_rob,
       order = "fit",
       shade = T)
title(main = "Forest Plot of Effects Mid-Intervention")
dev.off()

