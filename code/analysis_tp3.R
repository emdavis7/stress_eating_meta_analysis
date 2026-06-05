##### Analysis at timepoint 3 ######
bhl_meta %>% filter((mismatch_flag==F & follow_grp == 3)) %>% nrow() 

# All only have a single intervention.
third <- rma.mv(final_hedges_g, var_hedges_g, 
                    random = list(~ 1 | covidence_num/outcome_id),
                    tdist = T,
                    data = bhl_meta,
                    control=list(iter.max = 1000, rel.tol = 10^-10), 
                    subset =(mismatch_flag==F & follow_grp == 3),
                slab = full_id)
third_rob <- robust(third, cluster = covidence_num, clubSandwich = T)
third_rob

# I2, when restricted to the third timepoint:
W <- diag(1/third$vi)
X <- model.matrix(third)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
third_I2 <- 100*sum(third$sigma2) / 
  (sum(third$sigma2) + (third$k - third$p)/sum(diag(P)))
third_I2_comp <- 100*third$sigma2 / 
  (sum(third$sigma2) + (third$k - third$p)/sum(diag(P)))
third_I2

#### Pub bias check #############
png(filename="images/funnel_long.png", width=1200, height=1200, units = "px",
    family="sans",
    res = 200)
funnel(third, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot of Effects at and after 6 Months")
dev.off()

ranktest(third)
egger_third <- rma.mv(final_hedges_g, var_hedges_g, 
                      random = list(~ 1 | covidence_num/outcome_id),
                      mods = ~se_hedges_g,
                      tdist = T,
                      data = bhl_meta,
                      control=list(iter.max = 1000, rel.tol = 10^-10), 
                      subset =(mismatch_flag==F & follow_grp==3))
coef_test(egger_third, vcov = "CR2")


### Sensitivity analysis #####
# Do not do "include one more study" here, because this is just in tp_grp ==1

####### Forest plot ###########################################################
png(filename="images/forest_plot_third.png", width=1500, height=2100, units = "px",
    family="sans",
    res = 200)
forest(third_rob,
       order = "fit",
       shade = T)
title(main = "Forest Plot of Long-Term Effects (6 Months and Longer)")
dev.off()

