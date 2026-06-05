##### Analysis at timepoint 2 ######
bhl_meta %>% filter((mismatch_flag==F & follow_grp3 == 2)) %>% nrow()

second <- rma.mv(final_hedges_g, var_hedges_g, 
                random = list(~ 1 | covidence_num/intervention_id),
                tdist = T,
                data = bhl_meta,
                control=list(iter.max = 2000, rel.tol = 10^-10), 
                subset =(mismatch_flag==F & follow_grp == 2),
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
png(filename="images/funnel_sens_tp2_2.png", width=1200, height=1200, units = "px",
    family="sans",
    res = 200)
funnel(second, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot of Medium-Term Effects
(2 to 4 Months Post-Intervention)")
dev.off()

ranktest(second)
egger_second <- rma.mv(final_hedges_g, var_hedges_g, 
                      random = list(~ 1 | covidence_num/outcome_id),
                      mods = ~se_hedges_g,
                      tdist = T,
                      data = bhl_meta,
                      control=list(iter.max = 1000, rel.tol = 10^-10), 
                      subset =(mismatch_flag==F & follow_grp==2))
coef_test(egger_second, vcov = "CR2")


