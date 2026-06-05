#### Overall #########
# first, include all effects
overall_effect <- rma.mv(final_hedges_g, var_hedges_g, 
                         random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                       ~ TP_weeks |outcome_id),
                         struct = "CAR",
                         tdist = T,
                         data = bhl_meta,
                         control=list(iter.max = 1000, rel.tol = 10^-10), 
                         slab = full_id,
                         subset =(mismatch_flag==F))
robust(overall_effect, cluster = covidence_num, clubSandwich = T) 


# get I-sq
W <- diag(1/overall_effect$vi)
X <- model.matrix(overall_effect)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
overall_I2 <- 100*sum(overall_effect$sigma2) / 
  (sum(overall_effect$sigma2) + (overall_effect$k - overall_effect$p)/sum(diag(P)))
overall_I2_comps <- 100*overall_effect$sigma2 / 
  (sum(overall_effect$sigma2) + (overall_effect$k - overall_effect$p)/sum(diag(P)))
# Overall I^2 is 81% : percent of total variance due to heterogeneity.
# 67.6% of this is attributable to variation between studies.

png(filename="images/funnel.png", width=1200, height=1200, units = "px",
    family="sans",
    res = 200)
funnel(overall_effect, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot of Pooled Effects")
dev.off()


ranktest(overall_effect)

# Egger - run the same analysis, but include the standard error as
# a moderator.
egger_multi <- rma.mv(final_hedges_g, var_hedges_g, 
                             random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                           ~ TP_weeks |outcome_id),
                             struct = "CAR",
                      mods = ~se_hedges_g,
                             tdist = T,
                             data = bhl_meta,
                             control=list(iter.max = 1000, rel.tol = 10^-10), 
                             subset =(mismatch_flag==F))
coef_test(egger_multi, vcov = "CR2")

### Check for Influential studies #####
study_list <- unique(bhl_meta$covidence_num)
esize_vec <- c()
p_vec <- c()
i2_vec <- c()
for (i in study_list){
  loo_effect <- rma.mv(final_hedges_g, var_hedges_g, 
                           random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                         ~ TP_weeks |outcome_id),
                           struct = "CAR",
                           tdist = T,
                           data = bhl_meta,
                           control=list(iter.max = 1000, rel.tol = 10^-10), 
                           slab = full_id,
                           subset =(mismatch_flag==F & covidence_num != i))
  loo_robust <- robust(loo_effect, cluster = covidence_num, clubSandwich = T) 
  esize_vec <- c(esize_vec, loo_robust$b)
  p_vec <- c(p_vec, loo_robust$pval)
  
  # get I-sq
  W <- diag(1/loo_effect$vi)
  X <- model.matrix(loo_effect)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  loo_I2 <- 100*sum(loo_effect$sigma2) / 
    (sum(loo_effect$sigma2) + (loo_effect$k - loo_effect$p)/sum(diag(P)))
 i2_vec <- c(i2_vec, loo_I2)
}
min(i2_vec)
max(p_vec)
max(esize_vec)
min(esize_vec)

# hist(esize_vec)
# hist(p_vec)
# hist(i2_vec)
# study_list[which.min(i2_vec)]
# esize_vec[which.min(i2_vec)]
# plot(esize_vec, i2_vec)
#### End check #####











##### Moderator - design #########
## Exclude flag ###
bhl_meta$design <- as.factor(bhl_meta$design)
overall_m_effect <- rma.mv(final_hedges_g, var_hedges_g, 
                         random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                       ~ TP_weeks |outcome_id),
                         struct = "CAR",
                         mods = ~ design,
                         tdist = T,
                         data = bhl_meta,
                         control=list(iter.max = 1000, rel.tol = 10^-10), 
                         subset =(mismatch_flag==F))
robust(overall_m_effect, cluster = covidence_num, clubSandwich = T) 
predict(overall_m_effect, newmods = 1)

# get I-sq
W <- diag(1/overall_m_effect$vi)
X <- model.matrix(overall_m_effect)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
overall_m_I2 <- 100*sum(overall_m_effect$sigma2) / 
  (sum(overall_m_effect$sigma2) + (overall_m_effect$k - overall_effect$p)/sum(diag(P)))
overall_m_I2_comps <- 100*overall_m_effect$sigma2 / 
  (sum(overall_m_effect$sigma2) + (overall_m_effect$k - overall_m_effect$p)/sum(diag(P)))
overall_m_I2
# Overall I^2 is 75% : percent of total variance due to heterogeneity.
# 60% of this is attributable to variation between studies.

png(filename="images/funnel_design.png", width=1200, height=1200, units = "px",
    family="sans",
    res = 200)
funnel(overall_m_effect, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot of Pooled Effects
with Design as Moderator")
dev.off()

# Egger - run the same analysis, but include the standard error as
# a moderator.
egger_m_multi <- rma.mv(final_hedges_g, var_hedges_g, 
                      random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                    ~ TP_weeks |outcome_id),
                      struct = "CAR",
                      mods = ~se_hedges_g+design,
                      tdist = T,
                      data = bhl_meta,
                      control=list(iter.max = 1000, rel.tol = 10^-10), 
                      subset =(mismatch_flag==F))
coef_test(egger_m_multi, vcov = "CR2")


### Check for Influential studies #####
study_list <- unique(bhl_meta$covidence_num)
esize_vec1 <- c()
p_vec1 <- c()

esize_vec2 <- c()
p_vec2 <- c()
i2_vec <- c()
for (i in study_list){
  loo_effect <- rma.mv(final_hedges_g, var_hedges_g, 
                       random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                     ~ TP_weeks |outcome_id),
                       struct = "CAR",
                       mods = ~ design,
                       tdist = T,
                       data = bhl_meta,
                       control=list(iter.max = 1000, rel.tol = 10^-10), 
                       slab = full_id,
                       subset =(mismatch_flag==F & covidence_num != i))
  loo_robust <- robust(loo_effect, cluster = covidence_num, clubSandwich = T) 
  esize_vec1 <- c(esize_vec1, loo_robust$b[1])
  esize_vec2 <- c(esize_vec2, loo_robust$b[2])
  
  p_vec1 <- c(p_vec1, loo_robust$pval[1])
  p_vec2 <- c(p_vec2, loo_robust$pval[2])
  # get I-sq
  W <- diag(1/loo_effect$vi)
  X <- model.matrix(loo_effect)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  loo_I2 <- 100*sum(loo_effect$sigma2) / 
    (sum(loo_effect$sigma2) + (loo_effect$k - loo_effect$p)/sum(diag(P)))
  i2_vec <- c(i2_vec, loo_I2)
}
min(i2_vec)
max(p_vec1)
max(p_vec2)
max(esize_vec1)
min(esize_vec1)
max(esize_vec2)
min(esize_vec2)



##### Moderator - length of intervention #########
## Exclude flag ### 
overall_m_effect2 <- rma.mv(final_hedges_g, var_hedges_g, 
                           random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                         ~ TP_weeks |outcome_id),
                           struct = "CAR",
                           mods = ~ SILength,
                           tdist = T,
                           data = bhl_meta,
                           control=list(iter.max = 1000, rel.tol = 10^-10), 
                           subset =(mismatch_flag==F))
robust(overall_m_effect2, cluster = covidence_num, clubSandwich = T) 

# get I-sq
W <- diag(1/overall_m_effect2$vi)
X <- model.matrix(overall_m_effect2)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
overall_m2_I2 <- 100*sum(overall_m_effect2$sigma2) / 
  (sum(overall_m_effect2$sigma2) + (overall_m_effect2$k - overall_effect$p)/sum(diag(P)))
overall_m2_I2_comps <- 100*overall_m_effect2$sigma2 / 
  (sum(overall_m_effect2$sigma2) + (overall_m_effect2$k - overall_m_effect2$p)/sum(diag(P)))
overall_m2_I2
# Overall: 80.98%. Adding intervention length doesn't explain anything.

png(filename="images/funnel_intlength.png", width=1200, height=1200, units = "px",
    family="sans",
    res = 200)
funnel(overall_m_effect2, xlab = "Hedges' g", 
       digits = 2)
title(main = "Funnel Plot of Pooled Effects
with Intervention Length as Moderator")
dev.off()


# Egger - run the same analysis, but include the standard error as
# a moderator.
egger_m_multi <- rma.mv(final_hedges_g, var_hedges_g, 
                        random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                      ~ TP_weeks |outcome_id),
                        struct = "CAR",
                        mods = ~se_hedges_g+SILength,
                        tdist = T,
                        data = bhl_meta,
                        control=list(iter.max = 1000, rel.tol = 10^-10), 
                        subset =(mismatch_flag==F))
coef_test(egger_m_multi, vcov = "CR2")


### Check for Influential studies #####
study_list <- unique(bhl_meta$covidence_num)
esize_vec1 <- c()
p_vec1 <- c()

esize_vec2 <- c()
p_vec2 <- c()
i2_vec <- c()
for (i in study_list){
  loo_effect <- rma.mv(final_hedges_g, var_hedges_g, 
                       random = list(~ 1 | covidence_num/intervention_id/outcome_id,
                                     ~ TP_weeks |outcome_id),
                       struct = "CAR",
                       mods = ~ SILength,
                       tdist = T,
                       data = bhl_meta,
                       control=list(iter.max = 1000, rel.tol = 10^-10), 
                       slab = full_id,
                       subset =(mismatch_flag==F & covidence_num != i))
  loo_robust <- robust(loo_effect, cluster = covidence_num, clubSandwich = T) 
  esize_vec1 <- c(esize_vec1, loo_robust$b[1])
  esize_vec2 <- c(esize_vec2, loo_robust$b[2])
  
  p_vec1 <- c(p_vec1, loo_robust$pval[1])
  p_vec2 <- c(p_vec2, loo_robust$pval[2])
  # get I-sq
  W <- diag(1/loo_effect$vi)
  X <- model.matrix(loo_effect)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  loo_I2 <- 100*sum(loo_effect$sigma2) / 
    (sum(loo_effect$sigma2) + (loo_effect$k - loo_effect$p)/sum(diag(P)))
  i2_vec <- c(i2_vec, loo_I2)
}
min(i2_vec)
max(p_vec1)
max(p_vec2)
min(p_vec2)
max(esize_vec1)
min(esize_vec1)
max(esize_vec2)
min(esize_vec2)









####### Forest plot for Pooled ###########################################################
png(filename="images/forest_plot.png", width=1200, height=3200, units = "px",
    family="sans",
    res = 300)
 forest(overall_effect,
        order = "fit", 
        shade = T, cex = .5)
dev.off()






