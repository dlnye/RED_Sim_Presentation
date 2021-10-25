
####### PACKAGES ##### 

library(tidyverse)
library(brms)
library(AER)
library(rstanarm)

###### 
fit_brms <- function(data, measvar) {
  f1 <- bf(use_app ~ is_encouraged, family = binomial())
  f2 <- bf(as.formula(paste(measvar, "use_app", sep = "~")), set_rescor(rescor = FALSE))
  IV_brm_a <- brm(f1 + f2, data = data)
  return(IV_brm_a)
}

fit_stan_glm_wald <- function(data, measvar){
  
  itt_zt_sglm <- stan_glm(use_app ~ is_encouraged, data = data)
  itt_zy_sglm <- stan_glm(as.formula(paste(measvar, "is_encouraged", sep = "~")), data = data)
  wald_est <- coef(itt_zy_sglm)["is_encouraged"] / coef(itt_zt_sglm)["is_encouraged"]
  return(wald_est)
  
}



fit_iv_reg <- function(data, measvar, summary_only = TRUE){
  fm <-ivreg(as.formula(paste(measvar, "use_app", sep = "~")), instruments = ~ is_encouraged, data = data)
  if(summary_only){
    Est <- coef(fm)
    ci95 <- confint(fm)
    fm.summary <- cbind(Est, ci95)
    return(fm.summary)
  } else {
    
    return(fm)
    
  }
  
}



fit_2stage_glm <- function(data = d, measvar = "ebehave"){
  
  fit_a <- stan_glm(use_app ~ is_encouraged, data = data)
  data$use_app_hat <- fit_a$fitted
  fit_b <- stan_glm(as.formula(paste(measvar, "use_app_hat", sep = "~")), data = data)
  X_adj <- X <- model.matrix(fit_b)
  X_adj[, "use_app_hat"] <- data$use_app
  n <- nrow(X)
  p <- ncol(X)
  RMSE1 <- sqrt( sum( data[,measvar] - X %*% coef(fit_b)^2 ) / (n-p) )
  RMSE2 <- sqrt( sum( data[,measvar] - X_adj %*% coef(fit_b)^2 ) / (n-p) )
  se_adj <- se(fit_b)["use_app_hat"] * RMSE1 / RMSE2
  est <- coef(fit_b)[2]
  ci95 <- est + qnorm(c(lwr95ci=0.025, upr95ci=0.975))*se_adj
  out <- c(est, ci95)
  
  return(out)
  
}



####### RED SIM FUNCTION ##########
sim_red <- function(sample_size = 50, 
                    esave_est_app = 0.05,
                    esave_est_noapp = 0,
                    esave_sd = 0.05,
                    ebehave_rate_app = 0.15,
                    ebehave_rate_noapp = .1,
                    ebehave_sd = 0.05,
                    p_adopt_encouraged = 0.5,
                    p_adopt_noencourage = 0.1,
                    mvar = "esave",
                    output = "data") {
  
  encouraged <- c(0,1)
  is_encouraged <- sample(encouraged, size = sample_size, replace = TRUE)  # random assignment of encouragment
  
  # vectors for probabilities of app adoption given encouragement
  adopts <- c(0,1)
  adopt_probs_encouraged <- c(1-p_adopt_encouraged, p_adopt_encouraged)
  adopt_probs_noencouraged <- c(1-p_adopt_noencourage, p_adopt_noencourage)
  
  # simulate app adoption with different probs based on encouragement status
  use_app <- sapply(is_encouraged, function(x) 
    ifelse(x==1, 
           sample(adopts, size = 1, prob = adopt_probs_encouraged),
           sample(adopts, size = 1, prob = adopt_probs_noencouraged)))
  
  # vectors for probabilities of adopting energy-saving behaviors given use of SW app
  #ebehave_opts = c(0,1)
  #ebehave_probs_app <- c(1-ebehave_rate_app, ebehave_rate_app)
  #ebehave_probs_noapp <- c(1-ebehave_rate_noapp, ebehave_rate_noapp)
  
  # simulate energy behavior adoption with different probs based on encouragement status
  ebehave <- sapply(use_app, function(x) 
      ifelse(x==1, 
             rnorm(sample_size, 
                   mean = ebehave_rate_app,
                   sd = ebehave_sd ),
             rnorm(sample_size, 
                   mean = ebehave_rate_noapp,
                   sd = ebehave_sd)
      ))
    # ifelse(x==1, 
    #        sample(ebehave_opts, size = 1, prob = ebehave_probs_app),
    #        sample(ebehave_opts, size = 1, prob = ebehave_probs_noapp)))
    
  
  esave <- sapply(use_app, function(x)
    ifelse(x==1, 
           rnorm(sample_size, 
                 mean = esave_est_app,
                 sd = esave_sd ),
           rnorm(sample_size, 
                 mean = esave_est_noapp,
                 sd = esave_sd)
    )
  )
  
  
  
  df <- data.frame(id = 1:sample_size, 
                   is_encouraged = is_encouraged,
                   use_app = use_app,
                   ebehave = ebehave,
                   esave = esave)
  
  switch(ouput = output, 
         data = {
           return(df)
         },
         
         wald = {
           fit_stan_glm_wald(data = df, measvar = mvar)
         },
         
         iv_brms = {
           fit_brms(data = df, measvar = mvar)
         },
         
         iv_reg_ols = {
           fit_iv_reg(data = df, measvar = mvar) 
         },
         
         iv_reg_stan = {
           
           fit_2stage_glm(df, measvar = mvar)
         }
  )
}

#### RUN SIM GRID ####
run_sim_grid <-


###### RUN SIM FUCTION ########
run_sim <- function(reps = 100, 
                    ctrl =  list(sample_size = 50, 
                                 esave_est_app = 0.05,
                                 esave_est_noapp = 0,
                                 esave_sd = 0.05,
                                 ebehave_rate_app = 0.15,
                                 ebehave_rate_noapp = 0.1,
                                 ebehave_sd <- 0.05,
                                 p_adopt_encouraged = 0.5,
                                 p_adopt_noencourage = 0.1,
                                 mvar = "esave",
                                 output = "iv_reg_ols"))  {
  
  n_opts <- map_int(ctrl, length)
  if(sum(n_opts>1) > 1) stop("Function can only handle multiple options for one parameter at a time.")
  
  if(all(n_opts == 1)){
    
    out <- vector(mode = "list", length = reps)
    for(j in seq_along(out)){
      out[[j]] <- do.call(sim_red, args = as.list(ctrl))
    } 
    return(out)
    
  }else {
    
    (iter_var <- names(ctrl)[n_opts > 1])
    cat("Iterating on these levels of ",iter_var,":\n", sep = "")
    (iter_var_levels <- ctrl[iter_var][[1]])
    (out <- vector(mode = "list", length = length(iter_var_levels)))
    names(out) <- paste0(iter_var, ctrl[iter_var][[1]])
    cat(iter_var_levels, "\n")
    
    for(i in seq_along(iter_var_levels)){
      
      iv_level <- iter_var_levels[[i]]
      ctrl[iter_var][[1]] <- iv_level
      
      cat("working on ", iter_var, " iteration (",i," of ", length(out),")...\n", sep = "" )
      
      samp_out <- vector(mode = "list", length = reps)
      for(j in seq_along(samp_out)){
        samp_out[[j]] <- do.call(sim_red, args = as.list(ctrl))
        
      }
      out[[i]] <- samp_out
      
    }
   
    
    
    #n_opts <- purrr::map_int(ctrl, length)
    #(iter_var <- names(ctrl)[n_opts > 1])
    
    ctrl[iter_var][[1]] <- iter_var_levels
    (mvar <- ctrl$mvar[[1]])
    if(mvar == "esave") {
      true_val <-  ctrl$esave_est_app
    } else {
      true_val <- ctrl$ebehave_rate_app
    }
    
    # extract the point estimates for each iteration
    sim_est <- purrr::map(out, ~purrr::map_dbl(., ~.[[2]])) %>% 
      bind_cols() %>%
      gather(key = !!iter_var, value = !!mvar)%>%
      mutate(run = 1:n())
    
    # does the CI contain the true est behind the simulation?             
    trueval_inCI_est <- vector(mode = "list", length = length(out))
    names(trueval_inCI_est) <- names(out)
    if(length(true_val) == 1) true_val <- rep(true_val, length(out))
    for(i in seq_along(out)) {
      sim <- out[[i]]
      tval <- rep(true_val[[i]])
      x <- map_lgl(sim, ~tval >= .["use_app",2] & tval <= .["use_app", 3])
      trueval_inCI_est[[i]] <- x
    }  
    trueval_inCI_est <- trueval_inCI_est %>%
      bind_cols %>%
      gather(key = !!iter_var, value = "ci_contains_trueval")%>%
      mutate(run = 1:n())
    
    # does lower bound of CI exclude 0?
    alpha_est <- purrr::map(out, ~purrr::map_lgl(., ~.["use_app",2] > 0)) %>%
      bind_cols %>%
      gather(key = !!iter_var, value = "ci_excludes_zero") %>%
      mutate(run = 1:n())
    
    
    # what is the absolute width of the confidence interval?
    ci95width <- purrr::map(out, ~purrr::map_dbl(., ~.["use_app",3] - .["use_app", 2]))%>%
      bind_cols() %>%
      gather(key = !!iter_var, value = "ci95width") %>%
      mutate(run = 1:n())
    
    out_summary <- reduce(list(sim_est, trueval_inCI_est, alpha_est, ci95width), left_join, by = c("run", iter_var))
    
    out_summary_tbl <- out_summary%>%
      group_by_at(iter_var) %>%
      summarize(n_reps = n(), 
                detect_est_rate = mean(ci_contains_trueval, na.rm = TRUE),
                exclude_zero_rate = mean(ci_excludes_zero, na.rm = TRUE),
                ciwidth_mean = mean(ci95width, na.rm = TRUE),
                ciwidth_lwr = quantile(ci95width, 0.025, na.rm = TRUE),
                ciwidth_upr = quantile(ci95width, 0.975, na.rm = TRUE)
      )
    
    boot_ci_lwr <- out_summary %>% group_by_at(iter_var) %>% summarize_at(mvar,  ~quantile(.,0.025, na.rm = TRUE))
    boot_ci_upr <- out_summary %>% group_by_at(iter_var) %>% summarize_at(mvar,  ~quantile(.,0.975, na.rm = TRUE))  
    boot_ci <- reduce(list(boot_ci_lwr, boot_ci_upr), left_join, by = iter_var, suffix = c("_lwr_ci95", "_upr_ci95"))
    
    (out_summary_tbl <- out_summary_tbl %>% left_join(boot_ci, by = iter_var) %>%
        mutate_at(iter_var, ~str_remove(., !!iter_var)))
    
    
    (ctrl_df <- bind_rows(as.list(ctrl)) %>%
        mutate_at(iter_var, as.character))
    
    
    out_summary_tbl <- left_join(ctrl_df, out_summary_tbl, by = iter_var) %>%
      mutate_at(iter_var, as.numeric)
    
    
    all_results <- list(output_raw = out, 
                        output_data_frame = out_summary, 
                        summary_table = out_summary_tbl)
    
    
    
  }
  
  return(all_results)
  #return(out)
}
