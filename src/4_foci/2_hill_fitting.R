library(tidyverse)

# function for fitting hill plot ------------------------------------------
fitting_func <- function(X, Y, n_predict = 2, change_predict = 0.1){
  require(optimx)
  require(data.table)
  require(tidyverse)
  
  X <- X
  Y <- Y
  
  fun_estimate <- function (Y, X) { # Y = gene expression data of each dose or time-point
    resid <- function (param) {
      yhat <- param[1] + (param[2])*X^param[3] / (param[4] + X^param[3]) # Hill equation
      sum((Y-yhat)^2)	# return Residual sum of squares (RSS)
    }
    
    #omptimize
    initialp <- c(min(Y), max(Y), n_predict, change_predict)  #initial parameteres of parameter, 3rd argument = n, 4th argument = where the change is expected to happen
    RESULT <- optimx(par = initialp, fn = resid, control = list(all.methods = TRUE)) #you can constrain range of parameteres lower = or upper =
    RESULT %>% data.table(keep.rownames=TRUE) %>% 
      dplyr::filter(convcode==0) %>% # filter(min(p1,p2,p3,p4)>0) %>% fr(p4<12) filtering of the results if you want
      dplyr::filter(value == min(value)) -> RESULT_OPT
    return(RESULT_OPT)
  }
  
  # make hill function
  hillfunc <-function(X, param){
    param$p1 + (param$p2 - param$p1)*X^param$p3 / (param$p4^param$p3 + X^param$p3)
  } 
  
  # perform parameter estimation
  res <- fun_estimate(Y, X)
  
  fitted_df <- tibble(X = X, Y = Y, fitted_Y = hillfunc(X, res))
  
  # extract parameters and create df
  res_param <- tibble(param = c("n", "ka", "kd"), value = c(res$p3, res$p4, res$p4^res$p3))
  
  list(parameters = res,
       params_short = res_param,
       fitted = fitted_df) %>% return()
  
}


hillfunc <-function(X){
  param$p1 + param$p2*X^param$p3 / (param$p4 + X^param$p3)
}

# fitting to data -------------------------------------------
dose_dependent <- read_csv("data/foci/3_median_fitting.csv")

a <- fitting_func(X = c(0, 0.01, 0.1, 1, 10),
                  Y = dose_dependent %>%
                    filter(time == 20) %>%
                    select(-time) %>% as.numeric(),
                  change_predict = 0.1, n_predict = 1)


# plots -------------------------------------------------------------------
param <- a$parameters

# plot fitted and data points
ggplot(a$fitted, aes(x = X, y = Y)) + 
    stat_function(fun = hillfunc, size=3, color= "black") + 
    scale_x_continuous(trans = 'log10')
  
# n
param$p3
# Ka
param$p4
# RSS
sum((a$fitted$fitted_Y-a$fitted$Y)^2)




