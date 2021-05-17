
#' Monte carlo runs n_simulations and calls the callback function each time with the ... optional args
#'
#' @param n_simulations Number of simulations to run
#' @param callback Function to call each simulation

#' @export
monte_carlo <- function(n_simulations, callback, ...){
  simulations <- 1:n_simulations

  sapply(1:n_simulations, function(x){
    callback(...)
  })
}




#' Computes the Clopper/Pearon exact ci for a binomial success probability for x successes out of n trials with confidence coefficient conflev
#'
#' @param x Number of successes
#' @param n Number of trials
#' @param conflev Confidence level
#'
#' @return a list with lower bound, point estimate, and upper bound.
#' @export
exactci <- function(x,n,conflev = 0.95){
  alpha <- (1 - conflev)
  if (x == 0) {
    ll <- 0
    ul <- 1 - (alpha/2)^(1/n)
  }
  else if (x == n) {
    ll <- (alpha/2)^(1/n)
    ul <- 1
  }
  else {
    ll <- 1/(1 + (n - x + 1) / (x * qf(alpha/2, 2 * x, 2 * (n-x+1))))
    ul <- 1/(1 + (n - x) / ((x + 1) * qf(1-alpha/2, 2 * (x+1), 2 *
                                           (n-x))))
  }
  list(lower = ll,
       point_est = x/n,
       upper = ul)
}


#' Computes the Agresti-Coull `add 4' CI for x successes out of n trials with confidence coefficient conflev. Adds 2 successes and 4 trials.
#'
#' @param x Number of successes
#' @param n Number of trials
#' @param conflev Confidence level
#'
#' @return a list with lower bound, point estimate, and upper bound.
#' @export
add4ci <- function(x,n,conflev = 0.95){
  ptilde = (x+2)/(n+4)
  z = abs(qnorm((1-conflev)/2))
  stderr = sqrt(ptilde * (1-ptilde)/(n+4))
  ul = ptilde + z * stderr
  ll = ptilde - z * stderr
  if(ll < 0) ll = 0
  if(ul > 1) ul = 1
  list(lower = ll,
       point_est = x/n,
       upper = ul)
}



#' Square-and-Add or MOVER (Method of Variance Estimate Recovery) Confidence
#' Interval
#' Computes a (1 - alpha) * 100% confidence interval for the difference in two
#' independent proportions
#' Based on exact binomial (Clopper-Pearson) confidence intervals for single proportions
#' References:
#' (1) Robert G. Newcombe. Confidence Intervals for Proportions and Related Measures of Effect Size. CRC Press. 2013
#' (2) Zou, Huang, and Zhang. A note on confiddence interval estimation for a linear function of binomial proportions. Computational Statistics and Data Analysis. 2009
#'
#' @param x1 Number of successes in group 1
#' @param N1 Number of trials in group 1
#' @param x2 Number of successes in group 2
#' @param N2 Number of trials in group 2
#' @param confidence_level Confidence level
#'
#' @return a list with lower bound, point estimate, and upper bound.
#' @export
exact_ind_mover_ci <- function(x1, N1, x2, N2, confidence_level = 0.95) {
  # Input
  #      x1: number of successes in group 1
  #      N1: sample size of group 1
  #      x2: number of successes in group 2
  #      N2: sample size of group 2
  #      confidence_level: the confidence level of the confidence interval
  # Output
  #      vector with three elements:
  #        (1) lower confidence interval bound
  #        (2) point estimate of proportion difference
  #        (3) upper confidence interval bound
  if (any(is.na(c(x1, x2, N1, N2)))) {
    return(c(NA_real_, NA_real_, NA_real_))
  }

  alpha <- 1 - confidence_level
  phat1 <- x1 / N1
  phat2 <- x2 / N2
  p1_p2_diff <- phat1 - phat2

  # Use beta distribution quantiles to compute exact binomial confidence
  # interval of the proportion of successes in group 1
  lower_ci1 <- qbeta(alpha / 2, x1, N1 - x1 + 1)
  upper_ci1 <- qbeta(1 - (alpha / 2), x1 + 1, N1 - x1)

  # Use beta distribution quantiles to compute exact binomial confidence
  # interval of the proportion of successes in group 2
  lower_ci2 <- qbeta(alpha / 2, x2, N2 - x2 + 1)
  upper_ci2 <- qbeta(1 - (alpha / 2), x2 + 1, N2 - x2)

  # Use square and add approach to compute MOVER confidence interval bounds
  lower_mover <- p1_p2_diff - sqrt((phat1 - lower_ci1) ^ 2 +
                                     (upper_ci2 - phat2) ^ 2)
  upper_mover <- p1_p2_diff + sqrt((upper_ci1 - phat1) ^ 2 +
                                     (phat2 - lower_ci2) ^ 2)

  list(lower = lower_mover,
       point_est = p1_p2_diff,
       upper = upper_mover)
}




#'This function computes a confidence interval for a proportion. It is based on inverting the large-sample normal score test for the proportion.
#'
#' @param x Number of successes
#' @param n Number of trials
#' @param conflev Confidence level
#'
#' @return a list with lower bound, point estimate, and upper bound.
#' @export
scoreci <- function(x,n,conflev = 0.95){
  zalpha <- abs(qnorm((1-conflev)/2))
  phat <- x/n
  bound <- (zalpha*((phat*(1-phat)+(zalpha**2)/(4*n))/n)**(1/2))/
    (1+(zalpha**2)/n)
  midpnt <- (phat+(zalpha**2)/(2*n))/(1+(zalpha**2)/n)

  uplim <- round(midpnt + bound,digits=4)
  lowlim <- round(midpnt - bound,digits=4)

list(lower = lowlim,
     point_est = x/n,
     upper = uplim)
}



sim_CI<- function(method, x, n, nsims, conflev){
    x_sim <- rbinom(nsims,n,x/n)
    if(method == "Clopper-Pearson"){
    ci_alt<-  exactci(x = x_sim, n = n, conflev = conflev)
    }
    if(method == "Score"){
      ci_alt<-scoreci(x = x_sim, n = n, conflev = conflev)
    }
    if(method == "Agresti-Coull"){
      ci_alt<-add4ci(x_sim, n, conflev = conflev)
    }
    ci_alt
}
powers<-1:10000
for (i in 1:10000){
x1<-sim_CI(method = "Clopper-Pearson", x = .975*405, n = 405, nsims = 10000, conflev = 0.90)
powers[i]<-mean(x1$lower > .95)}
mean(powers)


runTest <- function(n, sens, p){
  x <- rbinom(n = 1, size = n, prob = sens)
  tmp <- binom.test(x, n, p = p, alt = "greater")
  tmp$conf.int[1] > p
}

get_power <- function(n_ci_simulations, n, sens, p) {

  mean(monte_carlo(n_simulations =  n_ci_simulations, callback = runTest,
            n = n,
            sens = sens,
            p = p), na.rm = TRUE)
}

get_many_powers <- function(n_power_simulations, n, sens, p, n_ci_simulations) {

  monte_carlo(n_simulations =  n_power_simulations,
              callback = get_power,
              n_ci_simulations = n_ci_simulations,
              n = n,
              sens = sens,
              p = p)
}

get_many_powers(n_power_simulations = 1000,
                n = 405,
                sens = .975,
                p = .95,
                n_ci_simulations = 1000)



calculate_power <- function(method, x, n, null_value, alternative = "two-sided", alpha = .05, nsims = 1000){
  conflev <- 1 - alpha
  simCIs<-sim_CI(method = method, x = x, n = n, nsims = nsims, conflev = conflev)
  if(alternative == "greater"){
  power<-mean(simCIs$lower > null_value, na.rm = T)
  }
  if(alternative == "less"){
  power<-mean(simCIs$upper < null_value, na.rm = T)
  }
  if(alternative == "two-sided"){
  power<-mean(simCIs$lower > null) + mean(simCIs$upper < null_value, na.rm = T)
  }
  power
}


#'This function performs a power analysis using the given method
#'
#' @param method Method of calculating confidence interval
#' @param effect_sizes Effect size(s) to examine
#' @param sample_sizes Sample sizes to examine
#'
#'
#' @return a list with lower bound, point estimate, and upper bound.
#' @export
power_calc_one_sample <- function(method, sample_p, sample_sizes, alternative, null_value, alpha = 0.05, nsims = 10000, nsims_power = 1000) {

variables<-
  purrr::cross_df(
    list(
      sample_p = sample_p,
      sample_sizes = sample_sizes,
      null_value = null_value,
      alternative = alternative
    )
  )
variables<- as.data.frame(variables)
# Maps the power simulation through every variable specified in "variables" above.


power_sims<- unlist(purrr::map(1:nrow(variables), ~
                                 mean(monte_carlo(n_simulations = nsims_power,
                                                  callback = calculate_power,
                                                  method = method,
                                                  x = variables$sample_p[.x]*variables$sample_sizes[.x],
                                                  n = variables$sample_sizes[.x],
                                                  alternative = alternative,
                                                  null_value = variables$null_value[.x],
                                                  alpha = alpha,
                                                  nsims = nsims), na.rm = TRUE)))


power_df <- variables %>% mutate(power = power_sims)
power_df
}





power_data<-power_calc_one_sample(method = "Score",
                      sample_p = .975,
                      sample_sizes = seq(300, 350, 5),
                      alternative = "greater",
                      null_value = .95,
                      alpha = 0.1,
                      nsims = 1000,
                      nsims_power = 10000)



ggplot(power_data, aes(x = sample_sizes, y = power))+
  geom_point()+
  geom_line()+
  mytheme()



powrange <- seq(0.4, .9, .1)

n <- sapply(powrange, function(i) power.prop.test(p1=0.975,p2=0.95,power=i,sig.level=0.10, alternative = "two.sided")$n)
plot(n, powrange, type="b",xlab="sample size",ylab="power")

power.prop.test(n = 405, p1 = 0.95, p2 = 0.975, sig.level = 0.05, alternative = "one.sided")

pwr.2p.test(h = .975-.95, n = 405, sig.level = 0.10, alternative = "two.sided")
pwr.p.test(h = .975-.95, n = 405, sig.level = 0.10, alternative = "two.sided")
binom.prop.test(.975, 405, .95)
binom.power(.975, n = 405, p = .95,alpha = 0.05, alternative = "greater", method = "exact" )
binom.power()
#' @import purrr














