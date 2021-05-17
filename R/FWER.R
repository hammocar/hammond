sens_spec<- data.frame(
  syndrome = c("T21", "T18", "T13", "MX", "22q11.2", "Aneuploidy Specificty", "22a Specificity", "No-call"),
  n = c(75, 10, 9, 44, 10, 500, 500, 648),
  prob = c(0.9925, 0.988, 0.999, 0.947, 0.92, 0.997, 0.997, 0.235),
  requirements = c(0.99, 0.95, 0.90, 0.90, 0.90, 0.999, 0.998, 0.04))

nsims<-10000

n_rejections<- vector(mode = "numeric", length = nsims)
alpha<- 0.05

for(i in 1:nsims){
  n_rejection<-0
  for(j in sens_spec$syndrome){
    random_sample <- rbinom(1,
                            size = sens_spec[sens_spec$syndrome == j, 'n'],
                            prob = sens_spec[sens_spec$syndrome == j, 'prob'])
    cp<-exactci(x = random_sample,
                n = sens_spec[sens_spec$syndrome == j, 'n'],
                conflev = 1-(alpha*2))

    if(j != "No-call"){
      if(cp$upper < sens_spec[sens_spec$syndrome == j, 'requirements']){
        n_rejection<-n_rejection + 1
      }
    }else{
      if(cp$lower > sens_spec[sens_spec$syndrome == j, 'requirements']){
        n_rejection <- n_rejection + 1
      }
    }
  }
  n_rejections[i]<-n_rejection
}

mean(n_rejections > 1)






##############################################################################################
##############################################################################################


reject_at_i <- function(observations, i, p1, p2){
  p1_sample <- rbinom(observations, 1, p1)
  p2_sample <- rbinom(observations, 1, p2)
  ( prop.test(c(sum(p1_sample[1:i]),sum(p2_sample[1:i])), c(i,i))$p.value ) < alpha
}


p_value_at_i <- function(observations, i, p1, p2){

  conversions_1 <- rbinom(observations, 1, p1)
  conversions_2 <- rbinom(observations, 1, p2)

  prop.test(c(sum(conversions_1[1:i]),sum(conversions_2[1:i])), c(i,i))$p.value
}


FWER <- function(ntrials, n_tests_per_trial, n_obs, p1, p2) {


  rejects <- 0

  for(i in 1:n_trials){
    # run the sim
    rejected.H0 <- monte_carlo(n_tests_per_trial,
                               callback=reject_at_i,
                               observations=n_obs,
                               i=n_obs,
                               p1 = p1,
                               p2 = p2
    )
    if(!is.na(table(rejected.H0)[2])) {
      rejects <- rejects + 1
    }
  }

  # Calculate FWER
  rejects/n_trials

}

FWER(1000, 4, 9720, .25,.25)



set.seed(1)



FWER_power <- function(ntrials, n_tests_per_trial, n_obs, p1, p2, mde) {

  res_bf <- c()
  res_holm <- c()

  for(i in 1:n_trials){
    null_true <- rbinom(1,1,prob = 0.5)
    #randomly setting the effect size to the minimum detectable effect in about half the cases.
    effect <- mde * null_true
    p2 <- (1+effect)*p1

    # run n_tests_per_trial
    p_values <- monte_carlo(n_tests_per_trial,
                            callback=p_value_at_i,
                            observations=n_obs,
                            i=n_obs,
                            p1=p1,
                            p2=p2
    )
    # Bonferroni: adjust the p-values and reject/accept
    reject_bf <- p.adjust(p_values, "bonferroni") <= alpha
    for(r in reject_bf){
      res_bf <- rbind(res_bf, c(r, null_true))
    }

    # Holm: adjust the p-values and reject/accept
    reject_holm <- p.adjust(sort(p_values), "holm") <= alpha
    for(r in reject_holm){
      res_holm <- rbind(res_holm, c(r, null_true))
    }

  }

  # the rows of the table represent the test result
  # while the columns represent the null truth
  table_bf <- table(Test=res_bf[,1], Null=res_bf[,2])
  table_holm <- table(Test=res_holm[,1], Null=res_holm[,2])

  # False positive rate
  fpr_bf <- table_bf['1','0']/sum(table_bf[,'0'])
  fpr_holm <- table_holm['1','0']/sum(table_holm[,'0'])

  print(paste0("FPR Bonferroni: ", round(fpr_bf,3), " FPR Holm: ", round(fpr_holm,3)))
  ## [1] "FPR Bonferroni: 0.029 FPR Holm: 0.029"
  # Power
  power_bf <- table_bf['1','1']/sum(table_bf[,'1'])
  power_holm <- table_holm['1','1']/sum(table_holm[,'1'])

  print(paste0("Power Bonferroni: ", round(power_bf,3), " Power Holm: ", round(power_holm,3)))
  ## [1] "Power Bonferroni: 0.713 Power Holm: 0.774"
  # Comparing the Power of Holm vs. Bonferroni
  print(paste0("Power Holm/Power Bonferroni: ", round(power_holm/power_bf,3)))
  ## [1] "Power Holm/Power Bonferroni: 1.084"


}
FWER_power(1000, 4, 9720, .25,.25, .05)

