#' Estimates qualitative reproducibility by calculating the concordance of each call with its within sample replicate majority call.
#'
#' @param Data Specifies the dataset to be used for the analysis.
#' @param SampleID Character string of the variable name that is used to specify samples.
#' @param Call Character string of the variable name specifying the call (e.g. 0/1 or positive/negative) for each replicate.
#' @return Returns the qualitative reproducibility estimate.
#' @examples
#' N_run = 6
#' N_reps_per_level_per_run = 30
#' N_Levels = 5
#' N_reps_per_run = N_reps_per_level_per_run * N_Levels
#' Run = rep(1:N_run,each=N_reps_per_run)
#' Intended = rep(rep(1:N_Levels,each=N_reps_per_level_per_run),N_run)
#' RunEffects = rep(rnorm(N_run,sd=0.2),each=N_reps_per_run)
#' Obs = RunEffects + rnorm(N_reps_per_run*N_run,mean=Intended,sd=0.3)
#' Noise = rnorm(N_reps_per_run*N_run,mean=Intended,sd=0.1)
#' Calls = ifelse(Obs+Noise>2.3,1,0)
#' Data = tibble(Run, Intended, Obs, Calls)
#' GetQualReprod(Data=Data, SampleID="Intended", Call="Calls")
#' @export

GetQualReprod = function(Data,SampleID,Call){
  #rename columns
  Data = Data %>% mutate("Samp" = get(SampleID))
  
  #calculate qualitative reproducibility 
  QualReprod = Data %>%  
    group_by(Samp) %>%
    summarize(MajorityCall=GetMode(get(Call))[1]) %>%
    right_join(Data,by=c("Samp")) %>%
    mutate("Concordant" = if_else(get(Call)==MajorityCall,1,0)) %>%
    group_by(1) %>%
    summarise("Perc_concordant"=sum(Concordant)/n()) %>%
    pull("Perc_concordant")
  
  return(QualReprod)
}



