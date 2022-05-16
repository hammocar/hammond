#' Fits variance component models using the VCA package and returns the precision profiles for each intended measureand level.
#'
#' @param Observed A character string containing the name of the response variable for the precision model.
#' @param VCvars A vector of character variable names for which variance components are to be fit. Levels should be specified assuming no nesting.
#' @param LevelVar A character string containing the name of the variable for the intended measureand levels. A random effects model will be fit for each unique value of LevelVar. 
#' @param Data Specifies the dataset to be used for the analysis.
#' @param Method Either "reml" or "anova" and is an argument passed to the fitVCA function.
#' @return Returns: a list of full AOV tables (one for each measureand level), vectors containing estimated repeatability and reproducibility CV, SD, and variance, and vectors containing the mean observed measureand and intended measureand levels.
#' @examples
#' N_run = 6
#' N_reps_per_level_per_run = 30
#' N_Levels = 5
#' N_reps_per_run = N_reps_per_level_per_run * N_Levels
#' Run = rep(1:N_run,each=N_reps_per_run)
#' Intended = rep(rep(1:N_Levels,each=N_reps_per_level_per_run),N_run)
#' RunEffects = rep(rnorm(N_run,sd=0.2),each=N_reps_per_run)
#' Obs = RunEffects + rnorm(N_reps_per_run*N_run,mean=Intended,sd=0.3)
#' Data = tibble(Run, Intended, Obs)
#' GetQuantPrecision(Observed="Obs", VCvars="Run", LevelVar="Intended", Data=Data)
#' @export
#' 

GetQuantPrecision = function(Observed,VCvars,LevelVar,Data,Method="reml"){
  #make sure 'Observed' is numeric and all VCvars are factors
  Data = Data %>% as_tibble() %>%
    mutate(across(VCvars,~as.factor(.x)),
           across(Observed,~as.numeric(.x)),
           across(LevelVar,~as.factor(.x)))
  #put together VCA formula
  Formula = as.formula(paste(Observed,paste(VCvars,collapse="+"),sep="~"))
  #create vector of numeric intended measureand levels
  Levels = Data %>% pull(get(LevelVar)) %>% table() %>% names() %>% as.numeric()
  #find the mean observed value for each intended level
  MeanObs = Data %>% 
    group_by(get(LevelVar)) %>%
    summarize("MeanObs"=mean(get(Observed))) %>% 
    pull(MeanObs)
  #initialize list/vectors for function output
  FullAOV=vector(mode="list",length=length(Levels))
  VCAfits=vector(mode="list",length=length(Levels))
  CVrepeat=rep(NA,length=length(Levels))
  CVreprod=rep(NA,length=length(Levels))
  VARrepeat=rep(NA,length=length(Levels))
  VARreprod=rep(NA,length=length(Levels))
  SDrepeat=rep(NA,length=length(Levels))
  SDreprod=rep(NA,length=length(Levels))
  #cycle through each intended level and calculated the desired output
  for(l in Levels){
    #create temporary dataset limited to a single intended level
    TempDat = Data %>% filter(get(LevelVar)==l) %>% as.data.frame()
    #fit VCA model
    Fit=fitVCA(Formula, Data=TempDat, method=Method)
    #save output
    VCAfits[[which(Levels==l)]] = Fit
    FullAOV[[which(Levels==l)]] = Fit$aov.tab
    CVrepeat[which(Levels==l)] = Fit$aov.tab["error","CV[%]"]
    CVreprod[which(Levels==l)] = Fit$aov.tab["total","CV[%]"]
    VARrepeat[which(Levels==l)] = Fit$aov.tab["error","VC"]
    VARreprod[which(Levels==l)] = Fit$aov.tab["total","VC"]
    SDrepeat[which(Levels==l)] = Fit$aov.tab["error","SD"]
    SDreprod[which(Levels==l)] = Fit$aov.tab["total","SD"]
  }
  
  return(list("AOV_tables"=FullAOV,"CVrepeat"=CVrepeat,"CVreprod"=CVreprod,
              "VARrepeat"=VARrepeat,"VARreprod"=VARreprod,
              "SDrepeat"=SDrepeat,"SDreprod"=SDreprod,
              "MeanObs"=MeanObs,"Levels"=Levels,"VCAfit"=VCAfits))
}



