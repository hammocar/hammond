#' Fits variance component models using the VCA package and returns the precision profiles for each intended measureand level.
#'
#' @param Observed A character string containing the name of the response variable for the precision model.
#' @param VCvars A vector of character variable names for which variance components are to be fit. Levels should be specified assuming no nesting.
#' @param LevelVar A character string containing the name of the variable for the intended measureand levels. A random effects model will be fit for each unique value of LevelVar. 
#' @param RunID Character string of the variable name specifying the run each replicate comes from.
#' @param SampleID Character string of the variable name that is used to specify samples.
#' @param Call Character string of the variable name specifying the call (e.g. 0/1 or positive/negative) for each replicate. Can be left as NULL if there is no qualitative component.
#' @param Data Specifies the dataset to be used for the analysis.
#' @return Returns: a list consisting of a summarized datavframe with estimated quantitative and qualitative repeatability and reproducibility CV, SD, and variance, as well as a data frame of precision profiles
#' @examples
#' @export
#' 
RnR <- function(Observed, 
                     VCvars, 
                     LevelVar, 
                     RunID,
                     SampleID,
                     Call=NULL, 
                     Data){
  
  ###Fit random effects models
  Out = GetQuantPrecision(Observed, VCvars, LevelVar, Data)
  #store variance component confidence limits
  OutInf = VCAinference(Out$VCAfit)
  RepeatLCLs = rep(NA,length=length(Out$VCAfit))
  RepeatUCLs = rep(NA,length=length(Out$VCAfit))
  ReprodLCLs = rep(NA,length=length(Out$VCAfit))
  ReprodUCLs = rep(NA,length=length(Out$VCAfit))
  for(i in 1:length(Out$VCAfit)){
    ReprodLCLs[i] = OutInf[[i]]$ConfInt$CV$TwoSided["total","LCL"]
    ReprodUCLs[i] = OutInf[[i]]$ConfInt$CV$TwoSided["total","UCL"]
    RepeatLCLs[i] = OutInf[[i]]$ConfInt$CV$TwoSided["error","LCL"]
    RepeatUCLs[i] = OutInf[[i]]$ConfInt$CV$TwoSided["error","UCL"]
  }
  
  ########fit smoothed precision profiles based on random effects output
  ReprodDFs = sapply(1:length(Out$AOV_tables),function(x){Out$AOV_tables[[x]]["total","DF"]})
  RepeatDFs = sapply(1:length(Out$AOV_tables),function(x){Out$AOV_tables[[x]]["error","DF"]})
  PrecProfData = data.frame(LevelVar = Out$Levels,"Mean"=Out$MeanObs,
                                "ReprodDF"=ReprodDFs,"ReprodVar"=Out$VARreprod,
                                "RepeatDF"=RepeatDFs,"RepeatVar"=Out$VARrepeat)
  ReprodProfFit = fit.vfp(PrecProfData,model.no=1:10,col.df="ReprodDF",col.var="ReprodVar")
  RepeatProfFit = fit.vfp(PrecProfData,model.no=1:10,col.df="RepeatDF",col.var="RepeatVar")
  
  ##generate data for gg plots of precision profile results
  PlotData = PlotQuantPrecision(Model=Out,PlotFitted=T,PlotData="include",
                                ReprodFit=ReprodProfFit,RepeatFit=RepeatProfFit)[[2]]
  
  #create some basic summary data for the output
  Summary_data<-Data %>% 
    group_by(get(LevelVar)) %>% 
    summarize(Observed_mean = mean(get(Observed)),
              n = n()) %>%
    rename(LevelVar=`get(LevelVar)`)
  
  #if 'Call' is specified calculate qualitative precision
  if(!is.null(Call)){
    QualReprod = sapply(Out$Levels,function(x){
      GetQualReprod(Data %>% filter(get(LevelVar) == x), 
                    SampleID, Call=Call)})
    QualRepeat = sapply(Out$Levels,function(x){
      GetQualRepeat(Data %>% filter(get(LevelVar) ==x), 
                    SampleID, RunID=RunID, Call=Call)})
    Out_Tab = data.frame("Reproducibility % CV" = round(Out$CVreprod,1),
                         "Reproducibility % CV, LCL" = round(ReprodLCLs,1),
                         "Reproducibility % CV, UCL" = round(ReprodUCLs,1),
                         "Repeatability % CV" = round(Out$CVrepeat,1),
                         "Repeatability % CV, LCL" = round(RepeatLCLs,1),
                         "Repeatability % CV, UCL" = round(RepeatUCLs,1),
                         "Qualitative Repeatability" = round(QualRepeat*100,1),
                         "Qualitative Reproducibility" = round(QualReprod*100,1))
  }else{
    #Precision Analysis Results
    Out_Tab = data.frame("Reproducibility % CV" = round(Out$CVreprod,1),
                         "Reproducibility % CV, LCL" = round(ReprodLCLs,1),
                         "Reproducibility % CV, UCL" = round(ReprodUCLs,1),
                         "Repeatability % CV" = round(Out$CVrepeat,1),
                         "Repeatability % CV, LCL" = round(RepeatLCLs,1),
                         "Repeatability % CV, UCL" = round(RepeatUCLs,1))
  }
  
  #build output and
  Out_List = list(Summary_RnR = cbind(Summary_data,Out_Tab),
                  Precision_Profile = PlotData,
                  Out = Out,
                  PrecProfData = PrecProfData,
                  ReprodProfFit = ReprodProfFit,
                  RepeatProfFit = RepeatProfFit
  )
  
  return(Out_List)
}
