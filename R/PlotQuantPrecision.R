#' Creates a precision profile plot based on modeling output.  
#'
#' @param Model The output of the GetQuantPrecision function.
#' @param PlotFitted Designates whether a fitted precision profile should also be plotted. If TRUE, ReprodFit and RepeatFit must be specified.
#' @param ReprodFit The output of the fit.vfp function containing a precision profile model for the quantitative reproducibility.
#' @param RepeatFit The output of the fit.vfp function containing a precision profile model for the quantitative repeatability.
#' @param PrecisType One of "Both", "Repeat", or "Reprod" and specifies which variance type to plot.
#' @param Yvar One of 'CV', 'VAR', or 'SD'. Specifies the plot's y-axis.
#' @param Xvar Either "Observed" or "Intended". Determines the x-axis.
#' @param logX TRUE or FALSE. Specifies whether the x-axis should be logged or not.
#' @param ShowLegend Specifies whether a legend should be printed with the plot or not.
#' @param Title The title of the output plot.
#' @param xLabel Label for the x-axis. Defaults to the specified Xvar.
#' @param Xticks Values at which x-axis ticks should be plotted. Defaults to match values of Xvar.
#' @return Returns a precision profile plot.
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
#' Mod = GetQuantPrecision(Observed="Obs", VCvars="Run", LevelVar="Intended", Data=Data)
#' PlotQuantPrecision(Model=Mod,logX=F,ShowLegend=T)
#' @export
#' 

PlotQuantPrecision = function(Model,PlotFitted=F,ReprodFit=NULL,RepeatFit=NULL,
                              PrecisType="Both",Yvar="CV",Xvar="Observed",
                              logX=T,ShowLegend=T,Title="Precision Profile",
                              xLabel=NULL,Xticks=NULL, PlotData = "exclude"){
  
  #vector of values to label with tick marks on x-axis
  if(is.null(Xticks)){
    if(Xvar=="Observed"){
      Xticks = round(Model$MeanObs,3)
    }
    if(Xvar=="Intended"){
      Xticks = Model$Levels
    }
    xLabel = Xvar
  }
  
  #build precision profile predicted curves
  if(PlotFitted){
    if(Xvar=="Observed"){NewDat = Model$MeanObs}
    if(Xvar=="Intended"){NewDat = Model$Levels}
    ReprodPredCV = sqrt(predict(ReprodFit,newdata=NewDat)$Fitted)/NewDat*100
    ReprodPredSD = sqrt(predict(ReprodFit,newdata=NewDat)$Fitted)
    ReprodPredVAR = predict(ReprodFit,newdata=NewDat)$Fitted
    RepeatPredCV = sqrt(predict(RepeatFit,newdata=NewDat)$Fitted)/NewDat*100
    RepeatPredSD = sqrt(predict(RepeatFit,newdata=NewDat)$Fitted)
    RepeatPredVAR = predict(RepeatFit,newdata=NewDat)$Fitted
  }else{
    ReprodPredCV = rep(NA,length(Model$MeanObs))
    ReprodPredSD = rep(NA,length(Model$MeanObs))
    ReprodPredVAR = rep(NA,length(Model$MeanObs))
    RepeatPredCV = rep(NA,length(Model$MeanObs))
    RepeatPredSD = rep(NA,length(Model$MeanObs))
    RepeatPredVAR = rep(NA,length(Model$MeanObs))
  }
  
  #build plotting dataframe
  PPPlotData =   tibble("Observed"=rep(Model$MeanObs,4),"Intended"=rep(Model$Levels,4),
                        "PrecType"=rep(c(rep("Repeatability",length(Model$MeanObs)),
                                         rep("Reproducibility",length(Model$MeanObs))),2),
                        "DataType"=c(rep("Mixed Model",length(Model$MeanObs)*2),
                                     rep("Smoothed",length(Model$MeanObs)*2)),
                        "CV"=c(Model$CVrepeat,Model$CVreprod,RepeatPredCV,ReprodPredCV),
                        "VAR"=c(Model$VARrepeat,Model$VARreprod,RepeatPredVAR,ReprodPredVAR),
                        "SD"=c(Model$SDrepeat,Model$SDreprod,RepeatPredSD,ReprodPredSD))
  #filter to only the data you want to plot
  if(PrecisType=="Repeat"){PPPlotData = PPPlotData %>% filter(PrecType=="Repeatability")}
  if(PrecisType=="Reprod"){PPPlotData = PPPlotData %>% filter(PrecType=="Reproducibility")}
  if(!PlotFitted){PPPlotData = PPPlotData %>% filter(DataType=="Observed")}
  
  #create clean y-axis label
  ylab=if(Yvar=="CV"){"% CV"}else{if(Yvar=="VAR"){"Variance"}else{if(Yvar=="SD"){"Standard Deviation"}}}
  #generate the plot
  if(PlotFitted){
    if(logX){
      #make plot w/out fitted precision profile, log x axis
      P = ggplot(PPPlotData, aes(x=log(get(Xvar)),y=get(Yvar),col=PrecType)) + 
        geom_point(aes(shape=DataType)) + geom_line(aes(linetype=DataType)) + 
        labs(x=xLabel,y=ylab,title=Title,col="Precision Type",linetype="Data Type",shape="Data Type") + 
        scale_x_continuous(breaks=log(Xticks), labels=Xticks) 
    }else{
      #make plot w/out fitted precision profile, non-log x axis
      P = ggplot(PPPlotData, aes(x=get(Xvar),y=get(Yvar),col=PrecType)) +
        geom_point(aes(shape=DataType)) + geom_line(aes(linetype=DataType)) + 
        labs(x=xLabel,y=ylab,title=Title,col="Precision Type",linetype="Data Type",shape="Data Type") + 
        scale_x_continuous(breaks=log(Xticks), labels=Xticks) 
    }
  }else{
    if(logX){
      #make plot with fitted precision profile, log x axis
      P = ggplot(PPPlotData, aes(x=log(get(Xvar)),y=get(Yvar),col=PrecType)) + 
        geom_point() + geom_line() + 
        labs(x=xLabel,y=ylab,title=Title,col="Precision Type") + 
        scale_x_continuous(breaks=log(Xticks), labels=Xticks) 
    }else{
      #make plot with fitted precision profile, non-log x axis
      P = ggplot(PPPlotData, aes(x=get(Xvar),y=get(Yvar),col=PrecType)) + 
        geom_point() + geom_line() + 
        labs(x=xLabel,y=ylab,title=Title,col="Precision Type") + 
        scale_x_continuous(breaks=log(Xticks), labels=Xticks) 
    }
  }
  
  #remove legend if desired
  if(!ShowLegend){P = P + theme(legend.position = "none")}
  
  ifelse(PlotData == "exclude",return(P), return(list(P, PPPlotData)))
}





