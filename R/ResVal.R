#' @title Validation Summary
#' @description  Displays the summary of the validation.
#' @param obs_c vector of observational climate data (rainfall) used for calibrating the model
#' @param mod_c vector of GCM/RCM climate data (rainfall) used for calibrating the model
#' @param obs_v vector of observational climate data (rainfall) used for validating the model
#' @param mod_v vector of GCM/RCM climate data (rainfall) used for validating the model
#' @param mod_fut vector of GCM/RCM future climate data (rainfall) need to be downscaled
#' @details
#'
#'1)	Dry-days correction / Defining threshold values
#'
#'    The relationship between the cumulative frequencies (thresholds) corresponding to the dry days of GCM/RCM data and that of the observational data is defined by a polynomial function given by;
#'
#'threshold_obs = (threshold_mod)^n
#'
#'n = ln(threshold_obs_c) / ln(threshold_mod_c)
#'
#'
#'2)	wet-days correction / Correcting the intensity of the GCM/RCM data
#'
#'Two parameter (shape and scale factors) gamma distribution function was used to model the frequency distributions of the rainfall data. The GCM/RCM rainfall above the threshold were corrected using unique correction factors for different cumulative frequencies.
#'
#'corrected_mod_fut = mod_fut * F-1(F.mod_fut, sh_obs_c,,sc_obs_c)/ F-1 (F.mod_fut,sh_mod_c,,sc_mod_c)
#'
#'where obs - observational data;  mod - GCM/RCM data; n - constant; c - calibration; v - validation; fut - future data; sh - shape factor; sc- scale factor; F. - cumulative density function and F-1 - inverse of cumulative density function
#' @export
#' @examples
#'
#' #subsetting dat_model
#'    mod_calibration=subset(data_model,(year==2003|year==2005|year==2007|year==2009|year==2011))
#'    mod_validation= subset(data_model,(year==2004|year==2006|year==2008|year==2010|year==2012))
#' #subsetting data_observation
#'    obs_calibration=subset(data_observation,(year==2003|year==2005|year==2007|year==2009|year==2011))
#'    obs_validation=subset(data_observation,(year==2004|year==2006|year==2008|year==2010|year==2012))
#' #creating the input vectors
#'    obs_c=obs_calibration$pr
#'    mod_c=mod_calibration$pr
#'    obs_v=obs_validation$pr
#'    mod_v=mod_validation$pr
#'    mod_fut= data_model_future$pr
#'
#'    ResVal(obs_c,mod_c,obs_v,mod_v,mod_fut)
#' @return NULL
ResVal= function(obs_c,mod_c,obs_v,mod_v,mod_fut) {

  m_obs_c=mean(obs_c)
  m_mod_c=mean(mod_c)
  m_obs_v=mean(obs_v)
  m_mod_v=mean(mod_v)


  s_obs_c=stats::sd(obs_c)
  s_mod_c=stats::sd(mod_c)
  s_obs_v=stats::sd(obs_v)
  s_mod_v=stats::sd(mod_v)


  sh_obs_c=(m_obs_c/s_obs_c)^2
  sh_mod_c=(m_mod_c/s_mod_c)^2
  sh_obs_v=(m_obs_v/s_obs_v)^2
  sh_mod_v=(m_mod_v/s_mod_v)^2


  sc_obs_c=s_obs_c^2/m_obs_c
  sc_mod_c=s_mod_c^2/m_mod_c
  sc_obs_v=s_obs_v^2/m_obs_v
  sc_mod_v=s_mod_v^2/m_mod_v



  ### finding thresholds for calibration data  (for both model and obsevartion)

  p=stats::ecdf(obs_c)
  thr_obs_c= p(0.0) # return the probabily at zero}
  threshold_obs_c= (round(p(0.0),3))*1000# return the probabily at zero in 2 digit                   # for calibration_obs

  p=stats::ecdf(mod_c)
  thr_mod_c= p(0.0)
  threshold_mod_c= (round(p(0.0),3))*1000                                                            # for calibration_mod


  p=stats::ecdf(obs_v)
  thr_obs_v= p(0.0)
  threshold_obs_v= (round(p(0.0),3))*1000                                                            # for validation_obs

  p=stats::ecdf(mod_v)
  thr_mod_v= p(0.0)
  threshold_mod_v= (round(p(0.0),3))*1000                                                            # for validation_mod


  n=log(thr_mod_c)/log(thr_obs_c)

  #---------------------------------------------------------------------------------------------------------------------------------------------------

  #----------------------------------------------------------------------------------------------------------------------------------------------

  ### Bias correction_validation

  xx_mod_v=mod_v
  ff_mod_v=c(1)
  crtd_mod_v=c(1)


  p=stats::ecdf(mod_v)
  thr= p(0.0)
  threshold=thr^(1/n)                   # this is new from the convensional bias correction

  #threshold                                           #  turn on/off to view it

  for (i in 1:length(xx_mod_v))
    ff_mod_v[i]= stats::pgamma(xx_mod_v[i],sh_mod_v,,sc_mod_v)

  #ff_mod_v                                            #  turn on/off to view it

  for (i in 1:length(ff_mod_v))
    if (ff_mod_v[i] >= threshold) {
      crtd_mod_v[i]= xx_mod_v[i] * stats::qgamma(ff_mod_v [i],sh_obs_c,,sc_obs_c) /stats::qgamma(ff_mod_v [i],sh_mod_c,,sc_mod_c)
    } else{
      crtd_mod_v[i]=0.0}


  #crtd_mod_v                                           #  turn on/off to view it

  #--------------------------------------------------------------------------------------------------------------------------------------------

  #-----------------------------------------------------------------------------------------------------------------------------------------------


  ### display results of the validation/summary


  # plotting /generating gamma curves for validation

  # modelling data using gamma distribution (above thresold)_ Observation
  # selecting threshold to plot thresholds (here, I have used the threshold of same data)


  p=stats::ecdf(obs_v)
  threshold_obs_v= (round(p(0.0),3))*1000
  #threshold_obs_v                               #  turn on/off to view it


  # modeling obs_val

  f_obs_v=c(1)
  x_obs_v=c(1)                                   # defining variables

  for (i in threshold_obs_v:999)
    f_obs_v[i-threshold_obs_v+1]=i/1000

  for (i in 1:length(f_obs_v))
    x_obs_v[i]=round((stats::qgamma(f_obs_v[i],sh_obs_v,,sc_obs_v)),1)

  #f_obs_v                                      #  turn on/off to view it
  #x_obs_v                                      #  turn on/off to view it

  p.plot= graphics::plot(x_obs_v,f_obs_v,xlim=c(0,50),type="l",col= "green")


  #x=x_obs_v
  x= obs_v
  p=stats::ecdf(x)
  graphics::lines(p,col="green")                          #  turn on/off to display actual data


  #-------------------------------------------------------------------------------------------


  #modelling data using gamma distribution (above thresold)_model
  # selecting threshold to plot thresholds (here, I have used the threshold of same data)

  p=stats::ecdf(mod_v)
  threshold_mod_v= round(p(0.0),3)*1000 # return the probabily at zero in 2 digit

  # threshold_mod_v                     #  turn on/off to view it

  f_mod_v=c(1)
  x_mod_v=c(1) # defining variables

  for (i in threshold_mod_v:999)
    f_mod_v[i-threshold_obs_v+1]=i/1000

  for (i in 1:length(f_mod_v))
    x_mod_v[i]=round((stats::qgamma(f_mod_v[i],sh_mod_v,,sc_mod_v)),1)

  #f_mod_v                               #  turn on/off to view it
  #x_mod_v                               #  turn on/off to view it

  graphics::lines(x_mod_v,f_mod_v,type="l",col= "red")

  #x=x_mod_v
  x= mod_v
  p=stats::ecdf(x)
  graphics::lines(p,col= "red")                      # fitting true values on to modeled curve, turn on/off to display actual data



  #-----------------------------------------------------------------------------------------



  # Modeling corrected data using gamma distribution
  # selecting threshold to plot thresholds (here, I have used the threshold of same data)

  p=stats::ecdf(crtd_mod_v)
  threshold_crtd_mod_v= round(p(0.0),3)*1000 # return the probabily at zero in 2 digit

  # threshold_crtd_mod_v                     #  turn on/off to view it


  #parameters
  m_crtd_mod_v= mean(crtd_mod_v)
  s_crtd_mod_v= stats::sd(crtd_mod_v)
  sh_crtd_mod_v= (m_crtd_mod_v/s_crtd_mod_v)^2
  sc_crtd_mod_v= s_crtd_mod_v^2/m_crtd_mod_v



  f_crtd_mod_v=c(1)
  x_crtd_mod_v=c(1) # defining variables

  for (i in threshold_crtd_mod_v:999)
    f_crtd_mod_v[i-threshold_obs_v+1]=i/1000

  for (i in 1:length(f_mod_v))
    x_crtd_mod_v[i]=round((stats::qgamma(f_crtd_mod_v[i],sh_crtd_mod_v,,sc_crtd_mod_v)),1)

  #f_crtd_mod_v                                            #  turn on/off to view it
  #x_crtd_mod_v                                            #  turn on/off to view it


  graphics::lines(x_crtd_mod_v,f_crtd_mod_v,type="l",col= "blue")


  #x=x_crtd_mod_v,f_crtd
  x= crtd_mod_v
  p=stats::ecdf(x)
  graphics::lines(p,col="blue")                                       # fitting true values on to modeled curve, turn on/off to display actual data




  #---------------------------------------------------------------------------------------------------------------------------------------------------


  #VALIDATION SUMMARY/RESULT



  ### creating th A matrix for validation index


  f=c(1)



  for (i in 1:1000)
    f[i]=(i-1)/1000

  #f
  #crtd_mod_v
  #obs_v$pr
  #mod_v$pr

  xx_obs_v=c(1)
  xx_mod_v=c(1)
  xx_crtd_mod_v=c(1)

  for(i in 1:length(f))   {
    xx_obs_v[i]=round((stats::qgamma(f[i], sh_obs_v,,sc_obs_v)),1)
    xx_mod_v[i]=round((stats::qgamma(f[i], sh_mod_v,,sc_mod_v)),1)
    xx_crtd_mod_v[i] =round((stats::qgamma(f[i], sh_crtd_mod_v,,sc_crtd_mod_v)),1)}

  f= matrix(f)
  xx_crtd_mod_v= matrix(xx_crtd_mod_v)
  xx_obs_v=matrix(xx_obs_v)
  xx_mod_v= matrix(xx_mod_v)

  A=cbind(f,xx_obs_v,xx_mod_v,xx_crtd_mod_v)

  #A
  #---------------------------------------------------------------------------------------------------------

  ### RMSV

  rmsv_obs_vs_mod=c(1)
  rmsv_obs_vs_crmod=c(1)


  for (i in 1:length(xx_obs_v))
    rmsv_obs_vs_mod[i]=((xx_mod_v[i]-xx_obs_v[i])^2)^.5

  for (i in 1:length(xx_obs_v))
    rmsv_obs_vs_crmod[i]=((xx_crtd_mod_v[i]-xx_obs_v[i])^2)^.5

  rmsv1=mean(rmsv_obs_vs_mod)
  rmsv2=mean(rmsv_obs_vs_crmod)

  #rmsv1
  #rmsv2      # has to be transformed into summmary

  #----------------------------------------------------------------------------------------------------


  ### gradient

  #plot (xx_obs_v,xx_crtd_mod_v,type = "l",col="blue")
  #lines(xx_obs_v,xx_obs_v,type = "l",col="green")
  #lines(xx_obs_v,xx_mod_v,type = "l",col="red")            #  turn on/off to view it


  graphics::plot (xx_obs_v,xx_obs_v,col="green")
  graphics::lines(xx_obs_v,xx_obs_v,type = "l",col="green")


  #lines (xx_obs_v,xx_crtd_mod_v,col="blue")
  fit1=stats::lm(xx_crtd_mod_v~xx_obs_v)
  graphics::abline(fit1, col="blue")
  #fit1                                                      #  turn on/off to view it
  # summary(fit1)                                           #  turn on/off to view it


  #lines (xx_obs_v,xx_mod_v)                                #  turn on/off to view it
  fit2=stats::lm(xx_mod_v~xx_obs_v)
  graphics::abline(fit2,col="red")
  #fit2
  #summary(fit2)                                            #  turn on/off to view it


  #coef(fit1)
  #coef(fit2)





  err=round((((threshold_obs_v- threshold_crtd_mod_v)^2)^.5)/threshold_obs_v*100,2)




  validation_output_list=list(percentage_error_in_the_thresholds =err,root_mean_squre_error_in_mod=rmsv1 ,root_mean_squre_error_in_crtdmod=rmsv2,gradient_of_obs_vs_mod_line=fit2$coefficients[2],gradient_of_obs_vs_crtdmod_line=fit1$coefficients[2])


  print(validation_output_list)
  return(NULL)
}


