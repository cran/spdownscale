#' @title Calibration Parameters
#' @description  Displays the shape factors, scale factors and the threshold values of the observation and GCM/RCM data set which ultimately define the model
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
#'    ParaCal(obs_c,mod_c,obs_v,mod_v,mod_fut)
#' @return NULL
ParaCal = function(obs_c,mod_c,obs_v,mod_v,mod_fut) {
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
  #returns
  shape_factors=matrix(c(sh_obs_c,sh_obs_v,sh_mod_c,sh_mod_v),nrow = 2)
  colnames(shape_factors)=c("obs","mod")
  rownames(shape_factors)=c("calibration","validation")
  shape_factors=round(shape_factors,4)
  scale_factors=matrix(c(sc_obs_c,sc_obs_v,sc_mod_c,sc_mod_v),nrow = 2)
  colnames(scale_factors)=c("obs","mod")
  rownames(scale_factors)=c("calibration","validation")
  scale_factors=round(scale_factors,4)
  thresholds=matrix(c(threshold_obs_c,threshold_obs_v,threshold_mod_c,threshold_mod_v),nrow = 2)
  colnames(thresholds)=c("obs","mod")
  rownames(thresholds)=c("calibration","validation")
  thresholds=round(thresholds/1000,4)
  n=n
  parameter_list = list(shape_factors=shape_factors, scale_factors=scale_factors, thresholds=thresholds,n=n)
  parameter_list
  print(parameter_list)
  return (NULL)
}

