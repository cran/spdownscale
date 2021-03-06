% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ResVal.R
\name{ResVal}
\alias{ResVal}
\title{Validation Summary}
\usage{
ResVal(obs_c, mod_c, obs_v, mod_v, mod_fut)
}
\arguments{
\item{obs_c}{vector of observational climate data (rainfall) used for calibrating the model}

\item{mod_c}{vector of GCM/RCM climate data (rainfall) used for calibrating the model}

\item{obs_v}{vector of observational climate data (rainfall) used for validating the model}

\item{mod_v}{vector of GCM/RCM climate data (rainfall) used for validating the model}

\item{mod_fut}{vector of GCM/RCM future climate data (rainfall) need to be downscaled}
}
\description{
Displays the summary of the validation.
}
\details{
1)	Dry-days correction / Defining threshold values

   The relationship between the cumulative frequencies (thresholds) corresponding to the dry days of GCM/RCM data and that of the observational data is defined by a polynomial function given by;

threshold_obs = (threshold_mod)^n

n = ln(threshold_obs_c) / ln(threshold_mod_c)


2)	wet-days correction / Correcting the intensity of the GCM/RCM data

Two parameter (shape and scale factors) gamma distribution function was used to model the frequency distributions of the rainfall data. The GCM/RCM rainfall above the threshold were corrected using unique correction factors for different cumulative frequencies.

corrected_mod_fut = mod_fut * F-1(F.mod_fut, sh_obs_c,,sc_obs_c)/ F-1 (F.mod_fut,sh_mod_c,,sc_mod_c)

where obs - observational data;  mod - GCM/RCM data; n - constant; c - calibration; v - validation; fut - future data; sh - shape factor; sc- scale factor; F. - cumulative density function and F-1 - inverse of cumulative density function
}
\examples{

#subsetting dat_model
   mod_calibration=subset(data_model,(year==2003|year==2005|year==2007|year==2009|year==2011))
   mod_validation= subset(data_model,(year==2004|year==2006|year==2008|year==2010|year==2012))
#subsetting data_observation
   obs_calibration=subset(data_observation,(year==2003|year==2005|year==2007|year==2009|year==2011))
   obs_validation=subset(data_observation,(year==2004|year==2006|year==2008|year==2010|year==2012))
#creating the input vectors
   obs_c=obs_calibration$pr
   mod_c=mod_calibration$pr
   obs_v=obs_validation$pr
   mod_v=mod_validation$pr
   mod_fut= data_model_future$pr

   ResVal(obs_c,mod_c,obs_v,mod_v,mod_fut)
}

