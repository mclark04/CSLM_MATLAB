# CSLM_MATLAB
The Canadian Small Lake Model, implemented in MATLAB for offline use

Written by M.G.Clark, for details contact dr.mg.clark@gmail.com
Current citation is in review, for now, please email me if you publish using this model so I can keep track of the users.

This is an early draft of the project to translate the Canadian Small Lakes Model from FORTRAN to MATLAB.  The original code was provided by Murray MacKay and is detailed in MacKay 2012, MacKay et al. 2016, and MacKay 2019 (see below).  It is currently working for both example files provided (Lake239.ini/Lake239.met and PitLake2019.ini/PitLake2019.met).  There is a reasonable fit between the MATLAB and the FORTRAN models.  More details will be provided in the near future.    Currently the code is broken up into a few sub routes and the main RunLake() function.  The subroutines are described as follows:

# These two are general functions called throughout the codebase:
###   XIT   
Exit function to crash safely, used in FORTRAN to log buggs, here it is just kept for compatabilty
###  EQNST 
Equation of state called by CLASSL and FREECONV, based on Farmer and Carmack, 1981

# These two are called from the main code (RunLake.f) (defined after the below functions)
###   CLASSI 
The machinery of CLASS to deal with atmosphereic inputs to the surface that are still useful for lakes, things like snow density and rain/snow partitioning.
###   CLASSL
This is the main lake model, everything is modelled here or subroutines are called from here.

# The rest of these are core functions called within CLASSL, I have listed them in exicution order
###   SNOWALBA 
Albedo calculations based on surface cover
###   TSLPREP
Preparing for temperature calculations based on class mechanisms defined in: CLASST.f/TPREP.f/TSPREP.f
###   TSOLVE 
Iterative surface temperature calculations for snow/soil from CLASS
###   DRCOEF 
Drag coefficent calculations used for surface stablity corrections from CLASS
###   FLXSURFZ
Not used, would be if the swtich ISLFD is overwritten and then it is used to determine surface stablity corrections.  - Graham has not translated this function
###   TLSPOST
Lake snow termperature and heat flux calculations for flux in/out of the water column.  Based on CLASS subroutines TSPOST.f and WPREP.f
###   SNOVAP 
Calculate sublimation from the snowpack.
###   TMELT 
Calculate the melting of the snowpack.
###   SNINFL 
Calculate rain infiltration into the snowpack
###   SNOWALBW 
Albedo adjustments post snowpack updates 
###   SNOADD 
Accumulation of snow on the ground
###   DRCOEFL 
Turbulant transfer into lake surface.
###   FREECONV 
Convective mixing in the lake's mixing layers
###   LKTRANS 
Light extinction coefficents, constants right now but will be updated in future verions of the CSLM
###   TSOLVL 
Compute the lake surface energy balance in the skin layer
###   MIXLYR 
Compute the mixing layer depth and TKE in the mixing layer from TKE  and Buoyancy
###   SCREENRH 
Caclulate the screen RH.  Discription states it does more, but it doesn't.  Either way MacKay Said it is unneeded
###   SLDIAG 
Not used inless ISLFD is overwritten, alternate surface stablity corrections
###   DIASURFZ 
Not used inless ISLFD is overwritten, alternate surface stablity  corrections
