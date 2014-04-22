;+
; NAME:  
;       AM_GETMETAL_LATE
; PURPOSE:
;        Compute metallicity of M4.5-M9.5 using an IRTF spectrum. Uses
;        the empirical formula described in Mann et al. (ApJ, in
;        review). This assumes that you have already put the spectrum in
;        vacuum and the rest frame of the star (correcting for RV
;        variations). 
;        Preferably you should use the same stars for the RV
;        correction as I did for consistency: For stars of M4.5 to
;        M6.5 we used the M5V template Gl 51, and for later-type stars
;        we used the M9V template LHS 2065. These templates can be
;        found at:
;        http://irtfweb.ifa.hawaii.edu/~spex/IRTF_Spectral_Library/
;        This software can be found at https://github.com/awmann/metal
;
;        Contact me at amann@astro.as.utexas.edu if you have problems
;        or suggested improvments. 
;
; Giving Credit:
;        If you make use of this code please cite Mann et al. (2014):
;        http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1403.5560
;        If you use the IRTF templates for RV correction make sure to
;        cite Cushing et al. (2005) and Rayner et al. (2009).
;
;
; CALLING SEQUENCE: 
;       AM_GETMETAL_LATE,lambda,spec,error,feh,feh_err,doerr=doerr
;
; INPUT PARAMETERS: 
;       lambda   Wavelength in microns
;       spec     Spectrum in erg/s/cm^2/A (the absolute units
;                are not important, just make sure the units are the same as
;                the error spectrum)
;       error    Error spectrum in erg/s/cm^2/A 
;
; OUTPUT PARAMETERS: 
;       feh      [Fe/H]
;       feh_err  error in [Fe/H], -99 if error spetrum not provided
;                (only based on SNR, does not include errors in calibration)
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       /doerr   set if you want to calculate monte carlo errors,
;       default = do errors (reccomended) 
;   
; EXAMPLES:
;       Download the example sample_late.fits
;       temp = mrdfits('sample_late.fits',0)
;       lambda = temp[*,0] 
;       spec = temp[*,1]
;       error = temp[*,2]
;       AM_GETMETAL_LATE,lambda,spec,error,feh,feh_Err 
;       print,feh,feh_Err
;
;       should print something like:
;       0.0036626064     0.061698792
;
;       The metallaciity of this star's primary is −0.01 ±
;       0.10, so everything is consistent.
;       
; RESTRICTIONS:
;       Users are cautioned that:
;       -This code is currently ONLY designed to work with IRTF
;       spectra. Use of this program on spectra of other instruments
;       will require correction factors due to differences in
;       resolution (among other things). I suggest collecting spectra
;       of >=10 stars from the binary sample to correct for systematic
;       offsets between the instruments.
;       -This program will only function if the wavelengths are given
;       in microns.
;       -This is only designed for M4.5-M9.5 (see am_getmetal.pro for
;       K5-M5).
;       -Reported errors ignore errors in the calibration itself. You
;       should add in a ~0.07 dex systematic calibration error
;       -This technique is untested for stars more metal-poor than
;       [Fe/H]=-0.60 or more rich than [Fe/H]=+0.5. Output outside
;       this range results should be treated with skepticism.
;
;
; DEPENDENCIES:
;
;
;       
; HISTORY: 
;       Routine written. A. Mann 11/18/2013
;       Squished bugs and added comments. A. Mann 01/20/2014
;       Formatted for GitHub 02/05/2014
;       Updated for publication 03/16/2014
;       
;-

;; just kind of a wrapper for pseudo_eqw
Pro metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,metal
  featurearr = dblarr(featurenum+1)
  for j = 0,featurenum-1 do begin
     getclosecontregion,lambda,spec,error,centers[j],continuum
     pseudo_eqw,lambda,spec,continuum,centers[j],widths[j],eqw
     featurearr[j] = eqw*10000.0
  endfor
  featurearr[j] = Te 
  metal = 0.0
  for k = 0,featurenum do metal+= featurearr[k]*coeffs[k] ;; herein lies the metallicity calibration!
  metal+=coeffs[k]
END


;; Meausre the equivalent width of a given feature.
;; Formula = sum(1-(Fl/Fc)*Dl)
PRO pseudo_eqw,lambda,flux,continuum,center,width,eqw
  eqrange = [center-width/2.0,center+width/2.0]
  featureloc = where(lambda ge eqrange[0] and lambda le eqrange[1])
  featurecont = continuum[featureloc]
  featureflux = flux[featureloc]
  dellambda = abs(lambda - shift(lambda,1))
  dellambda = dellambda[featureloc]
  eqw = total((1.0 - (featureflux/featurecont))*dellambda)
END


;; used for getting continuum in the optical
;; It grabs the continuum region just blueward and just redward of
;; center (the feature center).
;; WARNING: has some issues if continuum and feature overlap
;; Returns the continuum applicable to the given region
PRO getclosecontregion,blambda,bspec,berror,center,newcontinuum

  lambda = blambda
  spec = bspec
  error = berror

  continuum = getcontregions()
  leftcont = continuum[0,*]
  rightcont = continuum[1,*]
  
  lowloc = (where(leftcont - center lt 0))
  lowloc = lowloc[n_elements(lowloc)-1]
  highloc = (where(rightcont - center gt 0))[0]
  continuum = [[continuum[*,lowloc]],[continuum[*,highloc]]]

  contloc = dblarr(1)
  for j = 0,n_elements(continuum[0,*])-1 do contloc = [contloc,where(lambda ge continuum[0,j] and lambda le continuum[1,j])]
  shrink,contloc
  contlambda = lambda[contloc]
  contflux = spec[contloc]
  conterr = error[contloc]

  result = linfit(contlambda,contflux,measure_errors=conterr)
  newcontinuum = lambda*result[1] + result[0]

END

;; H2O-K2 (Rojas-Ayala et al.)
PRO geth20k2,newspec,lambda,h20
  loc1 = where(lambda gt 2.070 and lambda lt 2.090)
  loc2 = where(lambda gt 2.235 and lambda lt 2.255)
  loc3 = where(lambda gt 2.360 and lambda lt 2.380)
  h20 = (median(newspec[loc1]) / median(newspec[loc2])) / (median(newspec[loc2]) / median(newspec[loc3]))
END

PRO shrink, array
  if n_elements(array) gt 1 then array=array[1:n_elements(array)-1]
END

;; This program just holds all my favorite continuum regions.
FUNCTION getcontregions

  continuum =  1d0*[[0.32, 0.335], $
                    ;; B
                    [0.3520, 0.3565], [0.3660,0.3690], [0.3755,0.3785], [0.3865,0.3900], [0.4035,0.4080], [0.4135,.4180], [0.4425,0.4450], [0.4480, 0.4510], $
                    [0.4562,0.4575], [0.461,0.4625],  [0.468,0.4700], [0.4856,0.4865], [0.4895,0.491], [0.5050,0.5080], $
                    ;; R
                    [0.5269,0.5299], [.5660,.5675], [0.6586,0.6607], [0.7041,0.7049], [0.7390,0.7500], [0.810,0.816], [0.823,0.830], $
                    [0.8590,0.8620], [0.8890,0.8920], [0.91,0.912],[.9220,0.9255],[0.979,0.989], $
                    ;; J
                    [1.061,1.065], $
                    [1.126,1.130], [1.153,1.158], [1.189,1.193],[1.214,1.218],[1.225,1.23],[1.255,1.2634],[1.270,1.273],[1.295,1.297], $
                    [1.304,1.307],[1.3214,1.327], [1.409,1.415], [1.444,1.448], $
                    ;; H
                    [1.445,1.455], [1.4644,1.4710],[1.4921,1.4965],[1.5060,1.5090], [1.5190,1.522], $ ; [1.5396,1.5430],
                    [1.5920,1.5960],[1.623,1.631], [1.6935,1.6980],[1.7530,1.7570], $
                    ;; K
                    [1.8860,1.8900], [1.935,1.940], [1.961,1.97], [2.05,2.054], [2.08,2.0870], [2.1330,2.1351], [2.153,2.159], $
                    [2.167,2.172], [2.1940,2.1985], [2.213,2.219], [2.245,2.252], [2.2717,2.2781], $
                    [2.285,2.29], [2.305,2.3105], [2.360,2.364], [2.371,2.376], [2.395,2.405], $
                    ;; extra
                    [2.42,2.43] ]
  return,continuum

END


;; primary program
PRO am_getmetal_late,lambda,spec,error,feh,feh_err,doerr=doerr

  on_error,2

  basespec = spec ;; hold the base spectrum so we don't accidentally edit it

  ;; just some error handling and variable initialization
  if n_params() lt 4 then begin
     print,'Syntax -  GETMETAL,lambda,spec,error,feh,feh_err,doerr=doerr'
     return
  endif
  if (n_elements(lambda) ne n_elements(spec)) or (n_elements(spec) ne n_elements(error)) then begin
     print,'wavelength, spectrum, and error (if provided) arrays must have the same number of elements'
     feh = -99
     feh_err = -99
     stop
     return
  endif
  if n_elements(doerr) gt 1 or n_elements(doerr) eq 0 then doerr = 1 ;; default, do error analysis
  Centers = [2.2079, 2.2639999]                                      ;; Na, Ca
  Widths = [0.0068, 0.0059]                                          ;; in microns
  coeffs = [0.13116341, 0.21012097, -3.0679554, 1.3409050]           ;; Coefficients with full digits
  ;;coeffs = [0.131, 0.210, -3.07, 1.341]           ;; Coefficients from the paper (less sig figs but gives same standard deviation and r^2_ap)
  if n_elements(doerr) gt 1 or n_elements(doerr) eq 0 then doerr = 1 ;; default, do error analysis
  Centers = [2.2079, 2.264]                                          ;; Na, Ca
  Widths = [0.0068, 0.0059]                                          ;; in microns
  coeffs = [0.13116341, 0.21012097, -3.0679554, 1.3409050]           ;; Coefficients with all digits (not actually required)
  ;;coeffs = [0.131, 0.210, -3.07, 1.341] ;; Coefficients from the
  ;;paper (less sig figs but gives same standard deviation and
  ;;r^2_ap). Leave for testing!
  featurenum = 2                                   ;; 2 Features (Na and Ca), this is included to keep the program generic
  
  geth20k2,spec,lambda,Te                                             ;; H2O-K2 from Rojas-Ayala et al. (2012)
  metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh ;; actual metallicity calculation 
  feh_final = feh
  
  ;; error analysis (Measurement error only)
  if doerr eq 1 then begin
     nmonte = 500          ;; 100 is usually enough, decrease for better speed
     fehs = dblarr(nmonte) ;; variable initialization
     for i = 0,nmonte-1 do begin
        spec = basespec + error*randomn(seed,n_elements(basespec))          ;; add some noise to the spectrum
        geth20k2,spec,lambda,Te                                             ;; H2O-K2 from Rojas-Ayala et al. (2012)
        metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh ;; actual metallicity calculation 
        fehs[i] = feh
     endfor
     feh_err = robust_sigma(fehs) ;;
  endif else feh_Err = -99
  feh = feh_final              ;;

  spec = basespec ;; restore base spectrum

END
