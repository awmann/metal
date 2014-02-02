;+
; NAME:  
;       AM_GETMETAL
; PURPOSE:
;        Compute metallicity of K5-M4 using SNIFS or IRTF
;        spectrum. Uses the empirical formula described in Mann et
;        al. (2013 AJ 145 52). This assumes that you have
;        already put the spectrum in vacuum and the rest frame of the
;        star (correcting for RV variations). 
;        Assumes you have measured the Spectral Type via some other
;        method (I'll post something to do this eventually).
;
;        Contact me at amann@astro.as.utexas.edu if you have problems.
;
; Giving Credit:
;        If you make use of this code please cite the Mann et
;        al. (2013) paper. If you use the K-band relation, cite
;        Rojas-Ayala et al. (2012 ApJ 748 93) for the H2O-K2 index.
;        If you use the H-band relation, cite Terrien et al. (2012 ApJ
;        747 38) for the H2O-H index. If you use the
;        visible-wavelength relation, cite Hawley (2002 AJ 123 3409)
;        for use of the Color-1 relation. 
;
; CALLING SEQUENCE: 
;       AM_GETMETAL,lambda,spec,error,sptnum,band,feh,feh_err,mh,mh_err
;
; INPUT PARAMETERS: 
;       lambda   Wavelength in microns
;       spec     Spectrum in erg/s/cm^2/A (the absolute units
;                are not important, just make sure the units are the same as
;                error)
;       error    Error spectrum in erg/s/cm^2/A 
;       sptnum   Spectral subtypes past M0 (K5 = -2, M0 = 0, M2 = 2)
;       band     'Vis', 'J', 'H', 'K' specifying the wavelength regime
;                of interest
;
; OUTPUT PARAMETERS: 
;       feh      [Fe/H]
;       feh_err  error in [Fe/H], -1 if error spetrum not provided
;                (only based on SNR, does not include errors in calibration)
;       mh       [M/H]
;       mh_err   error in [M/H], -1 if error spetrum not provided
;                (only based on SNR, does not include errors in calibration)
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;       /doerr   set if you want to calculate monte carlo errors,
;       default = do errors (reccomended) 
;   
; EXAMPLES:
;
;       Download the example from
;       ifa.hawaii.edu~/amann/programs/sample.fits
;       temp = mrdfits('sample.fits',0)
;       lambda = temp[*,0] 
;       spec = temp[*,1]
;       error = temp[*,2]
;       AM_GETMETAL,lambda,spec,error,1.8,'Vis',feh,feh_Err,mh,mh_err 
;       print,feh,feh_Err,mh,mh_err
;       AM_GETMETAL,lambda,spec,error,1.8,'J',feh,feh_Err,mh,mh_err
;       print,feh,feh_Err,mh,mh_err 
;       AM_GETMETAL,lambda,spec,error,1.8,'H',feh,feh_Err,mh,mh_err 
;       print,feh,feh_Err,mh,mh_err 
;       AM_GETMETAL,lambda,spec,error,1.8,'K',feh,feh_Err,mh,mh_err 
;       print,feh,feh_Err,mh,mh_err
;
;       will print:
;     -0.19566243     0.056745980     -0.14632768     0.045831947
;     -0.23293379     0.056442321     -0.18807660     0.060947115
;     -0.18803976     0.057410258     -0.20218577     0.044678279
;     -0.19821758     0.032951810     -0.20726450     0.030075555
;
;     This is a good (somewhat atypical) range of consistency between
;     the 4 calibrations. J band in particular is sometimes an
;     outlier. I reccomend using H-band and K-band relations, which
;     are the most consistent/relaible/accurate. 
;       
;
; RESTRICTIONS:
;       -Users are cautioned that this code is currently ONLY designed
;       to work with IRTF and/or SNIFS spectra. Use of this program on
;       spectra of other instruments will require correction factors
;       due to differences in resolution (among other things). I
;       suggest collecting spectra of >=10 stars from the Mann et
;       al. (2013) binary sample to correct for systematic offsets
;       between the instruments.
;       -This program will only function if the wavelengths are given
;       in microns.
;       -Reported errors ignore errors in the calibration itself. You
;       should add the following (systematic) errors for a given
;       calibration: 
;       K-band: 0.08 dex
;       H-band: 0.09 dex
;       J-band: 0.12 dex
;       Vis: 0.10 dex
;       Do not use the (1-sigma) errors from the paper, as these
;       overestimate the error in the fit by attempting to account for
;       SNR errors.
;       -Errors are probably a bit lower for [Fe/H]>-0.5 and a bit
;       higher for [Fe/H]<-0.5. 
;       -This technique is untested past M5, earlier than K4, and for
;       stars more metal-poor than [Fe/H]=-1.0 or more rich than
;       [Fe/H]=+0.5. Output outside this range should not be trusted.
;
;
; DEPENDENCIES:
;
;       getcontregions.pro
;       (http://ifa.hawaii.edu/~amann/programs/getcontregions.pro)
;
;       
; HISTORY: 
;       Routine written. A. Mann 06/10/2013
;       Added more comments. A. Mann 07/01/2013
;       Added calibration errors to comments. A. Mann 10/01/2013
;-


Pro metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,metal

  featurearr = dblarr(featurenum+1)
  for j = 0,featurenum-1 do begin
     getclosecontregion,lambda,spec,error,centers[j],continuum
     pseudo_eqw,lambda,spec,continuum,centers[j],widths[j],eqw
     featurearr[j] = eqw*10000.0
  endfor
  featurearr[j] = Te 
  metal = 0.0
  for k = 0,featurenum do metal+= featurearr[k]*coeffs[k]
  metal+=coeffs[k]

END

;; This program is use to calculate band indices (as opposed to
;; equivlent widths).
PRO band_index,oldlambda,oldspec,center1,width1,center2,width2,index

  newlambda = oldlambda         ;generatearray(min(oldlambda),max(oldlambda),5d3);
  newspec = oldspec             ;interpol(oldspec,oldlambda,newlambda);

  lambda = newlambda
  spectrum = newspec
  
  flux1 = total(spectrum[where(lambda gt center1-width1/2.0 and lambda lt center1+width1/2.0)])
  flux2 = total(spectrum[where(lambda gt center2-width2/2.0 and lambda lt center2+width2/2.0)])

  index = (1. - (FLUX1/WIDTH1)/(FLUX2/WIDTH2))*WIDTH1
  index = (FLUX1/WIDTH1)/(FLUX2/WIDTH2)

END

PRO shrink, array
  if n_elements(array) gt 1 then array=array[1:n_elements(array)-1]
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


;; Meausre the equivalent width of a given feature.
PRO pseudo_eqw,lambda,flux,continuum,center,width,eqw
  eqrange = [center-width/2.0,center+width/2.0]
  featureloc = where(lambda ge eqrange[0] and lambda le eqrange[1])
  featurecont = continuum[featureloc]
  featureflux = flux[featureloc]
  dellambda = abs(lambda - shift(lambda,1))
  dellambda = dellambda[featureloc]
  eqw = total((1.0 - (featureflux/featurecont))*dellambda)
END

;; H2O-K1 (Covey et al.)
PRO geth20k1,newspec,lambda,h20
  loc1 = where(lambda gt 2.18 and lambda lt 2.20)
  loc2 = where(lambda gt 2.27 and lambda lt 2.29)
  loc3 = where(lambda gt 2.36 and lambda lt 2.38)
  h20 = (median(newspec[loc1]) / median(newspec[loc2])) / (median(newspec[loc2]) / median(newspec[loc3]))
END

;; H2O-K2 (Rojas-Ayala et al.)
PRO geth20k2,newspec,lambda,h20
  loc1 = where(lambda gt 2.070 and lambda lt 2.090)
  loc2 = where(lambda gt 2.235 and lambda lt 2.255)
  loc3 = where(lambda gt 2.360 and lambda lt 2.380)
  h20 = (median(newspec[loc1]) / median(newspec[loc2])) / (median(newspec[loc2]) / median(newspec[loc3]))
END

;; H2O-H (Terrien et al.)
PRO geth20h,newspec,lambda,h20h
  loc1 = where(lambda gt 1.595 and lambda lt 1.615)
  loc2 = where(lambda gt 1.680 and lambda lt 1.700)
  loc3 = where(lambda gt 1.760 and lambda lt 1.780)
  h20h = (median(newspec[loc1]) / median(newspec[loc2])) / (median(newspec[loc2]) / median(newspec[loc3]))
END

;; H2O-J (Mann et al.)
PRO geth20j,newspec,lambda,h20
  loc1 = where(lambda gt 1.210 and lambda lt 1.230)
  loc2 = where(lambda gt 1.313 and lambda lt 1.333)
  loc3 = where(lambda gt 1.331 and lambda lt 1.351)
  h20 = (median(newspec[loc1]) / median(newspec[loc2])) / (median(newspec[loc2]) / median(newspec[loc3]))
END


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


PRO am_getmetal,lambda,spec,error,sptnum,band,feh,feh_err,mh,mh_err,doerr=doerr

  on_error,2

  basespec = spec

  if n_params() lt 4 then begin
     print,'Syntax -  GETMETAL,lambda,spectrum,sptnum,band,feh,feh_err,mh,mh_err,error=error'
     return
  endif
  if (n_elements(lambda) ne n_elements(spec)) or (n_elements(spec) ne n_elements(error)) then begin
     print,'wavelength, spectrum, and error (if provided) arrays must have the same number of elements'
     feh = -99
     feh_err = -99
     mh = -99
     mh_err = -99
     return
  endif
  if n_elements(doerr) gt 1 or n_elements(doerr) eq 0 then doerr = 1 ;; default, do error analysis

  case band of
     'Vis': begin
        band_index,lambda,spec,.8950,.01,.74,.01,Te ;; calculate the Color-1 index
        if sptnum le 2.0 then begin
           ;; [Fe/H]
           Centers = [0.8208, 0.4648, 0.5608]    
           Widths = [0.0035, 0.0023, 0.0020]   
           coeffs = [0.53,0.26,-0.16,-0.784,-0.34] ;; Feature 7, 1, 2
           featurenum = 3
           metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

           ;; [M/H]
           Centers = [0.8208, 0.4648, 0.8684]      
           Widths = [0.0035, 0.0023, 0.0026]              
           coeffs = [0.38,0.21,0.29,-0.504,-0.79] ;; Feature 7, 1, 8
           featurenum = 3
           metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
        endif else begin
           ;; [Fe/H]
           Centers = [0.6416, 0.8684, 0.8208, 0.6118]      
           Widths = [0.0041, 0.0026, 0.0035, 0.0020]        
           coeffs = [-0.20,0.48,0.24,0.14,-0.204,-0.32] ;; Feature 5, 8, 7, 3
           featurenum = 4
           metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

           ;; [M/H]
           Centers = [0.6416, 0.6232, 0.7540]      
           Widths = [0.0041, 0.0020, 0.0020]       
           coeffs = [-0.065,-0.071,-0.3,0.719,-0.24] ;; Feature 5, 4, 6
           featurenum = 3
           metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
        endelse
     end
     'J': begin
        geth20j,spec,lambda,Te
        ;; [Fe/H]
        Centers = [1.2698, 1.1396, 1.3148, 1.3344]       
        Widths = [0.0098, 0.0026, 0.0050, 0.0023]    
        coeffs = [0.29, 0.21, 0.26, -0.26, -0.190, -1.03] ;; Feature 10, 9, 12, 13
        featurenum = 4
        metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

        ;; [M/H]
        Centers = [1.2698, 1.2908, 1.1396]       
        Widths = [0.0098, 0.0020, 0.0026]    
        coeffs = [0.32, 0.46, 0.076, 1.213, -1.97] ;; Feature 10, 11, 9
        featurenum = 3
        metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
     end
     'H': begin
        geth20h,spec,lambda,Te
        ;; [Fe/H]
        Centers = [1.6158, 1.4766, 1.7261]       
        Widths = [0.0023, 0.0041, 0.0032]    
        coeffs = [0.40, 0.51, -0.28, -1.460, 0.71] ;; Feature 17, 14, 18
        featurenum = 3
        metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

        ;; [M/H]
        Centers = [1.6158, 1.5172, 1.4836]       
        Widths = [0.0023, 0.0033, 0.0023]    
        coeffs = [0.38, 0.40, 0.41, .194, -0.76] ;; Feature 17, 16, 15
        featurenum = 3
        metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
     end
     'K': begin
        geth20k2,spec,lambda,Te
        ;; [Fe/H]
        Centers = [2.2079, 2.3844, 2.3242]       
        Widths = [0.0068, 0.0035, 0.0038]    
        coeffs = [0.19, 0.069, 0.083, 0.218, -1.55] ;; Feature 19, 22, 20
        featurenum = 3
        metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

        ;; [M/H]
        Centers = [2.2079, 2.3844, 2.3342]       
        Widths = [0.0068, 0.0035, 0.0035]    
        coeffs = [0.12, 0.086, 0.13, .245, -1.18] ;; Feature 19, 22, 21
        featurenum = 3
        metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
     end
     default: begin
        print,'band must be Vis, J, H, or K'
        return
     end
  endcase
  feh_final = feh
  mh_final = mh

  if doerr eq 1 then begin
     nmonte = 500
     
     fehs = dblarr(nmonte)
     mhs = dblarr(nmonte)

     for i = 0,nmonte-1 do begin
        spec = basespec + error*randomn(seed,n_elements(basespec))

        case band of
           'Vis': begin
              band_index,lambda,spec,.8950,.01,.74,.01,Te ;; calculate the Color-1 index
              if sptnum le 2.0 then begin
                 ;; [Fe/H]
                 Centers = [0.8208, 0.4648, 0.5608]    
                 Widths = [0.0035, 0.0023, 0.0020]   
                 coeffs = [0.53,0.26,-0.16,-0.784,-0.34] ;; Feature 7, 1, 2
                 featurenum = 3
                 metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

                 ;; [M/H]
                 Centers = [0.8208, 0.4648, 0.8684]      
                 Widths = [0.0035, 0.0023, 0.0026]              
                 coeffs = [0.38,0.21,0.29,-0.504,-0.79] ;; Feature 7, 1, 8
                 featurenum = 3
                 metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
              endif else begin
                 ;; [Fe/H]
                 Centers = [0.6416, 0.8684, 0.8208, 0.6118]      
                 Widths = [0.0041, 0.0026, 0.0035, 0.0020]        
                 coeffs = [-0.20,0.48,0.24,0.14,-0.204,-0.32] ;; Feature 5, 8, 7, 3
                 featurenum = 4
                 metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

                 ;; [M/H]
                 Centers = [0.6416, 0.6232, 0.7540]      
                 Widths = [0.0041, 0.0020, 0.0020]       
                 coeffs = [-0.065,-0.071,-0.3,0.719,-0.24] ;; Feature 5, 4, 6
                 featurenum = 3
                 metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
              endelse
           end
           'J': begin
              geth20j,spec,lambda,Te
              ;; [Fe/H]
              Centers = [1.2698, 1.1396, 1.3148, 1.3344]       
              Widths = [0.0098, 0.0026, 0.0050, 0.0023]    
              coeffs = [0.29, 0.21, 0.26, -0.26, -0.190, -1.03] ;; Feature 10, 9, 12, 13
              featurenum = 4
              metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

              ;; [M/H]
              Centers = [1.2698, 1.2908, 1.1396]       
              Widths = [0.0098, 0.0020, 0.0026]    
              coeffs = [0.32, 0.46, 0.076, 1.213, -1.97] ;; Feature 10, 11, 9
              featurenum = 3
              metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
           end
           'H': begin
              geth20h,spec,lambda,Te
              ;; [Fe/H]
              Centers = [1.6158, 1.4766, 1.7261]       
              Widths = [0.0023, 0.0041, 0.0032]    
              coeffs = [0.40, 0.51, -0.28, -1.460, 0.71] ;; Feature 17, 14, 18
              featurenum = 3
              metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

              ;; [M/H]
              Centers = [1.6158, 1.5172, 1.4836]       
              Widths = [0.0023, 0.0033, 0.0023]    
              coeffs = [0.38, 0.40, 0.41, .194, -0.76] ;; Feature 17, 16, 15
              featurenum = 3
              metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
           end
           'K': begin
              geth20k2,spec,lambda,Te
              ;; [Fe/H]
              Centers = [2.2079, 2.3844, 2.3242]       
              Widths = [0.0068, 0.0035, 0.0038]    
              coeffs = [0.19, 0.069, 0.083, 0.218, -1.55] ;; Feature 19, 22, 20
              featurenum = 3
              metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,feh

              ;; [M/H]
              Centers = [2.2079, 2.3844, 2.3342]       
              Widths = [0.0068, 0.0035, 0.0035]    
              coeffs = [0.12, 0.086, 0.13, .245, -1.18] ;; Feature 19, 22, 21
              featurenum = 3
              metalcalc,centers,widths,coeffs,featurenum,lambda,spec,error,te,mh
           end
           default: begin
              print,'band must be Vis, J, H, or K'
              return
           end
        endcase
        fehs[i] = feh
        mhs[i] = mh
     endfor
     feh_err = robust_sigma(fehs)
     mh_err = robust_sigma(mhs)
  endif else begin
     feh_err = -99
     mh_err = -99
  endelse

  mh = mh_final
  feh = feh_final

  spec = basespec

END
