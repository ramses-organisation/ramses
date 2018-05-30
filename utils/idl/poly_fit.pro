; $Id: poly_fit.pro,v 1.16 2003/05/16 17:27:36 chris Exp $
;+
; NAME:
;   POLY_FIT
;
; PURPOSE:
;   Perform a least-square polynomial fit with optional error
;estimates.
;
;   This routine uses matrix inversion.  A newer version of this
;routine,
;   SVDFIT, uses Singular Value Decomposition.  The SVD technique is
;more
;   flexible, but slower.
;
; CATEGORY:
;   Curve fitting.
;
; CALLING SEQUENCE:
;   Result = POLY_FIT(X, Y, Degree)
;
; INPUTS:
;   X:  The independent variable vector.
;
;   Y:  The dependent variable vector, should be same length as x.
;
;   Degree: The degree of the polynomial to fit.
;
; OUTPUTS:
;   POLY_FIT returns a vector of coefficients with a length of
;NDegree+1.
;
; KEYWORDS:
;   CHISQ:   Sum of squared errors divided by MEASURE_ERRORS if
;specified.
;
;   COVAR:   Covariance matrix of the coefficients.
;
;   DOUBLE:  if set, force computations to be in double precision.
;
;   MEASURE_ERRORS: Set this keyword to a vector containing standard
;       measurement errors for each point Y[i].  This vector must be
;the same
;       length as X and Y.
;
;     Note - For Gaussian errors (e.g. instrumental uncertainties),
;        MEASURE_ERRORS should be set to the standard
;        deviations of each point in Y. For Poisson or statistical
;weighting
;        MEASURE_ERRORS should be set to sqrt(Y).
;
;   SIGMA:   The 1-sigma error estimates of the returned parameters.
;
;     Note: if MEASURE_ERRORS is omitted, then you are assuming that
;           your model is correct. In this case, SIGMA is multiplied
;           by SQRT(CHISQ/(N-M)), where N is the number of points
;           in X and M is the number of terms in the fitting function.
;           See section 15.2 of Numerical Recipes in C (2nd ed) for
;details.
;
;   STATUS = Set this keyword to a named variable to receive the
;status
;          of the operation. Possible status values are:
;          0 for successful completion, 1 for a singular array (which
;          indicates that the inversion is invalid), and 2 which is a
;          warning that a small pivot element was used and that
;significant
;          accuracy was probably lost.
;
;    Note: if STATUS is not specified then any error messages will be
;output
;          to the screen.
;
;   YBAND:  1 standard deviation error estimate for each point.
;
;   YERROR: The standard error between YFIT and Y.
;
;   YFIT:   Vector of calculated Y's. These values have an error
;           of + or - YBAND.
;
; COMMON BLOCKS:
;   None.
;
; SIDE EFFECTS:
;   None.
;
; MODIFICATION HISTORY:
;   Written by: George Lawrence, LASP, University of Colorado,
;       December, 1981.
;
;   Adapted to VAX IDL by: David Stern, Jan, 1982.
;       Modified:    GGS, RSI, March 1996
;                    Corrected a condition which explicitly converted
;all
;                    internal variables to single-precision float.
;                    Added support for double-precision inputs.
;                    Added a check for singular array inversion.
;            SVP, RSI, June 1996
;                     Changed A to Corrm to match IDL5.0 docs.
;                    S. Lett, RSI, December 1997
;                     Changed inversion status check to check only for
;                     numerically singular matrix.
;                    S. Lett, RSI, March 1998
;                     Initialize local copy of the independent
;variable
;                     to be of type DOUBLE when working in double
;precision.
;       CT, RSI, March 2000: Changed to call POLYFITW.
;       CT, RSI, July-Aug 2000: Removed call to POLYFITW,
;                   added MEASURE_ERRORS keyword,
;                   added all other keywords (except DOUBLE),
;                   made output arguments obsolete.
;       CT, RSI, Jan 2003: Combine some vector math expressions,
;                   general code cleanup. About 50% faster.
;-

FUNCTION POLY_FIT, x, y, ndegree, $
    yfit_old, yband_old, yerror_old, corrm_old, $     ; obsolete arguments
    CHISQ=chisq, $
    COVAR=covar, $
    DOUBLE=double, $
    MEASURE_ERRORS=measure_errors, $
    SIGMA=sigma, $
    STATUS=status, $
    YBAND=yband, $
    YERROR=yerror, $
    YFIT=yfit

    COMPILE_OPT idl2

    ON_ERROR,2      ;RETURN TO CALLER IF ERROR

    n = N_ELEMENTS(x)
    IF (n NE N_ELEMENTS(y)) THEN MESSAGE, $
        'X and Y must have same number of elements.'
    m = ndegree + 1 ; # of elements in coeff vec

    double = (N_ELEMENTS(double) GT 0) ? KEYWORD_SET(double) : $
        (SIZE(x,/TNAME) EQ 'DOUBLE') OR (SIZE(y,/TNAME) EQ 'DOUBLE')

    haveMeasureError = (N_ELEMENTS(measure_errors) gt 0)
    sdev = 1d
    if (haveMeasureError) then $
        sdev *= measure_errors
    sdev2 = sdev^2

    haveYband = ARG_PRESENT(yband) || ARG_PRESENT(yband_old)

    ; construct work arrays
    covar = DBLARR(m,m) ; least square matrix, weighted matrix
    b = DBLARR(m)   ; will contain sum weights*y*x^j
    z = 1d          ; polynomial term (guarantees double precision calc)
    wy = DOUBLE(y)
    if (haveMeasureError) then $
        wy /= sdev2

    ; Just fill in some values in case we fail.
    yfit = !VALUES.D_NAN
    yband = !VALUES.D_NAN
    yerror = !VALUES.D_NAN


    covar[0,0] = haveMeasureError ?  TOTAL(1d/sdev2) : n
    b[0] = TOTAL(wy)


    FOR p = 1L,2*ndegree DO BEGIN   ; power loop
        z *= x  ; z is now x^p
        IF p LT m THEN b[p] = TOTAL(wy*z)   ; b is sum weights*y*x^j
        sum =  haveMeasureError ? TOTAL(z/sdev2) : TOTAL(z)
        for j = 0 > (p-ndegree), ndegree < p do $
            covar[j,p-j] = sum
     ENDFOR ; end of p loop, construction of covar and b


    covar = INVERT(TEMPORARY(covar), status)


    IF NOT ARG_PRESENT(status) THEN BEGIN
        CASE status OF
        1: MESSAGE, "Singular matrix detected."
        2: MESSAGE,/INFO, "Warning: Invert detected a small pivot element."
        ELSE:
     ENDCASE
     ENDIF

    if (status eq 1) then begin
        result = !VALUES.D_NAN
        goto, done
     endif


    result = (TEMPORARY(b) # covar)  ; construct coefficients


    ; compute optional output parameters.

    ; one-standard deviation error estimates, init
    yfit = result[ndegree]
    for k = ndegree-1L, 0, -1 do $
        yfit = result[k] + TEMPORARY(yfit)*x  ; sum basis vectors


    ; Vector of parameter errors.
    sigma = SQRT(ABS(covar[lindgen(M)*(M+1)]))


    if (haveMeasureError) then begin

        ; Only do this computation once.
        diff = (yfit - y)^2
        chisq = TOTAL(diff/sdev2)

        ; Experimental variance estimate, unbiased.
        var = (n gt m) ? TOTAL(diff)/(n-m) : 0d

     endif else begin

        ; No weighting, don't need to divide by sdev2.
        chisq = TOTAL((yfit - y)^2)

        ; Experimental variance estimate, unbiased.
        var = (n gt m) ? chisq/(n-m) : 0d

                                ; If MEASURE_ERRORS is omitted, then
                                ; you are assuming that your model
                                ; is correct. In this case, SIGMA is
                                ; multiplied by SQRT(chisq/(n-m)).
                                ; See section 15.2 of Numerical
                                ; Recipes in C (Second Edition) for
                                ; details.
        sigma *= SQRT(chisq/(n-m))

     endelse

    ; Overall fit error.
    yerror = SQRT(var)


    ; Only do this computation if user wants YBAND.
    if (haveYband) then begin
        z = REPLICATE(1d, n)
        yband = REPLICATE(covar[0,0], n)
        FOR p=1L,2*ndegree DO BEGIN ; compute correlated error estimates on y
            z *= x      ; z is now x^p
            sum = 0
            for j=0 > (p - ndegree), ndegree<p do $
                sum += covar[j,p-j]
            yband += sum * z  ; add in all the error sources
         ENDFOR  ; end of p loop

        yband *= var
        IF (MIN(yband) LT 0) OR (MIN(FINITE(yband)) EQ 0) THEN BEGIN
            status = 3
            IF NOT ARG_PRESENT(status) THEN MESSAGE, $
                'Undefined (NaN) error estimate encountered.'
         ENDIF ELSE yband = SQRT( TEMPORARY(yband) )
     endif

done:

    ; If necessary, convert all results to single precision.
    if (not double) then begin
        chisq = FLOAT(chisq)
        covar = FLOAT(covar)
        result = FLOAT(result)
        sigma = FLOAT(sigma)
        var = FLOAT(var)   ; needed for corrm_old below
        yfit = FLOAT(yfit)
        if (haveYband) then $
            yband = FLOAT(yband)
        yerror = FLOAT(yerror)
     endif

; fill in obsolete arguments, if necessary
    IF (N_PARAMS() GT 3) THEN BEGIN
        corrm_old = covar*var   ; convert to correlation matrix
        yerror_old = yerror
        yfit_old = yfit
        if (haveYband) then $
            yband_old = yband
     ENDIF

    RETURN, result
 END
