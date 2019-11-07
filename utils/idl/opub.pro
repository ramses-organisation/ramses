;+
; NAME: 
;				OPUB
;
; PURPOSE: 
;				set (open) system variables to produce a plot with
; 			publication standards. Characters are magnified and lines 
; 			are thickened so that a plot (with form=2 with ops.pro) can
; 			be reduced by about 50% and still be legible (as for eg. a
; 			one colum pplot in an APJ article)
; 			note: see cpub.pro to return to default values
;
; CALLING:
;				opub,CSIZE=csize,THICK=thick
;
; OPTIONAL INPUTS: 
;				csize: character size (default=1.4)
;       thick: line thickness (default=2.8)none
;	
; OUTPUT: 
;	
;
; HISTORY:
;	Written: A. Refregier 1996.
; re-Written: S. Pires 2003.
;-
;-------------------------------------------------------------------------------
pro opub,CSIZE=csize,THICK=thick

if not keyword_set(thick) then thick=3.2 ; thickness of lines and characters
if not keyword_set(csize) then csize=1.7
!p.charsize=csize
;!p.psym=csize
!p.charthick=thick
!p.thick=thick
!x.thick=thick
!y.thick=thick

return
end
