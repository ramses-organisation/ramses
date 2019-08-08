;+
; NAME: 
;				OPS
;
; PURPOSE: 
;				open a postscript file for plot output. (see cps.pro
; 			to close to .ps file and to print it)
;
; CALLING:
;				
;
; INPUTS: 
;				file: file name (default=idl.ps)
;     	form: page format:
;				0: top half
;				1: full page
;				2: square in the middle of the page (default)
;     	3: reduced in vertical direction
; 			landscape: use landscape mode encap: encapsulates postscript
;		  	color: produce color plot
;       nopub: do not use publication quality style
;       thick: line thickness (default=2.8)
;       csize: character size (default=1.4)
;       bits_per_pixel: number of bits for color table (default=4)v --- array
;	
; OUTPUT: 
;				any subsequent plot graphic output will be stored
; 			into the postscript file. 
;
; HISTORY:
; Feb 94 - Written by A. Refregier
; Nov 02 - Modified by Richard Massey to enable the bits_per_pixel option
;	re-Written: Sandrine Pires 2005.
;-
;-------------------------------------------------------------------------------
pro ops,FILE=file,FORM=form,LANDSCAPE=landscape,ENCAP=encap,$
  	COLOR=color,nopub=nopub,thick=thick,csize=csize,$
        bits_per_pixel=bits_per_pixel

if not keyword_set(nopub) then opub,thick=thick,csize=csize

set_plot,'ps'
if not keyword_set(file) then file='idl.ps'
if not keyword_set(encap) then encap=0
if not keyword_set(color) then color=0

if not keyword_set(form) then form=2
case form of
  0: device,filename=file,encapsul=encap,color=color,bits_per_pixel=bits_per_pixel
  1: device,filename=file,yoffset=5.,ysize=20.,encapsul=encap,$
	color=color,bits_per_pixel=bits_per_pixel
  2: device,filename=file,yoffset=5.3,ysize=17.,encapsul=encap,$
	color=color,bits_per_pixel=bits_per_pixel
  3: device,filename=file,yoffset=5.3,ysize=14.,encapsul=encap,$
	color=color,bits_per_pixel=bits_per_pixel
  else:	device,filename=file,encapsul=encap,color=color,bits_per_pixel=bits_per_pixel   ; (same as 0)
endcase 

if keyword_set(landscape) then device,/landscape,encapsul=encap,$
	color=color,bits_per_pixel=bits_per_pixel

return
end
