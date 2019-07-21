;+
; NAME: 
;				CPS
;
; PURPOSE: 
;				to save in format .ps
;
; CALLING:
;				cps,PRINTIT=print,PRINTER=printer,COLOR=color
;
; INPUTS: 
;				
;	
; OUTPUT: 
;	
;
; HISTORY:
;	Written: Sandrine Pires 2005.
;-
;-------------------------------------------------------------------------------
pro cps,PRINTIT=print,PRINTER=printer,COLOR=color
device,/close
set_plot,'x'
if not keyword_set(printer) then printer='hp4'
if keyword_set(color) then printer='color'
if keyword_set(print) then begin
  spawn,'lpr -P'+printer+' idl.ps'
  spawn,'lpq -P'+printer
endif
cpub

end
