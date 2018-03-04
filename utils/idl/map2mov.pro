pro map2mov,imin,imax,colt=colt,type=type,min=min,max=max

  if not keyword_set(colt) then colt=33
  loadct,colt
  tvlct,r,g,b,/get

  if not keyword_set(type) then type='dens'

  car=getcarnum(5000)

  if not keyword_set(min)then begin
     case type of
        'dens': min=0.01
        'temp': min=100.
        else: print, 'unkown type' 
     endcase     
  endif

  if not keyword_set(max)then begin
     case type of
        'dens': max=1e5
        'temp': max=1e6
        else: print, 'unkown type' 
     endcase
  endif

  for i=imin,imax do begin
     print,'====',i,', '+car(i)
     rd_im,d,file=type+'_'+car(i)+'.map'
     print,'Minmax',min(d),max(d)
     dd=d
     dd(0,0)=min
     dd(0,1)=max
     image=bytscl(alog10(dd),min=alog10(min),max=alog10(max))
     write_png,type+'_'+car(i)+'.png',image,r,g,b
  endfor

end
