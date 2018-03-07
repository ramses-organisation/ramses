pro map2mov,imin,imax,colt=colt,type=type,min=min,max=max

  if not keyword_set(colt) then colt=0
  loadct,colt
  tvlct,r,g,b,/get

  if not keyword_set(type) then type='dens'

  car=getcarnum(5000)

  if not keyword_set(min)then begin
     case type of
        'dm': min=1e-12
        'dens': min=0.01
        'temp': min=3.
        'stars': min=1e-13
        else: print, 'unkown type' 
     endcase     
  endif

  if not keyword_set(max)then begin
     case type of
        'dm': max=2e-7
        'dens': max=1e5
        'temp': max=1e7
        'stars': max=1e-7
        else: print, 'unkown type' 
     endcase
  endif

  for i=imin,imax do begin
     print,'====',i,', '+car(i)
     rd_im,d,file=type+'_'+car(i)+'.map',head=h
     print,'Minmax',min(d),max(d)
     dd=d
     dd(0,0)=min
     dd(0,1)=max
     image=bytscl(alog10(dd),min=alog10(min),max=alog10(max))
     time=string(h.t,format='(F5.3)')
     write_png,type+'_'+car(i)+'.png',image,r,g,b
     command='convert '+type+'_'+car(i)+'.png -font courier-bold -pointsize 20 -fill white -annotate +10+20 "time='+time+'" '+type+'_'+car(i)+'_annotated.png'
     spawn,command
;     print,command
  endfor

end
