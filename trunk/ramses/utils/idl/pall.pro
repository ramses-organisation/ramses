
pro pall,nin,typin

car=getcarnum(nin)
charin=car(nin-1)

rd_im,bx,file='blast_'+charin+'_'+typin+'_dirx.map'
rd_im,dx,file='dens_'+charin+'_'+typin+'_dirx.map'
rd_im,tx,file='temp_'+charin+'_'+typin+'_dirx.map'
rd_im,zx,file='met_'+charin+'_'+typin+'_dirx.map'
rd_im,sx,file='star_'+charin+'_'+typin+'_dirx.map'
window,0,xs=512,ys=512
window,1,xs=512,ys=512
window,2,xs=512,ys=512
window,3,xs=512,ys=512
window,4,xs=512,ys=512

loadct,33
wset,0 & mycontour,dx,/log,/tab,ncont=200
wset,1 & mycontour,tx,/log,/tab,ncont=200
wset,2 & mycontour,zx,/log,/tab,ncont=200
wset,3 & mycontour,bx,/log,/tab,ncont=200,min=1d-4
wset,4 & mycontour,sx,/log,/tab,ncont=200,min=1d-14

read,itest

rd_im,dy,file='dens_'+charin+'_'+typin+'_diry.map'
rd_im,ty,file='temp_'+charin+'_'+typin+'_diry.map'
rd_im,zy,file='met_'+charin+'_'+typin+'_diry.map'
rd_im,sy,file='star_'+charin+'_'+typin+'_diry.map'
rd_im,by,file='blast_'+charin+'_'+typin+'_diry.map'

wset,0 & mycontour,dy,/log,/tab,ncont=200
wset,1 & mycontour,ty,/log,/tab,ncont=200
wset,2 & mycontour,zy,/log,/tab,ncont=200
wset,3 & mycontour,by,/log,/tab,ncont=200,min=1d-4
wset,4 & mycontour,sy,/log,/tab,ncont=200,min=1d-14
read,itest

rd_im,dz,file='dens_'+charin+'_'+typin+'_dirz.map'
rd_im,tz,file='temp_'+charin+'_'+typin+'_dirz.map'
rd_im,zz,file='met_'+charin+'_'+typin+'_dirz.map'
rd_im,sz,file='star_'+charin+'_'+typin+'_dirz.map'
rd_im,bz,file='blast_'+charin+'_'+typin+'_dirz.map'

wset,0 & mycontour,dz,/log,/tab,ncont=200
wset,1 & mycontour,tz,/log,/tab,ncont=200
wset,2 & mycontour,zz,/log,/tab,ncont=200
wset,3 & mycontour,bz,/log,/tab,ncont=200,min=1d-4
wset,4 & mycontour,sz,/log,/tab,ncont=200,min=1d-14

end
