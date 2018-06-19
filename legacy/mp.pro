;----------------------------------------------------------------------
;	Frederic Masset   09/16/96
;	original procedures and idea : Annibal Hetem
;----------------------------------------------------------------------
;----------------------------------------------------------------------
pro majaxis
   common values,nonlyone,nfirst,nlast,firstrow
   common disk, dirname, number
   common planetnb, pn
   comm = 'wc '+dirname+strcompress('/orbit'+string(pn)+'.dat',/remove_all)
   spawn, comm, res
   lines=long(res)
   ll=lines(0)
   m = dblarr(6,ll-1)
   LoadColumn,nfirst,nlast,t,7
   tmin=min(t)/2.0/!PI
   tmax=max(t)/2.0/!PI
   openr,1,filepath(strcompress('orbit'+string(pn)+'.dat',/remove_all),root_dir='.',subdir=dirname)
   readf,1,m
   close,1
   t=m(0,*)/2.0/!PI
   a=m(2,*)
;   persoplot,t,a,xrange=[tmin,tmax],ytitle='!6Semi major-axis',xtitle='!6Time (orbits)',yr=[.9,1.1]
   persoplot,t,a,xrange=[tmin,tmax],ytitle='!6Semi major-axis',xtitle='!6Time (orbits)',/ynozero
   return
end

pro dermajaxis
   common values,nonlyone,nfirst,nlast,firstrow
   common disk, dirname, number
   common planetnb, pn
   comm = 'wc '+dirname+strcompress('/orbit'+string(pn)+'.dat',/remove_all)
   spawn, comm, res
   lines=long(res)
   ll=lines(0)
   m = dblarr(6,ll-1)
   LoadColumn,nfirst,nlast,t,7
   tmin=min(t)/2.0/!PI
   tmax=max(t)/2.0/!PI
   openr,1,filepath(strcompress('orbit'+string(pn)+'.dat',/remove_all),root_dir='.',subdir=dirname)
   readf,1,m
   close,1
   t=m(0,*)
   a=m(2,*)
   da=(a(1:ll-2)-a(0:ll-3))/(t(1:ll-2)-t(0:ll-3))
   tm=.5*(t(1:ll-2)+t(0:ll-3))
;   persoplot,tm/2.0/!PI,alog10(abs(da)+1e-9),xrange=[tmin,tmax],ytitle='!6da/dt',xtitle='!6Time (orbits)',/ynozero ;cuidadin
   persoplot,tm/2.0/!PI,da,xrange=[tmin,tmax],ytitle='!6da/dt',xtitle='!6Time (orbits)',/ynozero ;cuidadin
   return
end

pro eccentricity
   common values,nonlyone,nfirst,nlast,firstrow
   common disk, dirname, number
   common planetnb, pn
   comm = 'wc '+dirname+strcompress('/orbit'+string(pn)+'.dat',/remove_all)
   spawn, comm, res
   lines=long(res)
   ll=lines(0)
   m = dblarr(6,ll-1)
   LoadColumn,nfirst,nlast,t,7
   tmin=min(t)/2.0/!PI
   tmax=max(t)/2.0/!PI
   openr,1,filepath(strcompress('orbit'+string(pn)+'.dat',/remove_all),root_dir='.',subdir=dirname)
   readf,1,m
   close,1
   t=m(0,*)/2.0/!PI
   e=m(1,*)
   persoplot,t,e,xrange=[tmin,tmax],ytitle='!6Eccentricity',xtitle='!6Time (orbits)'
   return
end

pro persotv,image,off1,off2
   common share, numberscreens, indexscreen
   common imagesize, isize
   common tamponvideo, img, x0, y0
   common lastimagesize, sizex, sizey
   common positioning, xll, yll, xur, yur
   infosize = size(image)
   xsize = infosize(1)
   ysize = infosize(2)
   index = indexscreen
   ns = numberscreens
   if (n_elements(off1) EQ 0) then off1 = 0
   if (n_elements(off2) eq 0) then off2 = 0
   temp = image
   temp = congrid(image, isize/numberscreens, isize/numberscreens, cubic=-.5)
   offx = off1/numberscreens
   offx = offx+isize/ns*((index-1) mod ns)
   offy = off2/numberscreens
   offy = offy+isize/ns*(ns-1-floor((index-1)/ns))
   xll = offx
   yll = offy
   xur = offx+xsize/numberscreens
   yur = offy+ysize/numberscreens
   tvscl,temp,offx,offy
   chaine = strcompress ("min : "+string(min(temp))+"; max : "+string(max(temp)))
   txtmsg, chaine
   img = temp
   x0 = offx
   y0 = offy
   sizex = xsize/numberscreens
   sizey = ysize/numberscreens
   indexscreen = indexscreen + 1
   if (indexscreen gt ns * ns) then indexscreen = 1
END
;----------------------------------------------------------------------
pro persopolar,image,theta,rayons
 common values,nonlyone,nfirst,nlast,firstrow
 common zoompolar,xcenter,ycenter,zoomradius
 common share, numberscreens, indexscreen
 ns=numberscreens
 x0 = 1./ns*((indexscreen-1) mod ns)
 y0 = 1./ns*(ns-1-floor((indexscreen-1)/ns))
 x1 = x0+1./ns
 y1 = y0+1./ns
 if ((ns gt 1) and (indexscreen gt 1)) then begin
 polar_contour,image,theta,rayons,/fill,nlev=60,c_colors=indgen(60)*!D.TABLE_SIZE/60,$
               position=[x0,y0,x1,y1],/noerase,xstyle=1,ystyle=1,$
               xrange=[xcenter-zoomradius,xcenter+zoomradius],$
	       yrange=[ycenter-zoomradius,ycenter+zoomradius],$
	       levels=min(image)+(max(image)-min(image))*(findgen(60)/60.)^3
 endif else begin
 polar_contour,image,theta,rayons,/fill,nlev=60,c_colors=indgen(60)*!D.TABLE_SIZE/60,$
               position=[x0,y0,x1,y1],xstyle=1,ystyle=1,$
               xrange=[xcenter-zoomradius,xcenter+zoomradius],$
	       yrange=[ycenter-zoomradius,ycenter+zoomradius],$
	       levels=min(image)+(max(image)-min(image))*(findgen(60)/60.)^3
 endelse
 if (firstrow EQ 0) then indexscreen = indexscreen + 1
 if (indexscreen gt ns * ns) then indexscreen = 1
 return
END
;----------------------------------------------------------------------
pro getgasdim,nr,ns
relat_opnr,2,'dims.dat'
readf,2,a,b,c,d,e,f,nr,ns
close,2
END
;----------------------------------------------------------------------
pro relat_opnru,a,name,err
common disk, dirname, number
nomfichier = strcompress (dirname+name, /REMOVE_ALL)
print,nomfichier
openr,a,filepath(name, SUBDIR=dirname, ROOT_DIR='.'),ERROR=err
end
;----------------------------------------------------------------------
pro relat_opnr,a,name
common disk, dirname, number
nomfichier = strcompress (dirname+name, /REMOVE_ALL)
openr,a,filepath(name, SUBDIR=dirname, ROOT_DIR='.')
end
;-----------------------------------------------------------------------
pro persosmooth,v
common smoothing, width
print,width
if (width LE 1) then return
s = width
if (s GE ((size(v))(1))) then s = ((size(v))(1))
v = smooth(v,s,/edge_truncate)
return
end
;----------------------------------------------------------------------
pro persoplot, a, b, ynozero = ynz, psym=ps, xtitle=xt, ytitle=yt, title=ti,$
               xrange=xr, yrange=yr,charsize=chrsz
  
common values,nonlyone,nfirst,nlast,firstrow
common share, ns, index
if (n_elements(ynz) eq 0) then ynz=0
if (n_elements(ps)  eq 0) then ps= 0
if (n_elements(chrsz)  eq 0) then chrsz=1.0
if (n_elements(xt) eq  0) then xt=''
if (n_elements(yt) eq  0) then yt = ''
if (n_elements(ti) eq  0) then ti = ''
chrsz = chrsz/ns
persosmooth,b
case firstrow of
0 : begin
    xmin = 0.25/ns
    xmax = 0.95/ns
    ymin = 0.15/ns
    ymax = 0.9/ns
    dx = ((index-1) mod ns) * 1.0 / ns
    dy = (ns - 1 - floor((index-1)/ns)) * 1.0 / ns
    xmin = xmin + dx
    xmax = xmax + dx
    ymin = ymin + dy
    ymax = ymax + dy
    if (ns eq 1) then begin
    if (n_elements(xr) ne 0) then begin
    plot,a,b,ynozero=ynz,psym=ps,xtitle=xt,ytitle=yt,title=ti,pos=[xmin, ymin, xmax, ymax],xrange=xr, yrange=yr,charsize=chrsz, xstyle=1, font=17
    end
    if (n_elements(xr) eq 0) then begin
    plot,a,b,ynozero=ynz,psym=ps,xtitle=xt,ytitle=yt,title=ti,pos=[xmin, ymin, xmax, ymax],charsize=chrsz, xstyle=1, font=-1
    end
    end
    if (ns gt 1) then begin
    if (n_elements(xr) ne 0) then begin
    plot,a,b,ynozero=ynz,psym=ps,xtitle=xt,ytitle=yt,title=ti,pos=[xmin, ymin, xmax, ymax], xrange=xr, yrange=yr, /noerase,charsize=chrsz, xstyle=1, font=17
    end
    if (n_elements(xr) eq 0) then begin
    plot,a,b,ynozero=ynz,psym=ps,xtitle=xt,ytitle=yt,title=ti,pos=[xmin, ymin, xmax, ymax], /noerase,charsize=chrsz, xstyle=1, font=-1
    end
    end
    index = index+1
    if (index gt ns * ns) then index = 1
    end
1 : begin
    oplot,a,b,psym=ps
    end
endcase
end
;----------------------------------------------------------------------
;----------------------------------------------------------------------
pro persoplotlog, a, b, ynozero = ynz, psym=ps, xtitle=xt, ytitle=yt, title=ti,$
               xrange=xr, yrange=yr,charsize=chrsz
  
common values,nonlyone,nfirst,nlast,firstrow
common share, ns, index
if (n_elements(ynz) eq 0) then ynz=0
if (n_elements(ps)  eq 0) then ps= 0
if (n_elements(chrsz)  eq 0) then chrsz=1.0
if (n_elements(xt) eq  0) then xt=''
if (n_elements(yt) eq  0) then yt = ''
if (n_elements(ti) eq  0) then ti = ''
chrsz = chrsz/ns
persosmooth,b
case firstrow of
0 : begin
    xmin = 0.25/ns
    xmax = 0.95/ns
    ymin = 0.15/ns
    ymax = 0.9/ns
    dx = ((index-1) mod ns) * 1.0 / ns
    dy = (ns - 1 - floor((index-1)/ns)) * 1.0 / ns
    xmin = xmin + dx
    xmax = xmax + dx
    ymin = ymin + dy
    ymax = ymax + dy
    if (ns eq 1) then begin
    if (n_elements(xr) ne 0) then begin
    plot,a,b,ynozero=ynz,psym=ps,xtitle=xt,ytitle=yt,title=ti,pos=[xmin, ymin, xmax, ymax],xrange=xr, yrange=yr,charsize=chrsz, xstyle=1,/ylog
    end
    if (n_elements(xr) eq 0) then begin
    plot,a,b,ynozero=ynz,psym=ps,xtitle=xt,ytitle=yt,title=ti,pos=[xmin, ymin, xmax, ymax],charsize=chrsz, xstyle=1,/ylog
    end
    end
    if (ns gt 1) then begin
    if (n_elements(xr) ne 0) then begin
    plot,a,b,ynozero=ynz,psym=ps,xtitle=xt,ytitle=yt,title=ti,pos=[xmin, ymin, xmax, ymax], xrange=xr, yrange=yr, /noerase,charsize=chrsz, xstyle=1
    end
    if (n_elements(xr) eq 0) then begin
    plot,a,b,ynozero=ynz,psym=ps,xtitle=xt,ytitle=yt,title=ti,pos=[xmin, ymin, xmax, ymax], /noerase,charsize=chrsz, xstyle=1,/ylog
    end
    end
    index = index+1
    if (index gt ns * ns) then index = 1
    end
1 : begin
    oplot,a,b,psym=ps
    end
endcase
end
;----------------------------------------------------------------------
pro NumberOfSteps,nbsteps
  COMMON disk, dirname, number
  openr,1,filepath('dims.dat',SUBDIR=dirname,ROOT_DIR='.')
  readf,1,a,b,c,d,e,f
  close,1
  nbsteps = f
  end
;----------------------------------------------------------------------
pro getrayons,rmed,nr
  common disk, dirname, number
  radii=fltarr(nr+1)
  rmed = fltarr(nr)
  openr,1,filepath('used_rad.dat',SUBDIR=dirname,ROOT_DIR='.')
  readf,1,radii
  close,1
  rmed=radii(0:nr-1)+radii(1:nr)
  rmed=rmed*0.5
  return
end
;----------------------------------------------------------------------
pro convert_polar, piccart, picpol, nx, ny
  common polmode, polarchoice
  getgasdim,nr,ns
  if (polarchoice EQ 0) then begin
    picpol = piccart
    nx = ns
    ny = nr
    return
  endif
  image = piccart
  getrayons,rmed,nr
  theta=findgen(ns)/ns*2.0*!PI
  rmed=extrac(rmed,-1,nr+1)
  rmed(0)=2.0*rmed(1)-rmed(2)
  image=extrac(image,0,-1,ns,nr+1)
  persopolar,image,theta,rmed
  return
end
;-----------------------------------------------------------------------
pro gasvxvy, testnumber
  common values,nonlyone,nfirst,nlast,firstrow
  common disk, dirname, number
  common imagesize, isize
  common polmode, polarchoice
  common zoompolar,xcenter,ycenter,zoomradius
  common share, numberscreens, indexscreen
  ns=numberscreens
  x0 = 1./ns*((indexscreen-1) mod ns)
  y0 = 1./ns*(ns-1-floor((indexscreen-1)/ns))
  x1 = x0+1./ns
  y1 = y0+1./ns
  getgasdim,nr,ns
  vvr=dblarr(ns,nr)  ;cuidadin
  vvt=dblarr(ns,nr)  ;cuidadin
  filename = strcompress ('gasvtheta' $
		  + strtrim(string(nonlyone),2) + '.dat',/remove_all)
  openr,1,filepath(filename,SUBDIR=dirname,ROOT_DIR='.')
  readu,1,vvt
  close,1
  filename = strcompress ('gasvrad' $
		  + strtrim(string(nonlyone),2) + '.dat',/remove_all)
  openr,1,filepath(filename,SUBDIR=dirname,ROOT_DIR='.')
  readu,1,vvr
  close,1
  longueur = 1.0
  if (testnumber EQ 1) then begin
    filename = 'gasvtheta0.dat'
    openr,1,filepath(filename,SUBDIR=dirname,ROOT_DIR='.')
    vvt0=vvt
    readu,1,vvt0
    close,1
    filename = 'gasvrad0.dat'
    openr,1,filepath(filename,SUBDIR=dirname,ROOT_DIR='.')
    vvr0=vvr
    readu,1,vvr0
    close,1
    ;vvr=vvr-vvr0
    ;vvt=vvt-vvt0
  endif
  if (testnumber EQ 2) then begin
    longueur = 2.0
    timesteplookup,nonlyone,x,y,vx,vy,mass
    r=sqrt(x*x+y*y)
    vrp=(vx*x+vy*y)/r
    vtp=(-vx*y+vy*x)/r
    ;vvr=vvr-vrp
    ; NEW PART BELOW THIS COMMENT FM 020399
    getrayons,rmed,nr
    radius=rmed##replicate(1,ns)
    vttp = radius*vtp/r;
    ;vvt=vvt-vttp
    ; NEW PART ABOVE THIS COMMENT FM 020399
  endif
  getrayons,rmed,nr
  if polarchoice EQ 0 then begin
	vxc = congrid(vvt,31,31)
	vyc = congrid(vvr,31,31)
	xc  = findgen(31)/31.0*2.0*!PI
	yc  = findgen(31)/30.0*(rmed(nr-1)-rmed(0))+rmed(0)
	if ((indexscreen EQ 1) AND (firstrow EQ 0)) then begin
  		velovect,vxc,vyc,xc,yc,position=[x0,y0,x1,y1],length=longueur
	endif else begin
  		velovect,vxc,vyc,xc,yc,position=[x0,y0,x1,y1],/noerase,length=longueur
	endelse
  	indexscreen = indexscreen + 1
  	if (indexscreen gt numberscreens * numberscreens) then indexscreen = 1
  	return
  endif else begin
	vvr=extrac(vvr,0,0,ns,nr+1)
	vvr=0.5*(vvr(*,0:nr-1)+vvr(*,1:nr))
	vvt=0.5*(vvt+shift(vvt,-1,0))
	n_arrows=15
	angle=replicate(1,nr)##findgen(ns)/ns*2.0*!PI
	radius=rmed##replicate(1,ns)
	xx=radius*cos(angle)
	yy=radius*sin(angle)
	vx=vvr*cos(angle)-vvt*sin(angle)
	vy=vvr*sin(angle)+vvt*cos(angle)
	triangulate,xx,yy,tr,b
	vxc=trigrid(xx,yy,vx,tr,[zoomradius/n_arrows,zoomradius/n_arrows],$
	    [xcenter-zoomradius,ycenter-zoomradius,$
	     xcenter+zoomradius,ycenter+zoomradius])
	vyc=trigrid(xx,yy,vy,tr,[zoomradius/n_arrows,zoomradius/n_arrows],$
	    [xcenter-zoomradius,ycenter-zoomradius,$
	     xcenter+zoomradius,ycenter+zoomradius])
	xc = (findgen(2*n_arrows+1)/(n_arrows)-1.0)*zoomradius+xcenter
	yc = (findgen(2*n_arrows+1)/(n_arrows)-1.0)*zoomradius+ycenter
  	if ((indexscreen EQ 1) AND (firstrow EQ 0)) then begin
  		velovect,vxc,vyc,xc,yc,position=[x0,y0,x1,y1],$
		xrange=[xcenter-zoomradius,xcenter+zoomradius],$
		yrange=[ycenter-zoomradius,ycenter+zoomradius],length=longueur
  	endif else begin
  		velovect,vxc,vyc,xc,yc,position=[x0,y0,x1,y1],/noerase,$
		xrange=[xcenter-zoomradius,xcenter+zoomradius],$
		yrange=[ycenter-zoomradius,ycenter+zoomradius],length=longueur
  	endelse
  	indexscreen = indexscreen + 1
  	if (indexscreen gt numberscreens * numberscreens) then indexscreen = 1
  endelse
  return
END
;-----------------------------------------------------------------------
pro txtmsg,txt
  COMMON WTEXT4_Comm,WTEXT4_Id
  COMMON DRAW27_Comm,DRAW27_Id
  print,txt
  widget_control,WTEXT4_Id, set_value=string(txt)
  wset, draw27_id
end
;-----------------------------------------------------------------------
pro gasprofile,name,titrey
  common values,nonlyone,nfirst,nlast,firstrow
  COMMON disk, dirname, number
  fname = strcompress(name+strtrim(string(nonlyone),2)+'.dat',/remove_all)
  fname0 = strcompress(name+strtrim(string(0),2)+'.dat',/remove_all)
  getgasdim,nr,ns
  print,ns,nr
  f=dblarr(ns,nr)  ;cuidadin
  f0=f
  openr,1,filepath(fname,SUBDIR=dirname,ROOT_DIR='.')
  readu,1,f
  close,1
;  openr,1,filepath(fname0,SUBDIR=dirname,ROOT_DIR='.')
;  readu,1,f0
;  close,1
;  if (name EQ 'vtheta') then begin
     ; f=f-f0
;  endif
  rmed = fltarr(nr+1)
  getrayons,rmed,nr
  ;planetsector,nonlyone,plsec ; we do the following to avoid the planet
  plsec = 0.0
  f=shift(f,ns/4-plsec,0)
  redns = long(ns/2)
  redns=redns*2
  ;tvscl,f
;  f=f(redns:2*redns-1,*)
  ff=total(f,1)
  ;tvscl,f,200,0
  ff = ff/redns
  ;ff = f(200,*)
;  persoplotlog,rmed,ff,ynozero=1,xtitle='Radius',ytitle=titrey
 persoplot,rmed,ff,ynozero=1,xtitle='Radius',ytitle=titrey
;  persoplot,rmed,ff*rmed^1.5,ynozero=0,xtitle='Radius',ytitle=titrey
end
;-----------------------------------------------------------------------  
pro plot1
  common values,nonlyone,nfirst,nlast,firstrow
  common plot1d, nname1
  if nname1 EQ 0 then begin
    gasprofile, 'gasdens', 'Density'
    return
  endif
  if nname1 EQ 1 then begin
    gasprofile, 'gasvrad', 'Radial Velocity'
    return
  endif
  if nname1 EQ 2 then begin
    gasprofile, 'gasvtheta', 'Rotationnal Velocity'
    return
  endif
end
;-----------------------------------------------------------------------
pro namegasfile,nn,fname
  common values,nonlyone,nfirst,nlast,firstrow
  common gas, gaschoice
  namelist=['dens','dens','dens','dens','dens','vrad','vtheta','label']
  nameprefix = namelist(gaschoice)
  fname = strcompress('gas'+nameprefix+ strtrim(string(nn),2) + '.dat'$
  ,/remove_all)
end
;-----------------------------------------------------------------------
pro readgasfile,nomfich,y,nr,ns
  common values,nonlyone,nfirst,nlast,firstrow
  getgasdim,nr,ns
  print,'reading ',nomfich
  print,ns,nr
  y=dblarr(ns,nr)   ;cuidadin
  relat_opnru,1,nomfich,err
  if (err NE 0) then begin
    txtmsg, '*** Not available'
    return
  endif
  readu,1,y
  close,1
end
;-----------------------------------------------------------------------
pro animgas
  common values,nonlyone,nfirst,nlast,firstrow
  common imagesize,isize
  common gas, gaschoice
  common animframe, frameanim
  oldsize = isize
  isize = 193
  max = 0.0
  getgasdim,nr,ns
  nnr = round(ns/nr)*nr
  nfile = nlast-nfirst+1
;  spawn, 'rm -fR Movies/*'
  XINTERANIMATE, SET=[ns, nnr, nfile], /SHOWLOAD
  for i=nfirst,nlast do begin
    namegasfile,i,filename
    readgasfile,filename,yy,nr,ns
    IF (gaschoice EQ 1) then yy = alog(yy)-alog(min(yy))
    ;if (i EQ nfirst) then 
    max = max(yy)
    ;max = 0.001
    planetsector,i,j
    if (frameanim EQ 0) then begin
    	    yy = shift(yy,ns/2-j,0)
    endif
    Result = BYTSCL(yy,TOP=!D.TABLE_SIZE, MAX=max, min=0.0)
    picturetoplot = rebin(Result,ns,nnr)
    pp = picturetoplot
;    for j=0,ns-1 do begin
;        vec1 = pp(j,*)
;        vec2 = interpol(vec1,.1*(25.^(dindgen(nnr)/nnr)),.1+1.*dindgen(nnr)/nnr)
;        picturetoplot(j,*) = vec2
;    endfor
    ;XINTERANIMATE, FRAME = i-nfirst, image = congrid(Result,ns,nnr)
    XINTERANIMATE, FRAME = i-nfirst, image = picturetoplot
;    interstr = strcompress(string(i-nfirst),/remove_all)
;    astr = '000000'
;    strput, astr, interstr, 6-strlen(interstr)
;    filemov = "pic"+astr+".gif"
;    print,filemov
;    filetot = strcompress("Movies/"+filemov,/remove_all)
;    write_gif, filetot, picturetoplot
  endfor
  XINTERANIMATE
;  spawn, 'tar cvf movie.tar Movies/*.gif'
;  spawn, 'sstar'
  isize = oldsize
end
;-----------------------------------------------------------------------
PRO curseur
  common zoompolar,xcenter,ycenter,zoomradius
  cursor, x, y, 3, /data
  xcenter=x
  ycenter=y
  str = string(x,/print)+string(y,/print)
  txtmsg, str
END
;-----------------------------------------------------------------------
pro timesteplookup,timestep,x,y,vx,vy,mass
  COMMON disk, dirname, number
  common planetnb, pn
  common suffix,extension
  ncolumn = 8
  testvector = fltarr(ncolumn)
  testvector(0)=-1
  openr,1,filepath(strcompress('planet'+string(pn),/remove_all)+extension,SUBDIR=dirname,ROOT_DIR='.')
  while ((not eof(1)) and (testvector(0) ne  timestep)) do begin
    readf,1,testvector
  endwhile
  close,1
  x=testvector(1)
  y=testvector(2)
  vx=testvector(3)
  vy=testvector(4)
  mass=testvector(5)
  return
end
;-----------------------------------------------------------------------
pro PlanetSector,timestep,j
  getgasdim,nr,ns
  timesteplookup,timestep,x,y,vx,vy,mass
  j = atan(y,x)/2.0/!PI*ns
  return
end
;-----------------------------------------------------------------------
PRO hardcopy
  COMMON DRAW27_Comm, DRAW27_Id
  wset, draw27_id
  bb=tvrd ()
  set_plot,'ps'
  device,filename='protodisk.ps',/landscape,/color,$
	xsize=19,ysize=19,xoffset=1,yoffset=24.5
  tv,bb
  device, /close
  set_plot, 'X'
end
;-----------------------------------------------------------------------
PRO parameters
  COMMON disk, dirname, number
  nom = strcompress ('disc'+number+'.par', /REMOVE_ALL)
  titre='Parameters for disc number '+number
  xdisplayfile, filepath(nom, SUBDIR='.', ROOT_DIR='.'), title=titre
END
;-----------------------------------------------------------------------
PRO DEP3_Event, Event
common disk, dirname, number
common toggle, prefix
  WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev,GET_VALUE=dir_name
  CASE Ev OF 
  'FIELD4': BEGIN
    number = string(dir_name, /print)
    dirname = strcompress (prefix+'out'+dir_name)
    txtmsg, dirname
    widget_control, Event.top, /DESTROY
    END
  ENDCASE
END
;-----------------------------------------------------------------------
pro visugas, image_number
  common values,nonlyone,nfirst,nlast,firstrow
  common disk, dirname, number
  common imagesize, isize
  common polmode, polarchoice
  common representation, mode
  common share, numberscreens, indexscreen
  namegasfile,image_number,filename
  readgasfile,filename,yy,nr,ns
  convert_polar,yy,y,nx,ny
  if (polarchoice EQ 1) then return
  CASE mode OF
    0: BEGIN
      temp=y
      l1 = nx
      while ((l1+nx) LE isize) do l1 = l1 + nx
      l2 = ny
      while ((l2+ny) LE isize) do l2 = l2 + ny
      print,nx,ny,l1,l2
      persotv,temp
    END
    1: BEGIN
      temp = y
      shade_surf,temp,charsize=2.5
    END
  END
end
;-----------------------------------------------------------------------
pro WriteAGIFFile
  common values,nonlyone,nfirst,nlast,firstrow
  common disk, dirname, number
  common imagesize, isize
  common polmode, polarchoice
  common representation, mode
  common share, numberscreens, indexscreen
  namegasfile,nonlyone,filename
  readgasfile,filename,yy,nr,ns
  write_gif, 'mp.gif', bytscl(congrid(yy,384,384,/cubic))
end
;-----------------------------------------------------------------------
pro spectrum
  common wavenumber, wnb
  common positioning, xll, yll, xur, yur
  common disk, dirname, number
  common representation, mode
  common spectratype, sptype
  common spectraname, spname
  common zoomspectra, zoomsp
  fname=strcompress(spname+string(wnb)+".dat",/remove_all)
  relat_opnr,1,fname
  x=dblarr(2)
  i=lonarr(2)
  readu,1,i
  s=dblarr(i(0),i(1))
  readu,1,x
  print,i
  print,x
  readu,1,s
  close,1
  GetCorrectGrid, jfrac
  s= shift(s,0,wnb*(1.0/x(1.0))+0.5);cuidadin : for rotating frame only
  s = interpolate (s, jfrac, findgen(i(1))/zoomsp,/grid)
  if (sptype EQ 1) then s = alog(abs(s))
  if (sptype EQ 2) then s = sqrt(abs(s))
  if (mode EQ 1) then begin
	shade_surf,s
	return
  endif
  ;s= shift(s,0,i(1)/2);cuidadin 
  s = congrid (s, 680, 680, cubic=-0.5)
  persotv,s,50,50
  radii=fltarr(i(0)+1)
  openr,1,filepath('used_rad.dat',SUBDIR=dirname,ROOT_DIR='.')
  readf,1,radii
  close,1
  rmin = radii(0)
  rmax = radii(i(0))
  ray=findgen(1000)/1000.0*(rmax-rmin)+rmin
  omega=(ray+0.04)^(-1.5);cuidadin GFrame only
  kappa=omega
  plot,xrange=[rmin,rmax],yrange=[-x(1)/2.0/zoomsp,-x(1)/2.0/zoomsp+i(1)*x(1)/zoomsp],$
  ;plot,xrange=[rmin,rmax],yrange=[-x(1)/2.0-i(1)*x(1)/2.0,-x(1)/2.0+i(1)*x(1)-i(1)*x(1)/2.0],$
       findgen(10)*0.0,xstyle=1,ystyle=1,$
       xtitle="!6Radius", ytitle="!6Frequency !7x!6",$
       position=[xll,yll,xur,yur],/device,/noerase
  oplot,ray,omega*wnb
  oplot,ray,omega*wnb+kappa
  oplot,ray,omega*wnb-kappa
  ;oplot,ray,omega*(wnb+2),linestyle=1
  ;oplot,ray,omega*(wnb-2),linestyle=1
  return
  end
;-----------------------------------------------------------------------
pro spectrumcut
  common wavenumber, wnb
  common positioning, xll, yll, xur, yur
  common disk, dirname, number
  common representation, mode
  common spectratype, sptype
  common spectraname, spname
  common slider, fine, coarse
  fname=strcompress(spname+string(wnb)+".dat",/remove_all)
  relat_opnr,1,fname
  x=dblarr(2)
  i=lonarr(2)
  readu,1,i
  s=dblarr(i(0),i(1))
  readu,1,x
  readu,1,s
  close,1
  s= shift(s,0,wnb*(1.0/x(1.0))+0.5);cuidadin : for rotating frame only
  getgasdim,nr,ns
  getrayons,rmed,nr
  jfrac = 0.0

    r = fine/(100.0)*(rmed(nr-1)-rmed(0))+rmed(0)
    j = 0
    while (rmed(j) LT r) do j = j+1
    if (j EQ 0) then jfrac = 0
    if (j GT 0) then jfrac = j-1+(r-rmed(j-1))/(rmed(j)-rmed(j-1))

  s= shift(s,0,i(1)/2);cuidadin 

  s = interpolate (s, findgen(i(1))*0.0+jfrac, findgen(i(1)))
  if (sptype EQ 1) then s = alog(abs(s))
  if (sptype EQ 2) then s = sqrt(abs(s))
  persoplot,-x(1)*0.5-i(1)*x(1)/2.0+findgen(i(1))*x(1),s ;cuidadin rot.
;  persoplot,-x(1)*0.5+findgen(i(1))*x(1),s
  oplot,findgen(100)*0.0,sin(findgen(100)/100)*100.0,linestyle=1
  return
  end
;-----------------------------------------------------------------------
pro visuloggas, image_number
  common values,nonlyone,nfirst,nlast,firstrow
  common disk, dirname, number
  common representation, mode
  common imagesize, isize
  common polmode, polarchoice
  common share, numberscreens, indexscreen
  namegasfile,image_number,filename
  readgasfile,filename,yy,nr,ns
  yy = alog10 (yy)
  yy = yy-min(yy)
  convert_polar,yy,y,nx,ny
  if (polarchoice EQ 1) then return
  CASE mode OF
    0: BEGIN
      temp=y
      l1 = nx
      while ((l1+nx) LE isize) do l1 = l1 + nx
      l2 = ny
      while ((l2+ny) LE isize) do l2 = l2 + ny
      print,nx,ny,l1,l2
      persotv,temp
       END
    1: BEGIN
      temp = y
      shade_surf,temp,charsize=2.5
       END
   END
end
;-----------------------------------------------------------------------
PRO GasMenu, number
common gas, gaschoice
  CASE gaschoice OF
  0: BEGIN
	visugas, number
     END
  1: BEGIN
	visuloggas, number
     END
  2: BEGIN
        gasvxvy, 0
     END
  3: BEGIN
        gasvxvy, 1
     END
  4: BEGIN
        gasvxvy, 2
     END
  5: BEGIN
        visugas, number
     END
  6: BEGIN
        visugas, number
     END
  7: BEGIN
	visugas, number
     END
  ENDCASE
END
;-----------------------------------------------------------------------
PRO PDMENUShare_Event, Event
common share, numberscreens, indexscreen
COMMON DRAW27_Comm,DRAW27_Id
  CASE Event.Value OF 
  'Screen Share.Only one screen': BEGIN
	numberscreens = 1
	indexscreen = 1
    END
  'Screen Share.2x2 screen': BEGIN
	numberscreens = 2
	indexscreen = 1
    END
  'Screen Share.3x3 screen': BEGIN
	numberscreens = 3
	indexscreen = 1
    END
  'Screen Share.4x4 screen': BEGIN
	numberscreens = 4
	indexscreen = 1
    END
  ENDCASE
  wset, DRAW27_Id
  erase
  !p.charsize = 1.0/numberscreens
END
;-----------------------------------------------------------------------
pro LoadColumn,min,max,vector,columnnumber
  COMMON disk, dirname, number
  common planetnb, pn
  common suffix,extension
  ncolumn = 8
  testvector = fltarr(ncolumn)
  vector=fltarr(max-min+1)
  openr,1,filepath(strcompress('planet'+string(pn),/remove_all)+extension,SUBDIR=dirname,ROOT_DIR='.')
  while (not eof(1)) do begin
     readf,1,testvector
     if ((testvector(0) ge min) and (testvector(0) le max)) then begin
	vector(testvector(0)-min) = testvector(columnnumber)
     endif
  endwhile
  close,1
  return
end
;-----------------------------------------------------------------------
pro LoadMomCol,min,max,vector,columnnumber
  COMMON disk, dirname, number
  common suffix,extension
  ncolumn = 8
  testvector = fltarr(ncolumn)
  vector=fltarr(max-min)
  openr,1,filepath('momentum'+extension,SUBDIR=dirname,ROOT_DIR='.') 
  while (not eof(1)) do begin
     readf,1,testvector
     if ((testvector(0) gt min) and (testvector(0) le max)) then begin
	vector(testvector(0)-min-1) = testvector(columnnumber)
     endif
  endwhile
  close,1
  return
end
;-----------------------------------------------------------------------
pro Masscurve
  common values,nonlyone,nfirst,nlast,firstrow
  LoadColumn,nfirst,nlast,mass,5
  LoadColumn,nfirst,nlast,t,7
  persoplot,t/2.0/!PI,mass-mass(0),xtitle="Time (Orbits)",ytitle="Accreted Mass",title="Accreted Mass",charsize=1.0
  return
end
;-----------------------------------------------------------------------
pro LostMasscurve
  common values,nonlyone,nfirst,nlast,firstrow
  LoadColumn,nfirst,nlast,mass,6
  LoadColumn,nfirst,nlast,t,7
  persoplot,t/2.0/!PI,mass,xtitle="Orbits",ytitle="Disk lost mass at inner boundary",charsize=1.0
  return
end
;-----------------------------------------------------------------------
pro Accretioncurve
  common values,nonlyone,nfirst,nlast,firstrow
  LoadColumn,nfirst,nlast,mass,5
  LoadColumn,nfirst,nlast,t,7
  n=nlast-nfirst
  accrate = (mass(1:n)-mass(0:n-1))/(t(1:n)-t(0:n-1))
  t = 0.5*(t(1:n)+t(0:n-1))
;  persoplot,t/2.0/!PI,accrate*2.0*!PI,xtitle="Time (Orbits)",ytitle="Accretion Rate/Orbit", title="Mass Accretion Rate", charsize=1.0
  plot,t/2.0/!PI,accrate*2.0*!PI,xtitle="Time (Orbits)",ytitle="Accretion Rate/Orbit", title="Mass Accretion Rate", charsize=1.0,/ylog
  return
end
;-----------------------------------------------------------------------
pro Migration
  common values,nonlyone,nfirst,nlast,firstrow
  LoadColumn,nfirst,nlast,x,1
  LoadColumn,nfirst,nlast,y,2
  LoadColumn,nfirst,nlast,t,7
  persoplot,t/2.0/!PI,sqrt(x*x+y*y),xtitle="Time (Orbits)",ytitle="Distance",title="Planet - Star Separation", charsize=1.0,/ynozero
  return
end
;-----------------------------------------------------------------------
pro RocheLobe
  common values,nonlyone,nfirst,nlast,firstrow
  common polmode, polarchoice
  if (polarchoice EQ 0) then begin
	txtmsg, "Only in x-y polar mode"
	return
  end
  timesteplookup,nonlyone,x,y,vx,vy,mass
  print,x,y
  rroche=sqrt(x*x+y*y)*(mass/1.0/3.0)^(1.0/3.0)
  drawcircle, x, y, rroche, 0
  return
end
;-----------------------------------------------------------------------
pro drawcircle, xc, yc, radius, linestyle
  tt=findgen(10000)
  xx=xc+radius*cos(2.0*!PI*tt/10000.0)
  yy=yc+radius*sin(2.0*!PI*tt/10000.0)
  oplot,xx,yy,linestyle=linestyle
  return
end
;-----------------------------------------------------------------------
pro resonances
 common values,nonlyone,nfirst,nlast,firstrow
 timesteplookup,nonlyone,x,y,vx,vy,mass
 xstar = -x*mass ;only
 ystar = -y*mass ;approximate
 xstar = 0.0
 ystar = 0.0
 print,xstar,ystar
 dist = (x-xstar)*(x-xstar)+(y-ystar)*(y-ystar)
 dist = sqrt(dist)
 drawcircle,xstar,ystar,dist*(3.0/1.0)^(2.0/3.0),1
 drawcircle,xstar,ystar,dist*(4.0/2.0)^(2.0/3.0),1
 drawcircle,xstar,ystar,dist*(5.0/3.0)^(2.0/3.0),1
 drawcircle,xstar,ystar,dist*(2.0/1.0)^(2.0/3.0),0
 drawcircle,xstar,ystar,dist*(3.0/2.0)^(2.0/3.0),0
 drawcircle,xstar,ystar,dist*(4.0/3.0)^(2.0/3.0),0
return
end
;-----------------------------------------------------------------------
pro readvalue
  common tamponvideo, image, x0, y0
  COMMON DRAW27_Comm,DRAW27_Id
  wset, draw27_Id
  rdpix, image, x0, y0
end
;-----------------------------------------------------------------------
pro gaszoomplus
  common values,nonlyone,nfirst,nlast,firstrow
  common zoompolar,xcenter,ycenter,zoomradius
  zoomradius = zoomradius / 2.0
  txtmsg, 'Rmax : '+strtrim(string(zoomradius))
  gasmenu, nonlyone
end
;-----------------------------------------------------------------------
pro gaszoommoins
  common values,nonlyone,nfirst,nlast,firstrow
  common zoompolar,xcenter,ycenter,zoomradius
  zoomradius = zoomradius * 2.0
  txtmsg, 'Rmax : '+strtrim(string(zoomradius))
  gasmenu, nonlyone
end	
;-----------------------------------------------------------------------
pro gasdefault
  common values,nonlyone,nfirst,nlast,firstrow
  common zoompolar,xcenter,ycenter,zoomradius
  zoomradius = 3.0
  xcenter = 0.0
  ycenter = 0.0
  txtmsg, 'Rmax : '+strtrim(string(zoomradius))
  gasmenu, nonlyone
end	
;-----------------------------------------------------------------------
pro PlotDetail
  common values,nonlyone,nfirst,nlast,firstrow
  common zoompolar,xcenter,ycenter,zoomradius
  common polmode, polarchoice
  common modedetail, modet
  oldfirstrow = firstrow
  timesteplookup,nonlyone,x,y,vx,vy,mass
  rroche=sqrt(x*x+y*y)*(mass/1.0/3.0)^(1.0/3.0)
  xcenter = x
  ycenter = y
  zoomradius = rroche*2.5
  firstrow = 1  
  polarchoice = 1
  gasmenu,nonlyone
  gasvxvy,modet
  RocheLobe
  plotgrid
  firstrow=oldfirstrow
  return
end
;-----------------------------------------------------------------------
pro getcorrectgrid,jfrac
  getgasdim,nr,ns
  getrayons,rmed,nr
  jfrac = fltarr(nr)
  for i=0, nr-1 do begin
    r = i/(nr-1)*(rmed(nr-1)-rmed(0))+rmed(0)
    j = 0
    while (rmed(j) LT r) do j = j+1
    if (j EQ 0) then jfrac(i) = 0
    if (j GT 0) then jfrac(i) = j-1+(r-rmed(j-1))/(rmed(j)-rmed(j-1))
  endfor
  return
end
;-----------------------------------------------------------------------
pro plotgrid
  common disk, dirname, number
  getgasdim,nr,ns
  radii=fltarr(nr+1)
  openr,1,filepath('used_rad.dat',SUBDIR=dirname,ROOT_DIR='.')
  readf,1,radii
  close,1
  angle=findgen(3000)/3000.0*6.29
  cost = cos(angle)
  sint = sin(angle)
  for i=0, nr do begin
    xc=radii(i)*cost
    yc=radii(i)*sint
    oplot,xc,yc
  endfor
  for i=0, ns do begin
    a=2.0*i/ns*!PI+!PI/ns
    x0=radii(0)*cos(a)
    y0=radii(0)*sin(a)
    x1=radii(nr)*cos(a)
    y1=radii(nr)*sin(a)
    oplot,[x0,x1],[y0,y1]
  endfor
  return
end
;-----------------------------------------------------------------------
PRO PDMENU2_Event, Event
  common values,nonlyone,nfirst,nlast,firstrow
  LoadColumn,nfirst,nlast,t,7
  LoadColumn,nfirst,nlast,x,1
  LoadColumn,nfirst,nlast,y,2
  LoadColumn,nfirst,nlast,vx,3
  LoadColumn,nfirst,nlast,vy,4
  LoadColumn,nfirst,nlast,mp,5
  sig = mp*(vy*x-vx*y)
  dsig = (sig(1:nlast-nfirst)-sig(0:nlast-nfirst-1))
  dsig = dsig/(t(1:nlast-nfirst)-t(0:nlast-nfirst-1))
  halft = 0.5*(t(1:nlast-nfirst)+t(0:nlast-nfirst-1))
  LoadMomCol,nfirst,nlast,tt,1
  LoadMomCol,nfirst,nlast,dtmomacc,2
  LoadMomCol,nfirst,nlast,torqrings,3
  LoadMomCol,nfirst,nlast,torqinner,4
  LoadMomCol,nfirst,nlast,torqouter,5
  LoadMomCol,nfirst,nlast,torqorings,6
  LoadMomCol,nfirst,nlast,indirect,7
  CASE Event.Value OF 
  'Momemtum.Momentum(t)': BEGIN
    persoplot,t/2.0/!PI,sig,xtitle="Orbits",ytitle="Protoplanet Angular Momentum",charsize=2.0
    END
  'Momemtum.Derivatives.dMomentum/dt': BEGIN
    persoplot,halft/2.0/!PI,dsig,xtitle="Orbits",ytitle="dJ_planet/dt",charsize=2.0
    END
  'Momemtum.Derivatives.dMom/dt | Accretion': BEGIN
    persoplot,tt/2.0/!PI,dtmomacc,xtitle="Orbits",ytitle="dJ_planet/dt from accretion",charsize=2.0
    END
  'Momemtum.Derivatives.difference': BEGIN
    persoplot,tt/2.0/!PI,dsig-dtmomacc,xtitle="Orbits",ytitle="d(J-Jacc)_planet/dt from accretion",charsize=2.0
    END
  'Momemtum.Torques.Inner Disk': BEGIN
    persoplot,tt/2.0/!PI,torqinner,xtitle="Orbits",ytitle="Torque from inner disk",charsize=2.0
    END
  'Momemtum.Torques.Outer Disk': BEGIN
    persoplot,tt/2.0/!PI,torqouter,xtitle="Orbits",ytitle="Torque from outer disk",charsize=2.0,$
    yrange=[-max(abs(torqouter)),max(abs(torqouter))],xrange=[min(tt/2.0/!PI),max(tt/2.0/!PI)]
    END
  'Momemtum.Torques.Inner Rings': BEGIN
    persoplot,tt/2.0/!PI,torqrings,xtitle="Orbits",ytitle="Torque from inner rings",charsize=2.0
    END
  'Momemtum.Torques.Outer Rings': BEGIN
    persoplot,tt/2.0/!PI,torqorings,xtitle="Orbits",ytitle="Torque from outer rings",charsize=2.0
    END
  'Momemtum.Torques.Indirect': BEGIN
    persoplot,tt/2.0/!PI,indirect,xtitle="Orbits",ytitle="Torque from outer rings",charsize=2.0
    END
  'Momemtum.Torques.Total': BEGIN
    persoplot,tt/2.0/!PI,torqouter+torqinner,xtitle="Orbits",ytitle="Torque from disk",charsize=2.0
    END
  'Momemtum.Balance Residual': BEGIN
    persoplot,tt/2.0/!PI,dsig-dtmomacc-torqouter-torqinner-indirect,xtitle="Orbits",ytitle="Balance residual from angular momentum conservation",charsize=2.0
    END
  'Momemtum.Partial Balance': BEGIN
    persoplot,tt/2.0/!PI,dsig-dtmomacc-torqouter-torqinner,xtitle="Orbits",ytitle="Balance residual from angular momentum conservation",charsize=2.0
    END
  ENDCASE
END
;-----------------------------------------------------------------------
pro TroncateIndex, rayonmax, index
  getgasdim,nr,ns
  rmed = fltarr(nr+1)
  getrayons, rmed, nr
  j = 0
  for i=0, nr-1 do begin
    if (rmed(i) LT rayonmax) then j=i
  endfor
  if (j LT 3) then j=3
  index = j
return
end
;-----------------------------------------------------------------------
pro zoomdisplay
  COMMON DRAW27_Comm,DRAW27_Id
  wset, draw27_Id
  zoom
  end
;-----------------------------------------------------------------------
pro panels
set_plot,'ps'
device,filename='panels.ps'
AccretionCurve
Migration
MassCurve
device,/close
set_plot,'X'
return
end
;-----------------------------------------------------------------------
PRO PDMENUTRICKS_Event, Event
  common suffix,extension
  common modedetail, modet
  common wavenumber, wnb
  common spectratype, sptype
  common spectraname, spname
  common zoomspectra, zoomsp
  common animframe, frameanim
  common toggle, prefix
  common disk, dirname, number

  CASE Event.Value OF 
  'Tricks.Close Files': BEGIN
    close,1
    set_plot,'x'
    END
  'Tricks.Velocity field.V-VKep': BEGIN
    modet = 1
    END
  'Tricks.Velocity field.Omega-Omega_planet': BEGIN
    modet = 2
    END
  'Tricks.Spectra Aspect.Linear': BEGIN
    sptype = 0
    END
  'Tricks.Spectra Aspect.Log': BEGIN
    sptype = 1
    END
  'Tricks.Spectra Aspect.Sqrt': BEGIN
    sptype = 2
    END
  'Tricks.Spectra Zoom.1x': BEGIN
    zoomsp = 1.0
    END
  'Tricks.Spectra Zoom.2x': BEGIN
    zoomsp = 2.0
    END
  'Tricks.Spectra Zoom.4x': BEGIN
    zoomsp = 4.0
    END
  'Tricks.Spectra Zoom.8x': BEGIN
    zoomsp = 8.0
    END
  'Tricks.Spectra Zoom.16x': BEGIN
    zoomsp = 16.0
    END
  'Tricks.Spectra Type.Radius/Frequency': BEGIN
    spname = "spec"
    END
  'Tricks.Spectra Type.Radius/Time (cos)': BEGIN
    spname = "cos"
    END
  'Tricks.Spectra Type.Radius/Time (sin)': BEGIN
    spname = "sin"
    END
  'Tricks.extensions..dat': BEGIN
    extension='.dat'
    END
  'Tricks.extensions..old1': BEGIN
    extension='.old1'
    END
  'Tricks.extensions..old2': BEGIN
    extension='.old2'
    END
  'Tricks.extensions..old3': BEGIN
    extension='.old3'
    END
  'Tricks.extensions..old4': BEGIN
    extension='.old4'
    END
  'Tricks.extensions..old5': BEGIN
    extension='.old5'
    END
  'Tricks.extensions..old6': BEGIN
    extension='.old6'
    END
  'Tricks.extensions..old7': BEGIN
    extension='.old7'
    END
  'Tricks.extensions..old8': BEGIN
    extension='.old8'
    END
  'Tricks.extensions..old9': BEGIN
    extension='.old9'
    END
  'Tricks.Animations.Rotating frame': BEGIN
    frameanim = 0
    END
  'Tricks.Animations.Fix frame': BEGIN
    frameanim = 1
    END
  'Tricks.Wave Number.m = 0': BEGIN
    wnb = 0
    END
  'Tricks.Wave Number.m = 1': BEGIN
    wnb = 1
    END
  'Tricks.Wave Number.m = 2': BEGIN
    wnb = 2
    END
  'Tricks.Wave Number.m = 3': BEGIN
    wnb = 3
    END
  'Tricks.Wave Number.m = 4': BEGIN
    wnb = 4
    END
  'Tricks.Wave Number.m = 5': BEGIN
    wnb = 5
    END
  'Tricks.Wave Number.m = 6': BEGIN
    wnb = 6
    END
  'Tricks.Wave Number.m = 7': BEGIN
    wnb = 7
    END
  'Tricks.Wave Number.m = 8': BEGIN
    wnb = 8
    END
  'Tricks.Wave Number.m = 9': BEGIN
    wnb = 9
    END
  'Tricks.Resonances': BEGIN
    resonances
    END
  'Tricks.Toggle GRAND/Miro': BEGIN
    if (prefix EQ "") then begin
	prefix = "g"
    endif else begin
	prefix = ""
    endelse
    dirname = strcompress (prefix+'out'+string(number),/remove_all)
    txtmsg, dirname
    END
  ENDCASE
END
;-----------------------------------------------------------------------
PRO PDMENUMASS_Event, Event
  CASE Event.Value OF 
  'Mass.Planet Mass(t)': BEGIN
    MassCurve
    END
  'Mass.Lost Mass(t)': BEGIN
    LostMassCurve
    END
  'Mass.Accretion Rate (t)': BEGIN
    AccretionCurve
    END
  ENDCASE
END
;-----------------------------------------------------------------------
PRO PDMENUdiscnb_Event, Event
common disk, dirname, number
common toggle, prefix
  CASE Event.Value OF 
  'Disk Number.Number 1': BEGIN
	number=1
    END
  'Disk Number.Number 2': BEGIN
	number=2
    END
  'Disk Number.Number 3': BEGIN
	number=3
    END
  'Disk Number.Number 4': BEGIN
	number=4
    END
  'Disk Number.Number 5': BEGIN
	number=5
    END
  'Disk Number.Number 6': BEGIN
	number=6
    END
  'Disk Number.Number 7': BEGIN
	number=7
    END
  'Disk Number.Number 8': BEGIN
	number=8
    END
  'Disk Number.Number 9': BEGIN
	number=9
    END
  'Disk Number.Number 0': BEGIN
	number=0
    END
  ENDCASE
  dirname = strcompress (prefix+'out'+string(number),/remove_all)
  txtmsg, dirname
END
;-----------------------------------------------------------------------
PRO PDMENUplanetnb_Event, Event
common disk, dirname, number
  common planetnb, pn

  CASE Event.Value OF 
  'Planet Number.Number 1': BEGIN
	pn=1
    END
  'Planet Number.Number 2': BEGIN
	pn=2
    END
  'Planet Number.Number 3': BEGIN
	pn=3
    END
  'Planet Number.Number 4': BEGIN
	pn=4
    END
  'Planet Number.Number 5': BEGIN
	pn=5
    END
  'Planet Number.Number 6': BEGIN
	pn=6
    END
  'Planet Number.Number 7': BEGIN
	pn=7
    END
  'Planet Number.Number 8': BEGIN
	pn=8
    END
  'Planet Number.Number 9': BEGIN
	pn=9
    END
  'Planet Number.Number 0': BEGIN
	pn=0
    END
  ENDCASE
  txtmsg, 'planet nb : '+string(pn)
END
;-----------------------------------------------------------------------
PRO PDMENUmassprefix_Event, Event
common disk, dirname, number
common prefix, massp, aspectratiop, diskmass
  CASE Event.Value OF 
  'Mass prefix.120': BEGIN
	massp=120
    END
  'Mass prefix.60': BEGIN
	massp=60
    END
  'Mass prefix.30': BEGIN
	massp=30
    END
  'Mass prefix.15': BEGIN
	massp=15
    END
  'Mass prefix.10': BEGIN
	massp=10
    END
  'Mass prefix.5': BEGIN
	massp=5
    END
  'Mass prefix.3': BEGIN
	massp=3
    END
  'Mass prefix.2': BEGIN
	massp=2
    END
  'Mass prefix.1': BEGIN
	massp=1
    END
  ENDCASE
  dirname=strcompress('out_'+string(aspectratiop)+'_'+string(massp)+'_'+string(diskmass),/remove_all)
  txtmsg,dirname
END
;-----------------------------------------------------------------------
PRO PDMENUdiskmassprefix_Event, Event
common disk, dirname, number
common prefix, massp, aspectratiop, diskmass
  CASE Event.Value OF 
  'Disk Mass.1': BEGIN
	diskmass=1
    END
  'Disk Mass.2': BEGIN
	diskmass=2
    END
  'Disk Mass.4': BEGIN
	diskmass=4
    END
  'Disk Mass.8': BEGIN
	diskmass=8
    END
  ENDCASE
  dirname=strcompress('out_'+string(aspectratiop)+'_'+string(massp)+'_'+string(diskmass),/remove_all)
  txtmsg,dirname
END
;-----------------------------------------------------------------------
PRO PDMENUaspectratioprefix_Event, Event
common disk, dirname, number
common prefix, massp, aspectratiop, diskmass
  CASE Event.Value OF 
  'Aspect ratio prefix.10': BEGIN
	aspectratiop=10
    END
  'Aspect ratio prefix.7': BEGIN
	aspectratiop=7
    END
  'Aspect ratio prefix.6': BEGIN
	aspectratiop=6
    END
  'Aspect ratio prefix.5': BEGIN
	aspectratiop=5
    END
  'Aspect ratio prefix.4': BEGIN
	aspectratiop=4
    END
  'Aspect ratio prefix.3': BEGIN
	aspectratiop=3
    END
  'Aspect ratio prefix.2': BEGIN
	aspectratiop=2
    END
  'Aspect ratio prefix.1': BEGIN
	aspectratiop=1
    END
  ENDCASE
  dirname=strcompress('out_'+string(aspectratiop)+'_'+string(massp)+'_'+string(diskmass),/remove_all)
  txtmsg,dirname
END
;-----------------------------------------------------------------------
PRO MAIN13_Event,Event
  common values,nonlyone,nfirst,nlast,firstrow
  common tcut,lastcut
  common plot1d, nname1
  common disk, dirname, number
  common slider, fine, coarse
  common PointsZoomFactor, points_zoom
  common smoothing, width
  common gas, gaschoice
  common polmode, polarchoice
  common pointstype, both
  common representation, mode
  WIDGET_CONTROL,Event.Id,GET_UVALUE=Ev
  widget_control, /HourGlass
  CASE Ev OF 
  'BGROUP2' : BEGIN
      firstrow = Event.Value
      end
  'BUTTONDIR': BEGIN
      DEP3 = WIDGET_BASE(ROW=1, MAP=1, TITLE = 'Disc #',UVALUE='DEP3')
  FieldVal1100 = [ number ]
  FIELD4 = CW_FIELD( DEP3,VALUE=FieldVal1100, COLUMN=1, STRING=1, TITLE='Please enter disc number', $
      UVALUE='FIELD4', /RETURN_EVENTS)
  WIDGET_CONTROL,DEP3,/REALIZE
  XMANAGER,"dep3",DEP3
      END
  'SLIDER6': BEGIN
      fine = event.value
      nonlyone = fine+coarse*100
      txtmsg, strcompress('Output  : '+string(nonlyone),/remove_all)
      END
  'SLIDER6a': BEGIN
      coarse = event.value
      nonlyone = fine+coarse*100
      txtmsg, strcompress('Output  : '+string(nonlyone),/remove_all)
      END
  'Share': BEGIN 
	PDMENUShare_Event, Event
      END
  'SLIDER9': BEGIN
      nfirst = event.value
      END
  'SLIDER10': BEGIN
      nlast = event.value
      END
  'SLIDER11': BEGIN
      width = event.value
      END
  'Zoom': BEGIN
      zoomdisplay
	END
  'momentum': PDMENU2_Event, Event
  'BUTTON21': BEGIN
      xloadct
      END
  'BUTTON22': BEGIN
      Widget_Control,event.top,/Destroy
      END
  'Hardcopy': BEGIN
      hardcopy
      END
  'Panels': BEGIN
      Panels
      END
  'GIF': BEGIN
      WriteAGIFFile
      END
  'Read': BEGIN
      readvalue
      END
  'Parameters': BEGIN
      parameters
      END
  'e(t)': BEGIN
      eccentricity
      END
  'a(t)': BEGIN
      majaxis
      END
  'da/dt': BEGIN
      dermajaxis
      END
  'GASZOOMP': BEGIN
      gaszoomplus
      END
  'GASZOOMM': BEGIN
      gaszoommoins
      END
  'DEFAULT': BEGIN
      gasdefault
      END
  'PLANET': BEGIN
      RocheLobe
      END
  'MASS': BEGIN
      Masscurve
      END
  'LOSTMASS': BEGIN
      LostMasscurve
      END
  'RADIUS': BEGIN
      Migration
      END
  'ACCRETION': BEGIN
      Accretioncurve
      END
  'GRID': BEGIN
      PlotGrid
      END
  'DETAIL': BEGIN
      PlotDetail
      END
  'Spectrum': BEGIN
      Spectrum
      END
  'Cut': BEGIN
      SpectrumCut
      END
  'Cursor'    : BEGIN
      curseur
      END
  'LIST26': BEGIN
      nname1 = event.index
      END
  'LIST27': BEGIN
      gaschoice = event.index
      END
  'BUTTONGAS' : BEGIN
      gasmenu, nonlyone
      END
  'CURVE' : BEGIN
      plot1
      END
  'BUTTONGASANIM' : BEGIN
      animgas
      END
  'PDMENUTRICKS': BEGIN
      PDMENUTRICKS_Event, Event
      END
  'menumass': PDMENUMASS_Event, Event
  'discnb' : pdmenudiscnb_event, Event
  'planetnb' : pdmenuplanetnb_event, Event
  'menumassprefix' : pdmenumassprefix_event, Event
  'menuaspectratioprefix' : pdmenuaspectratioprefix_event, Event
  'menudiskmassprefix' : pdmenudiskmassprefix_event, Event
  'polarmode' : BEGIN
      polarchoice = event.value
      END
  'BMODE' : BEGIN
      mode = event.value
      END
  ENDCASE
  status
END

PRO STATUS
  common disk, dirname, number
  common suffix,extension
  common share, numberscreens, indexscreen
  common wavenumber, wnb
  COMMON WTEXT5_Comm,WTEXT5_Id
  common textprefix_Comm, textprefix_Id
  common toggle, prefix
  common planetnb, pn
  common prefix, massp, aspectratiop, diskmass
  x = (indexscreen-1) mod numberscreens
  y = (indexscreen-1) / numberscreens
  if (numberscreens EQ 1) then begin
	x = 0
	y = 0
  endif
  if (prefix EQ "") then begin 
		machine_name = "Miro"
  endif else begin 
	machine_name = "GRAND" 
  endelse
  widget_control,WTEXT5_Id, set_value='Disk number '+strcompress(string(number))+$
				' on '+machine_name+$
				'  ext.: '+extension+$
				'  screen  : '+strcompress(string(x+1)+','+string(y+1)+$
'     /   m = '+string(wnb)+' Planet '+string(pn))
  widget_control,TEXTPREFIX_Id, set_value=dirname
end


;-------------------------------------------------------------------------
PRO wid,GROUP=Group
  common values,nonlyone,nfirst,nlast,firstrow
  common wavenumber, wnb
  common tcut,lastcut
  common imagesize,isize
  common plot1d, nname1
  common disk, dirname, number
  common modedetail, modet
  common spectraname, spname
  common share, numberscreens, indexscreen
  common slider, fine, coarse
  common smoothing, width
  COMMON WTEXT4_Comm,WTEXT4_Id
  COMMON WTEXT5_Comm,WTEXT5_Id
  common textprefix_Comm, textprefix_Id
  common gas, gaschoice
  common polmode, polarchoice
  common representation, mode
  common spectratype, sptype
  common suffix,extension
  common zoompolar,xcenter,ycenter,zoomradius
  common zoomspectra, zoomsp
  common animframe, frameanim
  common planetnb, pn
  common toggle, prefix
  common prefix, massp, aspectratiop, diskmass
  prefix = ""
  dirname = prefix+'out1'
  IF N_ELEMENTS(Group) EQ 0 THEN GROUP=0
  NumberOfSteps, nframe
  nframe=nframe*100
  maxcoarse = floor(nframe/100)
  maxfine   = 99
;  isize = 745
  isize = 480
;  WIDGET_CONTROL,DEFAULT_FONT='10x20'
  WIDGET_CONTROL,DEFAULT_FONT='6x9'
  junk   = { CW_PDMENU_S, flags:0, name:'' }
  MAIN13 = WIDGET_BASE(GROUP_LEADER=Group,RESOURCE_NAME='protoplanets',$
      COLUMN=1, MAP=1, TITLE='FARGO output visualizer', UVALUE='MAIN13')
  Btnsradio200 = [ 'image', 'surface']
  Btnsradio400 = [ 'r - theta', 'x - y' ]
  Btns100 = ['Erase', 'Hold']

  ListVal190 = ['Gas Density', 'Gas Rad. Veloc.', 'Rot. Veloc.']
  ListVal200 = ['Gas Density', 'log Gas Dens.', 'Horiz. Veloc.',$
                'Residual Veloc.', 'Hor. Vel./Planet',$
	        'Radial Veloc.', 'Rot. Veloc.', 'Label']

  MenuDescShare = [ $
      { CW_PDMENU_S,       3, 'Screen Share' }, $ ;        0
        { CW_PDMENU_S,       0, 'Only one screen' }, $ ;        1
        { CW_PDMENU_S,       0, '2x2 screen' }, $ ;        2
        { CW_PDMENU_S,       0, '3x3 screen' }, $ ;        3 
        { CW_PDMENU_S,       2, '4x4 screen' } $  ;      4
  ]

  MenuDescMomentum = [ $
      { CW_PDMENU_S,       3, 'Momemtum' }, $ ;        0
        { CW_PDMENU_S,       0, 'Momentum(t)' }, $ ;        1
        { CW_PDMENU_S,       1, 'Derivatives' }, $ ;        2
          { CW_PDMENU_S,       0, 'dMomentum/dt' }, $ ;        3
          { CW_PDMENU_S,       0, 'dMom/dt | Accretion' }, $ ;        4
          { CW_PDMENU_S,       2, 'difference' }, $ ;        5
      { CW_PDMENU_S,       1, 'Torques' }, $ ;        6
        { CW_PDMENU_S,       0, 'Inner Disk' }, $ ;        7
        { CW_PDMENU_S,       0, 'Outer Disk' }, $ ;        8
        { CW_PDMENU_S,       0, 'Inner Rings' }, $ ;        9
        { CW_PDMENU_S,       0, 'Outer Rings' }, $ ;        10
        { CW_PDMENU_S,       0, 'Indirect' }, $ ;        10
        { CW_PDMENU_S,       2, 'Total' }, $ ;       11
      { CW_PDMENU_S,       0, 'Disk Lost Mass' }, $  ;     12
      { CW_PDMENU_S,       2, 'Balance Residual' } $  ;     12
  ]

  MenuTricks = [ $
      { CW_PDMENU_S,       3, 'Tricks' }, $ ;        0
        { CW_PDMENU_S,       1, 'extensions' }, $ ;        1
          { CW_PDMENU_S,       0, '.dat' }, $ ;        2
          { CW_PDMENU_S,       0, '.old1' }, $ ;        3
          { CW_PDMENU_S,       0, '.old2' }, $ ;        4
          { CW_PDMENU_S,       0, '.old3' }, $ ;        5
          { CW_PDMENU_S,       0, '.old4' }, $ ;        6
          { CW_PDMENU_S,       0, '.old5' }, $ ;        7
          { CW_PDMENU_S,       0, '.old6' }, $ ;        8
          { CW_PDMENU_S,       0, '.old7' }, $ ;        9
          { CW_PDMENU_S,       0, '.old8' }, $ ;       10
          { CW_PDMENU_S,       0, '.old9' }, $ ;       11
          { CW_PDMENU_S,       2, '.old0' }, $ ;       12
        { CW_PDMENU_S,       1, 'Wave Number' }, $ ;        1
          { CW_PDMENU_S,       0, 'm = 0' }, $ ;        2
          { CW_PDMENU_S,       0, 'm = 1' }, $ ;        3
          { CW_PDMENU_S,       0, 'm = 2' }, $ ;        4
          { CW_PDMENU_S,       0, 'm = 3' }, $ ;        5
          { CW_PDMENU_S,       0, 'm = 4' }, $ ;        6
          { CW_PDMENU_S,       0, 'm = 5' }, $ ;        7
          { CW_PDMENU_S,       0, 'm = 6' }, $ ;        8
          { CW_PDMENU_S,       0, 'm = 7' }, $ ;        9
          { CW_PDMENU_S,       0, 'm = 8' }, $ ;       10
          { CW_PDMENU_S,       2, 'm = 9' }, $ ;       11
        { CW_PDMENU_S,       1, 'Spectra Zoom' }, $ ;        1
          { CW_PDMENU_S,       0, '1x' }, $ ;        2
          { CW_PDMENU_S,       0, '2x' }, $ ;        3
          { CW_PDMENU_S,       0, '4x' }, $ ;        4
          { CW_PDMENU_S,       0, '8x' }, $ ;        5
          { CW_PDMENU_S,       2, '16x' }, $ ;        6
      { CW_PDMENU_S,       1, 'Velocity field' }, $ ;       13
        { CW_PDMENU_S,       0, 'V-VKep' }, $ ;       14
        { CW_PDMENU_S,       2, 'Omega-Omega_planet' }, $  ;     15
      { CW_PDMENU_S,       1, 'Animations' }, $ ;       13
        { CW_PDMENU_S,       0, 'Rotating frame' }, $ ;       14
        { CW_PDMENU_S,       2, 'Fix frame' }, $  ;     15
      { CW_PDMENU_S,       1, 'Spectra Aspect' }, $ ;       13
        { CW_PDMENU_S,       0, 'Linear' }, $ ;       14
        { CW_PDMENU_S,       0, 'Log' }, $ ;       14
        { CW_PDMENU_S,       2, 'Sqrt' }, $  ;     15
      { CW_PDMENU_S,       1, 'Spectra Type' }, $ ;       13
        { CW_PDMENU_S,       0, 'Radius/Frequency' }, $ ;       14
        { CW_PDMENU_S,       0, 'Radius/Time (cos)' }, $ ;       14
        { CW_PDMENU_S,       2, 'Radius/Time (sin)' }, $  ;     15
      { CW_PDMENU_S,       0, 'Resonances' }, $ ;       13 
      { CW_PDMENU_S,       0, 'Toggle GRAND/Miro' }, $ ;       13 
      { CW_PDMENU_S,       2, 'Close Files' } $ ;       13 
  ]
  MenuMass = [ $
      { CW_PDMENU_S,       3, 'Mass' }, $ ;        0
        { CW_PDMENU_S,       0, 'Planet Mass(t)' }, $ ;        1
        { CW_PDMENU_S,       0, 'Lost Mass(t)' }, $ ;        2
        { CW_PDMENU_S,       2, 'Accretion Rate (t)' } $  ;      3
  ]

  MenuDiscNumber = [ $
      { CW_PDMENU_S,       3, 'Disk Number' }, $ ;        0
        { CW_PDMENU_S,       0, 'Number 1' }, $ ;        1
        { CW_PDMENU_S,       0, 'Number 2' }, $ ;        2
        { CW_PDMENU_S,       0, 'Number 3' }, $ ;        3
        { CW_PDMENU_S,       0, 'Number 4' }, $ ;        4
        { CW_PDMENU_S,       0, 'Number 5' }, $ ;        5
        { CW_PDMENU_S,       0, 'Number 6' }, $ ;        6
        { CW_PDMENU_S,       0, 'Number 7' }, $ ;        7
        { CW_PDMENU_S,       0, 'Number 8' }, $ ;        8
        { CW_PDMENU_S,       0, 'Number 9' }, $ ;        9
        { CW_PDMENU_S,       2, 'Number 0' } $  ;     10

  ]

  MenuPlanetNumber = [ $
      { CW_PDMENU_S,       3, 'Planet Number' }, $ ;        0
        { CW_PDMENU_S,       0, 'Number 0' }, $ ;        9
        { CW_PDMENU_S,       0, 'Number 1' }, $ ;        1
        { CW_PDMENU_S,       0, 'Number 2' }, $ ;        2
        { CW_PDMENU_S,       0, 'Number 3' }, $ ;        3
        { CW_PDMENU_S,       0, 'Number 4' }, $ ;        4
        { CW_PDMENU_S,       0, 'Number 5' }, $ ;        5
        { CW_PDMENU_S,       0, 'Number 6' }, $ ;        6
        { CW_PDMENU_S,       0, 'Number 7' }, $ ;        7
        { CW_PDMENU_S,       0, 'Number 8' }, $ ;        8
        { CW_PDMENU_S,       2, 'Number 9' } $  ;     10

  ]

  MenuMassPrefix = [ $
      { CW_PDMENU_S,       3, 'Mass prefix' }, $ ;        0
        { CW_PDMENU_S,       0, '120' }, $ ;        9
        { CW_PDMENU_S,       0, '60' }, $ ;        9
        { CW_PDMENU_S,       0, '30' }, $ ;        1
        { CW_PDMENU_S,       0, '15' }, $ ;        2
        { CW_PDMENU_S,       0, '10' }, $ ;        3
        { CW_PDMENU_S,       0, '5' }, $ ;        4
        { CW_PDMENU_S,       0, '3' }, $ ;        5
        { CW_PDMENU_S,       0, '2' }, $ ;        6
        { CW_PDMENU_S,       2, '1' } $ ;        7
  ]

  MenuAspectRatioPrefix = [ $
      { CW_PDMENU_S,       3, 'Aspect ratio prefix' }, $ ;        0
        { CW_PDMENU_S,       0, '10' }, $ ;        9
        { CW_PDMENU_S,       0, '7' }, $ ;        1
        { CW_PDMENU_S,       0, '6' }, $ ;        2
        { CW_PDMENU_S,       0, '5' }, $ ;        3
        { CW_PDMENU_S,       0, '4' }, $ ;        4
        { CW_PDMENU_S,       0, '3' }, $ ;        5
        { CW_PDMENU_S,       0, '2' }, $ ;        6
        { CW_PDMENU_S,       2, '1' } $ ;        7
  ]

  MenuDiskMassPrefix = [ $
      { CW_PDMENU_S,       3, 'Disk Mass' }, $ ;        0
        { CW_PDMENU_S,       0, '8' }, $ ;        9
        { CW_PDMENU_S,       0, '4' }, $ ;        1
        { CW_PDMENU_S,       0, '2' }, $ ;        2
        { CW_PDMENU_S,       2, '1' } $ ;        7
  ]

  MAN13        = WIDGET_BASE(MAIN13,   ROW=1, MAP=1, UVALUE='B0')
  baseleft     = WIDGET_BASE(MAN13,    COLUMN=1, MAP=1, UVALUE='baseleft')
  base2ndleft  = WIDGET_BASE(MAN13,    COLUMN=1, MAP=1, UVALUE='baseleft')
  basewins     = WIDGET_BASE(MAN13,    COLUMN=1, MAP=1, UVALUE='basewins')
  baseprefix  = WIDGET_BASE(basewins, ROW=1, MAP=1, UVALUE='baseprefix', frame=3)
  basetopleft  = WIDGET_BASE(baseleft, COLUMN=1, MAP=1, UVALUE='basetopleft')
  basetimestep = WIDGET_BASE(baseleft, COLUMN=1, MAP=1, UVALUE='basetimestep', frame=3)
  baseaspect   = WIDGET_BASE(baseleft, COLUMN=1, MAP=1, UVALUE='baseaspect', frame=3)
  baseanim     = WIDGET_BASE(base2ndleft, COLUMN=1, MAP=1, UVALUE='baseanim', frame=3)
  basecurve    = WIDGET_BASE(base2ndleft, COLUMN=1, MAP=1, UVALUE='basecurve', frame=3)
  baseplot     = WIDGET_BASE(base2ndleft, COLUMN=1, MAP=1, UVALUE='baseplot', frame=3)
  basetasks    = WIDGET_BASE(MAIN13,    ROW=1, MAP=1, UVALUE='basetasks',frame=3)

  PDMENUNUMBER     = CW_PDMENU( basetopleft, MenuDiscNumber, /RETURN_FULL_NAME, UVALUE='discnb')
  ;BUTTONDIR        = WIDGET_BUTTON( basetopleft, UVALUE='BUTTONDIR', VALUE='Disc Number')
  BUTTONParameters = WIDGET_BUTTON (basetopleft, UVALUE='Parameters', VALUE='Parameters')
  WTEXT4           = WIDGET_TEXT( basetopleft, UVALUE='WTEXT4', YSIZE=1)
  ;LABEL11          = WIDGET_LABEL( basetimestep,UVALUE='',VALUE='time step selection')
  SLIDER6          = WIDGET_SLIDER( basetimestep, MAXIMUM=maxfine, MINIMUM=0,$
                     TITLE='Fine', UVALUE='SLIDER6', VALUE=0)
  SLIDER6a         = WIDGET_SLIDER( basetimestep, MAXIMUM=maxcoarse, MINIMUM=0,$
                     TITLE='x 100', UVALUE='SLIDER6a', VALUE=0)
  SLIDER11         = WIDGET_SLIDER ( basetimestep, MAXIMUM=201, MINIMUM=1,$
                     TITLE='Smoothing mask width', UVALUE='SLIDER11', VALUE=1)
  BGROUP2          = CW_BGROUP(baseaspect, Btns100, COLUMN=1, EXCLUSIVE=1, SET_VALUE=0, $
                     LABEL_TOP='Screen', UVALUE='BGROUP2')
  Bpolar           = CW_BGROUP( baseaspect, Btnsradio400, $
                     COLUMN=1, $
                     EXCLUSIVE=1, $
                     LABEL_TOP='Gas Polar Mode', $
                     UVALUE='polarmode', SET_VALUE=0) 
  BMODE            = CW_BGROUP( baseaspect, Btnsradio200, $
                     COLUMN=1, $
                     EXCLUSIVE=1, $
                     LABEL_TOP='Representation', $
                     UVALUE='BMODE', SET_VALUE=0) 
  LABEL11          = WIDGET_LABEL  ( baseanim, UVALUE='', VALUE='Animation limits')
  SLIDER9          = WIDGET_SLIDER ( baseanim, MAXIMUM=nframe, MINIMUM=0,$
  TITLE            = 'first', UVALUE='SLIDER9', VALUE=0)
  SLIDER10         = WIDGET_SLIDER ( baseanim, MAXIMUM=nframe, MINIMUM=1,$
                     TITLE='last', UVALUE='SLIDER10', VALUE=nframe)
  BUTTONGASANIM    = WIDGET_BUTTON( baseanim, UVALUE='BUTTONGASANIM', VALUE='Start')
  LABEL12          = WIDGET_LABEL ( basecurve, UVALUE='', VALUE='Curve Selection')
  LIST26           = WIDGET_LIST  ( basecurve,YSIZE=N_ELEMENTS(ListVal190),$
                     VALUE=ListVal190, UVALUE='LIST26')
  LABEL13          = WIDGET_LABEL ( baseplot, UVALUE='', VALUE='Plot Selection')
  LIST27           = WIDGET_LIST  ( baseplot,YSIZE=N_ELEMENTS(ListVal200),$
                     VALUE=ListVal200, UVALUE='LIST27')

  ButtonRead       = WIDGET_BUTTON(BASE2NDLEFT, UVALUE='Read', VALUE='Read Value') 
  ButtonZoom       = WIDGET_BUTTON(BASE2NDLEFT, UVALUE='Zoom', VALUE='Zoom') 
  BUTTONCursor     = WIDGET_BUTTON (BASE2NDLEFT, UVALUE='Cursor', VALUE='Cursor')
  BUTTON21         = WIDGET_BUTTON( BASE2NDLEFT, UVALUE='BUTTON21', VALUE='Palette')
  BUTTONHardCopy   = WIDGET_BUTTON (BASE2NDLEFT, UVALUE='Hardcopy', VALUE='Hardcopy')
  BUTTONPanels     = WIDGET_BUTTON (BASE2NDLEFT, UVALUE='Panels', VALUE='Panels')
  BUTTONGif        = WIDGET_BUTTON (BASE2NDLEFT, UVALUE='GIF', VALUE='GIF')
  PDMENUPlanetNb   = CW_PDMENU(base2ndleft, MenuPlanetNumber, /RETURN_FULL_NAME, UVALUE='planetnb')
  PDMENUShare      = CW_PDMENU( BASELEFT, MenuDescShare, /RETURN_FULL_NAME, UVALUE='Share')
  BUTTON22         = WIDGET_BUTTON( BASELEFT, UVALUE='BUTTON22', VALUE='quit')

;  PDMENU2          = CW_PDMENU(basetasks, MenuDescMomentum, /RETURN_FULL_NAME, UVALUE='momentum')
;  PDMENUMASS       = CW_PDMENU(basetasks, MenuMass, /RETURN_FULL_NAME, UVALUE='menumass')
  BUTTONMajAxis     = WIDGET_BUTTON (basetasks, UVALUE='a(t)', VALUE='a(t)')
  BUTTONdermajaxis     = WIDGET_BUTTON (basetasks, UVALUE='da/dt', VALUE='da/dt')
  BUTTONEcc     = WIDGET_BUTTON (basetasks, UVALUE='e(t)', VALUE='e(t)')
  PDMENUMASSPREFIX       = CW_PDMENU(baseprefix, MenuMassPrefix, /RETURN_FULL_NAME, UVALUE='menumassprefix')
  PDMENUASPECTRATIOPREFIX       = CW_PDMENU(baseprefix, MenuAspectRatioPrefix, /RETURN_FULL_NAME, UVALUE='menuaspectratioprefix')
  PDMENUDISKMASSPREFIX       = CW_PDMENU(baseprefix, MenuDiskMassPrefix, /RETURN_FULL_NAME, UVALUE='menudiskmassprefix')
  BUTTONCURVE      = WIDGET_BUTTON( basetasks, UVALUE='CURVE', VALUE='Curve Plot')
  BUTTONGAS        = WIDGET_BUTTON( basetasks, UVALUE='BUTTONGAS', VALUE='Gas Plot')
  BUTTONGasZoomP   = WIDGET_BUTTON (basetasks, UVALUE='GASZOOMP', VALUE='Zoom +')
  BUTTONGasZoomM   = WIDGET_BUTTON (basetasks, UVALUE='GASZOOMM', VALUE='Zoom -')
  BUTTONGasDefault = WIDGET_BUTTON (basetasks, UVALUE='DEFAULT', VALUE='Default')
  BUTTONPlanet     = WIDGET_BUTTON (basetasks, UVALUE='PLANET', VALUE='Planet')
  BUTTONRadius     = WIDGET_BUTTON (basetasks, UVALUE='RADIUS', VALUE='Radius(t)')
  ;BUTTONAccretion  = WIDGET_BUTTON (basetasks, UVALUE='ACCRETION', VALUE='Acc. Rate(t)')
  BUTTONGrid       = WIDGET_BUTTON (basetasks, UVALUE='GRID', VALUE='Grid')
  BUTTONDetail     = WIDGET_BUTTON (basetasks, UVALUE='DETAIL', VALUE='Detail')
  PDMENUTRICKS     = CW_PDMENU(basetasks, MenuTricks, /RETURN_FULL_NAME, UVALUE='PDMENUTRICKS')
  BUTTONSpectre     = WIDGET_BUTTON (basetasks, UVALUE='Spectrum', VALUE='Spectrum')
  BUTTONCut     = WIDGET_BUTTON (basetasks, UVALUE='Cut', VALUE='Cut')


  DRAW27           = WIDGET_DRAW( basewins, RETAIN=2, UVALUE='DRAW27',XSIZE=isize,YSIZE=isize)
  WTEXT5           = WIDGET_TEXT( basewins, UVALUE='WTEXT5', YSIZE=1)
  TEXTPREFIX           = WIDGET_TEXT( baseprefix, UVALUE='TESTPREFIX', YSIZE=1)

  widget_control,LIST26,SET_LIST_SELECT=0
  widget_control,LIST27,SET_LIST_SELECT=0
  WIDGET_CONTROL,MAIN13,/REALIZE
  WTEXT4_Id = WTEXT4
  WTEXT5_Id = WTEXT5
  textprefix_Id = TEXTPREFIX
  COMMON DRAW27_Comm,DRAW27_Id
  WIDGET_CONTROL,DRAW27,GET_VALUE=DRAW27_Id
  loadct,27	; Eos B
  xcenter = 0.0
  ycenter = 0.0
  zoomradius = 3.0
  nonlyone = 0
  nfirst = 0
  nlast = nframe
  nname1 = 0
  firstrow = 0
  number = '1' 
  mode = 0
  numberscreens = 1
  indexscreen = 1
  coarse = 0
  fine = 0
  gaschoice = 0
  polarchoice = 0
  width = 1
  extension='.dat'
  modet = 2
  wnb = 0
  sptype = 0
  zoomsp = 1.0
  frameanim = 0
  massp = 15
  aspectratiop = 4
  diskmass = 1
  spname = "spec"
  pn = 0
  status
  XMANAGER,'MAIN13',MAIN13
END

wid
END


