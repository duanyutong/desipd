pro pd_fpc,infile,imphot,indir=indir1,x0=x01,y0=y01,noprint=noprint1

; Input:
; 	infile (string) - input fits image name
;
; Output:
;	imphot (struct) - structure containing photometry params
; 
; Keywords:
; 	indir (string) - full path to data directory
;	x0,y0 (float) - center coordinate of bottom left fiber
;
; - Arjun (2016-08-05)

if n_params() lt 1 then begin
	print,'SYNTAX: pd_fpc,infile,imphot,indir=indir,x0=x0,y0=y0'
	print,'   where '
	print,'  infile = name of input FPC fits file'
	print,'  indir  = name of directory with input FITS file (default = ./)'
	print,'  imphot = name of output structure containing photometry'
	print,'  x0,y0 = pixel location of bottom left fiber (default = 1211,1067)'
	print,'  /print = print out the photometry'
	retall
endif

resolve_all, /continue_on_error, /quiet ; Avoid printing these while running

if keyword_set(x01) then x0 = x01 $
	else x0 = 1211.
if keyword_set(y01) then y0 = y01 $
	else y0 = 1067.
if (keyword_set(noprint1)) then noprint = noprint1 else noprint = 0
if (keyword_set(indir1)) then indir = indir1 $
    else indir = getenv('PDFPC_DATA')
if (NOT keyword_set(indir)) then $
    print, 'INDIR keyword and $PDFPC_DATA not set; default to current'

; Define the initial centroids of fiber images

fn0=['F1','F2','F3','F4']    ; fiber names
xc0=x0+[0.,1.,674.,681.] ; x-coord of fiber on FPC
yc0=y0+[0.,323.,0.,321.] ; y-coord of fiber on FPC
nfib = n_elements(xc0) ; number of fibers

; Define apertures

apr = [20.,25.,30.,35.,40.]  ; aperture radii in pixels
skyrad = [41.,60.]
nap = n_elements(apr)

; Define the output structure

imphot=create_struct( $
	'FILENAME','', $
	'EXPNUM',0L, $
	'DATEOBS','', $
	'MJDOBS',0.D0, $
	'FIBNAME','', $
	'APERDIA',fltarr(nap), $
	'FLUX',fltarr(nap), $
	'FLUXERR',fltarr(nap), $
	'APERPOSX',fltarr(nap), $
	'APERPOSY',fltarr(nap), $
	'BKGCTS',0., $
	'BKGERR',0.)
imphotj = imphot

; Read in the image

img=readfits(infile,hdr)

; define some image parameters

phpadu = 1.0 ; gain
badpix = [0.,32767]

display,img,/aspect,top=255
plotsym,0,3,color=c24(3),thick=3
oplot,xc0,yc0,psym=8

; Measure the photometry

for j=0,nfib-1 do begin
	aper, img, xc0[j], yc0[j], flux, errap, skyv, skyverr, phpadu, $
		apr, skyrad, badpix, /silent, /nan, /flux
	; allmags = allmags ; - 25.0 + zpt0 + 2.5*alog10(exptime)

	imphotj.filename=infile
	imphotj.dateobs = sxpar(hdr,'DATE-OBS')
	imphotj.mjdobs = sxpar(hdr,'MJD-OBS')
	imphotj.expnum = sxpar(hdr,'EXPNUM')
	imphotj.fibname=fn0[j]
	imphotj.aperdia=apr
	imphotj.flux = flux
	imphotj.fluxerr = errap
	imphotj.aperposx = xc0[j]
	imphotj.aperposy = yc0[j]
	imphotj.bkgcts = skyv
	imphotj.bkgerr = skyverr

	if j eq 0 then $
		imphot = imphotj $
	else $
		imphot = struct_append(imphot,imphotj)

endfor

if (NOT keyword_set(noprint)) then begin
	fmt = '(3x,a3,'+string(nap*2,format='(i2)')+'(x,f10.1))'
	fmt1 = '('+string(nap,format='(i1)')+'(f4.1,","))'
	print,'--------------------------'
	print,'FILE: '+infile
	print,'FIBER  FLUX ('+string(apr,format=fmt1)+')   ERRFLUX'
	for j=0,nfib-1 do begin
		print,imphot[j].fibname,imphot[j].flux,imphot[j].fluxerr,format=fmt
	endfor
endif

end
