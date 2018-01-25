;+
; NAME:
;   pd_fpc_continuous
;
; PURPOSE:
;   Continuously run pd_fpc
;
; CALLING SEQUENCE:
;   pd_fpc_continuous, [ expr1, topdir=, wtime=, $
;    /noloop, /verbose, /debug, _EXTRA= ]
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   expr       - Match expression for files; default to 'PROTODESI_FPC_????????.fits'
;   wtime      - Wait time before looking for new files; default to 10 sec
;   noloop     - Set to only read the files in the specified directory once,
;                then exit without waiting for new files to be written
;   _EXTRA     - Extra keywords to pass to PD_FPC, such as INDIR=, etc.
;
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   This code can be run as an infinite loop while observing at the telescope
;   to continuously update report the latest image summary information
;   be OBTYPE='object' and EXPTIME >= 30 sec.
;
;   It is simply a wrapper around PD_FPC to search for new files.
;
; EXAMPLES:
;
; BUGS:
;
; DATA FILES:
;   $TOPDIR/PROTODESI_FPC_????????.fits.fz
;
; REVISION HISTORY:
;   17-Aug-2016: A. Dey (Based on decstat_continuous, which was written by B. Nord)
;-
;------------------------------------------------------------------------------
; Get the sorted list of all raw data files; the last files are assumed
; to be the most recent.
function fpc_search, expr, topdir=topdir, nfile=nfile

   files = file_search(filepath(expr, root_dir=topdir), count=nfile)
   if (nfile GT 0) then files = files[sort(files)]
   return, files
end

;------------------------------------------------------------------------------
; run an open loop and call decstat on a new file whenever it appears
pro pd_fpc_continuous, expr1, topdir=topdir1, wtime=wtime1, $
 noloop=noloop, verbose=verbose1, debug=debug, _EXTRA=ExtraKeywords

   resolve_all, /continue_on_error, /quiet ; Avoid printing these while running

   ; set keywords
   if (keyword_set(expr1)) then expr = expr1 $
      else expr = 'PROTODESI_FPC_????????.fits'
   if (keyword_set(topdir1)) then topdir = topdir1 $
      else topdir = getenv('PDFPC_DATA')
   if (keyword_set(wtime1)) then wtime = (wtime1 > 1) $
      else wtime = 10
   if (keyword_set(verbose1)) then verbose = verbose1  $
      else verbose = 1
   if (NOT keyword_set(topdir)) then begin
      print, 'TOPDIR keyword or $PDFPC_DATA env must be set!'
      return
   endif

   nfile = 0L
   qdone = 0B
   file_temp = ''
   while (qdone EQ 0) do begin

      files = fpc_search(expr, topdir=topdir, nfile=nfile_new)

      for j=0., wtime, 1. do begin
         wait, 1
         if verbose then $
          print, 'PD_FPC: Last check '+systime()+string(13b),format='($,a,a1)'
      endfor

;      file_current = strmid(files[nfile_new-1],30,23)
      file_current = fileandpath(files[nfile_new-1]) ; strip directory name
;     file_current = repstr(file_current, '.fz', '') ; stripe '.fz' extension
      if (file_current NE file_temp) then begin
         print
         pd_fpc, file_current, _EXTRA=ExtraKeywords
      endif

      file_temp = file_current
      qdone = keyword_set(noloop)

   endwhile

   return
end

