PRO Nlparv5, patfile, outfile, lambda = lambda, dthresh = dthresh,searchWindow = searchWindow, $
  up1 = up1, up2=up2, nColScan = nColScan, nRowScan = nRowScan, $
  BackGroundSub=BackGroundSub, rescale = rescale,  $
  mask = mask, w0=w0, sigma_in=sigma_in, sigma_out=sigma_out,$
  RAMWarningoff = RAMWarningoff, blocksizein=blocksizein, threads = threads


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  This software was developed by employees of the US Naval Research Laboratory
; (NRL), an agency of the Federal Government. Pursuant to title 17 section 105 of
; the United States Code, works of NRL employees are not not subject to copyright
; protection, and this software is in the public domain. FiPy is an experimental
; system. NRL assumes no responsibility whatsoever for its use by other parties,
; and makes no guarantees, expressed or implied, about its quality, reliability,
; or any other characteristic. We would appreciate acknowledgment if the software
; is used.

; To the extent that NRL may hold copyright in countries other than the United
; States, you are hereby granted the non-exclusive irrevocable and unconditional
; right to print, publish, prepare derivative works and distribute this software,
; in any medium, or authorize others to do so on your behalf, on a royalty-free
; basis throughout the world.

; You may improve, modify, and create derivative works of the software or any
; portion of the software, and you may copy and distribute such modifications or
; works. Modified works should carry a notice stating that you changed the
; software and should note the date and nature of any such change. Please
; explicitly acknowledge the US Naval Research Laboratory as the original source.

; This software can be redistributed and/or modified freely provided that any
; derivative works bear some notice that they are derived from it, and any
; modified versions bear some notice that they have been modified.

 
; Author: David Rowenhorst; The US Naval Research Laboratory
; Date: 5 Dec 2018   
;
; Description: This is a wrapper code for the NLPAR algorithm
; This code will read in a up1/2 (version 3).
;
; Arguments (in IDL, arguements need to input in order):
;     patfile (input/required): A string giving the file location of an
;         EDAX up1 or up2 version 3 file that contains the patterns to be processed.
;
;     outfile (output/required): A string giving the file location of the output file.
;          No check will occur if a file is being overwritten.
;
; Keywords (in IDL keywords do not need to be in order):
;     lambda (input/default=0.6): a floating point number giving the magnitude of the
;          lambda smoothing factor.
;     dthresh (input/default=0.0): the normalized distance threshold for which two patterns are considered
;          identical to the pattern of interest.
;     searchWindow (input/default=3): value of the search window radius.  A value of 3 will
;          produce a search box of (2*(3)+1) x (2*(3)+1) == 7 x 7 search box.
;     up1/up1 (flag/required): Set one of these keywords to 1 to indicate that the
;          input file is either a up1 or up2 file.
;     nColScan (input/optional): A long integer of the number of columns in the EBSD scan. This is only used if a 
;           version 1 up1/2 is input.  Otherwise it is not required and ignored. 
;     nRowScan (input/optional): A long integer of the number of rows in the EBSD scan. This is only used if a 
;           version 1 up1/2 is input.  Otherwise it is not required and ignored. 
;     BackGroundSub (flag/default=0): optional flag that if set will do a background subtraction of
;          the patterns based on the "mean" pattern for the block.
;     rescale (flag/default=0): If set, the patterns will be rescaled to fit within the numberic range of the
;           output.
;     mask (input/optinal): an optinal input that is set to an array that is the size
;           of one pattern.  Values that are GT 0 will be considered in the NLPAR algorithm.
;           If not set to an valid array a default mask is created that is the standard cicular mask
;     w0 (output/optional): an optional output array that will be the size of the scan that contains
;           the weighting factor applied to the original pattern -- useful for detecting how much smoothing
;           was applied to pattern
;     sigma_in (input/optional): optinal input that is an array that is the same size of the scan that contains a
;           estimate of the standard deviation of the pattern noise at each point. If not provided a value will be
;           calculated using nearest neighbor values.
;     sigma_out (output/optional): output array that contains the results of the nearest neighbor sigma estimate
;           if sigma_in is not provided.
;     RAMWarningoff (flag, optional): This will turn the RAM warning off.  Not advised if you are running on a laptop.
;           If not set, the program will do a rough analysis of the parameters read in from the up1/2 files and will
;           ask the user if they want to continue if the RAM requirement is over 16GB.
;     blocksizein (input, default=150): this is the number of scan lines that will be processed at a time. In total, the
;           blocksize + 2*searchWindow will be read in -- potentially a large amount of RAM.
;     threads (input, default = 0 (all cpu cores)):  the number of threads that are availble for processing.  If
;           set to 0, all the cpu cores will be used.
;
;
;  Example input:
;
; nlparv5, '~/EBSDdata/scan1.up2', '~/EBSDdata/scan1nlpar.up1',  lambda = 1.10, dthresh = 0.0,searchWindow = 5, $
;  /up2, BackGroundSub=0,  mask = patmask, w0=weight0, sigma_out=sigma_NN,$
;  blocksizein=50, threads = 0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



  IF N_ELEMENTS(threads) EQ 0 THEN threads = 0
  IF threads LE 0 THEN threads = 0
  threads = LONG(threads)

  bitdepth = 0ll
  routine_info = routine_info('nlparv5', /source)
  path = FILE_DIRNAME(routine_info.path)
  library = "patternProcessingC.dylib"
  IF !VERSION.os_name EQ "linux" THEN library = "patternProcessingC_Linux64.so"
  library = path + PATH_SEP() + library
  IF KEYWORD_SET(up1) THEN bitdepth = 8ll
  IF KEYWORD_SET(up2) THEN bitdepth = 16ll
  IF bitdepth EQ 0 THEN BEGIN
    PRINT, 'Error: Keyword "UP1" or "UP2" must be set'
    RETURN
  ENDIF


  IF N_ELEMENTS(searchWindow) EQ 0 THEN searchWindow = 3
  IF N_ELEMENTS(searchWindow) EQ 1 THEN $
    searchWindow = REBIN(REFORM([searchWindow]), 2, /sam)


  IF N_ELEMENTS(lambda) EQ 0 THEN lambda = 0.6
  IF N_ELEMENTS(dthresh) EQ 0 THEN dthresh = 0.0


  bitdepthmax = 2ull^bitdepth

  typecode = (bitdepth EQ 8) + 12 * (bitdepth NE 8)


  IF KEYWORD_SET(up2) THEN onepat = Read_up2(patfile, read_interval=[0,0]) ELSE onepat = Read_up1(patfile, read_interval=[0,0])
  
  IF (*onepat.header).version GT 1 THEN BEGIN
    image_dim = LONG([(*onepat.header).cols, (*onepat.header).rows])
  ENDIF ELSE BEGIN
    IF ( (N_ELEMENTS(nColScan) GT 0) AND (N_ELEMENTS(nRowScan) GT 0) ) THEN BEGIN
      image_dim = LONG([nColScan[0], nRowScan[0]])
    ENDIF ELSE BEGIN
      Print, "Error: A version 1 up1/2 file has been read.  nColScan and nRowScan must be input."
      RETURN
    ENDELSE  
  ENDELSE
  w0 = FLTARR(image_dim, /nozero)
  sigma_out = FLTARR(image_dim, /nozero)
  nPoints = Product(image_dim, /int)

  FILE_COPY, patfile, outfile, /overwrite

  pat_dim = LONG([(*onepat.header).pattern_width, (*onepat.header).pattern_height])
  patstride = Product(pat_dim, /int);*bitdepth/8ll

  start_write0 = (*onepat.header).file_pos

  GET_LUN, lunout

  OPENU, lunout, outfile
  POINT_LUN, lunout, start_write0

  searchWindow <= image_dim
  IF N_ELEMENTS(blocksizein) EQ 0 THEN blocksize = 150 ELSE blocksize = blocksizein[0]
  blocksize += searchWindow[1]*2 ;number of rows to analyze at a time.
  blocksize <= image_dim[1]

  blockSizeByte = blockSize*image_dim[0]*pat_dim[0]*pat_dim[1]*(bitDepth/8)
  PRINT, "Approximate RAM usage: " + STRTRIM( 8.5*blocksizeByte/1.e9, 2) + "GB"

  IF ((8.5*blockSizeByte) GT 15e9) AND (~KEYWORD_SET(RAMwarningOff)) THEN BEGIN
    PRINT, "You are about to process a block of data that is going to be over 16GB in size."
    PRINT, "  Do you really want to continue? [Y/N]"
    Answer = ''
    READ, Answer, Prompt = 'Yes/[N]o:'
    answer = STRTRIM(answer, 2)
    IF answer EQ '' THEN answer = 'N'

    IF ~STRCMP(STRMID(answer, 0, 1) , 'Y', /fold_case) THEN BEGIN
      PRINT, "Returning.  Try a smaller block size."
      RETURN
    ENDIF

  ENDIF

  ;IF KEYWORD_SET(automask) THEN BEGIN
  IF N_ELEMENTS(mask) EQ 0 THEN BEGIN
    mask = SHIFT(Dist(pat_dim[0], pat_dim[1]), pat_dim[0]/2, pat_dim[1]/2) LT (MIN(pat_dim)/2)
  ENDIF

  IF N_ELEMENTS(mask) EQ 0 THEN maskindex = LINDGEN(patStride) ELSE maskindex = WHERE(mask GT 0, complement = notmaskindex)
  IF N_ELEMENTS(maskindex) LE 1 THEN maskindex = LINDGEN(patStride)
  IF N_ELEMENTS(notmaskindex) LE 1 THEN notmaskindex = -1
  nMask = LONG(N_ELEMENTS(maskindex))




  ;maskindex =  maskindex[0:*:10]

  ;newPatBlock = make_array(dimension = [pat_dim, image_dim[0], blocksize], type = typecode )
  patBuffer = MAKE_ARRAY(dimension = [pat_dim, image_dim[0], blocksize], type = typecode, /nozero )
  onePat = MAKE_ARRAY(dimension = [pat_dim], type = typecode )
  background = MAKE_ARRAY(dimension = [pat_dim], type = 4 )
  patBlock = MAKE_ARRAY(dimension = [nMask, image_dim[0], blocksize], type = 4, /nozero )
  nlparBlock = MAKE_ARRAY(dimension = [nMask, image_dim[0], blocksize], type = 4, /nozero )
  w0block = MAKE_ARRAY(dimension = [image_dim[0], blocksize], type = 4, /nozero )


  ; make sure our types are correct for the c-program
  cImage_dim = LONG([image_dim[0], blockSize])
  windowSize = LONG(searchWindow)
  pat_dim = LONG(pat_dim)
  nMask = LONG(N_ELEMENTS(maskindex))
  maskindex = LONG(maskindex)
  xstart = LONG(0)
  xend =  LONG(image_dim[0]-1)
  ystart =  LONG(0)
  IF blocksize LT image_dim[1] THEN $
    yend = LONG(blocksize-1 - windowSize[1]) ELSE $
    yend = LONG(image_dim[1]-1)

  readstart = 0l
  readend = readstart+blocksize-1


  xystartend = LONG([xstart, xend, ystart, yend])
  lambda = FLOAT(lambda)
  sigma = FLOAT(1.0)




  ii = 0
  WHILE (readend LE (image_dim[1]-1)) DO BEGIN
    PRINT, "Block Number:  ", STRTRIM(ii, 2)
    tic = SYSTIME(1)

    patInxStart = ULONG64(readstart)*image_dim[0]
    patInxEnd = ULONG64(readend+1)*image_dim[0]-1
    IF KEYWORD_SET(up2) THEN $
      patin = Read_up2(patfile, read_interval=[patInxStart, patInxEnd]) ELSE $
      patin = Read_up1(patfile, read_interval=[patInxStart, patInxEnd])

    nPatBlock = LONG((*patin.header).n_patterns)
    ;;;;;;;
    FOR q = 0, npatblock-1 DO BEGIN
      onePat[0] = (*patin.patterns)[q*patStride:(q+1)*patStride-1]
      patBlock[q*nMask:(q+1)*nMask-1] = FLOAT(onePat[maskindex])
    ENDFOR

    IF KEYWORD_SET(BackGroundSub) THEN BEGIN
      patdimbs = LONG([1,nMask])
      back = FLTARR(nMask)
      nPatBlock = LONG(nPatBlock)
      trash =  CALL_EXTERNAL(library, "cDJR_Background_wrap", $
        patBlock, patdimbs, nPatBlock, back)
    ENDIF
    ;print, "back done"



    ;; need a estimate of the sigma value at each pixel.  Making the assumption that each pixel is
    ; likely sourounded by other pixels with the same pattern.



    sig_search_kern = [[-1,-1], [-1,0], [-1,1], [0,-1], [0,1], [1, -1], [1,0], [1,1]]
    nsearch = N_ELEMENTS(sig_search_kern)/2

    IF N_ELEMENTS(sigma_in) EQ nPoints THEN BEGIN
      sigma = FLOAT(sigma_in[*, readstart:readend])
    ENDIF ELSE BEGIN

      ; s0 = FLTARR(image_dim[0], blocksize, nsearch)
      ; sigma = FLTARR(image_dim[0], blocksize, nsearch)
      ; FOR i=0, nsearch-1 DO BEGIN
      ;   temp = REFORM(patblock - SHIFT(patblock, 0, sig_search_kern[0,i], sig_search_kern[1,i]))
      ;   s0[*,*,i] = SQRT((TOTAL(temp*temp, 1))/(nMask*2.0))
      ; ENDFOR

      ; wh  = where(s0 lt 1.e-6, whcount) ; sometimes EDAX collects the same point twice, giving a bogus value
      ; if whcount gt 0 then s0[wh] = max(s0)+1
      ; sigma = min(s0, dim = 3)

      sigma = FLTARR(image_dim[0], blocksize)
      sigWin = LONG([1,1])
      xystartend_sig  = LONG([0, image_dim[0]-1, 0, blocksize-1])

      trash =  CALL_EXTERNAL(library, "cSigEst_wrap", $
        patBlock, nMask, cImage_dim, xystartend_sig, sigWin, sigma, threads)

    ENDELSE


    ;;;; Here is where NLPAR gets executed !!!!!!
    sigma = FLOAT(sigma)
    dthresh = FLOAT(dthresh)
    trash =  CALL_EXTERNAL(library, "cNLPARV5_wrap", $
      patBlock, nMask, cImage_dim, xystartend, windowSize, sigma, LAMBDA, dthresh, nlparBlock,w0block,  threads)
    ;   0         1         2           3         4          5       6       7        8          9        10



    ;; This is actaully a pretty slow way to rescale, but it does take care of outliers.
    ;; However, rescaling is probably better done during the indexing process.
    IF KEYWORD_SET(rescale) THEN BEGIN
      FOR q = 0, (image_dim[0]*blocksize)-1 DO BEGIN
        temppat = (nlparBlock[q*nMask:(q+1)*nMask-1])
        medianpat = MEDIAN(temppat)
        MAD = MEDIAN(ABS(temppat -medianpat))
        mi = ABS(0.6745*(temppat-medianpat)/MAD)
        whmi = WHERE(mi LE 3.5)
        meanpat = FLOAT(Mean(nlparblock[whmi], /double))
        stddpat = FLOAT(Stddev(nlparblock[whmi], /double))
        temppat -= (meanpat - 4*stddpat)
        temppat *= (bitdepthmax-1)/(8*stddpat)
        nlparBlock[q*nMask:(q+1)*nMask-1] = temppat
      ENDFOR
      ;      medianblck = median(nlparblock)
      ;      MAD = median(abs(nlparblock -medianblck))
      ;      mi = ABS(0.6745*(nlparblock-medianblck)/MAD)
      ;      whmi = where(mi le 3.5)
      ;      meanblock = float(mean(nlparblock[whmi], /double))
      ;      stddblock = float(stddev(nlparblock[whmi], /double))
      ;
      ;      nlparblock -= (meanblock - 3*stddblock)
      ;      nlparblock *= bitdepthmax/(6*stddblock)
      ;nlparblock *= bitdepthmax/MAX(nlparblock)
      ;IF N_ELEMENTS(notmaskindex) GT 1 THEN BEGIN
      ;  onepat[notmaskindex] = mean(nlparblock[0:nmask-1]); 0
    ENDIF
    ;ENDELSE

    nlparBlock >= 0
    nlparBlock <= (bitdepthmax-1)

    maskfiller = Flat(Mean(nlparblock, dim =1))
    ;IF N_ELEMENTS(notmaskindex) GT 1 THEN BEGIN
    ;  onepat[notmaskindex] = Mean(patblock);Mean(onepat[notmaskindex])
    ;ENDIF

    FOR q = 0, (image_dim[0]*blocksize)-1 DO BEGIN
      onePat[maskindex] = ROUND(nlparBlock[q*nMask:(q+1)*nMask-1])
      onePat[notmaskindex] = maskfiller[q]
      patBuffer[q*patStride:(q+1)*patStride-1]= onePat
    ENDFOR
    ;IF bitdepth EQ 8 THEN patBlockOut = BYTE(ROUND(nlparBlock)) ELSE patBlockOut = UINT(ROUND(nlparBlock))
    ;tvscl, background


    ;IF ii EQ 0 THEN writeStart = 0 ELSE writeStart = windowSize[1]+1
    w0[*, readstart+xystartend[2]:readstart+xystartend[3]] = w0block[*, xystartend[2]:xystartend[3]]
    sigma_out[*, readstart+xystartend[2]:readstart+xystartend[3]] = sigma[*, xystartend[2]:xystartend[3]]

    writeStart = xystartend[2]
    writeEnd = xystartend[3]

    patBlockOut = patBuffer[*,*,*, writeStart:writeEnd]
    WRITEU, lunout, patBlockOut


    IF readend EQ (image_dim[1]-1) THEN BEGIN ; set up to exit the While loop
      readend = (image_dim[1])
    ENDIF ELSE BEGIN

      readstart = (readstart)+xystartend[3]+1-windowsize[1]
      xystartend[2] = windowsize[1]

      readend = (readstart+blocksize-1)
      IF readend GE image_dim[1] THEN BEGIN
        readend = (image_dim[1] -1)
        blocksize = readend-readstart+1
        xystartend[3] = blocksize-1
        cImage_Dim[1] = blocksize
        patBlock = MAKE_ARRAY(dimension = [nMask, image_dim[0], blocksize], type = 4, /nozero )
        nlparBlock = MAKE_ARRAY(dimension = [nMask, image_dim[0], blocksize], type = 4, /nozero )
        w0block = MAKE_ARRAY(dimension = [image_dim[0], blocksize], type = 4, /nozero )
      ENDIF
    ENDELSE
    ii +=1
  ENDWHILE



  CLOSE, lunout

  FREE_LUN, lunout

END



