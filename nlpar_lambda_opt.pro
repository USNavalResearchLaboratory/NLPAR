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
; Description: This code will generate a set of estimates for the lambda value to use 
; in the NLPAR algorithm.  This can also be used to make calculation of the sigma values
; for each pattern based on nearest neighbor pattern differences.
;
; Arguments (in IDL, arguements need to input in order):
;     patfile (input/required): A string giving the file location of an
;         EDAX up1 or up2 version 3 file that contains the patterns to be processed.
;
;     
; Keywords (in IDL keywords do not need to be in order):
;      lambda (input/default=0.6): a floating point number giving the magnitude of the
;          lambda smoothing factor.
;      
;      up1/up1 (flag/required): Set one of these keywords to 1 to indicate that the
;          input file is either a up1 or up2 file.
;     target_weights (input/optional/default=[0.5, 0.375, 0.25]) : Set to an array that contains
;           the target weights to be optimized for.  Ideally these should be between 1/9 and 1.0, but 
;           the defaults work very well for all the structures that we have tried.     
;     lambda (output/optional): optional return of the lambda values associated with each target weight. 
;     BackGroundSub (flag/default=0): optional flag that if set will do a background subtraction of
;          the patterns based on the "mean" pattern for the block.
;      
;     mask (input/optinal): an optinal input that is set to an array that is the size
;           of one pattern.  Values that are GT 0 will be considered in the NLPAR algorithm.
;           If not set to an valid array a default mask is created that is the standard cicular mask
;
;     sigma_out (output/optional): output array that contains the results of the nearest neighbor sigma estimate
;           if sigma_in is not provided. 
;     RAMWarningoff (flag, optional): This will turn the RAM warning off.  Not advised if you are running on a laptop.
;           If not set, the program will do a rough analysis of the parameters read in from the up1/2 files and will
;           ask the user if they want to continue if the RAM requirement is over 16GB.
;     blocksizein (input, default=150): this is the number of scan lines that will be processed at a time. In total, the
;           blocksize + 2*searchWindow will be read in -- potentially a large amount of RAM.
;     
;
;
;  Example input:
;  nlpar_lambda_opt, './scan7.up2', /up2
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION LambOpt, r
; helper funtion for the optimization routine. 
  COMMON LamOpt_block, d, tW, dthresh0

  ;szd = size(d)

  dw = exp(-(d > dthresh0)/r[0]^2)
  w = total(dw,3)+1.0


  metric = mean(abs(tw-1.0/w))
  return, metric



END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; main program here

PRO nlpar_lambda_opt, patfile, sigma_out=sigma_out, $
  up1 = up1, up2=up2, nColScan = nColScan, nRowScan = nRowScan, $
  dthresh=dthresh,mask = mask, $
  RAMWarningoff = RAMWarningoff, blocksize = blocksize, BackGroundSub=BackGroundSub, $
  target_weights=target_weights,lambda=lambda

  COMMON LamOpt_block, d, tW, dthresh0
  ;  bitdepth = 0ll
  ;  path = FILE_DIRNAME(ROUTINE_FILEPATH("djr_nlparv5", /either))
  ;  library = "patterProcessingC.dylib"
  ;  IF !VERSION.os_name EQ "linux" THEN library = "patterProcessingC_Linux64.so"
  ;  library = path + PATH_SEP() + library
  
  routine_info = routine_info('nlpar_lambda_opt', /source)
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
  If n_elements(blocksize) eq 0 then blocksize = 151l else blocksize = long(blocksize)
  threads = 0
  If n_elements(target_weights) EQ 0 THEN target_weights = [0.5, 0.375, 0.25]
  IF N_ELEMENTS(dthresh) EQ 0 THEN dthresh = 0.0
  dthresh0 = dthresh
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

  
  ;w0 = fltarr(image_dim, /nozero)
  sigma_out = fltarr(image_dim, /nozero)
  nPoints = Product(image_dim, /int)


  pat_dim = LONG([(*onepat.header).pattern_width, (*onepat.header).pattern_height])
  patstride = Product(pat_dim, /int);*bitdepth/8ll

  start_write0 = (*onepat.header).file_pos

  onepat = 0 ; let IDL Auto GC take care of pointer references
  ;blocksize = 151l ;number of rows to analyze at a time.
  blocksize <= image_dim[1]

  blockSizeByte = blockSize*image_dim[0]*pat_dim[0]*pat_dim[1]*(bitDepth/8)
  ;print, blockSizeByte
  print, "Approximate RAM usage: " + strtrim( 7.6*blocksizeByte/1.e9, 2) + "GB"

  IF ((7.6*blockSizeByte) GT 15e9) AND (~KEYWORD_SET(RAMwarningOff)) THEN BEGIN
    PRINT, "You are about to process a block of data that is going to be over 15GB in size."
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

  nMask = LONG(N_ELEMENTS(maskindex))




  ;maskindex =  maskindex[0:*:10]

  ;newPatBlock = make_array(dimension = [pat_dim, image_dim[0], blocksize], type = typecode )
  patBuffer = MAKE_ARRAY(dimension = [pat_dim, image_dim[0], blocksize], type = typecode, /nozero )
  onePat = MAKE_ARRAY(dimension = [pat_dim], type = typecode )
  background = MAKE_ARRAY(dimension = [pat_dim], type = 4 )
  patBlock = MAKE_ARRAY(dimension = [nMask, image_dim[0], blocksize], type = 4, /nozero )
  ;nlparBlock = MAKE_ARRAY(dimension = [nMask, image_dim[0], blocksize], type = 4, /nozero )
  ;w0block = MAKE_ARRAY(dimension = [image_dim[0], blocksize], type = 4, /nozero )


  ; make sure our types are correct for the c-program
  cImage_dim = LONG([image_dim[0], blockSize])
  windowSize = [1,1];LONG(searchWindow)
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
  ;lambda = FLOAT(lambda)
  sigma = FLOAT(1.0)


  lam = list()

  ii = 0
  WHILE (readend LE (image_dim[1]-1)) DO BEGIN
    PRINT, "Block Number:  ", STRTRIM(ii, 2)
    tic = SYSTIME(1)

    patInxStart = ulong64(readstart)*image_dim[0]
    patInxEnd = ulong64(readend+1)*image_dim[0]-1
    IF KEYWORD_SET(up2) THEN $
      patin = Read_up2(patfile, read_interval=[patInxStart, patInxEnd]) ELSE $
      patin = Read_up1(patfile, read_interval=[patInxStart, patInxEnd])

    nPatBlock = long((*patin.header).n_patterns)
    ;;;;;;;
    FOR q = 0, npatblock-1 DO BEGIN
      onePat[0] = (*patin.patterns)[q*patStride:(q+1)*patStride-1]
      patBlock[q*nMask:(q+1)*nMask-1] = FLOAT(onePat[maskindex])
    ENDFOR

    IF KEYWORD_SET(BackGroundSub) THEN BEGIN
      patdim = long([1,nMask])
      back = fltarr(nMask)
      nPatBlock = LONG(nPatBlock)
      trash =  CALL_EXTERNAL(library, "cDJR_Background_wrap", $
        patBlock, patdim, nPatBlock, back)
    ENDIF


     sig_search_kern = [[-1,-1], [-1,0], [-1,1], [0,-1], [0,1], [1, -1], [1,0], [1,1]]
     nsearch = N_ELEMENTS(sig_search_kern)/2

    ;sigma = FLTARR(image_dim[0], blocksize)
    ;sigWin = long([1,1])
    ;xystartend_sig  = long([0, image_dim[0]-1, 0, blocksize-1])
    ;trash =  CALL_EXTERNAL(library, "cSigEst_wrap", $
    ;  patBlock, nMask, cImage_dim, xystartend_sig, sigWin, sigma, threads)

    ;sigma_out[*, readstart+xystartend[2]:readstart+xystartend[3]] = sigma[*, xystartend[2]:xystartend[3]]

    s0 = FLTARR(image_dim[0], blocksize, nsearch)
    sigma = FLTARR(image_dim[0], blocksize, nsearch)
    d = FLTARR(image_dim[0], blocksize, nsearch)
    FOR i=0, nsearch-1 DO BEGIN
       temp = REFORM(patblock - SHIFT(patblock, 0, sig_search_kern[0,i], sig_search_kern[1,i]))
       temp = total(temp^2,1)
       wh  = where(temp lt 1e-6, whcount) ; sometimes EDAX collects the same point twice, giving a bogus value
       if whcount gt 0 then temp[wh] = max(temp)+1
       d[*,*,i] =temp
       ;s0[*,*,i] = SQRT((temp)/(nMask*2.0))
     ENDFOR

    s0 = SQRT((d)/(nMask*2.0))
    sigma = min(s0, dim = 3)

    sigma_out[*, readstart+xystartend[2]:readstart+xystartend[3]] = sigma[*, xystartend[2]:xystartend[3]]
    sigma2 = sigma^2
      
    FOR i=0, nsearch-1 DO BEGIN
      d[*,*,i] -=  nmask*(sigma2+(shift(sigma2, sig_search_kern[0,i], sig_search_kern[1,i])))
      d[*,*,i] /= (sigma2+(shift(sigma2, sig_search_kern[0,i], sig_search_kern[1,i])))*sqrt(2.0*nmask)
    ENDFOR

    r1 = [1.0]
    stepin = [0.001]
    lblock = fltarr( n_elements(target_weights))
    For jj = 0, n_elements(target_weights)-1 Do Begin
      tW = target_weights[jj]
      l = Amoeba(1.e-5, function_name = 'LambOpt', p0 = r1, scale = stepin, function_value = fval, ncalls = ncalls, simplex=simplex)
      lblock[jj] = l
    EndFOr
    lam.add, lblock
    
    
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
  lambdatemp = fltarr(n_elements(target_weights), n_elements(lam))
  for ii = 0, n_elements(lam)-1 Do lambdatemp[*,ii] = lam[ii]

  ;lambda = mean(reform(transpose(lam.toarray()),n_elements(target_weights), ii ), dim = 2)
  lambda = mean(reform(lambdatemp,n_elements(target_weights), ii), dim=2)
  Print, ''
  Print, '****************'
  Print, "Guess at lambda"
  Print, 'Low: ' + String(lambda[0]) + '  High: '+ String(lambda[-1]) + '  Mean: ' + String(mean(lambda))
  ;print, lambda, mean(lambda)

END


