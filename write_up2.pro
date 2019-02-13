pro write_up2, filename, data, version=version
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

; a quick and simple EDAX up2 file writer.
; The default is to write a version 3 up2 file.  The user can option to output a
; version 1.  
; 
; 

  if n_elements(version) eq 0 then version=3

  if (version eq 3) then begin
    get_lun, lun
    openw, lun, filename

    writeu, lun, ulong(3); data.header.version
    writeu, lun, ulong((*data.header).pattern_width)
    writeu, lun, ulong((*data.header).pattern_height)
    writeu, lun, ulong(42); ulong(*(data.header).file_pos)
    writeu, lun, byte((*data.header).extra_patterns)
    writeu, lun, ulong((*data.header).cols)
    writeu, lun, ulong((*data.header).rows)
    writeu, lun, byte((*data.header).hex_flag)
    writeu, lun, double((*data.header).xstep)
    writeu, lun, double((*data.header).ystep)

    sz0 = size(*data.patterns)
    *data.patterns = reform(*data.patterns, $
            (*data.header).pattern_width, (*data.header).pattern_height, $
            (*data.header).n_patterns, /overwrite)
    ;patterns = data.patterns
    sz = size(*data.patterns)
    N = sz[3] ;data.rows*data.cols
    ;print, N
    ;patterns = reform(patterns, sz[1], sz[2], N)
    
    for i=0ll, N-1,1 Do Begin
      temp = uint((*data.patterns)[*,*,i])
      writeu, lun, temp
    endfor
    
    *data.patterns = REFORM(*data.patterns, sz0[1:sz0[0]], /over)
    
    close, lun
    free_lun, lun
    return
  endif

  IF (version EQ 1) THEN BEGIN

    GET_LUN, lun
    OPENW, lun, filename

    ;patterns = data.patterns
    sz = SIZE(*data.patterns)

    header = ULONARR(4)
    header[0] = [1, sz[1], sz[2], 16]

    WRITEU, lun, header

    N = *data.N_Patterns
    PRINT, N
    *data.patterns = REFORM(*data.patterns, sz[1], sz[2], N, /overwrite)
    temp = uintarr(sz[1], sz[2])
    
    ;temp[0,0] = mean(patterns[*,*,0:(1000< (N-1))], dim=3)
    ;writeu, lun, temp ; write out a garbage pattern
    FOR i=0ll, N-1,1 DO BEGIN
      temp[0,0] = uint((*data.patterns)[*,*,i])
      WRITEU, lun, temp
    ENDFOR
    *data.patterns = REFORM(*data.patterns, sz[1:sz[0]], /over)
    CLOSE, lun
    FREE_LUN, lun

    RETURN
  ENDIF

  print, "Error: Version ", strtrim(version,2), " is not supported"


end