function read_up1, filename, read_interval=read_interval
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
; a file reader for a .up1 from EDAX
; version 1 or 3
; Reeturns a structure with header and patterns in an array of [patdimX, patdimY, N]
;
; inputs:
;         filename: string that describes the file location of a up1 file
;         read_interval (optional input): 2 element vector indicating a pattern interval to start/stop, 
;                                         encoded as a 1-D index of the pattern numbers. 
;           For example, to read the first 5 patterns in the file: read_interval = [0,4].  
;           The default is to read all patterns in a file. 
; 
;                       
;         header (optional input/output): If empty, on the read of the file this will be assigned with the 
;                                         relevant header information for the .up1 in a structure.
;                                         On subsequent reads, if inputted, reading the header in the file 
;                                         will be skipped, and the data in the header structure is used instead.  
; 
;  version 1 will have dummy values for some of the header info

; open the file
  if n_elements(filename) eq 0 then filename = ''
  if filename eq '' then begin ; a way to get a dummy header 
    header = {version:3, pattern_width:0, pattern_height:0, $
      file_pos:42, extra_patterns:0, cols:0ll, n_patterns:0ll, $
      rows:0ll, hex_flag:0, xstep:1.0d, ystep:1.0d, $
      index_start:0ll, index_end:0ll}

    patterns =[0]
    return, {header:ptr_new(header, /no_copy), patterns:ptr_new(data, /no_copy)}
  endif
  
  
  file_info = file_info(filename)
  get_lun, lun
  openr, lun, filename

  ; read the first number to determine the version number
  version = ulong(0)
  readu, lun, version
  bitdepth = 8
  ; does the header contain any information?
  ; skip reading the header
  if n_elements(header) ne 0 then begin
    header=header
    N = header.n_patterns
    pat_width = header.pat_width
    pat_height = header.pat_height
    file_pos0 = header.file_pos
  endif else begin
    if version eq 1 then begin
      sz = ulonarr(3)
      readu, lun, sz
      pat_width = sz[0]
      pat_height = sz[1]
     ; N = Round(double(file_info.size-16ll)/(sz[0]*sz[1]))-1
      N = Round((file_info.size-16ll)/(sz[0]*sz[1]))
      ; Add placeholders for version 1
      extra_patterns = 0
      cols = 0
      rows = 0
      file_pos0 = 16l ;63
      hex_flag = 0
      step_size = [0.0d, 0.0d]

      
    endif else if version eq 3 or version eq 2 then begin 
      sz = ulonarr(3)
      readu, lun, sz
      pat_width = sz[0]
      pat_height = sz[1]
      file_pos0 = sz[2]

      extra_patterns = byte(1)
      readu, lun, extra_patterns

      cols_rows = ulonarr(2)
      readu, lun, cols_rows
      cols=cols_rows[0]
      rows=cols_rows[1]

      hex_flag = byte(1)
      readu, lun, hex_flag

      step_size = dblarr(2)
      readu, lun, step_size

      if hex_flag eq 0 then N = cols*rows
      if hex_flag eq 1 and extra_patterns eq 1 then N = cols*rows
      if hex_flag eq 1 and extra_patterns eq 0 then N = floor(cols*rows - rows / 2)    
    endif    
    
   header = {version:version, pattern_width:sz[0], pattern_height:sz[1], $
      file_pos:file_pos0, extra_patterns:extra_patterns, cols:cols, n_patterns:N, $
      rows:rows, hex_flag:hex_flag, xstep:step_size[0], ystep:step_size[1], $
      index_start:0, index_end:N-1}
    
  endelse

  ; now do the interval part
  if n_elements(read_interval) eq 0 then read_interval = [0,N-1]
  if n_elements(read_interval) eq 1 then read_interval = rebin(reform(read_interval), 2)
  nPatRead = read_interval[1] - read_interval[0] + 1
  header.n_patterns = nPatRead
  data = uintarr(pat_width, pat_height, nPatRead)
  header.index_start = read_interval[0]
  header.index_end = read_interval[1]
  stride = header.pattern_height*header.pattern_width*bitdepth/8ll
  file_pos = file_pos0+read_interval[0]*stride

  ; ready to start reading patterns
  point_lun, lun, file_pos
  
  line = bytarr(sz[0], sz[1])  
  for ii=read_interval[0],read_interval[1] do begin
    readu, lun, line
    data[*,*,ii-read_interval[0]] = line
  endfor
  
  close, lun
  free_lun, lun

  return,  {header:ptr_new(header, /no_copy), patterns:ptr_new(data, /no_copy)}
  
end


