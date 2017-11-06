pro setcolors, NAMES=names, VALUES=ctindx, $
               NPLOTCOLORS=nplotcol, TABLESIZE=size, $
               $
               DECOMPOSED=decomposed, $
               SYSTEM_VARIABLES=system_variables, $
               GETINFO=getinfo, $
               SILENT=silent, $
               TEST=test, $
               START=start, $
               NOCOLORS=nocolors,$
               $
               PLOTCOLORS=plotcolors, $
               R_PLOT=r_plot, $
               G_PLOT=g_plot, $
               B_PLOT=b_plot, $
               $
               PSEUDO256=pseudo256, $
               $
               _REF_EXTRA=extra
                                ;+
                                ; NAME:
                                ;       SETCOLORS
                                ;
                                ; PURPOSE:
                                ;       To set up device independent line plot colors. Designed to
                                ;       work correctly on both X Window and PostScript devices.
                                ;
                                ; CALLING SEQUENCE:
                                ;       SETCOLORS [, NAMES=string array][, VALUES=array][,
                                ;                 NPLOTCOLORS=variable][, TABLESIZE=variable][,
                                ;                 DECOMPOSED={0|1}][,/SYSTEM_VARIABLES][,
                                ;                 /GETINFO][,/SILENT][,/TEST][,/START][,
                                ;                 PLOTCOLORS=string
                                ;                 array, R_PLOT=array,
                                ;                 G_PLOT=array,
                                ;                 B_PLOT=array][,/PSEUDO256]
                                ;
                                ; INPUTS:
                                ;       None.
                                ;
                                ; OPTIONAL INPUTS:
                                ;       (You can provide your own set of line plot colors, but you
                                ;        must provide all four of the following keywords and they
                                ;        must all have the same size! (See EXAMPLE.) If omitted,
                                ;        a set of 12 basic colors is used.)
                                ;
                                ;       PLOTCOLORS = A string array of
                                ;       user-defined line plot color
                                ;       names.
                                ;
                                ;       R_PLOT = The RED values of line plot colors.
                                ;
                                ;       G_PLOT = The GREEN values of line plot colors.
                                ;
                                ;       B_PLOT = The BLUE values of line plot colors.
                                ;
                                ; OUTPUTS:
                                ;       None.
                                ;
                                ; OPTIONAL OUTPUTS:
                                ;       NAMES = Returns a string array of the line plot color names.
                                ;
                                ;       VALUES = Returns an array of
                                ;       the color table indices
                                ;       (8-bit) or
                                ;                24-bit integers
                                ;                corresponding to the
                                ;                line plot colors.
                                ;
                                ;       NPLOTCOLORS = Returns the
                                ;       number of line plot colors
                                ;       defined.
                                ;
                                ;       TABLESIZE = Returns the number of color table indices.
                                ;
                                ; KEYWORDS:
                                ;       DECOMPOSED = Set this keyword
                                ;       to explicitly use decomposed
                                ;       color
                                ;                    on 24-bit
                                ;                    machines. Has no
                                ;                    effect on devices
                                ;                    which
                                ;                    do not support
                                ;                    decomposed
                                ;                    color. Set this
                                ;                    keyword
                                ;                    to 0 to turn off
                                ;                    decomposed color
                                ;                    on 24-bit
                                ;                    devices.
                                ;
                                ;       /SYSTEM_VARIABLES: Define
                                ;       system variables with the name
                                ;       of each
                                ;                          color and
                                ;                          the value
                                ;                          of the
                                ;                          color table
                                ;                          index
                                ;                          or the
                                ;                          24-bit
                                ;                          integer
                                ;                          corresponding to each
                                ;                          color,
                                ;                          i.e. for
                                ;                          24-bit
                                ;                          color,
                                ;                          !orange=32767L.
                                ;                          This is
                                ;                          very useful
                                ;                          since
                                ;                          system
                                ;                          variables
                                ;                          have
                                ;                          global
                                ;                          scope,
                                ;                          therefore
                                ;                          these
                                ;                          colors can
                                ;                          be
                                ;                          used in any procedure once defined.
                                ;
                                ;       /GETINFO: Do not do anything to the color table, just print
                                ;                 out IDL color information.
                                ;
                                ;       /SILENT: Do not print out IDL
                                ;       color information. Default is
                                ;       to
                                ;                print info.
                                ;
                                ;       /TEST: Display the entire color table and the names of
                                ;              the line plot colors in their respective colors.
                                ;
                                ;       /START: Put the line plot colors at the start of the color
                                ;               table. Default is to store them at the top.
                                ;
                                ;       /NOCOLORS: Do not load line plot colors into color table!
                                ;                  Has no effect if color undecomposed.
                                ;
                                ;       /PSEUDO256: Default is to establish PseudoColor visual using
                                ;                   only the available color indices when the first
                                ;                   window is created.  This keyword forces IDL to
                                ;                   use all 256 colors.
                                ;
                                ; COMMON BLOCKS:
                                ;       COLORS: The IDL color common block.
                                ;
                                ; SIDE EFFECTS:
                                ;       The color table may be modified.  If a window is not open, a
                                ;       pixmap is created then
                                ;       destroyed to establish the
                                ;       visual class
                                ;       and depth.  If the /TEST keyword is set and the device is X
                                ;       Windows, a free window is created to display the color table
                                ;       and line plot colors.
                                ;
                                ; RESTRICTIONS:
                                ;       Only supported by IDL v5.2 or higher.
                                ;
                                ; PROCEDURES CALLED:
                                ;       STRETCH
                                ;
                                ; EXAMPLE:
                                ;       To view the color table and line plot names and colors:
                                ;
                                ;       IDL> setcolors, /test
                                ;
                                ;
                                ;       To define system variables with the names of each color
                                ;       and values of the corresponding color table index or 24-bit
                                ;       long integer:
                                ;
                                ;       IDL> setcolors, /SYSTEM_VARIABLES
                                ;
                                ;
                                ;       To use the line plot colors
                                ;       (set DECOMPOSED=1 if 24-bit
                                ;       display):
                                ;
                                ;       IDL> setcolors, NAMES=cnames, VALUES=cindx
                                ;       IDL> plot, findgen(30),
                                ;       co=cindx[(where(cnames eq
                                ;       'yellow'))[0]]
                                ;       *OR*
                                ;       IDL> plot, findgen(30), co=getcolor('yellow',cnames,cindx)
                                ;
                                ;
                                ;       If you'd like to display an image but you've added line plot
                                ;       colors to the color table, you need to scale the image to
                                ;       avoid using the plotcolors:
                                ;
                                ;       IDL> setcolors, NPLOTCOLORS=nplotcolors, TABLESIZE=tablesize
                                ;       IDL> tv, bytscl(image, top=tablesize-nplotcolors-1)
                                ;       If you load the line plot
                                ;       colors at the start of the
                                ;       color table:
                                ;
                                ;       IDL> setcolors,
                                ;       NPLOTCOLORS=nplotcolors,
                                ;       TABLESIZE=tablesize, /START
                                ;       IDL> tv, bytscl(image,
                                ;       top=tablesize-nplotcolors-1)+nplotcolors
                                ;
                                ;       (Do not use TVSCL!  "If
                                ;       color is important to
                                ;       you (and it almost
                                ;        always is), then you probably
                                ;        never want to use the TVScl
                                ;        command.
                                ;        Instead, you will want to
                                ;        scale your image data
                                ;        yourself, and use
                                ;        the TV command to display it." - David W. Fanning)
                                ;
                                ;
                                ;       To load your own line plot colors:
                                ;
                                ;       IDL> setcolors,
                                ;       PLOTCOLORS=['red','yellow','blue','green'], $
                                ;       IDL> R_PLOT=[255,255,0,0], G_PLOT=[0,255,0,255], $
                                ;       IDL> B_PLOT=[0,0,255,0], NAMES=cnames, VALUES=cindx
                                ;
                                ;       (Color 'WHITE' will be
                                ;       appended to the top of the
                                ;       color table
                                ;        to protect !p.color/!p.background.  If /START is set, color
                                ;        'BLACK' is appended to bottom of color table.)
                                ;
                                ;
                                ;       If you want to make a grayscale image and add colored
                                ;       annotation (when color decomposition is off):
                                ;
                                ;       IDL> setcolors,
                                ;       NPLOTCOLORS=nplotcolors,
                                ;       TABLESIZE=tablesize, $
                                ;       IDL> NAMES=cnames, VALUES=cindx
                                ;       IDL> tv, bytscl(image, top=tablesize-nplotcolors-1)
                                ;       IDL> xyouts, 200, 150, 'M31', /dev, $
                                ;       IDL> co=cindx[(where(cnames eq 'cyan'))[0]]
                                ;       *OR*
                                ;       IDL> xyouts, 200, 150, 'M31',
                                ;       /dev,
                                ;       co=getcolor('cyan',cnames,cindx)
                                ;
                                ;       WARNING!  After opening a PostScript file you should run
                                ;       SETCOLORS!  If you're using a 24-bit visual device,
                                ;       the system variables and color
                                ;       indices will be 24-bit
                                ;       longword
                                ;       integers (with values greater than 256) but PostScript is an
                                ;       8-bit device with only 256
                                ;       colors available, so right
                                ;       after opening
                                ;       a PostScript file, you should
                                ;       run SETCOLORS to re-establish
                                ;       the
                                ;       color system variables color
                                ;       indices. Likewise, when you
                                ;       close
                                ;       a PostScript image, you should
                                ;       run SETCOLORS again in order
                                ;       to
                                ;       re-establish the colors for X Windows.
                                ;
                                ; NOTES:
                                ;       Once a window is open, the
                                ;       visual class is established
                                ;       for the IDL
                                ;       session.  There is no way to
                                ;       change it once it is
                                ;       established.
                                ;       However, if no window
                                ;       is open (there's no way
                                ;       to tell if one has
                                ;       been opened in the past) we
                                ;       open a pixmap and delete it
                                ;       just to
                                ;       establish the visual class and depth.
                                ;
                                ;       I thought at first it would be
                                ;       convenient to store the color
                                ;       table
                                ;       indices of each color in a
                                ;       variable named after the
                                ;       corresponding
                                ;       color and storing each color
                                ;       variable in a common block.
                                ;       However, I
                                ;       was convinced by two very
                                ;       bright folks that common
                                ;       blocks should be
                                ;       avoided, especially for this use.
                                ;
                                ; RELATED PROCEDURES:
                                ;       GETCOLOR
                                ;
                                ; MODIFICATION HISTORY:
                                ;       Written by Tim Robishaw, Berkeley 13 Aug 2001
                                ;       Added /NOCOLORS keyword. TR, Berkeley 03 Oct 2001
                                ;       Added /SYSTEM_VARIABLES keyword. TR, 07 Feb 2002
                                ;       Added /PSEUDO256 keyword. TR, 15 Feb 2002
                                ;-

  on_error, 2

                                ; FIND IDL VERSION NUMBER... 5.2 OR HIGHER...
  if (float(!version.release) lt 5.2) then begin
     message,'This routine is only supported on IDL version 5.2 or higher',/INFO
     return
  endif

                                ; IF NOT IN X WINDOWS, DECOMPOSED DOESN'T MEAN ANYTHING!
  if (!d.name eq 'X') then begin

                                ; SET COLOR DECOMPOSITION EXPLICITLY?
     if (N_elements(DECOMPOSED) gt 0) then device, DECOMPOSED=decomposed

                                ; A WINDOW NEEDS TO HAVE BEEN CREATED
                                ; TO ESTABLISH THE VISUAL TYPE...
     if (!d.window eq -1) then begin
        window, /FREE, /PIXMAP, COLORS=256*keyword_set(PSEUDO256)
        wdelete, !d.window
     endif

                                ; WHAT VISUAL STATES HAVE BEEN SET?
     device, get_visual_name=visual, $
             get_visual_depth=depth, $
             get_decomposed=decomposed

                                ; CAN'T DISPLAY COLORS IF WE HAVE A GRAYSCALE VISUAL...
     if (strupcase(visual) eq 'GRAYSCALE') OR $
        (strupcase(visual) eq 'STATICGRAY') then begin
        message, 'Cannot set colors with '+strtrim(visual,2)+' visual!', /INFO
        if keyword_set(SILENT) $
        then return $
        else goto, output
     endif

  endif else decomposed=0

  size = !d.table_size          ; WHAT IS THE COLOR TABLE SIZE?
  ncol = !d.n_colors            ; HOW MANY COLORS ARE AVAILABLE?

                                ; DO YOU JUST WANT THE CURRENT VISUAL INFO...
  if keyword_set(GETINFO) then goto, output

                                ; THE IDL-DEFINED COLOR COMMON BLOCK...
  common colors, r_orig, g_orig, b_orig, r_curr, g_curr, b_curr

                                ; LOAD THE ORIGINAL COLOR TABLE...
  if (N_elements(r_orig) eq 0) then begin
     r_orig = bytscl(bindgen(size))
     g_orig = r_orig
     b_orig = r_orig
  endif
  tvlct, r_orig, g_orig, b_orig

                                ; DO WE USE PREDEFINED PLOTCOLORS OR USER-SPECIFIED...
  if (N_elements(PLOTCOLORS) eq 0) then begin

                                ; DO WE NOT WANT TO ADD LINE PLOT COLORS...
     if decomposed or (not decomposed AND not keyword_set(NOCOLORS)) then begin

                                ; NAMES OF THE PLOTCOLORS...
        names = [ 'black',  $
                  'red',    $
                  'orange', $
                  'green',  $
                  'forest', $
                  'yellow', $
                  'cyan',   $
                  'blue',   $
                  'magenta',$
                  'purple', $
                  'gray',   $
                  'white']

                                ; GENERATE RGB COLORS FOR COLORED LINE PLOTS...
                                ;              k  r   o  g   f  y  c  b  m   p   g  w
        r_plot = 255b*[0, 1,  1, 0,  0, 1, 0, 0, 1,  0,  0, 1]$
                 +byte([0, 0,  0, 0, 35, 0, 0, 0, 0,153,127, 0])
        g_plot = 255b*[0, 0,  0, 1,  0, 1, 1, 0, 0,  0,  0, 1]$
                 +byte([0, 0,127, 0,142, 0, 0, 0, 0, 50,127, 0])
        b_plot = 255b*[0, 0,  0, 0,  0, 0, 1, 1, 1,  0,  0, 1]$
                 +byte([0, 0,  0, 0, 35, 0, 0, 0, 0,205,127, 0])

     endif

  endif else begin

                                ; ONLY NEEDS TO BE DONE IF UNDECOMPOSED!!!!!
                                ; COLORS PROTECTED OTHERWISE!!!

                                ; MAKE SURE RIGHT # OF RGB & NAMES...
     if (N_elements(R_PLOT)+N_elements(G_PLOT)+N_elements(B_PLOT))/3. $
        ne N_elements(PLOTCOLORS) then $
           message, 'PLOTCOLORS, R_PLOT, G_PLOT, & B_PLOT must have same size!'

                                ; PLOTCOLORS HAD BETTER BE A STRING ARRAY...
     if (size(PLOTCOLORS,/TNAME) ne 'STRING') then $
        message, 'PLOTCOLORS must be of type STRING!'

                                ; IF YOU'RE SILLY AND PUT MORE THAN ONE VALUE OF A COLOR IN
                                ; THE TABLE... OH WELL, YOUR PROBLEM...
     if not keyword_set(START) then begin

                                ; FORCE COLOR 'WHITE' TO BE AT TOP OF COLOR TABLE...
        good = where(strupcase(plotcolors) ne 'WHITE', ngood)
        if (ngood ne 0) then begin
           plotcolors = plotcolors[good]
           names  = [plotcolors,'white']
           r_plot = [r_plot[good],size-1]
           g_plot = [g_plot[good],size-1]
           b_plot = [b_plot[good],size-1]
        endif else names = plotcolors

     endif else begin

                                ; FORCE COLOR 'BLACK' TO BE AT BOTTOM
                                ; OF COLOR TABLE...
        good = where(strupcase(plotcolors) ne 'BLACK', ngood)
        if (ngood ne 0) then begin
           plotcolors = plotcolors[good]
           names  = ['black',plotcolors]
           r_plot = [234,r_plot[good]]
           g_plot = [123,g_plot[good]]
           b_plot = [211,b_plot[good]]
        endif else names = plotcolors

     endelse
  endelse

                                ; HOW MANY PLOTCOLORS HAVE WE CREATED?
  nplotcol = N_elements(names)

                                ; DO YOU EVEN HAVE THIS MANY COLORS AVAILABLE???
                                ; AS CRAZY AS IT SOUNDS, THIS HAS HAPPENED!
  if (ncol le nplotcol) then begin
     message, 'Close some windows! Not enough available colors!', /info
     return
  endif

  if not decomposed then begin
                                ; WE HAVE AN 8-BIT DISPLAY, OR 24-BIT
                                ; WITH COLOR DECOMPOSITION OFF

                                ; STRETCH THE CURRENT COLOR TABLE OVER THE RANGE THAT
                                ; WON'T BE USED BY OUR PLOT COLORS...
                                ; THE CURRENT COLOR TABLE IS STORED IN "ORIG" COLOR VECTORS.
                                ; AFTER THE STRETCH, THE "CURR" COLOR VECTORS WILL CONTAIN
                                ; THE NEW COLOR TABLE...
     if not keyword_set(START) then begin
        start = size-nplotcol
        stretch, 0, start-1
     endif else begin
        start = 0
        stretch, nplotcol, size-1
     endelse

                                ; IF WE WANT PLOT COLORS, ADD THEM TO COLOR TABLE...
     if not keyword_set(NOCOLORS) then begin
                                ; THE COLOR TABLE INDICES OF THE LINE PLOT COLORS...
        ctindx = bindgen(nplotcol)+byte(start)

                                ; STORE THEM IN THE SPECIFIED ROWS OF
                                ; THE NEW COLOR TABLE...
        r_curr[ctindx] = r_plot
        g_curr[ctindx] = g_plot
        b_curr[ctindx] = b_plot

                                ; LOAD THE NEW COLOR TABLE...
        tvlct, r_curr, g_curr, b_curr
     endif

                                ; WE WANT TO LEAVE THE "ORIG" VALUES
                                ; AS THEY WERE WHEN SENT IN!
                                ; IF THEY WERE SET TO "CURR" VALUES
                                ; NOW, IF THIS ROUTINE WERE CALLED
                                ; AGAIN, THE PLOTCOLORS WOULD BE
                                ; STRETCHED INTO THE COLOR TABLE!

  endif $
                                ; OTHERWISE WE HAVE A 24-BIT DISPLAY...
  else ctindx = long(r_plot) + ishft(long(g_plot),8) + ishft(long(b_plot),16)

                                ; SET THESE COLOR NAMES AS SYSTEM VARIABLES...
  if keyword_set(SYSTEM_VARIABLES) then $
     for i = 0, N_elements(names)-1 do begin
     success = execute('defsysv, "!'+names[i]+'", long(ctindx[i]), 0')
     if not(success) $
     then message, '!'+strupcase(names[i])+' not changed.', /cont
  endfor

                                ; DISPLAY THE COLOR TABLE AND SEE IF THE COLORS ARE CORRECT...
  if keyword_set(TEST) then begin
     if (!d.name eq 'X') then begin
        window, /FREE, /PIXMAP
        wn = !d.window & wdelete, wn
        window, /FREE, YSIZE=2*size, XSIZE=size, $
                TITLE='IDL '+strtrim(wn,2)+': Plot Colors'
     endif
     if (!d.name eq 'PS') then $
        device, xsize=3, ysize=6, xoff=(8.5-3)/2., yoff=(11.-6)/2
     device, set_font='Helvetica Bold', /tt_font
     tv, matrix_multiply(intarr(size)+1,[indgen(size),indgen(size)])
     prefix = keyword_set(SYSTEM_VARIABLES) ? '!!' : ''
     for i = 0.0, nplotcol-1.0 do begin
        printcolor = 'xyouts, 0.45, '+strtrim((i+1)/(nplotcol+1),2)$
                     +', '''+prefix+names[i]+''', font=1, /norm, charsize=2.0, co='$
                     +strtrim(long(ctindx[i]),2)
        useless = execute(printcolor)
     endfor
  endif

                                ; OUTPUT THE COLOR TABLE INFO...
  output:
  if not keyword_set(SILENT) OR keyword_set(GETINFO) then begin
     print, 'Display Device  : ', !d.name, format='(%"\N",2A)'
     if (!d.name eq 'X') then begin
        print, 'Visual Class    : ', visual
        print, 'Visual Depth    : ', strtrim(depth,2)+'-Bit'
     endif
     print, 'Color Table Size: ', strtrim(!d.table_size,2)
     print, 'Number of Colors: ', strtrim(!d.n_colors,2)
     print, 'Decomposed Color: ', strtrim(decomposed,2), format='(2A,%"\N")'
  endif

end                             ; setcolors

