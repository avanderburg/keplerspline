function keplerspline2test, t, f, ndays=ndays, maxiter = maxiter, rms = rms, yplot = yplot, goodout = goodout, thisincludeinds = thisincludeinds



  t2 = (t - min(t)) / (max(t) - min(t))
  bksp = ndays / (max(t) - min(t))
  
  lastgood = findgen(n_elements(t2))
  if n_elements(thisincludeinds) ne 0 then lastgood = thisincludeinds
  
  for i = 0, maxiter do begin
  
  
    thist = t2(lastgood)
    thisf = f(lastgood)
    
    sset = bspline_iterfit(thist,thisf,nord=4,maxiter=0,bkspace=bksp)
    
    bg = bspline_valu(t2,sset)
    
    ;    if total(1 - finite(bg)) gt 0 then begin
    ;      stop
    ;      return, dblarr(n_elements(bg)) + 1d
    ;    endif
    
    junk = robust_mean(f - bg, 3, goodind = good)
    
    good = cgsetintersection(good, thisincludeinds)
    if n_elements(good) gt 1 then rms = stdev(bg(good) - f(good))
    if keyword_set(yplot) then begin
    
      if i le yplot then begin
        small = 0.7
        large = 1
        title = 'Iteration '+ trim(i + 1, 1)
        if i eq 0 then title = 'EPIC 60019950: Iteration '+ trim(i + 1, 1)
        if i eq 2 then xtitle = 'BJD - 2454833'
        
        plot, t, f, ps=8, /ynozero,/nodata, title = title,ytitle = 'Relative Flux', xtitle = xtitle
        oplot, t, f, ps=8, symsize = large
        oplot, t, f, ps=8, symsize = small, color = cgcolor('red')
        oplot, t(lastgood), f(lastgood), ps=8;, color = cgcolor('black'), symsize = large
        oplot, t, bg, thick = 15;, color = cgcolor('black')
        oplot, t, bg, thick = 12, color = cgcolor('orange')
      endif
    endif
    if n_elements(lastgood) eq n_elements(good) then begin
      if total(lastgood ne good) eq 0 then begin
        if keyword_set(yplot) then print, 'converged in ', i + 1 , ' iterations'
        break
      endif
    end
    
    
    lastgood = good
    
  end
  goodout = good
  return, bg
  
end
function keplersplinetest, tin, fin, ndays=ndays, maxiter = maxiter, rms = rms, yplot = yplot, $
  breakpin = breakpin, goodind = goodind, bpndays = bpndays, includeindsin = includeindsin
  
  
  t = tin
  f = fin
  difft = diff(t)
  resort = 0
  if keyword_set(breakpin) then breakp = breakpin
  
  if keyword_set(includeindsin) then includeinds = includeindsin
  if 1 - keyword_set(includeindsin) then includeinds = lindgen(n_elements(t))
  if 1 - keyword_set(ndays) then ndays = 1.5
  if 1 - keyword_set(bpndays) then bpndays = ndays
  
  if total(difft lt 0) ne 0 then begin

    ;print, 'FIX ME!!!'
    ;print, 'The time array is not sorted. Keplerspline will break. '

    resort = 1; so that later in the end we can check
    sor = sort(t)
    t = t[sor]
    f=f[sor]
    if n_elements(breakp) gt 0 then begin
      newbreakp = indgen(n_elements(breakp))
      for i = 0, n_elements(breakp)-1 do begin
        x = where(sor eq breakp[i])
        if x[0] ne -1 then newbreakp[i] = x[0]
      endfor
      breakp = newbreakp
    end
    
    ; do includeinds here too
    
    isincluded = intarr(n_elements(t))
    isincluded[includeinds] = 1
    isincluded = isincluded[sor]
    includeinds = where(isincluded)


    difft = diff(t)
    


  endif
  
  if n_elements(maxiter) eq 0 then maxiter = 5
  gaps = where(difft gt bpndays)
  if keyword_set(breakp) then begin
    if gaps[0] eq -1 then gaps = [breakp]
    if gaps[0] gt -1 then begin
      gaps = [gaps, breakp]
      gaps = gaps[UNIQ(gaps, SORT(gaps))]
    endif
  endif
  if gaps[0] eq -1 then begin
    s = keplerspline2test(t, f, ndays = ndays, maxiter = maxiter, yplot=yplot, good = good, thisincludeinds=includeinds) ; make sure we separate the continuum normalizations between gaps in data
    goodind = good
  endif
  if gaps[0] gt -1 then begin
    s = dblarr(n_elements(t))
    goodind = []
    gaps = [0,gaps,n_elements(t)]
    for i =1, n_elements(gaps)-1 do begin
      ;print, max(diff(k.t[gaps[i-1]:gaps[i]-1]))
      thisii = includeinds[where(includeinds ge gaps[i-1] and includeinds le gaps[i]-1)] - gaps[i-1]
      s[gaps[i-1]:gaps[i]-1] = keplerspline2test(t[gaps[i-1]:gaps[i]-1], f[gaps[i-1]:gaps[i]-1], ndays = ndays,$
         maxiter = maxiter,yplot=yplot, good = good, thisincludeinds = thisii)
      goodind = [goodind, good + gaps[i-1]]
    endfor
  end
  
  rms = stdev(s(goodind) - f(goodind))
  
  if resort then begin

    ssor = sort(sor)
    s = s[ssor]
    ;print, 'Fix Goodinds!'
    isgoodinds = intarr(n_elements(t))
    isgoodinds[goodind] = 1
    isgoodinds = isgoodinds[ssor]
    goodind = where(isgoodinds)

  endif
  
  return, s
end

