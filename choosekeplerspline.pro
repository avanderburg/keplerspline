function choosekeplerspline, t, f, range = range, num = num

  if 1 - keyword_set(num) then num = 20
  if 1 - keyword_set(range) then begin
    range = [.25, 20]
  endif
  
  ttot = max(t) - min(t)
  
  nds = logspace(range[0], range[1], num)
  ses = dblarr(n_elements(t), num)
  
  npars =  ttot / nds + 3 + 1 ; number of knots + order of spline + 1
  rmses = dblarr(num)
  ngood = rmses
  chisqs = rmses
  
  errors = median(abs(f[1:*] - f[0:n_elements(f)-1])) * 1.48 /sqrt(2) ; this is the point to point scatter
  
  errors *= 0.8
  
  for i = 0, num-1 do begin
    ses[*,i] = keplerspline(t, f, nd = nds[i], rms = rms, goodind = goodind, bpndays = 0.75)
    chisqs[i] = total((f(goodind) - ses[goodind,i])^2 / errors^2)
    rmses[i] = rms
    ngood[i] = n_elements(goodind)
    
  endfor
  
  bic = chisqs + npars * alog(ngood)
  redchisq = chisqs / (ngood - npars)
  closeall
  mbic = min(bic, ind)
  
  plot, nds, bic, ps = -8,/yno,/xl,/yl
  oplotvert, nds(ind), color = cgcolor('red')
  figure
  plot, t, f, ps = 3,/yno, xr = [800,1200]
  oplot, t, ses[*,ind], color = cgcolor('red'), ps = 3
  print, nds[ind]
  

  return, ses[*,ind]
  
  
  
end

a = getkoidata(351)

s = choosekeplerspline(a.t, a.f)

stop

end


