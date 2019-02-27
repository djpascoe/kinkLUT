pro load_kinklut, fname
  ; Loads kink damping profile lookup table data into common block
  ; David Pascoe 28/11/2018
  compile_opt idl2
  common kinklut, x1, y1, f1s, g1, triangles1, n_amp, xvals, xmax, version
  if n_elements(fname) eq 0 then fname='kinklut_v1_0.sav'
  restore,fname
  message,'kinkLUT version '+version+' loaded',/informational
end


function damping_lut, x, zeta, epsilon
  ; Damping profile for kink oscillations based on lookup table
  ; (linear inhomogeneous layer profile)
  ; x: normalised times for damping profile (t/P)
  ; zeta: density contrast ratio
  ; epsilon: normalised inhomogeneous layer width (l/R)
  ; David Pascoe 28/11/2018
  compile_opt idl2
  common kinklut, x1, y1, f1s, g1, triangles1, n_amp, xvals, xmax, version

  amps=dblarr(n_amp)
  for i=0,n_amp-1 do amps[i]=griddata(x1,y1,f1s[i],/linear, $
    xout=[zeta],yout=[epsilon],triangles=triangles1)
  ii=where(x le xmax,/null)
  d_lut=spline(xvals,amps,x[ii],/double)

  ii=where(x gt xmax,count)
  if count gt 0 then begin
    ; extrapolate damping profile assuming exponential asymptotic state
    extrap_grad=griddata(x1,y1,g1,/linear, $
      xout=[zeta],yout=[epsilon],triangles=triangles1)
    d_extrap=d_lut[-1]*exp(extrap_grad[0]*(x[ii]-x[ii[0]-1]))
    d_lut=[d_lut,d_extrap]
  endif

  return,d_lut
end


function damping_exp, x, zeta, epsilon
  ; Exponential damping profile for kink oscillations
  ; (thin boundary approximation and linear inhomogeneous layer profile)
  ; x: normalised times for damping profile (t/P)
  ; zeta: density contrast ratio
  ; epsilon: normalised inhomogeneous layer width (l/R)
  ; David Pascoe 28/11/2018
  compile_opt idl2

  kappa=(zeta - 1.d)/(zeta + 1.d)
  q_exp=4.d/(!dpi^2*kappa*epsilon)
  d_exp=exp(-x/q_exp)
  return,d_exp
end


function damping_gau, x, zeta, epsilon
  ; Gaussian damping profile for kink oscillations
  ; (thin boundary approximation and linear inhomogeneous layer profile)
  ; x: normalised times for damping profile (t/P)
  ; zeta: density contrast ratio
  ; epsilon: normalised inhomogeneous layer width (l/R)
  ; David Pascoe 28/11/2018
  compile_opt idl2

  kappa=(zeta - 1.d)/(zeta + 1.d)
  q_gau=2.d/(!dpi*kappa*sqrt(epsilon))
  d_gau=exp(-0.5d*(x/q_gau)^2)
  return,d_gau
end


function damping_gdp, x, zeta, epsilon
  ; General damping profile (GDP) for kink oscillations
  ; (thin boundary approximation and linear inhomogeneous layer profile)
  ; The GDP is a Gaussian profile up to a switch time and then exponential
  ; as proposed in A&A 551, A40 (2013) DOI: 10.1051/0004-6361/201220620
  ; x: normalised times for damping profile (t/P)
  ; zeta: density contrast ratio
  ; epsilon: normalised inhomogeneous layer width (l/R)
  ; David Pascoe 28/11/2018
  compile_opt idl2

  kappa=(zeta - 1.d)/(zeta + 1.d)
  q_gau=2.d/(!dpi*kappa*sqrt(epsilon))  ; tau_g / P
  q_exp=4.d/(!dpi^2*kappa*epsilon)  ; tau_d / P
  t_s=q_gau^2/q_exp  ; switch time normalised to period
  d_gdp=exp(-0.5d*(x/q_gau)^2)*(x le t_s) + $
    exp(-0.5d*(t_s/q_gau)^2)*exp(-(x-t_s)/q_exp)*(x gt t_s)
  return,d_gdp
end


pro kinklut_example
  ; Example profiles using kinkLUT and comparison to analytical profiles based
  ; on the thin boundary (TB) approximation. All profiles are for an inhomogeneous
  ; layer with a linear density profile. For further information see:
  ;  "Coronal loop seismology using standing kink oscillations with a lookup table"
  ;  David James Pascoe, Alan William Hood, Tom Van Doorsselaere
  ;  Frontiers in Astronomy and Space Sciences
  ; David Pascoe 28/11/2018
  compile_opt idl2

  load_kinklut

  period=5.d  ; minutes
  time=dindgen(200)/199.d*30.d  ; minutes
  x=time/period


  ; thin boundary example
  epsilon=0.1d
  zeta=6.d

  ; sinusoidal oscillation with kinkLUT damping profile
  y=sin(2.d*!dpi*x)*damping_lut(x,zeta,epsilon)

  window,0
  plot,time,y,yrange=[-1.1,1.1],title='!7e!X = 0.1, !7f!X = 6.0', $
    xtitle='time (mins)',ytitle='amplitude',/xstyle,/ystyle
  ; overplot TB damping envelopes for comparison
  oplot,time,-damping_exp(x,zeta,epsilon),linestyle=1
  oplot,time,-damping_gau(x,zeta,epsilon),linestyle=2
  oplot,time,damping_gdp(x,zeta,epsilon),linestyle=0


  ; thick boundary example
  epsilon=1.d
  zeta=1.5d

  ; sinusoidal oscillation with kinkLUT damping profile
  y=sin(2.d*!dpi*x)*damping_lut(x,zeta,epsilon)

  window,1
  plot,time,y,yrange=[-1.1,1.1],title='!7e!X = 1.0, !7f!X = 1.5', $
    xtitle='time (mins)',ytitle='amplitude',/xstyle,/ystyle
  ; overplot TB damping envelopes for comparison
  oplot,time,-damping_exp(x,zeta,epsilon),linestyle=1
  oplot,time,-damping_gau(x,zeta,epsilon),linestyle=2
  oplot,time,damping_gdp(x,zeta,epsilon),linestyle=0

end
