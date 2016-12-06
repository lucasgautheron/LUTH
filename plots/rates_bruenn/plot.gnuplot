set terminal epslatex
set output 'plot.tex'

set format xy "10^{%T}"
set logscale xy

set xlabel '$n_b$ [fm$^{-3}$]'
set ylabel '$\langle\lambda_{\mbox{ec}}\rangle$ [s$^{-1}$]'

plot '../../Noyaux/total_rates_15.res' u 2:4 t '15 $M_\odot$ A. Juodagalvis, K.Langanke et al. 2010' lc rgb 'red' w l, '' u 2:7 lc rgb 'red' t '$A \geq 20$, 15 $M_\odot$', '../../Noyaux/total_rates_25.res' u 2:4 t '25 $M_\odot$ A. Juodagalvis, K.Langanke et al. 2010' lc rgb 'blue' w l, '' u 2:7 lc rgb 'blue' t '$A \geq 20$, 25 $M_\odot$'
set output
