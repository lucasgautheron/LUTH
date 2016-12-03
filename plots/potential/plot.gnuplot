set terminal epslatex
set output 'plot.tex'

set format xy "10^{%T}"
set logscale xy

set xlabel '$\rho Y_e$ [g.cm$^{-3}$]'
set ylabel '$\mu_e$ [MeV]'

plot '../../Noyaux/potential.res' u 1:2 t 'Calculated' lc rgb 'red' w l, '../../Noyaux/states' u ($1*$2):5 lc rgb 'black' t 'A. Juodagalvis, K.Langanke et al. 2010'
set output
