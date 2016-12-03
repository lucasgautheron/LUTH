set terminal epslatex
set output 'plot.tex'

set format y "10^{%T}"
set logscale y

set xlabel '$Q$ [MeV]'
set ylabel '$\lambda_{\mbox{ec}}$ [s$^{-1}$]'

plot '../../Noyaux/rates.res' u 1:4 t '(A)' lc rgb 'red' w l, '../../Noyaux/rates.res' u 1:5 lc rgb 'black' t '(B)' w l
set output
