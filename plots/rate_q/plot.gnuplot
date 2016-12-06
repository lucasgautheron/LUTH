set terminal epslatex
set output 'plot.tex'

set format y "10^{%T}"
set logscale y

set xlabel '$Q$ [MeV]'
set ylabel '$\lambda_{\mbox{ec}}$ [s$^{-1}$]'

plot '../../Noyaux/rates.res' u 1:4 t '(A) analytical' lc rgb 'red' w l, '' u 1:5 lc rgb 'black' t '(B)' w l, '' u 1:2 t '(A) Bruenn' w l, '' u 1:3 t '(B) Bruenn' w l,
'../../Noyaux/mmc01.txt' u 
set output
