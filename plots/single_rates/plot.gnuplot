set terminal epslatex
set output 'plot.tex'

set format xy "10^{%T}"
set logscale xy

set xlabel '$\lambda_{\mbox{ec}}$ [s$^{-1}$] Langanke'
set ylabel '$\lambda_{\mbox{ec}}$ [s$^{-1}$] '

set xrange [1e-4:1e6]
set yrange [1e-4:1e6]

set pointsize 0.3333

plot '../../Noyaux/single_rates.res' u 5:(($3>0.1 && $3<=0.5)?$7:1/0) lc rgb 'green' t '0.1 MeV $\leq T \leq$ 0.5 MeV', '' u 5:(($3>0.5 && $3<=1)?$7:1/0) lc rgb 'blue' t '0.5 MeV $\leq T \leq$ 1 MeV', '' u 5:(($3>1)?$7:1/0) lc rgb 'red' t '1 MeV $\leq T$ ', x w l lc rgb 'black' notitle
set output
