reset
# 5in,3.75in = size
# solid = no dashed lines, only solid ones
# 10 = size of font to use
set term epslatex standalone solid color size 12cm,8cm 12

set output 'Fig5b.tex'
set xlabel '$h$'
set ylabel '$F$ [$k_{\mathrm{B}}T$]' offset 1

#set key left Left reverse spacing 1.3
set key off

set label '$-\beta\Delta G_{0\alpha} = 5.8$, $8.7$, $11.6$, $14.5$, $17.4$, $20.3$'  at 0.3,5
set arrow from 0.9,4 to 0.9,-5 lw 2 lc rgbcolor "black"

plot [0:2][-7:7] \
     'fig5b.txt' index 0 using 1:2 with lines lw 4 lc rgbcolor "black", \
     'fig5b.txt' index 1 using 1:2 with lines lw 4 lc rgbcolor "black", \
     'fig5b.txt' index 2 using 1:2 with lines lw 4 lc rgbcolor "black", \
     'fig5b.txt' index 3 using 1:2 with lines lw 4 lc rgbcolor "black", \
     'fig5b.txt' index 4 using 1:2 with lines lw 4 lc rgbcolor "black"

unset output
set terminal X

!epstopdf Fig5b-inc.eps
!pdflatex Fig5b.tex

!rm Fig5b.aux
!rm Fig5b.log
!rm Fig5b.tex
!rm Fig5b-inc.eps
!rm Fig5b-inc.pdf
!xpdf Fig5b.pdf&
