set terminal cairo size 15cm,7cm standalone color colortext header \
  "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'ratio_all.tex'

set view map

set multiplot 

set xrange[0:1]
set yrange[0:1]


set xlabel "from state"
set ylabel "to state"
set xtics 0.2
set ytics 0.2


set origin 0,0
set size 0.48,1


set cblabel "probability $p_{ij}$" offset 0.8,0
set label "a)" at graph -0.2,1
plot "../../data/Jaynes/T_2910.dat" matrix u ($1/60+1./120.):($2/60+1./120.):3 with image notitle
unset label


set origin 0.48,0
set size 0.48,1

set cbrange[-1:1]

set palette defined (-1 "red", 0 "black", 1 "red")
set label "b)" at graph -0.2,1
set cblabel "deviation $\\log_{10}  \\; \\frac{\\pi_i p_{ij}}{\\pi_j p_{ji}}  $" offset 1.3,0,0
plot "../../data/Jaynes/err_2910.dat" matrix u ($1/60+1./120.):($2/60+1./120.):(log10($3)) with image notitle

pause -1
