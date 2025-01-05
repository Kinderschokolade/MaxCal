
set terminal cairo size 15cm,5cm standalone color colortext header \
   "\\renewcommand\\familydefault{\\sfdefault}\\usepackage{cmbright}"

set output 'ISAW.tex'

set style line 2 lc rgb '#0060ad' lt 1 lw 5 pt 10 ps 0.5   #  blue
set style line 3 lc rgb '#dd181f' lt 1 lw 5 pt 5 ps .5   # red
set style line 4 lc rgb '#66A61E' lt 1 lw 5 pt 8 ps .5   #  green
set style line 6 lc rgb '#ff8800' lt 1 lw 5 pt 7 ps .5   #  orange
set style line 7 lc rgb '#9FAFDF' lt 1 lw 5 pt 8 ps 0.5   #royal blue
set style line 8 lc rgb '#1E90FF' lt 1 lw 5 pt 8 ps 0.5   #  slate blue
set style line 1 lc rgb '#F897C5' lt 1 lw 5 pt 8 ps 0.5   # taffy
set style line 5 lc rgb 'black' lt 1 lw 5 pt 5 ps 0.5   # black

set key outside right

set xlabel "$E / J$"
set ylabel "$P(E)$"
set xtics 30
set ytics 0.01


set xrange[-110:0]

plot "../../data/Jaynes/Ehist_4.dat" u (-$2):($3/50000) ls 2  w l title "Simulation$: T=4 J/k_{\\mathrm{B}} $",\
"../../data/Jaynes/Ehist_5.dat" u (-$2):($3/75000) w l ls 3 title "Simulation$: T=5 J/k_{\\mathrm{B}} $",\
"../../data/Jaynes/Ehist_8.dat" u (-$2):($3/75000) w l ls 4 title "Simulation$: T=8 J/k_{\\mathrm{B}} $",\
"../../data/Jaynes/Ehist_8.dat" u (-$2):($3*exp((1/4.-1/8.)*$2)/9187671) w p ls 2 pt 6 ps 0.3 title "$ T=8 J/k_{\\mathrm{B}} \\rightarrow  T=4 J/k_{\\mathrm{B}} $" ,\
"../../data/Jaynes/Ehist_8.dat" u (-$2):($3*exp((1/5.-1/8.)*$2)/1097800) w p ls 3 ps 0.3  title "$ T=8 J/k_{\\mathrm{B}} \\rightarrow  T=5 J/k_{\\mathrm{B}} $" ,\

pause -1


