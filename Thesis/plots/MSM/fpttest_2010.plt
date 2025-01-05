set logscale x

plot  "../../data/MSM/fpt_2010.dat" u 0:1:7 w err, \
 "" u 0:2:8 w err,\
 "" u 0:3:9 w err ,\
 "" u 0:4:10 w err ,\
 "" u 0:5:11 w err ,\
 "" u 0:6:12 w err ,\
 "../../data/MSM/markovana_200_2010.dat" u ($0*4):($2/4) w l ls 8 ,\
 "" u ($0*4):($3/4) w l ls 8,  "" u ($0*4):($4/4) w l ls 8 ,\
  "" u ($0*4):($6/4) w l ls 8 , "" u ($0*7):($3/4) w l ls 8,  "" u ($0*4):($8/4) w l ls 8


pause -1
