set terminal png size 680,600
set pm3d map
s=0.3
unset key
set rmargin 0
set lmargin 0
set bmargin 0
set tmargin 0
set yrange [-5.6 to 5.6]
set xrange [-11.6 to -0.4]
splot "rate_table.lst" using 1:2:3 with pm3d title 'rate', "field_table.lst" using 1:2:(0):(s*$3):(s*$5):(0) with vectors linecolor rgb "#000000" title 'B/|B|'
