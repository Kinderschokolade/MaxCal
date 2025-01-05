#set terminal svg 
#set output  "bell.svg"

set xrange[-0.5:2.5]
unset border
unset tics
unset key

U(x) = x*x
V(x) = 0.975 + 0.3*(x-1.5)*(x-1.5)

f = 0.3

S(x) = f*(x-1.8)*(x-1.8)

plot  x< 1 ? U(x) : V(x) lc "grey" lw 3 , x< 1 ? U(x)+S(x) : V(x)+S(x) lc "black" lw 3 


pause -1 
