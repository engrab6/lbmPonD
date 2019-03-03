set key left Left reverse
set ylabel "Rho"
set grid
plot [1400:2400][] \
"~/progs/lbmPonD/exact_3.dat" u 1:2 w l lw 3 lc 0 t "Exact",\
"~/progs/lbmPonD/lbm_test3_tau5.dat" u 1:2 i 0 w lp ps 0.5 lw 0.5 lc 1 t "tau=5 Nappr=3 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1.dat" u 1:2 i 0 w lp ps 0.5 lw 0.5 lc 2 t "tau=1 Nappr=3 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1_no5.dat" u 1:2 i 0 w lp ps 0.5 lw 0.5 lc 3 t "tau=1 Nappr=5 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1_d3.dat" u 1:2 i 0 w lp ps 0.5 lw 0.5 lc 4 t "tau=1 Nappr=3 dt=0.00033333333"
pause -1
set key center bottom Left reverse
set ylabel "Velocity"
plot [1400:2400][] \
"~/progs/lbmPonD/exact_3.dat" u 1:4 w l lw 3 lc 0 t "Exact",\
"~/progs/lbmPonD/lbm_test3_tau5.dat" u 1:($3*1000) i 0 w lp ps 0.5 lw 0.5 lc 1 t "tau=5 Nappr=3 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1.dat" u 1:($3*1000) i 0 w lp ps 0.5 lw 0.5 lc 2 t "tau=1 Nappr=3 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1_no5.dat" u 1:($3*1000) i 0 w lp ps 0.5 lw 0.5 lc 3 t "tau=1 Nappr=5 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1_d3.dat" u 1:($3*3000) i 0 w lp ps 0.5 lw 0.5 lc 4 t "tau=1 Nappr=3 dt=0.00033333333"
pause -1
set ylabel "T"
plot [1400:2400][] \
"~/progs/lbmPonD/exact_3.dat" u 1:($3/$2) w l lw 3 lc 0 t "Exact",\
"~/progs/lbmPonD/lbm_test3_tau5.dat" u 1:($4*1000000) i 0 w lp ps 0.5 lw 0.5 lc 1 t "tau=5 Nappr=3 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1.dat" u 1:($4*1000000) i 0 w lp ps 0.5 lw 0.5 lc 2 t "tau=1 Nappr=3 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1_no5.dat" u 1:($4*1000000) i 0 w lp ps 0.5 lw 0.5 lc 3 t "tau=1 Nappr=5 dt=0.001",\
"~/progs/lbmPonD/lbm_test3_tau1_d3.dat" u 1:($4*3000*3000) i 0 w lp ps 0.5 lw 0.5 lc 4 t "tau=1 Nappr=3 dt=0.00033333333"
pause -1

