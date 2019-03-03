set xlabel "CFsteps"
set ylabel "Average waiting time, clocks"
set grid
set xrange [:1000]
set logscale xy
set ytics 10,2,1000000
plot [][300:3000] \
"./sem_stats_titanX.dat" \
   u 6:8 i 0 t "Titan X SEM_CASE_0, morton, rank=5" w lp lw 1 lc 1,\
"" u 6:8 i 1 t "Titan X SEM_CASE_1, morton, rank=5" w lp lw 3 lc 1,\
"" u 6:8 i 1 t "Titan X SEM_CASE_2, morton, rank=5" w lp lw 5 lc 1,\
"" u 6:8 i 3 t "Titan X SEM_CASE_0, usual , rank=5" w lp lw 1 lc 1,\
"" u 6:8 i 4 t "Titan X SEM_CASE_1, usual , rank=5" w lp lw 3 lc 1,\
"" u 6:8 i 5 t "Titan X SEM_CASE_2, usual , rank=5" w lp lw 5 lc 1,\
"./sem_stats_gtx1080.dat" \
   u 6:8 i 0 t "gtx1080 SEM_CASE_0, morton, rank=5" w lp lw 1 lc 2,\
"" u 6:8 i 1 t "gtx1080 SEM_CASE_1, morton, rank=5" w lp lw 3 lc 2,\
"" u 6:8 i 1 t "gtx1080 SEM_CASE_2, morton, rank=5" w lp lw 5 lc 2,\
"" u 6:8 i 3 t "gtx1080 SEM_CASE_0, usual , rank=5" w lp lw 1 lc 2,\
"" u 6:8 i 4 t "gtx1080 SEM_CASE_1, usual , rank=5" w lp lw 3 lc 2,\
"" u 6:8 i 5 t "gtx1080 SEM_CASE_2, usual , rank=5" w lp lw 5 lc 2,\
"./sem_stats_v100.dat" \
   u 6:8 i 0 t "TeslaV100 SEM_CASE_0, morton, rank=5" w lp lw 1 lc 3,\
"" u 6:8 i 1 t "TeslaV100 SEM_CASE_1, morton, rank=5" w lp lw 3 lc 3,\
"" u 6:8 i 1 t "TeslaV100 SEM_CASE_2, morton, rank=5" w lp lw 5 lc 3,\
"" u 6:8 i 3 t "TeslaV100 SEM_CASE_0, usual , rank=5" w lp lw 1 lc 3,\
"" u 6:8 i 4 t "TeslaV100 SEM_CASE_1, usual , rank=5" w lp lw 3 lc 3,\
"" u 6:8 i 5 t "TeslaV100 SEM_CASE_2, usual , rank=5" w lp lw 5 lc 3

pause -1

plot [][300:3000] \
"./sem_stats_titanX.dat" \
   u 6:8 i 6  t "Titan X SEM_CASE_0, morton, rank=6" w lp lw 1 lc 1,\
"" u 6:8 i 7  t "Titan X SEM_CASE_1, morton, rank=6" w lp lw 3 lc 1,\
"" u 6:8 i 8  t "Titan X SEM_CASE_2, morton, rank=6" w lp lw 5 lc 1,\
"" u 6:8 i 9  t "Titan X SEM_CASE_0, usual , rank=6" w lp lw 1 lc 1,\
"" u 6:8 i 10 t "Titan X SEM_CASE_1, usual , rank=6" w lp lw 3 lc 1,\
"" u 6:8 i 11 t "Titan X SEM_CASE_2, usual , rank=6" w lp lw 5 lc 1,\
"./sem_stats_gtx1080.dat" \
   u 6:8 i 6  t "gtx1080 SEM_CASE_0, morton, rank=6" w lp lw 1 lc 2,\
"" u 6:8 i 7  t "gtx1080 SEM_CASE_1, morton, rank=6" w lp lw 3 lc 2,\
"" u 6:8 i 8  t "gtx1080 SEM_CASE_2, morton, rank=6" w lp lw 5 lc 2,\
"" u 6:8 i 9  t "gtx1080 SEM_CASE_0, usual , rank=6" w lp lw 1 lc 2,\
"" u 6:8 i 10 t "gtx1080 SEM_CASE_1, usual , rank=6" w lp lw 3 lc 2,\
"" u 6:8 i 11 t "gtx1080 SEM_CASE_2, usual , rank=6" w lp lw 5 lc 2,\
"./sem_stats_v100.dat" \
   u 6:8 i 6  t "TeslaV100 SEM_CASE_0, morton, rank=6" w lp lw 1 lc 3,\
"" u 6:8 i 7  t "TeslaV100 SEM_CASE_1, morton, rank=6" w lp lw 3 lc 3,\
"" u 6:8 i 8  t "TeslaV100 SEM_CASE_2, morton, rank=6" w lp lw 5 lc 3,\
"" u 6:8 i 9  t "TeslaV100 SEM_CASE_0, usual , rank=6" w lp lw 1 lc 3,\
"" u 6:8 i 10 t "TeslaV100 SEM_CASE_1, usual , rank=6" w lp lw 3 lc 3,\
"" u 6:8 i 11 t "TeslaV100 SEM_CASE_2, usual , rank=6" w lp lw 5 lc 3

pause -1
set ytics 0.00001,10,1000000
set ylabel "Average repeats"
plot [][] \
"./sem_stats_titanX.dat" \
   u 6:($9-1) i 0  t "Titan X SEM_CASE_0, morton, rank=5" w lp lw 1 lc 1,\
"" u 6:($9-1) i 1  t "Titan X SEM_CASE_1, morton, rank=5" w lp lw 3 lc 1,\
"" u 6:($9-1) i 1  t "Titan X SEM_CASE_2, morton, rank=5" w lp lw 5 lc 1,\
"" u 6:($9-1) i 3  t "Titan X SEM_CASE_0, usual , rank=5" w lp lw 1 lc 1,\
"" u 6:($9-1) i 4  t "Titan X SEM_CASE_1, usual , rank=5" w lp lw 3 lc 1,\
"" u 6:($9-1) i 5  t "Titan X SEM_CASE_2, usual , rank=5" w lp lw 5 lc 1,\
"./sem_stats_gtx1080.dat" \
   u 6:($9-1) i 0  t "gtx1080 SEM_CASE_0, morton, rank=5" w lp lw 1 lc 2,\
"" u 6:($9-1) i 1  t "gtx1080 SEM_CASE_1, morton, rank=5" w lp lw 3 lc 2,\
"" u 6:($9-1) i 1  t "gtx1080 SEM_CASE_2, morton, rank=5" w lp lw 5 lc 2,\
"" u 6:($9-1) i 3  t "gtx1080 SEM_CASE_0, usual , rank=5" w lp lw 1 lc 2,\
"" u 6:($9-1) i 4  t "gtx1080 SEM_CASE_1, usual , rank=5" w lp lw 3 lc 2,\
"" u 6:($9-1) i 5  t "gtx1080 SEM_CASE_2, usual , rank=5" w lp lw 5 lc 2,\
"./sem_stats_v100.dat" \
   u 6:($9-1) i 0  t "TeslaV100 SEM_CASE_0, morton, rank=5" w lp lw 1 lc 3,\
"" u 6:($9-1) i 1  t "TeslaV100 SEM_CASE_1, morton, rank=5" w lp lw 3 lc 3,\
"" u 6:($9-1) i 1  t "TeslaV100 SEM_CASE_2, morton, rank=5" w lp lw 5 lc 3,\
"" u 6:($9-1) i 3  t "TeslaV100 SEM_CASE_0, usual , rank=5" w lp lw 1 lc 3,\
"" u 6:($9-1) i 4  t "TeslaV100 SEM_CASE_1, usual , rank=5" w lp lw 3 lc 3,\
"" u 6:($9-1) i 5  t "TeslaV100 SEM_CASE_2, usual , rank=5" w lp lw 5 lc 3

pause -1

plot [][] \
"./sem_stats_titanX.dat" \
   u 6:($9-1) i 6  t "Titan X SEM_CASE_0, morton, rank=6" w lp lw 1 lc 1,\
"" u 6:($9-1) i 7  t "Titan X SEM_CASE_1, morton, rank=6" w lp lw 3 lc 1,\
"" u 6:($9-1) i 8  t "Titan X SEM_CASE_2, morton, rank=6" w lp lw 5 lc 1,\
"" u 6:($9-1) i 9  t "Titan X SEM_CASE_0, usual , rank=6" w lp lw 1 lc 1,\
"" u 6:($9-1) i 10 t "Titan X SEM_CASE_1, usual , rank=6" w lp lw 3 lc 1,\
"" u 6:($9-1) i 11 t "Titan X SEM_CASE_2, usual , rank=6" w lp lw 5 lc 1,\
"./sem_stats_gtx1080.dat" \
   u 6:($9-1) i 6  t "gtx1080 SEM_CASE_0, morton, rank=6" w lp lw 1 lc 2,\
"" u 6:($9-1) i 7  t "gtx1080 SEM_CASE_1, morton, rank=6" w lp lw 3 lc 2,\
"" u 6:($9-1) i 8  t "gtx1080 SEM_CASE_2, morton, rank=6" w lp lw 5 lc 2,\
"" u 6:($9-1) i 9  t "gtx1080 SEM_CASE_0, usual , rank=6" w lp lw 1 lc 2,\
"" u 6:($9-1) i 10 t "gtx1080 SEM_CASE_1, usual , rank=6" w lp lw 3 lc 2,\
"" u 6:($9-1) i 11 t "gtx1080 SEM_CASE_2, usual , rank=6" w lp lw 5 lc 2,\
"./sem_stats_v100.dat" \
   u 6:($9-1) i 6  t "TeslaV100 SEM_CASE_0, morton, rank=6" w lp lw 1 lc 3,\
"" u 6:($9-1) i 7  t "TeslaV100 SEM_CASE_1, morton, rank=6" w lp lw 3 lc 3,\
"" u 6:($9-1) i 8  t "TeslaV100 SEM_CASE_2, morton, rank=6" w lp lw 5 lc 3,\
"" u 6:($9-1) i 9  t "TeslaV100 SEM_CASE_0, usual , rank=6" w lp lw 1 lc 3,\
"" u 6:($9-1) i 10 t "TeslaV100 SEM_CASE_1, usual , rank=6" w lp lw 3 lc 3,\
"" u 6:($9-1) i 11 t "TeslaV100 SEM_CASE_2, usual , rank=6" w lp lw 5 lc 3

pause -1

unset logscale y
set ytics 1e7
set ylabel "Perf, CFsteps/sec"
plot [][] \
"./sem_stats_titanX.dat" \
   u 6:($10/512.*1e6) i 0  t "Titan X SEM_CASE_0, morton, rank=5" w lp lw 1 lc 1,\
"" u 6:($10/512.*1e6) i 1  t "Titan X SEM_CASE_1, morton, rank=5" w lp lw 3 lc 1,\
"" u 6:($10/512.*1e6) i 1  t "Titan X SEM_CASE_2, morton, rank=5" w lp lw 5 lc 1,\
"" u 6:($10/512.*1e6) i 3  t "Titan X SEM_CASE_0, usual , rank=5" w lp lw 1 lc 1,\
"" u 6:($10/512.*1e6) i 4  t "Titan X SEM_CASE_1, usual , rank=5" w lp lw 3 lc 1,\
"" u 6:($10/512.*1e6) i 5  t "Titan X SEM_CASE_2, usual , rank=5" w lp lw 5 lc 1,\
"./sem_stats_gtx1080.dat" \
   u 6:($10/512.*1e6) i 0  t "gtx1080 SEM_CASE_0, morton, rank=5" w lp lw 1 lc 2,\
"" u 6:($10/512.*1e6) i 1  t "gtx1080 SEM_CASE_1, morton, rank=5" w lp lw 3 lc 2,\
"" u 6:($10/512.*1e6) i 1  t "gtx1080 SEM_CASE_2, morton, rank=5" w lp lw 5 lc 2,\
"" u 6:($10/512.*1e6) i 3  t "gtx1080 SEM_CASE_0, usual , rank=5" w lp lw 1 lc 2,\
"" u 6:($10/512.*1e6) i 4  t "gtx1080 SEM_CASE_1, usual , rank=5" w lp lw 3 lc 2,\
"" u 6:($10/512.*1e6) i 5  t "gtx1080 SEM_CASE_2, usual , rank=5" w lp lw 5 lc 2,\
"./sem_stats_v100.dat" \
   u 6:($10/512.*1e6) i 0  t "TeslaV100 SEM_CASE_0, morton, rank=5" w lp lw 1 lc 3,\
"" u 6:($10/512.*1e6) i 1  t "TeslaV100 SEM_CASE_1, morton, rank=5" w lp lw 3 lc 3,\
"" u 6:($10/512.*1e6) i 1  t "TeslaV100 SEM_CASE_2, morton, rank=5" w lp lw 5 lc 3,\
"" u 6:($10/512.*1e6) i 3  t "TeslaV100 SEM_CASE_0, usual , rank=5" w lp lw 1 lc 3,\
"" u 6:($10/512.*1e6) i 4  t "TeslaV100 SEM_CASE_1, usual , rank=5" w lp lw 3 lc 3,\
"" u 6:($10/512.*1e6) i 5  t "TeslaV100 SEM_CASE_2, usual , rank=5" w lp lw 5 lc 3

pause -1

plot [][] \
"./sem_stats_titanX.dat" \
   u 6:($10/512.*1e6) i 6  t "Titan X SEM_CASE_0, morton, rank=6" w lp lw 1 lc 1,\
"" u 6:($10/512.*1e6) i 7  t "Titan X SEM_CASE_1, morton, rank=6" w lp lw 3 lc 1,\
"" u 6:($10/512.*1e6) i 8  t "Titan X SEM_CASE_2, morton, rank=6" w lp lw 5 lc 1,\
"" u 6:($10/512.*1e6) i 9  t "Titan X SEM_CASE_0, usual , rank=6" w lp lw 1 lc 1,\
"" u 6:($10/512.*1e6) i 10 t "Titan X SEM_CASE_1, usual , rank=6" w lp lw 3 lc 1,\
"" u 6:($10/512.*1e6) i 11 t "Titan X SEM_CASE_2, usual , rank=6" w lp lw 5 lc 1,\
"./sem_stats_gtx1080.dat" \
   u 6:($10/512.*1e6) i 6  t "gtx1080 SEM_CASE_0, morton, rank=6" w lp lw 1 lc 2,\
"" u 6:($10/512.*1e6) i 7  t "gtx1080 SEM_CASE_1, morton, rank=6" w lp lw 3 lc 2,\
"" u 6:($10/512.*1e6) i 8  t "gtx1080 SEM_CASE_2, morton, rank=6" w lp lw 5 lc 2,\
"" u 6:($10/512.*1e6) i 9  t "gtx1080 SEM_CASE_0, usual , rank=6" w lp lw 1 lc 2,\
"" u 6:($10/512.*1e6) i 10 t "gtx1080 SEM_CASE_1, usual , rank=6" w lp lw 3 lc 2,\
"" u 6:($10/512.*1e6) i 11 t "gtx1080 SEM_CASE_2, usual , rank=6" w lp lw 5 lc 2,\
"./sem_stats_v100.dat" \
   u 6:($10/512.*1e6) i 6  t "TeslaV100 SEM_CASE_0, morton, rank=6" w lp lw 1 lc 3,\
"" u 6:($10/512.*1e6) i 7  t "TeslaV100 SEM_CASE_1, morton, rank=6" w lp lw 3 lc 3,\
"" u 6:($10/512.*1e6) i 8  t "TeslaV100 SEM_CASE_2, morton, rank=6" w lp lw 5 lc 3,\
"" u 6:($10/512.*1e6) i 9  t "TeslaV100 SEM_CASE_0, usual , rank=6" w lp lw 1 lc 3,\
"" u 6:($10/512.*1e6) i 10 t "TeslaV100 SEM_CASE_1, usual , rank=6" w lp lw 3 lc 3,\
"" u 6:($10/512.*1e6) i 11 t "TeslaV100 SEM_CASE_2, usual , rank=6" w lp lw 5 lc 3


pause -1
