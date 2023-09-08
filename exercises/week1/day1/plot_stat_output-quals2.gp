
            set terminal png size 600,400 truecolor
            set output "plot_stat_output-quals2.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set multiplot
            
            set yrange [0:50]
            set ylabel "Quality"
            set xlabel "Cycle (fwd reads)"
            plot '-' using 1:2:3 with filledcurve lt 1 lc rgb "#cccccc" t '25-75th percentile' , '-' using 1:2 with lines lc rgb "#000000" t 'Median', '-' using 1:2 with lines lt 1 t 'Mean'
        1	15	30
2	15	30
3	15	30
4	16	30
5	17	30
6	17	30
7	18	30
8	19	30
9	21	31
10	22	31
11	22	31
12	23	31
13	24	31
14	24	31
15	24	31
16	24	31
17	25	31
18	25	31
19	25	31
20	25	31
21	25	31
22	25	31
23	24	31
24	25	31
25	24	31
26	24	31
27	23	31
28	23	31
29	20	30
30	20	30
31	19	30
32	19	30
33	18	30
34	18	30
35	17	30
36	18	30
end
1	25
2	25
3	26
4	26
5	26
6	26
7	26
8	27
9	28
10	29
11	28
12	29
13	29
14	29
15	29
16	29
17	29
18	29
19	29
20	29
21	29
22	29
23	29
24	29
25	29
26	29
27	29
28	29
29	28
30	28
31	28
32	28
33	28
34	28
35	27
36	27
end
1	23.64
2	23.51
3	23.61
4	24.01
5	24.22
6	24.32
7	24.57
8	24.96
9	26.57
10	26.85
11	26.96
12	27.41
13	27.75
14	27.77
15	28.00
16	28.00
17	28.08
18	28.03
19	28.06
20	28.19
21	28.19
22	28.25
23	28.08
24	28.15
25	27.87
26	27.53
27	27.44
28	27.22
29	25.88
30	25.50
31	25.30
32	25.21
33	25.08
34	24.71
35	24.65
36	24.78
end
