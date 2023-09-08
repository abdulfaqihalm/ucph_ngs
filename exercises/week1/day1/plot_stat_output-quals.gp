
            set terminal png size 600,400 truecolor
            set output "plot_stat_output-quals.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set ylabel "Average Quality"
            set xlabel "Cycle"
            set yrange [0:50]
            set title "output_aln.stat" noenhanced
            plot '-' using 1:2 with lines title 'Forward reads' 
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
