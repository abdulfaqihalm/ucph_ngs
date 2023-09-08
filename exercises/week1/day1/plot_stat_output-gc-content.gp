
            set terminal png size 600,400 truecolor
            set output "plot_stat_output-gc-content.png"
            set grid xtics ytics y2tics back lc rgb "#cccccc"
            set title "output_aln.stat" noenhanced
            set ylabel "Normalized Frequency"
            set xlabel "GC Content [%]"
            set yrange [0:1.1]
            set label sprintf("%.1f",39.95) at 39.95,1 front offset 1,0
            plot '-' smooth csplines with lines lc 1 title 'First fragments' 
        1	0.011433
4	0.007842
6	0.014417
9	0.024809
12	0.045263
15	0.080463
17	0.139695
20	0.224988
23	0.341128
26	0.477860
28	0.628020
31	0.775917
34	0.899630
37	0.981668
39	1.000000
42	0.928828
45	0.810910
48	0.683943
51	0.568194
54	0.471421
56	0.387053
59	0.307188
62	0.227196
65	0.153617
67	0.099623
70	0.061092
73	0.034338
76	0.018763
78	0.009134
81	0.004442
84	0.002177
87	0.001359
89	0.001110
92	0.000727
95	0.000377
98	0.000162
end
