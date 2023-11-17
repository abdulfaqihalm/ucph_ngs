times=5
for number in $(seq 1 $times)
do

        cafe5 -c 4 -i obp_all.hexapoda.tsv -t hexapoda.nwk -y hexapoda_2lambdamodel.nwk -o cafe_out_2l_r$number
done

