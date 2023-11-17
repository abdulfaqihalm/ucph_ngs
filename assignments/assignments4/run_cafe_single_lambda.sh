times=5
for number in $(seq 1 $times)
do

        cafe5 -c 4 -i obp_all.hexapoda.tsv -t hexapoda.nwk -o cafe_out_one_lambda_r$number
done

