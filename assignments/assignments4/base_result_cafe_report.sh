read -p "Enter prefix cafe output: " prefix
times=5
for number in $(seq 1 $times)
do
        head -2 ${prefix}${number}/Base_results.txt
done




