#for i in {1..5}
#do
#random_number=$WEST_RAND16
#echo $random_number
#done
random_number=$RANDOM
random_number1=$(echo $random_number + 1 | bc) 
echo $random_number
echo $random_number1
