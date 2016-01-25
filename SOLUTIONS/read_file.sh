while read col1 col2
do
        sum=$(($sum + $col1 ))
done < input.txt
echo $sum
