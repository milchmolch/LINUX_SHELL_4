SumLines() {  # Iterating over stdin 
        local sum=0
        local line=""
        while read line; do
                sum=$((${sum} + ${line}))
        done
        echo $sum 
}
cat $1 | SumLines
