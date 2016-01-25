# Use backticks ` to execute a command inside another
for i in `ls`
do
	echo $i
done

for i in `ls | grep foo | tr a-z A-Z`
do
	echo $i
done
