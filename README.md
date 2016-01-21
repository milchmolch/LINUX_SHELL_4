### University of Zurich
### URPP Evolution in Action
![URPP logo](Logo_URPP_kl2.png)

Stefan Wyder & Heidi Lischer

stefan.wyder@uzh.ch  
heidi.lischer@ieu.uzh.ch


## Bash scripting 4 - repetition and extension


The goal is to learn to write simple bash scripts in order to automate tasks  

Topic             |  
----------------- | --------------------------
Repetition on bash scripting: command structures (loops,if/else) | Heidi
Writing functions | Heidi
Writing clean code (10 min) | Stefan
Parallel jobs (10 min) | Stefan
Quick intro to cluster submission systems (10 min) | Heidi
Programming exercises | Heidi&Stefan



## Writing safe code in bash

Bash shell scripting is handy to automate tasks like file manipulation, program execution or printing text. It is often
used to set up a pipeline that executes a series of programs. However, bash syntax is often unintuitive
(e.g. string manipulation), difficult to read/remember and it does not provide libraries. 

If your script is longer than a few dozen or hundred lines of code then you should rather use a general-purpose programming language like Python or Perl.


Bash scripts are prone to errors. For [safer scripting](http://robertmuth.blogspot.ch/2012/08/better-bash-scripting-in-15-minutes.html) start every bash script with the following lines.

```
#!/bin/bash
set -o nounset
set -o errexit
```

This takes care of 2 very common errors
1. Referencing undefined variables (often due to typos)
2. Ignoring failing commands


### Debugging

To perform a syntax check/dry run of your bash script, run:
```
bash -n script.sh
```

If something does not work as expected, you can trace a script using the `-x` option, like so: 
```
bash -x script.sh
```
This will print each command (predeceded by a "+") before it is executed. Also try `-v` option (=`bash --verbose script.sh`).  
  

## Parallel jobs

[GNU Parallel](https://www.gnu.org/software/parallel/) makes analysis faster by parallelizing jobs.

- can handle any number of jobs
- even on remote computers


### Installation

Ubuntu

```
sudo apt-get install parallel
```

Mac using Homebrew

```
brew install parallel
```

### Examples

Run 10 FASTQC jobs in parallel 
```
ls *.fq | parallel  -j 12 "fastqc {} --outdir ."
```

Index your BAM files in parallel. Prints only the commands, doesn't run them (--dry-run)
```
ls *.bam | parallel --dry-run 'samtools index {}'
```

We can actually drop the ticks:
```
ls *.bam | parallel samtools index {}
```

There is an alternative way to run something on all BAM files:
```
parallel samtools index ::: *.bam
```

Running Freebayes in parallel - 12 chromosomes, each in a separate job
```
command="freebayes --ploidy 2 -f GENOME.fa BAMfile.bam"
seq 1 12 | awk '{print "chr"$1}' | parallel --keep-order -j 10 "$command -r {}" 
```

Running your own script in parallel (`--load 90` max load 90% on the computer)
```
nohup ls SINGLE_END/*.fastq.gz | parallel -j20 --load 90 ./Run_kallisto_SingleEnd.sh 2>log.kallisto & 
```

```
#!/bin/sh
# Runs kallisto on all Fastq Files in Single_End with arg $1
# for use with GNU parallel 
  
KALLISTO=~/APPL/KALLISTO/kallisto_linux-v0.42.1/kallisto
INDEX_FILE=KALLISTO_OUT/Athaliana_longestTranscript.idx
OUT_DIR=KALLISTO_OUT
  
FASTQ=$1
  
# Get sample name
filename=${FASTQ%%.fastq.gz}
filename=$(basename $filename)
  
# run kallisto quantification
$KALLISTO quant -i $INDEX_FILE --single -l 180 $FASTQ -o TEMP_${filename}
```


## Customizing

.bashrc
.bash_profile
```
alias grepcol="grep --color"
```

Try `grepcol` 

```
alias tarup="tar -zcf"
alias tardown="tar -zxf"
```

## Exercises

### Repetition

1. Explore the various ways to go wrong when you declare and use variables
```
a = "URPP Evolution"
a= "URPP Evolution"
a="URPP Evolution"
echo a
echo $a
```

Variables must start with a letter, most not contain spaces or punctuation marks  
Valid examples: BaseDir, my_project_dir  


2. Working with Variables
```
# Put the output of a command into a variable
a=$(whoami)
echo $a
free=$(df -h . | tail -1 | tr -s " " " " | cut -d" " -f4)
echo $free

a=4
b=$(( a+3 ))
echo $b

a="foo.bar"
echo ${a%%.bar}
echo ${a##foo.}
```


3. Functions

Use a function, if you use a code block more than once. 
```
ExtractBashComments() {
    egrep "^#"
} 
```
defines a function `ExtractBashComments`. It extracts comment lines from a script:  `cat test_safe_script.sh | ExtractBashComments`   

We can also define a function without an argument (the keyword `function` is optional, we put it for better readibility):

```
function today {
	echo "Today's date is: "
	date +"%A, %B %-d, %Y"
}
```
which you can run by `today`  


Another example with local variables:
```
SumLines() {  # Iterating over stdin 
	local sum=0
	local line=""
	while read line; do
		sum=$((${sum} + ${line}))
	done
	echo $sum
}
cat $1 | SumLines
```

We can now execute the script `bash SumLines.sh numbers.txt`


4. Write a loop that prints out chromosomes chr5 - chr9
```
for i in 5 6 7 8 9
do
	echo "chr"$i 
done
```

or alternatively (see also <nested commands> above):
```
for i in `seq 5 9`
do
        echo "chr"$i
done
``

5. Write a script that reads from a file


6. Modify the following script to make it safe 
```
#!/bin/bash

cutoff=0.05

echo $cutoff
# Referencing undefined variables (which default to "") 
echo $Cutoff
echo $unset_variable

# failing commands are ignored
cd NON_EXISTING_FOLDER
echo "last line"
```

7. Write a bash script that asks the user to enter a number 1-3 and prints out a text (Tip: use the `case` command) 

```
#!/bin/bash

echo -n "Enter a number between 1-3 > "
read choice
case $choice in
        1) echo "You selected 1."
           # some commands
           ;;
        2) echo "You selected 2."
           # some commands
           ;;   
        3) echo "You selected 3."
           # some commands
           ;;   
        *) echo "You did not enter a number between 1-3"
           ;;
esac
```



## Advanced topics

### Nested commands

Use backticks ` to execute a command inside another

```
for i in `ls`; do
        echo $i
done
```

One can also pipe multiple commands together
```
for i in `ls | grep foo | tr a-z A-Z` do
        echo $i
done
```

## Sources / Links

https://portal.tacc.utexas.edu/documents/13601/1080823/Shell+scripting+2014+eijkhout+%281%29.pdf/52353fc6-0fff-4efb-bfda-4a511827332e

https://portal.tacc.utexas.edu/documents/13601/1080823/LinuxIntro-20141009-eijkhout+%281%29.pdf/bcdcefad-47c5-4741-ab9f-c3380e63df93

http://explainshell.com/

https://github.com/swcarpentry/good-enough-practices-in-scientific-computing/blob/gh-pages/index.md

http://www.gnu.org/software/coreutils/manual/coreutils.html


- **GNU parallel**  
  http://www.gnu.org/software/parallel
  https://www.biostars.org/p/63816/

- **Safe Bash scripting***
  http://robertmuth.blogspot.ch/2012/08/better-bash-scripting-in-15-minutes.html

- **Tips & Tricks for using the shell on Mac OS**
  http://furbo.org/2014/09/03/the-terminal/

