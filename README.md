### University of Zurich
### URPP Evolution in Action
![URPP logo](Logo_URPP_kl2.png)

Stefan Wyder & Heidi Lischer

stefan.wyder@uzh.ch  
heidi.lischer@ieu.uzh.ch


## Bash scripting 4 - repetition and extension


The goal is to learn to write simple bash scripts in order to automate tasks  

Topic             | Person 
----------------- | --------------------------
Repetition on bash scripting: command structures (loops,if/else) | Heidi
Writing functions | Heidi
Writing safe code (10 min) | Stefan
Parallel jobs (10 min) | Stefan
Quick intro to cluster submission systems (10 min) | Heidi
Exercises | Heidi&Stefan
  
  
Heidi's parts: [Theory](URPP_Tutorial_BashScripting_2_HL.pdf) | [Exercises](Exercises_BashScripting_2_HL.pdf)


## Writing safe code in bash

Bash shell scripting is handy to automate tasks like file manipulation, program execution or printing text. It is often
used to set up a pipeline that executes a series of programs. However, bash syntax is often unintuitive
(e.g. string manipulation), difficult to read/remember and it does not provide libraries. 

If your script is longer than a few dozen or hundred lines of code then you should rather use a general-purpose programming language like Python or Perl.


Bash scripts are quite prone to errors. For [safer scripting](http://robertmuth.blogspot.ch/2012/08/better-bash-scripting-in-15-minutes.html) start every bash script with the following lines:

```
#!/bin/bash
set -o nounset
set -o errexit
```

This takes care of 2 very common errors:    
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
This will print each command (predeceded by a "+") before it is executed. Also try the `-v` option (=`bash --verbose script.sh`).  
  

## Parallel jobs

[GNU Parallel](https://www.gnu.org/software/parallel/) makes analysis faster by parallelizing jobs, i.e. by running multiple jobs simultaneously.

- can handle any number of jobs
- by default uses all CPUs on your machine (limit using the `-j` option)
- even on remote computers


### Installation

On Ubuntu:

```
sudo apt-get install parallel
```

On Mac using Homebrew:

```
brew install parallel
```

### Examples

Run 10 FASTQC jobs in parallel 
```
ls *.fq | parallel  -j 12 "fastqc {} --outdir ."
```

Index your BAM files in parallel. Here we only prints the commands, we don't run them (--dry-run)
```
ls *.bam | parallel --dry-run 'samtools index {}'
```

We can actually drop the ticks:
```
ls *.bam | parallel --dry-run samtools index {}
```

There is an alternative way to run something on all BAM files in the current directory:
```
parallel samtools index ::: *.bam
```

Running Freebayes in parallel - 12 chromosomes, each in a separate job (`--keep-order` keeps the output in input order):
```
command="freebayes --ploidy 2 -f GENOME.fa BAMfile.bam"
seq 1 12 | awk '{print "chr"$1}' | parallel --keep-order -j 10 "$command -r {}" 
```

Running your own script in parallel (`--load 90` max load 90% on the computer):
```
ls SINGLE_END/*.fastq.gz | parallel -j20 --load 90 ./Run_kallisto_SingleEnd.sh 2>log.kallisto
```

```
more Run_kallisto_SingleEnd.sh
#!/bin/bash
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

This will the script `Run_kallisto_Single` using 20 CPUs on all FASTQ files of the SINGLE_END subdirectory. Error messages will be
written to the file `log.kallisto`.


## (optional) Customizing the shell

### Define aliases

We can define an alias, i.e. a sort of shortcut which is easier to remember and/or saves some typing.

```
alias grepcol="grep --color"
```
defines a new command `grepcol` that color-marks the matching text. Try e.g. `ls | grepcol "\.sh"` 
  
The following aliases define aliases to create or extract tar archives.

```
alias tarup="tar -zcf"
alias tardown="tar -zxf"
```

An alias is only active in the current window. To make it permanent we have to add it to the hidden configuration files
`.bashrc` and `.bash_profile` located in your home directory.  
  
The two files have different roles. The `~/.bashrc` file is a script executed whenever a new terminal session is started in **interactive
mode** (if open a new terminal window or tab). In contrast, if you log in to session (asking for a password) the script `~/.bash_profile`
is executed (e.g. if you login via ssh).    
  
Ubuntu comes with well commented `.bashrc` and `.profile` files that already define some useful aliases. We modify `~/.bashrc` for:  
- Creating useful aliases (for example alias ll='ls -l').
- Setting new environment variables.

For adding more directories to PATH we add to `.profile` the following line:   
```
export PATH="/path/to/dir:$PATH"
```

~/.bashrc is actually only executed by the bash shell, other shells work differently.

On Mac OS the configuration also look different, they are called `~./bash_profile` and `~/.profile` 
  
Check both files in your home directory.


## Exercises

### Repetition - Iron Ration

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

a="test.bam"
echo ${a%.bam}   # shortest match from back - removes file extension
echo ${a%.*}     # removes any file extension
echo ${a#*.}     # shortest match from front
echo ${a:1:2}    # get substring
echo ${a/test/test2}	 # replace: ${string/pattern/replacement}

filename="Sample3.BWA.bam"
echo ${filename##*.}    # longest match from front - prints file extension even with multiple "."
echo ${filename%%.*}    # longest match from back - prints sample name
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
		sum=$(($sum + $line))
	done
	echo $sum
}
cat $1 | SumLines
```

  We can now execute the script by doing `bash SumLines.sh numbers.txt`


4. Write a script that prints out chromosomes chr5 - chr9  
  
  
5. Modify the following script to make it safe  
  
  
6. Write a bash script that asks the user to enter a number 1-3 and prints out a text (Tip: use the `case` command) 



## Advanced topics

### Nested commands

Use backticks ` to execute a command inside another

```
for i in `ls`
do
        echo $i
done
```

One can also pipe multiple commands together
```
for i in `ls | grep foo | tr a-z A-Z`
do
        echo $i
done
```

## Solutions

4. Write a script that prints out chromosomes chr5 - chr9

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
```

5. Modify the following script to make it safe 

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

6. Write a bash script that asks the user to enter a number 1-3 and prints out a text (Tip: use the `case` command) 

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

## Sources / Links

- Explain shell commands http://explainshell.com/

- List of commands  
  http://www.gnu.org/software/coreutils/manual/coreutils.html

- **GNU parallel**  
  http://www.gnu.org/software/parallel
  https://www.biostars.org/p/63816/

- **Safe Bash scripting***  
  http://robertmuth.blogspot.ch/2012/08/better-bash-scripting-in-15-minutes.html

- **Tips & Tricks for using the shell on Mac OS**  
  http://furbo.org/2014/09/03/the-terminal/

- Some examples are from  
  https://portal.tacc.utexas.edu/documents/13601/1080823/
