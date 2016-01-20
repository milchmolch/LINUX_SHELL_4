### University of Zurich
### URPP Evolution in Action
![URPP logo](Logo_URPP_kl2.png)

Stefan Wyder & Heidi Lischer

stefan.wyder@uzh.ch  
heidi.lischer@ieu.uzh.ch


## Bash scripting 3 - repetition and extension


The goal is to improve/learn to write simple bash scripts in order to automate tasks:

- repetition on bash scripting: command structures (loops, if else), writing functions (Heidi)
- programming exercises to write simple scripts (traversing through files, read in / rename / merge multiple files, run software on all files in a folder, etc)
- writing clean code (10 min) Stefan
- parallel jobs (10 min) Stefan
- quick intro to cluster submission system (10 min) Heidi





## Writing safe code in bash

Bash shell scripting is handy to automate tasks or to put together a pipeline that launches a series
of programs.

If your script is longer than a few hundred lines of code then rather use a scripting language like Python or Perl.


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

## Links

http://robertmuth.blogspot.ch/2012/08/better-bash-scripting-in-15-minutes.html

# Customizing

.bashrc
.bash_profile

alias grepcol="grep --color"
grepcol 

## Exercises

### Repetition

- 1 - Explore the various ways to go wrong when you declare and use variables
``
a = "URPP Evolution"
a= "URPP Evolution"
a="URPP Evolution"
echo a
echo $a
``

Variables must start with a letter, most not contain spaces or punctuation marks
Valid examples: BaseDir, my_project_dir


- 2 - Working with Variables
``
a=$(whoami)
echo $a

a=4
b=$(( a+3 ))
echo $b

a="foo.bar"
echo ${a%%.bar}
echo ${a##foo.}
``


function today {
	echo "Today's date is: "
	date +"%A, %B %-d, %Y"
}




## Sources / Links

https://portal.tacc.utexas.edu/documents/13601/1080823/Shell+scripting+2014+eijkhout+%281%29.pdf/52353fc6-0fff-4efb-bfda-4a511827332e

https://portal.tacc.utexas.edu/documents/13601/1080823/LinuxIntro-20141009-eijkhout+%281%29.pdf/bcdcefad-47c5-4741-ab9f-c3380e63df93

http://explainshell.com/


