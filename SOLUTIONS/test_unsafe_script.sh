#!/bin/bash

cutoff=0.05

echo $cutoff
# Referencing undefined variables (which default to "") 
echo $Cutoff
echo $unset_variable

# failing commands are ignored
cd NON_EXISTING_FOLDER
echo "last line"
