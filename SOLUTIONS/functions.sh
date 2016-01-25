#!/bin/bash

ExtractBashComments() {
    egrep "^#"
} 

cat test_safe_script.sh | ExtractBashComments | wc -l
