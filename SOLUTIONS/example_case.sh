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
