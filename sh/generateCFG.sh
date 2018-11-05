#!/bin/bash

#script to generate config files

if [ $# -ne 1 ]
then
	echo "Usage $0 sourcePathRelativeToCNFfolder"
	exit 1
fi

path=$1
input="../cnf/"$path
path=`echo $path | tr "/" .`
output="../cfg/"$path".txt"

> $output
for f in `ls -1 $input | grep "cnf"`
do
    sed -i -e '/^c/ d' $input"/"$f
	echo ${f%.cnf} >> $output
    rm $input"/"${f%.cnf}.cnf-e
done
