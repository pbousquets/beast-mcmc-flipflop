#!/bin/bash

usage="$0 folder file"

if [ $# -ne 2 ] || [ ! -f "$2" ] || [ ! -d "$1" ] || [ $1 == "-h" ] || [ $1 == "--help" ]
then
    echo -e $usage
    exit 1
fi

folder=$1
seed=$folder
file=$(readlink -e $2)

cd $folder

~/Downloads/FlipFlopBEASTv0.1/bin/beast -beagle_off -seed $seed $file 
~/Downloads/FlipFlopBEASTv0.1/bin/treeannotator -heights median -burnin 75000 $(basename $file | sed "s/.xml/.trees/") $(basename $file | sed "s/.xml/.tree/")
