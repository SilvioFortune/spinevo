#!/bin/bash

#INDIR="/e/ldata/users/sfortune/Dropbox/Apps/Overleaf/SpinEvolution/figures"
#INDIR="/home/moon/sfortune/spinevo/plots/thesis/"
INDIR="/e/ldata/users/sfortune/Dropbox/Apps/Overleaf/SpinEvolution/figures"
FILES="$INDIR/*.pdf"

echo 
echo 

echo 
for f in $FILES
do
    echo -e "\e[1A\e[K   Cropping $f"
    pdf-crop-margins -mo $f
done
echo -e "\e[1A\e[K Done Cropping!"
echo 

# clean-up

echo 
UNCROPPED="./*uncropped.pdf"
for u in $UNCROPPED
do
	echo -e "\e[1A\e[K   Deleting $u"
    rm $u
done
echo 

echo 