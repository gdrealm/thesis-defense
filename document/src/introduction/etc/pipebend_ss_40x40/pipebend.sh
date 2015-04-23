#!/bin/sh

# Use -xv at the end of previous line to enter debug mode in terminal

SCRIPTPATH=$FEMDOCROOT/scripts/         # Store the path to the scripts folder
export PATH=$PATH:$SCRIPTPATH           # export the path to the scripts folder
chmod u+x $FEMDOCROOT/scripts/mkInput   # make sure the script is executible

rm pipebend.so
rm pipebend.o
rm pipebend.exo
rm pipebend.e-s.*

mkdir tempDir

# Compile dynamically linked input files
echo mkInput $FEMDOCROOT pipebend
mkInput $FEMDOCROOT pipebend
echo
PATH=${PATH/$SCRIPTPATH}
export PATH
