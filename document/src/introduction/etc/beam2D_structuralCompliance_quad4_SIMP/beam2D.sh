#!/bin/sh

# Use -xv at the end of previous line to enter debug mode in terminal

SCRIPTPATH=$FEMDOCROOT/scripts/         # Store the path to the scripts folder
export PATH=$PATH:$SCRIPTPATH           # export the path to the scripts folder
chmod u+x $FEMDOCROOT/scripts/mkInput   # make sure the script is executible

echo ' ... Removing old files'
rm *.e-s.* *.data *.dat

mkdir tempDir

# Compile dynamically linked input files
echo
mkInput $FEMDOCROOT beam2D xdof opt
echo
PATH=${PATH/$SCRIPTPATH}
export PATH