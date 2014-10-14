#!/bin/bash

# release job from hold using qrls (q release)

# $1 job id to be removed from hold 
# $2 output file

cat blastout* > $2

echo "release job with jid = "$1 

qrls $1
