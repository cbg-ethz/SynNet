#!/bin/sh
# This function gets one column in the file which includes Ensemble IDs with version and will remove the versioning from it.
# Pej Arp 2015
#--------------------
inDLM=','
outDLM='\t'

awk -F $inDLM -v Col="$2" -v OFS=$outDLM '{
                           n = match($Col,/\./) ;
                           if (n != 0) $Col=substr($Col,1,n-1);
                           else $Col=$Col  # If you exclude this the output delimiter for these lines remain the same as input!
                           print $0
                       }' $1
