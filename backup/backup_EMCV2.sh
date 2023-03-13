#!/bin/bash

FOLDER=$1

mkdir $FOLDER
mkdir $FOLDER/EMC

## copies all files below 10Mb keeping folder structure, however not BACKUP folder (otherwise infinite recursion)
##find "/home/emil" -size -10M -exec cp -r --parents {} /emc/emil/BACKUP/${FOLDER} \;

echo "BE CAREFUL WILL COPY EVERYTHING IN TO BACKUP FOLDER INDEFINETELY"

#########################################################################################
## DOES INFINITE RECURSIONS WHERE CP BACKUP FOLDER INTO BACKUP FOLDER OF BACKUPED
## SO YOU WILL EVERYTHING YOU BACKED UP IN BACKUP FOLDER AND SO ON...
## THEREFORE WATCH THIS WHILE RUNNING AND STOP!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
###################################################################################


find /emc/emil/ -type f ! -path '/emc/emil/BACKUP/*' -size -10M -type f -exec cp -p --parents -t /emc/emil/BACKUP/${FOLDER}/EMC/ {} +


