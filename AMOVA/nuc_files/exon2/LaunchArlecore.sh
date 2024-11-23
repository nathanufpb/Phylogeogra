#!/bin/bash

#Laurent Excoffier November 2009
#
#The script will launch arlecore on all arlequin project files in turn
#It assumes that it is launched in a directory containing:
#         - a series of *.arp files to be analysed
#         - a settings file containing the settings specifying which 
#           computations are to be performed (usually obtained through the WinArl35.exe
#           graphical interface).

#Modify the following line to state which version of arlsumstat you are using
arlecore=arlecore3512_64bit
#Modify this if you wan to use another setting file
settingsFile=cmyc_pahsed_exon2.arp

if [ -f $settingsFile ]; then

	#This loop will analyse all files with the same settings file
	for file in *.arp
	do  	
		echo "Processing file $file"
		#Launch arlecore with the same settings for all files
		./$arlecore  $file $settingsFile	run_silent
	done

	#The following loop would be an alternative to perform specific computations on each file
	#assuming that there is a different setting file with the same name as the project file, but with 
	#the extension *.ars 

	#for file in *.arp
	#do  	
	#	echo "Processing file $file"
	#	#Launch arlecore with project specific settings
	#	./$arlecore  $file ${file%.*}.ars run_silent
	#done
else 
	echo "Settings file $settingsFile does not exist. Aborting script."
fi
