#!/bin/bash
cd ~/analysis_oindree/

h_end=261;
end_events=162000;
delta=1000;
for((phase_number=3;phase_number<4;phase_number++))
do
    if((phase_number==1 ))
    then phase_number=3
    fi
   
    for((filter_number=3;filter_number<4;filter_number++))
    do 
       #if((filter_number==4 ))
       #then filter_number=5
       #fi

       #if((filter_number==1 ))
       #then filter_number=3
       #fi

       if((filter_number==3 ))
       then phase_number=3
       fi
       
	h=155;
	while (( h < h_end ))
	do
	    qstat -u dailey.110 > output.txt
	    j=$(perl getnumber.pl)
	    
	    for((k=j; k<40; k+=40))
	    do
	    
		for((i=0; i<end_events; i+=delta))
		do
		
		FILENAMER="Run$h_$i";
		let "end_number=i+delta-1";
		let "start_number=i";
		qsub run90sample_pieces.sh -v RUNNUMBER=$h,FILENAME=$FILENAMER,FILTER=$filter_number,PHASE=$phase_number,SIMFILE=$start_number,start_number=$start_number,end_number=$end_number;
		sleep 2
		done
		
		let "h=h+1"  
		if((h >= h_end))
		
		then break
		fi
	    done
	    sleep 300
	done
	
#	if((filter_number==3 ))
#	then phase_number=0
#	fi
    done
done 