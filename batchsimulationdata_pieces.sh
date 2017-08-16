#!/bin/bash
cd ~/analysis_oindree/

h_end=2;
end_events=331000;
delta=3000;
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
       
	h=1;
	    for((i=0; i<end_events; i+=delta))
	    do
		
		FILENAMER="Run$h";
		let "end_number=i+delta-1";
		let "start_number=i";
		qsub runsimulationdata_pieces.sh -v RUNNUMBER=$h,FILENAME=$FILENAMER,FILTER=$filter_number,PHASE=$phase_number,SIMFILE=$start_number,start_number=$start_number,end_number=$end_number;

	    done
	
	    #sleep 20
	
	
#	if((filter_number==3 ))
#	then phase_number=0
#	fi
    done
done 