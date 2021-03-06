#!/bin/bash
cd ~/analysis_oindree/

h_end=262
for((phase_number=3;phase_number<4;phase_number++))
do
    if((phase_number==1 ))
    then phase_number=3
    fi
   
    for((filter_number=3;filter_number<4;filter_number++))
    do 
       if((filter_number==4 ))
       then filter_number=5
       fi

       if((filter_number==2 ))
       then filter_number=3
       fi
	h=12
	while (( h < h_end ))
	do
	    qstat -u dailey.110 > output.txt
	    j=$(perl getnumber.pl)
	    
	    for((i=j; i<40; i+=1))
	    do
	    
		FILENAMER="Run$h";
		qsub runfilterpulser.sh -v RUNNUMBER=$h,FILENAME=$FILENAMER,FILTER=$filter_number,PHASE=$phase_number;
		let "h=h+1"  
		if((h >= h_end))
		then break
		fi
		
	   
	    done
	
	    sleep 20
	
	done

    done
done 