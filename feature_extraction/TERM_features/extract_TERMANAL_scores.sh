#!/bin/bash
# for each design, take the abundance, design and structural score files that contain all the per-position information and extract it to a CSV.
# Reffer to "Tertiary structural propensities reveal fundamental sequence/structure relationships" on how to obtain these values and their meaning.
echo ",abd_w,abd_b,abd_av,dsc_w,dsc_b,dsc_av,ssc_b,ssc_w,ssc_av"
for i in $(ls termanal_digs/0????/*abd50 | cut -d'.' -f1);
	do echo\
	$(basename $i),\
	$(cat $i.abd50 | sort -k2 -n | awk '{print $2 }' | head -n1 ),\
	$(cat $i.abd50 | sort -k2 -n | awk '{print $2 }' | tail -n1),\
	$(cat $i.abd50 | awk '{print $2 }' | awk '{sum+=$1} END { print sum/NR}' ),\
	$(cat $i.dsc50 | sort -k2 -n | awk '{print $2 }' | head -n1 ),\
        $(cat $i.dsc50 | sort -k2 -n | awk '{print $2 }' | tail -n1),\
        $(cat $i.dsc50 | awk '{print $2 }' | awk '{sum+=$1} END { print sum/NR}' ),\
	$(cat $i.ssc50 | sort -k2 -n | awk '{print $2 }' | head -n1 ),\
        $(cat $i.ssc50 | sort -k2 -n | awk '{print $2 }' | tail -n1),\
        $(cat $i.ssc50 | awk '{print $2 }' | awk '{sum+=$1} END { print sum/NR}' );
done
