#!/bin/bash

cat score_tjfrags_* > fragment_met_folder.tab
for i in $(awk '{print $1}' fragment_met_folder.tab); do TESTPATH=$(ls -ltr $i | awk '{print $11}'); grep $i fragment_met_folder.tab | awk -v dt="$TESTPATH" 'BEGIN{FS=OFS=" "}{$1=dt}1' ; done > fragment_met_path.tab
