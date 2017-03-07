#!/bin/bash

lt1=longtables/1.canu.corcov500minrl500.aspb-e02.pball.ont2d.quiverfinal1.pilon2x
lt2=longtables/4.falcon.seed25.relaxed.pball.ontmol.quiverfinal1.pilon2x
names=../names.txt

bash ../convertLongTableToLatex.sh $names $lt1 $lt2
