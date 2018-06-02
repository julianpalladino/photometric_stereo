#!/bin/bash
 
for w in 1 5 10 15 20 30 50 100
do
 AVG=1 WINDOW=$w .././calibrate res/mate 2> /dev/null >> experimentacion/calibracion/luces_rgb_kbz/luces_mate_$w.txt
done
