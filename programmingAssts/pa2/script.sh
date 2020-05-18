#!/bin/bash

echo lets narrow in!

for j in 8, 10, 12, 14, 16, 18
do
    for i in 137, 260, 277, 313, 380, 420, 463, 521, 557, 580, 620, 769, 780, 820, 887, 1123;
    do
        ./strassen $j $i text
    done
done

echo i am finished