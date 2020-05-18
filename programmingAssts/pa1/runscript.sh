#!/bin/bash

echo hello, Mitz!

# for j in 128 256 512 1024 2048 4096 8192 
# do
#     for i in 1 2 3 4;
#     do
#         ./randmst 0 $j 5 $i
#     done
# done

# for j in 16384 32768 65536
# do
#     for i in 1 2 3 4;
#     do
#         ./randmst 0 $j 5 $i
#     done
# done


# for i in 1 2 3 4;
# do
#     ./randmst 0 131072 3 $i
# done


for i in 2 2;
do
    ./randmst 0 262144 1 $i
done