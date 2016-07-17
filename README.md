# MixTaR
MixTaR: software for de novo tandem repeat detection

Compilation
pour version gatb-core <= 1.1.0
g++  main.cpp -I../../../gatb-core-1.1.0-Linux/include -L../../../gatb-core-1.1.0-Linux/lib -lgatbcore -lhdf5 -ldl -lz -lpthread -fopenmp -O3 -o main.o

pour version gatb-core > 1.1.0
g++  main.cpp -I../../../gatb-core-1.2.1-bin-Linux/include -L../../../gatb-core-1.2.1-bin-Linux/lib -std=c++11 -lgatbcore -lhdf5 -ldl -lz -lpthread -fopenmp -O3 -o main.o

Executions instances de test
./main.o -k 17 -t 10 -s ../tests/pairedShortReads_1.fasta ../tests/pairedShortReads_2.fasta -c 20 -l ../tests/longReads.fasta -d 20 -r 2 -R 100 > testMini.txt 
