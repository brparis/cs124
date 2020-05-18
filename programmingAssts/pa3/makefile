# dual dependency of parition.cpp
all: partition.cpp maxheap.cpp maxheap.h
	g++ -c maxheap.cpp -o maxheap.o
	g++ -o partition partition.cpp maxheap.o

parition: partition.cpp, maxheap.o
	 g++ -o partition partition.cpp maxheap.o
     
# object file to link
maxheap: maxheap.cpp maxheap.h
	g++ -c maxheap.cpp -o maxheap.o

# removes the executable
clean:
	 rm partition
	 rm maxheap.o