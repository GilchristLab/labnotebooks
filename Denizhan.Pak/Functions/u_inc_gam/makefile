LINKS = -lm -lgsl -lgslcblas
INCLUDES = /usr/local/include
COMPILER = gcc
FLAGS = -I $(INCLUDES)

accurracy_test: accurracy_test.o mygamma.o mygamma.h
	$(COMPILER) -o accurracy_test accurracy_test.o mygamma.o $(FLAGS) $(LINKS)
accurracy_test.o: accurracy_test.c mygamma.h
	$(COMPILER) -c accurracy_test.c
mygamma.o: mygamma.c mygamma.h
	$(COMPILER) -c mygamma.c -lm
clean:
	rm -f accurracy_test accurracy_test.o mygamma.o
clear:
	rm -f Log*
c:
	rm -f Log*
	rm -f accurracy_test accurracy_test.o mygamma.o
