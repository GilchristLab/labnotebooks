LINKS = -lm -lgsl -lgslcblas
INCLUDES = /usr/local/include
COMPILER = gcc
FLAGS = -I $(INCLUDES)

test: test.o mygamma.o mygamma.h
	$(COMPILER) -o test test.o mygamma.o $(FLAGS) $(LINKS)
test.o: test.c mygamma.h
	$(COMPILER) -c test.c
mygamma.o: mygamma.c mygamma.h
	$(COMPILER) -c mygamma.c -lm
clean:
	rm -f test test.o mygamma.o
clear:
	rm -f Log*
c:
	rm -f Log*
	rm -f test test.o mygamma.o
