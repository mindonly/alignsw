
OBJS=

CC=gcc
CFLAGS= -std=c11 -g -Wall
CXX=g++
CXXFLAGS= -std=c++11 -g -Wall
LDLIBS=

EXE=align
RM=rm -f
ECHO=echo

$(P): $(OBJS)

all:	align test3
test1:
	@echo '*** TEST EX1 ***'
	./$(EXE) ex1_seq.txt ex1_unk.txt
	@echo
test2:
	@echo '*** TEST EX2 ***'
	./$(EXE) ex2_seq.txt ex2_unk.txt
	@echo
test3:
	@echo '*** TEST HIV ***'
	./$(EXE) HIV-1_db.fasta HIV-1_Polymerase.txt
	@echo
test:	test1
clean:
	$(RM) $(EXE)
