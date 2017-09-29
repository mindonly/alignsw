
OBJS=

CC=clang
CFLAGS= -std=c11 -g -Wall -O2
CXX=clang++
CXXFLAGS= -std=c++14 -g -Wall -O2
LDLIBS=

EXE=align
RM=rm -f
ECHO=echo

$(P): $(OBJS)

all:	align t1 t2

t1:
	@echo '*** TEST EX1 ***'
	./$(EXE) ex1_seq.txt ex1_unk.txt
	@echo
	@sleep 2

t2:
	@echo '*** TEST EX2 ***'
	./$(EXE) ex2_seq.txt ex2_unk.txt
	@echo
	@sleep 2

t3:
	@echo '*** TEST HIV ***'
	./$(EXE) HIV-1_db.fasta HIV-1_Polymerase.txt
	@echo

t:	t1 t2 t3

clean:
	$(RM) $(EXE)
