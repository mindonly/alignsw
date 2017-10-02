
OBJS=

CC=clang
CFLAGS= -std=c11 -g -Wall -O2
CXX=clang++
CXXFLAGS= -std=c++14 -g -Wall -O2
LDLIBS= -pthread

EXE=
RM=rm -f
ECHO=echo

$(P): $(OBJS)

all:	clear $(P) t1 t2

t1:
	@echo '************'
	@echo '* TEST EX1 *'
	@echo '************'
	./$(EXE) ex1_seq.txt ex1_unk.txt
	@echo

t2:
	@echo '************'
	@echo '* TEST EX2 *'
	@echo '************'
	./$(EXE) ex2_seq.txt ex2_unk.txt
	@echo

t3:
	@echo '************'
	@echo '* TEST HIV *'
	@echo '************'
	./$(EXE) HIV-1_db.fasta HIV-1_Polymerase.txt
	@echo

t4:
	@echo '****************'
	@echo '* TEST HIV (2) *'
	@echo '****************'
	./$(EXE) HIV-1_db.fasta.wolffe HIV-1_Polymerase.txt.wolffe
	@echo

t:	t1 t2 t3

clean:
	$(RM) $(P)

clear:
	@clear
