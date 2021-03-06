# flags:s
# FL: full projection matrices (V * V)

.SUFFIXES: .c .u
CC= gcc

####### TI sampler, IP ###########
CFLAGS = -Wall -I/home/hardyman/Documents/ebuka/gsl/include -L/home/hardyman/Documents/ebuka/gsl/lib -O3 -DFL -DIP -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF # sample only b

########################################################################################

LDFLAGS = -lm -lgsl -lgslcblas
#LDFLAGS = -lm -lgsl -lcblas
LOBJECTS= seqTM.o estimate.o gibbs.o est_gibbs.o inf_gibbs.o util.o corpus.o sort.o params.o stable.o stirling.o model.o sym_polya_fit.o arms.o sample_b.o sample_alpha.o hash.o
LSOURCE= seqTM.c estimate.c gibbs.c est_gibbs.c inf_gibbs.c util.c corpus.c sort.c params.c stable.c stirling.c model.c sym_polya_fit.c arms.c sample_b.c sample_alpha.c hash.c

all:	$(LOBJECTS)
	$(CC) $(CFLAGS) $(LOBJECTS) -o spyp $(LDFLAGS)

debug:	$(LOBJECTS)
	$(CC) $(CFLAGS) $(LOBJECTS) -o spyp $(LDFLAGS)

clean:
	-rm -f *.o
	-rm -f spyp

build:
	make all

TRSETTING = ./tr_setting.txt

##############################################################

DATASET = test  # folder name of dataset
VDATA = ../SPYP_Exp/data/$(DATASET)/elcap2.vo # vocabulary file

TDATA = ../SPYP_Exp/data/$(DATASET) # data file path


OFILE = ../SPYP_Exp/results/$(DATASET) # output folder


RESULT_DIR = ../SPYP_Exp/results/ # result directory

#############################################################
est:
	./spyp -e -n $(RESULT_DIR) -z 0 -f NULL -v $(VDATA) -t $(TDATA) -o $(OFILE) -s $(TRSETTING) #>$(OFILE)/test_$(RESULT_DIR).log
