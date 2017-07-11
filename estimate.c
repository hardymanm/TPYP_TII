/*
 * seqTM.c
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_math.h>
#include <libgen.h>

#include "seqTM.h"
#include "corpus.h"
#include "params.h"

extern dpyp_params PARAMS;

struct Args
{
	int do_ei;
	char* run_num;
	int heldout;
	int I;
	int K;
	int z;
	int seed;

	char* z_root;
	char* setting;
	char* data_name;
	char* v_data;
	char* mod_file;
	char* o_file;
	char* pmifile;
}Args;

/*
 * Command line processing
 */
void initialize_args()
{
	Args.do_ei = 0;
	Args.run_num = NULL; //1;
	Args.heldout = 1;
	Args.I = 0;
	Args.K = 0;
	Args.z = 0;
	Args.seed = -1;

	Args.setting = NULL;
	Args.data_name = NULL;
	Args.v_data = NULL;
	Args.mod_file = NULL;
	Args.o_file = NULL;
	Args.pmifile = NULL;
}

static const char* optString = "eipus:n:v:V:S:t:T:I:o:m:P:a:b:B:A:G:q:z:f:k:h?";

int main(int argc, char* argv[])
{
	int opt = 0;
	//gsl_vector* aa = gsl_vector_calloc()
	initialize_args();
	opt = getopt(argc, argv, optString);
	while(opt != -1)
	{
		switch(opt){
			case 'e': // for running estimation
				Args.do_ei = 1;
				break;
			case 'i': // for running inference
				Args.do_ei = -1;
				break;
			case 'p':
				Args.do_ei = 0;
				break;
			case 'u': // using held-out/left-to-right likelihood calculation
				Args.heldout = 0;
				break;
			case 'n':
				Args.run_num = optarg; //atoi(optarg);
				break;
			case 'I':
				Args.I = atoi(optarg);
				break;
			////////////////////////////////////////////
			case 'v': // for vocabulary file
				Args.v_data = optarg;
				break;
			case 't': // for training or testing data
				Args.data_name = optarg;
				break;
			case 's': // for training/testing setting file
				Args.setting = optarg;
				break;
			case 'S':
				Args.seed = atoi(optarg);
				break;
			//////////////////////////////////////////
			case 'o':
				Args.o_file = optarg;
				break;
			case 'm':
				Args.mod_file = optarg;
				break;
			case 'P':
				Args.pmifile = optarg;
				break;
			case 'z':
				Args.z = atoi(optarg);
				break;
			case 'f':
				Args.z_root = optarg;
				break;
			case 'h': // for command help
				printf("help!!!\n");
				break;
			default:
				break;
		}
		opt = getopt(argc, argv, optString);
	}

	vocabulary* v;
	Corpus *c_tr, *c_te;
	char str[BUFSIZ];
	char root[BUFSIZ];

	printf("COMMAND-LINE: %s ", basename(argv[0]));
	int loop;
	for (loop=1; loop < argc; loop++){
		printf("%s ", argv[loop]);
	}
	printf("\n\n");

	if(Args.do_ei == 1){
		printf("Sampling the Shadow PYP model ...\n");
		v = read_vocabulary(Args.v_data);
		sprintf(root, "%s", Args.o_file);
		mkdir(root, S_IRUSR | S_IWUSR | S_IXUSR);
		read_params(Args.setting);
		sprintf(root, "%s/r_%s", root, Args.run_num);
		mkdir(root, S_IRUSR | S_IWUSR | S_IXUSR);
		sprintf(str, "%s/README", root);
		write_params(str);
		print_params();
		printf("Reading training corpus ...\n");

		c_tr = (Corpus*)malloc(sizeof(Corpus));
		if(PARAMS.tr_percent < 1){
			c_te = (Corpus*)malloc(sizeof(Corpus));
		}else{
			c_te = NULL;
		}
		read_data(Args.data_name, PARAMS.train_data, PARAMS.I, c_tr, c_te, PARAMS.tr_percent);
		//save_local_vocabulary(ic, v, root);
		/*
		 * Training routing*/
		estimate(c_tr, c_te, v, root, Args.data_name, Args.z, Args.z_root);
		printf("End sampling the model ...\n");
		printf("Results saved at %s\n", root);
	}else if(Args.do_ei == 0){//calculate pmi, not used, ignore
		if(Args.mod_file == NULL){
			fprintf(stderr, "Error: enter trained model file using command -m !!!\n");
			exit(1);
		}
		printf("Calculating PMIs ...\n");
		read_pmi_params(Args.setting);
		calPMI(Args.mod_file, PARAMS.train_data, PARAMS.I, PARAMS.burn_in, Args.pmifile);
	}
	free_params(&PARAMS);
	return (0);
}
