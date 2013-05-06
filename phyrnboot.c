#include <pthread.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <string.h>
#include <memory.h>
#include <math.h>
#include <assert.h>

#include "phyrn_phylip.h"

#define SUCCESS 0
#define FAILURE 1
#define NAME_SIZE 1000
#define NUM_SIZE 100
#define PHYLIP_QNAME_SIZE 10

char *PROGRAM_NAME = "phyrnboot";		// program name
char *DEFAULT_OUTFILE_NAME = "phyrnboot_outfile.txt";	// default output file name

char outfile[NAME_SIZE];			// output file
char ShortQueryNames[] = "ls_identifier_list.txt";
char mifile[NAME_SIZE]; 			// microarray file
char** onames = NULL;
char** qnames = NULL;				// query names
FILE *ofp = NULL;				

int boot_or_jack = 0;				// 0 : running bootstrap (default), 1: running jacknife
int use_allzero_cols = 0;			// 0 : sample all zero cols (default), 1: not sample all zero cols
double** signalVal = NULL;			// signal value for each profile (e.g. hits, product, etc.) 

double fracsample = 100;	// fraction of sample
long reps = 100;		// # of replicates
longer seed;
long spp = 0;			// # of queries
long sites = 0;			// # of columns in microarray (# of profile) e.g. ph_pfam_3p
long allzero_sites = 0;		// # of columns having all zero values in microarray

long output_cols = 8;		// # of columns of PHYLIP distance output 

/* show usage for program users */
void usage(int status){
	if(status == SUCCESS){
		printf("Usage: %s [-m NxM matrix (.tab) file] [-q #_of_queries] [-r #_of_replicates] [-s random_number_seed] [-p fraction_of_sample (default: 100)] [-o output_file (default: phyrnboot_outfile.txt)] [-b 0 for bootstrap/1 for jacknife (default: 0)] [-n 0 for sampling all zero cols, 1 for NOT sampling all zero cols (default: 0)\n", PROGRAM_NAME);
		printf("	-h help\n");
		printf("	-m microarray file\n");
		printf("	-q # of queries\n");
		printf("	-r # of replicates (must be positive)\n");
		printf("	-s random number seed (must be odd)\n");
		printf("	-p fractin of sample in % [0 - 100]\n");
		printf("	-o output file name (default: phyrnboot_outfile.txt)\n");
		printf("	-b 0 for boostrap/1 for jacknife (default: boostrap)\n");
		printf("	-n 0 for sampling all zero cols, 1 for NOT sampling all zero cols (default: 0)\n");
		printf("e.g.)\n");
		printf("%s -m product1.txt -q 40 -r 5000 -s 51\n", PROGRAM_NAME);
		printf("%s -m product1.txt -q 40 -r 5000 -s 73 -p 80 -b 1\n", PROGRAM_NAME);	
	}else if(status == FAILURE){
		fprintf(stderr, "Try '%s -h' for usage information.\n", PROGRAM_NAME);
	}

	exit(status);
}

/* print program running info */
void printRunningInfo(){
        long ss = 0;
        printf("\n%s\n", PROGRAM_NAME);
        if(boot_or_jack == 0){
                printf("program option : bootstrap\n");
        }else{
                printf("program option: jacknife\n");
        }
        printf("# of queries : %ld\n", spp);
        printf("# of sites : %ld\n", sites);
        printf("# of all zero sites : %ld\n", allzero_sites);

        if(use_allzero_cols == 0){
                printf("sample all zero columns\n");
                ss = fracsample * sites;
        }else{
                printf("NOT sample all zero columns\n");
                ss = fracsample * (sites - allzero_sites);
        }
        printf("fraction of sample : %ld\%\%\n", (long)(fracsample * 100));
        printf("Max. # of sample size : %ld\n\n", ss);

}


void openOutfile(){
	if((ofp=fopen(outfile,"w"))==NULL){
		fprintf(stderr, "%s cannot be read.\n", outfile);
		exit(0);
	}	
}

void closeOutfile(){
	fclose(ofp);
}

void allocSignalVal(){
	long i = 0;
	
	signalVal = calloc(spp, sizeof(double *));

	for(i=0; i<spp; i++){
		signalVal[i] = (double *)calloc(sites, sizeof(double));	
	}

	qnames = calloc(spp, sizeof(char *));
	onames = calloc(spp, sizeof(char *));
	
	for(i=0; i<spp; i++){
		qnames[i] = (char *)calloc(1, PHYLIP_QNAME_SIZE);
		onames[i] = (char *)calloc(1, 250);
	}	
}

void freeSignalVal(){
	long i=0;
	for(i=0; i<spp; i++){
		free(signalVal[i]);
	}

	free(signalVal);

	for(i=0; i<spp; i++){
		free(qnames[i]);
		free(onames[i]);
	}
	
	free(qnames);
	free(onames);
}

/* read microarray data and calculate distance between each pair of queries */
void readMicroarrayData(){
	FILE *fp;
	FILE *sqn = NULL;
	char c;
	int k;
	long i, j, idx, p;
	char val[NUM_SIZE];

	if((fp=fopen(mifile,"r")) == NULL){
		fprintf(stderr, "%s cannot be read.\n", mifile);
		exit(0);
	}

	/* skip headers */
	for(i=0; i<1; i++){
		do{
			c = fgetc(fp);
		}while(c != '\t');
	}

	/* count # of profiles (# of columns) */
	do{
		c = fgetc(fp);
		if(c == '\t') sites++;	
	}while(c != '\n' && (int)c != 13);	
	sites++;

	/* skip the second line */	
	/*do{
		c = fgetc(fp);
	}while(c != '\n');*/

	allocSignalVal();
	
	if((sqn=fopen(ShortQueryNames,"w"))==NULL){
		fprintf(stderr, "%s cannot be read.\n", ShortQueryNames);
		exit(0);
	}

	for(i=0; i<spp; i++){
		k = 0;

		/*do{
			c = fgetc(fp);
		}while(c != '\t');

		qnames[i][k++] = fgetc(fp);*/
		
		//column1
		do{
			c = fgetc(fp);
			if (k==0 && (c == '\n' || (int)c == 13)) {
				c = fgetc(fp);
			}
			if(c != '\t') onames[i][k++] = c;
		}while(c != '\t');
		
		k = sprintf(qnames[i], "j%d", i+1);
		
		
		fprintf(sqn, "%s\t%s\n", qnames[i], onames[i]);
		
		while(k < PHYLIP_QNAME_SIZE){
			qnames[i][k++] = ' ';
		}
		
		/*for(j=0; j<2; j++){
			do{
				c = fgetc(fp);
			}while(c != '\t');
		}*/

		for(p=0; p<sites; p++){
			idx = 0;
			memset(val, 0, NUM_SIZE);	
			do{
				c = fgetc(fp);
				if(c != '\t' && c != '\n'  && (int)c != 13) val[idx] = c;
				idx++;
			}while(c != '\t' && c != '\n' && (int)c != 13);

			signalVal[i][p] = atof(val) / 10000;	
		}
	}	

	fclose(sqn);
	
	//printf("sites = %ld\n", sites);
	fclose(fp);

	
}

int isAllZeroCol(long profile_idx){
	int allzero = 1;
	long k;

	for(k=0; k<spp; k++){
		if(signalVal[k][profile_idx] > 0){
			allzero = 0;
			break; 	
		} 
	}

	return allzero;	
}


/* sets up weights by resampling data with boostrap */
void bootweights(long* weight){
	long j, i;
	long sample_size;

	if(use_allzero_cols == 1){
		sample_size = fracsample * (sites - allzero_sites); 	
	}else{
		sample_size = fracsample * sites;
	}

	for(i=0; i<sample_size; i++){
//		weight[i]++;

		j = (long)(sites * randum(seed));
		
		if(use_allzero_cols == 1){
			// if selecting all zero col, resample
			while(isAllZeroCol(j)){
				j = (long)(sites * randum(seed));
			}	
		}

		weight[j]++;	

	}
}

/* sets up weights by resampling data with jacknife */
void jackweights(long* weight){
	long j, i;
        long sample_size;

	if(use_allzero_cols == 1){
                sample_size = fracsample * (sites - allzero_sites);
        }else{
                sample_size = fracsample * sites;
        }

	i = 0;
	while(i < sample_size){
		j = (long)(sites * randum(seed));
		
		if(use_allzero_cols == 1){
			// if selecting all zero col, resample
                        while(isAllZeroCol(j)){
                                j = (long)(sites * randum(seed));
                        }
                }

		if(weight[j] == 0){
			weight[j] = 1;
			i++;
		}
	}
}

/* calculate eucledian distance between two sequences */
double getEucledianDistTwoSeqs(long q1, long q2, long* weight){
	long i, w;
	double dist = 0;

	for(i=0; i<sites; i++){
		w = weight[i];
		if(w > 0){
			dist += (pow((signalVal[q1][i] - signalVal[q2][i]), 2)*w); 
		}	
	}

	dist = sqrt(dist);

	return dist;
}

/* calculate eucledian distance using sampled data from bootstrap */
double** calculateEucledianDist(long* weight){
	long i, p, q, k;
	double distVal = 0;
	double** dist = NULL;

	dist = calloc(spp, sizeof(double*));
	for(i=0; i<spp; i++){
		dist[i] = (double*)calloc(spp, sizeof(double));
	}

	for(p=0; p<spp; p++){
		for(k=0; k<(p+1); k++){
			dist[p][k] = dist[k][p];
		}
		for(q=(p+1); q<spp; q++){
			distVal = getEucledianDistTwoSeqs(p, q, weight);
			dist[p][q] = distVal;	
		}
	}
/*
	for(p=0; p<spp; p++){
		for(q=0; q<spp; q++){
			printf("%lf\t", p, q, dist[p][q]);
		}
		printf("\n");
	}
*/
	return dist;
}

/* write distance to a file */
void writeDist(double** dist){
	long s, i, si, li, tm;

	fprintf(ofp, "   %ld\n", spp);

	for(s=0; s<spp; s++){
		fprintf(ofp, "%s", qnames[s]);
		si = 1;
        	li = output_cols;
		do{
			si--;
			li--;
			if(li > spp){
				li = spp;
			}
	
			for(i=si; i<li; i++){
				fprintf(ofp, " %.15lf", dist[s][i]);		
			}
			fprintf(ofp, "\n");

			if(li == spp) break;

			si += output_cols;
			li += output_cols;
		}while(1);
	}
	fprintf(ofp, "\n");	
} 

/* write distance to a file in MEGA format */
void writeDist_MEGA_FORMAT(double** dist){
        long s, i, si, li, tm;

	fprintf(ofp, "#mega\n");
	fprintf(ofp, "!TITLE  (%ld x %ld) distance matrix;\n", spp, spp);
	fprintf(ofp, "!Format DataType=distance;\n");
	fprintf(ofp, "!Description Converted for MEGA4;\n\n");

	for(s=0; s<spp; s++){
		fprintf(ofp, "#%s\n", qnames[s]);
	}
	fprintf(ofp, "\n\n");

	li = 1;
        for(s=1; s<spp; s++){
		for(i=0; i<li; i++){
               		fprintf(ofp, "   %lf", dist[s][i]);
                }
                fprintf(ofp, "\n");
		li++;
        }
}

 
void startBootDist(){
	double** dist = NULL;
	long* weight = NULL;
	long i;

	weight = (long*)calloc(sites, sizeof(long));
	memset(weight, 0, sites);

	/* sets up weights by resampling data */
	if(boot_or_jack == 0){
		bootweights(weight);
	}else{
		jackweights(weight);
	}

	/* calculate eucledian distance using sampled data from bootstrap */
	dist = calculateEucledianDist(weight);

	/* write distance to a file */
	writeDist(dist);
	//writeDist_MEGA_FORMAT(dist);
		
	free(weight);	
}

/* count # of columns having all zeros */
void countAllZeroCols(){
	long i, cnt;

	cnt = 0;
	for(i=0; i<sites; i++){
		if(isAllZeroCol(i)){
			cnt++;
		}	
	}

	allzero_sites = cnt;
}

/* start GDDA bootstrap */
void startGDDAbootstrap(){
	int rr = 0;

	/* read microarray data and calculate distance between each pair of queries */
	readMicroarrayData();

	/* count # of columns having all zeros */
	countAllZeroCols();

	/* print program running info */
	printRunningInfo();

	openOutfile();

	/* do bootstrapping # of replicates times*/
	for(rr=0; rr<reps; rr++){
		/* do bootstrapping and calculate distances */	
		startBootDist();
		printf("No.%d replicate done.\n", (rr+1));
	}

	closeOutfile();

	freeSignalVal();
}

/* main */
int main(int argc, char *argv[]){
        time_t s_time, f_time;
	int c, pflg = 0, oflg = 0, cmdformat = SUCCESS;
	long inseed, inseed0;
	extern char *optarg, optopt;
//	char mifile[NAME_SIZE];	// microarray file

        (void) time(&s_time);   // start time

	const char *optstring = "hm:p:r:s:q:o:b:n:";
	PROGRAM_NAME = argv[0];
	
	if(!argv[1]){
		fprintf(stderr, "%s: needs option and file name.\n", PROGRAM_NAME);
		usage(FAILURE);
	}
	
	while((c=getopt(argc, argv, optstring)) != EOF){
		switch(c){
			case 'h':
				usage(SUCCESS);
				break;
			case 'm':
				strcpy(mifile, optarg);					
				break;
			case 'q':
				spp = atol(optarg);
				break;
			case 'r':
				reps = atol(optarg);
				break;
			case 's':
				inseed = atol(optarg);
				break;
			case 'p':
				pflg++;
				fracsample = atof(optarg);
				break;
			case 'o':
				oflg++;
				strcpy(outfile, optarg);
				break;
			case 'b':
				boot_or_jack = atoi(optarg);
				break;
			case 'n':
				use_allzero_cols = atoi(optarg);
				break;
			case ':':
				fprintf(stderr, "%s: Option -%c requires an option-argument.\n",
					PROGRAM_NAME, optopt);
				break;
			case '?':
				fprintf(stderr,"%s: Unrecognized option.\n",PROGRAM_NAME);
				usage(FAILURE);
				break;			
		}
	}

	if(pflg && (fracsample <= 0.0 || fracsample > 100.0)){
		printf("Fraction of sample: out of range [0 - 100]\n");
		cmdformat = FAILURE;
		//usage(FAILURE);
	}else{
		fracsample = fracsample/100.0;
	}

	if(reps <= 0){
		printf("# of replicates: must be positive\n");
		cmdformat = FAILURE;
		//usage(FAILURE);
	}

	if(inseed > 0 && (inseed & 0x1)){
		initseed(&inseed, &inseed0, seed);
	}else{
		printf("Random number seed: must be odd\n");
		cmdformat = FAILURE;
		//usage(FAILURE);
	}
	
	if(boot_or_jack != 1 && boot_or_jack != 0){
                printf("Boot/Jacknife: must be 0 or 1\n");
                cmdformat = FAILURE;
        }

	if(use_allzero_cols != 1 && use_allzero_cols != 0){
                printf("Boot/Jacknife: must be 0 or 1\n");
                cmdformat = FAILURE;
        }

	if(cmdformat){
		usage(FAILURE);
	}

	if(oflg == 0){
		strcpy(outfile, DEFAULT_OUTFILE_NAME);
	}

	/* run GDDA bootstrap */
	startGDDAbootstrap();

        (void) time(&f_time);   // finish time

        printf("\nProgram compeleted!\nTotal time : %d seconds\n",(int)f_time-s_time);

}

