#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXWIDTH 10000
#define BUFFSIZE 50000

double *table; // Z-table to speed up lookup

int partition(double *set, int *groups, int lo, int hi){
	double pivot, dummy;
	int i,j, foo;

	pivot = set[hi];
	i = lo - 1;
	for(j = lo; j < hi; j++){
		if(set[j] < pivot){
			i++;
			dummy  = set[i];
			set[i] = set[j];
			set[j] = dummy;

			foo       = groups[i];
			groups[i] = groups[j];
			groups[j] = foo;
		}
	}
	dummy    = set[i+1];
	set[i+1] = set[hi];
	set[hi]  = dummy;

	foo         = groups[i+1];
	groups[i+1] = groups[hi];
	groups[hi]  = foo;
	return(i+1);
}

/* Usage: array, groups, 0, len-1 */
void quicksort(double *set, int *groups, int lo, int hi){
	int p;

	if(lo < hi){
		p = partition(set, groups, lo, hi);
		quicksort(set, groups, lo, p - 1);
		quicksort(set, groups, p + 1, hi);
	}
}

float mid(int lo, int hi){
	return( (hi - lo) / 2.0 + lo + 1.0);
}

/* Convert an ordered list into ranks */
int rank(double *set, int num, double *rank){
	int i, loindex, j, ranki;

	loindex = 0;
	ranki = 0;

	for(i = 1; i < num; i++){
		if(set[i-1] != set[i]){
			if(i-loindex > 1){
				for(j = loindex; j < i; j++){
					rank[ranki] = mid(loindex, i - 1);
					ranki++;
				}
			}
			else{
				rank[ranki] = i;
				ranki++;
			}
			loindex = i;
		}
	}

	if(i - loindex > 1){
		for(j = loindex; j < i; j++){
			rank[ranki] = mid(loindex, i - 1);
			ranki++;
		}
	}
	else{
		rank[ranki] = i;
		ranki++;
	}
	return(0);
}

/* Porbability density function of standard normal distribution */
double PDF(double x){
	return(1.0 / sqrt(2.0 * M_PI) * exp(-1.0*x*x/2));
}

/* Z-test calculation. We store results in table to speed up calculations */
void fillZtable(){
	double z, s = 0.0, prev = 0.0, actual;
	int i;

	table = malloc(sizeof(double) * 60000);

	for(z = -6.0, i = 0; z < 0.0; z += 0.0001, i++){
		actual = PDF(z);
		s += 0.0001 * ( (prev + actual) / 2.0);
		table[i] = s;
		prev = actual;
	}
}

/* Mann-Whitney test */
double mannwhitney(double *set, int *groups, int num) {
	double *r, sa = 0.0, sb = 0.0, numa = 0.0, numb = 0.0;
	double Ua, Ub, U, z, p;
	int i;

	r = malloc(sizeof(double) * num);

	quicksort(set, groups, 0, num-1);
        rank(set, num, r);

	for(i = 0; i < num; i++){
		if(groups[i] == 0){
			sa += r[i];
			numa += 1.0;
		}
		else{
			sb += r[i];
			numb += 1.0;
		}
	}
	if(numa == 0.0 || numb == 0.0){
		printf("WARNING:No two groups\n");
		return(2.0);
	}

	Ua = sa - (numa * (numa + 1.0) / 2.0);
	Ub = sb - (numb * (numb + 1.0) / 2.0);

	U = Ub;
	if(Ua < Ub){
		U = Ua;
	}

	z = (U - (numa * numb / 2)) / sqrt(numa*numb*(numa+numb+1)/12.0);
	p = table[(int)((z + 6.0) * 10000.0)];

	free(r);
	return(p);
}

void onegroup(FILE *grpfile, char *rowid, FILE *valuefile, FILE *output){
	char *buffer;
	int buffsize = BUFFSIZE;
	char *actrowid;
	char *value;
	int *groups;
	int width = 0;
	int count;
	double *set;
	double pvalue;
	int i;
	int linenum = -1;
	int prevperc = 0;
	int perc;
	char found = 0;

	buffer = malloc(sizeof(char) * buffsize);
	groups = malloc(sizeof(int) * MAXWIDTH);
	set    = malloc(sizeof(double) * MAXWIDTH);

	while(fgets(buffer, buffsize, grpfile) != NULL){
		actrowid = strtok(buffer, "\t");
		if(!strcmp(actrowid, rowid)){
			// extract group information
			while(1){
				value = strtok(NULL, "\t\n");
				if(value == NULL){
					break;
				}
				if(!strcmp(value, "1")){
					groups[width] = 1;
				}
				else{
					groups[width] = 0;
				}
				width++;
			};
			found = 1;
			break;
		}
	}

	if(!found){
		printf("WARNING:Row not found\n");
		free(buffer);
		free(groups);
		free(set);
		return;
	}

	while(fgets(buffer, buffsize, valuefile) != NULL){
		linenum++;
	}
	rewind(valuefile);

	fgets(buffer, buffsize, valuefile); // read header
	fprintf(output, "Gene\tPvalue\n");
	for(i = 0; i < linenum; i++){
		fgets(buffer, buffsize, valuefile);
		// Do MannWhitney
		actrowid = strtok(buffer, "\t");
		count = 0;
		while(1){
			value = strtok(NULL, "\t\n");
			if(value == NULL){
				break;
			}
			set[count] = atof(value);
			count++;
		}
		pvalue = mannwhitney(set, groups, width);
		fprintf(output, "%s\t%f\n", actrowid, pvalue);

		perc = i * 100 / linenum;
		if( perc % 10 == 0 && perc != prevperc){
			printf("MESSAGE:%d%%\n", perc);
			prevperc = perc;
		}
	}

	free(buffer);
	free(groups);
	free(set);
}

void onevalue(FILE *valuefile, char *rowid, FILE *grpfile, FILE *output){
	char *buffer;
	int buffsize = BUFFSIZE;
	char *actrowid;
	int count = 0;
	double *set;
	char *value;
	double pvalue;
	int *groups;
	int i;
	int linenum = -1;
	int perc;
	int prevperc = 0;
	char found = 0;

	buffer = malloc(sizeof(char) * buffsize);
	set    = malloc(sizeof(double) * MAXWIDTH);
	groups = malloc(sizeof(int) * MAXWIDTH);

	while( fgets(buffer, buffsize, valuefile) != NULL){
		actrowid = strtok(buffer, "\t");
		if(!strcmp(actrowid, rowid)){
			while(1){
				value = strtok(NULL, "\t\n");
				if(value == NULL){
					break;
				}
				set[count] = atof(value);
				count++;
			}
			found = 1;
			break;
		}
	}

	if(!found){
		printf("WARNING:Row not found\n");
		free(buffer);
		free(set);
		free(groups);
		return;
	}

	while(fgets(buffer, buffsize, grpfile) != NULL){
		linenum++;
	}
	rewind(grpfile);

	fgets(buffer, buffsize, grpfile); // read header
	fprintf(output, "Gene\tPvalue\n");
	for(i = 0; i < linenum; i++){
		fgets(buffer, buffsize, grpfile);
		actrowid = strtok(buffer, "\t");
		count = 0;
		while(1){
			value = strtok(NULL, "\t\n");
			if(value == NULL){
				break;
			}
			if(!strcmp(value, "1")){
				groups[count] = 1;
			}
			else{
				groups[count] = 0;
			}
			count++;
		}
		pvalue = mannwhitney(set, groups, count);
		fprintf(output, "%s\t%f\n", actrowid, pvalue);

		perc = i * 100 / linenum;
		if( perc % 10 == 0 && perc != prevperc){
			printf("MESSAGE:%d%%\n", perc);
			prevperc = perc;
		}
	}

	free(buffer);
	free(set);
	free(groups);
}

/*
 * Usage: powermw onegroup|onevalue rowid groupfile valuefile resultfile
 */
int main(int argc, char **argv){
	FILE *grpfile;
	FILE *valuefile;
	FILE *output;

	if(argc != 6){
		printf("Usage: powermw onegroup|onevalue rowid groupfile valuefile resultfile\n");
		return(0);
	}

	grpfile   = fopen(argv[3], "r");
	valuefile = fopen(argv[4], "r");
	output    = fopen(argv[5], "w");
	fillZtable();

	if(!strcmp(argv[1], "onegroup")){
		onegroup(grpfile, argv[2], valuefile, output);
	}
	else{
		onevalue(valuefile, argv[2], grpfile, output);
	}

	fclose(grpfile);
	fclose(valuefile);
	fclose(output);

	free(table);

	return(0);
}
