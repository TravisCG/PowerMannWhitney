#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXWIDTH 10000
#define BUFFSIZE 50000

double *table; // Z-table to speed up lookup
enum filtertype {none, include, exclude};

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
		s += 0.0001 * ( (prev + actual) / 2.0); // Numerical integral calculation
		table[i] = s * 2.0; // Two-tail p value
		prev = actual;
	}
}

/* Mann-Whitney test */
double mannwhitney(double *set, int *groups, int num, double *log2foldch) {
	double *r, sa = 0.0, sb = 0.0, numa = 0.0, numb = 0.0;
	double s1 = 0.0, s2 = 0.0;
	double Ua, Ub, U, z, p;
	int i;

	r = malloc(sizeof(double) * num);
	quicksort(set, groups, 0, num-1);

        rank(set, num, r);

	for(i = 0; i < num; i++){
		if(groups[i] == 0){
			sa += r[i];
			s1 += set[i];
			numa += 1.0;
			printf("%f:0\n", set[i]);
		}
		else{
			sb += r[i];
			s2 += set[i];
			numb += 1.0;
			printf("%f:1\n", set[i]);
		}
	}

	free(r);

	Ua = sa - (numa * (numa + 1.0) / 2.0);
	Ub = sb - (numb * (numb + 1.0) / 2.0);

	if(numa == 0.0 || numb == 0.0){
		// Avoiding division by zero, return value is non-significant
		*log2foldch = 0.0;
		return(1.0);
	}

	U = Ub;
	if(Ua < Ub){
		U = Ua;
	}

	s1 = s1 / numa;
	s2 = s2 / numb;

	*log2foldch = log2(s1 / s2);

	z = (U - (numa * numb / 2)) / sqrt(numa*numb*(numa+numb+1)/12.0);
	if(z < -6.0){
		z = -6.0;
	}
	if(z > 0.0){
		z = 0.0;
	}
	p = table[(int)((z + 6.0) * 10000.0)];

	return(p);
}

/* Just count the lines in a text file and revert the file pointer to the beginning */
int countlines(FILE *f, char *buffer, int buffsize){
	int linenum = -1;
	while(fgets(buffer, buffsize, f) != NULL){
		linenum++;
	}
	rewind(f);

	return(linenum);
}

/* Just count the lines in a text file and revert the file pointer to the beginning */
int countlineswithfilter(FILE *f, char *buffer, int buffsize, enum filtertype filt, char *rowid2, char *impcols){
	int linenum = -1, cols = 0, impcolnum = 0;
	char *fields;
	char found = 0;

	while(fgets(buffer, buffsize, f) != NULL){
		linenum++;
		if(filt != none){
			fields = strtok(buffer, "\t");
			if(!strcmp(fields, rowid2)){
				while(1){
					fields = strtok(NULL, "\t\n");
					if(fields == NULL){
						break;
					}
					if( (fields[0] == '1' && filt == include) || (fields[0] == '0' && filt == exclude)){
						impcols[cols] = 1;
						impcolnum++;
					}
					else{
						impcols[cols] = 0;
					}
					cols++;
				}
				found = 1;
			}
		}
	}
	rewind(f);

	if(filt != none && found == 0){
		printf("WARNING:Filter gene not found\n");
		return(-1);
	}
	if(filt != none && impcolnum < 20){
		printf("WARNING:Sample number is less than 20\n");
		return(-2);
	}

	return(linenum);
}

/* Simple progress inficator */
void progress(int i, int max){
	static int prevperc = 0;
	int perc;

	perc = i * 100 / max;
	if( perc % 10 == 0 && perc != prevperc){
		printf("MESSAGE:%d%%\n", perc);
		prevperc = perc;
	}
}

int parsegroup(int *groups){
	int count = 0;
	char *value;

	while(1){
		value = strtok(NULL, "\t\n");
		if(value == NULL){
			break;
		}
		if(value[0] == '1'){
			groups[count] = 1;
		}
		else{
			groups[count] = 0;
		}
		count++;
	}

	return(count);
}

int parsevalue(double *set){
	int count = 0;
	char *value;

	while(1){
		value = strtok(NULL, "\t\n");
		if(value == NULL){
			break;
		}
		set[count] = strtod(value, NULL);
		count++;
	}
	return(count);
}

int parseimpcol(char *cols, enum filtertype f){
	int count = 0;
	char *value;

	while(1){
		value = strtok(NULL, "\t\n");
		if(value == NULL){
			break;
		}
		if( (value[0] == '1' && f == include) || (value[0] == '0' && f == exclude)){
			cols[count] = 1;
		}
		else{
			cols[count] = 0;
		}
		count++;
	}
	return(count);
}

int reorder(double *set, int *groups, char *impcols, int width){
	int i, j = 0;
	double d1;
	char d2;

	for(i = 0; i < width; i++){
		if(impcols[i] == 1){
			d1 = set[i];
			set[i] = set[j];
			set[j] = d1;

			d2 = groups[i];
			groups[i] = groups[j];
			groups[j] = d2;
			j++;
		}
	}
	return(j);
}

void onegroup(FILE *grpfile, char *rowid, FILE *valuefile, FILE *output, enum filtertype f, char *rowid2){
	char *buffer, *importantcols;
	int buffsize = BUFFSIZE;
	char *actrowid;
	int *groups;
	int *cpygroups;
	int width = 0;
	double *set;
	double pvalue, log2fc = 100.0, bonferroni = 1.0;
	int i;
	int linenum, impcolcount = 0;
	char foundr1 = 0, foundr2 = 1;

	buffer    = malloc(sizeof(char) * buffsize);
	groups    = malloc(sizeof(int) * MAXWIDTH);
	set       = malloc(sizeof(double) * MAXWIDTH);
	cpygroups = malloc(sizeof(int) * MAXWIDTH);
	importantcols = malloc(sizeof(char) * MAXWIDTH);

	if(f != none){
		foundr2 = 0;
	}

	while(fgets(buffer, buffsize, grpfile) != NULL){
		actrowid = strtok(buffer, "\t");
		if(!strcmp(actrowid, rowid)){
			// extract group information
			width = parsegroup(groups);
			foundr1 = 1;
		}
		if(f != none && !strcmp(actrowid, rowid2)){
			impcolcount = parseimpcol(importantcols, f);
			foundr2 = 1;
		}
		if(foundr1 == 1 && foundr2 == 1){
			break;
		}
	}

	if(!foundr1 || !foundr2){
		printf("WARNING:Row not found\n");
		free(buffer);
		free(groups);
		free(set);
		free(cpygroups);
		free(importantcols);
		return;
	}

	if(f != none && impcolcount < 20){
		printf("WARNING:Number of samples less than 20\n");
		free(buffer);
		free(groups);
		free(set);
		free(cpygroups);
		free(importantcols);
		return;
	}

	linenum = countlines(valuefile, buffer, buffsize);

	fgets(buffer, buffsize, valuefile); // read header
	fprintf(output, "Gene\tlog2FC\tPvalue\tBonferroni\n");
	for(i = 0; i < linenum; i++){
		fgets(buffer, buffsize, valuefile);
		// Do MannWhitney
		actrowid = strtok(buffer, "\t");
		width = parsevalue(set);
		memcpy(cpygroups, groups, sizeof(int) * width);
		if(f != none){
			width = reorder(set, cpygroups, importantcols, width);
		}
		pvalue = mannwhitney(set, cpygroups, width, &log2fc);
		bonferroni = pvalue * linenum;
		if(bonferroni > 1.0) bonferroni = 1.0;
		fprintf(output, "%s\t%f\t%f\t%f\n", actrowid, log2fc, pvalue, bonferroni);
		progress(i, linenum);
	}

	free(buffer);
	free(groups);
	free(cpygroups);
	free(set);
	free(importantcols);
}

void onevalue(FILE *valuefile, char *rowid, FILE *grpfile, FILE *output, enum filtertype f, char *rowid2, int mutprevalence){
	char *buffer, *importantcols;
	int buffsize = BUFFSIZE;
	char *actrowid;
	int count = 0;
	double *set;
	double *cpyset;
	double pvalue, log2fc = 100.0, bonferroni = 1.0;
	int *groups;
	int i;
	int linenum = -1;
	char found = 0;

	buffer = malloc(sizeof(char) * buffsize);
	set    = malloc(sizeof(double) * MAXWIDTH);
	cpyset = malloc(sizeof(double) * MAXWIDTH);
	groups = malloc(sizeof(int) * MAXWIDTH);
	importantcols = malloc(sizeof(char) * MAXWIDTH);

	while( fgets(buffer, buffsize, valuefile) != NULL){
		actrowid = strtok(buffer, "\t");
		if(!strcmp(actrowid, rowid)){
			count = parsevalue(set);
			found = 1;
			break;
		}
	}

	if(!found){
		printf("WARNING:Row not found\n");
		free(buffer);
		free(set);
		free(cpyset);
		free(groups);
		free(importantcols);
		return;
	}

	linenum = countlineswithfilter(grpfile, buffer, buffsize, f, rowid2, importantcols);
	if(linenum < 0){
		free(buffer);
		free(set);
		free(cpyset);
		free(groups);
		free(importantcols);
		return;
	}

	fgets(buffer, buffsize, grpfile); // read header
	fprintf(output, "Gene\tlog2FC\tPvalue\tBonferroni\n");
	for(i = 0; i < linenum; i++){
		fgets(buffer, buffsize, grpfile);
		actrowid = strtok(buffer, "\t");
		count = parsegroup(groups);
		memcpy(cpyset, set, sizeof(double) * count);
		if(f != none){
			count = reorder(cpyset, groups, importantcols, count);
		}
		pvalue = mannwhitney(cpyset, groups, count, &log2fc);
		bonferroni = pvalue * linenum;
		if(bonferroni > 1.0) bonferroni = 1.0;
		fprintf(output, "%s\t%f\t%f\t%f\n", actrowid, log2fc, pvalue, bonferroni);
		progress(i, linenum);
	}

	free(buffer);
	free(set);
	free(cpyset);
	free(groups);
	free(importantcols);
}

/*
 * Usage: powermw onegroup|onevalue rowid groupfile valuefile resultfile
 */
int main(int argc, char **argv){
	FILE *grpfile;
	FILE *valuefile;
	FILE *output;
	char *rowid2;
	enum filtertype filter;

	if(argc < 6){
		printf("Usage: powermw onegroup|onevalue rowid groupfile valuefile resultfile include|exclude rowid2 mutprevalence\n");
		return(0);
	}

	grpfile   = fopen(argv[3], "r");
	valuefile = fopen(argv[4], "r");
	output    = fopen(argv[5], "w");
	if(argc >= 8 && !strcmp(argv[6], "exclude")){
		filter = exclude;
		rowid2 = argv[7];
	} 
	else if(argc >= 8 && !strcmp(argv[6], "include")){
		filter = include;
		rowid2 = argv[7];
	}
	else{
		filter = none;
		rowid2 = NULL;
	}


	fillZtable();

	if(!strcmp(argv[1], "onegroup")){
		onegroup(grpfile, argv[2], valuefile, output, filter, rowid2);
	}
	else{
		onevalue(valuefile, argv[2], grpfile, output, filter, rowid2, atoi(argv[8]));
	}

	fclose(grpfile);
	fclose(valuefile);
	fclose(output);

	free(table);

	return(0);
}
