#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAXWIDTH 10000
#define BUFFSIZE 50000
#define SHOW_NO_MUTPREV 0
#define SHOW_MUTPREV 1

double *table; // Z-table to speed up lookup
enum filtertype {none, include, exclude};

typedef struct Result {
	char genename[30];
	double FC;
	double p;
	double mutprev;
	double mutexp;
	double wtexp;
	double FDR;
} Result;

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

void swapres(Result *r, int i, int j){
	Result dummy;

	dummy = r[i];
	r[i]  = r[j];
	r[j]  = dummy;
}

int partres_desc(Result *r, int lo, int hi){
	double pivot;
	int i, j;

	pivot = r[hi].p;
	i = lo - 1;
	for(j = lo; j < hi; j++){
		if(r[j].p > pivot){
			i++;
			swapres(r, i, j);
		}
	}
	swapres(r, i + 1, hi);
	return(i+1);
}

int partres_asc(Result *r, int lo, int hi){
	double pivot;
	int i, j;

	pivot = r[hi].p;
	i = lo - 1;
	for(j = lo; j < hi; j++){
		if(r[j].p < pivot){
			i++;
			swapres(r, i, j);
		}
	}
	swapres(r, i + 1, hi);
	return(i+1);

}

void sortresult(Result *r, int lo, int hi, int order){
	int p;

	if(lo < hi){
		if(order == 0){
			p = partres_desc(r, lo, hi);
		}
		else{
			p = partres_asc(r, lo, hi);
		}
		sortresult(r, lo, p - 1, order);
		sortresult(r, p + 1, hi, order);
	}
}

float mid(int lo, int hi){
	return( (hi - lo) / 2.0 + lo + 1.0);
}

double min(double a, double b){
	return( a < b ? a : b);
}

double max(double a, double b){
	return( a < b ? b : a);
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

/* Print out the results in a pretty way */
void prettyprint(FILE *output, Result res, int show, double plimit){
	if(res.p >= plimit){
		return;
	}
	fprintf(output, "%s\t%.2f\t", res.genename, res.FC);
	if(res.p == table[0]){
		fprintf(output, "<");
	}
	if(res.p < 0.001){
		fprintf(output, "%.1e\t", res.p);
	}
	else{
		fprintf(output, "%.3g\t", res.p);
	}
	if(res.FDR < 0.001){
		fprintf(output, "%.2e", res.FDR);
	}
	else{
		fprintf(output, "%.3g", res.FDR);
	}
	if(show == SHOW_NO_MUTPREV){
		fprintf(output, "\n");
	}
	else{
		fprintf(output, "\t%.1f%%\t%.2f\t%.2f\n", res.mutprev * 100.0, res.mutexp, res.wtexp);
	}
}

/* Mann-Whitney test */
double mannwhitney(double *set, int *groups, int num, double *foldch, double *mutexp, double *wtexp) {
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
		}
		else{
			sb += r[i];
			s2 += set[i];
			numb += 1.0;
		}
	}

	free(r);

	Ua = sa - (numa * (numa + 1.0) / 2.0);
	Ub = sb - (numb * (numb + 1.0) / 2.0);

	if(numa == 0.0 || numb == 0.0){
		// Avoiding division by zero, return value is non-significant
		*foldch = 0.0;
		return(1.0);
	}

	U = Ub;
	if(Ua < Ub){
		U = Ua;
	}

	s1 = s1 / numa;
	s2 = s2 / numb;

	*foldch = s2 / s1;
	*mutexp = s2;
	*wtexp  = s1;

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
		printf("MESSAGE:%d\n", perc);
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

// Calculate the ratio of mutant samples
double calcmutprev(int *groups, int count, int *mut){
	int i, m = 0;

	for(i = 0; i < count; i++){
		if(groups[i] == 1){
			m++;
		}
	}

	*mut = m;
	return((double)m / (double)count);
}

void storeres(Result *r, char *rowid, double fc, double p, double m, double me, double we){
	strcpy(r->genename, rowid);
	r->FC      = fc;
	r->p       = p;
	r->mutprev = m;
	r->mutexp  = me;
	r->wtexp   = we;
}

void onegroup(FILE *grpfile, char *rowid, FILE *valuefile, FILE *output, enum filtertype f, char *rowid2, double foldlimit, double plimit){
	char *buffer, *importantcols;
	int buffsize = BUFFSIZE;
	char *actrowid;
	int *groups;
	int *cpygroups;
	int width = 0;
	double *set, mutexp = 0.0, wtexp = 0.0;
	double pvalue, fc = 100.0, alpha, minalpha;
	int i;
	int linenum, impcolcount = 0;
	char foundr1 = 0, foundr2 = 1;
	Result *res;
	int resnum = 0;

	buffer        = malloc(sizeof(char) * buffsize);
	groups        = malloc(sizeof(int) * MAXWIDTH);
	set           = malloc(sizeof(double) * MAXWIDTH);
	cpygroups     = malloc(sizeof(int) * MAXWIDTH);
	importantcols = malloc(sizeof(char) * MAXWIDTH);

	printf("MESSAGE:Measure input size\n");
	linenum = countlines(valuefile, buffer, buffsize);

	res = malloc(sizeof(Result) * linenum);

	if(f != none){
		foundr2 = 0;
	}

	printf("MESSAGE:Finding gene\n");
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
		goto end;
	}

	if(f != none && impcolcount < 20){
		printf("WARNING:Number of samples less than 20\n");
		goto end;
	}

	buffer = fgets(buffer, buffsize, valuefile); // read header
	fprintf(output, "Gene\tFold change\tP value\tFDR\n");
	for(i = 0; i < linenum; i++){
		buffer = fgets(buffer, buffsize, valuefile);
		// Do MannWhitney
		actrowid = strtok(buffer, "\t");
		width = parsevalue(set);
		memcpy(cpygroups, groups, sizeof(int) * width);
		if(f != none){
			width = reorder(set, cpygroups, importantcols, width);
		}
		pvalue = mannwhitney(set, cpygroups, width, &fc, &mutexp, &wtexp);
		// Don't ask what the hell it is. 
		if(foldlimit < 1.0){
			foldlimit = foldlimit + 1.0;
		}
		if(fc > max(foldlimit, 1.0 / foldlimit) || fc < min(foldlimit, 1.0 / foldlimit)){
			storeres(&res[resnum], actrowid, fc, pvalue, 0, mutexp, wtexp);
			resnum++;
		}
		progress(i, linenum);
	}

	sortresult(res, 0, resnum-1, 0);
	minalpha = 1.0;
	for(i = 0; i < resnum; i++){
		// Benjamini-Hochberg procedure
		alpha = (double)resnum / (double)(resnum - i) * res[i].p;
		if(alpha < minalpha){
			minalpha = alpha;
		}
		res[i].FDR = minalpha;
	}

	sortresult(res, 0, resnum - 1, 1);
	for(i = 0; i < resnum; i++){
		prettyprint(output, res[i], SHOW_NO_MUTPREV, plimit);
	}
end:
	free(res);
	free(buffer);
	free(groups);
	free(cpygroups);
	free(set);
	free(importantcols);
}

void onevalue(FILE *valuefile, char *rowid, FILE *grpfile, FILE *output, enum filtertype f, char *rowid2, double foldlimit, double plimit, double mutplimit){
	char *buffer, *importantcols;
	int buffsize = BUFFSIZE;
	char *actrowid;
	int count = 0;
	double *set, mutexp = 0.0, wtexp = 0.0;
	double *cpyset;
	double pvalue, fc = 0.0, mutprev, alpha, minalpha;
	int *groups;
	int i, mutnum;
	int linenum = -1;
	char found = 0;
	Result *res;
	int resnum = 0;

	buffer = malloc(sizeof(char) * buffsize);
	set    = malloc(sizeof(double) * MAXWIDTH);
	cpyset = malloc(sizeof(double) * MAXWIDTH);
	groups = malloc(sizeof(int) * MAXWIDTH);
	importantcols = malloc(sizeof(char) * MAXWIDTH);

	printf("MESSAGE:Measure input size\n");
	linenum = countlineswithfilter(grpfile, buffer, buffsize, f, rowid2, importantcols);
	if(linenum < 0){
		goto end;
	}
	res = malloc(sizeof(Result) * linenum);

	printf("MESSAGE:Finding gene\n");
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
		goto end;
	}

	buffer = fgets(buffer, buffsize, grpfile); // read header
	fprintf(output, "Gene\tFold change\tP value\tFDR\tMutation Prevalence\tAverage mutant expression\tAverage wild type expression\n");
	for(i = 0; i < linenum; i++){
		buffer = fgets(buffer, buffsize, grpfile);
		actrowid = strtok(buffer, "\t");
		count = parsegroup(groups);
		memcpy(cpyset, set, sizeof(double) * count);
		if(f != none){
			count = reorder(cpyset, groups, importantcols, count);
		}

		mutprev = calcmutprev(groups, count, &mutnum);
		if(mutprev < mutplimit){
			continue;
		}

		pvalue = mannwhitney(cpyset, groups, count, &fc, &mutexp, &wtexp);
		if(foldlimit < 1.0){
			foldlimit = foldlimit + 1.0;
		}
		if(fc > max(foldlimit, 1.0 / foldlimit) || fc < min(foldlimit, 1.0 / foldlimit)){
			storeres(&res[resnum], actrowid, fc, pvalue, mutprev, mutexp, wtexp);
			resnum++;
		}
		progress(i, linenum);
	}

	sortresult(res, 0, resnum-1, 0);
	minalpha = 1.0;

	for(i = 0; i < resnum; i++){
		// Benjamini-Hochberg procedure
		alpha = (double)resnum / (double)(resnum - i) * res[i].p;
		if(alpha < minalpha){
			minalpha = alpha;
		}
		res[i].FDR = minalpha;
	}

	sortresult(res, 0, resnum - 1, 1);
	for(i = 0; i < resnum; i++){
		prettyprint(output, res[i], SHOW_MUTPREV, plimit);
	}
end:
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
	char *rowid1 = NULL, *rowid2 = NULL;
	enum filtertype filter = none;
	int i;

	char *type = NULL;
	char *grpfilename = NULL;
	char *valuefilename = NULL, *outname = NULL;
	double mutplim = 0.01, foldlimit = 0.5, plimit = 0.05;

	for(i = 1; i < argc; i++){
		if(!strcmp(argv[i], "-t")){
			type = argv[i+1];
		}
		if(!strcmp(argv[i], "-r")){
			rowid1 = argv[i+1];
		}
		if(!strcmp(argv[i], "-g")){
			grpfilename = argv[i+1];
		}
		if(!strcmp(argv[i], "-v")){
			valuefilename = argv[i+1];
		}
		if(!strcmp(argv[i], "-o")){
			outname = argv[i+1];
		}
		if(!strcmp(argv[i], "-f")){
			if(!strcmp(argv[i+1], "include")){
				filter = include;
			}
			else if(!strcmp(argv[i+1], "exclude")){
				filter = exclude;
			}
		}
		if(!strcmp(argv[i], "-b")){
			rowid2 = argv[i+1];
		}
		if(!strcmp(argv[i], "-h")){
			printf("Usage: powermw -t onegroup|onevalue -r rowid -g groupfile -v valuefile -o resultfile -f include|exclude|none -b rowid2 -m 0.1 -l 0.5 -p 0.05\n");
			return(0);
		}
		if(!strcmp(argv[i], "-m")){
			mutplim = atof(argv[i+1]);
		}
		if(!strcmp(argv[i], "-l")){
			foldlimit = atof(argv[i+1]);
		}
		if(!strcmp(argv[i], "-p")){
			plimit = atof(argv[i+1]);
		}
	}

	grpfile   = fopen(grpfilename, "r");
	valuefile = fopen(valuefilename, "r");
	output    = fopen(outname, "w");

	fillZtable();

	if(!strcmp(type, "onegroup")){
		onegroup(grpfile, rowid1, valuefile, output, filter, rowid2, foldlimit, plimit);
	}
	else{
		onevalue(valuefile, rowid1, grpfile, output, filter, rowid2, foldlimit, plimit, mutplim);
	}

	fclose(grpfile);
	fclose(valuefile);
	fclose(output);

	free(table);

	return(0);
}
