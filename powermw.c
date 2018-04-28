#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

double mannwhitney(double *set, int *groups, int num) {
	double *r, sa = 0.0, sb = 0.0, numa = 0.0, numb = 0.0;
	double Ua, Ub, U, z;
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

	Ua = sa - (numa * (numa + 1.0) / 2.0);
	Ub = sb - (numb * (numb + 1.0) / 2.0);

	U = Ub;
	if(Ua < Ub){
		U = Ua;
	}

	z = (U - (numa * numb / 2)) / sqrt(numa*numb*(numa+numb+1)/12.0);

	free(r);
	return(z);
}

void onegroup(FILE *grpfile, char *colid, FILE *valuefile){
	while(fgets(buffer, buffsize, grpfile) != NULL){
		if(strncmp(buffer, colid) >= 0){
			// We found the actual group entry
		}
	}

	while(fgets(buffer, buffsize, valuefile) != NULL){
		// Do MannWhitney
	}
}

int main(int argc, char **argv){
	double *data;
	int  *groups;
	int i;
	FILE *grpfile;
	FILE *valuefile;
	char *buffer;

	grpfile   = fopen("", "r");
	valuefile = fopen("", "r");

	for(i = 0; i < 60000; i++){
		mannwhitney(data, groups, 31);
	}

	fclose(grpfile);
	fclose(valuefile);

	return(0);
}
