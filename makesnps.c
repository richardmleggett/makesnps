/*----------------------------------------------------------------------*
 * File:    makesnpss.c                                                 *
 *                                                                      *
 * Purpose: Make a copy of a reference with SNPs in it.                 *
 *                                                                      *
 * Author:  Richard Leggett                                             *
 *          The Genome Analysis Centre                                  *
 *          Norwich Research Park, Colney, Norwich, NR4 7UH, UK         *
 *          richard.leggett@bbsrc.ac.uk                                 * 
 *                                                                      *
 * History: 25-May-11: RML: Created                                     *
 *----------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <stdint.h>
#include <math.h>

/*----------------------------------------------------------------------*
 * Constants                                                            *
 *----------------------------------------------------------------------*/
#define MAX_SNPS 10000

/*----------------------------------------------------------------------*
 * Global variables                                                     *
 *----------------------------------------------------------------------*/
char* input_filename = 0;
char* output_filename = 0;
char* snp_list_filename = 0;
char* output_id = 0;
int n_snps_to_insert = 1000;
int min_distance_between = 100;
char* reference = 0;
int reference_size = 0;
unsigned int snp_position[MAX_SNPS];
int n_snps = 0;
int column_width = 70;
char* default_id = "makesnps";

/*----------------------------------------------------------------------*
 * Function: is_nucleotide                                              *
 * Purpose:  Check if char is A, C, G or T                              *
 * Params:   c = character to check                                     *
 * Returns:  1 if nucleotide, 0 otherwise                               *
 *----------------------------------------------------------------------*/
int is_nucleotide(char c)
{
	char u = toupper(c);
	
	if ((u == 'A') || (u == 'C') || (u == 'G') || (u == 'T')) {
		return 1;
	}
	
	return 0;
}

/*----------------------------------------------------------------------*
 * Function: read_reference                                             *
 * Purpose:  Read the reference sequence                                *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void read_reference()
{
    FILE* fp;
	long file_size;
	int index = 0;
	char line[1024];
	char c;
	
	fp = fopen(input_filename, "r");
	if (!fp) {
		printf("Error: can't open reference.\n");
		exit(1);
	}
	
	fseek(fp, 0L, SEEK_END);
	file_size = ftell(fp);
	fseek(fp, 0L, SEEK_SET);
	
	reference = malloc(file_size);
	if (!reference) {
		printf("Error: can't get space to store reference.\n");
		fclose(fp);
		exit(1);
	}
	
	if (!fgets(line, 1024, fp)) {
		printf("Error: Couldn't get header line.\n");
		fclose(fp);
		exit(1);
	}
	
	if (line[0] != '>') {
		printf("Error: File should begin with FASTA header.\n");
		fclose(fp);
		exit(1);
	}
	
	while (!feof(fp)) {
		c = fgetc(fp);
		
		if (c == '>') {
			printf("Error: File should only have 1 sequence in it.\n");
			fclose(fp);
			exit(1);
		} else if (is_nucleotide(c)) {
			reference[index++] = toupper(c);
		}
	}
	
	reference_size = index;
	printf("Reference read... %d nucleotides.\n", reference_size);
	
	fclose(fp);
}

/*----------------------------------------------------------------------*
 * Function: positiokn_ok                                               *
 * Purpose:  Check if SNP position is ok (ie. not too close to another. *
 * Params:   p = proposed position                                      *
 * Returns:  1 = ok, 0 = not ok                                         *
 *----------------------------------------------------------------------*/
int position_ok(int p)
{
	int i;
	int ok = 0;
	
	if (n_snps == 0) {
		ok = 1;
	} else {	
		for (i=0; i<n_snps; i++) {
			if (snp_position[i] > p) {
				break;
			}
		}
		
		if (i == n_snps) {
			if ((p - snp_position[n_snps-1]) >= min_distance_between) {
				ok = 1;
			}
		} else if (i == 0) {
			if ((snp_position[0] - p) >= min_distance_between) {
				ok = 1;
			}
		} else {
			if (((snp_position[i] - p) >= min_distance_between) && ((p - snp_position[i-1]) >= min_distance_between)) {
				ok = 1;
			}
		}
	}
		
	return ok;
}

/*----------------------------------------------------------------------*
 * Function: compare                                                    *
 * Purpose:  Integer comparison function, used by qsort.                *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
int compare(const void* a, const void* b)
{
	return (*(int*)a - *(int*)b);
}

/*----------------------------------------------------------------------*
 * Function: make_snps                                                  *
 * Purpose:  Make SNPs                                                  *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void make_snps(void)
{
	int i, p;
	
	printf("Making SNPs...\n");
	
	for (i = 0; i<n_snps_to_insert; i++) {
		do {
			p = random() % reference_size;
		} while (!position_ok(p));
		
		snp_position[n_snps++] = p;	
		qsort(snp_position, n_snps, sizeof(int), compare);
	}
}

/*----------------------------------------------------------------------*
 * Function: display_snp_stats                                          *
 * Purpose:  Display SNP statistics                                     *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void display_snp_stats(void)
{
	int n_bins = 10;
	int bins[n_bins];
	int bin_size = reference_size / n_bins;
	int i, b;
		
	for (i=0; i<n_bins; i++) {
		bins[i] = 0;
	}
	
	for (i=0; i<n_snps; i++) {
		b = snp_position[i] / bin_size;
		bins[b]++;
	}
	
	printf("Distribution of SNPs:\n");
	for (i=0; i<n_bins; i++) {
		printf("Bin %d (%10d to %10d) \t: %d\n", i, i*bin_size, ((i+1)*bin_size)-1, bins[i]);
	}
	
	for (i=1; i<n_snps; i++) {
		int d = snp_position[i] - snp_position[i-1];
		if (d < min_distance_between) {
			printf("d: %i\n", d);
		}
	}
}

/*----------------------------------------------------------------------*
 * Function: make_snp                                                   *
 * Purpose:  Make SNP at position and write details to SNP position     *
 *           file.                                                      *
 * Params:   position = position of SNP                                 *
 *           current_n = reference nucleotide                           *
 *           csv_fp -> SNP position file handle                         *
 * Returns:  SNP nucleotide char                                        *
 *----------------------------------------------------------------------*/
char make_snp(int position, char current_n, FILE* csv_fp)
{
	char new_n;
	char nucleotides[] = {'A','C','G','T'};
	
	do {
		new_n = nucleotides[random() % 4];
	} while (new_n == current_n);
	
	fprintf(csv_fp, "%i,%c,%c\n", position, current_n, new_n);
	
	return new_n;
}

/*----------------------------------------------------------------------*
 * Function: write_output_files                                         *
 * Purpose:  Write output files                                         *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void write_output_files(void)
{
	FILE* out_fp;
	FILE* csv_fp;
	int p = 0;
	int c = 0;
	int snp_pos = 0;
	
	printf("Writing output files...\n");
	out_fp = fopen(output_filename, "w");
	if (!out_fp) {
		printf("Error: Can't open file %s.\n", output_filename);
		exit(1);
	}

	csv_fp = fopen(snp_list_filename, "w");
	if (!csv_fp) {
		printf("Error: Can't open file %s.\n", snp_list_filename);
		exit(1);
	}
	
	
	fprintf(out_fp, ">%s\n", output_id);
	fprintf(csv_fp, "Position,Reference,SNP\n");
	
	for (p=0; p<reference_size; p++) {
		int nucleotide = reference[p];
		
		if ((snp_pos < n_snps) && (p == snp_position[snp_pos])) {
			nucleotide = make_snp(snp_position[snp_pos], nucleotide, csv_fp);
			snp_pos++;
		}
		
		fputc(nucleotide, out_fp);

		c++;
		if (c == column_width) {
			fputc('\n', out_fp);
			c = 0;
		}
	}
	
	fclose(out_fp);	
	fclose(csv_fp);
}

/*----------------------------------------------------------------------*
 * Function: parse_string                                               *
 * Purpose:  Return a string parameter (command line parsing)           *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 *           i -> pointer to argument number counter                    *
 * Returns:  a string and updates i                                     *
 *----------------------------------------------------------------------*/
char* parse_string(int argc, char* argv[], int* i)
{
    char* token = 0;
    
    if (strlen(argv[*i]) > 2)
        token = argv[*i] + 2;
    else if (*i < (argc-1)) {
        *i = *i + 1;
        token = argv[*i];
    }
    
    if ((token) && (token[0] == '-'))
        token = 0;
    
    return token;
}

/*----------------------------------------------------------------------*
 * Function: parse_int                                                  *
 * Purpose:  Return an integer parameter (command line parsing)         *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 *           i -> pointer to argument number counter                    *
 * Returns:  a number and updates i                                     *
 *----------------------------------------------------------------------*/
int parse_int(int argc, char* argv[], int* i)
{
    int v = -1;
    char* token = parse_string(argc, argv, i);
    
    if (token)
        v = atoi(token);
    
    return v;
}

/*----------------------------------------------------------------------*
 * Function: parse_command_line_args                                    *
 * Purpose:  Deal with command line arguments                           *
 * Params:   argc = number of arguments                                 *
 *           argv -> argument array                                     *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void parse_command_line_args(int argc, char* argv[])
{
    int i = 1;
	
    if (argc < 4)
    {
        printf("Syntax: makesnps [-i filename] [-o filename] [-c filename] [options]\n");
        printf("where [-i filename] specifies the name of a reference genome in FASTA format.\n");
        printf("      [-o filename] specifies the name of an output FASTA file.\n");
		printf("      [-c filename] specifies the name of a CSV file to output containing SNP positions.\n");
		printf("      [-s id] specifies the output sequence id (default '%s').\n", default_id);
		printf("      [-n int] specifies the number of SNPs to insert (default %d).\n", n_snps_to_insert);
		printf("      [-m int] specifies the minimum distance between SNPs (default %d).\n", min_distance_between);
		printf("      [-w int] specifies the column width of the output file (default %d).\n", column_width);
        printf("\n");
        exit(1);
    }
	
    while(i < argc)
    {
        char* parameter = argv[i];
        if (parameter[0] == '-')
        {
            switch (parameter[1])
            {
                case 'c':
                    snp_list_filename = parse_string(argc, argv, &i);
                    break;
                case 'i':
                    input_filename = parse_string(argc, argv, &i);
                    break;
                case 'm':
                    min_distance_between = parse_int(argc, argv, &i);
                    break;
                case 'n':
                    n_snps_to_insert = parse_int(argc, argv, &i);
					if ((n_snps_to_insert < 1) || (n_snps_to_insert > MAX_SNPS)) {
						printf("Error: number of SNPs must be between 1 and %d.\n", MAX_SNPS);
						exit(1);
					}
                    break;
                case 'o':
                    output_filename = parse_string(argc, argv, &i);
                    break;
                case 's':
                    output_id = parse_string(argc, argv, &i);
                    break;
				case 'w':
					column_width = parse_int(argc, argv, &i);
					break;
                default:
                    printf("Error: Invalid parameter %c\n", parameter[1]);
                    exit(1);
            }
        }
        
        i++;
    }
    
    if (!input_filename) {
        printf("Error: You must specify an input file.\n");
        exit(1);
    }   

    if (!output_filename) {
        printf("Error: You must specify an output file.\n");
        exit(1);
    }
	
	if (!snp_list_filename) {
		printf("Error: You must specify a SNP list filename.\n");
	}
	
	if (!output_id) {
		default_id = output_id;
	}
}

/*----------------------------------------------------------------------*
 * Function: display_parameters                                         *
 * Purpose:  Display command line parameters                            *
 * Params:   None                                                       *
 * Returns:  None                                                       *
 *----------------------------------------------------------------------*/
void display_parameters(void) 
{
	printf("      Input filename: %s\n", input_filename);
	printf("     Output filename: %s\n", output_filename);
	printf("  Output sequence ID: %s\n", output_id);
	printf("   SNP list filename: %s\n", snp_list_filename);
	printf("      Number of SNPs: %d\n", n_snps_to_insert);
	printf("Min distance between: %d\n", min_distance_between);
	printf("        Column width: %d\n\n", column_width);
}

/*----------------------------------------------------------------------*
 * Program entry                                                        *
 *----------------------------------------------------------------------*/
int main (int argc, char* argv[])
{
	printf("\nmakesnps - create copy of genome with SNPs inserted.\n\n");
	parse_command_line_args(argc, argv);
	display_parameters();	
	time_t now = time(NULL);
	srandom(now);
	read_reference();
	make_snps();
	display_snp_stats();
	write_output_files();
	printf("Finished.\n");
	
	return 0;
}
