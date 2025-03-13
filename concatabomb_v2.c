#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <signal.h>

#ifndef TRUE
#define TRUE 1
#endif
#ifndef FALSE
#define FALSE 0
#endif
#ifndef TEMP
#define TEMP 2
#endif


void clean_exit(int x);
void clear_memory(void);
int comment(FILE *file);
int read_nexus(FILE *nexusfile);
int read_fasta(FILE *infile1);
int assign_arrays(void);
int read_sequences(FILE *infile1);
void removeChar(char *str, char garbage);


char **species_names = '\0', **sequences = '\0', *concatabomination ='\0',  characters[100];
int  numseqs = 0, aln_len =0, *incompatibilities='\0', *incompatibilities2='\0', calculate=TRUE, num_chars=0, chars_found=0, **char_matrix='\0', tot_incompats=0, *orig_assign='\0', current_incompats=0; 


/*
character matrix combinations:

character matches in the following order:

		0	1
	0	0	0	
	1	0	0 

Will contain a 1 if that character combination exists in the pairwise comparison.
in the arrray char_matrix, for each position of the alignment, we use a binary value to represent the what it contains. 
They are in the following order:

[0,0] [0,1] [1,0] [1,1]

in the above example the coding would be: 0000 which in binary is equal to (0x2^3)+(0x2^2)+(0X2^1)+(0X2^0) = 0

In the case of
	0	1
0	0	1	
1	1	0 

This would be equal to: 0110 whic in binary is equal to (0x2^3)+(1x2^2)+(1X2^1)+(0X2^0) = 6;

These values will be used to represent what already exists in the matrix. 

We are only interested in identifying the cases where there are already 3 1's in the matrix, they are:

1110 = (1x2^3)+(1x2^2)+(1X2^1)+(0X2^0) = 16
1101 = (1x2^3)+(1x2^2)+(0X2^1)+(1X2^0) = 13
1011 = (1x2^3)+(0x2^2)+(1X2^1)+(1X2^0) = 11 
0111 = (0x2^3)+(1x2^2)+(1X2^1)+(1X2^0) = 7

For all concatabominations, we only have to recalculate the incompatibilities for those sites that have binary values of 16, 13, 11 and 7.

We can calculate in advance the specific combinations that will make each of these cases incompatible (so we only need to know for the concatabominated sequences):
They are:

for values of: 
16: 1->1 will make is incompatible
13: 1->0 will make is incompatible
11: 0->1 will make is incompatible
7 : 0->0 will make is incompatible


LASTLY:

we don't have to look at combinations of characters on the concatabominated sequence, where both characters have come from the same inital sequence,
only those combinationa of characters where the charaters are originally from either of the concatabominated sequences


*/


int main(int argc, char *argv[])
    {
    FILE *infile1 = '\0', *outfile = '\0', *infile2='\0', *outfile2='\0';
    char filename[100], c = '\0', *taxa1 = '\0', *taxa2 = '\0';
    int i = 0, j = 0, k=0, l=0, compat[2][2], taxon1num=-1, taxon2num=-1, num=0, error=FALSE,  found=FALSE;
    if(argc < 3)
        {
        printf("\n\nconcatabomb version 2 \
		\n\nThis software takes as input a fasta formatted character matrix and a list of taxa to concatabominate\
		\nIt calculates the pairwise incompatibilities for the original matrix and then for all concatabominations specified\
		\n\
		\nThis currently is only implementd for matrices containing '0','1','?' and '-' \
		\nwhere '?' amd '-' represents missing data and is not considered when determining compatibility\
		\n\nUsage: concatabomb FastaFile ConcatabominationFile <nocalculate>\nwhere FastaFile contains a character matrix in fasta format\
		\nConcatabominationFile is a file containg the list of taxa to be concatabominated in turn\
		\nnocalculate is an option if specified will cause the software to only print the concatabominated nexus files to disk and not calculate compatibility\
		\nThis is to retain compatibility with the version of the pipline that uses compass (i.e. for multistate characters)\n\n");
		
        exit(1);
        }
    taxa1 = malloc(1000*sizeof(char));
    taxa2 = malloc(1000*sizeof(char));
	taxa1[0] = '\0';
	taxa2[0] = '\0';
	characters[0] = '\0';

	if(argc == 4)
		{
		if(strcmp(argv[3], "nocalculate") == 0)
			calculate=FALSE;
		}
		
	/* open file  */    
    infile1 = fopen(argv[1], "r");
    
	filename[0] = '\0';
	strcpy(filename, argv[1]);
	strcat(filename, ".abomscores");
	outfile2=fopen(filename, "w");

	error=read_sequences(infile1);
	if(error == TRUE)
		{
		clean_exit(432);
		}
    /*
	*
	*
	*/

	fclose(infile1);
	if(calculate && strcmp(characters, "01") != 0 && strcmp(characters, "10") != 0)
		{
		fprintf(stderr, "ERROR: Concatabomb can only calculate pairwise compatabilities for binary characters, Characters found are: %s\n", characters);
		clean_exit(-1);
		}
	

	
	/* Now calculate the pairwise compatibilities between all the characters of the original matrix */
	
	if(calculate)
		{
		for(i=0; i<aln_len; i++)
			{
			for(j=i+1; j<aln_len; j++)
				{
				if(i!=j)
					{
					compat[0][0] = compat[0][1] = compat[1][0] = compat[1][1] = 0;
					for(k=0; k<numseqs; k++)
						{
						if(sequences[k][i] != '?' && sequences[k][i] != '-' &&  sequences[k][j] != '?' && sequences[k][j] != '-')
							compat[sequences[k][i] - '0'][sequences[k][j] - '0'] = 1;
						}
					char_matrix[i][j] = char_matrix[j][i] = ((compat[0][0]*8) + (compat[0][1]*4) + (compat[1][0]*2) + (compat[1][1]*1));
					if(char_matrix[i][j] == 15)
						{	
						incompatibilities[i]++;
						incompatibilities[j]++;
						}
					}
				}
			}
		

	/* Print out the number of incompatibilities for each character for the original matrix */
		filename[0] = '\0';
		strcpy(filename, argv[1]);
		strcat(filename, ".IncompatCounts");
		outfile=fopen(filename, "w");
		
		for(i=0; i<aln_len; i++)
			{
			tot_incompats+=incompatibilities[i];
			fprintf(outfile, "%d\n", incompatibilities[i]);	
			}
		
		fclose(outfile);
		}	
	
/* Now read in the list of species to concatabominate and for each calculate if they add any more incompatibilities */

	infile2=fopen(argv[2], "r");
	if(!infile2)
		{
		fprintf(stderr, "ERROR: File named: %s does not exist\n", argv[2]);
		error=TRUE;
		}
		
	if(!error)
		{
		num=1;
		while(!feof(infile2))
			{
			/*read in the names of the taxa to be concatabominated */
			fscanf(infile2, "%s\t%s\n", taxa1, taxa2);
			taxon1num=taxon2num=-1;
			for(i=0; i<numseqs; i++)
				{
				if(strcmp(taxa1, species_names[i]) == 0)
					taxon1num=i;
				if(strcmp(taxa2, species_names[i]) == 0)
					taxon2num=i;
				}
			if(taxon1num == -1)
				{
				fprintf(stderr,"Error taxon name %s on line %d in concatabomination file does not exist in the input matrix\n", taxa1, num);
				error=TRUE;
				}
			if(taxon2num == -1)
				{
				fprintf(stderr,"Error taxon name %s on line %d in concatabomination file does not exist in the input matrix\n", taxa2, num);
				error=TRUE;
				}
			if(taxon2num == taxon1num && taxon1num != -1)
				{
				fprintf(stderr,"Error the same taxa (%s) was specified twice on line %d\n", taxa2, num);
				error=TRUE;
				}			
			if(error)
				fprintf(stderr,"Skipping concatabomination specified in line %d\n", num);
		
			if(!error)
				{
				/* do the concatabomination of the two taxa */
				for(i=0; i<aln_len; i++)
					{
					if(sequences[taxon1num][i] == '?' || sequences[taxon1num][i] == '-' )
						{
						if(sequences[taxon2num][i] == '?' || sequences[taxon2num][i] == '-' )
							{
							orig_assign[i] = 3;
							concatabomination[i] = '?';
							}
						else
							{
							orig_assign[i] = 2;
							concatabomination[i] = sequences[taxon2num][i];
							}
						}
					else
						{
						if(sequences[taxon2num][i] == '?' || sequences[taxon2num][i] == '-' )
							{
							orig_assign[i] = 1;
							concatabomination[i] = sequences[taxon1num][i];
							}
						else
							{
							if(sequences[taxon1num][i] != sequences[taxon2num][i])
								{
								fprintf(stderr,"ERROR: A character conflict exists between taxa %s and %s at character %d (%c versus %c)\n", taxa1, taxa2, i, sequences[taxon1num][i], sequences[taxon2num][i]);
								fprintf(stderr,"Skipping this taxon combination and continuing....\n");
								error=TRUE;
								}
							else
								{
								orig_assign[i] = 3;
								concatabomination[i] = sequences[taxon1num][i];
								}
							}
						}
					if(error) i=aln_len;
					}
				if(calculate) /* if we are to calculate the compatibility scores */
					{
					if(!error) /* now calculate if there are any new conflicts created by this concatabomination */
						{
						for(i=0; i<aln_len; i++) incompatibilities2[i]=incompatibilities[i]; /* Start by assigning the number ob incompatitibilities for this character the same as calculated from the original matrix */

						for(i=0; i<aln_len; i++)
							{
							if(orig_assign[i] != 3) /* anything with a three means that it was the same character in both taxa before concataominating so the calculation for the iginal matrix is correct */
								{
								for(j=i+1; j<aln_len; j++)
									{
									if(orig_assign[j] != 3) /* anything with a three means that it was the same character in both taxa before concataominating so the calculation for the iginal matrix is correct */
										{
										if(orig_assign[i] != orig_assign[j]) /* If the have the same number (1 or 2) then both sites are from the same original taxa and the calculation for the original matrix is correct */
											{
											switch(char_matrix[i][j]) 
												{
												case 14: /* 1->1 will make it incompatible */
													if(concatabomination[i] == '1' && concatabomination[j] == '1') /* This makes these characters incompatible */
														{
														incompatibilities2[i]++;
														incompatibilities2[j]++;
														}
													break;
												case 13: /* 1->0 will make it incompatible */
													if(concatabomination[i] == '1' && concatabomination[j] == '0') /* This makes these characters incompatible */
														{
														incompatibilities2[i]++;
														incompatibilities2[j]++;
														}						
													break;
												case 11: /* 0->1 will make it incompatible */
													if(concatabomination[i] == '0' && concatabomination[j] == '1') /* This makes these characters incompatible */
														{
														incompatibilities2[i]++;
														incompatibilities2[j]++;
														}
													break;
												case 7: /* 0->0 will make it incompatible */
													if(concatabomination[i] == '0' && concatabomination[j] == '0') /* This makes these characters incompatible */
														{
														incompatibilities2[i]++;
														incompatibilities2[j]++;
														}
													break;
												default : /* Nothing will make this incompatible so don't calculate */
													break;
												}
											}
										}
									}
								}
							}
						
						/* Print out incompatibility scores for this taxa concatibomination */
						/*filename[0] = '\0';
						strcpy(filename, argv[1]);
						strcat(filename, ".");
						strcat(filename, taxa1);
						strcat(filename, ".");
						strcat(filename, taxa2);
						strcat(filename, ".IncompatCounts");
						outfile=fopen(filename, "w"); */
						current_incompats=0;
						for(i=0; i<aln_len; i++)
							{
							current_incompats+=incompatibilities2[i];
							/*fprintf(outfile, "%d\n", incompatibilities2[i]);	*/
							}
						
						/*fclose(outfile);*/					
						
						
						
						fprintf(outfile2, "%s&%s %d %d\n", taxa1, taxa2, tot_incompats, current_incompats);
						fprintf(stdout, "%s&%s %d %d\n", taxa1, taxa2, tot_incompats, current_incompats);
						}
					else
						{
						error=FALSE; /* reset this for the next combination of taxa to concatabominate */
						}
					}
				else	/* Just print the concatabominated version of the nexus file to disk */
					{
					filename[0]='\0';
					strcpy(filename, argv[1]);
					strcat(filename, ".abom.");
					strcat(filename, taxa1);
					strcat(filename, "-");
					strcat(filename, taxa2);
					outfile=fopen(filename, "w");
					fprintf(outfile, "#NEXUS\n[Baum-Ragan coding ]\n\n\nBegin data;\nDimensions ntax=%d  nchar=%d;\nFormat datatype=standard missing=? gap=- symbols=\"",numseqs-1, aln_len);
					for(i=0; i<chars_found; i++)
						fprintf(outfile, "%c", characters[i]);
					fprintf(outfile, "\";\nMatrix\n");
					
					for(i=0; i<numseqs; i++)
						{
						if(i != taxon1num && i != taxon2num)
							{
							fprintf(outfile,"\'%s\' ", species_names[i]);
							for(j=0; j<aln_len; j++)
								fprintf(outfile, "%c", sequences[i][j]);
							fprintf(outfile, "\n");
							}
						}
					fprintf(outfile,"\'%s&%s\' ", taxa1, taxa2);
					for(j=0; j<aln_len; j++)
						{
						fprintf(outfile, "%c", concatabomination[j]);
						}
					fprintf(outfile, "\n  ;\nEnd;\n");
					
					fclose(outfile);

					}
				}
			num++;
			}
		fclose(infile2);
		}
	if(taxa1 != '\0') free(taxa1);
	if(taxa2 != '\0') free(taxa2);

	clear_memory();
    fclose(outfile2);
    }
    



void clean_exit(int x)
    {
    if(x > 0) fprintf(stderr,"Error: out of memory at %d\n", x); 
	clear_memory();
    exit(0);
    }
	
void clear_memory(void)
	{
	int i, j;
	

	if(species_names != '\0')
		{
		for(i=0; i<numseqs; i++)
			{
			if(species_names[i] != '\0')
				{
				free(species_names[i]);
				species_names[i] = '\0';
				}
			}
		free(species_names);
		species_names = '\0';
		}


	if(sequences != '\0')
		{
		for(j=0; j<numseqs; j++)
			{
			if(sequences[j] != '\0')
				{
				free(sequences[j]);
				sequences[j] = '\0';
				}
			}
		free(sequences);
		sequences = '\0';
		}
	if(incompatibilities != '\0')
		{
		free(incompatibilities);
		incompatibilities = '\0';
		}
	if(incompatibilities2 != '\0')
		{
		free(incompatibilities2);
		incompatibilities2 = '\0';
		}
		
		
	if(char_matrix != '\0')
		{
		for(i=0; i<aln_len; i++)
			{
			if(char_matrix[i] != '\0')
				{
				free(char_matrix[i]);
				char_matrix[i] = '\0';
				}
			}
			free(char_matrix);
			char_matrix = '\0';
		}
		
	if(concatabomination != '\0')
		{
		free(concatabomination);
		concatabomination='\0';
		}
	}

int read_sequences(FILE *infile1)
	{
	int i, j, error=FALSE;
	char c;
	
	/* read through the input file and determine the number of taxa and the length of the character matrix */
	aln_len=0; numseqs=0;
	
	/* First Determine if this is a fasta or Nexus file */
	
	while((c=getc(infile1)) == ' ' || c == '\n' || c == '\t' || c == '\r' );
	switch(c)
		{
		case '#':
			/* nexus format */
			rewind(infile1);
			error= read_nexus(infile1);
			break;
		case '>':
			/* Fasta Format */
			rewind(infile1);
			error = read_fasta(infile1);
			
			break;
		default:
			error=TRUE;
			printf("ERROR: Unknown input file type\n");
			break;
		}

	/* Check data read in */
	/*printf("numtaxa = %d\tnumchars=%d\n", numseqs, aln_len);
	for(i=0; i<numseqs; i++)
		{
		printf ("%s\t%s\n", species_names[i], sequences[i]);
		}
	*/
	return(error);
	}


int assign_arrays(void)
	{
	int i, j, k, error=FALSE;
	
	/* assign the arry to hold the names of the species  */
	species_names = malloc(numseqs*sizeof(char*));
    if(species_names == '\0') clean_exit(1);
	
	for(i=0; i<numseqs; i++)
		{
		species_names[i] = malloc(1000*sizeof(char));
		if(!species_names[i]) clean_exit(23);
		species_names[i][0] = '\0';
		}

	/* Assign the array to hold all the sequences of characters per taxon */
	sequences = malloc(numseqs*sizeof(char*));
    if(sequences == '\0') clean_exit(2);

	for(j=0; j<numseqs; j++)
		{
		sequences[j] = malloc(aln_len+1*sizeof(char));
		if(sequences[j] == '\0') clean_exit(3);
		for(k=0; k<aln_len+1; k++) 
			{
			sequences[j][0]= '\0';
			}
		}

	concatabomination=malloc(aln_len*sizeof(char));
	if(!concatabomination) clean_exit(111);
	for(i=0; i<aln_len; i++)
		{
		concatabomination[i] = '\0';
		}

	orig_assign=malloc(aln_len*sizeof(int));
	for(i=0; i<aln_len; i++)
		{
		orig_assign[i] = 0;
		}
	if(!orig_assign) clean_exit(112);
	

	if(calculate)
		{
		/* this records the number of incompatibilites with each character */
		incompatibilities=malloc(aln_len*sizeof(int));
		if(!incompatibilities) clean_exit(47);
		for(i=0; i<aln_len; i++)
			incompatibilities[i]=0;

	/* this records the number of incompatibilites with each character */
		incompatibilities2=malloc(aln_len*sizeof(int));
		if(!incompatibilities2) clean_exit(48);
		for(i=0; i<aln_len; i++)
			incompatibilities2[i]=0;		
			
		/* char_matrix records in binary format the character combinations for each pairwise comparison of characters */
		
		char_matrix=malloc(aln_len*sizeof(int*));
		if(char_matrix == '\0') clean_exit(99);
		for(i=0; i<aln_len; i++)
			{
			char_matrix[i]=malloc(aln_len*sizeof(int));
			if(char_matrix[i] == '\0') clean_exit(101);
				for(j=0; j<aln_len; j++)
				{
				char_matrix[i][j] = 0;
				}
			}
			
		}
		
	/*printf("Done assignment\n");*/
	return(error);
	}
	
	
	
int read_fasta(FILE *infile1)
	{
	int i, j,k, error=FALSE, found=FALSE;
	char c;
	
	while(!feof(infile1))
		{
		c = getc(infile1);	
		if(c == '>' ) 
			{ 
			numseqs++;
			if( aln_len == 0)
				{
				while(!feof(infile1) && (c=getc(infile1)) != '\n'); /* go to the end of the name */
				while(!feof(infile1) && (c=getc(infile1)) != '>')
					{
					if(c != ' ' && c != '\t' && c != '\n' && c != '\r') aln_len++;
					}
				if(c == '>') numseqs++;
				}
			}
		}
					
	rewind(infile1);
	
	/*printf("number of taxa = %d\nnumber of characters = %d\n", numseqs, aln_len);*/
	
	
	error=assign_arrays();
	if(error == TRUE)
		{
		printf("Error assigning arrays, maybe out of memory\n");
		clean_exit(555);
		}
		
	/**** READ IN THE SEQUENCES  ***/
    i=-1; j=0;
	while(!feof(infile1))
		{
		if((c = getc(infile1)) == '>' ) 
			{
			if(i>=0)
				{
				sequences[i][j] = '\0';
				}
			i++;
			j=0;
			while(!feof(infile1) && (c=getc(infile1)) != '\n' && c != '\r')
				{
				if(c != '\n' && c != '\r') 
					{
					species_names[i][j] = c;
					j++;
					}
				}
			species_names[i][j] = '\0';
			j=0;

			}
		else
			{
			if(!feof(infile1) && c != '>' && c != ' ' && c != '\t' && c != '\n' && c != '\r' && c != '\t' )
				{
				sequences[i][j] =c;
				if(!feof(infile1) && c != '?' && c != '-')
					{

					found=FALSE;
					for(k=0; k<chars_found; k++)
						{
						if(c==characters[k]) found=TRUE;
						}
					if(!found)
						{
						characters[chars_found]=c;
						characters[chars_found+1]='\0';
						chars_found++;
						}
					found=FALSE;
					}
				
				j++;
				}
			}
		}
	sequences[i][j] = '\0';
	return(error);
	}
	
int read_nexus(FILE *nexusfile)
	{
	char c, begin[6] = {'b', 'e', 'g', 'i', 'n', '\0'}, tree[5] = {'t','r','e','e','\0'};
	int error = FALSE, i, j, k, l, found = FALSE;
	char *string = '\0', *token ='\0', *res='\0', *string1 = '\0';
	
	string=malloc(1000*sizeof(char));
	string1=malloc(1000*sizeof(char));
	string[0] = '\0';
	string1[0] = '\0';
	while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r' ) && !feof(nexusfile));

	while(!feof(nexusfile) && !error)
		{
		while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r' ) && !feof(nexusfile));
		switch(tolower(c))
			{
			case '[':
				
				comment(nexusfile);
				break;
			
			case 'b':   /* this is the beginning of a block **/
				for(i=1; i<5; i++)
					if((c = tolower(getc(nexusfile))) != begin[i]) error = TRUE;
				if(!error)
					{
					while(((c = getc(nexusfile)) == ' ' || c == '\t' || c == '\n' || c == '\r' ) && !feof(nexusfile));
					i = 0;
					do {	/* find the type of block it is **/
						string[i] = tolower(c);
						i++;
						c = getc(nexusfile);
						}while(c != ';' && c != ' ');
					string[i] = '\0';

					if(strcmp(string, "data") == 0)
						{
						while(strcmp(string, "end;") != 0 && strcmp(string, "endblock;") != 0 && !feof(nexusfile) && !error)
							{
							fscanf(nexusfile, "%s", string);
							for(i=0; i<strlen(string); i++) string[i]=tolower(string[i]);
							if(strstr(string, "ntax") != '\0') /* Extract the number of taxa */
								{
								token='\0';
								token = strtok(string, "=");
								token = strtok(NULL, "=");
								strcpy(string1, token);
								removeChar(string1, ';');
								numseqs=atoi(string1);
								}
							if(strstr(string, "nchar") != '\0') /* Extract the number of characters */
								{
								token='\0';
								token = strtok(string, "=");
								token = strtok(NULL, "=");
								strcpy(string1, token);
								removeChar(string1, ';');
								aln_len=atoi(string1);
								}

							if(strstr(string, "matrix") != '\0') /* Start reading the character matrix */
								{
								if(aln_len == 0 || numseqs == 0)
									{
									error = TRUE;
									printf("Error: numtaxa or numcharacters is not defined properly\n");
									}
								else
									{
									assign_arrays(); /* assign the arrrays */
									i=-1; j=0;

									while(!feof(nexusfile) && (c=getc(nexusfile)) != ';')
										{
										/* Next read in the matrix */
										if(c=='\'') /*read in the name */
											{
											i++;
											j=0;
											while((c=getc(nexusfile)) != '\'')
												{
												if(c != ' ' && c != '\t') 
													{
													species_names[i][j] = c;
													j++;
													}
												}
											species_names[i][j] = '\0';
											j=0;
											while((c=getc(nexusfile)) != '\n' && c != '\r') /* read in the characters */
												{
												if(c != ' ' && c != '\t') 
													{
													sequences[i][j] =c;
													j++;
													}
												
												if(!feof(nexusfile) && c != '?' && c != '-' && c!=' ' && c !='\t') /* record the characters */
													{

													found=FALSE;
													for(k=0; k<chars_found; k++)
														{
														if(c==characters[k]) found=TRUE;
														}
													if(!found)
														{
														characters[chars_found]=c;
														characters[chars_found+1]='\0';
														chars_found++;
														}
													
													}
												}
											sequences[i][j] = '\0';
											j=0;
											}
										}
									}
								}
							}
						}
					else
						{
						
							/*** skip to the end of the block -- it must contain trees or commands or something else **/
							strcpy(string, "");
							while(strcmp(string, "end;") != 0 && strcmp(string, "endblock;") != 0 && !feof(nexusfile) && !error)
								{
								i=0;
								while((c == ' ' || c == '\n' || c == '\r' || c == '\t' || c == ';') && !feof(nexusfile)) c = getc(nexusfile);
								while(c != ';' && c != ' ' && c != '\n' && c != '\r' && c != '\t' && !feof(nexusfile))
									{
									string[i] = tolower(c);
									if(i < 9997)i++;
									c = tolower(getc(nexusfile));
									}
								string[i] = ';';
								string[i+1] = '\0';
								}
							
						}	
					}
				else
					{
					/*** what to do if it doesn't say begin? ***/
					error = TRUE;
					}
				break;
				
			default:
				break;


			}
		}
	
/*	fclose(nexusfile); */
	free(string);
	return(error);
	}



int comment(FILE *file)
	{
	char c;
	c = getc(file);
	while(c != ']' && !feof(file))
		{
		if(c == '[') comment(file);
		c = getc(file);
		}
	if(!feof(file)) return(FALSE);
	else return(TRUE);
	}


void removeChar(char *str, char garbage) {

    char *src, *dst;
    for (src = dst = str; *src != '\0'; src++) {
        *dst = *src;
        if (*dst != garbage) dst++;
    }
    *dst = '\0';
}



	
