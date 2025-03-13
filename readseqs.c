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

void clear_memory(void);
void clean_exit(int x);
int comment(FILE *file);
int read_nexus(FILE *nexusfile);
int read_fasta(FILE *infile1);
int assign_arrays(void);
int read_sequences(FILE *infile1);
void removeChar(char *str, char garbage);


char **species_names = '\0', **sequences = '\0', *concatabomination ='\0',  characters[100];
int  numseqs = 0, aln_len =0, *incompatibilities='\0', *incompatibilities2='\0', calculate=TRUE, num_chars=0,  **char_matrix='\0', tot_incompats=0, *orig_assign='\0', current_incompats=0; 


int main(int argc, char *argv[])
    {
    FILE *infile1 ='\0';
    

    /* open file  */    
    infile1 = fopen(argv[1], "r");
    read_sequences(infile1);

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
	printf("numtaxa = %d\tnumchars=%d\n", numseqs, aln_len);
	for(i=0; i<numseqs; i++)
		{
		printf ("%s\t%s\n", species_names[i], sequences[i]);
		}

	clear_memory();	

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
					for(k=0; k<num_chars; k++)
						{
						if(c==characters[k]) found=TRUE;
						}
					if(!found)
						{
						characters[num_chars]=c;
						characters[num_chars+1]='\0';
						num_chars++;
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
	int num_taxa =0, num_chars =0;
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
							if(strstr(string, "ntax") != '\0') /* Extract the number of taxa */
								{
								token='\0';
								token = strtok(string, "=");
								token = strtok(NULL, "=");
								strcpy(string1, token);
								removeChar(string1, ';');
								num_taxa=atoi(string1);
								}
							if(strstr(string, "nchar") != '\0') /* Extract the number of characters */
								{
								token='\0';
								token = strtok(string, "=");
								token = strtok(NULL, "=");
								strcpy(string1, token);
								removeChar(string1, ';');
								num_chars=atoi(string1);
								}

							if(strstr(string, "matrix") != '\0') /* Start reading the character matrix */
								{
								if(num_chars == 0 || num_taxa == 0)
									{
									error = TRUE;
									printf("Error: numtaxa or numcharacters is not defined properly\n");
									}
								else
									{
									aln_len=num_chars;
									numseqs=num_taxa;
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
												if(!feof(nexusfile) && c != '?' && c != '-') /* record the characters */
													{
													found=FALSE;
													for(k=0; k<num_chars; k++)
														{
														if(c==characters[k]) found=TRUE;
														}
													if(!found)
														{
														characters[num_chars]=c;
														characters[num_chars+1]='\0';
														num_chars++;
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


void removeChar(char *str, char garbage) {

    char *src, *dst;
    for (src = dst = str; *src != '\0'; src++) {
        *dst = *src;
        if (*dst != garbage) dst++;
    }
    *dst = '\0';
}





