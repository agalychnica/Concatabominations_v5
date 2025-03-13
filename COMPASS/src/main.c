/******************/
/*Included headers*/
/******************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <ctype.h>

/*********/
/*Defines*/
/*********/

#define Maxtax 500
//#define Maxchar 15000
#define Maxchar 3500000
#define Maxstates 30
#define Maxiterations 10000


/*********************/
/*Function prototypes*/
/*********************/

int  doExecuteFile(void);
int tests(void);
//int getfiles(void);
int Matread(void);
void Doreadcompassblock(long spot);
//void bootstrap(void);
//void boildown(void);
//void CCSR(void);
//void NDEV(void);
void LeQuesne(void);
void iterations(void);
void LeQuesnefast(void);
void iterationsfast(void);
void comparetype (int i, int h, int iteration);
void comparetypefast (int i, int h, int iteration);
void compare(int i, int h, int j, int iteration, int fast);
void randomise(int i);
void ressler(int i, char tempchar[Maxtax], char newchar[], int tempcount);
int countstates(char character[Maxtax]);
void resetmatrix(int chr1, int chr2, int h, int iteration, int fast);
void delnoncorners(void);
int delloop(int chr1, int chr2, int h);
void ecs(int start, int x, int y, int statechar1, int statechar2);
void why(int start, int x, int y, int statechar1, int statechar2);
//void printcompat(void);
void boildown(void);
void boildownfast(void);
void delchar(int worst, int p, int uninf);
void delcharfast(int worst);
//void deltax(int tax);
int finduninf(void);
int removeuninf(void);
void Removebrakets(void);
void CalculateTestValues (void);
void DoTestOutput (void);
void Domask(char name[50], char cur);
void Docharset(char name[50], char cur);
void writeplot(char name[5]);
void writematrix(char name[5]);
int ran(int k);
void Flush( void );
void DoError(char *message, char usage);
void Usage( void );
void readdata(int number);
void DoInteractive(int executeok);
void printoptions(void);
void getfilename(int executeok);
void printheading(void);
//void Clique (void);
//void Toss (int i, /*int nic, int cs[Maxchar][Maxchar],*/ int ce[Maxchar], int ne[Maxchar], int clique[Maxchar]);

//char *CFStringToCString( CFStringRef input );






/******************/
/*Global variables*/
/******************/

FILE* ftemp;
FILE* ftempb;
FILE* fout;
FILE* flog;
FILE* maskfile;
FILE* charsetfile;
FILE* fdata;


int gtaxa =0;
char taxonnames[Maxtax][100];
int gchars = 0;
int alignlength = 0;
int gnostates;
char missing = '?';
char gap = '?';
int nouninf;
int noequiv;
int noconst;
char gstate[Maxstates];
//char command = 'y';
char inputfile[500];
char outputfile[520];
int *charnum;
int *orignum;
int taxnum[Maxtax];
char comparison[Maxtax];
//float LQP[Maxchar];
char **gMatrix=NULL;
char **tempMatrix=NULL;
char *keep;
int ***compat;
float taxaremoved[Maxtax][2];
int *actualincompat;
int permutationincompat;
float *Teststatistic;
int *tot;
float *tot2;
int *greater;
int *type;
int *same;

int tns;
int xmax;
int ymax;
int xaxis[Maxstates];
int yaxis[Maxstates];
int d[Maxstates][Maxstates];

int *removedlastround;
int noremovedlast=0;
int *included;
int gfrom;
int gto;
//int iterationincompats[Maxchar][Maxiterations];
int **iterationincompats;
char string[100];
char mask='n';
char charset='n';
int nodatafiles=0;
char datafile[1000][100];
float cutoff=0;
float gstopvalue=0;

//int nic=0; //number of chars in clique

//char boiltype;

char loop;
time_t startseed;
time_t charseed[Maxchar];


int currentnumber=0;
int currentchar=0;

//char uninffound = 'y';

//char gAnalysisdone = 'n';
char gCompattype = 'C';
//char gJactaxtype = 'A';
int gPermutations = 0;
//int gjacPermutations = 0;
int gJacpermtype = 0;
int gFuzzy = 0;
int gMissingconstant = 1;
int gPrintUninf = 1;
int gBoildown = 0;
//char gBoot='n';
//char executedyn = 'n';
//char resultsyn='n';
//int gFast = 0;
//int gRemove = 0;
float gCharset =0;
char interactive='n';

//char gSearchtype='H';
//char gAdditiontype='S';
//int gBootreps=100;
//int gHeuristicrands=10;
//int gSeed=0;
//int gLE=0;
//char IncMaxtrees='n';
//int gMaxtrees;
//char PAUPinfile[50];

char ftempopen = 'n';
char newexecutionyn = 'n';




int  doExecuteFile(void)
{

    int readok = 1;

    if (ftempopen=='y'){fclose(ftemp);}

    if (interactive=='n')printf("\nExecuting File \"%s\"...\n\n", inputfile);



    if ((ftempb=fopen(inputfile, "r"))==NULL){
        printf("Cannot find file \"%s\"", inputfile);
        printf("\n...Execution terminated without completion\n");
//        executedyn='n';
        return readok;
//        printf("%s", string);
//        DoError(string, 'n');
    }
    fseek (ftempb, 0, SEEK_SET);
    
    if (gto!=0){
        sprintf(outputfile, "%s.%d.%d.log", inputfile, gfrom+1, gto);
    }
    else {
        sprintf(outputfile, "%s.%d.end.log", inputfile, gfrom+1);
    }
    flog=fopen(outputfile, "w");    
    if (mask=='y'){
        if (gto>=0){
            sprintf(outputfile, "%s.%d.%d.msk", inputfile, gfrom+1, gto);
        }
        else {
            sprintf(outputfile, "%s.%d.end.msk", inputfile, gfrom+1);
        }
        maskfile=fopen(outputfile, "w");
    }
    
    if (charset=='y'){
        if (gto>=0){
            sprintf(outputfile, "%s.%d.%d.charset", inputfile, gfrom+1, gto);
        }
        else {
            sprintf(outputfile, "%s.%d.end.charset", inputfile, gfrom+1);
        }
        charsetfile=fopen(outputfile, "w");
    }
    
    fprintf(flog,"\nCOMPASS (Compatibility Analysis Site Stripping)\n\nv1.0 (2009) written by Simon Harris, The Wellcome Trust Sanger Institute, Hinxton, UK\n\n");
    fprintf(flog,"\nExecuting File \"%s\"...\n\n", inputfile);
    newexecutionyn = 'y';
    ftemp = tmpfile();
    Removebrakets();
    fclose(ftempb);
    ftempopen='y';
    readok = Matread();
    //uninffound='n';

    if (readok==0){
        if (interactive=='n')printf("\n...File successfully executed\n\n");
        fprintf(flog,"\n...File successfully executed\n");
        readok=0;
    }
    else{
        printf("\n...Execution terminated without completion\n\n");
        fprintf(flog,"\n...Execution terminated without completion\n");
//        executedyn='n';
        return readok;
    }
//    resultsyn='n';


    return readok;
}


int tests( void )
{

    int p;
    int i;
    int noconst=0;
    
    
//Initialisation
for (p=0; p<gchars; p++)
{	//Stores character numbers in charnum field
        charnum[p] = p+1;
        included[p]=1;
        for (i=0; i<gPermutations; i++){
            iterationincompats[p][i]=0;
        }
}
charnum[gchars+1]='\0';


for (p=0; p<gtaxa; p++)
{	/*Stores taxon numbers in taxnum field*/
taxnum[p] = p+1;
}
taxnum[gtaxa+1]='\0';

if (gFuzzy==1){fprintf(flog,"\n\nFuzzy ");printf("Fuzzy ");}
else if (gFuzzy==0){fprintf(flog,"\n\n");}

if (gCompattype=='C'){fprintf(flog,"CCSR");printf("CCSR");}
else if (gCompattype=='N'){fprintf(flog,"NDev");printf("NDev");}
else if (gCompattype=='L'){fprintf(flog,"LQP");printf("LQP");}
else if (gCompattype=='A'){fprintf(flog,"Number of Incompatibilities");printf("Number of Incompatibilities");}


if ((gCompattype=='C')||(gCompattype=='N')||(gCompattype=='L')) 
{fprintf(flog," analysis of \"%s\":\n\t%d random permutations", inputfile, gPermutations);
fprintf(flog,"\n\tstarting seed = %ld", startseed);
printf(" analysis of \"%s\":\n\t%d random permutations", inputfile, gPermutations);
printf("\n\tstarting seed = %ld", startseed);


    if (gMissingconstant == 1){
        fprintf(flog,"\n\tMissing data points held constant");printf("\n\tMissing data points held constant");}
    else if (gMissingconstant == 0){
        fprintf(flog,"\n\tMissing data points not held constant");printf("\n\tMissing data points not held constant");}
}
else if (gCompattype=='A'){ fprintf(flog," in matrix: \"%s\"", inputfile);}

//if (j==0){ffprintf(flog,fout, "\n\nCompatibility matrix of \"%s\":\n", inputfile);}
if (gCompattype!='Q'){
    noconst = finduninf();
    
    
    if (mask=='y'){
        sprintf(string, "Constant");
        Domask(string, 'c');
        sprintf(string, "Uninformative");
        Domask(string, 'u');
    }
    if (charset=='y'){
        sprintf(string, "charset Constant =");
        Docharset(string, 'c');
        sprintf(string, "charset Uninformative =");
        Docharset(string, 'u');
    }


    fprintf(flog,"\n\t%d characters are parsimony uninformative", nouninf);
    fprintf(flog,"\n\t%d of which are constant", noconst);
    fprintf(flog,"\n\t%d characters are equivalent to others\n", noequiv);
    
    printf("\n\t%d characters are parsimony uninformative", nouninf);
    printf("\n\t%d of which are constant", noconst);
    printf("\n\t%d characters are equivalent to others\n", noequiv);
   // uninffound = 'y';
}

if (gto-gfrom<gchars){
    fprintf(flog,"\tAnalysing characters %d to %d (%d characters)\n", gfrom+1, gto, gto-gfrom);
    printf("\tAnalysing characters %d to %d (%d characters)\n\n", gfrom+1, gto, gto-gfrom);
}
else {
    fprintf(flog,"\tAnalysing %d characters\n", gchars-(noequiv+nouninf));
    printf("\tAnalysing %d characters\n\n", gchars-(noequiv+nouninf));
}

	
for (p=0;p<gchars;p++){
       charseed[p]=startseed+charnum[p];
              }
//fout=fopen("test.out", "w");
if (nodatafiles==0){
    LeQuesne();
}
else {
    for (p=0;p<nodatafiles;p++){
        readdata(p);
    }
    //writematrix();
}



if (gCompattype!='Q'){
    CalculateTestValues();
    
    writematrix("out");
    
    writeplot("");


    if ((gBoildown==1)||((gBoildown==1)&&(gCompattype!='L')))
    {
        DoTestOutput();
        
        boildown();
    	DoTestOutput();
        writeplot("_boil");
    }
    else if (gBoildown==2)
    {
        DoTestOutput();
        boildownfast();
    	DoTestOutput();
        writeplot("_boil");
    }
    else DoTestOutput();
}




/*if ((gBoot=='y')&&(exittest=='n')){
    printf("\n\nend;");
    fprintf(flog,"\n\nThe file \"%s\" should now be executed in PAUP.\nNOTE: This file MUST be in the same folder as your data file", PAUPinfile);
    fprintf(flog," (\"%s\"). Results from the PAUP analysis will be saved in the same folder in a tree file named \"%s.tre\"", inputfile, PAUPinfile);
    fprintf(flog,". This tree file should then be read back into Boildown using the Analyses->Boildown Bootstrap->Read PAUP Output command.");
}*/

return 1;
}


int Matread(void){
    /*Read nexus input file*/

    int a;
    int i;
    int j;
    char yn= 'y';
    char temp[100];
    //char tempb[100];
    long spot=0;
    long oldspot;
    char found[5] = "nnnnn";
    int readok=0;
    char datatype='s';
    int gcharsnew;

    for (i=0; i<20; i++){
        gstate[i] = '\0';
    }

    fscanf(ftemp, " %100s", temp);


    if (strncasecmp(temp, "#nexus", 6)!=0){
        readok=1;
        DoError("Not a Valid NEXUS File. Expecting #NEXUS at start of file", 'n');
        return readok;
    }



    while (found[3]=='n'){

        oldspot = spot;
        spot = ftell(ftemp);

        if (spot!=oldspot){

            fscanf(ftemp, " %100s", temp);

            if (strncasecmp(temp, "ntax", 4)==0){
                fseek (ftemp, spot, SEEK_SET);
                fscanf (ftemp, " %[^=]=", &temp[0]);
                fscanf (ftemp, " %d", &gtaxa);

                if((gtaxa>Maxtax) || (gtaxa<3)){
                    readok=1;
                    sprintf(string,"Invalid number of taxa defined in input file. Must be between 3 and %d", Maxtax);
                    DoError(string, 'n');
                    return readok;
                }
                found[0] = 'y';
                fprintf(flog,"\tNumber of taxa = %d\n", gtaxa);
                if (interactive=='n')printf("\tNumber of taxa = %d\n", gtaxa);
            }

            else if (strncasecmp(temp, "nchar", 5)==0){
                fseek (ftemp, spot, SEEK_SET);
                fscanf (ftemp, " %[^=]=", &temp[0]);
                fscanf (ftemp, " %d", &gchars);
                if((gchars>Maxchar) || (gchars<3)){
                    readok=1;
                    sprintf(string,"Error encountered during execution\n\tInvalid number of characters defined in input file. Must be between 3 and %d\n", Maxchar);
                    DoError(string, 'n');
                    return readok;
                }
                found[1] = 'y';
                if (gfrom>gchars){
                    DoError("From parameter must be less than number of characters in input file", 'n');
                }
                else if (gfrom<0){
                    gfrom=0;
                }
                fprintf(flog,"\tNumber of characters = %d\n", gchars);
                if (interactive=='n')printf("\tNumber of characters = %d\n", gchars);
                if (gto>gchars){
                    printf("\nError: To parameter greater than number of chars. Resetting to %d\n\n", gchars);
                    gto=gchars;
                }
                else if (gto==0){
                    gto=gchars;
                }
            }

            else if (strncasecmp(temp, "symbols", 7)==0){
                gnostates = 0;
                fseek (ftemp, spot, SEEK_SET);
                fscanf (ftemp, " %[^=]=", &temp[0]);
                fscanf (ftemp, " \"");
                fscanf(ftemp, " %c", &temp[0]);
                while (temp[0] !='"'){
                    if (gnostates>(Maxstates-1)){
                        readok=1;
                        sprintf(string,"Error encountered during execution\n\tToo many character states in state definition. Maximum = %d\n", Maxstates);
                        DoError(string, 'n');
                        return readok;
                    }
                    gstate[gnostates] = temp[0];
                    for (i=0; i<gnostates; i++){
                        if (gstate[i] == gstate[gnostates]){
                            readok=1;
                            DoError("Repeated state in state definition", 'n');
                            return readok;
                        }
                    }
                    gnostates++;
                    fscanf(ftemp, " %c", &temp[0]);
                }

                fprintf(flog,"\t%d states \(%s)\n", gnostates, gstate);
                if (interactive=='n')printf("\t%d states \(%s)\n", gnostates, gstate);
                found[2] = 'y';
            }

            else if (strncasecmp(temp, "datatype", 7)==0){
                gnostates = 0;
                fseek (ftemp, spot, SEEK_SET);
                fscanf (ftemp, " %[^=]=", &temp[0]);
                fscanf(ftemp, " %s", temp);
                if (strncasecmp(temp, "nucleotide", 10)==0)
                {
                    gnostates=5;
                    strcpy(gstate,"acgtu");
                    datatype='n';
                }
                else if (strncasecmp(temp, "dna", 3)==0)
                {
                    gnostates=5;
                    strcpy(gstate,"acgtu");
                    datatype='d';
                }
                else if (strncasecmp(temp, "rna", 3)==0)
                {
                    gnostates=5;
                    strcpy(gstate,"acgu");
                    datatype='r';
                }
                else if (strncasecmp(temp, "protein", 7)==0)
                {
                    gnostates=20;
                    strcpy(gstate,"arndceqghilkmfpstwyv");
                    datatype='p';
                }
                else if (strncasecmp(temp, "standard", 8)==0)
                {
                    datatype='s';
                }
                else {
                    readok=1;
                    DoError("Illegal datatype defined", 'n');
                    return readok;
                }

                fprintf(flog,"\tDatatype = %s\n", temp);
                if (interactive=='n')printf("\tDatatype = %s\n", temp);
                found[2] = 'y';
            }


            else if (strncasecmp(temp, "missing", 7)==0){
                fseek (ftemp, spot, SEEK_SET);
                fscanf (ftemp, " %[^=]=", &temp[0]);

                fscanf (ftemp, " %c", &missing);

                fprintf(flog,"\tMissing = %c\n", missing);
                if (interactive=='n')printf("\tMissing = %c\n", missing);
            }

            else if (strncasecmp(temp, "gap", 3)==0){
                fseek (ftemp, spot, SEEK_SET);
                found[5] = 'y';
                fscanf (ftemp, " %[^=]=", &temp[0]);

                fscanf (ftemp, " %c", &gap);

                fprintf(flog,"\tGap = %c\n", gap);
                if (interactive=='n')printf("\tGap = %c\n", gap);
            }



            else if (strncasecmp(temp, "matrix", 6)==0){

                if (found [0]!='y'){
                    readok=1;
                    DoError("Could not find number of taxa in input file", 'n');
                    return readok;
                }
                else if (found [1]!='y'){
                    readok=1;
                    DoError("Could not find number of characters in input file", 'n');
                    return readok;
                }
                else if (found [2]!='y'){
                    gstate[0]='0';
                    gstate[1]='1';
                    gnostates=2;
                }
                found[3] = 'y';
            }
        
        else if (strncasecmp(temp, "begin", 5)==0){
            fseek (ftemp, spot, SEEK_SET);
            //printf ("HERE1");
            fscanf (ftemp, " %s", temp);
            if (strncasecmp(temp, "compass", 7)==0){
                Doreadcompassblock(ftell(ftemp));
            }
        }
        }
        else {
            readok=1;
            DoError("Could not find matrix in input file", 'n');
            return readok;
        }



}

/*Allocate Matrix Memory*/
tempMatrix = malloc((gchars+1) * sizeof(char *));
if (tempMatrix == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory for matrix", 'n');
} 
for(i = 0; i < (gchars+1); i++){
    tempMatrix[i] = malloc((gtaxa+1) * sizeof(char));
    if (tempMatrix[i] == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
    DoError("Couldn't allocate memory for matrix", 'n');
} 
}


/*read in matrix*/
for (i=0; i<gtaxa; i++){
    fscanf (ftemp, " %99s", temp); /*input taxon name to 'taxonnames'*/
            for (a=0; a<30; a++){
                taxonnames[i][a]=temp[a];
            }
            for (j=0; j<gchars; j++){
                fscanf(ftemp, " %c", &tempMatrix[j][i]); /*input character state to 'matemp'*/
                       tempMatrix[j][i]=tolower(tempMatrix[j][i]);
                       if ((datatype=='d')||(datatype=='n')||(datatype=='r')){
                           if((tempMatrix[j][i]=='r')||(tempMatrix[j][i]=='m')||(tempMatrix[j][i]=='w')||(tempMatrix[j][i]=='y')||(tempMatrix[j][i]=='s')||(tempMatrix[j][i]=='k')||(tempMatrix[j][i]=='v')||(tempMatrix[j][i]=='h')||(tempMatrix[j][i]=='d')||(tempMatrix[j][i]=='b')||(tempMatrix[j][i]=='n')){
                               tempMatrix[j][i]=missing;
                           }
                       }
                if((tempMatrix[j][i] != missing)&&(tempMatrix[j][i]!= '?')&&(tempMatrix[j][i]!= gap)){
                    yn = 'n';
                    for(a=0; a<gnostates; a++){
                        if(tempMatrix[j][i] == gstate[a]){
                            yn = 'y';
                            break;
                        }
                    }
                    if (yn != 'y'){
                        readok=1;
                        sprintf(string,"Unidentified character state in matrix: %c taxon %d, character %d\n", tempMatrix[j][i], i+1, j+1);
                        DoError(string, 'n');
                        return readok;
                    }
                }
            }
}


fscanf (ftemp, " %c", &temp[0]);
if (temp[0] != ';'){
    readok=1;
    DoError("Expecting \';\' after matrix", 'n');
    return readok;
}


orignum = malloc((gchars+1) * sizeof(int));
if (orignum == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory for charnum", 'n');
}


gcharsnew=removeuninf();

alignlength=gchars;


/*Allocate Matrix Memory*/
gMatrix = malloc((gcharsnew+1) * sizeof(char *));
if (gMatrix == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory for matrix", 'n');
} 
for(i = 0; i < (gcharsnew+1); i++){
    gMatrix[i] = malloc((gtaxa+1) * sizeof(char));
    if (gMatrix[i] == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
    DoError("Couldn't allocate memory for matrix", 'n');
} 
}


//produce reduced matrix
for (i=0; i<gtaxa; i++){
	a=0;
	for (j=0; j<gchars; j++){
		if (keep[j]=='y'){
			   orignum[a]=j+1;
	           gMatrix[a][i]=tempMatrix[j][i];
	           a++;
		}
		}
}


printf ("\nRemoving %d sites that are uninformative for memory reasons\n", gchars-gcharsnew);



gchars=gcharsnew;
gto=gchars;

free(keep);

//free allocated memory
for (i=0; i<gchars; i++){
    free(tempMatrix[i]);
}
free(tempMatrix);
/*Allocate compatibility matrix Memory*/
compat = malloc((gchars+1) * sizeof(int **));
if (compat == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
for(i = 0; i < (gchars+1); i++){
    compat[i] = malloc((gchars+1) * sizeof(int *));
    if (compat[i] == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
    DoError("Couldn't allocate memory\n", 'n');
    }
    for(j = 0; j < (gchars+1); j++){
        compat[i][j] = malloc(2 * sizeof(int));
        if (compat[i][j] == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
    DoError("Couldn't allocate memory\n", 'n');
        }
    }
}

/*Allocate iterationincompats Memory*/
iterationincompats = malloc((gchars+1) * sizeof(int *));
if (iterationincompats == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
for(i = 0; i < (gchars+1); i++){
    iterationincompats[i] = malloc((gPermutations+1) * sizeof(int));
    if (iterationincompats[i] == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
    DoError("Couldn't allocate memory\n", 'n');
} 
}

/*Allocate Other Memory*/
charnum = malloc((gchars+1) * sizeof(int));
if (charnum == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
included = malloc((gchars+1) * sizeof(int));
if (included == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
same = malloc((gchars+1) * sizeof(int));
if (same == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
actualincompat = malloc((gchars+1) * sizeof(int));
if (actualincompat == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
Teststatistic = malloc((gchars+1) * sizeof(float));
if (Teststatistic == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
tot = malloc((gchars+1) * sizeof(int));
if (tot == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
tot2 = malloc((gchars+1) * sizeof(float));
if (tot2 == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
greater = malloc((gchars+1) * sizeof(int));
if (greater == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
type = malloc((gchars+1) * sizeof(int));
if (type == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 
removedlastround = malloc((gchars+1) * sizeof(int));
if (removedlastround == NULL) {
   /* Memory could not be allocated, so print an error and exit. */
   DoError("Couldn't allocate memory\n", 'n');
} 



//fclose(ftemp);

return readok;
}


void Doreadcompassblock(long spot){

printf ("HERE!!!");

fseek (ftemp, spot, SEEK_SET);

char temp[100];
int exclude;


fscanf(ftemp, " %100s", temp);

if (strncasecmp(temp, "end", 3)==0){
        return;
    }
else if (strncasecmp(temp, "exclude", 7)==0){
        printf ("\nExcluding");
        fscanf (ftemp, " %[^=]=", &temp[0]);
        
        fscanf (ftemp, " %d", &exclude);
        
        if ((exclude > 0)&&(exclude<=gchars)){
            printf (" %d", exclude);
            delcharfast(exclude-1);
        }
        else {
            return;
        }
        
    }

}



/*void bootstrap(void){

    int i;
    int j;
    int randchar;



    for (i=0; i<gchars; i++){
        randchar=ran(gchars);
        strcpy(gMatrix[i], tempMatrix[randchar]);
    }

    fpaup = fopen(paupfiles, "w");


    fprintf(fpaup, "#NEXUS\n\n");
    fprintf(fpaup, "Begin data;\n");
    fprintf(fpaup, "Dimensions Ntax=%d Nchar=%d;\n", gtaxa, gchars);
    fprintf(fpaup, "Format symbols=\"%s\" Missing=%c Gap=%c;\n", gstate, missing, gap);
    fprintf(fpaup, "Matrix");

    for (i=0; i<gtaxa; i++){
        fprintf(fpaup, "\n%s\t", taxonnames[i]);
        for (j=0; j<gchars; j++){
            if (gMatrix[j][i]=='?'){fprintf(fpaup, "%c", missing);}
            else{fprintf(fpaup, "%c", gMatrix[j][i]);}
        }
    }

    fprintf(fpaup, ";\nEnd;");

    fclose(fpaup);


}*/



void LeQuesne(void){


    int i;
    int k;
    int l;
//    int timeb;
    float count=0.0;
    int tot=0;

//Clear matrices

    for (l=gfrom; l<gto; l++){
        for (k=gfrom; k<gto; k++){
            compat[l][k][0] = 0;
            compat[l][k][1] = 0;
        }
    }   /*clear compat matrix*/
for(l=0; l<gtaxa+1; l++){
    comparison[l] = '\0';
}


//timeb=clock();

//If CCSR selected do it quickly

if (gCompattype=='A'){
    for (i=gfrom; i<gto; i++){ /*for all characters compare with all other characters*/
        
        if ((included[i]==1)&&(type[i]==0)){

            for(l=0; l<gtaxa; l++){
                comparison[l] = gMatrix[i][l];}
            actualincompat[i]=0;
            comparetype(i, 0, -1);   /*call comparison with other characters*/
            count++;
            printf("%.2f%% complete\r", (count/(gchars-(nouninf+noequiv)))*100);
            //printf(".");
            //if (count==50){
             //   count=0;
              //  tot=tot+50;
               // printf(" \[%.2f characters complete]\r", float(tot)/gchars-(nouninf+noequiv));
                
           // }
            fflush(stdout);
        }
    }
}

//If there are permutations to run, call iterations function

if (gPermutations>0){
    iterations();
}

}



void iterations(void){

    int i;
    int k;
    int l;
//    int timeb;
    float count=0.0;
    int tot=0;


    currentnumber=0;
    currentchar=0;
//    timeb=clock();

//clear matrices

    for (i=gfrom; i<gto; i++)
    {
        if (included[i]==1){
            greater[i] = 0;
            actualincompat[i]=0; //clears actualincompat
        }
    }

    for (i=gfrom; i<gto; i++)
    {

//        currentchar=i;
//        startseed=time(NULL);
//        printf("\n%d) %ld ", charnum[i], startseed+charnum[i]);
        srand( charseed[i] ); //sets random seed so it can be recalled
//        fprintf(fout, "\ni=%d seed=%ld\n", i, charseed[i]);
//        printf("\n%ld ", startseed);
        

        for(l=0; l<gtaxa; l++){
            comparison[l] = gMatrix[i][l];}

        comparetype(i, 0, -1);
        
        if ((included[i]==1)&&(type[i]==0)){


            k=0;
            while(k<gPermutations)
            {


                randomise (i);  /*call function to randomise character*/
                permutationincompat=0;
                comparetype(i, 1, k);   /*call comparison with other characters*/

                if((actualincompat[i]>=permutationincompat) && (actualincompat[i]>0))
                {
                    greater[i]++;
                }
                k++;
            }
            count++;
            
            printf("%.2f%% complete\r", (count/(gchars-(nouninf+noequiv)))*100);
            //printf(".");
            //if (count==50){
            //    count=0;
            //    tot=tot+50;
            //    printf(" \[%d characters complete]\n", tot);
            //}
            fflush(stdout);
            }//if same[i]==-1
    }
}

void LeQuesnefast(void){


    int i;
    int l;


/*    for (l=0; l<gchars; l++){
        for (k=0; k<gchars; k++){
            compat[l][k][0] = 0;
            compat[l][k][1] = 0;
        }
    }*/   /*clear compat matrix*/
for(l=0; l<gtaxa+1; l++){
    comparison[l] = '\0';
}


if (gCompattype=='A'){
    for (i=gfrom; i<gto; i++){ /*for all characters compare with all other character*/
        
        if ((included[i]==1)&&(type[i]==0)){

            for(l=0; l<gtaxa; l++){
                comparison[l] = gMatrix[i][l];}
            actualincompat[i]=0;
            comparetype(i, 0, -1);   /*call comparison with other characters*/
//            printf("%d\n", i+1);


        }

        //		addconflicts(0);  /*call function to add conflicts together*/

    }
}

if (gPermutations>0){
    iterationsfast();
}
}



void iterationsfast(void){

    int i;
    int iteration;
    int l;
    /*    int numberupdate;
    int currentnumber=0;
    int currentchar=0;*/


    currentnumber=0;
    currentchar=0;
    for (i=gfrom; i<gto; i++)
    {
        if (included[i]==1){
            greater[i] = 0;
//            actualincompat[i]=0; //clears actualincompat
        }
    }

    for (i=gfrom; i<gto; i++)
    if (actualincompat[i]>0)//if character has no incompatibilities don't bother recalculating in boildown
    {
    {
        currentchar=i;
//        startseed=time(NULL);
//        printf("\n%d) %ld ", charnum[i], startseed+charnum[i]);
//        fprintf(fout, "\ni=%d seed=%ld\n", i, charseed[i]);
        srand( charseed[i]);
//        printf("\n%ld ", startseed);
        

        for(l=0; l<gtaxa; l++){
            comparison[l] = gMatrix[i][l];}

        comparetypefast(i, 0, -1);
        if ((included[i]==1)&&(type[i]==0)){
            

            iteration=0;
            while(iteration<gPermutations)
            {
                randomise (i);  /*call function to randomise character*/
                permutationincompat=0;
                comparetypefast(i, 1, iteration);   /*call comparison with other characters*/

                if((actualincompat[i]>=iterationincompats[i][iteration]) && (actualincompat[i]>0))
                {
                    greater[i]++;
                }
                iteration++;
            }
            }//if same[i]==-1
            }
    }
}




void randomise(int i){

    int l;
    int n;
    char tempchar[Maxtax];
    char newchar[Maxtax];
    int tempcount = 0;


    if (gMissingconstant == 1)
    {
        for(l=0; l<gtaxa; l++)
        {
            if((gMatrix[i][l] != missing)&&(gMatrix[i][l]!= '?')&&(gMatrix[i][l]!= gap))
            {
                comparison[l] = '\0';
                newchar[tempcount] = '\0';
                tempchar[tempcount]=gMatrix[i][l];
                tempcount++;
            }
            else
            {
                comparison[l] = missing;
            }
        }
    }
    else
    {
        for(l=0; l<gtaxa; l++)
        {
            comparison[l] = '\0';
            tempchar[tempcount]=gMatrix[i][l];
            tempcount++;
        }
    }
    tempchar[tempcount]='\0';
    newchar[tempcount]='\0';

    ressler(i, tempchar, newchar, tempcount);


    n=0;
    for(l=0; l<gtaxa; l++)
    {

        if (comparison[l]=='\0')
        {

//            m=ran(tempcount);
            comparison[l]=newchar[n];
            n++;

//            for (n=m; n<tempcount; n++)
//            {
//                newchar[n]=newchar[n+1];
//            }
//            tempcount--;
        }
    }
}

void ressler(int i, char tempchar[Maxtax], char newchar[], int tempcount){

int n;
int m;
//printf (" %d", tempcount);

        m=ran(tempcount);
//        fprintf(fout, " %d", m);

        newchar[tempcount-1]=tempchar[m];
//        printf(" %c", tempchar[m]);

        for (n=m; n<tempcount; n++)
        {
            tempchar[n]=tempchar[n+1];
        }
        tempcount--;
        if (tempcount>0)
            {
                ressler(i, tempchar, newchar, tempcount);
            }
        else{
            return;
            newchar[tempcount]=tempchar[tempcount];
//            printf (" %s %c", newchar, newchar[0]);
        }


}


void comparetype (int i, int h, int iteration){

    int j;
//    int numberupdate;

//    numberupdate=(((gchars*gchars)+gchars)/2)-gchars;
//    numberupdate=(numberupdate*(gPermutations+1));
    //    time=clock();

    if ((gCompattype=='C')||(gCompattype=='N')||(gCompattype=='A'))
    {
        for (j=i+1; j<gchars; j++)
        {        //for each character
            compare (i, h, j, iteration, 0);
        }
    }
else if (gCompattype=='L')
{
    for (j=0; j<gchars; j++)
    {        //for each character
        if (included[j]==1){
            compare (i, h, j, iteration, 0);
            }
//        if ((actualincompat[i]<permutationincompat)&&(h==1)){
//            return;
//        }
    }
    }
else {DoError("Invalid compatibility type", 'y');}

}


void comparetypefast (int i, int h, int iteration){//fast version for testing

    int j;
    int numberupdate;

    numberupdate=(((gchars*gchars)+gchars)/2)-gchars;
    numberupdate=(numberupdate*(gPermutations+1));
    //    time=clock();

    if ((gCompattype=='C')||(gCompattype=='N')||(gCompattype=='A'))
    {
    for (j=0; j<noremovedlast; j++)
        {        /*for each character removed last round*/
            compare (i, h, removedlastround[j], iteration, 1); //is the 0 correct???
        }
    }
else if (gCompattype=='L')
{
    for (j=0; j<noremovedlast; j++)
    {        /*for each character removed last time around*/
        compare (i, h, removedlastround[j], iteration, 1);
//        if ((actualincompat[i]<permutationincompat)&&(h==1)){
//            return;
//        }
    }
}
else {DoError("Invalid compatibility type", 'y');}


}


void compare (int i, int h, int j, int iteration, int fast){

    int k;
    int l;
    int m;

    if ((gCompattype=='C')||(gCompattype=='N')||(gCompattype=='A')){

        if (currentnumber<currentchar){
            currentnumber=currentchar;
            //currentchar++;

        }
    }


    tns=0;
    if (j!=i)
    {
        /*except the character in question*/
        for (k=0; k<gnostates; k++)
        {        /*clear state comparison matrix and other variables*/
            for (l=0; l<gnostates; l++)
            {
                d[k][l] = 0;  /*0 in matrix = combination of states not present*/
            }
            xaxis[k]=0;
            yaxis[k]=0;
        }
        xmax=0;
        ymax=0;

        for (m=0; m<gtaxa; m++)
        {            /*fill in state comparison matrix*/

            for (k=0; k<gnostates; k++)
            {
                for (l=0; l<gnostates; l++)
                {
                    /*for each state combination checks to see if it is present in a character pair. If so, it adds 1 to that cell in combination matrix. If the cell was originally 0, adds 1 to the count of incompatibilities in each axis and the total number of incompatible cells(tns). Also increases the xmax and ymax to increase the size of the combination matrix examined*/
                    if ((comparison[m] == gstate[k]) && (gMatrix[j][m] == gstate[l]))
                    {

                        if (d[k][l] == 0)
                        {
                            if (k>(xmax-1))
                            {
                                xmax=(k+1);
                            }
                            if (l>(ymax-1))
                            {
                                ymax=(l+1);
                            }
                            xaxis[k]++;
                            yaxis[l]++;
                            tns++;
                        }
                        d[k][l]++;    /*>0 in matrix = combination of states present*/

                    }
                }
            }
        }
    }

if (tns>3)
{  	         /*If no of incompatibilities in matrix < 4, chars cannot be incompatible*/

resetmatrix( i, j, h, iteration, fast);  /*Else call resetmatrix function to check for incompatibility*/

}



//    }
// exit:;
}




int countstates(char character[Maxtax]){

    int i;
    int j;
    int k = 0;

    for (i=0; i<gnostates; i++){
        for (j=0; j<gtaxa; j++){

            if (character[j] == gstate[i]){
                k++;
                break;
            }

        }
    }

    return k;
}




void resetmatrix( int chr1, int chr2, int h, int iteration, int fast){

    /*int i;
    int j;*/
    int minincomp = 0;

    

    delnoncorners(); /*call delnoncorners to delete any cells in comparison matrix that are not corners of a loop. If corners remain after this, then characters are incompatible*/
    
//    fprintf (fout, "\nchr1=%d, chr2=%d, h=%d, iteration=%d, fast=%d, tns=%d\n", chr1, chr2, h, iteration, fast, tns);

    if (tns>3) /*if tns is<4 characters cannot be incompatible*/
    {
        if (gFuzzy==1) /*if fuzzy compatibility is selected*/
        {
            minincomp = delloop(chr1, chr2, h); /*call delloop function to break loops. Returns the mininum number of taxa that had to be removed to break loops*/
            if (h==0) /*if h=0, i.e. if this is the original, unpermuted data*/
            {
                
                if (gCompattype!='L'){
                    compat[chr2][chr1][h] = minincomp; /*put minincomp value into actual compat matrix*/	actualincompat[chr2] = actualincompat[chr1]+minincomp;
                }
                
                if (fast==0){
                    actualincompat[chr1] = actualincompat[chr1]+minincomp;
                    compat[chr1][chr2][h] = minincomp; /*put minincomp value into actual compat matrix*/
                }
                else {
                    actualincompat[chr1] = actualincompat[chr1]-minincomp;
                    compat[chr1][chr2][h] = 0;
                }
                
            }
            else if (h==1) /*if this is a permutation of the data*/
            {
                
                if (gCompattype!='L'){
                    compat[chr2][chr1][h] = (compat[chr2][chr1][h]+minincomp); /*add minincomp value to permutation compat matrix*/
                }
                if (fast==0){
                    compat[chr1][chr2][h] = (compat[chr1][chr2][h]+minincomp); /*add minincomp value to permutation compat matrix*/
                }
                else {
                    compat[chr1][chr2][h] = (compat[chr1][chr2][h]-minincomp);
                }
                permutationincompat=permutationincompat+minincomp;
                if (iteration>-1){
                    if (fast==0){
                        iterationincompats[chr1][iteration]=iterationincompats[chr1][iteration]+minincomp;
                    }
                    else{
                        iterationincompats[chr1][iteration]=iterationincompats[chr1][iteration]-minincomp;
                    }
                }
            }
        }
        else if (gFuzzy==0) /*if fuzzy is not selected*/
        {
            if (h==0)
            {
                 /*set actual compat matrix cell to 1*/
                if (gCompattype!='L'){
                    compat[chr2][chr1][h] = 1; /*set actual compat matrix cell to 1*/
                    actualincompat[chr2]++;
                }
                
                if (fast==0){
                    compat[chr1][chr2][h] = 1;
                    actualincompat[chr1]++;
                }
                else {
                    compat[chr1][chr2][h] = 0;
                    actualincompat[chr1]--;
                }
                
            }
            else if (h==1)
            {
                
                if (fast==0){
                    compat[chr1][chr2][h]++; /*add 1 to iterations compat matrix cell*/
                }
                else{
                    compat[chr1][chr2][h]--;
                }
                if (gCompattype!='L'){
                    compat[chr2][chr1][h]++; /*add 1 to iterations compat matrix cell Is this correct???*/
                }
                permutationincompat++;
                if (iteration>-1){
                    if (fast==0){
                        iterationincompats[chr1][iteration]++;
                    }
                    else{
                        iterationincompats[chr1][iteration]--;
                    }
                }
            }
        }	
    }
// fprintf (fout, "actual incompat=%d iterationincompat=%d\n", actualincompat[chr1], iterationincompats[chr1][iteration]);
}


void delnoncorners(void){

    int i;
    int j;
    char noones;

    do{
        noones='y';
        for (i=0; i<gnostates; i++){
            if (xaxis[i]==1){
                xaxis[i]--;
                for (j=0; j<gnostates; j++){		
                    if (d[i][j]>0){
                        d[i][j]=0;
                        tns--;
                        yaxis[j]--;
                        break;
                    }
                }
            }
        }

        for (j=0; j<gnostates; j++){
            if (yaxis[j]==1){
                noones='n';
                yaxis[j]--;
                for (i=0; i<gnostates; i++){		
                    if (d[i][j]>0){
                        d[i][j]=0;
                        tns--;
                        xaxis[i]--;
                        break;
                    }
                }
            }
        }

    }
    while ((noones!='y')&&(tns>3));

}



int delloop(int chr1, int chr2, int h){

    int i;
    int j;
    int k=0;
    int minvalue;
    //int maxvalue;
    int minincomp = 0;
    int stai;
    int staj;
    /*char noones;*/
    int statechar1=0;
    int statechar2=0;
    int e[Maxstates][Maxstates];

    do
    {	/*finds point with smallest number of taxa*/

        stai=0;
        staj=0;

        for (i=0; i<gnostates; i++)
        {
            if (xaxis[i]>0)
            {
                stai++;
            }
            if (yaxis[i])
            {
                staj++;
            }
        }

        minvalue=Maxtax;
        for (i=0; i<gnostates; i++)
        {
            for (j=0; j<gnostates; j++)
            {
                if ((d[i][j]<minvalue)&&(d[i][j]>0))
                {
                    minvalue=d[i][j];
                }
            }
        }

        for (i=0; i<gnostates; i++)
        {		/*Deletes first point with smallest number of taxa*/
            for (j=0; j<gnostates; j++)
            {

                if (d[i][j]==minvalue)
                {
                    d[i][j]=0;
                    statechar1=i;
                    statechar2=j;
                    tns--;
                    xaxis[i]--;
                    yaxis[j]--;
                    goto next;
                }
            }
        }

next:
loop='y';
k=0;
if ((stai>3)&&(staj>3)&&(xaxis[statechar1]>1)&&(yaxis[statechar2]>1))
{ /*To sort out problem of points that are corners, but not in loops*/
loop='n';
for (i=0; i<gnostates; i++)
{
    for (j=0; j<gnostates; j++)
    {
        e[i][j]=d[i][j];		/*e matrix saves original matrix*/
    }
}


for (j=0; j<gnostates; j++)
{
    d[statechar1][j]=0;
}

i=0;

do 
{
    j=0;
    d[statechar1][statechar2]=1;
    do 
    {		
        if(e[statechar1][k]>0)
        {
            d[statechar1][k]=1;
            j=1;		
        }	
        k++;
    } while((k<gnostates)&&(j==0));

    ecs(0, statechar1, k, statechar1, statechar2);
    i++;
} while((i<xaxis[statechar1])&&(loop!='y'));


for (i=0; i<gnostates; i++)
{
    for (j=0; j<gnostates; j++)
    {
        d[i][j]=e[i][j];			/*reset d matrix from e*/
    }
}
/*printf("\n%d %d %c", charnum[chr1], charnum[chr2], loop);*/
}

/*for (i=0; i<xmax; i++){
printf("%d ", xaxis[i]);}
printf("\ny = ");
for (i=0; i<ymax; i++){
    printf("%d ", yaxis[i]);}
for (i=0; i<xmax; i++){
    printf("\n");
    for (j=0; j<ymax; j++){
        printf("%d ", d[i][j]);}}*/


delnoncorners();

if (loop=='y')
{
    minincomp=minincomp+minvalue;
    for (i=0; i<gtaxa; i++)
    {
        if ((comparison[i]==gstate[statechar1]) && (gMatrix[chr2][i]==gstate[statechar2]))
        {
            taxaremoved[i][h]++;
        }
    }
}
    } while (tns>3);

return minincomp;
}



void ecs(int start, int x, int y, int statechar1, int statechar2){ /*with why checks if a character pair forms a loop*/

int i;

d[x][y]=0;

for (i=start; i<gnostates; i++){

    if ((d[i][y]>0)){
        if ((i==statechar1)&&(y==statechar2)){
            loop='y';
        }
        if (loop!='y'){
            why(0, i, y, statechar1, statechar2);
        }
    }
}

}


void why(int start, int x,int y, int statechar1, int statechar2){

    int i;

    for (i=start; i<gnostates; i++){

        if (d[x][i]>0){
            if ((x==statechar1)&&(i==statechar2)){
                loop='y';
            }
            if (loop!='y'){
                ecs(0, x, i, statechar1, statechar2);
            }
        }
    }

}







void boildown(void){

    int i;
    int noexcluded=0;
    int prevexcluded;
    int j=0;
    int z=0;
    int nummin=0;
    float maxincompat = 0;
    float minincompat = 0;
    float matrixincompat;


    if (gCompattype=='L')
    {
        if (gstopvalue==0){
            minincompat=(gPermutations+1);
            minincompat=1/minincompat;
            }
        else{
            minincompat=gstopvalue;
            }
    }
    
    fprintf(flog,"\n\nBoildown:");
    printf("\n\nBoiling Down:\n");

    //fprintf(fout, "\ninclude all; cleartrees; hsearch addseq=random nreps=10 rseed=1; savetrees file=%dexcluded.tre append;\n", noexcluded);

    do
    {
    
        if (z!=0)
        {
            //fprintf(fout, "\nexclude");
            fprintf(flog,"\n%d) Character(s) removed: ", z);
            printf("\n%d) Removing character(s): ", z);

 //           if (gBoot=='y'){
 //               printf("\nExclude");
 //           }

            prevexcluded=noexcluded;
            j=0;
            for (i=0; i<gchars; i++)
            {
                if (Teststatistic[i]==maxincompat && included[i]==1)
                {
                    fprintf(flog," %d", orignum[i]);
                    printf(" %d", orignum[i]);
//                    if (gBoot=='y'){
//                        printf(" %d", charnum[i-j]);
//                    }
                    delcharfast(i);
                    j++;
                    noexcluded++;
                }
/*                if ((gCompattype=='L')&&(gRemove==1)&&(Teststatistic[i-j]==minincompat))
                {
                    if (gBoot=='y'){
                        printf(" %d", charnum[i-j]);
                    }
                    delchar(i-j, charnum[i-j]-1, 0);
                    j++;
                    noexcluded++;
                    nummin++;
                }*/
            }

            if (gCompattype=='C'){fprintf(flog," at a CCSR");}
            else if (gCompattype=='N'){fprintf(flog," at an NDev");}
            else if (gCompattype=='L'){fprintf(flog," at an LQP");}
            else if (gCompattype=='A'){fprintf(flog," with an incompatibility");}

            if (gCompattype=='A') {fprintf(flog," value of %.0f.", maxincompat);}
            else {fprintf(flog," value of %.3f.", maxincompat);}
            
            if (gCompattype=='C'){printf(" at a CCSR");}
            else if (gCompattype=='N'){printf(" at an NDev");}
            else if (gCompattype=='L'){printf(" at an LQP");}
            else if (gCompattype=='A'){printf(" with an incompatibility");}

            if (gCompattype=='A') {printf(" value of %.0f.", maxincompat);}
            else {printf(" value of %.3f.", maxincompat);}

            if (gCompattype=='C'){fprintf(flog," Matrix CCSR = %.3f.\n", matrixincompat);}
            else if (gCompattype=='N'){fprintf(flog," Matrix NDev = %.3f.\n", matrixincompat);}
            else if (gCompattype=='L'){fprintf(flog," Matrix LQP = %.3f.\n", matrixincompat);}
            else if (gCompattype=='A'){fprintf(flog," Matrix average incompatibility = %.3f.\n", matrixincompat);}

/*            if (gBoot=='y')
            {
                printf("; Bootstrap NReps=%d GrpFreq=No ", gBootreps);
                if (gLE==1){
                    printf("Keepall=yes ");
                }
                

                switch (gSearchtype)
                {
                    case 'H':
                        printf("Search=Heuristic/ AddSeq=");
                        switch (gAdditiontype)
                        {
                            case 'S':
                                printf("Simple");
                                break;
                            case 'C':
                                printf("Closest");
                                break;
                            case 'A':
                                printf("AsIs");
                                break;
                            case 'R':
                                printf("Random NReps=%d", gHeuristicrands);
                                if (gSeed!=0){
                                    printf(" RSeed=%d", gSeed);
                                }
                                break;
                        }
                        break;
                    case 'B':
                        printf("Search=BandB");
                        break;
                }

                printf("; Savetrees From=1 To=1 File=%s.tre SaveBootP=Brlens Append=yes;", PAUPinfile);

            }*/

         
            
//            for (i=0; i<Maxchar; i++){
//            printf(" %d", included[i]);}
            
            //fprintf(fout, ";\ncleartrees; hsearch addseq=random nreps=10 rseed=1;");

            //for (i=prevexcluded; i<noexcluded; i++)
            //{
            //    fprintf(fout, " savetrees file=%dexcluded.tre append;", i+1);
            //}

//            gchars=gchars;

            if (gCompattype=='L')
            {
//                startseed=time(NULL);
//                srand( startseed );
//                printf("\n%ld ", startseed);
                LeQuesne();
            }

            CalculateTestValues();
//            DoTestOutput();//REMOVE!!!
        }

        maxincompat=0;
        matrixincompat=0;
        maxincompat=0;
        matrixincompat=0;
        for (i=gfrom; i<gto; i++)
        {
            if (Teststatistic[i]>maxincompat && included[i]==1)
            {
                maxincompat=Teststatistic[i];
            }
            if (included[i]==1){
                matrixincompat=matrixincompat+Teststatistic[i];
                //printf("\n stat%f %d %d", Teststatistic[i], charnum[i], type[i]);
            }
        }

        matrixincompat=(matrixincompat+(minincompat*nummin));
        matrixincompat=(matrixincompat/(gchars+nummin));

        z++;
    }while ((maxincompat>0)&&(gchars>0)&&(maxincompat>minincompat));

    
    if (gCompattype=='L')
    {
        fprintf(flog,"\nAll remaining characters are compatible at the %.3f level\n", minincompat);
    }
    else {fprintf(flog,"\nAll remaining characters are compatible\n");}
    
}


void boildownfast(void){

    int i;
    int noexcluded=0;
    int prevexcluded;
    int j=0;
    int z=0;
    float maxincompat = 0;
    float minincompat = 0;
    float matrixincompat;
    int cutvalue=cutoff;

//    for (i=0; i<gchars; i++)
//    {
//        removedlastround[i]=0;
//    }


    if (gCompattype=='L')
    {
        if (gstopvalue==0){
            minincompat=(gPermutations+1);
            minincompat=1/minincompat;
            }
        else{
            minincompat=gstopvalue;
            }
    }
    
    fprintf(flog,"\n\nBoildown:");
    printf("\n\nBoiling Down:\n");

    do
    {
    
        if (z!=0)
        {
     //       printf ("\ni=%d actual=%d, iterations %d %d %d %d %d %d %d %d %d %d\n", 50, actualincompat[50], iterationincompats[50][0], iterationincompats[50][1], iterationincompats[50][2], iterationincompats[50][3], iterationincompats[50][4], iterationincompats[50][5], iterationincompats[50][6], iterationincompats[50][7], iterationincompats[50][8], iterationincompats[50][9]);
            fprintf(flog,"\n%d) Character(s) removed: ", z);
            printf("\n%d) Removing character(s): ", z);

            prevexcluded=noexcluded;
            j=0;
            for (i=gfrom; i<gto; i++)
            {
                //printf (" %d %f %f\n", charnum[i], Teststatistic[i], maxincompat);
                if ((Teststatistic[i]==maxincompat) && (included[i]==1))
                {
     //               if (charnum[i]==51){
    //                    printf ("\ni=%d actual=%d, iterations %d %d %d %d %d %d %d %d %d %d\n", i, actualincompat[i], iterationincompats[i][0], iterationincompats[i][1], iterationincompats[i][2], iterationincompats[i][3], iterationincompats[i][4], iterationincompats[i][5], iterationincompats[i][6], iterationincompats[i][7], iterationincompats[i][8], iterationincompats[i][9]);
     //               }
                    fprintf(flog," %d", orignum[i]);
                    printf(" %d", orignum[i]);
                    removedlastround[j]=charnum[i]-1;
                    delcharfast(i);
                    j++;
                    noexcluded++;
                }

            }
            noremovedlast=j;
            if (gCompattype=='C'){fprintf(flog," at a CCSR");}
            else if (gCompattype=='N'){fprintf(flog," at an NDev");}
            else if (gCompattype=='L'){fprintf(flog," at an LQP");}
            else if (gCompattype=='A'){fprintf(flog," with an incompatibility");}

            if (gCompattype=='A') {fprintf(flog," value of %.0f.", maxincompat);}
            else {fprintf(flog," value of %.3f.", maxincompat);}
            
            if (gCompattype=='C'){printf(" at a CCSR");}
            else if (gCompattype=='N'){printf(" at an NDev");}
            else if (gCompattype=='L'){printf(" at an LQP");}
            else if (gCompattype=='A'){printf(" with an incompatibility");}

            if (gCompattype=='A') {printf(" value of %.0f.", maxincompat);}
            else {printf(" value of %.3f.", maxincompat);}

            if (gCompattype=='C'){fprintf(flog," Matrix CCSR = %.3f.", matrixincompat);}
            else if (gCompattype=='N'){fprintf(flog," Matrix NDev = %.3f.", matrixincompat);}
            else if (gCompattype=='L'){fprintf(flog," Matrix LQP = %.3f.", matrixincompat);}
            else if (gCompattype=='A'){fprintf(flog," Matrix average incompatibility = %.3f.", matrixincompat);}
            
            if (gCompattype=='L')
            {
                LeQuesnefast();
            }

            CalculateTestValues();
            
   //         printf ("\n%d %d %d %f\n", gchars, nouninf, noexcluded, ((float)(noexcluded)/((float)(gchars)-(nouninf)))*100);
//            DoTestOutput();//REMOVE!!
            writematrix("boil");//REMOVE!!

//        gchars=gchars;
        }
        
        if ((maxincompat>cutvalue) && (mask=='y' || charset=='y')){
            
            if (mask=='y'){
                if (maxincompat!=0){
                    sprintf(string, "Incompatibility_%.3f =", maxincompat);
                    Domask(string, 'r');
                }
            }
            if (charset=='y'){
                if (maxincompat!=0){
                    sprintf(string, "Incompatibility_%.3f =", maxincompat);
                    Docharset(string, 'r');
                }
            }
            cutvalue=cutvalue+cutoff;
        }


        maxincompat=0;
        matrixincompat=0;
        for (i=gfrom; i<gto; i++)
        {
            if (Teststatistic[i]>maxincompat && included[i]==1)
            {
                maxincompat=Teststatistic[i];
            }
            if (included[i]==1){
                matrixincompat=matrixincompat+Teststatistic[i];
                //printf("\n stat%f %d %d", Teststatistic[i], charnum[i], type[i]);
            }
        }
        
        
        matrixincompat=(matrixincompat/(gchars-noexcluded));

        z++;
    }while ((maxincompat>0)&&(gchars>0)&&(maxincompat>minincompat));

    
    if (gCompattype=='L')
    {
        fprintf(flog,"\nAll remaining characters are compatible at the %.3f level\n", minincompat);
        printf("\nAll remaining characters are compatible at the %.3f level\n", minincompat);
    }
    else {fprintf(flog,"\nAll remaining characters are compatible\n");
        printf("\nAll remaining characters are compatible\n");}

    if (charset=='y'){
        sprintf(string, "Incompatibility_%.3f =", maxincompat);
        Docharset(string, 'r');
    }
    if (mask=='y'){
        sprintf(string, "Incompatibility_%.3f =", maxincompat);
        Domask(string, 'r');
    }

}


void delchar (int worst, int p, int uninf){/*Check code for compat since extra matrix added!*/

int i;
int j;


//if ((type[p]==0)||(type[p]==3)){type[p] = 4;}

delcharfast(worst);
for (i=worst; charnum[i] != '\0'; i++){   /*Deletes worst character & number of that character*/
     charnum[i] = charnum[i+1];
     same[i] = same[i+1];
     type[i] = type[i+1];
     for (j=0; j<gtaxa; j++)
     {
         gMatrix[i][j] = gMatrix[i+1][j];
     }
    
     if (uninf==0)
     {
         Teststatistic[i] = Teststatistic[i+1];
         /*tot2[i] = tot2[i+1];
         tot[i] = tot[i+1];*/
         for (j=0; j<gchars; j++){
             compat[i][j][0] = compat[i+1][j][0];
         }
         for (j=0; j<gchars; j++){
             compat[j][i][0] = compat[j][i+1][0];
         }
         for (j=0; j<gchars; j++){
             compat[i][j][1] = compat[i+1][j][1];
         }
         for (j=0; j<gchars; j++){
             compat[j][i][1] = compat[j][i+1][1];
         }
     }		
}

}

void delcharfast (int worst){/*Check code for compat since extra matrix added!*/

included[worst] = 0;

}




int finduninf(void)
{
    int i;
    int j;
    int p=0;
    int k;
    int l;
    int m;
    int nostate[Maxstates];
    char yn = 'y';


        nouninf=0;
        noequiv=0;
    
//for (i=0; i<gchars; i++){type[i]=0;}

    for (i=0; i<gchars; i++){//for each character

        same[i]=-1;

        type[i]=1;//set all characters to type 1 (uninformative) as default
        for (l=0; l<gnostates; l++){
            nostate[l] = 0;
        }
        k = countstates(gMatrix[i]);//count the number of states in the character
        if (k<=1){
            type[i]=2;//if the character has one or less states, set it to type 2 (constant)

        }
        else{

            for (l=0; l<gnostates; l++){
                for (m=0; m<gtaxa; m++){
                    if (gMatrix[i][m] == gstate[l]){
                        nostate[l]++;//for each state count the number of taxa with that state
                    }
                }
            }

            for (l=0; l<gnostates; l++){
                if (nostate[l] > 1){//check to see if more than one state is present in at least 2 taxa
                    for (m=(l+1); m<gnostates; m++){
                        if (nostate[m] > 1){
                            type[i]=0;//if so set characters to state 0 (informative)
                            break;
                        }
                    }
                    break;
                }
            }
        }
        
        
        if (type[i]==0){//if the character is informative, check to see if it is the same as any of the preceeding characters
            for (j=0; j<i; j++){
                yn = 'y';
                for (k=0; k<gtaxa; k++){
                    if (gMatrix[i][k]!=gMatrix[j][k]){
                        yn = 'n';
                        break;
                    }
                }
                if (yn == 'y'){
                    
                    same[i]=j;//if it is, link the two together using the same array
                    type[i]=3;//and change the type of the current cahracter to 3 (same as a previous char)
                    noequiv++;
                    break;
                }
            }
        }
    }

    for (i=0; i<gchars; i++){
        if ((type[i]==1)||(type[i]==2)){//if characters are uninformative 

            if (type[i]==2){p++;}

            delcharfast(i);//delete the character
           // compatiblechars[nouninf]=i+1;
            nouninf++;
        }
        
    }

//    gchars = gchars-nouninf;

    //HideControl (UninfCharsMessage);
    //ShowControl (EquivCharsMessage);

    

   /* for (i=1; i<gchars; i++){
        if (exittest=='y'){
            return p;
        }
        if(i==0){time=clock();}
        progressbarvalue=((i*100)/gchars);
        SetControlValue (ProgressBar, progressbarvalue);
        DrawOneControl (ProgressBar);
        timeb=clock();
        if(timeb-time>20){
            gotEvent = WaitNextEvent(everyEvent, &eventStructure, 0, NULL);
            time=clock();
        }
        for (j=0; j<i; j++){
            yn = 'y';
            for (k=0; k<gtaxa; k++){
                if (gMatrix[i][k]!=gMatrix[j][k]){
                    yn = 'n';
                    break;
                }
            }
            if (yn == 'y'){
            
                tempint=i-nouninf;
                same[charnum[i]]=charnum[j];
                
                tempint=charnum[j];
                type[charnum[i-nouninf]-1]=3;
                break;
            }
        }
    }*/
    //HideWindow (gUninfProgressWindow);

return p;
}




int removeuninf(void)
{
    int i;
    int j;
    int p=0;
    int k;
    int l;
    int m;
    int nostate[Maxstates];
    char yn = 'y';
    int gcharsnew=0;

	keep = malloc((gchars+1) * sizeof(char));
	if (keep == NULL) {
	   /* Memory could not be allocated, so print an error and exit. */
	   DoError("Couldn't allocate memory for keep\n", 'n');
	}
	for (i=0; i<gchars; i++){
	keep[i]='n';
	}

        nouninf=0;
        noequiv=0;
    
//for (i=0; i<gchars; i++){type[i]=0;}

    for (i=0; i<gchars; i++){//for each character
        for (l=0; l<gnostates; l++){
            nostate[l] = 0;
        }
        k = countstates(tempMatrix[i]);//count the number of states in the character
        if (k>1){

            for (l=0; l<gnostates; l++){
                for (m=0; m<gtaxa; m++){
                    if (tempMatrix[i][m] == gstate[l]){
                        nostate[l]++;//for each state count the number of taxa with that state
                    }
                }
            }

            for (l=0; l<gnostates; l++){
                if (nostate[l] > 1){//check to see if more than one state is present in at least 2 taxa
                    for (m=(l+1); m<gnostates; m++){
                        if (nostate[m] > 1){
                            gcharsnew++;
                            keep[i]='y';
                            break;
                        }
                    }
                    break;
                }
            }
        }
        
        

    }

return gcharsnew;
}






void Removebrakets(void){
    char temp;

    while (fscanf(ftempb, "%c", &temp)!=EOF) 
    {
        if ((temp != '[')&&(temp != '(')&&(temp != '{')){
            fprintf(ftemp, "%c", temp);
        }
        else if (temp == '[') {
            do {
                fscanf(ftempb, "%c", &temp);
            } while (temp != ']');
        }
        else if (temp == '(') {
            printf("%d", temp);
            do {
                fscanf(ftempb, "%c", &temp);
            } while (temp != ')');
        }
        else if (temp == '{') {
            fprintf(ftemp, "%c", missing);
            do {
                fscanf(ftempb, "%c", &temp);
            } while (temp != '}');
        }
    }

    fseek (ftemp, 0, SEEK_SET);

    /*  for (i=0; i<10; i++){
        fscanf (ftemp, "%s", tempstring);
    printf("%s", tempstring);
    }

fseek (ftemp, 0, SEEK_SET);*/
}



void CalculateTestValues (void)
{

    int i=0;
    int j;
	//printf("\n\nhere\n\n");
	
	for (i=0; i<gchars; i++)
    {
    	if (same[i]>=0)
        {
    	
    		for (j=0; j<gchars; j++){
    			compat[j][i][0]=compat[j][same[i]][0];
    			compat[j][i][1]=compat[j][same[i]][1];
    			compat[i][j][0]=compat[same[i]][j][0];
    			compat[i][j][1]=compat[same[i]][j][1];
            }
    	
    	
        }
    }
	
	
    for (i=0; i<gchars; i++)
    {
    	//if (included[m]==1){
        //if(type[m]==0)
        //printf("\nhere!!!%d\n", m);
        if (same[i]==-1)
        {
            if (included[i]==1){
            tot[i]=0;
            for (j=0; j<gchars; j++)
            {
                tot[i]=(tot[i]+compat[i][j][0]);
            }

            tot2[i]=0;
            for (j=0; j<gchars; j++){
                tot2[i]=(tot2[i]+compat[j][i][1]);
            }
            tot2[i]=(tot2[i]/gPermutations);
			
	            if (gCompattype=='C'){
	                Teststatistic[i]=(tot[i]/tot2[i]);
	            }
	            else if (gCompattype=='N')
	            {
	                if (tot2[i]<tot[i]){
	                    Teststatistic[i]=(tot2[i]-tot[i]+0.5)/(sqrt(tot2[i]*((((gchars-nouninf)*(gchars-(nouninf+1))/2)-tot[i]))/((gchars-nouninf)*(gchars-(nouninf+1))/2)));
	                }
	                else{
	                    Teststatistic[i]=(tot2[i]-tot[i]-0.5)/(sqrt(tot2[i]*( ( ((gchars-nouninf)*(gchars-(nouninf+1))/2) -tot[i])) /((gchars-nouninf)*(gchars-(nouninf+1))/2)));
	                }
	            }
	            else if (gCompattype=='L')
	            {
	                Teststatistic[i]=(((float) (greater[i]+1)) / (float) (gPermutations+1));
	            }
	            else if (gCompattype=='A')
	            {
	                Teststatistic[i]=(tot[i]);
	            }
			}
        }
       // else if(type[m]==3)
       else
        {
            tot[i]=tot[same[i]];
            tot2[i]=tot2[same[i]];
            Teststatistic[i]=Teststatistic[same[i]];
            //break;

        }
    	//}
    	
    }
}



void DoTestOutput (void)
{
    int i;

    int m;


    fprintf(flog,"\n\nResults:");
    fprintf(flog,"\nChar\t");
    if (gCompattype!='A'){ fprintf(flog,"\t Obs\t\t ");
        if ((gCompattype!='L')/*||(gFast==0)*/){ fprintf(flog,"Exp\t\t");}
    }

    if (gCompattype=='C')
    {
        fprintf(flog,"CCSR");
    }
    else if (gCompattype=='N')
    {
        fprintf(flog,"NDev");
    }
    else if (gCompattype=='L')
    {
        fprintf(flog,"LQP");
    }
    else if (gCompattype=='A')
    {
        fprintf(flog,"No. Incompatibilities");
    }
    i=0;
    
    for (m=0; m<gchars; m++)
    {
        //if ((type[m]==0)||(type[m]==3)||(gPrintUninf==1))
        if ((type[m]==0)||(type[m]==3))
        {
            fprintf(flog,"\n%4d:", orignum[m]);
                
            if ((type[m]==0)&&(included[m]==1)){//If character is informative and unique and included

                fprintf(flog,"\t%4d", tot[m]);
                if (gCompattype!='A'){
                    if ((gCompattype!='L'/*)||(gFast==0*/)){
                        fprintf(flog,"\t%6.1f", tot2[m]);
                    }
                    fprintf(flog,"\t%5.3f", Teststatistic[m]);
                }
                i++;
            }
            else if(type[m]==1){//If character is uninformative
                fprintf(flog,"\tUNINFORMATIVE CHARACTER\t0");
            }
            else if(type[m]==2){//If character is constant
                fprintf(flog,"\tCONSTANT CHARACTER\t0");
            }
            else if((type[m]==3)&&(included[m]==1))//If character is equivalent
            {

                fprintf(flog,"\tEQUIVALENT TO CHARACTER %d", orignum[same[m]]);
                if (gCompattype!='A')
                {
                    fprintf(flog,"\t%.3f", Teststatistic[m]);
                }
                else
                {
                    fprintf(flog,"\t%d", tot[m]);
                }
                i++;
            }
            else {
                fprintf(flog,"\tRemoved by boildown\t%.3f", Teststatistic[m]);
            }
//            else {DoError("Inapplicable character type");}
//            printf("\t%d", included[m]);
            
        }
        //else {i++;}
    } 

 fprintf(flog,"\n");


//    if (gCharset>0){
//        writecharset();
//    }
    
}

void Domask(char name[50], char cur)
{
int i;

fprintf(maskfile, "%s ", name);

if (cur=='c'){
    for (i=0; i<gchars; i++){
        if (type[i]==2){
            fprintf (maskfile,"*");
        }
        else {
            fprintf (maskfile,"-");
        }
    }
}
else if (cur=='u'){
    for (i=0; i<gchars; i++){
        if ((type[i]==1)||(type[i]==2)){
            fprintf (maskfile,"*");
        }
        else {
            fprintf (maskfile,"-");
        }
    }
}
else{
    for (i=0; i<gchars; i++){
        if (included[i]==0){
            fprintf (maskfile,"-");
        }
        else {
            fprintf (maskfile,"*");
        }
    }

}

fprintf (maskfile, "\n");

}

void Docharset(char name[50], char cur)
{
int i;

fprintf(charsetfile, "%s", name);

if (cur=='c'){
    for (i=0; i<gchars; i++){
        if (type[i]==2){
            fprintf (charsetfile," %d", charnum[i]);
        }
    }
}
else if (cur=='u'){
    for (i=0; i<gchars; i++){
        if ((type[i]==1)||(type[i]==2)){
            fprintf (charsetfile," %d", charnum[i]);
        }
    }
}
else{
    for (i=0; i<gchars; i++){
        if (included[i]==0){
            fprintf (charsetfile," %d", charnum[i]);
        }
    }

}



fprintf (charsetfile,";\n");
    
}

/*void Docharset(void)
{

    float top=1;
    float bottom;
    int i;
    int count=0;
    int first=0;

        do{

            bottom=top-gCharset;
            
            if (bottom<0){bottom=0;}
            
            count=0;
            first=0;
            for (i=0; i<gchars; i++){
                if ((Teststatistic[i]<=top)&&(Teststatistic[i]>bottom)){
                    first++;
                    count=1;
                    if ((count==1)&&(first==1)){
                        printf("\ncharset %.3f_to_%.3f =", bottom, top);
                    }

                    printf(" %d", charnum[i]);
                }
            }
            if (count==1){
                printf(";");
            }
            top=bottom;

        }while(bottom!=0);  

}*/


void writeplot(char name[5]){
    int i;
    int j;
    
    
	sprintf(outputfile, "%s.%c%s.plot", inputfile, gCompattype, name);
	
	fout=fopen(outputfile, "w");
	j=0;
	for (i=0; i<gchars; i++){
		j++;
		//printf("\n%d %d", i, orignum[i]);
		while (j<orignum[i]){
			fprintf(fout,"0\n");
			j++;	
		}
		if (gCompattype!='A')
		{
			fprintf(fout,"%.3f\n", Teststatistic[i]);
		}
		else
		{
			fprintf(fout,"%d\n", tot[i]);
		}
	}
	while (j<alignlength){
			fprintf(fout,"0\n");
			j++;	
		}
	fclose(fout);
}



void writematrix(char name[5]){
    int i;
    int j;
    int k;
    int l;
    
    if (gto>=0){
        sprintf(outputfile, "%s.%d.%d.%s", inputfile, gfrom+1, gto, name);
    }
    else {
        sprintf(outputfile, "%s.%d.end.%s", inputfile, gfrom+1, name);
    }
    
    fout=fopen(outputfile, "w");
    
    //fprintf(fout,"%d %d %d", gfrom, gto, gPermutations);
    
    //fprintf(fout,"\n\n");
    
    for (i=gfrom; i<gto; i++){
    	
    	if ((type[i]==0)||(type[i]==3)){
    		if (type[i]==3){
    			k=same[i];
    		}
    		else{
    			k=i;
    		}
    		//printf("\n%d %d %d %d", i, k,compat[i][0][0],compat[k][0][0]);
	        for (j=0; j<gchars; j++){
	        	
	    		if ((type[j]==0)||(type[j]==3)){
	    			if (type[j]==3){
	    				l=same[j];
		    		}
		    		else{
		    			l=j;
		    		}
		    		
	            	fprintf(fout," %d",compat[k][l][0]);
	    		}
	        }
	        fprintf(fout,"\n");
    	}
    }
    

    fclose(fout);
    
    return;
    
    for (i=gfrom; i<gto; i++){
        for (j=0; j<gchars; j++){
            fprintf(fout," %d",compat[i][j][0]);
        }
        fprintf(fout,"\n");
    }
    
        fprintf(fout,"\n\n");
    
    for (i=gfrom; i<gto; i++){
    //	if (type[i]==0){
	        for (j=0; j<gchars; j++){
	        	
    //			if (type[j]==0){
	            	fprintf(fout," %d",compat[i][j][1]);
    			//}
	        }
    	    fprintf(fout,"\n");
    	//}
    }
    //return;
    fprintf(fout,"\n\n");
    
    for (i=gfrom; i<gto; i++){
    fprintf(fout, " %d", greater[i]);
    }
   
   
    fprintf(fout,"\n\n");
     
    for (i=gfrom; i<gto; i++){
        fprintf(fout," %d",actualincompat[i]);
        for (j=0; j<gPermutations; j++){
            fprintf(fout," %d",iterationincompats[i][j]);
        }
        fprintf(fout,"\n");
    }
    
    fprintf(fout,"\n\n");
    
    for (i=gfrom; i<gto; i++){
        fprintf(fout, " %d", included[i]);
    }
    
    fprintf(fout,"\n\n");
    
    for (i=gfrom; i<gto; i++){
        fprintf(fout, " %d", type[i]);
    }
    fclose(fout);
}


int ran(int k)
{

    double  x = RAND_MAX + 1.0;
    int y;
    y = rand()*(k/x);
//    printf("%d  ", y);

    return y;
}


void Flush( void )
{
    while ( getchar() != '\n' )
        ;
}


void DoError( char *message, char usage)
{
    printf( "\n!!!Error: %s!!!\n\n", message );
    if (usage=='y')
        {
            Usage();
        }
    exit( 0 );
}

void Usage( void )
{
    printf("Usage:\n\n\t-i\tinput file in nexus format\n\t-t\tcompatibility type\t\tA = Number of pairwise Incompatibilities\n\t\t\t\t\t\tC = CCSR (coefficient of character-state randomness)\n\t\t\t\t\t\tL = LQP (LeQuesne Probability)\n\t\t\t\t\t\tN = Normal Deviate");
    printf("\n\t-p\tnumber of permutations\t\t0 - 10,000\n\t-b\trun boildown\t\t\t0 = no\n\t\t(default = 0)\t\t\t1 = y\n\t\t\t\t\t\t2 = fast boildown\n\t-f\tfuzzy analysis\t\t\t0 = no\n\t\t(default = 0)\t\t\t1 = y");
    printf("\n\t-s\tstart character (optional)\n\t-e\tend character (optional)\n\t-r\tstarting seed (optional, default = random)");
    printf("\n\t-m\tcreate mask file (optional)\n\t-c\tcreate charset file (optional)\n\t-x\tLQP cutoff for charset and mask (optional, default = 1)\n\t-I\tLaunch interactive mode\n\n");
}


void readdata(int number){

int start;
int end;
int perms;
int i;
int j;

printf("\nReading Data From \"%s\"...", datafile[number]);

if ((fdata=fopen(datafile[number], "r"))==NULL){
    sprintf(string,"Cannot open file \"%s\"", datafile[number]);
//    printf("%s", string);
    DoError(string, 'n');
}

fscanf(fdata, " %d %d %d", &start, &end, &perms);

if (gPermutations==0){
    gPermutations=perms;
    }
else if (gPermutations!=perms){
    DoError("Datasets have differents numbers of permutations", 'n');
}

for (i=start; i<end; i++){
    for (j=0; j<gchars; j++){
        fscanf(fdata," %d",&compat[i][j][0]);
    }
}

for (i=start; i<end; i++){
    for (j=0; j<gchars; j++){
        fscanf(fdata," %d",&compat[i][j][1]);
    }
}

for (i=start; i<end; i++){
    fscanf(fdata, " %d", &greater[i]);
}

for (i=start; i<end; i++){
    fscanf(fdata," %d",&actualincompat[i]);
    for (j=0; j<perms; j++){
        fscanf(fdata," %d",&iterationincompats[i][j]);
    }
}


    for (i=start; i<end; i++){
        fscanf(fdata, " %d", &included[i]);
    }
    
    for (i=start; i<end; i++){
        fscanf(fdata, " %d", &type[i]);
    }

fclose(fdata);

}


int main(int argc, char *argv[])
{

    struct tm *systime;
    
    clock_t c;
    long days;
    long hours;
    long mins;
    long secs;
    int i;
    int j;
    int executeok=1;
    
    gto=0;
    gfrom=-1;
    inputfile[0]='\0';
    outputfile[0]='\0';
    datafile[0][0]='\0';
    gCompattype='L';
    gPermutations=100;
    gFuzzy=0;
    gBoildown=0;


    startseed = time(NULL);
    systime = localtime (&startseed);
    printheading();


//c=clock();

//startseed=1160728179;


//input = boildown inputfile compattype permutations fuzzy boildown > outputfile

//check that inputfile and compattype are defined

if (argc==1){
    interactive='y';
}

for (i = 1; i < argc; i++) {

	/* Check for a switch (leading "-"). */

	if (argv[i][0] == '-') {

	    /* Use the next character to decide what to do. */

	    switch (argv[i][1]) {
            
                case 'm':	mask='y'; //
				break;
                
                case 'c':	charset='y'; //
				break;
                                
                case 'x':	cutoff=atof(argv[i+1]); //
                                if ((cutoff>1)||(cutoff<=0)){
                                    DoError("Print frequency must be a percentage of characters removed", 'y');
                                }
				break;

		case 'i':	strcpy(inputfile,argv[i+1]); //read inputfile name and execute the file. Quit if file cannot be opened or executed
				break;
                
                case 'I':       interactive='y';
                                break;
                                
                case 'o':       strcpy(outputfile,argv[i+1]); //read inputfile name and execute the file. Quit if file cannot be opened or executed
				break;

		case 't':	gCompattype = toupper(argv[i+1][0]); //read compattype and check it is legal
                                if ((gCompattype!='C')&&(gCompattype!='A')&&(gCompattype!='L')&&(gCompattype!='N'))
                                {
                                    DoError("Illegal compatibility type selected", 'y');
                                }
				break;

		case 'p':	gPermutations = atoi(argv[i+1]); //extract number of permutations and check it adheres to program rules
                                if ((gPermutations<1)||(gPermutations>10000)){
                                    DoError("Number of permutations must be between 1 and 10,000", 'y');
                                }
				break;

		case 'b':	if ((argv[i+1][0]=='0')||(argv[i+1][0]=='1')||(argv[i+1][0]=='2')){ //extract boildown parameter if present and check it is '0' or '1'. If not present set to 0
                                    gBoildown = atoi(argv[i+1]);
                                }
                                else {
                                       DoError("Boildown parameter must be 0, 1 or 2", 'y');
                                }
				break;
                
                case 'f':	gFuzzy=1; //Set fuzzy to 1
				break;

		case 's':	if (atoi(argv[i+1])>1){
                                    gfrom = atoi(argv[i+1])-1;
                                }
                                else if ((atoi(argv[i+1])==0)||(atoi(argv[i+1])==1)){
                                    gfrom = 0;
                                }
                                else {
                                    DoError("Start parameter must be a positive integer", 'y');
                                }
				break;
                
                case 'e':	if (atoi(argv[i+1])>0){
                                    gto = atoi(argv[i+1]);
                                }
                                else{
                                    DoError("End parameter must be a positive integer", 'y');
                                }
				break;

		case 'r':	startseed = atoi(argv[i+1]);
				break;
                
                case 'h':	Usage();
                                exit(0);
				break;
                                
                case 'd':	j=i;
                                nodatafiles=0;

                                while (argv[j+1]&&argv[j+1][0]!='-'){
                                    strcpy(datafile[nodatafiles],argv[j+1]);
                                    j++;
                                    nodatafiles++;
                                }
				break;
                                
                case 'z':	gstopvalue=atof(argv[i+1]); //
                                if ((gCompattype=='L')&&((gstopvalue>1)||(gstopvalue<0))){
                                    DoError("Stop value must be between 0 and 1. 0 indicates stop when no more incompatibility in matrix", 'y');
                                }
				break;
                        
                default:        sprintf(string,"Unknown operator: %c", argv[i][1]);
                                DoError(string,'y');

	    }
	}
    }

//printf("\n%s\n", datafile[0]);

srand( startseed );

if (inputfile[0]!='\0'){
    executeok=doExecuteFile();
    if (executeok==1){
        interactive='y';
    }
}
else {
    printf("\nNo input file selected");
    interactive='y';
}

if (interactive=='y'){
    DoInteractive(executeok);
}

//if (inputfile[0]=='\0'){strcpy(inputfile,"/Users/simonharris/Desktop/concat_nuc+eocyte.nex");gPermutations=1000;}//{DoError("Expecting input file", 'y');}

if ((gto<gfrom+1)&&(gto!=-1)){
    DoError("End parameter must be greater than start parameter", 'y');
}



//extract number of permutations and check it adheres to program rules

if (gCompattype!='A'){
    if (argv[3]==NULL){
        DoError("Must define number of permutations for compatibility types C, L and N", 'y');
    }
}
else {gPermutations=0;}



/*if ((argc<=9)||(gCompattype!='L')){
    gCharset=0;
}
else{
    if (atof(argv[9])<=1){
        gCharset=atof(argv[9]);
    }
    else{
        DoError("Charset parameter must be between 0 and 1", 'y');
    }
}*/


//run test

tests();

//free allocated memory
for (i=0; i<gchars; i++){
    free(gMatrix[i]);
}
free(gMatrix);

for (i=0; i<gchars; i++){
    for (j=0; j<2; j++){
        free(compat[i][j]);
    }
    free(compat[i]);
}
free(compat);

for (i=0; i<gchars; i++){
    free(iterationincompats[i]);
}
free(iterationincompats);

free(included);
free(same);
free(actualincompat);
free(Teststatistic);
free(tot);
free(tot2);
free(greater);
free(type);
free(removedlastround);
free(charnum);

printf("\n\nAnalysis complete. Time taken: ");

c = clock();

secs=c/CLOCKS_PER_SEC;
mins=0;
hours=0;
days=0;

mins=secs/60;
hours=mins/60;
days=hours/24;
secs=secs%60;
mins=mins%60;
hours=hours%24;


/*if (secs>60){
    mins=floor(secs/60);
    secs=secs-(mins*60);
}
if (mins>60){
    hours=floor(mins/60);
    mins=mins-(hours*60);
}
if (hours>24){
    days=floor(hours/24);
    hours=hours-(days*24);
}*/

if (floor(days)>0){printf("%ld days ", days);}
{printf("%.2ld:%.2ld:%.2ld\n", hours, mins, secs);}

printf("\n");

fclose(flog);

return 1;
}

void DoInteractive(int executeok){

if (executeok==1){getfilename(executeok);printheading();}

printoptions();



//exit(1);


}

void printoptions(void){

char option;
char temp[500];

printf("\n\ni: Input file:\t\t\t");
if (inputfile[0]=='\0'){
    printf("None");
}
else {
    printf("%s\n\t\t\t\t%d taxa\n\t\t\t\t%d characters", inputfile, gtaxa, gchars);
    
}
/*printf("\no: Output file:\t\t\t");
if (outputfile[0]=='\0'){
    printf("None");
}
else {
    printf("%s", outputfile);
}*/
printf("\nt: Compatibility type:\t\t");
if (gCompattype=='A'){
    printf("Number of incompatibilities");
}
else if (gCompattype=='C'){
    printf("Coefficient of Character Starte Randomness (CCSR)");
}
else if (gCompattype=='L'){
    printf("LeQuesne Probability (LQP)");
}
else if (gCompattype=='N'){
    printf("Normal Deviate (NDev)");
}


if (gCompattype!='A'){
    printf("\np: Number of permutations:\t%d", gPermutations);
    if (gPermutations>0){
        printf("\nr: Starting seed:\t\t%ld", startseed);
    }
}
printf("\nf: Fuzzy compatibility:\t\t");
if (gFuzzy==0){
    printf("No");
}
else {
    printf("Yes");
}
printf("\nb: Boildown:\t\t\t");
if (gBoildown==0){
    printf("No");
}
else if (gBoildown==1){
    printf("Yes (slow)");
}
else {
    printf("Yes (fast)");
}
printf("\ns: Characters to analyse:\t");
if((gfrom==0)&&(gto==gchars)){
    printf("All");
}
else{
    printf("%d to %d", gfrom, gto);
}
printf("\nz: Stop value:\t\t\t%g",gstopvalue);
printf("\nm: Print mask file:\t\t");
if (mask=='n'){
    printf("No");
}
else {
    printf("Yes");
}
printf("\nc: Print character set file:\t");
if (charset=='n'){
    printf("No");
}
else {
    printf("Yes");
}

printf("\nh: Show help");
printf("\nq: Quit");


printf("\n\nChoose an option to change settings or type y to run analysis: ");
fflush(stdout);
option=getchar();
getchar();

if (tolower(option)=='q') exit(0);



//for (i = 1; i < argc; i++) {

	/* Check for a switch (leading "-"). */

//	if (argv[i][0] == '-') {

	    /* Use the next character to decide what to do. */

	    switch (tolower(option)) {
            
                case 'm':	if (mask=='y')mask='n'; 
                                else mask='y';
				break;
                
                case 'c':	if (charset=='y')charset='n';
                                else charset='y';//
				break;
                                
/*                case 'x':	cutoff=atof(argv[i+1]); //
                                if ((cutoff>1)||(cutoff<=0)){
                                    DoError("Print frequency must be a percentage of characters removed", 'y');
                                }
				break;*/

		case 'i':	getfilename(1); //read inputfile name and execute the file. Quit if file cannot be opened or executed
				break;
                                
/*                case 'o':       strcpy(outputfile,argv[i+1]); //read inputfile name and execute the file. Quit if file cannot be opened or executed
				break;
*/
		case 't':	if (gCompattype=='A') gCompattype='C'; //Set boildown to 0/1/2
                                else if (gCompattype=='C'){ gCompattype='L';
                                    if (gstopvalue>1)gstopvalue=0;
                                    }
                                else if (gCompattype=='L') gCompattype='N';
                                else gCompattype='A';
				break;

		case 'p':       if (gCompattype!='A'){
                                    gPermutations=0;
                                    while ((gPermutations<1)||(gPermutations>10000)){
                                        printf("\nEnter number of permutations: ");
                                        fflush(stdout);
                                        fgets(temp,500,stdin);
                                        if (temp[strlen(temp)-1] == '\n') temp[strlen(temp)-1] = '\0';
                                        gPermutations=atoi(temp);
                                
                                        if ((gPermutations<1)||(gPermutations>10000)){
                                            printf("Number of permutations must be an integer between 1 and 10,000");
                                        }
                                    }
                                }
                                else printf("Unknown choice: %c", option);
				break;
    
		case 'b':       if (gBoildown==0) gBoildown=1; //Set boildown to 0/1/2
                                else if (gBoildown==1) gBoildown=2;
                                else gBoildown=0;
				break;
            
                case 'f':	if (gFuzzy==0) gFuzzy=1; //Set fuzzy to 0/1
                                else gFuzzy=0;
				break;

		case 's':	gfrom=-1;
                                while ((gfrom<0)||(gfrom>(gchars-1))){
                                    printf("\nEnter site number to analyse from: ");
                                    fflush(stdout);
                                    fgets(temp,500,stdin);
                                    if (temp[strlen(temp)-1] == '\n') temp[strlen(temp)-1] = '\0';
                                    gfrom=atoi(temp);
                                
                                    if ((gfrom<0)||(gfrom>(gchars-1))){
                                        printf("start number must be an integer between 0 and %d", gchars-1);
                                    }
                                }
                                gto=-1;
                                while ((gto<(gfrom+1))||(gto>gchars)){
                                    printf("\nEnter site number to analyse to: ");
                                    fflush(stdout);
                                    fgets(temp,500,stdin);
                                    if (temp[strlen(temp)-1] == '\n') temp[strlen(temp)-1] = '\0';
                                    gto=atoi(temp);
                                
                                    if ((gto<(gfrom+1))||(gto>gchars)){
                                        printf("End number must be an integer between %d and %d", gfrom+1, gchars);
                                    }
                                }
				break;

		case 'r':	if (gCompattype!='A'){
                                    startseed=0;
                                    while ((startseed<1)||(startseed>1000000000)){
                                        printf("\nEnter starting seed: ");
                                        fflush(stdout);
                                        fgets(temp,500,stdin);
                                        if (temp[strlen(temp)-1] == '\n') temp[strlen(temp)-1] = '\0';
                                        startseed=atol(temp);
                                
                                        if ((startseed<1)||(startseed>1000000000)){
                                            printf("Starting seed must be an integer between 1 and 1,000,000,000");
                                        }
                                    }
                                }
                                else printf("Unknown choice: %c", option);
				break;
                
                case 'h':	Usage();
                                exit(0);
				break;
  /*                              
                case 'd':	j=i;
                                nodatafiles=0;

                                while (argv[j+1]&&argv[j+1][0]!='-'){
                                    strcpy(datafile[nodatafiles],argv[j+1]);
                                    j++;
                                    nodatafiles++;
                                }
				break;
   */                                
                case 'z':	gstopvalue=-1;
                                if (gCompattype=='L'){
                                while ((gstopvalue<0)||(gstopvalue>1)){
                                    printf("\nEnter value to stop analysis: ");
                                    fflush(stdout);
                                    fgets(temp,500,stdin);
                                    if (temp[strlen(temp)-1] == '\n') temp[strlen(temp)-1] = '\0';
                                    gstopvalue=atof(temp);
                                
                                    if ((gstopvalue<0)||(gstopvalue>1)){
                                        printf("Stopping value must be between 0 and 1 for LQP analysis");
                                    }
                                }
                                }
                                else {
                                while (gstopvalue<0){
                                    printf("\nEnter value to stop analysis: ");
                                    fflush(stdout);
                                    fgets(temp,500,stdin);
                                    if (temp[strlen(temp)-1] == '\n') temp[strlen(temp)-1] = '\0';
                                    gstopvalue=atof(temp);
                                
                                    if (gstopvalue<0){
                                        printf("Stopping value must be 0 or greater");
                                    }
                                }
                                
                                }
				break;
                     
                case 'y':       return;
                                break;
                default:        printf("Unknown choice: %c", option);

	    }

printheading();
printoptions();


}

void getfilename(int executeok){

while (executeok==1){
    printf("\n\nPlease enter the name of your input file (or q to quit): ");
    fflush(stdout);
    fgets(inputfile,500,stdin);
    if (inputfile[strlen(inputfile)-1] == '\n') inputfile[strlen(inputfile)-1] = '\0';
    if (tolower(inputfile[0])=='q'){
        exit(0);
    }
    if (inputfile[0]!='\0'){
        executeok=doExecuteFile();
    }
}

}


void printheading(void){

system("clear");
    
printf("\nCOMPASS (Compatibility Analysis Site Stripping)\n\nv1.0 (2009) written by Simon Harris, The Wellcome Trust Sanger Institute, Hinxton, UK\n\n");

printf(asctime(localtime (&startseed)));

}
