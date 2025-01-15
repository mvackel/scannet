/* 
    O SCANNET visa simplificar a utilização do método PSN (Protein Similarity Networks)
	a partir da unificação dos seus passos.
	
    Copyright (C) 2015 Austeclino Magalhães Barros Júnior

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	Modified: 2025-01 by Marcos Ackel
 */
 
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h> 
#include <time.h>
#include <stdbool.h>

#ifdef WINDOWS
    #include <direct.h>
    #define GetCurrentDir _getcwd
#else
    #include <unistd.h>
    #define GetCurrentDir getcwd
#endif

typedef struct CN{
    int** adjacencyMatrix;
    int** neighborhoodMatrix;
    int** colorsMatrix;
    float** dendrogramGraphic;
    int numberOfCommunities;
    int qtdOfEdges;
    int numberRemovedEdges;
    int numberDendrogramLines;
    int* realIndex;
    int* realIndexDendrogramGraphic;
} complexNetwork;

typedef struct C_EDGE{
    int V1;
    int V2;
}edge;

// struct for matrix indexing
typedef struct {
	unsigned int row;
	unsigned int col;
} index_t;

// Prototypes
void loadSimilarityMatrix();
void loadSequenceList(char[], char**);
void createComplexNetworks();
complexNetwork createComplexNetwork();
void executeNga(complexNetwork);
void dijkstra (int, int**, int**, int);
void processSimilarityMatrix(int**);
void createAdjacencyMatrix(int, int);
int** createDynamicMatrix(int, int);
float** createDynamicMatrixFloat(int, int);
int* createDynamicArray(int);
void createNeighborhoodMatrix(complexNetwork);
void DestroyDynamicMatrix(int**, int );
void clearMatrix(int**);
void clearMatrixBySize(int**, int );
bool isThereEdge(int**);
void calculateBetweenness(int, complexNetwork, float**);
void removeEdgesWithBiggerBetweenness(index_t[], float**, int**);
void processGraphicDendrogram(index_t[],complexNetwork*, int**);
void startDendrogram(complexNetwork*);
void copyMatrixOfFloat(int**, float**);
void copyMatrix(int**, int**);
void showRemovedEdges(unsigned long totalRemoved, time_t start);
void conveterArrayInMatrix(int*, int**, int, int);
void conveterMatrixInArray(int**, int*, int , int );
int findInArray(int*, int, int);
void sortGraphicDendrograma(complexNetwork*);
void createColorsMatrix(complexNetwork*);
void saveMatrizInFile();
void saveMatrizInFileFloat();
int findRealIndex(int* realIndex, int);
void adjustNeighborhoodMatrix(index_t[], int**, int**, int);
int countNumberLinesOfFile(char[30]);
void plotGraphic();
void setupOutputNames(char*, int, const char*); //char[10]);
void saveNetwork(int**, char[30]);
void orderSequences(complexNetwork*);
void saveArrayInFile(char**, char[], int);
void saveGraphicDendrogramInFile(float*);
char** createDynamicCharMatrix(int);
complexNetwork selectBetterNetworkToFindCommunities(FILE*);
void calculateDistances(int);
int calculateDiameter();
void showtime();
void showElapsed(const char* text, const time_t start);		// show elapsed time from 'start' in format HH:mm:ss
float calculateDistance();
void showRunningClock(void);

//Global vars

int sizeSimilarityMatrix = 1;
float** similarityMatrix;
int NCutPoints = 101;
int N_EDGES = 100;
int indexMelhorRede = 0;
int selectedNetwork = 0;; 
float biggerDistance = 0;
float distances[101];
int i_distance = 0;
int qtdDendrogramLines = 0;
complexNetwork *complexNetworks;
edge *removedEdges;
int qtdNodeWithEdgeRemoved =0;
char current_path_entrada[1000];
char current_path_saida[1000];
int rank;
int size;
int optionLoadSequenceFile = 0;
char** sequenceArray;

time_t start_t_RC, end_t_RC, tempo_total_RC, start_t_Distancia, end_t_Distancia, 
       tempo_total_Distancia, start_t_NGA, end_t_NGA, tempo_total_NGA, start_t_MC, 
       end_t_MC, tempo_total_MC, start_t_TOTAL, end_t_TOTAL, tempo_total_TOTAL ;

// Not used:
//double start_t_Paralelo, end_t_Paralelo, total_t_Paralelo;
//double * tempos_paralelo;
char filename[31];
char sequenceFilename[21];
FILE *arqLog;
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------

int main(int argc, char** argv) {
        complexNetwork betterNetwork;
        int option = 1;       
		
		printf(" Copyright (C) 2015 Austeclino Magalhaes Barros Junior");
		printf("\n\n This program is free software: you can redistribute it and/or modify");
		printf("\n it under the terms of the GNU General Public License as published by");
		printf("\n the Free Software Foundation, either version 3 of the License, or");
		printf("\n (at your option) any later version.");
     	printf("\n\n This program is distributed in the hope that it will be useful,");
		printf("\n but WITHOUT ANY WARRANTY; without even the implied warranty of");
		printf("\n MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the");
		printf("\n GNU General Public License for more details.");
		printf("\n\n You should have received a copy of the GNU General Public License");
		printf("\n along with this program.  If not, see <http://www.gnu.org/licenses/>.");
		
		printf("\n\n\n --------------------------- SCANNET v1.1 ------------------------------ \n\n");
		 
        while (option == 1 || option == 2){
            
            if(option == 1){
                
                printf("\n\n Enter the similarity file name (Ex: input.txt):  ");
                scanf("%s", filename);

                sizeSimilarityMatrix = countNumberLinesOfFile(filename);
				getcwd(current_path_entrada, 255);
               
                #ifdef WINDOWS
					strcat(current_path_entrada, "\\");
				#else
					strcat(current_path_entrada, "/");
				#endif
				
                similarityMatrix = createDynamicMatrixFloat(sizeSimilarityMatrix, sizeSimilarityMatrix);
                loadSimilarityMatrix(filename, similarityMatrix);   
                start_t_TOTAL = time(NULL);
                optionLoadSequenceFile = 0;
				
                while((optionLoadSequenceFile != 1)&&(optionLoadSequenceFile != 2)){
                    printf("\n\n Do you want to load the initial order of proteins (1-Yes/2-No):  ");
                    scanf("%d", &optionLoadSequenceFile);

                    if(optionLoadSequenceFile == 1){
                         printf("\n\n Enter the order of protein file name (Ex: initialOrder.csv):  ");
                         scanf("%s", sequenceFilename);
                         sequenceArray = createDynamicCharMatrix(sizeSimilarityMatrix);
                         loadSequenceList(sequenceFilename, sequenceArray);   
                    }else if((optionLoadSequenceFile != 1)&&(optionLoadSequenceFile != 2)){
                         printf("\n Invalid Option  ");
						 continue;
                    }                 
                }
                
				char path[1000];
                strcpy(path,current_path_entrada);
                char log[30] = "log_";
                setupOutputNames(log, 0, ".txt");
                strcat(path,log);
                arqLog = fopen(path, "wt");
                fprintf(arqLog,"\n\n ------- Summary execution --------- \n");
				
                start_t_RC = time(NULL);
                //printf("\n\n The complex networks are being created, please wait... \n");
                createComplexNetworks();  // time consuming function
		        end_t_RC = time(NULL);
                tempo_total_RC = end_t_RC - start_t_RC;
            }
			betterNetwork = selectBetterNetworkToFindCommunities(arqLog);  
		    printf("\n Number of edges: %d ", betterNetwork.qtdOfEdges);
			int i = 0;
            for (i = 0; i < sizeSimilarityMatrix; i++){ 
                betterNetwork.realIndex[i] = i;
                betterNetwork.realIndexDendrogramGraphic[i] = i;
            }
			char networkName[30] = "network_";
            setupOutputNames(networkName, 1, ".csv");
            saveNetwork(betterNetwork.neighborhoodMatrix, networkName);
	        copyMatrix(betterNetwork.neighborhoodMatrix, betterNetwork.colorsMatrix);
            printf("\n\n The Girvan-Newman algorithm is running, please wait ...");
			showtime("Start NGA");
			printf("\n");
            executeNga(betterNetwork);			// very time consuming function
            end_t_TOTAL = time(NULL);

			for (i = 0; i < sizeSimilarityMatrix; i++){ 
               betterNetwork.dendrogramGraphic[1][i] = betterNetwork.realIndexDendrogramGraphic[i]+1;     
            }
			saveGraphicDendrogramInFile(betterNetwork.dendrogramGraphic[1]);
            
            char mColorName[30] = "mColors_";
            setupOutputNames(mColorName, 1, ".dat");
            saveMatrizInFile(betterNetwork.colorsMatrix, mColorName, sizeSimilarityMatrix, sizeSimilarityMatrix );

            if(optionLoadSequenceFile == 1){
                char sequenceOutput[30] = "orderedSequence_";
                setupOutputNames(sequenceOutput, 1, ".csv");
                saveArrayInFile(sequenceArray, sequenceOutput,sizeSimilarityMatrix);
            }
            
			tempo_total_TOTAL = tempo_total_RC + tempo_total_NGA;
            showtime("End of execution");  
            fprintf(arqLog,"\n\n The selected optimal network has %d edges", betterNetwork.qtdOfEdges);
            fprintf(arqLog,"\n Time spent on the creation of 101 complex networks (s): %ld", tempo_total_RC); 
            fprintf(arqLog,"\n Time spent on NGA (s): %ld", tempo_total_NGA);
            fprintf(arqLog,"\n Running time (s): %ld \n", tempo_total_TOTAL);
            fprintf(arqLog,"\n --------------------------------------------------- \n\n\n");
            
            printf("\n\n What do you want to do now?  ");
            printf("\n\n 1 - Analyze new sequence  ");
            printf("\n 2 - Analyze this sequence with other similarity value ");
            printf("\n 3 - Exit \n");
            scanf("%d", &option);
            if(option == 1 || option == 3){
				 free(sequenceArray);
                 free(similarityMatrix);
				 free(complexNetworks);
                 fclose(arqLog);
            }
        } 
	
        return 0;
}
//-----------------------------------------------------------------------------

void saveNetwork(int** matrix , char networkName[30]){
      char path[1000];
      int i, j;
	  FILE *arqNetwork;
      strcpy(path,current_path_entrada);
      strcat(path,networkName);
      arqNetwork = fopen(path, "wt");
      fprintf(arqNetwork,"source;target;type\n");
      
       for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for(j = (i+1); j < sizeSimilarityMatrix; j++){
           if(matrix[i][j] == 1){
               fprintf(arqNetwork,"%d;%d;undirected\n",i+1,j+1);
           } 
        } 
    }

    fclose(arqNetwork);
    printf("\n\n The file %s was successfully saved! \n", networkName);
}
//-----------------------------------------------------------------------------

void setupOutputNames(char * nameOfFile, int concatNetwork, const char *format){ //[10] ){
    
    int sizeFileName = strlen (filename);
    char fileNameWithoutExtension[30] = "";
    char buffer [33];
    
    int i = 0;
    for (i = 0; i < sizeFileName; i++){
        if(filename[i] == '.'){
            break;
        }else{
            sprintf(buffer, "%c", filename[i]);
            strcat(fileNameWithoutExtension, buffer);
        }
    }
    strcat(nameOfFile, fileNameWithoutExtension); 
    
    if(concatNetwork == 1){
        strcat(nameOfFile, "_");
        char index[30];
		char timeString[30];
        sprintf(index, "%d", selectedNetwork);
	    strcat(nameOfFile, index);
    }
    strcat(nameOfFile, format);
}
//-----------------------------------------------------------------------------

int countNumberLinesOfFile(char nameOfFile[30]){   
    FILE *fileptr;
    int numberLines = 1;
    int chr;
	int chrBefore;
    fileptr = fopen(nameOfFile, "r");
    chr = getc(fileptr);
    while ((chr != EOF) && (chr != -1)){
        if (chr == '\n' && chrBefore != '\n'){
            numberLines = numberLines + 1;
        }
		chrBefore = chr;
        chr = getc(fileptr);
    }
	if (chrBefore == '\n'){
		numberLines = numberLines - 1;
    }
    fclose(fileptr);
    return numberLines;
}
//-----------------------------------------------------------------------------

void conveterArrayInMatrix(int* array, int** matrix, int n, int m){
    int i = 0;  
    int j = 0;
    int index = 0;
    int size = n*m;
    for (index = 0; index < size; index++){ 
        matrix[i][j] = array[index];
        
        if(j == (m-1) ){
            j = 0;
            i++;
        }else{
            j++;
        } 
    }
}
//-----------------------------------------------------------------------------

void conveterMatrixInArray(int** matrix, int* array, int n, int m){
    int i = 0;  
    int j = 0;
    int index = 0;
    
    for (i = 0; i < n; i++){ 
        for (j = 0; j < m; j++){ 
            array[index] = matrix[i][j];
            index++;
        }
    }
}
//-----------------------------------------------------------------------------

void showtime(char texto[30]){
    time_t currentTime;
    struct tm *timeinfo;
    currentTime = time(NULL);
    timeinfo = localtime(&currentTime);
    printf("\n %s: %02d:%02d:%02d \n",texto, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec);   
}
//-----------------------------------------------------------------------------

// display the text and elapsed time from start in the format HH:mm:ss
// Doesn't include any newline.
void showElapsed(const char *text, const time_t start)
{
	time_t now = time(NULL);
	time_t elapsed = now - start;
	printf("%s %02ld:%02ld:%02ld", text, elapsed/3600, (elapsed/60)%60, elapsed%60);
}
//-----------------------------------------------------------------------------

int** createDynamicMatrix(int sizeLine, int sizeColumn){
    int** matrix = malloc(sizeLine * sizeof(int*));    
    int i;
    if(matrix == NULL){
         printf("\n  Insufficient memory. \n");
		 system("pause");
         exit(1);
    }
    for(i=0; i < sizeLine; i++){
        matrix[i] = calloc(sizeColumn, sizeof(int));
        if(matrix[i] == NULL){
            printf("\n Insufficient memory. Total = %d, Linha = %d \n", sizeLine, i);
			system("pause");
            exit(1);
        }
    }
    
    return matrix;
}
//-----------------------------------------------------------------------------

int* createDynamicArray(int size){
    int* matrix = calloc(size, sizeof(int));     
  
    if(matrix == NULL){
         printf("\n Insufficient memory. \n");
         exit(1);
    }
    
    return matrix;
}
//-----------------------------------------------------------------------------

char** createDynamicCharMatrix(int size){
    char** matrix = malloc(size*sizeof(char*));
    int i;
  
    if(matrix == NULL){
         printf("\n Insufficient memory. \n");
         exit(1);
    }
        
    for(i=0; i < size; i++){
        matrix[i] = calloc(1000, sizeof(int));
        if(matrix[i] == NULL){
            printf("\n Insufficient memory. \n");
            exit(1);
        }
    }
    
    return matrix;
}
//-----------------------------------------------------------------------------

void reiniciaMatriz(int sizeLine, int sizeColumn, int** matrix){
    int i, j;
    
    for(i=0; i < sizeLine; i++){
        for(j=0; j < sizeColumn; j++){
            matrix[i][j] = 0;
        }
    }
}
//-----------------------------------------------------------------------------

void createAdAndNM(int sizeLine, int sizeColumn, int** adjacencyMatrix, int ** neighborhoodMatrix){
    int i;
    
    for(i=0; i < sizeLine; i++){
        adjacencyMatrix[i] = calloc(sizeColumn, sizeof(int));
        neighborhoodMatrix[i] = calloc(sizeColumn, sizeof(int));
    }
}
//-----------------------------------------------------------------------------

float** createDynamicMatrixFloat(int sizeLine, int sizeColumn){

	float ** matrix = malloc(sizeLine * sizeof(float*));    
    int i;

    if(matrix == NULL){
         printf("\n Insufficient memory. \n");
		 system("pause");
         exit(1);
    }

    for(i=0; i < sizeLine; i++){
        matrix[i] = calloc(sizeColumn, sizeof(float));
         if(matrix[i] == NULL){
			printf("\n Insufficient memory. Total = %d, Linha = %d \n", sizeLine, i);
			system("pause");
            exit(1);
         }
    }

    return matrix;
}
//-----------------------------------------------------------------------------

void DestroyDynamicMatrix(int** matrix, int sizeLine){
    int i;
    
    for(i=0; i < sizeLine; i++){
        free(matrix[i]); 
    }
    free(matrix);
}
//-----------------------------------------------------------------------------

void DestroyDynamicMatrixFloat(float** matrix, int sizeLine){
    int i;
    
    for(i=0; i < sizeLine; i++){
        free(matrix[i]); 
    }
    free(matrix);
}
//-----------------------------------------------------------------------------

void loadSimilarityMatrix(char nameOfFile[30], float** matrix){
    FILE *arqBLAST; 
    char *endLinha; 
    int size = sizeSimilarityMatrix * 4;
    char linha[size];

    char *token = NULL;
    float n = 0.0;
    int indexLine = 0;
    int indexCollumn = 0;
  
    char path[1000];
    strcpy(path,current_path_entrada);  
    strcat(path,nameOfFile);
    
    arqBLAST = fopen(path, "rt");
    
    if (arqBLAST == NULL){
        printf("\n\n\n The file with the similarity matrix can not be opened \n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    while (indexLine < sizeSimilarityMatrix) {
        endLinha = fgets(linha, size, arqBLAST);
        if (endLinha){ 
            token= strtok(linha, " " );
            n = atof(token);
            matrix[indexLine][indexCollumn] = n;
            indexCollumn++;
            
            while (indexCollumn < sizeSimilarityMatrix) {
                 token= strtok(NULL, " " ); 
                 n = atof(token);
                 matrix[indexLine][indexCollumn] = n;
                 indexCollumn++;
            }             
        }
        indexCollumn = 0;
        indexLine++;
    }
         
   fclose(arqBLAST);
} 
//-----------------------------------------------------------------------------

void loadSequenceList(char nameOfFile[30], char** array){
    FILE *arqSequence; 
    char ch;
    int n = 0;
    int indexLine = 0;
  
    char path[1000];
    strcpy(path,current_path_entrada);  
    strcat(path,nameOfFile);
    
    arqSequence = fopen(path, "rt");
    
    if (arqSequence == NULL){
        printf("\n\n\n The file with the sequence list can not be opened \n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    while (indexLine < sizeSimilarityMatrix) {
        ch = getc(arqSequence);
       
	if(ch == '\n'){
            indexLine++;
            n = 0;  
        }else{
            array[indexLine][n] = ch;
            n++;
        }       
    }
         
   fclose(arqSequence);
} 
//-----------------------------------------------------------------------------

void saveMatrizInFile(int** matriz, char nameOfFile[30], int n, int m){         
    FILE *arqMATRIZ;
    int w,k,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
   
    if (arqMATRIZ == NULL){
        printf("\n Problems in the creation of the similarity matrix file.\n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    for (w = 0; w < n; w++){
        for (k = 0; k < m; k++){
            if(matriz[w][k] < 10){
                result = fprintf(arqMATRIZ,"%d   ", matriz[w][k]);
            }else{
                result = fprintf(arqMATRIZ,"%d  ", matriz[w][k]);
            }
            if (result == EOF) printf("\n Error in file recording\n");
        }

        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("\n Error in the recording of the similarity matrix file \n");
    }

    fclose(arqMATRIZ);

    printf("\n\n\n The file %s was successfully saved! \n", nameOfFile);
}
//-----------------------------------------------------------------------------
    
void saveMatrizInFileFloat(float** matriz, char nameOfFile[30], int n, int m){         
    FILE *arqMATRIZ;
    int w,k,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
    
    if (arqMATRIZ == NULL){ 
        printf("\n Problems in the creation of the similarity matrix file.\n");
        system("pause");
        exit (EXIT_FAILURE);
    } 

    for (w = 0; w < n; w++){
        for (k = 0; k < m; k++){
            result = fprintf(arqMATRIZ,"%.2f ", matriz[w][k]);
            if (result == EOF) printf("Error in file recording\n");
        }

        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("\n Error in the recording of the similarity matrix file\n");
    }

    fclose(arqMATRIZ);

    printf("\n\n The file %s was successfully saved! \n", nameOfFile);
}
//-----------------------------------------------------------------------------

void saveArrayInFileFloat(float* array, char nameOfFile[30], int size){         
    FILE *arqMATRIZ;
    int w,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "wt");
   
    if (arqMATRIZ == NULL){
        printf("\n Problems in the creation of the distances file \n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    for (w = 0; w < size; w++){
        result = fprintf(arqMATRIZ,"%.2f", array[w]);
        result = fprintf(arqMATRIZ,"\n");
        if (result == EOF) printf("\n Error in file recording\n");
    }

    fclose(arqMATRIZ);

    printf("\n\n\n The file %s was successfully saved! \n", nameOfFile);
}
//-----------------------------------------------------------------------------

void saveGraphicDendrogramInFile(float* array){     
    char nameOfFile[50] = "mDendrogramGraphic_";
    setupOutputNames(nameOfFile, 1, ".dat");
	int size = sizeSimilarityMatrix;

    FILE *arqMATRIZ;
    int w,result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqMATRIZ = fopen(path, "a");
   
    if (arqMATRIZ == NULL){
        printf("\n Problems in the creation of the graphic dendrogram file \n");
        system("pause");
        exit (EXIT_FAILURE);
    }
	
	result = fprintf(arqMATRIZ,"%.2f ", array[sizeSimilarityMatrix]);
    if (result == EOF) printf("\n Error in file recording\n");
    
	for (w = 0; w < size; w++){
        result = fprintf(arqMATRIZ,"%.2f ", array[w]);
        if (result == EOF) printf("\n Error in file recording\n");
    }
	
    result = fprintf(arqMATRIZ,"\n");
	if (result == EOF) printf("\n Error in file recording\n");
    fclose(arqMATRIZ);
}
//-----------------------------------------------------------------------------

void saveArrayInFile(char** array, char nameOfFile[30], int size){         
    FILE *arqSequence;
    int i, result = 0;
    char path[1000];
    strcpy(path,current_path_entrada);
    strcat(path,nameOfFile);
    arqSequence = fopen(path, "wt");
   
    if (arqSequence == NULL){
        printf("\n Problems in the creation of the sequence file \n");
        system("pause");
        exit (EXIT_FAILURE);
    }

    for (i = 0; i < size; i++){
        for (i = 0; i < size; i++){
            result = fprintf(arqSequence,"%s", array[i]);
            result = fprintf(arqSequence,"\n");
            if (result == EOF) printf("\n Error in file recording \n");
        }
    }

    fclose(arqSequence);
    printf("\n\n\n The file %s was successfully saved! \n", nameOfFile);
}
//-----------------------------------------------------------------------------

void createComplexNetworks(){
    complexNetworks = (complexNetwork *)malloc(2 * sizeof(complexNetwork));

    if(complexNetworks == NULL){
         printf("\n Error in memory allocation. \n");
         exit(1);
    }
    
	int i = 0;
	complexNetworks[0] = createComplexNetwork();
	createAdjacencyMatrix(i,0);
	createNeighborhoodMatrix(complexNetworks[0]);
    
	printf("\n\n");
    for (i = 1;  i < NCutPoints; i++){
		printf("\r >>> Creating Complex Networks: %d/%d", i+1, NCutPoints);
		fflush(NULL);
	    complexNetworks[1] = createComplexNetwork();;
        createAdjacencyMatrix(i,1);
        createNeighborhoodMatrix(complexNetworks[1]);
		calculateDistances(i-1);
		complexNetworks[0] = complexNetworks[1];
    }
	
	char mDistanceGraphic[30] = "mDistanceGraphic_";
    setupOutputNames(mDistanceGraphic, 0, ".dat");
    saveArrayInFileFloat(distances, mDistanceGraphic ,i_distance);
}	
//-----------------------------------------------------------------------------

complexNetwork createComplexNetwork(){
    complexNetwork complexNetwork; 
	int i = 0;
	complexNetwork.adjacencyMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);
	complexNetwork.neighborhoodMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);
	complexNetwork.numberOfCommunities = 1;
	complexNetwork.numberRemovedEdges = 0;
	complexNetwork.numberDendrogramLines = 0;
	complexNetwork.qtdOfEdges = 0;
	
    return complexNetwork;
}   
//-----------------------------------------------------------------------------

void createAdjacencyMatrix(int cutPoint, int index){
    int i,j;   
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for(j = (i+1); j < sizeSimilarityMatrix; j++){
			if(similarityMatrix[i][j] >= cutPoint){
		        complexNetworks[index].adjacencyMatrix[i][j] = 1;
                complexNetworks[index].adjacencyMatrix[j][i] = 1;
                complexNetworks[index].qtdOfEdges++;
            }
        } 
    }
}
//-----------------------------------------------------------------------------

void createNeighborhoodMatrix(complexNetwork complexNetwork){
    int i;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){
         dijkstra(i, complexNetwork.adjacencyMatrix, complexNetwork.neighborhoodMatrix, sizeSimilarityMatrix);
    }
}
//-----------------------------------------------------------------------------

void dijkstra (int Vi, int** adjacencyMatrix, int** neighborhoodMatrix, int maxSize)
{       
    int dis[maxSize];
    char vis[maxSize];
    memset (vis, 0, sizeof (vis));
    memset (dis, 0x7f, sizeof (dis));
    dis[Vi] = 0;

    while (1){
        int i, n = -1;
        for (i = 0; i < maxSize; i++){
            if (!vis[i] && (n < 0 || dis[i] < dis[n])){
                    n = i;
            }        
        }

        if (n < 0){
           break;
        }

        vis[n] = 1;
        for (i = 0; i < maxSize; i++){
            if (adjacencyMatrix[n][i] && dis[i] > dis[n] + adjacencyMatrix[n][i]){
                dis[i] = dis[n] + adjacencyMatrix[n][i];
                neighborhoodMatrix[Vi][i] = dis[i];
                neighborhoodMatrix[i][Vi] = dis[i];
            }        
        }        
    }
}
//-----------------------------------------------------------------------------

void calculateDistances(int i){
	  float distance = 0;
      complexNetwork betterNetwork;
	  distance =  calculateDistance(complexNetworks[0], complexNetworks[1]);  
	  distances[i_distance] = distance;
	  i_distance = i_distance + 1;
	  if(distance > biggerDistance){
		  biggerDistance = distance;
		  indexMelhorRede = i;
	  }
}
//-----------------------------------------------------------------------------

complexNetwork selectBetterNetworkToFindCommunities(FILE* arqLog){
	int i = 0;     
    plotGraphic(); 
	for (i = 0; i < NCutPoints-1; i++){
		 fprintf(arqLog,"\n Distance between the network %2d and %2d = %f", i, i+1, distances[i]);
         printf("\n Distance between the network %d and %d = %f", i, i+1, distances[i]);
	}		
      
	fprintf(arqLog, "\n\n The similarity value %d has the largest distance: %f  ", indexMelhorRede,  biggerDistance);
	printf("\n\n The similarity value %d has the largest distance: %f  ", indexMelhorRede,  biggerDistance);  
    printf("\n\n\n Select the similarity value you want to use:  ");
    scanf("%d", &selectedNetwork);
      
    fprintf(arqLog, "\n\n Similarity value selected: %d  ", selectedNetwork);
	complexNetworks[1] = createComplexNetwork();
	createAdjacencyMatrix(selectedNetwork, 1);
    createNeighborhoodMatrix(complexNetworks[1]);  
           
	complexNetworks[1].colorsMatrix = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);      
	complexNetworks[1].dendrogramGraphic = createDynamicMatrixFloat(2, sizeSimilarityMatrix+1);
	complexNetworks[1].realIndex = createDynamicArray(sizeSimilarityMatrix);
	complexNetworks[1].realIndexDendrogramGraphic = createDynamicArray(sizeSimilarityMatrix);
	
    return complexNetworks[1];
}
//-----------------------------------------------------------------------------

void plotGraphic(){

    char path[1024];
    char syscommand[1024];
    strcpy(path,current_path_entrada);
    char mDistanceGraphic[100] = "mDistanceGraphic_";
    setupOutputNames(mDistanceGraphic, 0, ".dat");
    strcat(path,mDistanceGraphic);
	
    #ifdef WIN32
        sprintf(syscommand, "start gnuplot -persist -e \"plot '%s' with lines title 'SCANNET - Distances chart'; set xlabel 'Similarity Value'; set ylabel 'Distance'; set xtics 2; set ytics 0.2; pause -1 'Press enter to exit. . .';\"\"",path);
	#else
		sprintf(syscommand, "gnuplot -p -e \"plot '%s' with lines title 'SCANNET - Distances chart';set xlabel 'Similarity Value'; set ylabel 'Distance'; set xtics 2; set ytics 0.2;\"",path);
	#endif
	
    system(syscommand);
}
//-----------------------------------------------------------------------------

float calculateDistance(complexNetwork networkOne, complexNetwork networkTwo){
    int i = 0;
    int j = 0;
    float t1 = 1;
    float t2 = (sizeSimilarityMatrix * sizeSimilarityMatrix);
    float a = t1/t2;
    float r = 0;
    float total = 0;
     
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = (i+1); j < sizeSimilarityMatrix; j++){ 
            r = (networkOne.neighborhoodMatrix[i][j]) - (networkTwo.neighborhoodMatrix[i][j]);
            r = r * r;
            total = total + r;        
        }
    }
    
    return a * total;
}
//-----------------------------------------------------------------------------

int calculateDiameter(complexNetwork network){
    int i = 0;
    int j = 0;
    int bigger = 0;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
         for (j = (i+1); j < sizeSimilarityMatrix; j++){ 
             if(network.neighborhoodMatrix[i][j] > bigger){
                 bigger = network.neighborhoodMatrix[i][j];
             }   
         }
    }

    return bigger;
}
//-----------------------------------------------------------------------------

void executeNga(complexNetwork complexNetwork){
    start_t_NGA = time(NULL); //getTime();
    
    float** betweennessMatrix = createDynamicMatrixFloat(sizeSimilarityMatrix, sizeSimilarityMatrix);
    int** neighborhoodMatrixPrevious = createDynamicMatrix(sizeSimilarityMatrix, sizeSimilarityMatrix);
    int diameter = 0, i = 1;
    int parcialArestasRemovidas = 0, totalArestasRemovidas = 0; 
    startDendrogram(&complexNetwork);

    while( isThereEdge(complexNetwork.adjacencyMatrix) ){
        diameter = calculateDiameter(complexNetwork);
        calculateBetweenness(diameter, complexNetwork, betweennessMatrix); 
        //int nodeWithEdgeRemoved[i][2];
        index_t nodeWithEdgeRemoved[i];
        qtdNodeWithEdgeRemoved = 0;
        copyMatrix(complexNetwork.neighborhoodMatrix ,neighborhoodMatrixPrevious);
        removeEdgesWithBiggerBetweenness(nodeWithEdgeRemoved, betweennessMatrix, complexNetwork.adjacencyMatrix); 		
        complexNetwork.numberRemovedEdges++;  
        adjustNeighborhoodMatrix(nodeWithEdgeRemoved, complexNetwork.adjacencyMatrix, complexNetwork.neighborhoodMatrix, diameter);
        processGraphicDendrogram(nodeWithEdgeRemoved, &complexNetwork, neighborhoodMatrixPrevious); 		

		parcialArestasRemovidas += qtdNodeWithEdgeRemoved;
    	totalArestasRemovidas += qtdNodeWithEdgeRemoved;

		showRunningClock();
		if (parcialArestasRemovidas >= N_EDGES) {
			showRemovedEdges(totalArestasRemovidas, start_t_NGA);
			parcialArestasRemovidas = 0;
		}
    }
	printf("\n");
	saveGraphicDendrogramInFile(complexNetwork.dendrogramGraphic[complexNetwork.numberDendrogramLines-1]);

    end_t_NGA = time(NULL); //getTime();
    tempo_total_NGA = end_t_NGA - start_t_NGA;

    start_t_MC = time(NULL);
    createColorsMatrix(&complexNetwork);
    end_t_MC = time(NULL);
    tempo_total_MC = end_t_MC - start_t_MC;
    
    DestroyDynamicMatrixFloat(betweennessMatrix, sizeSimilarityMatrix);
    DestroyDynamicMatrix(neighborhoodMatrixPrevious, sizeSimilarityMatrix);
    
    qtdDendrogramLines = complexNetwork.numberDendrogramLines;
    if(optionLoadSequenceFile == 1){
        orderSequences(&complexNetwork);
    }
}
//-----------------------------------------------------------------------------

void startDendrogram(complexNetwork *complexNetwork){
    int i,j, index = 0;
    int index2 = 0;
    int* newCommunityItens = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix); 
    int* newCommunityItens2 = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);  
    int pontoMaximo = 0;
    float dividendo = 2.0;
    int realIndex = 0;
    int hasNewCommunity = 0;
    int repeatLine = 1;
	int t=0;
    for (i = 0; i < sizeSimilarityMatrix; i++){
        (*complexNetwork).dendrogramGraphic[0][i] = (sizeSimilarityMatrix+1) / dividendo;
    }
    (*complexNetwork).numberDendrogramLines = 1;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){
        if(findInArray(newCommunityItens2, i, index2) == -1){
            index = 0;
            newCommunityItens[index] = i;
            newCommunityItens2[index2] = i;
            realIndex = findRealIndex((*complexNetwork).realIndex, i);
            index2++;
            index++;
                           
             for (j = 0; j < sizeSimilarityMatrix; j++){
                 if( (*complexNetwork).neighborhoodMatrix[i][j] != 0 ){
                      hasNewCommunity = 1;
                      if( (findInArray(newCommunityItens, j, index) == -1) ){
                           realIndex = findRealIndex((*complexNetwork).realIndex, j);
                           newCommunityItens[index] = j;
                           index++;
                           newCommunityItens2[index2] = j;
                           index2++;
                      }
                 }
            }
            
            if(hasNewCommunity && repeatLine){
                repeatLine = 0;
                for (j = 0; j < sizeSimilarityMatrix+1; j++){
                    (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][j] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines - 1][j];
                }
				saveGraphicDendrogramInFile((*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1]);
            }
            
            for (j = 0; j < index; j++){
               realIndex = findRealIndex((*complexNetwork).realIndexDendrogramGraphic, newCommunityItens[j]);
               (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][realIndex] = pontoMaximo + ((index+1) / dividendo);  
            }
         
            (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][sizeSimilarityMatrix] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1][sizeSimilarityMatrix] + 1;
         	
            sortGraphicDendrograma(complexNetwork);
            
            pontoMaximo = pontoMaximo + index;
 
            (*complexNetwork).numberOfCommunities++;   			
        } 
    }
	
   	saveGraphicDendrogramInFile((*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1]);
	for (j = 0; j < sizeSimilarityMatrix+1; j++){
      (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines - 1][j] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][j];
    }
	
    (*complexNetwork).numberOfCommunities--;
    free(newCommunityItens);
    free(newCommunityItens2);
}
//-----------------------------------------------------------------------------

int findRealIndex(int* realIndex, int virtualIndex){
    int j;
    for (j = 0; j < sizeSimilarityMatrix; j++){
        if(realIndex[j] == virtualIndex ){
            return j;
        }
    }  
    return -1;
}
//-----------------------------------------------------------------------------

void processGraphicDendrogram(index_t nodeWithEdgeRemoved[1], complexNetwork *complexNetwork, int** neighborhoodMatrixPrevious){
    int i,j, index = 0, index2 = 0, index3 = 0;
    int* newCommunityItens = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);
    int* itensOldCommunity = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);  
    int* AuxCommunity = createDynamicArray(sizeSimilarityMatrix*sizeSimilarityMatrix);  
    
    int hasNewCommunity = 0;
    int finded1 = 0; 
    int finded2 = 0; 
    
    int indiceJ = nodeWithEdgeRemoved[0].row;
    int indiceI = nodeWithEdgeRemoved[0].col;
    
    float dividendo = 2.0;
    int realIndex = 0;

    if( (*complexNetwork).neighborhoodMatrix[indiceI][indiceJ] == 0) {
        hasNewCommunity = 1;
        newCommunityItens[index] = indiceJ;
        index++;
        itensOldCommunity[index2] = indiceI;
        index2++; 
        for (j = 0; j < sizeSimilarityMatrix; j++){
            if(((*complexNetwork).neighborhoodMatrix[indiceI][j] == 0) && (neighborhoodMatrixPrevious[indiceI][j] != 0) ){
                 if( (findInArray(newCommunityItens, j, index) == -1)  && ((*complexNetwork).neighborhoodMatrix[indiceJ][j] != 0) ){
                     newCommunityItens[index] = j;
                     index++;
                 }
            }else if ((*complexNetwork).neighborhoodMatrix[indiceI][j] != 0) {
                if (findInArray(itensOldCommunity, j, index2) == -1) {
                    itensOldCommunity[index2] = j;
                    index2++;
                }
            }
        } 
    }
    
    if(index2 > index){
        for (j = 0; j < index2; j++){
            AuxCommunity[j] = itensOldCommunity[j];
        }
        
        for (j = 0; j < index; j++){
            itensOldCommunity[j] = newCommunityItens[j];
        }
        
        for (j = 0; j < index2; j++){
            newCommunityItens[j] = AuxCommunity[j];
        } 
        index3 = index2;
        index2 = index;
        index = index3;
    }
    
    
    for (j = 0; j < sizeSimilarityMatrix; j++){
        (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][j] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1][j];             
    }
    (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][sizeSimilarityMatrix] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1][sizeSimilarityMatrix] ;             

    saveGraphicDendrogramInFile((*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1]);
    
	i = (*complexNetwork).numberDendrogramLines;
        
    if(hasNewCommunity == 1){
        float positionBefore = 0.0;
        (*complexNetwork).numberOfCommunities = (*complexNetwork).numberOfCommunities + 1;
        
        for (j = 0; j < sizeSimilarityMatrix; j++){
            finded1 = findInArray(newCommunityItens, j, index);
            finded2 = findInArray(itensOldCommunity, j, index2);
            realIndex = findRealIndex((*complexNetwork).realIndexDendrogramGraphic, j);
            positionBefore = (*complexNetwork).dendrogramGraphic[i-1][realIndex];
            
            if(finded2 != -1){    
                (*complexNetwork).dendrogramGraphic[i][realIndex] = positionBefore - (index/dividendo);
            }else if(finded1 != -1){
                (*complexNetwork).dendrogramGraphic[i][realIndex] = positionBefore + (index2/dividendo);
            }else{
                (*complexNetwork).dendrogramGraphic[i][realIndex] = (*complexNetwork).dendrogramGraphic[i-1][realIndex];             
            }                    
        }
   
        sortGraphicDendrograma(complexNetwork);
        
        (*complexNetwork).dendrogramGraphic[i][sizeSimilarityMatrix] = (*complexNetwork).dendrogramGraphic[i-1][sizeSimilarityMatrix] + 1; 

	    saveGraphicDendrogramInFile((*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1]);
		for (j = 0; j < sizeSimilarityMatrix+1; j++){
			(*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines-1][j] = (*complexNetwork).dendrogramGraphic[(*complexNetwork).numberDendrogramLines][j];             
		}
    }  
            
    free(newCommunityItens);
    free(AuxCommunity);
    free(itensOldCommunity);   
}
//-----------------------------------------------------------------------------

void sortGraphicDendrograma(complexNetwork *complexNetwork){
    int i = (*complexNetwork).numberDendrogramLines;
    int j, l, m, auxIndice, realIndexL, realIndexM = 0;
    float aux = 0.0;
    
    for (j = (sizeSimilarityMatrix-1); j > -1; j--){
        m = j; 
        for (l = (j+1); l < sizeSimilarityMatrix; l++){ 
             if((*complexNetwork).dendrogramGraphic[i][l] < (*complexNetwork).dendrogramGraphic[i][m]){
                 aux = (*complexNetwork).dendrogramGraphic[i][m];
                 (*complexNetwork).dendrogramGraphic[i][m] = (*complexNetwork).dendrogramGraphic[i][l];
                 (*complexNetwork).dendrogramGraphic[i][l] = aux;
                 auxIndice = (*complexNetwork).realIndexDendrogramGraphic[m];
                 (*complexNetwork).realIndexDendrogramGraphic[m] =  (*complexNetwork).realIndexDendrogramGraphic[l];
                 (*complexNetwork).realIndexDendrogramGraphic[l] = auxIndice;
                 realIndexL = (*complexNetwork).realIndexDendrogramGraphic[l];                          
                 realIndexM = (*complexNetwork).realIndexDendrogramGraphic[m];
                 m++;
             }
         }
    }
}
//-----------------------------------------------------------------------------

void createColorsMatrix(complexNetwork *complexNetwork){
    int index = 0, index1, index2 = 0, i;
    int value = 0;
    int * arrayToSort = createDynamicArray(sizeSimilarityMatrix);
  
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        arrayToSort[i] = i;
    } 

    for (index = 0; index < sizeSimilarityMatrix; index++){
        index1 =(*complexNetwork).realIndex[index];
        if(index1 != index){       
            
            index1 = findRealIndex(arrayToSort, index1);
            value = arrayToSort[index1];
            arrayToSort[index1] = arrayToSort[index];
            arrayToSort[index] = value;
            
            for (index2 = 0; index2 < sizeSimilarityMatrix; index2++){              
                value = (*complexNetwork).colorsMatrix[index2][index1];
                (*complexNetwork).colorsMatrix[index2][index1] = (*complexNetwork).colorsMatrix[index2][index];
                (*complexNetwork).colorsMatrix[index2][index] = value;
            }

            for (index2 = 0; index2 < sizeSimilarityMatrix; index2++){              
                value = (*complexNetwork).colorsMatrix[index1][index2];
                (*complexNetwork).colorsMatrix[index1][index2] = (*complexNetwork).colorsMatrix[index][index2];
                (*complexNetwork).colorsMatrix[index][index2] = value;
            }          
        }     
    }  
}
//-----------------------------------------------------------------------------

void orderSequences(complexNetwork *complexNetwork){
    int index = 0, index1, i, valueInt;
    char* value = 0;
    
    int * arrayToSort = createDynamicArray(sizeSimilarityMatrix);
  
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        arrayToSort[i] = i;
    } 
     
    for (index = 0; index < sizeSimilarityMatrix; index++){
        index1 = (*complexNetwork).realIndexDendrogramGraphic[index];
        if(index1 != index){     
            index1 = findRealIndex(arrayToSort, index1);
            
            valueInt = arrayToSort[index1];
            arrayToSort[index1] = arrayToSort[index];
            arrayToSort[index] = valueInt;
            
            value = sequenceArray[index];
            sequenceArray[index] = sequenceArray[index1];
            sequenceArray[index1] = value;  
        }     
    }  
}
//-----------------------------------------------------------------------------

int findInArray(int* array, int item, int sizeOfArray){    
     int i;
     for (i = 0; i < sizeOfArray; i++){ 
         if(array[i] == item){
             return item;
         }
     }
     return -1;
}
//-----------------------------------------------------------------------------

void showRunningClock(void)
{
	// running bar:
	static int cnt = 0;
	const char aSign[] = {'|', '/', '-', '\\'};
	static int idx = 0;
	if (cnt % 20 == 0) {
		idx = (idx < sizeof(aSign)-1)? idx + 1 : 0;
		printf("\r %c", aSign[idx]); fflush(NULL);
	}
	cnt = (cnt < 10)? cnt+1 : 0;
}
//-----------------------------------------------------------------------------

void showRemovedEdges(unsigned long totalRemoved, time_t start)
{
	char sBuf[64] = {0};
	sprintf(sBuf, "%s %d  ", "\r   Edges removed:", totalRemoved);
	showElapsed(sBuf, start);
}
//-----------------------------------------------------------------------------

bool isThereEdge(int** adjacencyMatrix){
    int i,j;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = 0; j < sizeSimilarityMatrix; j++){
            if(adjacencyMatrix[i][j] == 1){
                return true;
            }
        }
     }

    return false;
}
//-----------------------------------------------------------------------------

void removeEdgesWithBiggerBetweenness(index_t nodeWithEdgeRemoved[1], float** betweennessMatrix, int** adjacencyMatrix){
    int i,j;
    float biggerValue = 0;
    index_t *already;		//already[sizeSimilarityMatrix*200][2];
    int sizeOfAlready = 0;
    int indexToRemove;

    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i+1; j < sizeSimilarityMatrix; j++){
            if(betweennessMatrix[i][j] > biggerValue){
               biggerValue = betweennessMatrix[i][j];
            }
        }
    }
	
	// count the number of edges to remove:
	sizeOfAlready = 0;
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i+1; j < sizeSimilarityMatrix; j++){
            if(betweennessMatrix[i][j] == biggerValue){  
                 sizeOfAlready++;
            }
        }
    }  
	// allocate the necessary memory to the 'already' matrix:
    already = (index_t *)malloc( (sizeSimilarityMatrix * sizeOfAlready) * sizeof(index_t));
    if(already == NULL){
         printf("\n ERROR IN CREATING THE REMOVED EDGES INDEX MATRIX. PROBLEM IN MEMORY ALLOCATION. \n");
         exit(1);
    }
	// fill the 'already' matrix with the indexes of the edges to be removed:
	sizeOfAlready = 0; 
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i+1; j < sizeSimilarityMatrix; j++){
            if(betweennessMatrix[i][j] == biggerValue){  
                 already[sizeOfAlready].row = i;
                 already[sizeOfAlready].col = j;
                 sizeOfAlready++;
            }
        }
    }  

    srand(time(NULL));
    // srand(100);  For TESTS only! VERIRY....
    indexToRemove = rand() % sizeOfAlready;
    i = already[indexToRemove].row;
    j = already[indexToRemove].col;
    adjacencyMatrix[i][j] = 0;
    adjacencyMatrix[j][i] = 0;
    nodeWithEdgeRemoved[qtdNodeWithEdgeRemoved].row = i;
    nodeWithEdgeRemoved[qtdNodeWithEdgeRemoved].col = j;
    qtdNodeWithEdgeRemoved++;
}
//-----------------------------------------------------------------------------

void adjustNeighborhoodMatrix(index_t nodeWithEdgeRemoved[1], int** adjacencyMatrix, int** neighborhoodMatrix, int diameter){
    int indiceI, indiceJ, p, i, d, v1,v2, pathSize, item;
    int vertex[sizeSimilarityMatrix];
    int sizeVertex = 0, removeEdge = 1, sizeAux = 0, sizeEdges = 0, hasPath = 0;
    int aux, finded;
    
    indiceI = nodeWithEdgeRemoved[0].row;
    indiceJ = nodeWithEdgeRemoved[0].col;
    neighborhoodMatrix[indiceI][indiceJ] = 0;
    neighborhoodMatrix[indiceJ][indiceI] = 0;
    vertex[sizeVertex] = indiceI;
    sizeVertex++;
    vertex[sizeVertex] = indiceJ;
    sizeVertex++;

    removedEdges = (edge *)malloc( (sizeSimilarityMatrix * 500) * sizeof(edge));
    
    if(removedEdges == NULL){
         printf("\n ERROR IN CREATING THE EDGES LIST. PROBLEM IN MEMORY ALLOCATION. \n");
         exit(1);
    }
    removedEdges[sizeEdges].V1 = indiceI;
    removedEdges[sizeEdges].V2 = indiceJ;
    sizeEdges++;
       
    for (d = 1; d < diameter; d++){
        sizeAux = 0;
        for (i = 0; i < sizeVertex; i++){ 
            item = vertex[i];
            for (v1 = 0; v1 < sizeSimilarityMatrix; v1++){
                if(neighborhoodMatrix[item][v1] == (d+1) ){
                    for (v2 = 0; v2 < sizeSimilarityMatrix; v2++){
                         removeEdge = 1;
                         if( ( (neighborhoodMatrix[item][v2] == 1)&&(neighborhoodMatrix[v1][v2] == d) ) || ( (neighborhoodMatrix[v1][v2] == 1)&&(neighborhoodMatrix[item][v2] == d) ) )  {
                             removeEdge = 0;
                             break;
                         }
                     }
                     if(removeEdge){
                         neighborhoodMatrix[item][v1] = 0;
                         neighborhoodMatrix[v1][item] = 0;
                         finded = 0;
                        
                         for(aux = 0; aux < (sizeVertex-1)+sizeAux; aux++){
                             if(vertex[aux] == v1){
                                 finded = 1;
                                 break;
                             }
                         }

                         if(finded == 0){
                             sizeAux++;
                             vertex[(sizeVertex-1)+sizeAux] = v1;
                         }
                      
                         if(item < v1){
                              removedEdges[sizeEdges].V1 = item;
                              removedEdges[sizeEdges].V2 = v1;
                         }else{
                              removedEdges[sizeEdges].V1 = v1;
                              removedEdges[sizeEdges].V2 = item;
                         }
                         sizeEdges++;   
                     } 
                }
            }
        }
        sizeVertex = sizeVertex + sizeAux;
    }   

    for (d = 2; d <= (2*diameter); d++){
        pathSize = 0;
        for (i = 0; i < sizeEdges; i++){ 
            v1 =  removedEdges[i].V1;
            v2 =  removedEdges[i].V2;            
            if(neighborhoodMatrix[v1][v2] == 0 ){
                pathSize = 1;
                hasPath = 0;
                while( (hasPath==0) && (pathSize < d) ){
                    for (p = 0; p < sizeSimilarityMatrix; p++){                        
                         if( ( (neighborhoodMatrix[v1][p] == pathSize)&&(neighborhoodMatrix[v2][p] == (d-pathSize)) ) || ( (neighborhoodMatrix[v1][p] == (d-pathSize))&&(neighborhoodMatrix[v2][p] == pathSize) ) )  {
                             hasPath = 1;
                             break;
                         }                        
                    }
                 
                    if(hasPath){
                        neighborhoodMatrix[v1][v2] = d;
                        neighborhoodMatrix[v2][v1] = d;
                    }
                    
                    pathSize++;
                }             
            }
        }

        if(pathSize == 0){
            break;
        }
    }           
    free(removedEdges);
}
//-----------------------------------------------------------------------------

void calculateBetweenness(int diameter, complexNetwork complexNetwork, float** betweennessMatrix){
    int i, j, l, k, m;
    copyMatrixOfFloat(complexNetwork.adjacencyMatrix, betweennessMatrix);
    int** D = complexNetwork.neighborhoodMatrix;
    
    int requirementOK[sizeSimilarityMatrix];
    memset(requirementOK, 0x7f, sizeof(requirementOK));
    int reqCount = 0;
    
    for (l = diameter; l > 1; l--){
        for (i = 0; i < sizeSimilarityMatrix; i++){
              for (j = (i+1); j < sizeSimilarityMatrix; j++){
                  
                  if( D[i][j] == l ){
                      reqCount = 0;
                      for (k = 0; k < sizeSimilarityMatrix; k++){
                           if( ( (D[i][k] == 1) && (D[j][k] == (l-1))  ) || ( (D[i][k] == (l-1)) && (D[j][k] == 1) ) ){
                               requirementOK[reqCount] = k;
                               reqCount++;
                           }
                       }
                       
                       for (m = 0; m < reqCount; m++){
                           k = requirementOK[m];
                           betweennessMatrix[i][k] = betweennessMatrix[i][k] + ((betweennessMatrix[i][j] + 1)/reqCount);
                           betweennessMatrix[j][k] = betweennessMatrix[j][k] + ((betweennessMatrix[i][j] + 1)/reqCount);
                           betweennessMatrix[k][i] = betweennessMatrix[i][k];
                           betweennessMatrix[k][j] = betweennessMatrix[j][k]; 
                       }
                  }
              }
         }
    }

    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i; j < sizeSimilarityMatrix; j++){
            if (complexNetwork.adjacencyMatrix[i][j] == 0){
                 betweennessMatrix[i][j] = 0;
                 betweennessMatrix[j][i] = 0;
            }         
        }
    }
}
//-----------------------------------------------------------------------------

void copyMatrixOfFloat(int** sourceMatrix, float** targetMatrix){
    int i,j;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i; j < sizeSimilarityMatrix; j++){
            targetMatrix[i][j] = (float)sourceMatrix[i][j];
            targetMatrix[j][i] = (float)sourceMatrix[j][i];
        }
    }   
}
//-----------------------------------------------------------------------------

void copyMatrix(int** sourceMatrix, int** targetMatrix){
    int i,j;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i+1; j < sizeSimilarityMatrix; j++){
            targetMatrix[i][j] = sourceMatrix[i][j];
            targetMatrix[j][i] = sourceMatrix[j][i];
        }
    }    
}
//-----------------------------------------------------------------------------

void clearMatrix(int** matrix){
    int i,j;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        for (j = i; j < sizeSimilarityMatrix; j++){
            matrix[i][j] = 0;
            matrix[j][i] = 0;
        }
    }    
}
//-----------------------------------------------------------------------------

void clearMatrixBySize(int** matrix, int size){
    int i,j;
    
    for (i = 0; i < size; i++){ 
        for (j = i; j < size; j++){
            matrix[i][j] = 0;
            matrix[j][i] = 0;
        }
    }
}
//-----------------------------------------------------------------------------

void clearVetor(int* vetor){
    int i;
    
    for (i = 0; i < sizeSimilarityMatrix; i++){ 
        vetor[i] = 0;
    }    
}
//-----------------------------------------------------------------------------
//-----------------------------------------------------------------------------