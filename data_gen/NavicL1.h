char ***createStrArr(int m, int n, int stringLength);
char ***loadStrArr(char *filename, int m, int n, int stringLength);
void writeStrArr(char *filename, char ***matrix, int m, int n, int stringLength);
void printStrArr(char ***matrix, int m, int n);
void octalToBinaryDigit(char octalDigit, char* binaryResult) ;
char ***OctToBin(char ***arr,int m,int n,int l);
int* shift(int *register_array, int feedback, int length) ;
int* genNavicCaCode_Data(int sv ,char ***SV_L1_Data_r0,char ***SV_L1_Data_r1,char ***SV_L1_Data_C);
int* genNavicCaTable_Data(float samplingFreq, int codeLength, float codeFreqBasis, int* satId, int satIdLength, char ***r0,char ***r1,char ***c);
int *genNavicCaCode_Pilot(int sv,char ***SV_L1_Pilot_r0,char ***SV_L1_Pilot_r1,char ***SV_L1_Pilot_C);
int* genNavicCaTable_Pilot(float samplingFreq, int codeLength, float codeFreqBasis, int* satId, int satIdLength, char ***r0,char ***r1 ,char ***c);
int* genNavicCaCode_Pilot_Overlay(int sv,char ***SV_L1_Pilot_Overlay_r0,char ***SV_L1_Pilot_Overlay_r1);
int* genNavicCaTable_Pilot_Overlay(float samplingFreq, int codeLength, float codeFreqBasis, int* satId, int satIdLength,char ***r0,char ***r1);
char ***concateZero(char ***buff,int r,int c,int n);
int *packbits(char ***nbuff,int r,int l,int n);
//dddouble* genNavicCaTable_Data(samplingFreq, codeLength, codeFreqBasis, satId );


// function to convert octal to binary
char ***OctToBin(char ***arr,int m,int n,int l)
{
	l=strlen(arr[0][0]);
	char ***binaryResult=createStrArr(m,n,3*l);
	char ***new_binary=createStrArr(m,n,3*l);
    //binaryResult[i][j][0] = '\0';  // Ensure the string is empty initially
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
    	    binaryResult[i][j][0] = '\0';  // Ensure the string is empty initially
	    //binaryResult[i][j][0] = (int)arr[i][j][0];
        //	strcat(binaryResult, arr[i][j][0]);
	    //binaryResult[i][j][3*l+1] = '\0';
            for(int k=1;k<l;k++)
	    { //binaryResult[i][j][0] = (int)arr[i][j][0];
		octalToBinaryDigit((int)arr[i][j][k],binaryResult[i][j]);
	    }
    new_binary[i][j][0]=arr[i][j][0];
    new_binary[i][j][55]='\0';
    strcat(new_binary[i][j],binaryResult[i][j]);
        }
    }

/*    for(int i=0;i<m;i++)
    {
	for(int j=0;j<n;j++)
	{
		for(int k=0;k<l;k++)
		{
		char ***temp=createStrArr(m,n,)
		}
	}
    }
*/
    return new_binary;
    //return binaryResult;
}



void octalToBinaryDigit(char octalDigit, char* binaryResult) {
    switch (octalDigit) {
        case '0': strcat(binaryResult, "000"); break;
        case '1': strcat(binaryResult, "001"); break;
        case '2': strcat(binaryResult, "010"); break;
        case '3': strcat(binaryResult, "011"); break;
        case '4': strcat(binaryResult, "100"); break;
        case '5': strcat(binaryResult, "101"); break;
        case '6': strcat(binaryResult, "110"); break;
        case '7': strcat(binaryResult, "111"); break;
        case '\0': break;
        //default: continue;
        default: printf(" %d ",octalDigit);
        //default: printf("Invalid octal digit");
    }
}




// Function to create a array of strings
char ***createStrArr(int m, int n, int stringLength)
{
    int i, j;
    char ***matrix;

    // Allocate memory for the matrix
    matrix = (char ***)malloc(m * sizeof(char **));
    for (i = 0; i < m; i++)
    {
        matrix[i] = (char **)malloc(n * sizeof(char *));
        for (j = 0; j < n; j++)
        {
            matrix[i][j] = (char *)malloc((stringLength + 1) * sizeof(char)); // +1 for null terminator
        }
    }

    return matrix;
}

// Function to read a array of strings from a file
char ***loadStrArr(char *filename, int m, int n, int stringLength)
{
    FILE *fp;
    char ***matrix;
    int i, j;

    matrix = createStrArr(m, n, stringLength);
    fp = fopen(filename, "r");

    if (fp == NULL)
    {
        perror("Error opening file");
        exit(1);
    }

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            fscanf(fp, "%s", matrix[i][j]);
        }
    }

    fclose(fp);

    return matrix;
}

// Function to write a 3D array of strings to a file
void writeStrArr(char *filename, char ***matrix, int m, int n, int stringLength)
{
    FILE *fp;
    int i, j;

    fp = fopen(filename, "w");

    if (fp == NULL)
    {
        perror("Error opening file");
        exit(1);
    }

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            fprintf(fp, "%s ", matrix[i][j]);
        }
        fprintf(fp, "\n"); // Move to the next line after each row
    }

    fclose(fp);
}


// Function to print a array of strings
void printStrArr(char ***matrix, int m, int n)
{
    int i, j;

    for (i = 0; i < m; i++)
    {
        for (j = 0; j < n; j++)
        {
            printf("%s ", matrix[i][j]);
           // printf("%d ", j);
           // printf("\n ");
        }
        printf("\n");
    }
    printf("\n");
}


int* shift(int *register_array, int feedback, int length) {
    for (int i = 0; i < length-1; i++) {
        register_array[i] = register_array[i+1];
    }
    register_array[length-1] = feedback;
}


int* genNavicCaCode_Data(int sv ,char ***SV_L1_Data_r0,char ***SV_L1_Data_r1,char ***SV_L1_Data_C) {
    if (sv < 0 || sv > 63) {
        printf("Error: PRN ID out of bounds!\n");
        return NULL;
    } 
    
        int r0[55],r1[55],C[5];
	for(int i=0;i<55;i++)
	{
		r0[i]=SV_L1_Data_r0[0][sv][i]-'0';
		r1[i]=SV_L1_Data_r1[0][sv][i]-'0';
		if(i<5)
		{
			C[i]=SV_L1_Data_C[0][sv][i]-'0';
		}
        }

   int codeLength = 10230;
         int* cad = malloc( codeLength* sizeof(int));
        for (int j = 0; j < codeLength; j++) {
            int a = ((r0[50] ^ r0[45]) ^ r0[40]) ^ ((r0[20] ^ r0[10]) ^ (r0[5] ^ r0[0]));
     
            int sigma2A = ((r0[50] ^ r0[45]) ^ r0[40]) & ((r0[20] ^ r0[10]) ^ (r0[5] ^ r0[0]));
            int sigma2B = ((r0[50] ^ r0[45]) & r0[40]) ^ ((r0[20] ^ r0[10]) & (r0[5] ^ r0[0]));
            int sigma2C = (r0[50] & r0[45]) ^ ((r0[20] & r0[10]) ^ (r0[5] & r0[0]));
            int sigma2 = (sigma2A ^ sigma2B) ^ sigma2C;
            int temp = ((r0[40] ^ r0[35]) ^ (r0[30] ^ r0[25])) ^ (r0[15] ^ r0[0]);
            int R1A = sigma2 ^ temp;
            int R1B = ((r1[50] ^ r1[45]) ^ (r1[40] ^ r1[20])) ^ ((r1[10] ^ r1[5]) ^ (r1[0]));
            int b = R1A ^ R1B;
            int z = r1[0] ^ C[0];
            shift(C, C[0], 5);
            shift(r1, b, 55);
            shift(r0, a, 55);
            cad[j] = z;
        }
//	printf("cad=%d\n",*cad);
        return cad ;
         }







//function to upsample the Data PRN sequence generated to required sampling rate
int* genNavicCaTable_Data(float samplingFreq, int codeLength, float codeFreqBasis, int* satId, int satIdLength, char ***r0,char ***r1,char ***c) {
    float samplingPeriod = 1.0 / samplingFreq;  
    int sampleCount = (int)(samplingFreq / (codeFreqBasis / codeLength));
    float* indexArr = (float*)malloc(sampleCount * sizeof(float));
    int* result = (int*)malloc(sampleCount * satIdLength * sizeof(int));

    if (indexArr == NULL || result == NULL) {
        // Handle memory allocation failure
        free(indexArr);
        free(result);
        return NULL;
    }

    // Precompute indexArr values outside the loop
    for (int i = 0; i < sampleCount; i++) {
        indexArr[i] = i * samplingPeriod * codeFreqBasis;
    }

    FILE* file_data = fopen("code_table_data.txt", "w");
    if (file_data == NULL) {
        printf("Failed to open the file for writing.\n");
        free(indexArr);
        free(result);
        return NULL;
    }

    for (int i = 0; i < satIdLength; i++) {
        int* satPrn = genNavicCaCode_Data(satId[i],r0,r1,c);
        if (satPrn == NULL) {
            // Handle error in genNavicCaCode_Data function
            fclose(file_data);
            free(indexArr);
            free(result);
            return NULL;
        }

        for (int j = 0; j < sampleCount; j++) {
            result[i * sampleCount + j] = satPrn[(int)indexArr[j]];
            fprintf(file_data, "%d\n", result[i * sampleCount + j]);
        }

        free(satPrn); // Free the memory allocated in genNavicCaCode_Data function
    }

    fclose(file_data);
    free(indexArr);
    printf("d_result=%d\n",*result);
    
    return result;

}



int *genNavicCaCode_Pilot(int sv,char ***SV_L1_Pilot_r0,char ***SV_L1_Pilot_r1,char ***SV_L1_Pilot_C) {
    if (sv < 0 || sv > 63) {
        printf("Error: PRN ID out of bounds!\n");
        return NULL;
    }
       int r0_p[55],r1_p[55],C_p[5];
	 for(int i=0;i<55;i++)
	{
		r0_p[i]=SV_L1_Pilot_r0[0][sv][i]-'0';
		r1_p[i]=SV_L1_Pilot_r1[0][sv][i]-'0';
		if(i<5)
		{
			C_p[i]=SV_L1_Pilot_C[0][sv][i]-'0';
		}
        }

        int codeLength = 10230;
    	  int* cap = malloc(codeLength * sizeof(int));

    // Initialize r0_p, r1_p, and C_p based on SV_L1_Pilot_r0, SV_L1_Pilot_r1, and SV_L1_Pilot_C arrays

    for (int k = 0; k < codeLength; k++) {
        int r_p = (r0_p[50] ^ r0_p[45]) ^ (r0_p[40]) ^ ((r0_p[20] ^ r0_p[10]) ^ (r0_p[5] ^ r0_p[0]));
        int sigma2A = (r0_p[50] ^ r0_p[45] ^ r0_p[40]) & (r0_p[20] ^ r0_p[10] ^ r0_p[5] ^ r0_p[0]);
        int sigma2B = ((r0_p[50] ^ r0_p[45]) & (r0_p[40])) ^ ((r0_p[20] ^ r0_p[10]) & (r0_p[5] ^ r0_p[0]));
        int sigma2C = (r0_p[50] & r0_p[45]) ^ (r0_p[20] & r0_p[10]) ^ (r0_p[5] & r0_p[0]);
        int sigma2 = sigma2A ^ sigma2B ^ sigma2C;
        int temp = r0_p[40] ^ r0_p[35] ^ r0_p[30] ^ r0_p[25] ^ r0_p[15] ^ r0_p[0];
        int R1A = sigma2 ^ temp;
        int R1B = ((r1_p[50] ^ r1_p[45]) ^ (r1_p[40] ^ r1_p[20])) ^ ((r1_p[10] ^ r1_p[5]) ^ (r1_p[0]));
        int r3_p = R1A ^ R1B;
        int z1 = (r1_p[0] + C_p[0]) % 2;

        shift(r0_p, r_p, 55);
        shift(r1_p, r3_p, 55);
        shift(C_p, C_p[0], 5);

        cap[k] = z1;
    }
   // printf("cap=%d\n",*cap);
    return cap;
}



int* genNavicCaTable_Pilot(float samplingFreq, int codeLength, float codeFreqBasis, int* satId, int satIdLength, char ***r0,char ***r1 ,char ***c) {
    float samplingPeriod = 1.0 / samplingFreq;  
    int sampleCount = (int)(samplingFreq / (codeFreqBasis / codeLength));
    float* indexArr = (float*)malloc(sampleCount * sizeof(float));
    int* result = (int*)malloc(sampleCount * satIdLength * sizeof(int));

    if (indexArr == NULL || result == NULL) {
        // Handle memory allocation failure
        free(indexArr);
        free(result);
        return NULL;
    }

    // Precompute indexArr values outside the loop
    for (int i = 0; i < sampleCount; i++) {
        indexArr[i] = i * samplingPeriod * codeFreqBasis;
    }

    FILE* file_data = fopen("code_table_pilot.txt", "w");
    if (file_data == NULL) {
        printf("Failed to open the file for writing.\n");
        free(indexArr);
        free(result);
        return NULL;
    }

    for (int i = 0; i < satIdLength; i++) {
        int* satPrn = genNavicCaCode_Pilot(satId[i],r0,r1,c);
        if (satPrn == NULL) {
            // Handle error in genNavicCaCode_Pilot function
            fclose(file_data);
            free(indexArr);
            free(result);
            return NULL;
        }

        for (int j = 0; j < sampleCount; j++) {
            result[i * sampleCount + j] = satPrn[(int)indexArr[j]];
            fprintf(file_data, "%d\n", result[i * sampleCount + j]);
        }

        free(satPrn); // Free the memory allocated in genNavicCaCode_Pilot function
    }

    fclose(file_data);
    free(indexArr);
    printf("p_result%d\n",*result);
    return result;
}




int* genNavicCaCode_Pilot_Overlay(int sv,char ***SV_L1_Pilot_Overlay_r0,char ***SV_L1_Pilot_Overlay_r1) {
    // Check if sv is within the valid range
    if (sv < 0 || sv > 63) {
        printf("Error: PRN ID out of bounds!\n");
        return NULL;
    }

    // Initialize registers
    int r0_pl[10], r1_pl[10];
    for (int i = 0; i < 10; i++) {
        r0_pl[i] = SV_L1_Pilot_Overlay_r0[0][sv][i] - '0';
        r1_pl[i] = SV_L1_Pilot_Overlay_r1[0][sv][i] - '0';
    }

    int* ca = (int*)malloc(1800 * sizeof(int)); // Allocate memory for the PRN sequence

    // Generate Pilot Overlay PRN sequence
    for (int l = 0; l < 1800; l++) {
        int r_pl = (r0_pl[5] ^ r0_pl[2]) ^ (r0_pl[1] ^ r0_pl[0]);
        // Compute σ2A
        int sigma2A = (r0_pl[5] ^ r0_pl[2]) & (r0_pl[1] ^ r0_pl[0]);
        // Compute σ2B
        int sigma2B = (r0_pl[5] & r0_pl[2]) ^ (r0_pl[1] & r0_pl[0]);
        // Compute σ2
        int sigma2 = sigma2A ^ sigma2B;
        // Compute r1A
        int temp = r0_pl[6] ^ r0_pl[3] ^ r0_pl[2] ^ r0_pl[0];
        int R1A = sigma2 ^ temp;
        // Compute r1B
        int R1B = r1_pl[5] ^ r1_pl[2] ^ r1_pl[1] ^ r1_pl[0];
        int r3_pl = R1A ^ R1B;
        int z2 = r1_pl[0];

        shift(r0_pl, r_pl, 10);
        shift(r1_pl, r3_pl, 10);

        // Modulo 2 add and store in the PRN sequence array
        ca[l] = z2;
    }
//printf("ca%d\n",*ca);
    return ca;
}



int* genNavicCaTable_Pilot_Overlay(float samplingFreq, int codeLength, float codeFreqBasis, int* satId, int satIdLength,char ***r0,char ***r1) {
    float samplingPeriod = 1.0 / samplingFreq; 
    int sampleCount = (int)(samplingFreq / (codeFreqBasis / codeLength));

    float* indexArr = (float*)malloc(sampleCount * sizeof(float));
    int* result = (int*)malloc(sampleCount * satIdLength * sizeof(int));

    if (indexArr == NULL || result == NULL) {
        // Handle memory allocation failure
        free(indexArr);
        free(result);
        return NULL;
    }

    // Precompute indexArr values outside the loop
    for (int i = 0; i < sampleCount; i++) {
        indexArr[i] = i * samplingPeriod * codeFreqBasis;
    }
  FILE* file_data = fopen("code_table_pilot_overlay.txt", "w");
    if (file_data == NULL) {
        printf("Failed to open the file for writing.\n");
        free(indexArr);
        free(result);
        return NULL;
    }

    for (int i = 0; i < satIdLength; i++) {
        int* satPrnOverlay = genNavicCaCode_Pilot_Overlay(satId[i],r0,r1);
        if (satPrnOverlay == NULL) {
            // Handle error in genNavicCaCode_Pilot_Overlay function
            fclose(file_data);
            free(indexArr);
            free(result);
            return NULL;
        }

        for (int j = 0; j < sampleCount; j++) {
            result[i * sampleCount + j] = satPrnOverlay[(int)indexArr[j]];
            fprintf(file_data, "%d\n", result[i * sampleCount + j]);
        }

        free(satPrnOverlay); // Free the memory allocated in genNavicCaCode_Pilot_Overlay function
    }

    fclose(file_data);
    free(indexArr);
printf("ol_result%d\n",*result);
    return result;
}



///function for padding n zeros before buff  
char ***concateZero(char ***buff,int r,int c,int n)
{

char ***newbuff=createStrArr(r, c+n, 1);


    for (int i = 0; i < r; i++)
    {
        for (int j = 0; j < n+c; j++)
        {
    	    newbuff[i][j][0]= '\0';  // Ensure the string is empty initially
	if(j<n)	
	    newbuff[i][j][0] = '0';
	else	
	    newbuff[i][j][0] = buff[i][j-n][0];

	}
    }	
/*
    for (int i = 0; i < r; i++)
    {
        for (int j = n; j < c+n; j++)
        {
    	    newbuff[i][j][0]= '\0';  // Ensure the string is empty initially
		newbuff[i][j][0] = buff[i][j-n][0];
		printf(",");

	}
    }
*/
return newbuff;
}




/*
int *packbits(char ***bufft,int r,int l,int n)
{

int k= l/n;
//int *packed_data = (int  *)malloc( k * sizeof(int));

for (int i=0;i<k;i++)
{	int dec=0;
	int t=(i+1)*n-1,ind=0;
for(int j=t;j>=t-n-1;j--)
{
dec=dec+ (*bufft[0][j]-'0')*pow(2,7-ind);
ind++;
}
//packed_data[i]=dec;
printf("%d\n",dec);
}


//return packed_data;
}
*/
