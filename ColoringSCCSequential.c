#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"


int main(int argc, char* argv[]){
    FILE *f;
    MM_typecode matcode;
    int  ret_code;
    int cols, rows, nz;
    int i;

    if (argc < 2){
        printf("Where is the file, bro?\n");
        exit(1);
    }
    else if ((f = fopen(argv[1], "r")) == NULL){
        printf("It's like you are not even trying, bro\n");
        exit(1);
    }

    if(mm_read_banner(f, &matcode) != 0){
        printf("Now, now, it's not your fault, bro\n");
        exit(1);
    }

    if ((ret_code = mm_read_mtx_crd_size(f, &cols, &rows, &nz)) !=0)
        exit(1);
    
    //printf("Lines: %d\nColumns: %d\nNon-zero: %d\n", M, N, nz);

    //CSR data structure
    int value;
    int* CSR_ROW_INDEX = (int *) malloc((rows + 1) * sizeof(int));
    int* CSR_COL_INDEX = (int *) malloc(nz * sizeof(int));
    int readRow;

    for (i = 0; i < rows + 1; i++)
        CSR_ROW_INDEX[i] = 0;

    for (i = 0; i < nz; i++) {
        fscanf(f, "%d %d %d\n", &value, &readRow, &CSR_COL_INDEX[i]);
        CSR_COL_INDEX[i]--;
        CSR_ROW_INDEX[readRow]++;
    }

    for (i = 0; i < rows; i++)
        CSR_ROW_INDEX[i + 1] += CSR_ROW_INDEX[i];

    //reset file pointer at start of data
    fseek(f, 0, SEEK_SET);
    ret_code = mm_read_mtx_crd_size(f, &cols, &rows, &nz);

    //CSC data structure
    int* CSC_ROW_INDEX = (int *) malloc(nz * sizeof(int));
    int* CSC_COL_INDEX = (int *) malloc((cols + 1) * sizeof(int));
    int readCol;

    for (i = 0; i < nz; i++)
        CSC_COL_INDEX[i] = 0;
    
    for (i = 0; i < nz; i++){
        fscanf(f, "%d %d %d\n", &value, &CSC_ROW_INDEX[i], &readCol);
        CSC_ROW_INDEX[i]--;
        CSC_COL_INDEX[readCol]++;
    }

    for (i = 0; i < cols; i++)
        CSC_COL_INDEX[i + 1] += CSC_COL_INDEX[i];


    free(CSR_ROW_INDEX);
    free(CSR_COL_INDEX);

    return 0;
}