#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include "mmio.h"

int p = 4;

typedef struct
{
    int *vertices;
    int *rowsOutgoingEdges;
    int *colsOutgoingEdges;
    int *colors;
    int u;
    int *colorChanges;
} ColorArgs;

typedef struct
{
    int *colors;
    int *rowsIncomingEdges;
    int *colsIncomingEdges;
    int u;
    int numVertices;
    int *vertices;
    int *SSCIDs;
} BFSArgs;

void *BFS(void *arguments)
{
    BFSArgs *args = (BFSArgs *)arguments;
    int i = -1, j, size = 0;
    int *queue = (int *)calloc(args->numVertices, sizeof(int));
    int *visited = (int *)calloc(args->numVertices, sizeof(int));
    int v;
    queue[++i] = args->u;
    visited[args->u] = 1;

    while (i != -1)
    {
        v = queue[i--];
        args->vertices[v] = 0;
        args->SSCIDs[v] = args->u + 1;
        size++;
        for (j = args->colsIncomingEdges[v]; j < args->colsIncomingEdges[v + 1]; j++)
        {
            if (args->colors[args->rowsIncomingEdges[j]] == args->u + 1 && visited[args->rowsIncomingEdges[j]] != 1)
            {
                queue[++i] = args->rowsIncomingEdges[j];
                visited[args->rowsIncomingEdges[j]] = 1;
            }
        }
    }

    free(queue);
    free(visited);

    int *result = (int *)malloc(sizeof(int));
    *result = size;

    return (void *)result;
}

void *propagateColors(void *arguments)
{
    ColorArgs *args = (ColorArgs *)arguments;
    int i;
    for (i = args->rowsOutgoingEdges[args->u]; i < args->rowsOutgoingEdges[args->u + 1]; i++)
    {
        if (args->colors[args->u] > args->colors[args->colsOutgoingEdges[i]] && args->vertices[args->colsOutgoingEdges[i]] != 0)
        {
            args->colors[args->colsOutgoingEdges[i]] = args->colors[args->u];
            args->colorChanges[0] = 1;
        }
    }
}

int *coloringSCC(int *rowsOutgoingEdges, int *colsOutgoingEdges, int *rowsIncomingEdges, int *colsIncomingEdges, int numVertices)
{
    int *colors = (int *)malloc(numVertices * sizeof(int));
    int *vertices = (int *)malloc(numVertices * sizeof(int));
    int *SCCIDs = (int *)calloc(numVertices, sizeof(int));

    int *verticesRemoved;
    pthread_t *threads = (pthread_t *)malloc(p * sizeof(pthread_t));
    int *openThreads = (int *)calloc(p, sizeof(int));
    ColorArgs *colorArgs = (ColorArgs *)malloc(p * sizeof(ColorArgs));
    BFSArgs *bfsArgs = (BFSArgs *)malloc(p * sizeof(BFSArgs));

    int verticesRemaining = numVertices;
    int i, j;
    int colorChange;

    for (i = 0; i < numVertices; i++)
        vertices[i] = i + 1;

    while (verticesRemaining != 0)
    {
        for (i = 0; i < numVertices; i++)
        {
            colors[i] = vertices[i];
        }

        do
        {
            int *colorChanges = (int *)calloc(p, sizeof(int));
            colorChange = 0;
            i = 0;

            while (i < numVertices)
            {
                for (j = 0; j < p; j++)
                {
                    while (vertices[i] == 0 && i < numVertices)
                        i++;
                    if (i == numVertices)
                        break;

                    colorArgs[j].vertices = vertices;
                    colorArgs[j].rowsOutgoingEdges = rowsOutgoingEdges;
                    colorArgs[j].colsOutgoingEdges = colsOutgoingEdges;
                    colorArgs[j].colors = colors;
                    colorArgs[j].u = i;
                    colorArgs[j].colorChanges = &colorChanges[j];
                    pthread_create(&threads[j], NULL, propagateColors, (void *)&colorArgs[j]);
                    openThreads[j] = 1;

                    i++;
                }
                for (j = 0; j < p; j++)
                {
                    if (openThreads[j] == 1)
                    {
                        pthread_join(threads[j], NULL);
                        openThreads[j] = 0;
                    }
                }
            }

            for (i = 0; i < p; i++)
            {
                if (colorChanges[i] == 1)
                {
                    colorChange = 1;
                    break;
                }
            }
            free(colorChanges);
        } while (colorChange);

        while (i < numVertices)
        {
            for (j = 0; j < p; j++)
            {
                while (colors[i] != i + 1 && i < numVertices)
                    i++;
                if (i == numVertices)
                    break;

                bfsArgs[j].colors = colors;
                bfsArgs[j].rowsIncomingEdges = rowsIncomingEdges;
                bfsArgs[j].colsIncomingEdges = colsIncomingEdges;
                bfsArgs[j].u = i;
                bfsArgs[j].numVertices = numVertices;
                bfsArgs[j].vertices = vertices;
                bfsArgs[j].SSCIDs = SCCIDs;
                pthread_create(&threads[j], NULL, BFS, (void *)&bfsArgs[j]);

                i++;
            }

            for (j = 0; j < p; j++)
            {
                if (pthread_join(threads[j], (void **)&verticesRemoved) == 0)
                    verticesRemaining = verticesRemaining - *verticesRemoved;
            }
        }
    }
    return SCCIDs;
}

int trimmingCSC(int **colsPointers, int **rowsIndexes, int size, int nz);

int main(int argc, char *argv[])
{
    FILE *f;
    MM_typecode matcode;
    int ret_code;
    int cols, rows, nz;
    int i, j;

    if (argc < 2)
    {
        printf("Where is the file, bro?\n");
        exit(1);
    }
    else if ((f = fopen(argv[1], "r")) == NULL)
    {
        printf("It's like you are not even trying, bro\n");
        exit(1);
    }

    if (mm_read_banner(f, &matcode) != 0)
    {
        printf("Now, now, it's not your fault, bro\n");
        exit(1);
    }

    if ((ret_code = mm_read_mtx_crd_size(f, &rows, &cols, &nz)) != 0)
        exit(1);

    printf("Rows: %d\nColumns: %d\nNon-zero: %d\n", rows, cols, nz);

    // CSC data structure
    int *CSC_COL_INDEX = (int *)malloc((cols + 1) * sizeof(int));
    int *CSC_ROW_INDEX = (int *)malloc(nz * sizeof(int));
    int readCol;

    for (i = 0; i < cols + 1; i++)
        CSC_COL_INDEX[i] = 0;

    if (matcode[2] == 'I' || matcode[2] == 'R')
    {
        double value;
        for (i = 0; i < nz; i++)
        {
            fscanf(f, "%d %d %lf\n", &CSC_ROW_INDEX[i], &readCol, &value);
            CSC_ROW_INDEX[i]--;
            CSC_COL_INDEX[readCol]++;
        }
    }
    else
    {
        for (i = 0; i < nz; i++)
        {
            fscanf(f, "%d %d\n", &CSC_ROW_INDEX[i], &readCol);
            CSC_ROW_INDEX[i]--;
            CSC_COL_INDEX[readCol]++;
        }
    }

    for (i = 0; i < cols; i++)
        CSC_COL_INDEX[i + 1] += CSC_COL_INDEX[i];

    int trimmed = trimmingCSC(&CSC_COL_INDEX, &CSC_ROW_INDEX, cols, nz);
    cols = cols - trimmed;
    rows = cols;
    nz = CSC_COL_INDEX[cols];

    // CSR data structure
    int *CSR_COL_INDEX = (int *)malloc(nz * sizeof(int));
    int *CSR_ROW_INDEX = (int *)malloc((rows + 1) * sizeof(int));
    int row, dest, temp, last = 0, cumsum = 0;

    for (i = 0; i < rows + 1; i++)
    {
        CSR_ROW_INDEX[i] = 0;
    }

    // CONVERT CSC TO CSR
    for (i = 0; i < nz; i++)
        CSR_ROW_INDEX[CSC_ROW_INDEX[i]]++;

    for (i = 0; i < rows; i++)
    {
        temp = CSR_ROW_INDEX[i];
        CSR_ROW_INDEX[i] = cumsum;
        cumsum += temp;
    }
    CSR_ROW_INDEX[rows] = nz;

    for (i = 0; i < cols; i++)
    {
        for (j = CSC_COL_INDEX[i]; j < CSC_COL_INDEX[i + 1]; j++)
        {
            row = CSC_ROW_INDEX[j];
            dest = CSR_ROW_INDEX[row];
            CSR_COL_INDEX[dest] = i;

            CSR_ROW_INDEX[row]++;
        }
    }

    for (i = 0; i < rows + 1; i++)
    {
        temp = CSR_ROW_INDEX[i];
        CSR_ROW_INDEX[i] = last;
        last = temp;
    }

    fclose(f);

    // find SCCs
    struct timeval start, end;

    gettimeofday(&start, NULL);
    int *SCCIDs = coloringSCC(CSR_ROW_INDEX, CSR_COL_INDEX, CSC_ROW_INDEX, CSC_COL_INDEX, rows);
    gettimeofday(&end, NULL);

    int *showTimes = (int *)calloc(rows + 1, sizeof(int));
    int CC = 0, TRIVIAL = 0;

    for (i = 0; i < rows; i++)
        showTimes[SCCIDs[i]]++;
    for (i = 1; i < rows + 1; i++)
    {
        if (showTimes[i] == 1)
            TRIVIAL++;
        else if (showTimes[i] > 1)
            CC++;
    }

    printf("SCCs: %d, Trivial: %d\n", CC, TRIVIAL + trimmed);

    long seconds = end.tv_sec - start.tv_sec;
    long microseconds = end.tv_usec - start.tv_usec;
    double elapsed = seconds + microseconds * 1e-6;
    printf("Time taken: %f\n", elapsed);

    return 0;
}

int trimmingCSC(int **colsPointers, int **rowsIndexes, int cols, int nz)
{
    int i, j, k, flag;
    int numColsRemoved = 0, newSizeCols = 0, newSizeRows = 0;
    int *tempRows = (int *)malloc(nz * sizeof(int));
    int *colsRemoved = (int *)calloc(cols, sizeof(int));
    int *numElements = (int *)calloc(cols, sizeof(int));

    // note columns that are going to be removed
    for (i = 0; i < cols; i++)
    {
        if (colsPointers[0][i] == colsPointers[0][i + 1])
        {
            colsRemoved[numColsRemoved++] = i;
            numElements[i] = -1;
        }
    }

    int *tempCols = (int *)malloc((cols + 1 - numColsRemoved) * sizeof(int));

    // note num and which elements each new column will have
    for (i = 0; i < cols; i++)
    {
        for (j = colsPointers[0][i]; j < colsPointers[0][i + 1]; j++)
        {
            flag = 0;
            for (k = 0; k < numColsRemoved; k++)
            {
                if (rowsIndexes[0][j] == colsRemoved[k])
                    flag = 1;
            }
            if (flag == 0)
            {
                tempRows[newSizeRows++] = rowsIndexes[0][j];
                numElements[i]++;
            }
        }
    }

    // make the new CSC structure on temps
    // temp CSC_COL_INDEX
    for (i = 0; i < cols; i++)
    {
        if (numElements[i] != -1)
            tempCols[++newSizeCols] = numElements[i];
    }

    for (i = 0; i < newSizeCols; i++)
        tempCols[i + 1] += tempCols[i];
    tempCols[0] = 0;
    newSizeCols++;

    // fix indexes of temp CSC_ROW_INDEX
    int indexSubstract;
    for (i = 0; i < newSizeRows; i++)
    {
        indexSubstract = 0;
        for (j = 0; j < numColsRemoved; j++)
        {
            if (tempRows[i] > colsRemoved[j])
                indexSubstract++;
        }
        tempRows[i] -= indexSubstract;
    }

    free(colsPointers[0]);
    colsPointers[0] = (int *)malloc(newSizeCols * sizeof(int));
    free(rowsIndexes[0]);
    rowsIndexes[0] = (int *)malloc(newSizeRows * sizeof(int));

    for (i = 0; i < newSizeCols; i++)
        colsPointers[0][i] = tempCols[i];

    for (i = 0; i < newSizeRows; i++)
        rowsIndexes[0][i] = tempRows[i];

    free(tempRows);
    free(tempCols);
    return numColsRemoved;
}
