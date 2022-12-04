#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mmio.h"

void BFS(int *colors, int *rowsIncomingEdges, int *colsIncomingEdges, int u, int numVertices, int *vertices, int *SSCIDs)
{
    int i = -1, j;
    int *queue = (int *)calloc(numVertices, sizeof(int));
    int *visited = (int *)calloc(numVertices, sizeof(int));
    int v;
    queue[++i] = u;
    visited[u] = 1;

    while (i != -1)
    {
        v = queue[i--];
        vertices[v] = 0;
        SSCIDs[v] = u + 1;
        for (j = colsIncomingEdges[v]; j < colsIncomingEdges[v + 1]; j++)
        {
            if (colors[rowsIncomingEdges[j]] == u + 1 && visited[rowsIncomingEdges[j]] != 1)
            {
                queue[++i] = rowsIncomingEdges[j];
                visited[rowsIncomingEdges[j]] = 1;
            }
        }
    }

    free(queue);
    free(visited);
}

// coloring algorithm
int *coloringSCC(int *rowsOutgoingEdges, int *colsOutgoingEdges, int *rowsIncomingEdges, int *colsIncomingEdges, int numVertices)
{
    int *colors = (int *)malloc(numVertices * sizeof(int));
    int *vertices = (int *)malloc(numVertices * sizeof(int));
    int *verticesRemaining = (int *)malloc(numVertices * sizeof(int));
    int *SCCIDs = (int *)malloc(numVertices * sizeof(int));
    int *uniqueColors = (int *)malloc(numVertices * sizeof(int));

    int numUnique;
    int numVerticesRemaining = numVertices;
    int colorChange;
    int i, j;

    for (i = 0; i < numVertices; i++)
    {
        vertices[i] = i + 1;
        verticesRemaining[i] = i + 1;
    }

    while (numVerticesRemaining != 0)
    {
        for (i = 0; i < numVertices; i++)
        {
            colors[i] = vertices[i];
        }

        do
        {
            colorChange = 0;
            for (i = 0; i < numVerticesRemaining; i++)
            {
                for (j = rowsOutgoingEdges[verticesRemaining[i] - 1]; j < rowsOutgoingEdges[verticesRemaining[i]]; j++)
                {
                    if (colors[verticesRemaining[i] - 1] > colors[colsOutgoingEdges[j]] && vertices[colsOutgoingEdges[j]] != 0)
                    {
                        colors[colsOutgoingEdges[j]] = colors[verticesRemaining[i] - 1];
                        colorChange = 1;
                    }
                }
            }
        } while (colorChange);

        numUnique = 0;
        for (i = 0; i < numVertices; i++)
            if (colors[i] == i + 1)
                uniqueColors[numUnique++] = i;

        for (i = 0; i < numUnique; i++)
            BFS(colors, rowsIncomingEdges, colsIncomingEdges, uniqueColors[i], numVertices, vertices, SCCIDs);

        numVerticesRemaining = 0;
        for (i = 0; i < numVertices; i++)
            if (vertices[i] != 0)
                verticesRemaining[numVerticesRemaining++] = vertices[i];
    }

    free(colors);
    free(vertices);
    free(verticesRemaining);
    free(uniqueColors);
    return SCCIDs;
}

int trimming(int **indexes, int **nonzeros, int size, int nz);

void convert(int *indexA, int *nzA, int *indexB, int *nzB, int size, int nz);

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

    printf("--SEQUENTIAL--\n");
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

    fclose(f);

    int *CSR_ROW_INDEX = (int *)malloc((rows + 1) * sizeof(int));
    int *CSR_COL_INDEX = (int *)malloc(nz * sizeof(int));

    struct timeval start, end;
    long seconds, microseconds;
    double elapsed;
    int trimmed = 0, removed = 0;

    gettimeofday(&start, NULL);

    // trim indegree zero
    removed = trimming(&CSC_COL_INDEX, &CSC_ROW_INDEX, cols, nz);
    trimmed += removed;
    cols = cols - removed;
    rows = cols;
    nz = CSC_COL_INDEX[cols];

    // trim outdegree zero
    CSR_COL_INDEX = (int *)realloc(CSR_COL_INDEX, nz * sizeof(int));
    convert(CSC_COL_INDEX, CSC_ROW_INDEX, CSR_ROW_INDEX, CSR_COL_INDEX, cols, nz);

    removed = trimming(&CSR_ROW_INDEX, &CSR_COL_INDEX, cols, nz);
    trimmed += removed;
    rows = rows - removed;
    cols = rows;
    nz = CSR_ROW_INDEX[rows];

    // remake CSC data structure
    CSC_ROW_INDEX = (int *)realloc(CSC_ROW_INDEX, nz * sizeof(int));
    convert(CSR_ROW_INDEX, CSR_COL_INDEX, CSC_COL_INDEX, CSC_ROW_INDEX, rows, nz);
    gettimeofday(&end, NULL);

    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    elapsed = seconds + microseconds * 1e-6;
    printf("Time taken for trimming: %f\n", elapsed);

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

    seconds = end.tv_sec - start.tv_sec;
    microseconds = end.tv_usec - start.tv_usec;
    elapsed = seconds + microseconds * 1e-6;
    printf("Time taken: %f\n", elapsed);

    return 0;
}

// converts CSC/CSR(A) to CSR/CSC(B)
void convert(int *indexA, int *nzA, int *indexB, int *nzB, int size, int nz)
{
    int i, j, row, dest, temp, last = 0, cumsum = 0;

    for (i = 0; i < size + 1; i++)
        indexB[i] = 0;

    for (i = 0; i < nz; i++)
        indexB[nzA[i]]++;

    for (i = 0; i < size; i++)
    {
        temp = indexB[i];
        indexB[i] = cumsum;
        cumsum += temp;
    }
    indexB[size] = nz;

    for (i = 0; i < size; i++)
    {
        for (j = indexA[i]; j < indexA[i + 1]; j++)
        {
            row = nzA[j];
            dest = indexB[row];

            nzB[dest] = i;

            indexB[row]++;
        }
    }

    for (i = 0; i < size + 1; i++)
    {
        temp = indexB[i];
        indexB[i] = last;
        last = temp;
    }
}

int trimming(int **indexes, int **nonZeros, int size, int nz)
{
    int *numEdges = (int *)malloc(size * sizeof(int));
    int *tempIndexes = (int *)malloc((size + 1) * sizeof(int)); // plus 1 because of the size = cols and not actual size of indexes[0]
    int newNZ = 0;
    int newSize = 0;
    int i, j;

    // find number of edges of each vertice, if it is zero make it -1
    // that is because we want to remove only vertices that had from start 0 in/outdegree
    // and not because of the trimming
    for (i = 0; i < size; i++)
    {
        if (indexes[0][i] != indexes[0][i + 1])
            numEdges[i] = indexes[0][i + 1] - indexes[0][i];
        else
            numEdges[i] = -1;
    }

    // erase edges from/to zero in/outdegree vertices
    for (i = 0; i < size; i++)
    {
        for (j = indexes[0][i]; j < indexes[0][i + 1]; j++)
        {
            if (numEdges[nonZeros[0][j]] == -1)
            {
                nonZeros[0][j] = -1;
                numEdges[i]--;
            }
        }
    }

    // fix new num of edges of each vertice
    indexes[0][0] = 0;
    for (i = 0; i < size; i++)
    {
        if (numEdges[i] != -1)
            indexes[0][i + 1] = indexes[0][i] + numEdges[i];
        else
            indexes[0][i + 1] = indexes[0][i];
    }

    // remove vertices with zero in/outdegree
    for (i = 0; i < size; i++)
    {
        if (numEdges[i] != -1)
            tempIndexes[newSize++] = indexes[0][i];
    }
    tempIndexes[newSize++] = indexes[0][size];

    // keep edges that are only from non-zero in/outdegree vertices
    int *tempNonzeros = (int *)malloc(indexes[0][size] * sizeof(int));
    for (i = 0; i < nz; i++)
    {
        if (nonZeros[0][i] != -1)
            tempNonzeros[newNZ++] = nonZeros[0][i];
    }

    // we have to fix the nonZeros array because of the vertices we removed
    int removed = 0;
    for (i = 0; i < size; i++)
    {
        if (numEdges[i] != -1)
            numEdges[i] = removed;
        else
            removed++;
    }

    for (i = 0; i < newNZ; i++)
    {
        tempNonzeros[i] -= numEdges[tempNonzeros[i]];
    }

    // transfer new structure to arguments
    indexes[0] = (int *)realloc(indexes[0], newSize * sizeof(int));
    nonZeros[0] = (int *)realloc(nonZeros[0], newNZ * sizeof(int));

    for (i = 0; i < newSize; i++)
        indexes[0][i] = tempIndexes[i];

    for (i = 0; i < newNZ; i++)
        nonZeros[0][i] = tempNonzeros[i];

    free(tempIndexes);
    free(tempNonzeros);
    free(numEdges);
    return (size + 1 - newSize);
}