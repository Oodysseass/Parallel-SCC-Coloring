#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cilk/cilk.h>
#include "mmio.h"

int CC = 0;
int trivial = 0;

int BFS(int *SCVc, int *colors, int *rowsIncomingEdges, int *colsIncomingEdges, int u, int numVertices)
{
    int i = -1, j, size = 0;
    int* queue = (int *) calloc(numVertices, sizeof(int));
    int* visited = (int *) calloc(numVertices, sizeof(int));
    int v;
    queue[++i] = u;
    visited[u] = 1;
    SCVc[size++] = u;

    while (i != -1) {
        v = queue[i--];
        for (j = colsIncomingEdges[v]; j < colsIncomingEdges[v + 1]; j++){
            if (colors[rowsIncomingEdges[j]] == u + 1 && visited[rowsIncomingEdges[j]] != 1) {
                queue[++i] = rowsIncomingEdges[j];
                SCVc[size++] = rowsIncomingEdges[j];
                visited[rowsIncomingEdges[j]] = 1;
            }
        }
    }

    free(queue);
    free(visited);

    return size;
}

// coloring algorithm
int *coloringSCC(int *rowsOutgoingEdges, int *colsOutgoingEdges, int *rowsIncomingEdges, int *colsIncomingEdges, int numVertices)
{
    int *colors = (int *)malloc(numVertices * sizeof(int));
    int *vertices = (int *)malloc(numVertices * sizeof(int));
    int *SCCIDs = (int *)malloc(numVertices * sizeof(int));

    int verticesRemaining = numVertices;
    int colorChange;
    int i, j, id = 1, size;

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
            colorChange = 0;
            cilk_for(int k = 0; k < numVertices; k++)
            {
                if (vertices[k] == 0) // if vertices[k] == 0, vertice with id = k + 1 has been removed
                    continue;
                for (j = rowsOutgoingEdges[k]; j < rowsOutgoingEdges[k + 1]; j++)
                {
                    if (colors[k] > colors[colsOutgoingEdges[j]] && vertices[colsOutgoingEdges[j]] != 0)
                    {
                        colors[colsOutgoingEdges[j]] = colors[k];
                        colorChange = 1;
                    }
                }
            }
        } while (colorChange);

        cilk_for(int k = 0; k < numVertices; k++)
        {
            // check only vertices that kept their original color
            if (colors[k] != k + 1)
                continue;

            int* SCVc = (int *) calloc(numVertices, sizeof(int));
            size = BFS(SCVc, colors, rowsIncomingEdges, colsIncomingEdges, k, numVertices);
            if (size > 1)
                CC++;
            else
                trivial++;

            for (j = 0; j < size; j++)
            {
                SCCIDs[SCVc[j]] = id;   // vertices of the SCVc have the same id
                vertices[SCVc[j]] = 0;  // remove vertices of the SCVc from initial V
                verticesRemaining--;
            }

            id++;
            free(SCVc);
        }
    }

    return SCCIDs;
}


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


    // CSC data structure
    int* CSC_COL_INDEX = (int *) malloc((cols + 1) * sizeof(int));
    int* CSC_ROW_INDEX = (int *) malloc(nz * sizeof(int));
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
    } else {
        for (i = 0; i < nz; i++)
        {
            fscanf(f, "%d %d\n", &CSC_ROW_INDEX[i], &readCol);
            CSC_ROW_INDEX[i]--;
            CSC_COL_INDEX[readCol]++;
        }
    }

    for (i = 0; i < cols; i++)
        CSC_COL_INDEX[i + 1] += CSC_COL_INDEX[i];

    // CSR data structure
    int* CSR_COL_INDEX = (int *) malloc(nz * sizeof(int));
    int* CSR_ROW_INDEX = (int *) malloc((rows + 1) * sizeof(int));
    int row, dest, temp, last = 0, cumsum = 0;

    for (i = 0; i < rows + 1; i++){
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

    printf("Rows: %d\nColumns: %d\nNon-zero: %d\n", rows, cols, nz);

    clock_t time;

    time = clock();
    int *SCCIDs = coloringSCC(CSR_ROW_INDEX, CSR_COL_INDEX, CSC_ROW_INDEX, CSC_COL_INDEX, rows);
    time = clock() - time;

    printf("SCCs: %d, Trivial: %d\n", CC, trivial);
    printf("Time taken: %f\n", ((double) time) / CLOCKS_PER_SEC);

    return 0;
}
