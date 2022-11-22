#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cilk/cilk.h>
#include "mmio.h"


int BFS(int *colors, int *rowsIncomingEdges, int *colsIncomingEdges, int u, int numVertices, int *vertices, int *SSCIDs)
{
    int i = -1, j, size = 0;
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
        size++;
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

    return size;
}

int *coloringSCC(int *rowsOutgoingEdges, int *colsOutgoingEdges, int *rowsIncomingEdges, int *colsIncomingEdges, int numVertices)
{
    int *colors = (int *)malloc(numVertices * sizeof(int));
    int *vertices = (int *)malloc(numVertices * sizeof(int));
    int *SCCIDs = (int *)calloc(numVertices, sizeof(int));
    int *verticesRemoved = (int *)calloc(numVertices, sizeof(int));

    int verticesRemaining = numVertices;
    int i;
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
            int *colorChanges = (int *)calloc(numVertices, sizeof(int));
            colorChange = 0;
            cilk_for(int j = 0; j < numVertices; j++)
            {
                if (vertices[j] == 0) // if vertices[j] == 0, vertice with id = j + 1 has been removed
                    continue;

                for (int k = rowsOutgoingEdges[j]; k < rowsOutgoingEdges[j + 1]; k++)
                {
                    if (colors[j] > colors[colsOutgoingEdges[k]] && vertices[colsOutgoingEdges[k]] != 0)
                    {
                        colors[colsOutgoingEdges[k]] = colors[j];
                        colorChanges[j] = 1;
                    }
                }
            }
            for (i = 0; i < numVertices; i++){
                if (colorChanges[i] == 1){
                    colorChange = 1;
                    break;
                }
            }
            free(colorChanges);
        } while (colorChange);

        cilk_for(int j = 0; j < numVertices; j++)
        {
            // check only vertices that kept their original color
            if (colors[j] != j + 1)
                continue;

            int size = BFS(colors, rowsIncomingEdges, colsIncomingEdges, j, numVertices, vertices, SCCIDs);
            verticesRemoved[j] += size;
        }
        for (i = 0; i < numVertices; i++)
        {
            verticesRemaining -= verticesRemoved[i];
            verticesRemoved[i] = 0;
            if (verticesRemaining == 0)
                break;
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

    printf("Rows: %d\nColumns: %d\nNon-zero: %d\n", rows, cols, nz);

    clock_t time;

    time = clock();
    int *SCCIDs = coloringSCC(CSR_ROW_INDEX, CSR_COL_INDEX, CSC_ROW_INDEX, CSC_COL_INDEX, rows);
    time = clock() - time;

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


    printf("SCCs: %d, Trivial: %d\n", CC, TRIVIAL);
    printf("Time taken: %f\n", ((double) time) / CLOCKS_PER_SEC);

    return 0;
}
