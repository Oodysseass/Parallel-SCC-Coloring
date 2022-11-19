#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"

int findSize(int *Vc, int *rowsIncomingEdges, int *colsIncomingEdges, int u)
{
    int size = 0, i;
    for (i = colsIncomingEdges[u]; i < colsIncomingEdges[u + 1]; i++) // all vertices for which there is an edge <v, u>
    {
        if (Vc[rowsIncomingEdges[i]] == 1) // vertices for which there is an edge <v, u> AND have the same color as u
            size++;
    }
    return size;
}

void BFS(int *Vc, int *rowsIncomingEdges, int *colsIncomingEdges, int u, int *SCVc)
{
    int i, j = 0;

    for (i = colsIncomingEdges[u]; i < colsIncomingEdges[u + 1]; i++)
    {
        if (Vc[rowsIncomingEdges[i]] == 1)
        {
            SCVc[j++] = rowsIncomingEdges[i];
        }
    }
}

// Coloring Algorithm
int *ColoringSCC(int *rowsOutgoingEdges, int *colsOutgoingEdges, int *rowsIncomingEdges, int *colsIncomingEdges, int numVertices)
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
            for (i = 0; i < numVertices; i++)
            {
                if (vertices[i] == 0) // if vertices[i] == 0, vertice with id = i + 1 has been removed
                    continue;
                for (j = rowsOutgoingEdges[i]; j < rowsOutgoingEdges[i + 1]; j++)
                {
                    if (colors[i] > colors[colsOutgoingEdges[j]] && vertices[colsOutgoingEdges[j]] != 0)
                    {
                        colors[colsOutgoingEdges[j]] = colors[i];
                        colorChange = 1;
                    }
                }
            }
        } while (colorChange);

        for (i = 0; i < numVertices; i++)
        {
            // check only vertices that kept their original color
            if (colors[i] != i + 1)
                continue;

            int *Vc = (int *)calloc(numVertices, sizeof(int)); // worst case scenario, all vertices have the same color
            for (j = 0; j < numVertices; j++)
            {
                if (colors[j] == colors[i])
                    Vc[j] = 1;
            }

            size = findSize(Vc, rowsIncomingEdges, colsIncomingEdges, i);
            int SCVc[size + 1];
            if (size != 0)
                BFS(Vc, rowsIncomingEdges, colsIncomingEdges, i, &SCVc[0]);
            SCVc[size] = i;

            for (j = 0; j < size + 1; j++)
            {
                SCCIDs[SCVc[j]] = id; // vertices of the SCVc have the same id
                vertices[SCVc[j]] = 0; // remove vertices of the SCVc from initial V
                verticesRemaining--;
            }
            id++;
            free(Vc);
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

    // printf("Lines: %d\nColumns: %d\nNon-zero: %d\n", rows, cols, nz);

    // CSR data structure
    int value;
    int CSR_ROW_INDEX[rows + 1];
    int CSR_COL_INDEX[nz];
    int readRow;

    for (i = 0; i < rows + 1; i++)
        CSR_ROW_INDEX[i] = 0;

    for (i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %d\n", &value, &readRow, &CSR_COL_INDEX[i]);
        CSR_COL_INDEX[i]--;
        CSR_ROW_INDEX[readRow]++;
    }

    for (i = 0; i < rows; i++)
        CSR_ROW_INDEX[i + 1] += CSR_ROW_INDEX[i];

    // CSC data structure
    int CSC_ROW_INDEX[nz];
    int CSC_COL_INDEX[cols + 1];
    int col, dest, temp, last = 0, cumsum = 0;

    for (i = 0; i < nz; i++)
        CSC_COL_INDEX[i] = 0;
    // CONVERT CSR TO CSC
    for (i = 0; i < nz; i++)
        CSC_COL_INDEX[CSR_COL_INDEX[i]]++;

    for (i = 0; i < cols; i++)
    {
        temp = CSC_COL_INDEX[i];
        CSC_COL_INDEX[i] = cumsum;
        cumsum += temp;
    }
    CSC_COL_INDEX[cols] = nz;

    for (i = 0; i < rows; i++)
    {
        for (j = CSR_ROW_INDEX[i]; j < CSR_ROW_INDEX[i + 1]; j++)
        {
            col = CSR_COL_INDEX[j];
            dest = CSC_COL_INDEX[col];

            CSC_ROW_INDEX[dest] = i;

            CSC_COL_INDEX[col]++;
        }
    }

    for (i = 0; i < cols + 1; i++)
    {
        temp = CSC_COL_INDEX[i];
        CSC_COL_INDEX[i] = last;
        last = temp;
    }

    int *SCCIDs = ColoringSCC(CSR_ROW_INDEX, CSR_COL_INDEX, CSC_ROW_INDEX, CSC_COL_INDEX, rows);

/*     int* showTimes = (int *) calloc(rows, sizeof(int));

    for (i = 0; i < rows; i++) {
        showTimes[SCCIDs[i]]++;
    }
    int sum = 0;
    for (i = 0; i < rows; i++) {
        if (showTimes[i] != 1)
            sum++;
    }
    printf("%d\n", sum); */

    return 0;
}
