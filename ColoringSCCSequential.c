#include <stdio.h>
#include <stdlib.h>
#include "mmio.h"

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
                printf("neighboooooor %d %d\n", j, v);
                queue[++i] = colors[rowsIncomingEdges[j]];
                SCVc[size++] = rowsIncomingEdges[j];
                visited[rowsIncomingEdges[j]] = 1;
            }
        }
    }

    free(queue);
    free(visited);

    return size;
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

    //printf("Lines: %d\nColumns: %d\nNon-zero: %d\n", rows, cols, nz);

    // CSC data structure
    int value;
    int* CSC_COL_INDEX = (int *) malloc((cols + 1) * sizeof(int));
    int* CSC_ROW_INDEX = (int *) malloc(nz * sizeof(int));
    int readCol;

    for (i = 0; i < cols + 1; i++)
        CSC_COL_INDEX[i] = 0;

    for (i = 0; i < nz; i++)
    {
        fscanf(f, "%d %d %d\n", &CSC_ROW_INDEX[i], &readCol, &value);
        CSC_ROW_INDEX[i]--;
        CSC_COL_INDEX[readCol]++;
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

    int *SCCIDs = coloringSCC(CSR_ROW_INDEX, CSR_COL_INDEX, CSC_ROW_INDEX, CSC_COL_INDEX, rows);
/*
    int* showTimes = (int *) calloc(rows + 1, sizeof(int));

    for (i = 0; i < rows + 1; i++) {
        showTimes[SCCIDs[i]]++;
    }
    int sum = 0, trivial = 0;
    for (i = 0; i < rows; i++) {
        if (showTimes[i] > 1)
            sum++;
        else if (showTimes[i] == 1)
            trivial++;
    }
    printf("%d %d\n", sum, trivial); */

    return 0;
}
