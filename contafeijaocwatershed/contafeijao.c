/*=============================================================
 *           UNIFAL = Universidade Federal de Alfenas.
 *             BACHARELADO EM CIENCIA DA COMPUTACAO.
 * Trabalho..: Contador de Feijoes
 * Professor.: Luiz Eduardo da Silva
 * Aluno.....: Matheus Reis de Lima
 * Data......: 18/06/2024
 *=============================================================*/

// Libs
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Luizinho
#include "utils/imagelib.h"

// Inicio
image limiarizar_pbm(image In, int limiar)
{
    image Out = img_clone(In);
    for (int i = 0; i < In->nc * In->nr; i++)
    {
        if (In->px[i] >= limiar)
        {
            Out->px[i] = 0;
        }
        else
        {
            Out->px[i] = 1;
        }
    }
    return Out;
}

int find(int parent[], int i)
{
    while (parent[i] != i)
        i = parent[i];
    return i;
}

void Union(int parent[], int i, int j)
{
    int x = find(parent, i);
    int y = find(parent, j);
    parent[y] = x;
}

int quantidadeRotulos(int *parent, int numLabel)
{
    int qtdcomponentes = 0;
    int pertence;
    for (int i = 1; i <= numLabel; i++)
    {
        pertence = 0;
        for (int j = 0; j <= numLabel; j++)
        {
            if (i == parent[j])
            {
                pertence = 1;
                break;
            }
        }
        if (pertence)
        {
            qtdcomponentes++;
        }
    }
    return qtdcomponentes;
}

int label(image In)
{
    int nr = In->nr;
    int nc = In->nc;
    int *px = In->px;
    int numLabel = 0;
    int *parent = (int *)malloc(nc * nr * sizeof(int));
    for (int i = 0; i < nr * nc; i++)
    {
        parent[i] = i;
    }
    for (int i = 1; i < nr; i++)
        for (int j = 1; j < nc; j++)
        {
            int p = px[i * nc + j];
            int r = px[(i - 1) * nc + j];
            int t = px[i * nc + j - 1];
            if (p != 0)
            {
                if (r == 0 && t == 0)
                {
                    px[i * nc + j] = ++numLabel;
                }
                    
                if (r != 0 && t == 0)
                {
                    px[i * nc + j] = r;
                }
                    
                if (r == 0 && t != 0)
                {
                    px[i * nc + j] = t;
                }

                if (r != 0 && t != 0 && t == r)
                {
                    px[i * nc + j] = r;
                    // printf("Teste");
                }

                if (r != 0 && t != 0 && t != r)
                {
                    // printf("Teste");
                    px[i * nc + j] = t;
                    Union(parent, r, t);
                }
            }
        }
    for (int i = 0; i < nr * nc; i++)
    {
        In->px[i] = find(parent, In->px[i]);
    }
        
    In->ml = numLabel;
    return quantidadeRotulos(parent, numLabel);
}

void msg(char *s)
{
    printf("Uso:  %s  nome-imagem[.pgm] \n", s);
    printf("nome-imagem[.pgm] é o nome do arquivo da imagem \n");
    printf("beans[] \n");
    exit(1);
}

/*------------------------------------------------
 * Fila de prioridade
 *------------------------------------------------*/
typedef struct no *ptno;
typedef struct no
{
    int i, j;
    ptno next;
} no;

/*-----------------------------------------------
 *  init   *Q      new     init         *Q
 *  [a:]->[b:]-+   [c:]  ->[a:]->[b:]->[c:]-+
 *   ^---------+             ^--------------+
 *----------------------------------------------*/
void insQ(ptno *Q, int i, int j)
{
    ptno new = malloc(sizeof(no));
    new->i = i;
    new->j = j;
    if (!(*Q))
        new->next = new;
    else
    {
        new->next = (*Q)->next;
        (*Q)->next = new;
    }
    *Q = new;
}

/*-----------------------------------------------
 *  init         *Q                    *Q
 *  [a:]->[b:]->[c:]-+    [a:]  [b:]->[c:]-+
 *    ^--------------+            ^--------+
 *----------------------------------------------*/
void remQ(ptno *Q, int *i, int *j)
{
    ptno init = (*Q)->next;
    *i = init->i;
    *j = init->j;
    if (*Q == init)
        *Q = NULL;
    else
        (*Q)->next = init->next;
    free(init);
}

int isEmpty(ptno *Q)
{
    return *Q == NULL;
}

void initQPrior(ptno *QPrior, int mn)
{
    int i;
    for (i = 0; i < mn; i++)
        QPrior[i] = NULL;
}

void insert(ptno *QPrior, int i, int j, int p)
{
    insQ(QPrior + p, i, j);
}

int pop(ptno *QPrior, int *i, int *j, int *maxPrior, int mn)
{
    while (*maxPrior < mn && isEmpty(QPrior + *maxPrior))
        (*maxPrior)++;
    if (*maxPrior == mn)
        return 1;
    remQ(QPrior + *maxPrior, i, j);
    return 0;
}

/*------------------------------------------------
 * Calcula o gradient da imagem PGM
 *------------------------------------------------*/
image gradient(image In, int raio)
{
    int i, j, y, x, max, min;
    int nl = In->nr, nc = In->nc, mn = In->ml;
    image Out = img_clone(In);
    for (i = 0; i < nl; i++)
        for (j = 0; j < nc; j++)
        {
            max = -1;
            min = mn + 1;
            for (y = -raio; y <= raio; y++)     // -1 0 1
                for (x = -raio; x <= raio; x++) // -1 0 1
                {
                    int pi = i + y;
                    int pj = j + x;
                    if (pi >= 0 && pi < nl && pj >= 0 && pj < nc)
                    {
                        if (abs(x) + abs(y) <= raio && In->px[pi * nc + pj] > max)
                            max = In->px[pi * nc + pj];
                        if (abs(x) + abs(y) <= raio && In->px[pi * nc + pj] < min)
                            min = In->px[pi * nc + pj];
                    }
                }
            Out->px[i * nc + j] = max - min;
        }
    return Out;
}

/*------------------------------------------------
 * Calcula o Watershed da imagem
 *------------------------------------------------*/
image watershed(image In, int y, int x)
{
    int i, j, k, maxPrior = 0, stop = 0;
    int nl = In->nr, nc = In->nc, mn = In->ml;
    ptno qPrior[mn];
    image mark = img_clone(In);
    image Out = img_clone(In);
    enum
    {
        NONE,
        QUEUE,
        WSHED,
        MARK1,
        MARK2,
    };
    struct
    {
        int i, j;
    } n4[4] =
        {0, 1,   // vizinho da direita
         1, 0,   // vizinho de baixo
         0, -1,  // vizinho da esquerda
         -1, 0}; // vizinho de cima

    initQPrior(qPrior, mn);

    // Inicialização dos marcadores
    //-----------------------------
    for (i = 0; i < nl * nc; i++)
        mark->px[i] = NONE;

    // MARK1 = (y, x) - centro da região
    int raio = 10;
    for (i = -raio; i <= raio; i++)
        for (j = -raio; j <= raio; j++)
        {
            int pi = i + y;
            int pj = j + x;
            if (abs(i) + abs(j) <= raio)
                mark->px[pi * nc + pj] = MARK1;
        }

    // mark[y * nc + x] = MARK1;

    // MARK2 = borda da imagem
    for (i = 0; i < nl; i++)
    {
        mark->px[i * nc] = MARK2;
        mark->px[i * nc + nc - 1] = MARK2;
    }
    for (j = 0; j < nc; j++)
    {
        mark->px[j] = MARK2;
        mark->px[(nl - 1) * nc + j] = MARK2;
    }
    // Inicialização
    //--------------
    // Todos os vizinhos dos marcadores são colocados na fila
    // A prioridade é o nível de cinza do pixel
    for (i = 1; i < nl - 1; i++)
        for (j = 1; j < nc - 1; j++)
            if (mark->px[i * nc + j] == NONE)
            {
                int isAdj = 0;
                for (k = 0; k < 4; k++)
                {
                    int pi = i + n4[k].i;
                    int pj = j + n4[k].j;
                    int m = mark->px[pi * nc + pj];
                    if (m == MARK1 || m == MARK2)
                        isAdj = 1;
                }
                if (isAdj)
                {
                    mark->px[i * nc + j] = QUEUE;
                    insert(qPrior, i, j, In->px[i * nc + j]);
                }
            }
    // Crescimento dos Marcadores
    //---------------------------
    // Enquanto a fila de prioridade não está vazia faça
    // 1. retira um pixel da fila
    // 2. Se o pixel eh vizinho de UMA marca, recebe essa marca.
    //    Seus vizinhos sem rotulo são colocados na fila
    // 3. Se o pixel eh vizinho de DUAS marcas, ele eh WATERSHED
    while (!stop)
    {
        int m = NONE;
        int isWshed = 0;
        stop = pop(qPrior, &i, &j, &maxPrior, mn);
        if (!stop)
        {
            for (k = 0; k < 4; k++)
            {
                int pi = i + n4[k].i;
                int pj = j + n4[k].j;
                if (pi >= 0 && pi < nl && pj >= 0 && pj < nc)
                {
                    int mAdj = mark->px[pi * nc + pj];
                    if (mAdj == MARK1 || mAdj == MARK2)
                    {
                        if (m == NONE)
                            m = mAdj;
                        else if (m != mAdj)
                            isWshed = 1;
                    }
                }
            }
            if (isWshed)
                mark->px[i * nc + j] = WSHED;
            else
            {
                mark->px[i * nc + j] = m;
                for (k = 0; k < 4; k++)
                {
                    int pi = i + n4[k].i;
                    int pj = j + n4[k].j;
                    if (pi >= 0 && pi < nl && pj >= 0 && pj < nc)
                        if (mark->px[pi * nc + pj] == NONE)
                        {
                            int prior, px;
                            mark->px[pi * nc + pj] = QUEUE;
                            px = In->px[pi * nc + pj];
                            prior = (px < maxPrior) ? maxPrior : px;
                            insert(qPrior, pi, pj, prior);
                        }
                }
            }
        }
    }
    for (i = 0; i < nl * nc; i++)
        Out->px[i] = (mark->px[i] == WSHED) ? 0 : 255;
    img_free(mark);
    return Out;
}

int main(int argc, char *argv[])
{
    int nc, nr, ml, tp = GRAY, qComponentes;
    char nameIn[100], nameOut[100], cmd[110];
    image In, Out, watershedOut;

    if (argc < 2)
    {
        msg(argv[0]);
    }

    img_name(argv[1], nameIn, nameOut, tp, tp - 1);

    In = img_get(nameIn, tp);
    Out = limiarizar_pbm(In, 94);

    // Aplicar o Watershed
    watershedOut = watershed(Out, In->nr / 2, In->nc / 2);

    // Contar os componentes
    qComponentes = label(watershedOut);

    // Salvando a imagem do Watershed
    img_put(watershedOut, "watershed_output.pgm", tp - 1);

    printf("#componentes= %2.0d\n", qComponentes);

    img_free(In);
    img_free(Out);
    img_free(watershedOut);

    return 0;
}


