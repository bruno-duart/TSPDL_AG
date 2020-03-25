#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "grafos.h"

typedef struct {
    int *harbor;
    int distance;
}Solution;

//variáveis globais
Graph *G;
size_t CONT_GER;
int DIM, PSIZE, MAX_ITER, PERC_MUT, OPT_VAL;
int *DEMAND, *DRAFT;

int* ini_array();
void print_arr(int *arr);
Solution* new_solution();
void free_solution(Solution *S);
bool is_Solution(int *harbor);
Solution* random_swap(Solution* individuo);
int fitness(int *S);
int isIn(int cidade, int *filho);
int indexOf(Solution **Arr, int value);
void order1Crossover(Solution** Pai, int **filho);
void merge(Solution **Arr, int start, int middle, int end);
void mergeSort(Solution **Arr, int start, int end);

Solution* construcao();
int AlgGenetico();

int main(){
    clock_t tempo1, tempo2;
    srand(time(NULL));
    //Inicialização da matriz e dos vetores
    scanf("%d", &DIM);
    G = New_Graph(DIM);
    DEMAND = ini_array(DIM);
    DRAFT = ini_array(DIM);
    scanf("%d", &OPT_VAL);
    //parâmetros do algoritmo genético
    PSIZE = DIM * 2;
    MAX_ITER = 100;
    PERC_MUT = 5;


    tempo1 = clock();
    //sei la

    long int result, best, cbest, media[2];
    int qt_iter = 1;

    for(int i=1; i <= 20; i++){
        PERC_MUT = i;
        media[0] = media[1] = 0;
        best = 0;
        for(int j=0; j < qt_iter; j++){
            result = AlgGenetico();
            media[0] += result;
            media[1] += (long int) CONT_GER;
            if(result < best || !best){
                best = result;
                cbest = 1;
            }
            else if(result == best)
                cbest++;
        }
        media[0] /= qt_iter;
        media[1] /= qt_iter;
        cbest *= (100.0 / qt_iter);
        printf("[ %2i %% ]  MediaRes = %5ld   MediaCont = %5ld   "\
                "BestRes = %5ld ( %2i %%)  DesvMed = %2.2f  DesvMenor = %2.2f\n", i, media[0], 
                media[1], best, (int) cbest,(media[0]-OPT_VAL)*100.0/(OPT_VAL), (best-OPT_VAL)*100.0/(OPT_VAL));
    }
    //printf("\nCONT_GER = %ld\n\n", (long int) CONT_GER);
    
    free_Graph(G);
    free(DEMAND);
    free(DRAFT);
    
    return 0;
}

Solution* new_solution(){
    Solution* sol = malloc(sizeof(Solution));
    sol->harbor = malloc(sizeof(int) * DIM);
    sol->distance = __INT16_MAX__;

    for(int i=0; i < DIM; i++)
        sol->harbor[i] = 0;
    return sol;
}

void free_solution(Solution *S){
    free(S->harbor);
    free(S);
}

int* ini_array(){
    int *arr = malloc(sizeof(int) * DIM);
    int value;
    for(int i=0; i < DIM; i++)
        scanf("%d", &arr[i]);
    return arr;
}

void print_arr(int *arr){
    for(int i=0; i < DIM; i++)
        printf(" %d ",arr[i]); //erro com valor não inicializado, por conta do 'i'
    printf("\n");
}

bool is_Solution(int *harbor){
    int weight = G->V-1;
    int demand[DIM];

    //copia a demanda
    for(int i=0; i < G->V; i++)
        demand[i] = DEMAND[i];
    
    for(int i=0; i < G->V - 1; i++){
        //verifica o indice dos portos
        if(harbor[i] < 0)
            return false;
        //verifica se o porto de origem está sendo visitado antes da hora
        if(harbor[i] == 0)
            return false;
        //verifica a condição do calado
        if(DRAFT[harbor[i]] < weight)
            return false;
        //verifica se a demanda já foi atendida
        if(!demand[harbor[i]])
            return false;
        //senão, atende à demanda e diminui o peso
        demand[harbor[i]] = 0;
        weight--;
    }
    //verifica se o navio retorna ao porto de origem
    return (harbor[DIM-1] == 0);
}

int fitness(int *S){
    int distance = G->adj[0][ S[0] ];

    for(int i = 1; i < G->V; i++){
        distance += G->adj[ S[i-1] ][ S[i] ];
    }
    return distance;
}

int isIn(int cidade, int *filho){
    for(int i = 0; i < G->V; i++)
        if(filho[i] == cidade)
            return 1;
    return 0;
}

int indexOf(Solution **Arr, int value){
    for(int i = 0; i < 2 * G->V; i++)
        if(Arr[i]->distance == value)
            return i;
    return -1;
}

void order1Crossover(Solution** Pai, int **filho){
    int  index, cidade;    
    int i, j, k;
    int n[2];

    //sorteio da região (intervalo) de corte
    i = rand() % (G->V - 2) + 1;
    do
        j = rand() % (G->V - 2) + 1;
    while(i == j);

    if(j < i){
        k = i;
        i = j;
        j = k;
    }
    //preenchimento do segmento originado dos pais
    for(k = 0; k < G->V; k++){
        filho[0][k] = (k >= i && k <= j ? Pai[0]->harbor[k] : -1);
        filho[1][k] = (k >= i && k <= j ? Pai[1]->harbor[k] : -1);
    }
    
    n[0] = n[1] = G->V - 1;

    for(k = G->V-1; k >= 0; k--){
        if(k == j)      //salto do intervalo de análise
            k = i-1;
        if(n[0] == j)     //salto do índice de cópia 0
            n[0] = i-1;
        if(n[1] == j)     //salto do índice de cópia 1
            n[1] = i-1;
        //preenchimento dos trechos a direita de j e a esquerda de i do filho 0
        if(!isIn(Pai[1]->harbor[k], filho[0]) && filho[0][ n[0] ]==-1)
            filho[0][ n[0]-- ] = Pai[1]->harbor[k];
        //preenchimento dos trechos a direita de j e a esquerda de i do filho 1
        if(!isIn(Pai[0]->harbor[k], filho[1]) && filho[1][ n[1] ]==-1)
            filho[1][ n[1]-- ] = Pai[0]->harbor[k];
    }
    //caso falte posições no filho a serem preenchidas
    for(int m=0; m<2; m++)
        if(n[m] >= 0)
            for(k = G->V-1; k >= 0; k--){
                if(filho[m][k] == -1){
                    for(int l = G->V-1; l >= 0; l--){
                        if(!isIn(Pai[1-m]->harbor[l], filho[m])){
                            filho[m][k] = Pai[1-m]->harbor[l];
                            break;
                        }
                    }
                }
            }
            
    //teste de erro
    if(isIn(-1, filho[0]) || isIn(-1, filho[1])){
        printf("[ %i %i ]\n",i,j);
        print_arr(Pai[0]->harbor);
        print_arr(Pai[1]->harbor);
        printf("\n");
        print_arr(filho[0]);
        print_arr(filho[1]);
        printf("----------------------------------------------\n");
    }
}

bool indice(int size, int value, int *arr){
    for(int i = 0; i < size; i++)
        if(arr[i] == value)
            return true;
    return false;
}

void prova_do_ortiz(Solution **p, int nFighter, Solution **pais){
    //Retorna os índices dos pais para cruzamento
    int visitado[PSIZE];
    int candidatos[nFighter];
    int idx, idx1, idx2, nCandidato = nFighter, i;
    //inicialização
    for(i = 0; i < PSIZE; i++)
        visitado[i] = 0;
    //seleção de candidatos
    for(i = 0; i < nFighter; i++){
        do
            idx = rand() % PSIZE;
        while(visitado[idx]);
        candidatos[i] = idx;
        visitado[idx] = 1;
    }
    //reset
    for(i = 0; i < PSIZE; i++)
        visitado[i] = 0;
    //torneio até sobrarem apenas 2 universitários fudidos
    while(nFighter > 2){
        //seleção aleatória de 2 candidatos ainda
        //do
            idx1 = rand() % (nCandidato);
       // while(visitado[candidatos[idx1]]);
        if(visitado[candidatos[idx1]])
            while(visitado[candidatos[ (++idx1)%nCandidato ]]);
        idx1 = idx1 % nCandidato;
        do{
            idx2 = rand() % (nCandidato);
            if(visitado[candidatos[idx2]])
                while(visitado[candidatos[ (++idx2)%nCandidato ]]);
            idx2 = idx2 % nCandidato;
        }while(idx1 == idx2);
        //escolha do perdedor (pior candidato)
        idx = (p[candidatos[idx1]]->distance > p[candidatos[idx2]]->distance ? idx1 : idx2);
        //eliminação do perdedor
        visitado[candidatos[idx]] = 1;
        nFighter--;
    }
    //atribuição dos selecionados ao vetor de pais
    for(i = 0, idx = 0; i < nCandidato && idx < 2; i++)
        if(!visitado[candidatos[i]]){
            for(int j = 0; j < G->V; j++)
                pais[idx]->harbor[j] = p[candidatos[i]]->harbor[j];
            pais[idx]->distance = p[candidatos[i]]->distance;
            idx++;
        }
}

void print_population(Solution** p, int size_p){
    for(int i = 0; i < size_p ; i++){
        printf("[Cromossomo %2d - Distância %4d]  : ", i, p[i]->distance);
        print_arr(p[i]->harbor);
    }
}

Solution* construcao(){
    /*Constrói uma solução viável para o conjunto de soluções relativo à 
    população inicial do algoritmo genético. A construção é realizada
    com o método guloso com características aleatórias, através do qual
    dois portos são sorteados aleatoriamente para comporem o início da solução, 
    e em seguida, aplica-se o método guloso para obter uma solução viável.*/

    Solution *novaSolucao = new_solution(G->V);
    int port1, port2, weight = G->V-1, menor, position, k;
    int *demand = malloc(sizeof(int) * G->V);
    
    //Cópia do array de demandas, permitindo que seja possível alterá-lo sem preocupações
    for(int i=0; i<G->V; i++)
        demand[i] = DEMAND[i];
    
    //Seleção dos dois portos iniciais. Port1 e port2 devem ser diferentes entre si
    //e devem ser viáveis.
    do{
        port1 = rand() % (G->V-1) + 1;
    }while(DRAFT[port1] < (weight));
    weight--;

    do{
        port2 = rand() % (G->V-1) + 1;
    }while(DRAFT[port2] < weight || port2 == port1);
    weight--;

    //Inicialização da solução, atualização da distância percorrida
    novaSolucao->harbor[0] = port1;
    novaSolucao->harbor[1] = port2;
    novaSolucao->distance = G->adj[0][port1] + G->adj[port1][port2];
    demand[port1] = demand[port2] = 0;

    //printf("Porto 1: %d\nPorto 2: %d\n", port1, port2);
    k = port2;
    //printf("Weight: %d\n",weight);

    //Etapa gulosa da construção.
    for(int i = 2; i < G->V; i++){
        menor = __INT16_MAX__;
        for(int j = 0; j < G->V; j++){
            if(G->adj[k][j] < menor && demand[j] == 1 && weight <= DRAFT[j]){
                position = j;
                menor = G->adj[k][j];
            }
        }        
        if(menor != __INT16_MAX__){
            //printf("K = %d\n", position);
            k = position;
            novaSolucao->harbor[i] = k;
            weight--;
            novaSolucao->distance += menor;
            demand[k] = 0;
        }
       // print_arr(G->V, novaSolucao->harbor);
    }
    novaSolucao->distance += G->adj[k][0];
    //print_arr(G->V, novaSolucao->harbor);
    //printf("Weight: %d\n",weight);
    //exit(1);

    free(demand);
    return novaSolucao;    
}

void mutacao(Solution *P){
    /*Responsável por realizar mutação nos cromossomos, permitindo maior
    variabilidade no espaço de busca. Sorteiam-se dois índices e trocam-se os dois
    de locus. Verifica-se a viabilidade da solução, e caso positivo, aceita-se a troca.
    Caso contrário, retorna-se ao estado original e refaz-se o sorteio*/
    int ind1, ind2, auxTroca;

    do{
        ind1 = rand() % (G->V - 1);
        do{
            ind2 = rand() % (G->V - 1);
        }while(ind1 == ind2);

        auxTroca = P->harbor[ind1];
        P->harbor[ind1] = P->harbor[ind2];
        P->harbor[ind2] = auxTroca;

        if(is_Solution(P->harbor))
            break;
        
        //Caso a troca gere uma solução não viável, desfaz-se a troca.
        auxTroca = P->harbor[ind1];
        P->harbor[ind1] = P->harbor[ind2];
        P->harbor[ind2] = auxTroca;
    }while(1);

    P->distance = fitness(P->harbor);
}

void copiar(Solution *S, int *solucao){
    /*Função utilizada para efetuar a cópia do conteúdo dos ponteiros de 
    Solution*/
    for(int i = 0; i < G->V; i++)
        S->harbor[i] = solucao[i];
    S->distance = fitness(solucao);
}

void merge(Solution **Arr, int start, int middle, int end){
    Solution **temp = malloc(sizeof(Solution) * (end - start + 1)); 
    int i = start, j = middle + 1, k = 0;

    while(i <= middle && j <= end){
        if(Arr[i]->distance <= Arr[j]->distance){
            temp[k] = Arr[i];
            k++;
            i++;
        }
        else{
            temp[k] = Arr[j];
            k++;
            j++;
        }
    }
    while(i <= middle){
        temp[k] = Arr[i];
        k++;
        i++;
    }
    while(j <= end){
        temp[k] = Arr[j];
        k++;
        j++;
    }
    for(int i = start; i <= end; i++)
        Arr[i] = temp[i - start];
    
    //free(temp);
}

void mergeSort(Solution **Arr, int start, int end){
    if (start < end){
        int middle = start + (end - start)/2;
        mergeSort(Arr, start, middle);
        mergeSort(Arr, middle+1, end);
        merge(Arr, start, middle, end);
    }
}

int AlgGenetico(){
    int nFighter;
    int taxaMutacao;
    int i, j, **bebes = malloc(sizeof(int*) * 2), numFilhos;
    int *mutantes, gerAtual = 0, ultMelhor = 0;
    int index[4];
    Solution **population = malloc(sizeof(Solution*) * PSIZE);
    Solution **pais = malloc(sizeof(Solution*) * 2);
    Solution **filhos = malloc(sizeof(Solution*) * PSIZE);
    Solution *Melhor = new_solution();
    
    pais[0] = new_solution();
    pais[1] = new_solution();
    bebes[0] = malloc(sizeof(int) * DIM);
    bebes[1] = malloc(sizeof(int) * DIM);

    for(i = 0; i < PSIZE; i++){
        population[i] = construcao();
        if(i == 0 || population[i]->distance < Melhor->distance)
            copiar(Melhor, population[i]->harbor);
        filhos[i] = new_solution();
    }

    //início do ciclo
    CONT_GER = 0;
    while((gerAtual - ultMelhor) < MAX_ITER){
        if((gerAtual - ultMelhor) == (MAX_ITER/2)){
            for(i = 0; i < PSIZE; i++)
                copiar(population[i], random_swap(population[i])->harbor);
        }
        //Cruzamento
        numFilhos = 0;
        for(i = 0; i < DIM; i++)
            if((rand() % 100) < 95){ // Probabilidade de ocorrer 'crossover'
                //Seleção dos pais
                nFighter = rand() % (DIM/2) + (DIM/4);
                prova_do_ortiz(population, nFighter, pais);

                do
                    order1Crossover(pais, bebes);
                while(!is_Solution(bebes[0]) || !is_Solution(bebes[1]));

                copiar(filhos[numFilhos++], bebes[0]);
                copiar(filhos[numFilhos++], bebes[1]);
            }
        
        //Mutação (taxa entre 0 e 20%)
        taxaMutacao = ((gerAtual - ultMelhor) * PERC_MUT) / MAX_ITER;
        for(i = 0; i < numFilhos; i++){
            if((rand() % 100) < taxaMutacao){
                i = rand() % numFilhos;
                mutacao(filhos[i]);
            }
        }

        while(numFilhos < PSIZE){ //Preenche o conjunto de cromossomos com pais aleatórios
            i = rand() % PSIZE;
            copiar(filhos[numFilhos++], population[i]->harbor);
        }

        //0 e 1 para os melhores filhos
        //2 e 3 para os piores pais
        index[0] = index[2] = 0; 
        index[1] = index[3] = 1;

        if(filhos[ index[0] ]->distance > filhos[ index[1] ] ->distance){
            i = index[0]; 
            index[0] = index[1];
            index[1] = i;
        }

        for(i = 2; i < PSIZE; i++){
            if(filhos[i]->distance < filhos[ index[0] ]->distance){
                index[1] = index[0];
                index[0] = i;
            } else if(filhos[i]->distance < filhos[ index[1] ]->distance){
                index[1] = i;
            }
        }

        if(population[ index[2] ]->distance < population[ index[3] ] ->distance){
            i = index[2]; 
            index[2] = index[3];
            index[3] = i;
        }
        for(i = 2; i < PSIZE; i++){
            if(population[i]->distance > population[ index[2] ]->distance){
                index[3] = index[2];
                index[2] = i;
            } else if(population[i]->distance > population[ index[3] ]->distance ){
                index[3] = i;
            }
        }
        /*
        for(i = 0; i < PSIZE; i++){ //Atualiza as gerações
            copiar(population[i], filhos[i]->harbor);
            
            if(population[i]->distance < Melhor->distance){ //Atualiza a melhor solução corrente
                copiar(Melhor, population[i]->harbor);
                ultMelhor = gerAtual;
            }
        }*/

        if(filhos[ index[0] ]->distance < population[ index[2] ]->distance){
            copiar(population[ index[2] ], filhos[ index[0] ]->harbor);
            if(filhos[ index[1] ]->distance < population[ index[3] ]->distance)
                copiar(population[ index[3] ], filhos[ index[1] ]->harbor);
        }else if(filhos[ index[1] ]->distance < population[ index[2] ]->distance){
            copiar(population[ index[2] ], filhos[ index[1] ]->harbor);
        }else if(filhos[ index[0] ]->distance < population[ index[3] ]->distance){
            copiar(population[ index[3] ], filhos[ index[0] ]->harbor);
        }else if(filhos[ index[1] ]->distance < population[ index[3] ]->distance){
            copiar(population[ index[3] ], filhos[ index[1] ]->harbor);
        }
        /*
        copiar(population[ index[2] ], filhos[ index[0] ]->harbor);
        copiar(population[ index[3] ], filhos[ index[1] ]->harbor);
        */
        for(i = 0; i < PSIZE; i++){ //Atualiza as gerações
           
            if(population[i]->distance < Melhor->distance){ //Atualiza a melhor solução corrente
                copiar(Melhor, population[i]->harbor);
                ultMelhor = gerAtual;
            }
        }
        gerAtual++;
        CONT_GER++;
        //print_arr(Melhor->harbor);
    }

    //Liberação de memória alocada
    for(i = 0; i < PSIZE; i++){
        free_solution(population[i]);
        free_solution(filhos[i]);
    }
    
    free(bebes[0]);
    free(bebes[1]);
    free(bebes);
    free_solution(pais[0]);
    free_solution(pais[1]);
    free(pais);
    free(population);
    free(filhos);

    int res = Melhor->distance;
    free_solution(Melhor);
    return res;
}

Solution* random_swap(Solution* individuo){
    int index_1, index_2, aux, distance_i;
    int *copy = malloc(sizeof(int)*(G->V));
    Solution *s = malloc(sizeof(Solution));
    s->harbor = malloc(sizeof(int) * (G->V));
    s->distance = individuo->distance;
    
    for(int i=0; i < 100; i++){
        //gera uma cópia da solução de entrada
        for(int k=0; k < G->V; k++)
            copy[k] = individuo->harbor[k];
        //permanece no loop até encontrar uma solução válida
        while(1){
            index_1 = rand() % (G->V-1);
            do{
                index_2 = rand() % (G->V-1);
            }while(index_1 == index_2);

            aux = copy[index_2];
            copy[index_2] = copy[index_1];
            copy[index_1] = aux;

            if(is_Solution(copy))
                break;

            aux = copy[index_2];
            copy[index_2] = copy[index_1];
            copy[index_1] = aux;
        }
        //calcula o custo
        
        distance_i = G->adj[0][copy[0]];
        for(int j=1; j < G->V; j++)
            distance_i += G->adj[ copy[j-1] ][ copy[j] ];

        if(distance_i < s->distance){
            s->distance = distance_i;
            for(int k=0; k < G->V; k++)
                s->harbor[k] = copy[k];
        }
    }
    free(copy);
    return s;
}