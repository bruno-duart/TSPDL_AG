#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "grafos.h"

typedef struct {
    int *solution;
    int distance;
}Solution;

int* ini_array(int dimension);
void print_arr(int dimension, int *arr);
Solution* new_solution(int vsize);
void free_solution(Solution *S);
Solution* greed_method(Graph *G, int *demand, int *draft);
bool is_Solution(Graph *G, int *demand, int *draft, int *solution);
Solution* random_swap(Graph *G, int *demand, int *draft, int* solution, int distance, int param);
Solution* random_solution(Graph *G, int *demand, int *draft);
int fitness(Graph *G, int *S);
int isIn(Graph *G, int cidade, int *filho);
int indexOf(Graph *G, Solution **Arr, int value);
void order1Crossover(Graph *G, Solution** Pai, int **filho);
Solution genetico(Graph *G, int *demand, int *draft);

void merge(Solution **Arr, int start, int middle, int end);
void mergeSort(Solution **Arr, int start, int end);
int torneio(Graph *G, Solution **population);
Solution* construcao(Graph *G, int *demand, int *draft);
Solution AlgGenetico(Graph *G, int *demand, int *draft);

int main(){    
    //Inicialização da matriz e dos vetores
    clock_t tempo1, tempo2;
    int dimension, value, distance, *solution;
    Solution  melhor;
    scanf("%d", &dimension);
    Graph *G = New_Graph(dimension);
    //Graph_print(G);    
    int *demand, *draft;
    int weight = G->V - 1;
    demand = ini_array(dimension);
    draft = ini_array(dimension);
         
    srand(time(NULL));
    tempo1 = clock();
    /* melhor =*/ AlgGenetico(G, demand, draft);
    
    /*
    tempo2 = clock();
    printf("%d\t\t%f\n",melhor.distance, (double)(tempo2 - tempo1)/CLOCKS_PER_SEC);*/
    
    free_Graph(G);
    free(demand);
    free(draft);
    
    return 0;
}

Solution* new_solution(int vsize){
    Solution* sol = malloc(sizeof(Solution));
    sol->solution = malloc(sizeof(int) * vsize);
    sol->distance = 0;
    return sol;
}

void free_solution(Solution *S){
    free(S->solution);
    free(S);
}

int* ini_array(int dimension){
    int *arr = malloc(sizeof(int) * dimension);
    int value;
    for(int i=0; i < dimension; i++){
        scanf("%d", &value);
        arr[i] = value;
    }
    return arr;
}

void print_arr(int dimension, int *arr){
    for(int i=0; i < dimension; i++)
        printf(" %d ",arr[i]); //erro com valor não inicializado, por conta do 'i'
    printf("\n");
}

Solution* greed_method(Graph *G, int *demand, int *draft){
    Solution *s = malloc(sizeof(Solution));
    s->solution = malloc(sizeof(int)*G->V);
    s->distance = 0;
    int menor, weight = G->V-1;
    int position = 0, k=0;
    int *demand_1 = malloc(sizeof(int)*G->V);
    for(int i=0; i<G->V; i++)
        demand_1[i] = demand[i];

    for(int i = 0; i < G->V; i++){
        menor = 9999;
        for(int j = 0; j < G->V; j++){
            if(G->adj[k][j] < menor && demand_1[j] == 1 && weight <= draft[j]){
                position = j;
                menor = G->adj[k][j];
            }
        }
        if(menor != 9999){
            k = position;
            s->solution[i] = k;
            weight--;
            s->distance += menor;
            demand_1[k] = 0;
        }
    }
    s->distance += G->adj[k][0];
    return s;
}

bool is_Solution(Graph *G, int *demand, int *draft, int *solution){
    int weight = G->V-1;
    int *demand_1 = malloc(sizeof(int)*G->V);
    for(int i=0; i<G->V; i++)
        demand_1[i] = demand[i];
    
    for(int i=0; i < G->V; i++){
        if(solution[i] == 0 && i != (G->V-1)){
            free(demand_1);
            return false;
        }
        if(draft[solution[i]] < (weight-i)){
            free(demand_1);
            return false;
        }
        demand_1[solution[i]] = 0;
    }
    for(int i=0; i < G->V; i++)
        if(demand_1[i] != 0){
            free(demand_1);
            return false;
        }
    free(demand_1);
    return true;
}

Solution* random_solution(Graph *G, int *demand, int *draft){
    Solution *s = malloc(sizeof(Solution));
    s->solution = malloc(sizeof(int)*G->V);
    s->distance = 0;
    s->solution[G->V-1] = 0;
    
    do{
        for(int i=0; i<G->V-1; i++)
        {
            s->solution[i] = rand() % G->V;
            for(int j=0; j<i; j++)
                if(s->solution[j] == s->solution[i])
                {
                    i--;
                    break;
                }
        }
    }while(!is_Solution(G, demand, draft, s->solution));

    s->distance += G->adj[0][s->solution[0]];
        for(int j=1; j < G->V; j++)
            s->distance += G->adj[s->solution[j-1]][s->solution[j]];
    
    return s;
}

Solution* random_swap(Graph *G, int *demand, int *draft, int* solution, int distance, int param){
    int index_1, index_2, aux, distance_i;
    int *copy = malloc(sizeof(int)*(G->V));
    Solution *s = malloc(sizeof(Solution));
    s->solution = malloc(sizeof(int) * (G->V));
    s->distance = distance;
    
    for(int i=0; i < 100; i++){
        for(int k=0; k < G->V; k++)
            copy[k] = solution[k];
        do{
            index_1 = rand() % (G->V-1);
            do{
                index_2 = rand() % (G->V-1);
            }while(index_1 == index_2);

            aux = copy[index_2];
            copy[index_2] = copy[index_1];
            copy[index_1] = aux;
        }while(!is_Solution(G, demand, draft, copy));
        
        distance_i = G->adj[0][copy[0]];
        for(int j=1; j < G->V; j++)
            distance_i += G->adj[ copy[j-1] ][ copy[j] ];

        if(distance_i < s->distance){
            s->distance = distance_i;
            for(int k=0; k < G->V; k++)
                s->solution[k] = copy[k];
        }
    }
    free(copy);
    return s;
}
int fitness(Graph *G, int *S){
    int distance = G->adj[0][S[0]];

    for(int i = 1; i < G->V; i++){
        distance += G->adj[S[i-1]][S[i]];
    }
    return distance;
}

int isIn(Graph *G, int cidade, int *filho){
    for(int i = 0; i < G->V; i++)
        if(filho[i] == cidade)
            return 1;
    return 0;
}

int indexOf(Graph *G, Solution **Arr, int value){
    for(int i = 0; i < 2 * G->V; i++)
        if(Arr[i]->distance == value)
            return i;
    return -1;
}

void order1Crossover(Graph *G,  Solution** Pai, int **filho){
    int  index, cidade;    
    int i, j, k, n1, n0;

    i = rand() % (G->V - 1);
    do{
        j = rand() % (G->V - 1);
    }while(i == j);

    if(j < i){
        k = i;
        i = j;
        j = k;
    }
    
    //preenchimento do segmento originado dos pais
    for(k = 0; k < G->V; k++){
        filho[0][k] = (k >= i && k <= j ? Pai[0]->solution[k] : -1);
        filho[1][k] = (k >= i && k <= j ? Pai[1]->solution[k] : -1);
    }
    
    n0 = n1 = G->V-1;
    for(k = G->V-1; k >= 0; k--){
        if(k == j)
            k = i-1;
        if(n0 == j)
            n0 = i-1;
        if(n1 == j)
            n1 = i-1;
        //preenchimento dos trechos a esquerda de j e a direita de i do filho 0
        if(!isIn(G,Pai[1]->solution[k], filho[0]) && filho[0][n0]==-1)
            filho[0][n0--] = Pai[1]->solution[k];
        //preenchimento dos trechos a esquerda de j e a direita de i do filho 1
        if(!isIn(G,Pai[0]->solution[k], filho[1]) && filho[1][n1]==-1)
            filho[1][n1--] = Pai[0]->solution[k];
    }
    
    if(n0 > 0)
        for(k = G->V-1; k >= 0; k--){
            if(filho[0][k] == -1){
                for(int l = G->V-1; l >= 0; l--){
                    if(!isIn(G,Pai[1]->solution[l], filho[0])){
                        filho[0][k] = Pai[1]->solution[l];
                        break;
                    }
                }
            }
        }
    
    if(n1 > 0)
        for(k = G->V-1; k >= 0; k--){
            if(filho[1][k] == -1){
                for(int l = G->V-1; l >= 0; l--){
                    if(!isIn(G,Pai[0]->solution[l], filho[1])){
                        filho[1][k] = Pai[0]->solution[l];
                        break;
                    }
                }
            }
        }
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

void shuffle(Graph *G, Solution **Arr, int nTrocas){
    Solution *aux;
    int idx1, idx2;

    for(int i = 0; i < nTrocas; i++){
        idx1 = rand() % (2 * G->V - 1);
        do{
            idx2 = rand() % (2 * G->V - 1);
        }while(idx1 == idx2);
        aux = Arr[idx1];
        Arr[idx1] = Arr[idx2];
        Arr[idx2] = aux;
    }
}

bool indice(int size, int value, int *arr){
    for(int i = 0; i < size; i++)
        if(arr[i] == value)
            return true;
    return false;
}

int torneio(Graph *G, Solution **p){
    int index1, index2, idx, numCandidatos = 2*G->V;
    Solution *looser;

    while(numCandidatos > G->V){
        index1 = rand() % (numCandidatos);
        do{
            index2 = rand() % (numCandidatos);
        }while(index1 == index2);

        idx = (p[index1]->distance > p[index2]->distance ? index1 : index2);
        looser = p[idx];

        for(int j = idx; j < 2*G->V - 1; j++)
            p[j] = p[j+1];
        p[2*G->V-1] = looser;
        numCandidatos--;
    }
    return numCandidatos;
}

void torneio2(Graph *G, Solution **p, int nFighter, Solution **pais){
    //Retorna os índices dos pais para cruzamento
    int *visitado = malloc(sizeof(int) * 2 * G->V);
    int *candidatos = malloc(sizeof(int) * nFighter);
    int idx, idx1, idx2, nCandidato = nFighter, i;
    //inicialização
    for(i = 0; i < 2*G->V; i++)
        visitado[i] = 0;
    //seleção de candidatos
    for(i = 0; i < nFighter; i++){
        do{
            idx = rand() % (2 * G->V);
        }while(visitado[idx]);
        candidatos[i] = idx;
        visitado[idx] = 1;
    }
    //reset
    for(i = 0; i < 2*G->V; i++)
        visitado[i] = 0;
    //torneio até sobrarem apenas 1 família tradicional brasileira
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
    idx = 0;
    for(i = 0; i < nCandidato && idx < 2; i++){
        if(!visitado[candidatos[i]]){
            for(int j = 0; j < G->V; j++)
                pais[idx]->solution[j] = p[candidatos[i]]->solution[j];
            pais[idx]->distance = p[candidatos[i]]->distance;
            idx++;
        }
    }
    
    free(visitado);
    free(candidatos);
}

void print_population(Solution** p, int size_p, int size_v){
        for(int i = 0; i < size_p ; i++){
        printf("[Cromossomo %2d - Distância %4d]  : ", i, p[i]->distance);
        print_arr(size_v, p[i]->solution);
    }
}

Solution* construcao(Graph *G, int *demand, int *draft){
    /*Constrói uma solução viável para o conjunto de soluções relativo à 
    população inicial do algoritmo genético. A construção é realizada
    com o método guloso com características aleatórias, através do qual
    dois portos são sorteados aleatoriamente para comporem o início da solução, 
    e em seguida, aplica-se o método guloso para obter uma solução viável.*/

    Solution *novaSolucao = new_solution(G->V);
    int port1, port2, weight = G->V-1, menor, position, k;
    int *demand_1 = malloc(sizeof(int) * G->V);
    
    //Cópia do array de demandas, permitindo que seja possível alterá-lo sem preocupações
    for(int i=0; i<G->V; i++)
        demand_1[i] = demand[i];

    
    //Seleção dos dois portos iniciais. Port1 e port2 devem ser diferentes entre si
    //e devem ser viáveis.
    do{
        port1 = rand() % (G->V-1) + 1;
    }while(draft[port1] < (weight));

    do{
        port2 = rand() % (G->V-1) + 1;
    }while(draft[port2] < (weight - 1) || port2 == port1);

    //Inicialização da solução, atualização da distância percorrida
    novaSolucao->solution[0] = port1;
    novaSolucao->solution[1] = port2;
    novaSolucao->distance = G->adj[0][port1] + G->adj[port1][port2];
    demand_1[port1] = demand_1[port2] = 0;
    weight -= 2;

    //printf("Porto 1: %d\nPorto 2: %d\n", port1, port2);
    k = port2;
    //printf("Weight: %d\n",weight);

    //Etapa gulosa da construção.
    for(int i = 2; i < G->V; i++){
        menor = __INT16_MAX__;
        for(int j = 0; j < G->V; j++){
            if(G->adj[k][j] < menor && demand_1[j] == 1 && weight <= draft[j]){
                position = j;
                menor = G->adj[k][j];
            }
        }        
        if(menor != __INT16_MAX__){
            //printf("K = %d\n", position);
            //printf("%d\n",i);
            k = position;
            novaSolucao->solution[i] = k;
            weight--;
            novaSolucao->distance += menor;
            demand_1[k] = 0;
        }
       // print_arr(G->V, novaSolucao->solution);
    }
    novaSolucao->distance += G->adj[k][0];
    //printf("Weight: %d\n",weight);
    //exit(1);

    free(demand_1);
    return novaSolucao;    
}

void mutacao(Graph *G, int *demand, int *draft, Solution *P){
    int ind1, ind2, auxTroca;

    do{
        ind1 = rand() % (G->V - 1);
        do{
            ind2 = rand() % (G->V - 1);
        }while(ind1 == ind2);

        auxTroca = P->solution[ind1];
        P->solution[ind1] = P->solution[ind2];
        P->solution[ind2] = auxTroca;

        if(is_Solution(G, demand, draft, P->solution))
            break;
        
        //Caso a troca gere uma solução não viável, desfaz-se a troca.
        auxTroca = P->solution[ind1];
        P->solution[ind1] = P->solution[ind2];
        P->solution[ind2] = auxTroca;
    }while(1);

    P->distance = fitness(G, P->solution);
}

void copiar(Graph *G, Solution *S, int *solucao){
    for(int i = 0; i < G->V; i++)
        S->solution[i] = solucao[i];
    S->distance = fitness(G, solucao);
}

Solution AlgGenetico(Graph *G, int *demand, int *draft){
    int sizep = 2 * G->V, sizer = G->V;
    int i, **bebes = malloc(sizeof(int*) * 2), numFilhos;
    int *mutantes, gerAtual = 0, ultMelhor = 0;
    Solution **population = malloc(sizeof(Solution*) * sizep);
    Solution **pais = malloc(sizeof(Solution*) * 2);
    Solution **filhos = malloc(sizeof(Solution*) * sizep);
    Solution *Melhor = new_solution(sizer);
    
    pais[0] = new_solution(sizer);
    pais[1] = new_solution(sizer);
    bebes[0] = malloc(sizeof(int) * sizer);
    bebes[1] = malloc(sizeof(int) * sizer);

    //print_arr(sizer, draft);
   
    for(i = 0; i < sizep; i++){
        population[i] = construcao(G, demand, draft);
        if(i == 0 || population[i]->distance < Melhor->distance)
            copiar(G, Melhor, population[i]->solution);
        filhos[i] = new_solution(sizer);
    }
        
    //print_population(population, sizep, sizer);   
    printf("Menor: %d\n", Melhor->distance);
    numFilhos = 0;
    
    while((gerAtual - ultMelhor) < 1){
        //Cruzamento
        for(i = 0; i < sizer; i++){
            if((rand() % 100) < 95){ // Probabilidade de ocorrer 'crossover'

                //Seleção dos pais
                //printf("Pais selecionados\n");
                torneio2(G, population, 6, pais);

                //print_population(pais, 2, sizer);            
                do{
                    order1Crossover(G, pais, bebes);
                }while(!is_Solution(G, demand, draft, bebes[0]) || !is_Solution(G, demand, draft, bebes[1]));

                // printf("Filho 1: Distância: %d ",fitness(G,bebes[0]));
                // print_arr(sizer, bebes[0]);
                // printf("Filho 2: Distância: %d ",fitness(G,bebes[1]));
                // print_arr(sizer, bebes[1]);

                copiar(G, filhos[numFilhos++], bebes[0]);
                copiar(G, filhos[numFilhos++], bebes[1]);
            }
        }
        
        //Mutação
        for(i = 0; i < numFilhos; i++){
            if((rand() % 100) < 5 ){
                i = rand() % numFilhos;
                // printf("Antes da mutação: \nDistância: %d: ", population[i]->distance);
                // print_arr(sizer, population[i]->solution);    
                mutacao(G, demand, draft, filhos[i]);
                // printf("Depois da mutação: \nDistância: %d: ", population[i]->distance);
                // print_arr(sizer, population[i]->solution);    
                // mutacao(G, demand, draft, population[i]);
            }
        }
        //printf("Antes Filhos\n");
        //print_population(filhos, numFilhos, sizer);

        while(numFilhos < sizep){
            i = rand() % sizep;
            copiar(G, filhos[numFilhos++], population[i]->solution);
        }
        
        for(i = 0; i < sizep; i++){
            copiar(G, population[i], filhos[i]->solution);
            
            if(population[i]->distance < Melhor->distance){
                copiar(G, Melhor, population[i]->solution);
                ultMelhor = gerAtual;
            }
        }
        
        gerAtual++;
    }
    printf("Menor: %d\n", Melhor->distance);
    //printf("Depois Filhos\n");
    //print_population(filhos, numFilhos, sizer);

    //Liberação de memória alocada
    for(i = 0; i < sizep; i++){
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
}


// Primeira tentativa de implementar algoritmos genéticos. 
// Não houve sucesso. Parte das funções serão reaproveitadas.

/* Solution genetico(Graph *G, int *demand, int *draft){
    Solution *population[2*G->V];
    Solution *pais[2];
    Solution *filhos[2*G->V];
    Solution melhor;
    int iterAtual = 0, ultMelhor = 0, index1, index2, numBebes, tamPopulacao;
    int auxTroca, prob, nFighter;
    int *bebes[2];
    bool flag = false;

    //bebes = malloc(sizeof(int*) * 2);
    bebes[0] = malloc(sizeof(int) * G->V);
    bebes[1] = malloc(sizeof(int) * G->V);

    for(int i=0; i< G->V; i++){
        bebes[0][i] = 0;
        bebes[1][i] = 0;
    }

    pais[0] = new_solution(G->V);
    pais[1] = new_solution(G->V);

    // Geração de população inicial de soluções
    for(int i = 0; i < 2 * G->V; i++){
        population[i] = random_solution(G, demand, draft);
        //printf("Solução %d gerada\n",i);
        if(i == 0 || population[i]->distance < melhor.distance)
            melhor = *population[i];
    }    
   // printf("Melhor distância: %d\n",melhor.distance);
    for(int i = 0; i < 2 * G->V; i++)
        filhos[i] = new_solution(G->V);

    while((iterAtual - ultMelhor)<100){
       
        tamPopulacao = torneio(G,population);
        numBebes = 0;
        //printf("\nIteração %d\n", iterAtual);

        for(int k = 0; k < G->V/2 && numBebes < G->V; k++){
            // Crossover
            prob = rand() % 100;
            
            if(prob < 98){
                    
                do{
                    index1 = rand() % tamPopulacao;
                    do{
                        index2 = rand() % tamPopulacao;
                    }while(index1 == index2);

                    
                    pais[0]->solution = population[index1]->solution;
                    pais[1]->solution = population[index2]->solution;
                
                    order1Crossover(G, pais, bebes);
                }while(!is_Solution(G,demand,draft, bebes[0]) || !is_Solution(G,demand,draft,bebes[1]));
                numBebes += 2;

                for(int i = 0; i < G->V; i++){
                    filhos[numBebes-2]->solution[i] = bebes[0][i];
                    filhos[numBebes-1]->solution[i] = bebes[1][i];
                }
                
                filhos[numBebes-2]->distance = fitness(G,bebes[0]);
                filhos[numBebes-1]->distance = fitness(G,bebes[1]);
            }
            
        }

        for(int i = 0; i < numBebes; i++)
            population[tamPopulacao + i] = filhos[i];

        for(int i = 0; i < 2*G->V; i++){
            prob = rand() % 100;
            if(prob < 5){

                index1 = rand() % numBebes;
                do{
                
                    index2 = rand() % (G->V-2);

                    auxTroca = population[index1]->solution[index2];
                    population[index1]->solution[index2] = population[index1]->solution[index2-1];
                    population[index1]->solution[index2-1] = auxTroca;

                    population[index1]->distance = fitness(G, population[index1]->solution);

                    if(is_Solution(G, demand, draft, population[index1]->solution))
                        flag = false;
                    
                    else{
                        auxTroca = population[index1]->solution[index2];
                        population[index1]->solution[index2] = population[index1]->solution[index2-1];
                        population[index1]->solution[index2-1] = auxTroca;

                        flag = true;
                        population[index1]->distance = fitness(G, population[index1]->solution);
                    }
                }while(flag);
            }
        } 
            
        for(int i = 0; i < 2 * G->V; i++){
           if(population[i]->distance < melhor.distance){
                melhor = *population[i];
                ultMelhor = iterAtual;
           }
        }
        iterAtual++;
    }
   
    // Liberação de memória alocada
    
    free(pais[0]->solution);
    free(pais[1]->solution);
    free(pais[0]);
    free(pais[1]);
    free(bebes[0]);
    free(bebes[1]);
    for(int i = 0; i < 2 * G->V; i++){
        free(filhos[i]->solution);
        free(filhos[i]);
        //printf("%p\n", population[i]->solution);
        //free(population[i]->solution);
        //free(population[i]);
    }
    //free(population);
    //free(population[0]);
    
    return melhor;
}   */

