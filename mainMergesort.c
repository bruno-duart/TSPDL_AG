#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <string.h>
#include "grafos.h"

typedef struct {
    int *harbor;
    int distance;
}Solution;

//variáveis globais
Graph *G;
size_t CONT_GER;
int DIM, PSIZE, MAX_ITER, NUM_TESTES, PERC_MUT, OPT_VAL;
int *DEMAND, *DRAFT;

//cabeçalhos de função
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
void fixed_swap(Solution* s);
void merge(Solution **Arr, int start, int middle, int end);
void mergeSort(Solution **Arr, int start, int end);
void shuffle(Solution **Arr, int nTrocas);
Solution* construcao();
int AlgGenetico();

int main(){
    clock_t tempo, menortempo;
    srand(time(NULL));
    //Inicialização dos valores
    scanf("%d", &DIM);
    G = New_Graph(DIM);
    DEMAND = ini_array(DIM);
    DRAFT = ini_array(DIM);
    scanf("%d", &OPT_VAL);
    PSIZE = DIM * 2;
    MAX_ITER = 20;
    NUM_TESTES = 100;

    printf("Ótimo Conhecido: %d\n", OPT_VAL);
    printf("Número de Repetições do AG: %d\n", NUM_TESTES);
    printf("Máximo de Gerações sem Melhora: %d\n", MAX_ITER);
    long int result, best, cbest, media[3];
    
    //variação da taxa de mutação
    for(int i=0; i <= 100; i+=5){
        PERC_MUT = i;
        media[0] = media[1] = media[2] = 0;
        best = 0;
        //ciclo de execuções
        for(int j=0; j < NUM_TESTES; j++){
            tempo = clock();
            result = AlgGenetico();
            tempo = clock() - tempo;
            media[0] += result;
            media[1] += (long int) CONT_GER;
            media[2] += (long int) tempo;
            if(result < best || !best){
                best = result;
                cbest = 1;
                menortempo = tempo;
            }
            else if(result == best)
                cbest++;
        }
        //cálculo e impressão dos resultados
        media[0] /= NUM_TESTES;
        media[1] /= NUM_TESTES;
        media[2] /= NUM_TESTES;
        cbest *= 100.0 / (float) NUM_TESTES;
        printf("[%3i %% ]  MediaRes = %5ld   MediaCont = %4ld   "\
                "BestRes = %5ld (%3i %%)  ErroMed = %5.2f  ErroMenor = %5.2f  Tempo = %lf  MenTempo = %lf\n", i, media[0], 
                (media[1]-MAX_ITER), best, (int) cbest,(media[0]-OPT_VAL)*100.0/(OPT_VAL), (best-OPT_VAL)*100.0/(OPT_VAL),
                (double) media[2]/CLOCKS_PER_SEC, (double) menortempo/CLOCKS_PER_SEC);
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

    //copia o vetor de demanda
    for(int i=0; i < G->V; i++)
        demand[i] = DEMAND[i];
    
    for(int i=0; i < G->V - 1; i++){
        //verifica o indice dos portos
        if(harbor[i] < 0)
            return false;
        //verifica se o porto de origem está sendo marcado antes da hora
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

    for(int i = 1; i < G->V; i++)
        distance += G->adj[ S[i-1] ][ S[i] ];
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
    //preenchimento do intervalo com o segmento dos pais
    /*
                 [0]            [i]     [j]             [k]
                
    	Pai A:     a0  a1 ... ai-1 ai ... aj aj+1 ... ak-1 ak
    	Pai B:     b0  b1 ... bi-1 bi ... bj bj+1 ... bk-1 bk
    	          ____________________________________________
    	          
    	Filho A:   -1  -1 ...  -1  ai ... aj  -1  ...  -1  -1
    	Filho B:   -1  -1 ...  -1  bi ... bj  -1  ...  -1  -1
     
    */
    for(k = 0; k < G->V; k++){
        filho[0][k] = (k >= i && k <= j ? Pai[0]->harbor[k] : -1);
        filho[1][k] = (k >= i && k <= j ? Pai[1]->harbor[k] : -1);
    }
    
    //índice apontador para o fim da rota: construção de trás para frente
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
    //caso haja posições nos filhos a serem preenchidas
    for(int m=0; m<2; m++){
        if(n[m] < 0)
        		continue;
        	//se o índice ainda nao chegou no início, há valores faltando
         for(k = G->V-1; k >= 0; k--){
             if(filho[m][k] == -1){
             	//para cada região não preenchida, encontra algum valor faltante e insere
                 for(int l = G->V-1; l >= 0; l--)
                     if(!isIn(Pai[1-m]->harbor[l], filho[m])){
                         filho[m][k] = Pai[1-m]->harbor[l];
                         break;
                     }
             }
         }
    }
    return;        
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

void prova_do_ortiz(Solution **jovem, int nFighter, Solution **pai){//famoso torneio
	/*
		recebe uma população de jovens;
		sorteia uma quantidade n de candidatos;
		seleciona 2 pais dentre os candidatos para cruzamento;
	*/
    bool marcado[PSIZE];
    int candidato[nFighter];
    int idx, idx1, idx2, d1, d2;
    int nCandidatos = nFighter, i, j;
    //inicialização
    for(i = 0; i < PSIZE; i++)
        marcado[i] = false;
    //seleção de candidato
    for(i = 0; i < nFighter; i++){
        do
            idx = rand() % PSIZE;
        while(marcado[idx]);
        candidato[i] = idx;
        marcado[idx] = true;
    }
    //reset
    for(i = 0; i < PSIZE; i++)
        marcado[i] = false;
    //torneio até sobrarem apenas 2 universitários com chance de não reprovar
    while(nFighter > 2){
        //seleção aleatória do primeiro candidato
        idx1 = rand() % (nCandidatos);
        //avanço até encontrar o primeiro disponível (não marcado)
        if(marcado[candidato[idx1]])
            while(marcado[candidato[ (++idx1) % nCandidatos ]]);
        idx1 = idx1 % nCandidatos;
        //repetição dos passos acima para encontrar o segundo candidato
        do{
            idx2 = rand() % (nCandidatos);
            if(marcado[candidato[idx2]])
                while(marcado[candidato[ (++idx2)%nCandidatos ]]);
            idx2 = idx2 % nCandidatos;
        }while(idx1 == idx2);
        //eliminação do perdedor (pior candidato)
        d1 = jovem[ candidato[idx1] ]->distance;
        d2 = jovem[ candidato[idx2] ]->distance;
        idx = (d1 > d2 ? idx1 : idx2);
        marcado[candidato[idx]] = true;
        nFighter--;
    }
    //atribuição dos selecionados ao vetor de pais
    for(i = 0, idx = 0; i < nCandidatos && idx < 2; i++)
        if(!marcado[ candidato[i] ]){
            for(j = 0; j < G->V; j++)
                pai[idx]->harbor[j] = jovem[ candidato[i] ]->harbor[j];
            pai[idx]->distance = jovem[ candidato[i] ]->distance;
            idx++;
        }
}

void print_populacao(Solution** p, int size_p){
    for(int i = 0; i < size_p ; i++){
        printf("[Cromossomo %2d - Distância %4d]  : ", i, p[i]->distance);
        print_arr(p[i]->harbor);
    }
}

Solution* construcao(){
    /*
    	Constrói uma solução viável para o conjunto de soluções relativo à 
    	população inicial do algoritmo genético. A construção é realizada
   	com o método guloso com características aleatórias, através do qual
    	dois portos são sorteados aleatoriamente para comporem o início da solução, 
    	e em seguida, aplica-se o método guloso para obter uma solução viável.
    */
    Solution *novaSolucao = new_solution(G->V);
    int port1, port2, weight = G->V-1, menor, position, k;
    int *demand = malloc(sizeof(int) * G->V);
    
    //Cópia do array de demandas, permitindo que seja possível alterá-lo sem preocupações
    for(int i=0; i<G->V; i++)
        demand[i] = DEMAND[i];
    
    //Seleção dos dois portos iniciais. Port1 e port2 devem ser diferentes entre si
    //e devem ser viáveis.
    do
        port1 = rand() % (G->V-1) + 1;
    while(DRAFT[port1] < (weight));
    weight--;

    do
        port2 = rand() % (G->V-1) + 1;
    while(DRAFT[port2] < weight || port2 == port1);
    weight--;

    //Inicialização da solução, atualização da distância percorrida
    novaSolucao->harbor[0] = port1;
    novaSolucao->harbor[1] = port2;
    novaSolucao->distance = G->adj[0][port1] + G->adj[port1][port2];
    demand[port1] = demand[port2] = 0;
    k = port2;

    //Etapa gulosa da construção.
    for(int i = 2; i < G->V; i++){
        menor = __INT16_MAX__;
        for(int j = 0; j < G->V; j++)
            if(G->adj[k][j] < menor && demand[j] == 1 && weight <= DRAFT[j]){
                position = j;
                menor = G->adj[k][j];
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
    /*
		 Responsável por realizar mutação nos cromossomos, permitindo maior
		 variabilidade no espaço de busca. Sorteiam-se dois índices e trocam-se os dois
		 de locus. Verifica-se a viabilidade da solução, e caso positivo, aceita-se a troca.
	 	 Caso contrário, retorna-se ao estado original e refaz-se o sorteio
    */
    int ind1, ind2, auxTroca;

	//loop enquanto não gerar uma solução válida
    while(1){
        //sorteio
        ind1 = rand() % (G->V - 1);
        do
          ind2 = rand() % (G->V - 1);
        while(ind1 == ind2);
	     //troca dos portos sorteados
        auxTroca = P->harbor[ind1];
        P->harbor[ind1] = P->harbor[ind2];
        P->harbor[ind2] = auxTroca;
		  //condição de saída
        if(is_Solution(P->harbor))
            break;
        //Caso a troca gere uma solução não viável, desfaz-se a troca.
        auxTroca = P->harbor[ind1];
        P->harbor[ind1] = P->harbor[ind2];
        P->harbor[ind2] = auxTroca;
    }
    P->distance = fitness(P->harbor);
}

void copiar(Solution *S, int *solucao){
    /*
    	Função utilizada para efetuar a cópia do conteúdo 
    	dos ponteiros de Solution
    */
    for(int i = 0; i < G->V; i++)
        S->harbor[i] = solucao[i];
    S->distance = fitness(solucao);
}

void selectSubstitute(Solution **populacao, Solution **filhos){
    //posições 0 e 1 para os melhores filhos
    //posições 2 e 3 para os piores pais
    int index[4], i;

    index[0] = index[2] = 0; 
    index[1] = index[3] = 1;

    /*Verifica-se quais os melhores filhos obtidos no cruzamento atual
     *e armazena-os nas posições 0 e 1 do array index*/
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
    /*Verifica-se quais os piores elementos da geração atual
     *e armazena-os nas posições 2 e 3 do array index*/
    if(populacao[ index[2] ]->distance < populacao[ index[3] ] ->distance){
        i = index[2];
        index[2] = index[3];
        index[3] = i;
    }
    for(i = 2; i < PSIZE; i++){
        if(populacao[i]->distance > populacao[ index[2] ]->distance){
            index[3] = index[2];
            index[2] = i;
        } else if(populacao[i]->distance > populacao[ index[3] ]->distance ){
            index[3] = i;
        }
    }
    /*Verifica as possibilidades de inserção dos filhos na geração dos pais
     *de acordo com a ordem do valor da fitness entre os avaliados. */
    if(filhos[ index[0] ]->distance < populacao[ index[3] ]->distance){
        copiar(populacao[ index[3] ], filhos[ index[0] ]->harbor);
        if(filhos[ index[1] ]->distance < populacao[ index[2] ]->distance)//O segundo melhor filho
            copiar(populacao[ index[2] ], filhos[ index[1] ]->harbor);
    }else if(filhos[ index[1] ]->distance < populacao[ index[2] ]->distance){
        copiar(populacao[ index[2] ], filhos[ index[1] ]->harbor);
    }else if(filhos[ index[0] ]->distance < populacao[ index[2] ]->distance){
        copiar(populacao[ index[2] ], filhos[ index[0] ]->harbor);
    }else if(filhos[ index[1] ]->distance < populacao[ index[3] ]->distance){
        copiar(populacao[ index[3] ], filhos[ index[1] ]->harbor);
    }
}

void selectMergesort(Solution **populacao, Solution **filhos, int numFilhos){
    int ind = PSIZE - 1;

    mergeSort(populacao, 0, PSIZE - 1);
    mergeSort(filhos, 0, numFilhos - 1);

    for(int i = PSIZE - 1; i >= 0; i--)
        if(filhos[i]->distance < populacao[ind]->distance)
            if((rand() % 100) < 30)
                copiar(populacao[ind--], filhos[i]->harbor);   
    
    shuffle(populacao, PSIZE);
}

void array_swap(Solution** arr, int i, int j)
{
	Solution* aux = arr[i];
	arr[i] = arr[j];
	arr[j] = aux;
}

void exame(Solution** pais, Solution** filhos, int cota)
{
	float avg = 0;
	int i, j, x, aux;
	int size = PSIZE * 2;
	int nbytes = PSIZE*sizeof(Solution*);
	Solution* p[size];
	
	//preenchimento
	memcpy(p, filhos, nbytes);
	memcpy(p+PSIZE, pais, nbytes);
	
	//calculo da media
	for(i=0; i<size; i++)
		avg += p[i]->distance;
	avg /= (float) size;
	
	//seleção do pivô [x]
	x=0;
	for(i=1; i<size; i++)
		if(abs(p[i]->distance - avg) < abs(p[x]->distance - avg))
			x = i;
	
	//aplicação da cota
	for(i=0; i < cota; i++)
	{
		aux = i;
		//seleção do melhor
		for(j = i+1; j < size; j++)
			if(p[j]->distance < p[aux]->distance)
				aux = j;
		array_swap(p, i, aux);
	}
	
	//pivoteamento
	array_swap(p, i, x);
	x=i++;
	j=size-1;
	while(i < j)
	{
		while(p[i]->distance < p[x]->distance && i < j)
			i++;
		while(p[j]->distance > p[x]->distance && j > i)
			j--;
		if(i==j)
			break;
		array_swap(p,i,j);
		i++;
	}
	if(p[i]->distance > p[x]->distance)
		i--;
	array_swap(p, i, x);
	
	//restituição da população de entrada
	memcpy(pais, p, nbytes);
	memcpy(filhos, p+PSIZE, nbytes);
}

int AlgGenetico(){
    //variáveis
    int i, j, numFilhos, nFighter, taxaMutacao;
    int gerAtual = 0, ultMelhor = 0;
    int *bebes[2];
    
    Solution **populacao = malloc(sizeof(Solution*) * PSIZE);
    Solution **filhos = malloc(sizeof(Solution*) * PSIZE);
    Solution **pais = malloc(sizeof(Solution*) * 2);
    Solution *melhor = new_solution();
    
    //inicialização
    pais[0] = new_solution();
    pais[1] = new_solution();
    bebes[0] = malloc(sizeof(int) * DIM);
    bebes[1] = malloc(sizeof(int) * DIM);
    
    for(i = 0; i < PSIZE; i++){
        populacao[i] = construcao();
        if(!i || populacao[i]->distance < melhor->distance)
            copiar(melhor, populacao[i]->harbor);
        filhos[i] = new_solution();
    }

    //execução da busca local
    for(i = 0; i < PSIZE; i++){
        copiar(populacao[i], populacao[i]->harbor);
        fixed_swap(populacao[i]);
    }

    //início do ciclo
    while((gerAtual - ultMelhor) < MAX_ITER){
                
        //Cruzamento
        numFilhos = 0;
        for(i = 0; i < DIM; i++)
            if((rand() % 100) < 95){ // Probabilidade de ocorrer 'crossover'
                //Seleção dos pais
                nFighter = rand() % (DIM/2) + (DIM/4);
                prova_do_ortiz(populacao, nFighter, pais);
                //geração de soluções válidas
                do
                    order1Crossover(pais, bebes);
                while(!is_Solution(bebes[0]) || !is_Solution(bebes[1]));
                //atribuição ao vetor de filhos
                copiar(filhos[numFilhos++], bebes[0]);
                copiar(filhos[numFilhos++], bebes[1]);
            }
        
        //Preenchimento complementar do conjunto de filhos
        while(numFilhos < PSIZE){
            i = rand() % PSIZE;
            copiar(filhos[numFilhos++], populacao[i]->harbor);
        }

        //Seleção dos agentes da população original e dos filhos a serem substituídos
        //selectSubstitute(populacao, filhos);
        //selectMergesort(populacao, filhos, numFilhos);
        exame(populacao, filhos, (int) PSIZE * 0.1);
        
        //Mutação
        taxaMutacao = ((gerAtual - ultMelhor) * PERC_MUT) / MAX_ITER;
        for(i = 0; i < PSIZE; i++)
            if((rand() % 100) < PERC_MUT)
                mutacao(populacao[i]);
        
        //execução da busca local
        if((gerAtual - ultMelhor) % 3 == 0)
            for(i = 0; i < PSIZE; i++){
                copiar(populacao[i], populacao[i]->harbor);
                fixed_swap(populacao[i]);
            }

        //Atualização da melhor solução corrente e avanço de geração
        for(i = 0; i < PSIZE; i++)
            if(populacao[i]->distance < melhor->distance){
                copiar(melhor, populacao[i]->harbor);
                ultMelhor = gerAtual;
            }
        gerAtual++;
    }
    
    ultMelhor = melhor->distance;
    CONT_GER = gerAtual;
    
    //Liberação de memória alocada
    for(i = 0; i < PSIZE; i++){
        free_solution(populacao[i]);
        free_solution(filhos[i]);
    }
    
    free_solution(melhor);
    free_solution(pais[0]);
    free_solution(pais[1]);
    free(pais);
    free(bebes[0]);
    free(bebes[1]);
    free(populacao);
    free(filhos);

    return ultMelhor;
}

void fixed_swap(Solution* s){
    int aux, copy[DIM];
    
   //gera uma cópia da solução de entrada
   for(int k=0; k < G->V; k++)
       copy[k] = s->harbor[k];
    
    for(int i=0; i < DIM-1; i++)
        for(int j=i+1; j < DIM; j++){
            //gera uma solução
            aux = copy[j];
            copy[j] = copy[i];
            copy[i] = aux;

            if(is_Solution(copy) && fitness(copy) < s->distance)
                copiar(s, copy);
            
            //desfaz o movimento
            aux = copy[j];
            copy[j] = copy[i];
            copy[i] = aux;
        }
}

void merge(Solution **Arr, int start, int middle, int end){
    Solution *temp[(end - start + 1)]; 
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
}

void mergeSort(Solution **Arr, int start, int end){
    if (start < end){
        int middle = start + (end - start)/2;
        mergeSort(Arr, start, middle);
        mergeSort(Arr, middle+1, end);
        merge(Arr, start, middle, end);
    }
}

void shuffle(Solution **Arr, int nTrocas){
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