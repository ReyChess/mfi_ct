/*
 ============================================================================
 Name        : MaximalesDFS.c
 Author      : Reyder Cruz de la Osa
 Version     : 1.0
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#pragma GCC target("sse4.2")
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/timeb.h>
#include <nmmintrin.h>

#define MAX_ITEMS 15						// Máxima longitud permitida por defecto, para un itemset.
#define POW_DOS_MAX_ITEMS (1 << MAX_ITEMS)  // Dos elevado a la 15, para la cantidad de combinaciones.
#define MAX_UNIQUE_ITEMS 100000	    		// Máxima cantidad de ítems, por defecto.
#define MAX_TRANSACTIONS 1000000			// Máxima cantidad de transacciones, por defecto.
#define MAX_MAXIMAL 200000					// Máxima cantidad de maximales, por defecto.
#define BLOOM_FILTER_SIZE 1000000			// Tamaño del filtro BLOOM.



// Definición de la estructura Itemset
typedef struct {
	unsigned int items[MAX_ITEMS];
    unsigned int sup;
    char size;
} Itemset;

// Definición de la estructura Transaction
typedef struct
{
   unsigned int count;
   unsigned int* items;
} Transaction;

//Definición de la estructura MaximalList
typedef struct {
    Itemset sets[MAX_MAXIMAL];
    unsigned int count;
} MaximalList;

//Definición de la estructura Bloom
typedef struct {
    Itemset** sets;
    unsigned int sets_count;
    char flat;
} Bloom;

// Definición de la estructura Matriz Binaria
typedef struct
{
   unsigned int* bti_data;
   unsigned int** bti;
   unsigned int btiCountRow;
   unsigned int btiCountColumns;
   char btiItemsCountLastRow;
} Bti;

// Definición de la estructura UnoFrecuentes
typedef struct
{
	unsigned int* items;
	unsigned int items_count;
	unsigned int* items_support;
	unsigned int* map1_N;
} UnoFrecuentes;

// Definición de la estructura Block
typedef struct
{
   unsigned int value;
   unsigned int number;
} Block;



typedef struct {
    Block* blocks;
    unsigned int count;
} BlockList;



struct timeb T1, T2;
UnoFrecuentes* uno_frecuentes;
MaximalList* maximal_list[MAX_ITEMS + 1];
Bloom* bloom;
unsigned int* tempHash;
Transaction* transactions;
unsigned int transactions_count;



/* Activa los bits de la matriz compactada "b". Recorre las transacciones y para cada ítem,
 * si es frecuente (map1_N[item] > 0), activa el bit correspondiente en "b".
 */
void btiItem(const char* fileName, Bti* b, unsigned int* map1_N)
{
   for(unsigned int i = 0; i < transactions_count; i++)
   {
      for(unsigned int j = 0; j < transactions[i].count; j++)
      {
    	  if(map1_N[transactions[i].items[j]] > 0)
    	      b->bti[map1_N[transactions[i].items[j]] - 1][i / 32] |= 1U << (i % 32);
      }
   }
}

// Lee el conjunto de datos desde un fichero de nombre "fileName"
void read_transactions(const char* fileName, float min_support, Bti* b) {


	FILE* file = fopen(fileName, "r");
    if (!file) {
        perror("Error abriendo el fichero");
        exit(EXIT_FAILURE);
    }

    unsigned int unique_item_count = 0;
    unsigned int* unique_items = (unsigned int*) malloc(MAX_UNIQUE_ITEMS * sizeof(unsigned int));
    unsigned int* unique_supports = (unsigned int*) calloc(MAX_UNIQUE_ITEMS, sizeof(unsigned int));
    char* buffer = (char*) malloc(15000 * sizeof(char));
    char* ptr;
    char* eof = fgets(buffer,15000, file);
    unsigned int item, max_item = 0, max_support = 0;
    unsigned int* unique_check = (unsigned int*) calloc(MAX_UNIQUE_ITEMS, sizeof(unsigned int));     // bandera para saber si un ítem existe en la BD
    unsigned int** support_inv;													// matrix de soporte invertida para ordenar por soporte en o(n), n la cantidad de items
    unsigned int* support_inv_count;
    unsigned int* support_inv_size;

    /* se lee la BD
     * se almacenan los ítems únicos en unique_items y los soportes en unique_supports
     * se ignora el ítem con valor mayor a MAX_UNIQUE_ITEMS (retorno silencioso)
     * se almacena el valor máximo de los items en max_item y el valor máximo del soporte en max_support
     * se almacena la cantidad de transacciones en transaction_count.
     * se trunca la lectura de la BD si se intenta leer más de MAX_TRANSACTIONS transacciones (se emite un mensaje)
      * */
    transactions_count = 0;

    while (eof != NULL) {
    	transactions[transactions_count].items = (unsigned int*) malloc (MAX_UNIQUE_ITEMS * sizeof(unsigned int));
    	transactions[transactions_count].count = 0;
    	ptr = buffer;
    	while((*ptr) & 32)
        {
    	    item = 0;
    	    while((*ptr) & 16)
    	        item = item * 10 + (*ptr++) - 48;

    	    if((*ptr) & 32) ptr++;

    	    if(item < MAX_UNIQUE_ITEMS)
    	    {
    	        if(item > max_item)
    	        	max_item = item;

    	        if (unique_check[item] == 0)
    	        {
    	            unique_items[unique_item_count++] = item;
    	            unique_check[item] = 1;
    	        }

    	        unique_supports[item]++;
    	        if(unique_supports[item] > max_support)
    	           max_support = unique_supports[item];

    	        transactions[transactions_count].items[transactions[transactions_count].count++] = item;
    	    }
    	}

    	transactions[transactions_count].items = (unsigned int*) realloc (transactions[transactions_count].items, transactions[transactions_count].count * sizeof(unsigned int));
    	transactions_count++;
        if(transactions_count == MAX_TRANSACTIONS)
        {
            printf("Couldn't process all the transactions. Please modify the constant MAX_TRANSACTIONS.\n");
        	eof = NULL;
        }
        else eof = fgets(buffer,15000,file);
    }

    free(buffer);

    uno_frecuentes = (UnoFrecuentes*) malloc(sizeof(UnoFrecuentes));
    uno_frecuentes->items = (unsigned int*) malloc(unique_item_count * sizeof(unsigned int));
    uno_frecuentes->items_support = (unsigned int*) malloc(unique_item_count * sizeof(unsigned int));
    support_inv = (unsigned int**) malloc((max_support + 1) * sizeof(unsigned int*));
    support_inv_count = (unsigned int*) calloc(max_support + 1, sizeof(unsigned int));
    support_inv_size = (unsigned int*) calloc(max_support + 1, sizeof(unsigned int));
    for(int i = ceil(min_support * transactions_count); i <= max_support; i++)
    {
    	support_inv[i] = (unsigned int*) malloc(100 * sizeof(unsigned int));
    	support_inv_size[i] = 100;
    }
    uno_frecuentes->items_count = 0;

    // se almacena para cada valor de soporte, mayor o igual al umbral de soporte, los ítems que comparten ese soporte
    unsigned int sup_temp;
    for(int i = 0; i < unique_item_count; i++)
    {
    	sup_temp = unique_supports[unique_items[i]];
    	if(((float) sup_temp / transactions_count) >= min_support)
    	{
    		if(support_inv_count[sup_temp] == support_inv_size[sup_temp])
    		{
    			support_inv[sup_temp] = (unsigned int*) realloc(support_inv[sup_temp], (support_inv_size[sup_temp] + 100) * sizeof(unsigned int));
		        support_inv_size[sup_temp] += 100;

    		}
    		support_inv[sup_temp][support_inv_count[sup_temp]++] = unique_items[i];
    	}
    }

    // se almacenan los N uno frecuentes, ordenados por su soporte y se mapean, de 1 a N
    uno_frecuentes->map1_N = (unsigned int*) calloc((max_item + 1), sizeof(unsigned int));

    for(int i = ceil(min_support * transactions_count); i <= max_support; i++)
    {
    	if(support_inv_count[i] > 0)
    		for(int j = 0; j < support_inv_count[i]; j++)
    		{
    			uno_frecuentes->items[uno_frecuentes->items_count] = support_inv[i][j];
    			uno_frecuentes->map1_N[uno_frecuentes->items[uno_frecuentes->items_count]] = uno_frecuentes->items_count + 1;
    			uno_frecuentes->items_support[uno_frecuentes->items_count++] = i;
    		}
    }

    uno_frecuentes->items = (unsigned int*) realloc(uno_frecuentes->items, uno_frecuentes->items_count * sizeof(unsigned int));
    uno_frecuentes->items_support = (unsigned int*) realloc(uno_frecuentes->items_support, uno_frecuentes->items_count * sizeof(unsigned int));

    // se inicializa la matriz que almacena las transacciones a nivel de bits
    b->btiItemsCountLastRow = transactions_count % 32;
    b->btiCountRow = (b->btiItemsCountLastRow == 0) ? transactions_count/32 : transactions_count/32 + 1;

    b->bti_data = (unsigned int*) calloc(uno_frecuentes->items_count * b->btiCountRow, sizeof(unsigned int));
    b->bti = (unsigned int**) malloc(uno_frecuentes->items_count * sizeof(unsigned int*));
    for(int i = 0; i < uno_frecuentes->items_count; i++)
        b->bti[i] = &b->bti_data[i * b->btiCountRow];

    b->btiCountColumns = uno_frecuentes->items_count;

    btiItem(fileName,b, uno_frecuentes->map1_N);

    // se libera espacio
    free(unique_check);
    free(unique_items);
    free(unique_supports);
    /*for(int i = ceil(min_support * transactions_count); i <= max_support; i++)
        free(support_inv[i]);*/
    free(support_inv);
    free(support_inv_count);
    free(support_inv_size);
    fclose(file);
}

// Verifica si un conjunto es subconjunto de otro
int is_subset(unsigned int* subset, int subset_count, unsigned int* superset, int superset_count) {
	unsigned int low = 0;
	for (unsigned int i = 0; i < subset_count; i++)
	{
	    unsigned int high = superset_count - 1;
	    unsigned int found = 0;
	    while (low <= high)
	    {
	       unsigned int mid = (low + high) / 2;
	       if(uno_frecuentes->map1_N[superset[mid]] == uno_frecuentes->map1_N[subset[i]])
	       {
	          found = 1;
	          low = mid + 1;
	          break;
	       }
	       if(uno_frecuentes->map1_N[superset[mid]] < uno_frecuentes->map1_N[subset[i]])
	           low = mid + 1;
	       else
	           high = mid - 1;
	    }

	    if (!found) return 0;
	}
	return 1;
}


// Verifica si un conjunto es maximal
int is_maximal(unsigned int* itemset, int itemset_count) {
	unsigned long int h = 1;
	unsigned int seed = 19;
	unsigned int key;
	for (unsigned int i = 0; i < itemset_count; i++)
	    h = ((h * seed) + itemset[i]) % BLOOM_FILTER_SIZE;

	key = h;
	if(bloom[key].flat == 0)
		return 1;

	for(unsigned int j = 0; j < bloom[key].sets_count; j++)
		if (bloom[key].sets[j]->size > itemset_count && is_subset(itemset, itemset_count, bloom[key].sets[j]->items, bloom[key].sets[j]->size))
		    return 0; // No es maximal porque es subconjunto de otro conjunto maximal
	return 1;
}

// Función recursiva para encontrar patrones frecuentes maximales
void find_maximal_frequent(Bti* b, unsigned int* current_itemset, int current_count, BlockList* current_block, int start,
		float min_support, unsigned int current_support, unsigned int transaction_count, unsigned int* maximal_total) {

	BlockList new_block = {NULL, 0};
	unsigned int temp;
	unsigned int new_support = 0;
	unsigned int tail_pruning = 0;

	// Explorar extensiones del conjunto actual
    int has_frequent_superset = 0;

    for (int i = start; !tail_pruning && (i < uno_frecuentes->items_count) && (*maximal_total < MAX_MAXIMAL) && (current_count < MAX_ITEMS); i++)
    {
    	current_itemset[current_count] = uno_frecuentes->items[i];

        // calcular el soporte y crear el nuevo block
        new_support = 0;
        if((current_count + 1) == 1)
        {
        	new_support = uno_frecuentes->items_support[i];
        } else if((current_count + 1) == 2)
               {
        	       new_block.blocks = malloc(b->btiCountRow * sizeof(Block));
        	       unsigned int total_blocks = 0;

       	    	   for(unsigned int z = 0; z < b->btiCountRow; z++)
       	    	   {
        	           temp = b->bti[uno_frecuentes->map1_N[current_itemset[0]] - 1][z] & b->bti[uno_frecuentes->map1_N[current_itemset[1]] - 1][z];
        	    	   if(temp > 0)
        	    	   {
        	    		   new_block.blocks[total_blocks] = (Block){temp, z};
        	    		   total_blocks++;
        	    	   }
        	    	   new_support += _mm_popcnt_u32(temp);
        	       }

       	    	   new_block.count = total_blocks;
               }
               else
               {
            	   new_block.blocks = malloc(current_block->count * sizeof(Block));
            	   unsigned int valid_count = 0;

        	       for(int z = 0; z < current_block->count; z++ )
                   {
        	    	   temp = current_block->blocks[z].value & b->bti[uno_frecuentes->map1_N[current_itemset[current_count]] - 1][current_block->blocks[z].number];
        	           if(temp > 0)
        	           {
        	               new_support += _mm_popcnt_u32(temp);
        	               new_block.blocks[valid_count] = (Block){temp, current_block->blocks[z].number};
        	               valid_count++;
        	          }
                   }
        	       new_block.count = valid_count;
                }
        // calcular el soporte y crear el nuevo block

        if (new_support >= min_support)
        {
        	if(new_support == current_support)
         		tail_pruning  = 1;


            has_frequent_superset = 1; // Tiene un superconjunto frecuente
            find_maximal_frequent(b, current_itemset, current_count + 1, &new_block, i + 1, min_support, new_support, transaction_count, maximal_total);
            if(new_block.blocks != 0)
            	free(new_block.blocks);
        }
        else if(new_block.blocks != 0)
        	free(new_block.blocks);
    }

    // Si no tiene superconjuntos frecuentes, verificar si es maximal
    if (!has_frequent_superset && is_maximal(current_itemset, current_count)) {
        if(*maximal_total < MAX_MAXIMAL)
        {
        	memcpy(maximal_list[current_count]->sets[maximal_list[current_count]->count].items, current_itemset, current_count * sizeof(unsigned int));
        	maximal_list[current_count]->sets[maximal_list[current_count]->count].size = current_count;
        	maximal_list[current_count]->sets[maximal_list[current_count]->count].sup = current_support;
        	(*maximal_total)++;


            //con bloom
        	unsigned int combination[MAX_ITEMS - 1];
        	tempHash[0] = 1;
        	unsigned int tempHash_count = 1;
        	unsigned int it = 0;
        	unsigned long int h;
        	unsigned int key;
        	unsigned int seed = 19;
        	for(unsigned int r = 1; r < current_count; r++)
        	{
        		h = tempHash[it];
                key = (h * seed + current_itemset[r - 1]) % BLOOM_FILTER_SIZE;

                tempHash[tempHash_count++] = key;
                bloom[key].flat = 1;
                if(bloom[key].sets_count % 100 == 0)
                    bloom[key].sets = (Itemset**) realloc(bloom[key].sets, (bloom[key].sets_count + 100) * sizeof(Itemset*));
                bloom[key].sets[bloom[key].sets_count++] = &(maximal_list[current_count]->sets[maximal_list[current_count]->count]);

                // Inicializar la primera combinación
        		for (int i = 0; i < r; i++)
        		    combination[i] = i;

        	    while (1)
        	    {
        	        // Generar la siguiente combinación
        	        int i = r - 1;
        	        while (i >= 0 && combination[i] == current_count - r + i)
        	        {
        	            i--;
        	            it++;
        	        }


        	        if (i < 0)
        	            break; // Todas las combinaciones generadas

        	        combination[i]++;
        	        for (int j = i + 1; j < r; j++)
        	            combination[j] = combination[i] + j - i;

        	        h = tempHash[it];
        	        key = (h * seed + current_itemset[combination[r - 1]]) % BLOOM_FILTER_SIZE;

        	        tempHash[tempHash_count++] = key;
        	        bloom[key].flat = 1;
        	        if(bloom[key].sets_count % 100 == 0)
        	            bloom[key].sets = (Itemset**) realloc(bloom[key].sets, (bloom[key].sets_count + 100) * sizeof(Itemset*));
        	        bloom[key].sets[bloom[key].sets_count++] = &(maximal_list[current_count]->sets[maximal_list[current_count]->count]);
        	    }
        	}
            //con bloom

        	maximal_list[current_count]->count++;
        }
    }
}


// Función principal para encontrar patrones frecuentes maximales
void find_maximal_frequent_patterns(float min_support, Bti* b, unsigned int transaction_count) {
	unsigned int current_itemset[MAX_ITEMS];
    bloom = (Bloom*) calloc(BLOOM_FILTER_SIZE, sizeof(Bloom));

    for(int i = 1; i <= MAX_ITEMS; i++)
    {
    	maximal_list[i] = (MaximalList*) malloc(sizeof(MaximalList));
    	maximal_list[i]->count = 0;
    }

    unsigned int maximalTotal = 0;

    find_maximal_frequent(b, &(current_itemset[0]), 0, 0, 0, min_support, 0, transaction_count, &maximalTotal);

    printf("Total of Maximals %d\n", maximalTotal);
}

int main(int argc, char *argv[]) {
	ftime( &T1 );
	if (argc != 3) {
        fprintf(stderr, "Usage: %s <filename> <min_support>\n", argv[0]);
        return EXIT_FAILURE;
    }

    transactions = (Transaction*) malloc (MAX_TRANSACTIONS * sizeof(Transaction));



    Bti* b = (Bti*) malloc(sizeof(Bti));

    const char *filename = argv[1];
    float min_support = atof(argv[2]);

    if (min_support < 0 || min_support > 1) {
        fprintf(stderr, "Min support must be between 0 and 1\n");
        return EXIT_FAILURE;
    }

    FILE * ficheroSalida = 0;
    ficheroSalida = fopen("MFI_CTOut.txt","a");

    //Se usa en el llenado del filtro bloom
    tempHash = (unsigned int*) calloc(POW_DOS_MAX_ITEMS, sizeof(unsigned int));

    read_transactions(filename, min_support, b);
    min_support = min_support * transactions_count;



    // liberando espacio
    unsigned int id = 1;

    for(int i = 1; i <= MAX_ITEMS; i++)
    {
   	   for(int j = 0; j < maximal_list[i]->count; j++)
       {
           if(maximal_list[i]->sets[j].size > 0)
    	   {
    	    	fprintf(ficheroSalida, "%d -> ", id++);
    	    	for(int z = 0; z < maximal_list[i]->sets[j].size; z++)
    	    		fprintf(ficheroSalida, "%d ", maximal_list[i]->sets[j].items[z]);
    	    	fprintf(ficheroSalida, ": %d \n", maximal_list[i]->sets[j].sup);

    	   }
       }

    }


    free(tempHash);
    free(bloom);


    free(uno_frecuentes->items);
    free(uno_frecuentes->items_support);
    free(uno_frecuentes->map1_N);
    free(uno_frecuentes);

    free(b->bti_data);
    free(b->bti);
    free(b);

    /*for(int i = 0; i < transactions_count; i++)
       free(transactions[i].items);*/
    free(transactions);

    ftime( &T2 );
    int milisecs = ( ( T2.time - T1.time ) * 1000 ) + T2.millitm - T1.millitm;
    printf( "The search took ,\t %d sec, %d millisec\n", milisecs / 1000, milisecs % 1000 );
    fclose(ficheroSalida);


    return EXIT_SUCCESS;
}
