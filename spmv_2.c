/*

* This program executes spmv with COO and CSR format and a matrix with density 0.2
* Generates the matrix directly in the support structures for COO format by extracting a value different from zero with probability 0.8
* From the COO data structures we then create the CSR data structures
* All the data structures are dynamically allocated
* Tested with matrix from: 2x2 -- 200x200

*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define DEBUG 1

#define MAX_FLOAT_VALUE 10
#define ROWS 4
#define COLS 4
#define DENSITY 0.2
#define CODIFICA "COO"
#define ITERATIONS 6


/******************************************************************* DECLARING FUNCTIONS *********************************************************************/

void spmv_coo();
void spmv_csr();
void spmv_ell();

float *generate_x();
void print_info();
void print_y();
void print_debug();
void switch_float_array();
void switch_int_array();
void print_matrix();
float fun();
void print_ptr_indices_data_row();
void print_float(float f);
int count_max_nz_per_row(int *ptr);
void reset_y();

int number_of_nz = 0;
int ptr_iterator = 0;
int expected_nz = ROWS*COLS*DENSITY;
/******************************************************************** MAIN *************************************************************************************/

int main(int argc, char *argv[]){
	printf("\nInitializing all variables...\n");
	
	int rows = ROWS;
	int cols = COLS;
	int data_iterator = 0;
	int row_iterator = 0;
	int indices_iterator = 0;
	int max_nz_per_row = 0;
	int nz_per_row = 0;
	
	time_t t;
	srand((unsigned) time(&t));
	float extract;

	printf("\nAllocating data structures for spmv...\n");

	int *indices = (int*) malloc(expected_nz * sizeof(int));
	float *data = (float*) malloc(expected_nz * sizeof(float));
	int *row = (int*) malloc(expected_nz * sizeof(int));
	int *ptr = (int*) calloc((expected_nz) * sizeof(int));
	
	
	if (row == NULL) {
    	printf("Error: Memory not allocated.\n");
   		exit(0);
	}
	if (data == NULL) {
    	printf("Error: Memory not allocated.\n");
   		exit(0);
	}
	if (indices == NULL) {
    	printf("Error: Memory not allocated.\n");
   		exit(0);
	}





  //init structures
  for (int i = 0; i < number_of_nz; i++){
  	indices[i] = 0;
  	data[i] = 0;
  	row[i] = 0;
  	ptr[i] = 0;
  }

	ptr[0] = 0;
	ptr_iterator = 1;

	printf("\nStarting to populate data structures...\n");

	//generate data, row and indices...
	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){

			extract = fun();

			if (extract != 0){
				number_of_nz++;
				if(DEBUG == 1){
					printf("value extracted: ");
					print_float(extract);
					printf("\n");
				}
				data[data_iterator] = extract;
				row[row_iterator] = i;
				indices[indices_iterator] = j;
				data_iterator++;
				row_iterator++;
				indices_iterator++;
				if (number_of_nz == expected_nz){
					printf("Matrix with right density generated!\n");
					break;
				}
			}
		}
		if (number_of_nz == expected_nz){
			break;
		}
	}

	//generating ptr
	for (int i = 0; i < number_of_nz; i++){
		ptr[row[i] + 1]++;
	}
	for (int i = 0; i < rows; i++){
		ptr[i + 1] += ptr[i];
  }


	//generate dense array of float x
	float *x = generate_x(cols);

	//generate empty array of float y
	float *y = (float*) malloc(cols * sizeof(float));


	if(strcmp(CODIFICA,"COO")==0){
		if (DEBUG){
			print_ptr_indices_data_row(ptr, indices, data, row);
		}
		printf("\nStarting multiplication with COO format...\n");

		asm volatile("nop\nnop\nnop\nnop\n");
		for (int i=0; i < ITERATIONS; i++){
			printf("Iteration: %d\n", i+1);
			//setting y to 0
			reset_y(y);
			spmv_coo(row, indices, data, number_of_nz, x, y);
		}
		asm volatile("nop\nnop\nnop\nnop\n");
		free(row);
		free(indices);
		free(data);
	}

	else if (strcmp(CODIFICA,"CSR")==0){
		if (DEBUG){
			print_ptr_indices_data_row(ptr, indices, data, row);
		}
		printf("\nStarting multiplication with CSR format...\n");
		for (int i = 0; i < ITERATIONS; i++){
			asm volatile("nop\nnop\nnop\n");
			reset_y(y);
			spmv_csr(ptr, indices, data, number_of_nz, x, y);
		}
	}

	else if (strcmp(CODIFICA,"ELL")==0){
		int max_row_elements = count_max_nz_per_row(ptr);

		int *data_ell = (int *) malloc(max_row_elements*sizeof(float));
		int *indices_ell = (int *) malloc(max_row_elements*sizeof(int));

		printf("\nStarting multiplication with ELL format...\n");

		//spmv_ell(indices_EL2, data_EL2, max_nz_per_row, rows, x, y);
	}

	if (DEBUG == 1){
		print_y(y, cols);
	}

	free(x);
	free(y);

	printf("\nProgram successfully terminated.\n\n");
	return 0;
}

/****************************************************************** FUNCTIONS ********************************************************************************/

void spmv_ell(int *indices, float *data, int max_nz_per_row, int rows, float *x, float *y) {
  asm volatile("starting_computation_ell:\n");
  int i, j;
  for(j = 0; j <  max_nz_per_row; j++){
    for(i = 0; i < rows; i++){
    	//printf("%d\n", max_nz_per_row);
      y[i] += data[j * rows + i] * x[indices[j * rows + i]];
      //printf("y[%d] = data[%d * %d + %d] * x[indices[%d * %d + %d] = %f\n",i,j,max_nz_per_row,i,j,max_nz_per_row,i, y[i]);

    }
  }
  asm volatile("ending_computation_ell:\n");
}

void spmv_coo(int *row, int *indices, float *data, int number_of_nz, float *x, float *y) {
	asm volatile("nop\nnop\nnop\n");
	asm volatile("starting_computation_coo:\n");
	//asm volatile("nop\nnop\nnop\n");
	int i;
	int index_x;
	int index_y;

	float y_index_y;
	float data_i;
	float x_index_x;
	float product_data_x;
	float final_sum;


	for(i = 0; i < number_of_nz ; i++){
		asm volatile("saving_index_x:\n");
		index_x = indices[i];

		asm volatile("saving_index_y:\n");
		index_y = row[i];
		
		asm volatile("Saving_y_index_y_in_variable:\n");
		y_index_y = y[index_y];

		asm volatile("Saving_array_data:\n");
		data_i = data[i];

		asm volatile("Saving_array_x:\n");
		x_index_x = x[index_x];

		asm volatile("compute_product:\n");
		product_data_x = data_i * x_index_x;

		asm volatile("compute_sum:\n");
		final_sum = y_index_y + product_data_x;

		asm volatile("Final_assignment:\n");
		y[index_y] = final_sum;
		//printf("y[%d] = %f\n", i, y[i]);
		asm volatile("incrementing_i_loop:\n");

	}
  asm volatile("ending_computation_coo:\n");
  //asm volatile("nop\nnop\nnop\nnop\n");
}


void spmv_csr(int *ptr, int *indices, float *data, int number_of_nz, float *x, float *y) {
	asm volatile("nop\nnop\nnop\n");
	asm volatile("starting_computation_csr:\n");
	int i, j;
  	float temp;
  	float data_j;
  	int indices_j;
  	float x_indices_j;
  	float mul;

  	asm volatile("starting_outer_loop:\n");
  	for(i = 0; i < number_of_nz ; i++){
  		asm volatile("save_index_y_in_temp:\n");
    	temp = y[i];
    	asm volatile("inner_loop:\n");
    	for(j = ptr[i]; j < ptr[i+1]; j++){
    		asm volatile("saving_data_j:\n");
			data_j = data[j];
			asm volatile("saving_indices_j:\n");
			indices_j = indices[j];
			asm volatile("saving_x_indices_j:\n");
			x_indices_j = x[indices_j];
			asm volatile("data_j_mul_x_indices_j:\n");
			mul = (data_j * x_indices_j);
			asm volatile("temp_plus_mul:\n");
			temp = temp + mul;
			asm volatile("incrementing_inner_loop:\n");
		}
    	asm volatile("saving_new_temp:\n");
    	y[i] = temp;
    	asm volatile("incrementing_outer_loop:\n");
  	}
  	asm volatile("ending_computation_csr:\n");
  	asm volatile("nop\nnop\nnop\nnop\n");
  
}


/****************************************************************** SUPPORT FUNCTIONS ********************************************************************************/


float *generate_x(int C){
	printf("\nGenerating random dense array x...\n");
	srand(time(NULL));
	float r;
	int i;
	int x_int;
	float x_float;
	float x_cut;
	float small_float;
	int small_int;


	float *array = (float*) malloc(C * sizeof(float));
	for (i=0; i<C; i++){
		array[i] = (float)rand()/(float)(RAND_MAX/MAX_FLOAT_VALUE);
		x_int = (int) array[i];
		x_float = array[i];
		x_cut = (float) x_int;
		small_float = x_float - x_int;
		small_float = small_float * 100;
		small_int = (int) small_float;
		if (DEBUG == 1){
			printf("generated... x[%d] = %d.%d\n",i, x_int, small_int);
		}
	}
	return array;
}


void print_y(float *y, int length){
	float y_float;
	float small_float;
	float y_cut;
	int small_int;
	int y_int;
	printf("\nARRAY y: \n");
	for (int i = 0; i < length; i++){
		y_int = (int) y[i];
		y_float = y[i];
		y_cut = (float) y_int;
		small_float = y_float - y_int;
		small_float = small_float * 100;
		small_int = (int) small_float;
		printf("y[%d] = %d.%d\n", i, y_int, small_int);
	}
}

void print_float(float f){
	float y_float;
	float small_float;
	float y_cut;
	int small_int;
	int y_int;

	y_int = (int) f;
	y_float = f;
	y_cut = (float) y_int;
	small_float = y_float - y_int;
	small_float = small_float * 100;
	small_int = (int) small_float;
	printf("%d.%d", y_int, small_int);
}


float fun(){
	int n = rand() % 5; 
	float val;
	if (n != 0){
		return 0;
	}
	else{
		val = (float)rand()/(float)(RAND_MAX/MAX_FLOAT_VALUE);
		return val;
	}
}

void print_ptr_indices_data_row(int *ptr, int *indices, float *data, int *row){
	printf("\nSupport structures:\n");
	printf("row = [ ");
	for (int i = 0; i < number_of_nz; i++){
		printf("%d ", row[i]);
	}
	printf("]\n");
	printf("ptr = [ ");
	for (int i = 0; i < ROWS+1; i++){
		printf("%d ", ptr[i]);
	}
	printf("]\n");
	printf("indices = [ ");
	for (int i = 0; i < number_of_nz; i++){
		printf("%d ", indices[i]);
	}
	printf("]\n");
	printf("data = [ ");
	for (int i = 0; i < number_of_nz; i++){
		print_float(data[i]);
		printf(" ");
	}
	printf("]\n");
}

int count_max_nz_per_row(int *ptr){
	int val_1 = 0; 
	int val_2 = 0;
	int j;
	int temp = 0;
	int max = 0;
	//ptr dim is ROWS + 1
	for (int i = 0; i < ROWS; i++){
		j = i + 1;
		val_1 = ptr[i];
		val_2 = ptr[j];
		temp = val_2 - val_1;
		if (temp > max){
			max = temp;
		}
	}
	return max;
}

void reset_y(float *y){
	for (int i=0; i < ROWS; i++){
		y[i] = 0;
	}
}