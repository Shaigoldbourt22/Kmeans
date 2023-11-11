#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "spkmeans.h"
#include <unistd.h>

int main(int argc , char* argv[]) {
    int goal;
    char* file_name;
    double** points;
    int d, n;
    if(argc != 3){
        printf("Invalid Input!");
        exit(1);
    }
    if(strcmp(argv[1], "wam")==0){
        goal=2;
    }
    if(strcmp(argv[1], "ddg")==0){
        goal=3;
    }
    if(strcmp(argv[1], "lnorm")==0){
        goal=4;
    }
    if(strcmp(argv[1], "spk")==0){
        goal=1;
    }
    if(strcmp(argv[1], "jacobi")==0){
        goal=5;
    }
    file_name = argv[2];
    points = read_file(file_name);
    n = calc_num_points(file_name);
    d = calc_point_size(file_name);

    choose_algo(goal, points, NULL, NULL, n, d, 0);
    return 0;
}


int calc_num_points(char* file_name){
    FILE* file;
    char ch;
    int total_vec = 0;

    file = fopen(file_name, "r");
    for(ch = fgetc(file) ; ch != EOF ; ch = fgetc(file)) {
        if (ch == '\n') {
            total_vec += 1;
        }
    }

    fclose(file);
    return total_vec;
}


int calc_point_size(char* file_name){
    FILE *file;
    char c;
    int size_vec = 0;
    file = fopen(file_name, "r");
    for(c = fgetc(file) ; c != EOF ; c = fgetc(file)) {
        if (c == '\n')
            break;
        if(c == ',')
            size_vec++;
    }
    return size_vec + 1;
}


double** read_file(char* file_name){
    FILE *file;
    double a;
    int total_vec, size_vec, i, j;
    double** p;
    a = 0;

    total_vec= calc_num_points(file_name);
    size_vec= calc_point_size(file_name);

    p = malloc_two_d_arr(total_vec , size_vec);
    file = fopen(file_name, "r");
    for(i = 0 ; i < total_vec ; i++){
        for(j = 0 ; j < size_vec ; j++){
            if(fscanf(file, "%lf", &a) == 1){
                p[i][j] = (double)a;
                fgetc(file);
            }
        }
    }
    fclose(file);
    return p;
}

int compare(double x , double y){
    return (x > y) - (y > x);
}


int choose_algo(int goal ,double** points, double** centroid, double** T, int n ,int d,int k) {

    int i, j, z;

    double sum;
    double *arr_eig_val;
    double **matrix, **matrix1, **matrix2, **matrix3, **mat_T;

    if(k > n || k < 0){
        printf("Invalid Input!");
        exit(1);
    }

    if (goal == 1 && k > 0) {
        k_means(points, centroid, k, d, n);
        print_two_d_vec(centroid, k, d);
    }
    else {
        switch (goal) {
            case 1: {
                matrix = create_weighted_adjacency_mat(points, n, d);
                matrix1 = create_degree_mat(matrix, n);
                matrix2 = create_norm_laplacian_mat(matrix1, matrix, n);
                matrix3 = jacobi_calc(matrix2, n);


                /* matrix 3 : first line - eignvalues , under that - eignvectors as colomns */
                arr_eig_val = (double*) calloc(n, sizeof(double));
                if (arr_eig_val == NULL) {
                    printf("An Error Has Occurred");
                    exit(1);
                }

                for (i = 0; i < n; i++)
                    arr_eig_val[i] = matrix3[0][i];

                /* sort the matrix by the eignvalues */
                selectionSort(arr_eig_val, n);

                /* k -  calc_k - by the algorithm  - send sorted eignvalues */
                k = calc_k(arr_eig_val, n);
                if (k == 1){
                    printf("An Error Has Occurred");
                    exit(1);
                }

                mat_T = malloc_two_d_arr(n, k);
                for (i = 0; i < k; i++) {
                    for (j = 0; j < n; j++) {
                        if (matrix3[0][j] == arr_eig_val[i]) {
                            for (z = 0; z < n; z++)
                                mat_T[z][i] = matrix3[z + 1][j];
                            matrix3[0][j] = arr_eig_val[n - 1] + 1;
                            break;
                        }
                    }
                }
                /*at that point the eigen values in the first row in matrix3 are mixed*/
                /* find T  - take the matrix of the k first vector (sorted by the eignvalues) and normalize the rows!! */
                for (i = 0; i < n; i++) {
                    sum = 0;
                    for (j = 0; j < k; j++)
                        sum += mat_T[i][j] * mat_T[i][j];
                    sum = pow(sum, 0.5);
                    if(sum == 0)
                        continue;
                    for (j = 0; j < k; j++)
                        mat_T[i][j] = mat_T[i][j] / sum;
                }
                free(arr_eig_val);
                for (i = 0; i < n; i++)
                    for (j = 0; j < k; j++)
                        T[i][j] = mat_T[i][j];
                break;
            }
            case 2: {
                matrix = create_weighted_adjacency_mat(points, n, d);
                print_two_d_vec(matrix, n, n);
                free_two_d_arr(matrix, n);
                break;
            }
            case 3: {
                matrix = create_weighted_adjacency_mat(points, n, d);
                matrix1 = create_degree_mat(matrix, n);
                print_two_d_vec(matrix1, n, n);
                free_two_d_arr(matrix, n);
                free_two_d_arr(matrix1, n);
                break;
            }
            case 4: {
                matrix = create_weighted_adjacency_mat(points, n, d);
                matrix1 = create_degree_mat(matrix, n);
                matrix2 = create_norm_laplacian_mat(matrix1, matrix, n);
                print_two_d_vec(matrix2, n, n);
                free_two_d_arr(matrix, n);
                free_two_d_arr(matrix1, n);
                free_two_d_arr(matrix2, n);
                break;
            }
            case 5: {
                matrix = jacobi_calc(points, n);
                print_two_d_vec(matrix, n + 1, n);
                free(matrix);
                break;
            }
        }
    }

    return k;
}


double** create_weighted_adjacency_mat(double** points, int n, int d){
    int i, j;
    double w;
    double ** adj_mat= malloc_two_d_arr(n,n);

    for(i=0; i<n; i++){
        for(j = i + 1; j < n; j++){
            w = exp(-(norm_oc2(points[i],points[j], d))/2);
            adj_mat[i][j] = w;
            adj_mat[j][i] = w;
        }
    }
    for(i = 0; i < n; i++)
        adj_mat[i][i]=0;

    return adj_mat;
}


double norm_oc2(double* point1, double* point2, int d) {
    double sum;
    int i;

    sum = 0;
    for (i = 0; i < d; i++) {
        sum += pow((point2[i] - point1[i]), 2);
    }
    return pow(sum, 0.5);
}


double** create_degree_mat(double** weighted_adjacency_mat, int n) {
    int i, j;
    double sum;
    double **degree_mat = malloc_two_d_arr(n, n);

    for (i = 0; i < n; i++) {
        sum = 0;
        for (j = 0; j < n; j++)
            sum += weighted_adjacency_mat[i][j];
        degree_mat[i][i] = sum;
    }
    return degree_mat;
}


double** create_norm_laplacian_mat(double** degree_mat, double** w_mat, int n){
    double **power_d, **temp, **res;
    int i, j;
    power_d = power_minus_half(degree_mat, n);
    temp = malloc_two_d_arr(n,n);
    multiply_mat(power_d, w_mat, temp, n);
    res = malloc_two_d_arr(n,n);
    multiply_mat(temp, power_d, res,n);
    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            if(i == j)
                res[i][j]=1-res[i][j];
            else
                res[i][j]=0-res[i][j];
        }
    }
    free_two_d_arr(power_d, n);
    free_two_d_arr(temp, n);
    return res;
}


void multiply_mat(double** mat_1, double** mat_2, double** res_mat, int n) {
    int i, j, t;
    double sum;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            sum = 0;
            for (t = 0; t < n; t++)
                sum += mat_1[i][t] * mat_2[t][j];
            res_mat[i][j] = sum;
        }
    }
}


double** power_minus_half(double** mat, int n){
    int i;
    double** res_mat= malloc_two_d_arr(n,n);
    for(i = 0; i < n; i++){
        res_mat[i][i]= pow(mat[i][i], -0.5);
    }

    return res_mat;

}


/* return double** - (n+1) - first row - the n eignvalues and  in the next rows - the n eignvectors
 as columns corresponding to the order of the eignvalues
 sym_mat = norm_lap */
double** jacobi_calc(double** sym_mat, int n){
    int i, j, z;
    double **p, **temp, **A_prime, **temp_eig_vectors, **eig_vectors, **transpose_mat, **result_matrix;
    eig_vectors = malloc_two_d_arr(n,n);
    temp_eig_vectors = malloc_two_d_arr(n,n);
    p = malloc_two_d_arr(n,n);
    temp = malloc_two_d_arr(n,n);
    A_prime = malloc_two_d_arr(n,n);
    transpose_mat = malloc_two_d_arr(n,n);
    result_matrix = malloc_two_d_arr(n+1,n);

    for(i=0; i<n; i++)
        eig_vectors[i][i]=1;

    for(z = 0; z < 100; z++) {
        calc_p(sym_mat, p, n);
        multiply_mat(eig_vectors, p ,temp_eig_vectors, n);
        for(i=0; i<n; i++)
            for (j=0; j<n; j++)
                eig_vectors[i][j] = temp_eig_vectors[i][j];
        transpose(p,transpose_mat,n);
        multiply_mat(transpose_mat, sym_mat, temp, n);
        multiply_mat(temp, p ,A_prime, n);
        if ((compute_diff(sym_mat, n) - compute_diff(A_prime, n)) <= pow(10, -5))
            break;
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                sym_mat[i][j] = A_prime[i][j];
    }
    for(i=0; i<n; i++)
        result_matrix[0][i] = A_prime[i][i];

    for (i=1; i<n+1; i++)
        for (j=0; j<n; j++)
            result_matrix[i][j] = eig_vectors[i-1][j]; 

    free_two_d_arr(p, n);
    free_two_d_arr(A_prime, n);
    free_two_d_arr(temp, n);
    free_two_d_arr(transpose_mat, n);
    free_two_d_arr(eig_vectors, n);

    return result_matrix;
}


void calc_p(double** sym_mat, double** p, int n){
    double t, c, s, theta;
    int sign, i, j, row, col;
    for(i=0; i<n; i++)
        for(j=0; j<n; j++)
            if(i==j)
                p[i][j]=1;
            else
                p[i][j]=0;

    row = find_indexes(sym_mat, n)[0];
    col = find_indexes(sym_mat, n)[1];
    theta=(sym_mat[col][col]-sym_mat[row][row])/(2*sym_mat[row][col]);
    sign = (theta >= 0) ? 1 : -1;
    t = (theta == 0) ? 1: sign/(fabs(theta) + sqrt(pow(theta, 2) + 1));


    c = 1/(sqrt(pow(t,2)+1));
    s=t*c;
    p[row][row] = c;
    p[col][col] = c;
    p[row][col] = s;
    p[col][row] = -s;
}


int calc_k(double* eig_val, int n){
    int i, max_index = 0;
    double n_double, max = -1;
    n_double=(double) n;
    for(i=0; i< floor(n_double/2)-1; i++)
        if(max<fabs(eig_val[i]-eig_val[i+1])){
            max = fabs(eig_val[i]-eig_val[i+1]);
            max_index = i;
        }
    return max_index+1;
}


double compute_diff(double** mat, int n){
    int i, j;
    double sum = 0;
    for(i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            if (i != j)
                sum += pow(mat[i][j], 2);

    return sum;
}


/*sym_mat = A_mat*/
int* find_indexes(double** sym_mat, int n){
    double max_val = -1;
    int i, j, row=0 ,col=0;
    int* ret_arr;

    for(i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if(i != j && fabs(sym_mat[i][j]) > max_val){
                max_val = fabs(sym_mat[i][j]);
                row = i;
                col = j;
            }
        }
    }
    ret_arr = (int*)calloc(2,sizeof (int));
    ret_arr[0]=row;
    ret_arr[1]=col;
    return ret_arr;
}


void transpose(double** mat, double** trans_mat, int n) {
    int i, j;
    for (i = 0; i < n; i++)
        for (j = 0; j < n; j++)
            trans_mat[j][i]=mat[i][j];
}


double** malloc_two_d_arr(int rows, int cols) {
    double **arr;
    int i;

    arr = (double**)malloc(rows * sizeof(double*));

    if(arr == NULL){
        printf("An Error Has Occurred");
        exit(1);
    }

    for (i = 0; i < rows; i++) {
        arr[i] = (double *) calloc(cols, sizeof(double));
        if(arr[i] == NULL){
            printf("An Error Has Occurred");
            exit(1);
        }
    }

    return arr;
}


void free_two_d_arr(double** matrix, int rows) {
    int i;
    for (i=0; i<rows; i++)
        free(matrix[i]);
    free(matrix);
}


void print_one_d_vec(double* arr, int size){
    int i;
    for(i = 0; i<size; i++){
        if(arr[i] > -0.5/1000)
            printf("%.4f", fabs(arr[i]));
        else
            printf("%0.4f", arr[i]);
        if (i != size - 1)
            printf(",");
    }
}


void print_two_d_vec(double** arr, int rows, int cols){
    int i;
    for(i = 0; i<rows; i++){
        print_one_d_vec(arr[i], cols);
        printf("\n");
    }
}


/* k_means implimentation from HW2 */
static void k_means(double** points, double** centroids, int k, int d, int n){
    int i, j;
    double** good_centroids;
    good_centroids=loop_calc_centroids(centroids,points, k, d, n);
    for(i=0; i<k; i++) {
        for (j = 0; j < d; j++)
            centroids[i][j] = good_centroids[i][j];
        free(good_centroids[i]);
    }
    free(good_centroids);
}


static double** loop_calc_centroids(double** centroids, double** points, int k, int d, int n){
    int iter = 0;
    double** new_centroids, **calced_centroids;
    new_centroids = calc_centroids(centroids,points, k, d, n);
    while((smaller_then_eps(new_centroids,centroids,k,d) == 0) & (iter < 300)){
        iter+=1;

        calced_centroids = calc_centroids(new_centroids,points, k, d, n);
        centroids=new_centroids;
        new_centroids=calced_centroids;
    }

    return new_centroids;

}


static double** calc_centroids(double** centroids, double** points, int k, int d, int n){
    double** new_centroids;
    double* cent_num_points;
    int index_of_close_cent, i, j, t, z, p;
    double min_range, range;

    new_centroids = malloc_two_d_arr(k,d);
    cent_num_points = (double*) malloc(k * sizeof(double));
    for(i=0; i<n; i++) {
        index_of_close_cent = 0;
        min_range = norm_oc(centroids[0], points[i], d);
        for (j = 0; j < k; j++) {
            range = norm_oc(points[i], centroids[j], d);
            if (range < min_range) {
                min_range = range;
                index_of_close_cent = j;
            }
        }
        cent_num_points[index_of_close_cent] += 1;
        for (t = 0; t < d; t++) {
            new_centroids[index_of_close_cent][t] += points[i][t];
        }
    }
    for (z = 0; z < k; z++) {
        for (p = 0; p < d; p++) {
            new_centroids[z][p] = new_centroids[z][p] / cent_num_points[z];
        }
    }

    return new_centroids;
}


static double norm_oc(double* point1, double* point2, int d) {
    double sum;
    int i;

    sum = 0;
    for (i = 0; i < d; i++) {
        sum += pow((point2[i] - point1[i]), 2);
    }
    return pow(sum, 0.5);
}


static int smaller_then_eps(double** centroids, double** new_centroids, int size, int d){
    int i;
    for(i=0; i<size; i++){
        if(norm_oc(centroids[i],new_centroids[i], d) > 0.001){
            return 0;
        }
    }
    return 1;
}

void selectionSort(double* arr, int n){
        int i, j, min_ind;
        double temp, min;
        for (i = 0; i < n - 1; i++) {
            min = arr[i];
            min_ind = i;
            for (j = i + 1; j < n; j++)
                if (arr[j] < min){
                    min = arr[j];
                    min_ind = j;
                }
            temp = arr[i];
            arr[i] = min;
            arr[min_ind] = temp;
        }
}

