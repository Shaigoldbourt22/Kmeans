#ifndef SHAI_OREN_PROJ_SPKMEANS_H
#define SHAI_OREN_PROJ_SPKMEANS_H
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <assert.h>


enum Goal {spk = 1 , wam = 2, ddg = 3 , lnorm = 4, jacobi = 5, kmeans = 6};
int calc_k(double* eig_val, int n);
int compare(double x , double y);
int calc_num_points(char* file_name);
int calc_point_size(char* file_name);
int* find_indexes(double** sym_mat, int n);
double norm_oc2(double* point1, double* point2, int d);
double compute_diff(double** mat, int n);
double** jacobi_calc(double** sym_mat, int n);
double** create_degree_mat(double** weighted_adjacency_mat, int n);
double** power_minus_half(double** mat, int n);
double** create_weighted_adjacency_mat(double** points, int n, int d);
double** create_norm_laplacian_mat(double** degree_mat, double** w_mat, int n);
double** read_file(char* file_name);
double** malloc_two_d_arr(int rows, int cols);
static void k_means(double** points, double** centroids, int k, int d, int n);
static double** loop_calc_centroids(double** centroids, double** points, int k, int d, int n);
static double** calc_centroids(double** centroids, double** points, int k, int d, int n);
static double norm_oc(double* point1, double* point2, int d);
static int smaller_then_eps(double** centroids, double** new_centroids, int size, int d);
int choose_algo(int goal ,double** points, double** centroid, double** T, int n ,int d,int k);
void multiply_mat(double** mat_1, double** mat_2, double** res_mat, int n);
void calc_p(double** sym_mat, double** p, int n);
void transpose(double** mat, double** trans_mat, int n);
void free_two_d_arr(double** matrix, int rows);
void print_one_d_vec(double* arr, int size);
void print_two_d_vec(double** arr, int rows, int cols);
void selectionSort(double* arr, int n);
#endif

