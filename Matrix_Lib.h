#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <stdarg.h>

#ifndef __MATRIX_H__
#define __MATRIX_H__

extern double toGrade(double radiant);
extern void print_vector(int len, double* v);
extern void copy_vector(int len, double* v, double* v_cp);
extern void copy_matrix(int row, int col, double (*matrix_1)[col], double (*matrix_2)[col]);
extern void printMatrix(int row, int col, double (*matrix)[col]);
extern void printMatrix_pointer(int row, int col, double** matrix);
extern void transposed_matrix(int row, int col, double (*matrix)[col], double (*matrix_T)[col]);
extern int* getDegree(int row, int col);
extern double sarrus_rule(int row, int col, double (*matrix)[col]);
extern double getDeterminant(int row, int col, double (*matrix)[col], bool sarrus);
extern void reduceMatrix(int row, int col, int delete_row, int delete_col, double (*matrix)[col], double (*new_matrix)[col-1]);
extern void sum_vector(int len, double* v1, double* v2, double* v_sum);
extern void prod_s_vector(int len, double* v, double scale, double* v_prod);
extern double norm(double* v, int len);
extern double angle(double* v1, double* v2, int len, char* g_fun);
extern double s_product(double* v1, double* v2, int len);
extern double* v_product(double* v1, double* v2, int len);
extern void prod_s(int row, int col, double s, double (*matrix)[col], double (*m_p)[col]);
extern void sum_matrix(int row, int col, double (*matrix1)[col], double (*matrix2)[col], double (*m_s)[col]);
extern void matrix_prod(int row_a, int col_a, int row_b, int col_b, double (*m_a)[col_a], double (*m_b)[col_b], double (*m_p)[col_b]);
extern bool isMatrixScale(int row, int col, double (*matrix)[col]);
extern bool isMatrixReverse(int row, int col, double (*matrix_a)[col], double (*matrix_b)[col]);
extern void getReverseMatrix(int row, int col, double (*matrix)[col], double (*reverse_matrix)[col]);
extern int getRank(int row, int col, double (*matrix)[col], bool show_gauss_matrix);
extern bool isIndipendent(int len_vector, int args, ...);
extern void show(int row, int col, double (*matrix)[col], bool print_matrix, bool sarrus);

#endif
