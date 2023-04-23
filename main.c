#include "Matrix_Lib.h"

double toGrade(double radiant){
    return radiant * (180/M_PI);
}

void print_vector(int len, double* v){
    printf("\n");
    for(int i=0;i<len;i++){
        printf(" %lf ", *(v+i));
    }
    printf("\n");
}

void copy_vector(int len, double* v, double* v_cp){
    for(int i=0;i<len;i++){
        *(v_cp+i) = *(v+i);
    }
}

void copy_matrix(int row, int col, double (*matrix_1)[col], double (*matrix_2)[col]){
    for(int i=0;i<row;i++){
        for(int y=0;y<col;y++){
            matrix_2[i][y] = matrix_1[i][y];
        }
    }
}

void printMatrix(int row, int col, double (*matrix)[col]){
    printf("\n");
    for(int i=0;i<row;i++){
        for(int j=0;j<col;j++){
           printf(" %4lf ",matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printMatrix_pointer(int row, int col, double** matrix){
    for(int i=0;i<row;i++){
        for(int j=0;j<col;j++){
           printf(" %4lf ",matrix[i][j]);
        }
        printf("\n");
    }
}

void transposed_matrix(int row, int col, double (*matrix)[col], double (*matrix_T)[col]){
    for(int i=0;i<row;i++){
        for(int y=0;y<col;y++){
            matrix_T[i][y] = matrix[y][i];
        }
    }
}

int* getDegree(int row, int col){
    int* degree = (int*) malloc(sizeof(int)*2);

    *(degree) = row;
    *(degree+1) = col;

    return degree;
}

int* searchZero(int row, int col, double (*matrix)[col]){
    int index_row = 0;
    int row_zero = 0;
    for(int i=0;i<row;i++){
        int row_zero_temp = 0;
        for(int j=0;j<col;j++){
            if(matrix[i][j]==0) row_zero_temp++;
        }
        if(row_zero_temp>row_zero){
            row_zero = row_zero_temp;
            index_row = i;
        }
    }

    int index_col = 0;
    int col_zero = 0;
    for(int j=0;j<col;j++){
        int col_zero_temp = 0;
        for(int i=0;i<row;i++){
            if(matrix[i][j]==0) col_zero_temp++;
        }
        if(col_zero_temp>col_zero){
            col_zero = col_zero_temp;
            index_col = j;
        }
    }

    int* arr_for_laplace = (int*) malloc(2 * sizeof(int));
    if(row_zero>=col_zero){
        *(arr_for_laplace) = index_row;
        *(arr_for_laplace+1) = 0;
    }else{
        *(arr_for_laplace) = index_col;
        *(arr_for_laplace+1) = 1;
    }

    //if(arr_for_laplace[1]==0) printf("\nRiga con più zeri : %d, numero di zeri : %d",arr_for_laplace[0], row_zero);
    //else printf("\nColonna con più zeri : %d, numero di zeri : %d",arr_for_laplace[0], col_zero);
    return arr_for_laplace;
}

double** resize_sarrus_matrix(int row, int col, double (*matrix)[col]){
    double** sarrus_matrix = (double**) malloc(row * sizeof(double*));
    for (int i=0;i<row;i++) sarrus_matrix[i] = (double*) malloc((col+2) * sizeof(double));

    for(int i=0;i<row;i++){
        for(int j=0;j<col+2;j++){
            if(j>=3) sarrus_matrix[i][j] = matrix[i][j-3];
            else sarrus_matrix[i][j] = matrix[i][j];
        }
    }
    return sarrus_matrix;
}

double sarrus_rule(int row, int col, double (*matrix)[col]){
    double** sarrus_matrix = resize_sarrus_matrix(row, col, matrix);

    double result = 0, value = 1;
    int y = 0;


    for(int i=0;i<row;i++){
        y = 0;
        value = 1;
        for(int j=0;j<col;j++){
            value *= sarrus_matrix[j][i+y];
            y++;
        }
        result += value;
    }

    for(int i=0;i<row;i++){
        y = 0;
        value = 1;
        for(int j=col-1;j>=0;j--){
            value *= sarrus_matrix[j][i+y];
            y++;
        }
        result -= value;
    }

    for(int i=0;i<row;i++) free(sarrus_matrix[i]);
    free(sarrus_matrix);

    return result;
}

void reduceMatrix(int row, int col, int delete_row, int delete_col, double (*matrix)[col], double (*new_matrix)[col-1]){
    int x = 0, y = 0;
    for(int i=0;i<row;i++){
        if(i != delete_row){
            for(int j=0;j<col;j++){
                if(j != delete_col){
                    new_matrix[x][y++] = matrix[i][j];
                }
            }
            x++;
            y = 0;
        }
    }
}

double getDeterminant(int row, int col, double (*matrix)[col], bool sarrus){
    if(row != col){
       printf("\nThe matrix is not square!");
       return 0;
    }
    else if(row == 1) return matrix[0][0];
    else if(row == 2) return ((matrix[0][0]*matrix[1][1])-(matrix[1][0]*matrix[0][1]));
    else if(row == 3 && sarrus) return sarrus_rule(row, col, matrix);
    else{
        double determinant = 0;
        int* major_zeros = searchZero(row, col, matrix);
        for(int i=0;i<row;i++){
            if(major_zeros[1]==1){
                double new_matrix[row-1][col-1];
                reduceMatrix(row, col, i, major_zeros[0], matrix, new_matrix);
                determinant += (matrix[i][major_zeros[0]] * pow(-1,((i+1)+(major_zeros[0]+1))) * getDeterminant(row-1, col-1, new_matrix, false));
            }else{
                double new_matrix[row-1][col-1];
                reduceMatrix(row, col, major_zeros[0], i, matrix, new_matrix);
                determinant += (matrix[major_zeros[0]][i] * pow(-1,((i+1)+(major_zeros[0]+1))) * getDeterminant(row-1, col-1, new_matrix, false));
            }
        }
        free(major_zeros);
        return determinant;
    }
    return 0;
}

void sum_vector(int len, double* v1, double* v2, double* v_sum){
    for(int i=0;i<len;i++){
        *(v_sum+i) = (*(v1+i) + *(v2+i));
    }
}

void prod_s_vector(int len, double* v, double scale, double* v_prod){
    for(int i=0;i<len;i++){
        *(v_prod+i) = (*(v+i)*scale);
    }
}

double norm(double* v, int len){
    double norm = 0;
    for(int i=0;i<len;i++){
        norm += pow(*(v+i), 2);
    }
    return sqrt(norm);
}

double angle(double* v1, double* v2, int len, char* g_fun){
    double product = 0, angle = 0;

    for(int i=0;i<len;i++){
        product += (*(v1+i)) * (*(v2+i));
    }
    double norm_p = norm(v1, len) * norm(v2, len);
    angle = (product/norm_p);

    if(strcmp(g_fun, "sin") == 0) return toGrade(asin(angle));
    else return toGrade(acos(angle));

    return 0;
}

double s_product(double* v1, double* v2, int len){
    double product = 0;

    for(int i=0;i<len;i++){
        product += (*(v1+i)) * (*(v2+i));
    }

    return product;
}

double* v_product(double* v1, double* v2, int len){
    if(len != 3) return NULL;
    double* v_prod = (double*) malloc(sizeof(double));
    double m_vectors[len][2];

    for(int i=0;i<len;i++){
        m_vectors[i][0] = *(v1+i);
        m_vectors[i][1] = *(v2+i);
    }

    int j = 0;
    for(int i=len-1;i>=0;i--){
        for(int y=i-1;y>=0;y--){
            double matrix_d[2][2];
            matrix_d[0][0] = m_vectors[i][0];
            matrix_d[0][1] = m_vectors[i][1];
            matrix_d[1][0] = m_vectors[y][0];
            matrix_d[1][1] = m_vectors[y][1];

            printMatrix(2, 2, matrix_d);
            *(v_prod+j) = (j == 0 || j == 2) ? -(getDeterminant(2, 2, matrix_d, false)) : getDeterminant(2, 2, matrix_d, false);
            j++;
        }
    }

    return v_prod;
}

void sum_matrix(int row, int col, double (*matrix1)[col], double (*matrix2)[col], double (*m_s)[col]){
    for(int i=0;i<row;i++){
        for(int y=0;y<col;y++){
            m_s[i][y] = (matrix1[i][y] + matrix2[i][y]);
        }
    }
}

void prod_s(int row, int col, double s, double (*matrix)[col], double (*m_p)[col]){
    for(int i=0;i<row;i++){
        for(int y=0;y<col;y++){
            m_p[i][y] = s*matrix[i][y];
        }
    }
}

void matrix_prod(int row_a, int col_a, int row_b, int col_b, double (*m_a)[col_a], double (*m_b)[col_b], double (*m_p)[col_b]){
    if(col_a != row_b){
        printf("\nShape error with %d != %d", col_a, row_b);
    }

    for(int i=0;i<row_a;i++){
        for(int y=0;y<col_b;y++){
            double result = 0;
            for(int j=0;j<col_a;j++){
                result += m_a[i][j]*m_b[j][y];
            }
            m_p[i][y] = result;
        }
    }
}

bool isMatrixScale(int row, int col, double (*matrix)[col]){
    for(int i=0;i<row-1;i++){
        for(int y=0;y<col;y++){
            if(matrix[i][y] != 0){
                for(int j=i+1;j<row;j++){
                    for(int z=y;z>=0;z--){
                        if(matrix[j][z] != 0) return false;
                    }
                }
                break;
            }
        }
    }
    return true;
}

bool isMatrixReverse(int row, int col, double (*matrix_a)[col], double (*matrix_b)[col]){
    if(row != col) return false;
    double matrix_p[row][col];
    matrix_prod(row, col, row, col, matrix_a, matrix_b, matrix_p);

    for(int i=0;i<row;i++){
        for(int y=0;y<col;y++){
            if(i == y && matrix_p[i][y] != 1) return false;
            else if( i != y && matrix_p[i][y] != 0) return false;
        }
    }
    return true;
}

void filter_reverse_matrix(int row, int col, double (*matrix)[col]){
    for(int i=0;i<row;i++){
        for(int y=0;y<col;y++){
            if(matrix[i][y] == 0) matrix[i][y] = 0;
        }
    }
}

void getReverseMatrix(int row, int col, double (*matrix)[col], double (*reverse_matrix)[col]){
    if(row != col){
        printf("\nrow != col");
        return;
    }
    double determinant = getDeterminant(row, col, matrix, false);
    if(determinant == 0){
        printf("\ndeterminant = 0");
        return;
    }

    double ac_matrix[row][col];
    for(int i=0;i<row;i++){
        for(int y=0;y<col;y++){
            double reduced_matrix[row-1][col-1];
            reduceMatrix(row, col, i, y, matrix, reduced_matrix);
            ac_matrix[i][y] = (pow(-1,(i+y)) * getDeterminant(row-1, col-1, reduced_matrix, false));
        }
    }
    double ac_matrix_T[row][col];
    transposed_matrix(row, col, ac_matrix, ac_matrix_T);

    prod_s(row, col, (1/determinant), ac_matrix_T, reverse_matrix);
    filter_reverse_matrix(row, col, reverse_matrix);
}

int count_zero(int len, double* row){
    int num_zero = 0;
    for(int i=0;i<len;i++){
        if(*(row+i) == 0) num_zero++;
        else return num_zero;
    }
    return num_zero;
}

void swap_row(int row1, int row2, int col, double (*matrix)[col]){
    double temp_arr[col];
    for(int i=0;i<col;i++){
        *(temp_arr+i) = matrix[row1][i];
    }
    for(int i=0;i<col;i++){
        matrix[row1][i] = matrix[row2][i];
        matrix[row2][i] = temp_arr[i];
    }
}

void fix_row(int row, int col, double (*matrix)[col]){
    for(int i=0;i<row-1;i++){
        for(int y=i;y<col-1;y++){
            double row1[col], row2[col];
            for(int j=0;j<col;j++){
                *(row1+j) = matrix[i][j];
                *(row2+j) = matrix[y+1][j];
            }
            if(count_zero(col, row1) > count_zero(col, row2)){
                swap_row(i, y+1, col, matrix);
            }
        }
    }
}
//Todo : control
void gauss_jordan(int row, int col, double (*matrix)[col], double (*gauss_matrix)[col], bool show_gauss_matrix){
    copy_matrix(row, col, matrix, gauss_matrix);
    fix_row(row, col, gauss_matrix);
    if(isMatrixScale(row, col, gauss_matrix)){
        if(show_gauss_matrix) printMatrix(row, col, gauss_matrix);
        return;
    }
    else{
        for(int i=0;i<row-1;i++){
            for(int y=0;y<col;y++){
                if(gauss_matrix[i][y] != 0){
                    for(int j=row-1;j>i;j--){
                        for(int z=y;z>=0;z--){
                            if(gauss_matrix[j][z] != 0){
                                double scale = (-(gauss_matrix[j][z])/gauss_matrix[i][y]);
                                double row1[col], row2[col], row_p[col], row_s[col];

                                copy_vector(col, gauss_matrix[i], row1);
                                copy_vector(col, gauss_matrix[j], row2);

                                prod_s_vector(col, row1, scale, row_p);
                                sum_vector(col, row2, row_p, row_s);

                                for(int k=0;k<col;k++){
                                    gauss_matrix[j][k] = *(row_s+k);
                                }
                            }
                        }
                    }
                    fix_row(row, col, gauss_matrix);
                    if(isMatrixScale(row, col, gauss_matrix)){
                        if(show_gauss_matrix) printMatrix(row, col, gauss_matrix);
                        return;
                    }
                    break;
                }
            }
        }
    }
    if(show_gauss_matrix) printMatrix(row, col, gauss_matrix);
}

int getRank(int row, int col, double (*matrix)[col], bool show_gauss_matrix){
    double gauss_matrix[row][col];
    gauss_jordan(row, col, matrix, gauss_matrix, show_gauss_matrix);

    int n_pivot = 0;
    for(int i=0;i<row;i++){
        for(int y=0;y<col;y++){
            if(gauss_matrix[i][y] != 0){
                n_pivot++;
                break;
            }
        }
    }
    return n_pivot;
}

bool isIndipendent(int len_vector, int args, ...){
    va_list list_args;
    double matrix_v[args][len_vector];

    va_start(list_args, args);
    for(int i=0;i<args;i++){
        double* vector = (double*) va_arg(list_args, double*);
        for(int y=0;y<len_vector;y++){
            matrix_v[i][y] = *(vector+y);
        }
    }
    int rank = getRank(args, len_vector, matrix_v, false);

    va_end(list_args);
    printf("\nSubspace Spanned = R^%d\n", rank);
    if(rank == args) return true;
    else return false;
}

void show(int row, int col, double (*matrix)[col], bool print_matrix, bool sarrus){
    if(print_matrix) printMatrix(row, col, matrix);
    int* degree = getDegree(row, col);
    printf("\nDegree : (%d, %d)", *(degree), *(degree+1));
    printf("\nDeterminant : %lf", getDeterminant(row, col, matrix, sarrus));
}

