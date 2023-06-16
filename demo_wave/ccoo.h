typedef struct{
  int nrow, ncol;//matrix size m-rows, n-columns
  int nnz;//number of non-zeros
  int *row_ind;//row indices
  int *col_ind;//column indices
  complex *val;//values of the matrix A
} ccoo_t;//COOrdinate format

