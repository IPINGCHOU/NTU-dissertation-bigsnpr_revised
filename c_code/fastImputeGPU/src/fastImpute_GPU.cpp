#include <Rcpp.h>
#include <xgboost/c_api.h>
#include <xgboost_R.h>
#include <dmlc/omp.h>
#include <dmlc/logging.h>
#include <dmlc/base.h>
#include <stdio.h>
#include <vector>
#include <string>
#include <utility>
#include <cstring>
#include <cstdio>
#include <iostream>
#include <sstream>

#define safe_xgboost(call) {                                            \
int err = (call);                                                       \
if (err != 0) {                                                         \
  fprintf(stderr, "%s:%d: error in %s: %s\n", __FILE__, __LINE__, #call, XGBGetLastError()); \
  exit(1);                                                              \
}                                                                       \
}

using namespace Rcpp;

inline SEXP getAttrib(SEXP a, SEXP b) {
  return Rf_getAttrib(a, b);
}

std::vector<float> R_to_float_2D(SEXP mat, size_t nrow, size_t ncol)
{
  // SEXP dim = getAttrib(mat, R_DimSymbol);
  // size_t nrow = static_cast<size_t>(INTEGER(dim)[0]);
  // size_t ncol = static_cast<size_t>(INTEGER(dim)[1]);
  const bool is_int = TYPEOF(mat) == INTSXP;
  double *din;
  int *iin;
  if (is_int) {
    iin = INTEGER(mat);
  } else {
    din = REAL(mat);
  }
  std::vector<float> data(nrow * ncol);
  #pragma omp parallel for schedule(static)
    for (dmlc::omp_ulong i = 0; i < nrow; ++i) {
      for (size_t j = 0; j < ncol; ++j) {
        data[i * ncol +j] = is_int ? static_cast<float>(iin[i + nrow * j]) : din[i + nrow * j];
      }
    }
  return data;
}

std::vector<float> R_to_float_1D(SEXP mat, size_t v_len)
{
  const bool is_int = TYPEOF(mat) == INTSXP;
  double *din;
  int *iin;
  if (is_int) {
    iin = INTEGER(mat);
  } else {
    din = REAL(mat);
  }
  std::vector<float> data(v_len);
  #pragma omp parallel for schedule(static)
    for (dmlc::omp_ulong i = 0; i < v_len; ++i) {
      data[i] = is_int ? static_cast<float>(iin[i]) : din[i];
    }
  return data;
}

float get_base_score(std::vector<float> data)
{
  float average = accumulate( data.begin(), data.end(), 0.0) / data.size();
  average = std::max(average, (float)(1e-7));
  average = std::min(average, (float)(1-1e-7));
  return average;
}

std::tuple<size_t, size_t> get_mat_info(SEXP mat) {
  SEXP dim = getAttrib(mat, R_DimSymbol);
  size_t nrow = static_cast<size_t>(INTEGER(dim)[0]);
  size_t ncol = static_cast<size_t>(INTEGER(dim)[1]);
  return  std::make_tuple(nrow, ncol);
}

size_t get_vec_size(SEXP x) {
  return Rf_xlength(x);
}

void float_to_NumericVector(const float* x, int n, double* y) {
  for(int i = 0; i < n; ++i) y[i] = x[i];
}

// [[Rcpp::export]]
List cpp_xgboost_gpu(SEXP train, SEXP train_labels, SEXP valid, SEXP test, SEXP mode) {
  
  // initialize parameters
  int silent = 0;
  int n_print = 10;
  std::string input_mode = as<std::string>(mode); 
  std::string cpu_string = "cpu";
  std::string gpu_string = "gpu";

  std::vector<float> ftrain, ftrain_labels, fvalid, ftest;
  float base_score;
  size_t train_nrow, train_ncol;
  double missing = NA_REAL;

  // get nrow and ncol from train_df
  // transfer SEXP train, train_lables std::vector<float>
  std::tie(train_nrow, train_ncol) = get_mat_info(train);
  size_t v_size = get_vec_size(train_labels);
  ftrain = R_to_float_2D(train, train_nrow, train_ncol);
  ftrain_labels = R_to_float_1D(train_labels, v_size);

  // creat DMatrixHandle for train
  DMatrixHandle dtrain;
  safe_xgboost(XGDMatrixCreateFromMat(dmlc::BeginPtr(ftrain), train_nrow, train_ncol, missing, &dtrain));
  safe_xgboost(XGDMatrixSetFloatInfo(dtrain, "label", dmlc::BeginPtr(ftrain_labels), train_nrow));
  base_score = get_base_score(ftrain_labels);

  // settings
  BoosterHandle booster;
  DMatrixHandle create_dtrain[1] = {dtrain};
  safe_xgboost(XGBoosterCreate(create_dtrain, 1, &booster));
  safe_xgboost(XGBoosterSetParam(booster, "objective", "binary:logistic"));

  if(input_mode == gpu_string)
  {
    safe_xgboost(XGBoosterSetParam(booster, "tree_method", "gpu_hist"));
    safe_xgboost(XGBoosterSetParam(booster, "gpu_id", "0"));
  } else if(input_mode == cpu_string) {
    safe_xgboost(XGBoosterSetParam(booster, "tree_method", "hist"));
  } else {
    printf("What's your mode? keyin cpu or gpu (lower case)");
    exit(1);
  }


  // base score
  char buffer[64];
  snprintf(buffer, sizeof buffer, "%f", base_score);
  safe_xgboost(XGBoosterSetParam(booster, "base_score", buffer));
  safe_xgboost(XGBoosterSetParam(booster, "max_depth", "4"));
  safe_xgboost(XGBoosterSetParam(booster, "n_thread", "1"));
  safe_xgboost(XGBoosterSetParam(booster, "verbosity", silent ? "0" : "1"));
  safe_xgboost(XGBoosterSetParam(booster, "save_period", "0"));

  // update
  
  int n_rounds = 10;
  for(int i = 0; i < n_rounds; ++i)
  {
    safe_xgboost(XGBoosterUpdateOneIter(booster, i, dtrain));
  }
  
  // // predict
  size_t valid_nrow, valid_ncol, test_nrow, test_ncol;
  std::tie(valid_nrow, valid_ncol) = get_mat_info(valid);
  std::tie(test_nrow, test_ncol) = get_mat_info(test);

  DMatrixHandle dvalid, dtest;
  fvalid = R_to_float_2D(valid, valid_nrow, valid_ncol);
  ftest = R_to_float_2D(test, test_nrow, test_ncol);
  
  safe_xgboost(XGDMatrixCreateFromMat(dmlc::BeginPtr(fvalid), valid_nrow, valid_ncol, missing, &dvalid));
  safe_xgboost(XGDMatrixCreateFromMat(dmlc::BeginPtr(ftest), test_nrow, test_ncol, missing, &dtest));
  
  bst_ulong valid_out_len = 0, test_out_len = 0;
  const float* valid_result = NULL;
  const float* test_result = NULL;

  safe_xgboost(XGBoosterPredict(booster, dvalid, 0, 0, 0, &valid_out_len, &valid_result));
  NumericVector return_valid(valid_out_len);
  float_to_NumericVector(valid_result, valid_out_len, REAL(return_valid));
  safe_xgboost(XGBoosterPredict(booster, dtest, 0, 0, 0, &test_out_len, &test_result));
  NumericVector return_test(test_out_len);
  float_to_NumericVector(test_result, test_out_len, REAL(return_test));

  safe_xgboost(XGBoosterFree(booster));
  safe_xgboost(XGDMatrixFree(dtrain));
  safe_xgboost(XGDMatrixFree(dvalid));
  safe_xgboost(XGDMatrixFree(dtest));

  return List::create(Named("pred2") = return_valid,
                      Named("pred") = return_test);
}