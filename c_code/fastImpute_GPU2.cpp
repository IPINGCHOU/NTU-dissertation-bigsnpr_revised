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
using namespace Rcpp;

inline SEXP getAttrib(SEXP a, SEXP b) {
  return Rf_getAttrib(a, b);
}

void R_to_handle(SEXP mat, int mode, DMatrixHandle out)
{
    SEXP dim = getAttrib(mat, R_DimSymbol);
    size_t nrow = static_cast<size_t>(INTEGER(dim)[0]);
    size_t ncol = static_cast<size_t>(INTEGER(dim)[1]);
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

    if(mode == 1)
    {
        XGDMatrixCreateFromMat(dmlc::BeginPtr(data), nrow, ncol, -1, &out);
    } else if (mode == 2) {
        XGDMatrixSetFloatInfo(out, "label", dmlc::BeginPtr(data), nrow);
    }
}


// [[Rcpp::export]]
List cpp_xgboost_gpu(SEXP train, SEXP train_labels, CharacterVector base_score, SEXP valid, SEXP test) {
    int silent = 0;
    int use_gpu = 0;

    DMatrixHandle dtrain;
    // XGDMatrixCreateFromMat((float *)train, nrow, ncol, -1, &dtrain[0]);
    R_to_handle(train, 1, dtrain);
    R_to_handle(train_labels, 2, dtrain);


    // // read back the labels, just a sanity check
    // bst_ulong bst_result;
    // const float *out_floats;
    // XGDMatrixGetFloatInfo(dtrain[0], "label" , &bst_result, &out_floats);
    // for (unsigned int i=0;i<bst_result;i++)
    //     std::cout << "label[" << i << "]=" << out_floats[i] << std::endl;

    // // settings
    BoosterHandle booster;
    XGBoosterSetParam(booster, "tree_method", use_gpu ? "gpu_hist" : "hist");
    XGBoosterSetParam(booster, "gpu_id", "0");

    DMatrixHandle create_dtrain[1] = {dtrain};
    XGBoosterCreate(create_dtrain, 1, &booster);
    XGBoosterSetParam(booster, "objective", "binary:logistic");
    // // base score ?
    // XGBoosterSetParam(booster, "base_score", base_score);
    XGBoosterSetParam(booster, "max_depth", "4");
    XGBoosterSetParam(booster, "n_thread", "1");
    XGBoosterSetParam(booster, "verbosity", silent ? "0" : "1");
    XGBoosterSetParam(booster, "save_period", "0");

    // // update
    int n_rounds = 10;
    for(int i = 0; i < n_rounds; ++i)
    {
        XGBoosterUpdateOneIter(booster, i, dtrain);
    }

    // // predict
    DMatrixHandle dvalid, dtest;
    R_to_handle(valid, 1, dvalid);
    R_to_handle(test, 1, dtest);

    bst_ulong valid_out_len = 0, test_out_len = 0;
    const float* valid_result = NULL;
    const float* test_result = NULL;
    XGBoosterPredict(booster, dvalid, 0, 0, 0, &valid_out_len, &valid_result);
    XGBoosterPredict(booster, dtest, 0, 0, 0, &test_out_len, &test_result);

    XGBoosterFree(booster);
    XGDMatrixFree(dtrain);
    XGDMatrixFree(dvalid);
    XGDMatrixFree(dtest);

    NumericVector return_valid = *valid_result;
    NumericVector return_test = *test_result;

    return List::create(Named("valid") = return_valid,
                        Named("inputation") = return_test);
}