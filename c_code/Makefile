SRC=fastImpute_GPU2.cpp
TGT=fastImpute_GPU2.o

cc=g++
CFLAGS ?=-O2

XGBOOST_ROOT ?=/home/bigsnpr_medicine/c_code/xgboost
R_INCLUDE_DIR = -I"/usr/local/lib/R/site-library/Rcpp/include" -I"/home/bigsnpr_medicine/c_code" -I/usr/share/R/include -I/usr/local/lib/R/site-library/Rcpp/include -I/usr/local/include
XGBOOST_INCLUDE_DIR=-I$(XGBOOST_ROOT)/include -I$(XGBOOST_ROOT)/dmlc-core/include -I$(XGBOOST_ROOT)/rabit/include -I$(XGBOOST_ROOT)/R-package/src
LIB_DIR=-L$(XGBOOST_ROOT)/lib

build: $(TGT)

$(TGT): $(SRC) Makefile
	$(cc) -std=gnu++11 -I"/usr/share/R/include"-DNDEBUG $(R_INCLUDE_DIR) $(XGBOOST_INCLUDE_DIR) $(LIB_DIR) -fpic -g -c $(SRC) -o $(TGT) -lxgboost
