# include <stdio.h>
# include <stdlib.h>
# include <xgboost/c_api.h>

int main(int argc, char** argv)
{
	int silent = 0;
	int use_gpu = 1;
	
	DMatrixHandle dtrain, dtest;
	XGDMatrixCreateFromFile("agaricus.txt.train", silent, &dtrain);
	XGDMatrixCreateFromFile("agaricus.txt.test", silent, &dtest);

	BoosterHandle booster;
	DMatrixHandle eval_dmats[2] = {dtrain, dtest};
	XGBoosterCreate(eval_dmats, 2, &booster);

	XGBoosterSetParam(booster, "tree_method", use_gpu ? "gpu_hist" : "hist");

	if (use_gpu)
	{
		XGBoosterSetParam(booster, "gpu_id", "0");
	} else {
		XGBoosterSetParam(booster, "gpu_id", "-1");
	}

	
	XGBoosterSetParam(booster, "objective", "binary:logistic");
	XGBoosterSetParam(booster, "min_child_weight", "1");
	XGBoosterSetParam(booster, "gamma", "0.1");
	XGBoosterSetParam(booster, "max_depth", "3");
	XGBoosterSetParam(booster, "verbosity", silent ? "0" : "1");

	int n_trees = 10;
	const char* eval_names[2] = {"train", "test"};
	const char* eval_result = NULL;
	for(int i = 0; i < n_trees; ++i)
	{
		XGBoosterUpdateOneIter(booster, i, dtrain);
		XGBoosterEvalOneIter(booster, i, eval_dmats, eval_names, 2, &eval_result);
		printf("%s/n", eval_result);
	}

	bst_ulong out_len = 0;
	const float* out_result = NULL;
	int n_print = 10;

	XGBoosterPredict(booster, dtest, 0, 0, 0, &out_len, &out_result);
	printf("y_pred: ");
	for(int i = 0; i < n_print ; ++i)
	{
		printf("%1.4f ", out_result[i]);
	}

	printf("\n");

	XGBoosterFree(booster);
    XGDMatrixFree(dtrain);
	XGDMatrixFree(dtest);
	return 0;

}
