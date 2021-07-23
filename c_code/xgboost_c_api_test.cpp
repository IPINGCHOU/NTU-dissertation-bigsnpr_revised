// Template.cpp by kevinptt 20170429
# include <bits/stdc++.h>
//# include </opt/xgboost/include/xgboost/c_api.h>
# include <xgboost/c_api.h>
# include <stdio.h>

const int cols=3, rows=5;
float train[rows][cols];
float train_labels[rows];

int main()
{
	for(int i = 0; i < rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			train[i][j] = (i+1) * (j+1);
		}
	}

	for(int i = 0; i < rows; i++)
	{
		train_labels[i] = 1+i*i*i;
	}
	
	DMatrixHandle h_train[1];
	XGDMatrixCreateFromMat((float * ) train, rows, cols, -1, &h_train[0]);
	XGDMatrixSetFloatInfo(h_train[0], "label", train_labels, rows);

	bst_ulong bst_result;
	const float *out_floats;
	XGDMatrixGetFloatInfo(h_train[0], "label", &bst_result, &out_floats);
	for(unsigned int i = 0; i < bst_result; i++)
	{
		std::cout << "label[" << i << "]=" << out_floats[i] << std::endl;
	}

	BoosterHandle h_booster;
	XGBoosterCreate(h_train, 1, &h_booster);
	XGBoosterSetParam(h_booster, "booster", "gbtree");
	XGBoosterSetParam(h_booster, "objective", "reg:linear");
	XGBoosterSetParam(h_booster, "max_depth", "5");
	XGBoosterSetParam(h_booster, "eta", "0.1");
	XGBoosterSetParam(h_booster, "min_child_weight", "1");
	XGBoosterSetParam(h_booster, "subsample", "0.5");
	XGBoosterSetParam(h_booster, "colsample_bytree", "1");
	XGBoosterSetParam(h_booster, "num_parallel_tree", "1");

	for(int iter = 0; iter < 200 ; iter++)
	{
		XGBoosterUpdateOneIter(h_booster, iter, h_train[0]);
	}
	const int sample_rows = 5; 
	float test[sample_rows][cols];
	for(int i = 0; i < sample_rows; i++)
	{
		for(int j = 0; j < cols; j++)
		{
			test[i][j] = (i+1) * (j+1);
		}
	}

	DMatrixHandle h_test;
	XGDMatrixCreateFromMat((float *) test, sample_rows, cols, -1 , &h_test);
	bst_ulong out_len;
	const float* f;
	XGBoosterPredict(h_booster, h_test, 0, 0, 0, &out_len, &f);
	
	for(unsigned int i = 0; i < out_len; i++)
	{
		std::cout << "prediction[" << i << "]=" << f[i] << std::endl;
	}

	XGDMatrixFree(h_train[0]);
	XGDMatrixFree(h_test);
	XGBoosterFree(h_booster);

	
	return 0;
}

