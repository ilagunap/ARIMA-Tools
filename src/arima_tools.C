/*
 * arima_tools.cpp
 *
 *  Created on: Jun 9, 2010
 *      Author: Ignacio Laguna
 *     Contact: ilaguna@purdue.edu
 */

#include "arima_tools.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <iostream>

using namespace std;

/*void ARIMA_predict(const ARIMAModel &model, unsigned int n,
		const vector<double> &obs, Predictions &pred)
{
	size_t ma_params = model.ma_params.size();
	size_t obs_size;
	vector<double> obs_trans;
	Series dseries;

	// Differentiate series
	if (model.d > 0) {
		dseries.diff(obs, obs_trans, model.d);
	} else {
		obs_trans = obs;
	}

	obs_size = obs_trans.size();

	// Initialize error values
	vector<double> errors(ma_params,0);

	// Create whole window of values: observed + predicted
	vector<double> win(obs_trans);
	for (size_t i=0; i < n; i++)
		win.push_back(0.0);

	if (ma_params > 0) {
		// Make predictions when q > 0
		size_t start = model.ar_params.size();
		for (size_t i=start; i < win.size(); i++) {

			double new_val = 0;
			double error_now = 0;
			double phi_factor = 1.0;

			// Do the AR part
			if (model.ar_params.size() > 0) {
				for (size_t j=0; j < model.ar_params.size(); j++) {
					new_val += model.ar_params[j] * win[i-1-j];
					phi_factor -= model.ar_params[j];
				}
			}

			// Do the MA part
			for (unsigned int j=0; j < ma_params; j++)
				new_val += model.ma_params[j] * errors[ma_params-1-j];

			// Add intersect
			new_val += model.intercept * phi_factor;

			// Update predictions
			if (i >= obs_size) {
				win[i] = new_val;
				error_now = 0;
			} else {
				error_now = win[i] - new_val;
			}

			// Update errors
			for (size_t i=0; i < ma_params - 1; i++)
				errors[i] = errors[i+1];

			errors[ma_params-1] = error_now;
		}
	} else if (ma_params == 0) {
		// Here q = 0, so we don't need the previous error values
		for (size_t i = obs_size; i < win.size(); i++) {
			double new_val = 0;
			double phi_factor = 1.0;

			// Do the AR part
			if (model.ar_params.size() > 0) {
				for (size_t j=0; j < model.ar_params.size(); j++) {
					new_val += model.ar_params[j] * win[i-1-j];
					phi_factor -= model.ar_params[j];
				}
			}

			// Add intersect
			new_val += model.intercept * phi_factor;

			// Update predictions
			win[i] = new_val;
		}
	}

	// Update prediction structure
	pred.values.clear();
	pred.error.clear();
	for (size_t i=0; i < n; i++) {
		pred.values.push_back(win[obs_size+i]);
		pred.error.push_back(1.96 * sqrt(model.sigma_2));
	}

	// Undo difference only if d > 0
	if (model.d > 0) {
		vector<double> tmp(pred.values);
		dseries.undo_diff(tmp, pred.values);
	}
}*/

Predictions ARIMAForecaster::forecast(const ARIMAModel &model, unsigned int n,
			const vector<double> &obs)
{
	Predictions pred;

	size_t ma_params = model.ma_params.size();
	size_t obs_size;
	vector<double> obs_trans;
	Series dseries;

	// Differentiate series
	if (model.d > 0) {
		dseries.diff(obs, obs_trans, model.d);
	} else {
		obs_trans = obs;
	}

	obs_size = obs_trans.size();

	// Initialize error values
	vector<double> errors(ma_params,0);

	// Create whole window of values: observed + predicted
	vector<double> win(obs_trans);
	for (size_t i=0; i < n; i++)
		win.push_back(0.0);

	if (ma_params > 0) {
		// Make predictions when q > 0
		size_t start = model.ar_params.size();
		for (size_t i=start; i < win.size(); i++) {

			double new_val = 0;
			double error_now = 0;
			double phi_factor = 1.0;

			// Do the AR part
			if (model.ar_params.size() > 0) {
				for (size_t j=0; j < model.ar_params.size(); j++) {
					new_val += model.ar_params[j] * win[i-1-j];
					phi_factor -= model.ar_params[j];
				}
			}

			// Do the MA part
			for (unsigned int j=0; j < ma_params; j++)
				new_val += model.ma_params[j] * errors[ma_params-1-j];

			// Add intersect
			new_val += model.intercept * phi_factor;

			// Update predictions
			if (i >= obs_size) {
				win[i] = new_val;
				error_now = 0;
			} else {
				error_now = win[i] - new_val;
			}

			// Update errors
			for (size_t i=0; i < ma_params - 1; i++)
				errors[i] = errors[i+1];

			errors[ma_params-1] = error_now;
		}
	} else if (ma_params == 0) {
		// Here q = 0, so we don't need the previous error values
		for (size_t i = obs_size; i < win.size(); i++) {
			double new_val = 0;
			double phi_factor = 1.0;

			// Do the AR part
			if (model.ar_params.size() > 0) {
				for (size_t j=0; j < model.ar_params.size(); j++) {
					new_val += model.ar_params[j] * win[i-1-j];
					phi_factor -= model.ar_params[j];
				}
			}

			// Add intersect
			new_val += model.intercept * phi_factor;

			// Update predictions
			win[i] = new_val;
		}
	}

	// Update prediction structure
	pred.values.clear();
	pred.error.clear();
	for (size_t i=0; i < n; i++) {
		pred.values.push_back(win[obs_size+i]);
		pred.error.push_back(1.96 * sqrt(model.sigma_2));
	}

	// Undo difference only if d > 0
	if (model.d > 0) {
		vector<double> tmp(pred.values);
		dseries.undo_diff(tmp, pred.values);
	}

	return pred;
}

void Series::diff(const vector<double> &old_series,
		vector<double> &new_series, unsigned int d)
{
	unsigned int s = old_series.size();
	diff_d = d;

	if (d >= s) {
		printf("[ARIMA_TOOLS] Cannot differentiate series\n");
		fflush(stdout);
		exit(-1);
	}

	// Save original window
	vector<double> tmp;
	for (unsigned int i=0; i < s; i++)
		tmp.push_back(old_series[i]);
	dseries.win.push_back(tmp);

	// Differentiate d times
	for (unsigned int i=0; i < d; i++) {
		tmp.clear();
		for (unsigned int j=0; j < (dseries.win[i].size() - 1); j++) {
			tmp.push_back(dseries.win[i][j+1] - dseries.win[i][j]);
		}
		dseries.win.push_back(tmp);
	}

	new_series = dseries.win[d];
}

void Series::print(void)
{
	cout << "Printing series ..." << endl;
	for (unsigned int i=0; i < dseries.win.size(); i++) {
		for (unsigned int j=0; j < dseries.win[i].size(); j++) {
			cout << dseries.win[i][j] << " ";
		}
		cout << endl;
	}
}

void Series::undo_diff(vector<double> trans_predict, vector<double> &real_predict)
{
	if (diff_d <= 0) {
		printf("[ARIMA_TOOLS] Cannot undo differentiated series\n");
		fflush(stdout);
		exit(-1);
	}

	unsigned int size = trans_predict.size();

	real_predict.clear();
	for (unsigned int i=0; i < size; i++)
		real_predict.push_back(0);

	// Undo d-differentiated series
	for (unsigned int i=0; i < diff_d; i++) {
		unsigned int s = dseries.win.size();

		for (unsigned int k=0; k < size; k++) {
			double prev = 0;
			if (k == 0) prev = dseries.win[s-i-2].back();
			else prev = real_predict[k-1];

			real_predict[k] = prev + trans_predict[k];
		}

		trans_predict = real_predict;
	}
}
