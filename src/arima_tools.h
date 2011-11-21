/*
 * arima_tools.h
 *
 *  Created on: Jun 8, 2010
 *      Author: Ignacio Laguna
 *     Contact: ilaguna@purdue.edu
 */

#ifndef ARIMA_TOOLS_H_
#define ARIMA_TOOLS_H_

#include <vector>

typedef struct ARIMAModel_ {
	std::vector<double> ar_params;
	std::vector<double> ma_params;
	unsigned int d;
	double intercept;
	double sigma_2;
} ARIMAModel;

typedef struct Predictions_ {
	std::vector<double> values;
	std::vector<double> error;
} Predictions;

typedef struct DiffSeries_ {
	std::vector< std::vector<double> > win;
} DiffSeries;

/*
 * Difference and Undo-difference a series
 */
class Series {
	DiffSeries dseries;
	unsigned int diff_d;
public:
	void diff(const std::vector<double> &old_series, std::vector<double> &new_series, unsigned int d);
	void undo_diff(std::vector<double> trans_predict, std::vector<double> &real_predict);
	void print(void);
};

/*
 * Predicts 'n' values using an ARIMA model
 */
void ARIMA_predict(const ARIMAModel &model, unsigned int n,
		const std::vector<double> &obs, Predictions &pred);

#endif /* ARIMA_TOOLS_H_ */
