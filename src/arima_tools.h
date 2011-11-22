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

/**
 * ARIMA(p,d,q) model structure.
 * Defines the parameters (p,d,q) of an ARIMA model (already estimated).
 */
typedef struct {
	/** AR (autoregressive) part parameters */
	std::vector<double> ar_params;

	/** MA (moving average) part parameters */
	std::vector<double> ma_params;

	/** Integrated part */
	unsigned int d;

	/** Intercept (constant) */
	double intercept;

	/** Error obtained in the training part */
	double sigma_2;
} ARIMAModel;

/**
 * Vector of forecasted values.
 * Contains a vector for the forecast values and another vector for error values.
 */
typedef struct {
	/** Forecasted values */
	std::vector<double> values;

	/** Error values */
	std::vector<double> error;
} Predictions;

/**
 * Forecast 'n' values using a predefined ARIMA model.
 * @param model An ARIMAModel structure.
 * @param n Number of forecasts ahead in time.
 * @param obs Observations already seen.
 * @param pred Structure where forecasts will be saved in.
 */
//void ARIMA_predict(const ARIMAModel &model, unsigned int n,
		//const std::vector<double> &obs, Predictions &pred);


/**
 * Forecast observations using a predefined ARIMA model.
 */
class ARIMAForecaster {
public:
	/**
	 * Forecast 'n' values.
	 * @param model An ARIMAModel structure.
	 * @param n Number of forecasts ahead in time.
	 * @param obs Observations already seen.
	 * @return Predictions structure with forecasts and errors.
	 */
	static Predictions forecast(const ARIMAModel &model, unsigned int n,
			const std::vector<double> &obs);
};

/**
 * Vector with differentiated series.
 */
typedef struct {
	/** Matrix of values */
	std::vector< std::vector<double> > win;
} DiffSeries;

/**
 * Differentiate and undo-difference a series.
 */
class Series {
	DiffSeries dseries;
	unsigned int diff_d;
public:
	/**
	 * Differentiate a series.
	 * @param old_series Original series.
	 * @param new_series New (differentiated) series.
	 * @param d Degree.
	 */
	void diff(const std::vector<double> &old_series,
			std::vector<double> &new_series, unsigned int d);

	/**
	 * Undo-difference a series.
	 * @param trans_predict Converted series.
	 * @param real_predict Original series.
	 */
	void undo_diff(std::vector<double> trans_predict,
			std::vector<double> &real_predict);

	/**
	 * Print series.
	 */
	void print();
};

#endif /* ARIMA_TOOLS_H_ */
