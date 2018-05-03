#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
    */

	  VectorXd rmse(4);
		rmse << 0,0,0,0;

		// assert the following inputs:
		// estimation vector size should not be zero
		// estimation vector size should equal ground truth vector size
		if(estimations.size() != ground_truth.size()
				|| estimations.size() == 0){
			std::cout << "Invalid estimation or ground_truth data" <<std::endl;
			return rmse;
		}

		//accumulate squared difference
		for(uint8_t i=0; i < estimations.size(); ++i){
			VectorXd difference = estimations[i] - ground_truth[i];

			//coefficient-wise multiplication
			difference = difference.array()*difference.array();
			rmse += difference;
		}

		//find the mean and squared root of mean
		rmse = rmse/estimations.size();
		rmse = rmse.array().sqrt();
	  return rmse;

}
