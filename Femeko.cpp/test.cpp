#include "src/femeko.h"

int main(int argc, char const *argv[])
{
	
	std::vector<double> array = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
	Eigen::MatrixXd p = Eigen::Map<Eigen::MatrixXd> (array.data(), 3, 2);

	std::cout << p  << std::endl; // correct format

	// You could do it manually by
	Eigen::MatrixXd pManual = Eigen::MatrixXd::Zero(3, 2);
	int n = 0;
	for(int nd = 0; nd<2; nd++){ // loop over the nodes
		for(int i = 0; i<3; i++){ // Loop over x, y, z
			pManual(i, nd) = array[n];
			n++;
		}
	}

	std::cout << pManual << std::endl; // Also correct format

	return 0;
}