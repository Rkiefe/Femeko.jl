#include <iostream>
#include "../extern/eigen/Eigen/Dense"

int main(int argc, char const *argv[])
{
	
	std::vector<double> array = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
	Eigen::MatrixXd p = Eigen::Map<Eigen::MatrixXd> (array.data(), 3, 2);

	std::cout << p; // correct format

	return 0;
}