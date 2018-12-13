#ifndef ICP_H
#define ICP_H
#include <eigen3/Eigen/Dense>
#include "../KDTree/kd_tree.h"
using namespace Eigen;

struct Transformation{
	Matrix<float,3,3> rotation_matrix;
	Matrix<float,3,1> translation;
};

struct ICP_Result{
	Matrix<float,4,4> g;
	MatrixXf p;
	float error;
	MatrixXf target;
};


class ICP{
public:
	MatrixXf *p0_ICP,*p1_ICP,*SigmaS_ICP;
	KDTree *tree_M_sampled_ICP;

	ICP(MatrixXf*, MatrixXf*, KDTree*, MatrixXf*);

	ICP_Result compute(int);

	Transformation find_rigid_transform(MatrixXf*, MatrixXf*);
	Matrix<float,4,4> icp_Rt_to_matrix(Transformation);
};
#endif