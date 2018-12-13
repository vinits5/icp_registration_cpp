#include "ICP.h"
#include <iostream>
#include <eigen3/Eigen/Dense>
#include "../KDTree/kd_tree.h"
#include "../KDTree/type_defs.h"
using namespace Eigen;
using namespace std;

// Class Constructor.
ICP::ICP(MatrixXf *p0, MatrixXf *p1, KDTree *tree_M_sampled, MatrixXf *SigmaS){
	// Arguments->
	// *p0 		Sampled Model Points (Nx3)
	// *p1 		Sensor Points 		 (500x3)

	p0_ICP = p0; p1_ICP = p1; tree_M_sampled_ICP = tree_M_sampled; SigmaS_ICP = SigmaS;
	cout<<"ICP Class"<<endl;
}

// Find Rotation & Translation between sensor & target.
Transformation ICP::find_rigid_transform(MatrixXf *sensor, MatrixXf *target){
	MatrixXf A = *sensor;		// Sensor Data (500x3)
	MatrixXf B = *target;		// Target Model Data (500x3)
	VectorXf centroid_A = A.colwise().mean();			// Centroid of Sensor data. (3x1) (uA)
	VectorXf centroid_B = B.colwise().mean();			// Centroid of Target Model data. (3x1) (uB)
	
	A.transposeInPlace(); B.transposeInPlace();			// Convert data to (3x500)

	A = A.colwise()-centroid_A;							// Subtract mean from Sensor Data (3x500)
	B = B.colwise()-centroid_B;							// Subtract mean from Model Data (3x500)

	MatrixXf H = A*B.transpose();						// Matrix to find SVD (3x3)

	JacobiSVD<MatrixXf> svd(H, ComputeThinU | ComputeThinV);	// Computes U,Vt.T,diagonal of Sigma Matrices.
	Transformation T;
	T.rotation_matrix = svd.matrixV()*svd.matrixU().transpose();	// Compute Rotation Matrix (3x3)	(R = U*V.trans)
	
	if (T.rotation_matrix.determinant() < 0){
		// Special Reflection cases.
		// all rows and last column.
		MatrixXf temp = svd.matrixV();
		temp(0,2)=-1*temp(0,2); temp(1,2)=-1*temp(1,2); temp(2,2)=-1*temp(2,2);
		T.rotation_matrix = temp*svd.matrixU().transpose();
	}

	T.translation = centroid_B - T.rotation_matrix*centroid_A;			// Compute translation (3x1)		(t = uB-R*uA)
	return T;											// Return the structure containing Rotation & translation.
}

Matrix<float,4,4> ICP::icp_Rt_to_matrix(Transformation T){
	Matrix<float,4,4> transformation;
	transformation.block(0,0,3,3)=T.rotation_matrix;
	transformation(3,3)=1; transformation(3,0)=0; transformation(3,1)=0; transformation(3,2)=0; 
	transformation.block(0,3,3,1)=T.translation;
	return transformation;
}

ICP_Result ICP::compute(int iterations){
	float ftol = 1.0e-7;								// Define a tolerance limit.
	int dim_k = p0_ICP->cols();							// Store number of points in M_sampled matrix. 
	MatrixXf p = *p1_ICP;								// Copy the Sensor matrix to another matrix. (500x3)
	
	Matrix<float,4,4> g = MatrixXf::Identity(dim_k+1,dim_k+1);	// An identity matrix (4x4)
	vector<Matrix<float,4,4>> g_series;					// A vector to store transformation matrices (Nx4x4)
	g_series.push_back(g);								// Push the identity matrix to vector.

	ICP_Result result;

	// Iterations for ICP
	for(int itr=0; itr<iterations; itr++){
		PointCloud sensor_ptCloud(3,p.rows());					// Point Cloud data structure for Sensor data. (3x500)
		for(int i=0; i<p.rows(); i++){
			sensor_ptCloud.col(i)=p.row(i).cast<long double>();	// Convert MatrixXf Sensor data to PointCloud.
		}

		// Nearest Neighbour Search for Sensor Points from KDTree of Sampled Model Points.
		// target_ptCloud--> x,y,z,distance (Neighbouring Points from Sampled Model Points)
		MatrixXld target_ptCloud = kd_search_targets(&sensor_ptCloud, *tree_M_sampled_ICP);		// targets (neighbouring points from KDTree of Model points) (4x500)
		
		MatrixXf target(target_ptCloud.cols(),3);			// Find target points from target point cloud.	(500x3)
		for(int i=0; i<target_ptCloud.cols(); i++){
			target.row(i)=target_ptCloud.block(0,i,3,1).cast<float>().transpose();	// conversion from long double to float.
		}
		result.target = target;

		// Find the Rotation and Translation between Sensor Data & Target Data (Nearest Neighbour from Model Pts).
		Transformation T = ICP::find_rigid_transform(&p, &target);

		// Apply transformation to Sensor data.
		MatrixXf new_p = T.rotation_matrix*p.transpose();		// Apply rotation to sensor pt cloud (3x500)
		new_p = new_p.colwise()+T.translation;					// Apply translation to rotated sensor data (3x500)
		new_p.transposeInPlace();								// Transpose the sensor data (500x3)

		MatrixXf error_mat = (p-new_p).cwiseAbs();				// Take L1 norm b/w new sensor data & previous sensor data.
		result.error = error_mat.sum();
		if(error_mat.sum()<ftol){								// If below thershold break the iterations.
			break;
		}

		// Update for next iteration.
		p = new_p;												// Copy the New Sensor data to previous variable
		Matrix<float,4,4> dg = ICP::icp_Rt_to_matrix(T);		// Generate Transformation Matrix.
		g = dg*g;												// Apply it to previous transformations.
		g_series.push_back(g);									// Store the transformation in a vector.
		cout<<"Index No: "<<itr<<" & Error: "<<error_mat.sum()<<endl;
	}

	result.g = g;
	result.p = p;
	return result;
}