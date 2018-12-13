#include <iostream>
#include <string>
#include <eigen3/Eigen/Dense>
#include <fstream>
#include <cstdio>
#include <algorithm>
#include "KDTree/kd_tree.h"
#include "ICP/ICP.h"

using namespace std;
using namespace Eigen;

// Function Definitions
void data_shape(string, int *, int *);
void icp_test(MatrixXf *, MatrixXf *, MatrixXf *, KDTree *, MatrixXf *, KDTree *, MatrixXf *, MatrixXf *, int *, int *, int *);
void printMatrixSize(MatrixXf *);
MatrixXf read_file(string);
void sampleModelPoints(MatrixXf *, MatrixXf *, int *);
void write_data(MatrixXf *, string);

// Function Body defined here (arranged in alphabetical order):

// Used to calculate the columns and rows in .txt file.
void data_shape(string file_name, int *row_size, int *col_size){
	// Args.
	// file_name: name of file in string.
	// row_size: memory location to store no. of rows. (values seperated by ',')
	// col_size: memory location to store no. of cols. (total no of lines in a file)

	ifstream file(file_name);						// Open a file
	string x;
	getline(file,x,'\n');							// Read first line and store in x.
	*col_size = count(x.begin(),x.end(),',');		// count number of commas appeared in first line.
	*col_size = *col_size+1;						// no of columns = no of commas in first line + 1
	*row_size = *row_size+1;						// add 1 to no of rows as first line is read.

	while(!file.eof()){								// till the end of file.
		*row_size= *row_size+1;						// keep adding 1 to no of rows.
		getline(file,x,'\n');						// read each line.
	}
	file.close();									// close the file.
}

// ICP Function
void icp_test(MatrixXf *V, MatrixXf *S, MatrixXf *M, KDTree *tree_M, MatrixXf *M_sampled, KDTree *tree_M_sampled, MatrixXf *SigmaS, 
	MatrixXf *F, int *points_per_face, int *ICP_triangle_proj_switch, int *ICP_or_GICP_switch){
	// Some basic conversions to manage the matrix dimensions.
	V->transposeInPlace();				// N1x3 (vertices)
	S->transposeInPlace();				// 500x3 (sensor pts)
	M->transposeInPlace();				// N2x3 (model pts)
	M_sampled->transposeInPlace();		// Nix3 (sampled model pts)

	// Decide which class to be used.
	if (*ICP_triangle_proj_switch == 1){
		cout<<"Use Triangle Projection!"<<endl;
	}
	else{
		// Normal ICP method
		ICP icp = ICP(M_sampled, S, tree_M_sampled, SigmaS);	// Create object of ICP class.
		ICP_Result result = icp.compute(20);					// Call the compute method with 20 iterations.
		write_data(&result.p,"output.txt");
		MatrixXf transformation_data = result.g;
		write_data(&transformation_data,"transformation.txt");
	}
}

// Just to check the size of given matrix.
void printMatrixSize(MatrixXf *mat){
	// Just to save memory give the address of matrix.
	// mat->rows()		mat: Matrix class in Eigen library & rows() is a function in that class. 
	// calling function of a pointer needs arrow(->) operator.
	cout<<"Rows: "<<mat->rows()<<" & "<<"Cols: "<<mat->cols()<<endl;
}

// Function to read given file and return data as a Matrix in Eigen library.
MatrixXf read_file(string file_name){
	// Args.
	// file_name: name of file in string.
	// Output
	// MatrixXf: matrix data type in eigen library with the data in given file.

	int row_size=0,col_size=0;						// to define shape of matrix.
	data_shape(file_name,&row_size,&col_size);		// to find shape of matrix.
	MatrixXf data(row_size,col_size);				// define a matrix.
	ifstream file(file_name);						// read the file.
	string x;										// to store the string from the file.
	for(int i=0; i<row_size; i++){					// loop to read rows.
		for(int j=0; j<col_size; j++){				// loop to read cols.
			if(j!=(col_size-1)){					// delimiter is ',' if not the last column.	
				getline(file,x,',');				// store the string in x.
			}
			else{									// delimiter is '\n' if it is the last column.
				getline(file,x,'\n');				// store the string in x.
			}
			data(i,j) = stof(x);					// convert string to float value in matrix.
		}
	}
	file.close();									// close the file.
	return data;									// return the matrix.
}

// Sample model points after certain interval.
void sampleModelPoints(MatrixXf *M, MatrixXf *M_sampled, int *num_sampled_model_points){
	// Args.
	// M:				 			Model points (3xN)
	// M_sampled: 					Sampled Model Points (3xN1)
	// num_sampled_model_points: 	Interval to sample points from M.
	
	// #######################
	// Stopped here. Error: coeff is just a read-only method. Cannot assign new value to an element of matrix using it.
	// #######################
	for(int i=0; i<M_sampled->cols(); i++){
		M_sampled->block(0,i,1,1) = M->block(0,i*(*num_sampled_model_points),1,1);
		M_sampled->block(1,i,1,1) = M->block(1,i*(*num_sampled_model_points),1,1);
		M_sampled->block(2,i,1,1) = M->block(2,i*(*num_sampled_model_points),1,1);
	}
}

// Write Registered Point Cloud to a text file.
void write_data(MatrixXf *data, string file_name){
	// Args.
	// data			Registered Sensor Data. (Nx3)

	ofstream file(file_name);					// Open a file to write the data.
	int rows = data->rows(), cols = data->cols();	// Store the rows and columns in a point cloud.
	for(int i=0; i<rows; i++){					// Loop to write rows.
		for(int j=0; j<cols; j++){				// Loop to write cols.
			file<<data->coeff(i,j);	// Write each co-ordinate in a point.
			if(j<cols-1){
				file<<",";						// Add delimiter "," after each co-ordinate.
			}
		}
		if(i<rows-1){
			file<<"\n";							// Add delimiter "\n" after each point.
		}
	}
	file.close();
}


// Main code:
int main(){
	// Read all the files.
	MatrixXf S_given = read_file("armadillo20deg/S.txt");				// shape: (Nx3)
	MatrixXf V_given = read_file("armadillo20deg/M.txt");				// shape: (N1x3)
	MatrixXf M_given = read_file("armadillo20deg/Msampled.txt");		// shape: (N2x3)
	MatrixXf F = read_file("armadillo20deg/F.txt");					// Tejas stored F_given to F without having operation on it. (binary faces)
	MatrixXf gt_given = read_file("armadillo20deg/gt.txt");			// shape: (4x4)
	MatrixXf B_given = read_file("armadillo20deg/B.txt");			// shape: (N3x6)
	MatrixXf SigmaS_given = read_file("armadillo20deg/SigmaS.txt");	// shape: (N3x6)	(square root of covariance matrix)

	// Convert the data as per requirement.
	MatrixXf V = V_given.transpose();									// shape: (3xN1)	(Vertices)
	MatrixXf S = S_given.block(0,0,500,S_given.cols()).transpose();		// shape: (3x500)	(sensor pts.)
	MatrixXf M = M_given.transpose();									// shape: (3xN2)	(model pts.)

	// Matrix F is in data read process.

	// Ns: 500, Nm: no of model pts, Nv: no of vertices, Nf: rows in binary faces.
	int Ns = S.cols(), Nm = M.cols(), Nv = V.cols(), Nf = F.rows();		// Store the number of points in each data.
	int num_sampled_model_points = Nm/300;								// Step value to sample model pts from mesh data.
	MatrixXf M_sampled(M.rows(),(Nm/num_sampled_model_points)+1);		// Matrix to store the sampled model pts. (3xNi)
	sampleModelPoints(&M, &M_sampled, &num_sampled_model_points);		// Sample the model points from M matrix.

	int ICP_or_GICP_switch = 1; 	// 1 for ICP & 2 for GICP
	int ICP_triangle_proj_switch = 0; 	// 1 to use triangle projection 0 to not
	int show_plot_switch = 1; 

	int points_per_face = M.cols()/F.rows();	// Points per face
	cout<<"Points Per Face: "<<points_per_face<<endl;

	KDTree tree_M=NULL,tree_M_sampled=NULL;
	
	// Create KD tree of Matrix M_sampled_tree as tree_M_sampled.
	for(int i=0; i<M_sampled.cols(); i++){					// M_sampled (3xNi)
		// 1. Find ith column from M_sampled & cast it to Matrix<long double,3,1>
		// 2. Call insert function from kd_tree.h for each point.
		// 3. Arguments for insert-> i] 3D point as Matrix<long double,3,1>
								//  ii] Memory location of KDTree.

		insert(M_sampled.col(i).cast<long double>(),&tree_M_sampled);
	}

	// Create KD tree of Matrix M named as tree_M.
	for(int i=0; i<M.cols(); i++){
		insert(M.col(i).cast<long double>(),&tree_M);
	}

	// Call icp_test function.
	icp_test(&V, &S, &M, &tree_M, &M_sampled, &tree_M_sampled, &SigmaS_given, &F, &points_per_face, &ICP_triangle_proj_switch, &ICP_or_GICP_switch);

	return 0;
}