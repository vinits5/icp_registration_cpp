/*
 * File Header:
 * This file contains functions for performing KDTree Search.
 */

#ifndef KD_TREE
#define KD_TREE

#include <eigen3/Eigen/Dense>
#include "type_defs.h"
#include <vector>

using std::string;
using std::vector;

struct KDNode;
typedef struct KDNode *KDTree;
typedef struct KDNormalNode *KDNormalTree;

struct KDNode{
	Vector3ld value;
	KDTree left;
	KDTree right;
};

struct KDNormalNode{
	Vector3ld value;
	KDNormalTree left;
	KDNormalTree right;
	int index;	// Index of the value where it was originally in the ptcldFixed
};

struct KdResult{
	PointCloud pc;
	PointCloud pr;
	long double res;
};

struct KDNormalResult{
	long double res1;  // res1 = mean of all the point distances calculated
	long double res2;  // res2 = mean of all the normal distances calculated
	PointCloud pc;
	PointCloud pr;
	PointCloud normalc;
	PointCloud normalr;
};

/* insert:
 * 		Input: point (to be inserted into the tree), kd-tree (can't be NULL), 
 		Return: None. Modify the tree in place by inserting the point into the tree
 */
void insert(Vector3ld point, KDTree *T);

/*
 * kd_search:
 		Input: target point cloud, kd-tree, inlier ratio, Xreg from last iteration
 			   to transform points 
		Return: pc = set of all closest points
				pr = set of all target points in corresponding order with pc
 				res = mean of the sum of all the distances calculated
 */
MatrixXld kd_search_targets(PointCloud *targets, KDTree T);

// This function frees the tree
void free_tree(KDTree T);

#endif
