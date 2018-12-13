/*
 * File Header for kd_tree.cpp:
 * 		This file contains functions for performing kdtree search on 3-D points 
 		with or without the aid of normals. 

 		Search steps:
 			transform the moving dataset
			kd-search for closest points in fixed dataset 
			return search result

 */
#define NOMINMAX
#include <limits>
#include "kd_tree.h"
#include <iostream>
// #include "registration_tools.h"
// #include "compute_transformed_points.h"

using namespace Eigen;
using namespace std;

void call_error(string msg){
	cout << "Error: " << msg << endl;
	exit(1);
}

/* find_distance:
 * 		Input: two points in Vector3ld type
 *      Return: distance between two points
 */
long double find_distance(Vector3ld point1, Vector3ld point2) {
	return (point1 - point2).norm();
}

/* insert_helper:
 * 		Input: point (to be inserted into the tree), kd-tree (can't be NULL), 
 			   level that the point should be sorted on
 		Return: None. Modify the tree in place by inserting the point into the tree
 */
void insert_helper(Vector3ld point, KDTree *T, int level) {
	// Right now the tree only works for x, y, z point
	if (level < 0 || level > 2) 
		call_error("Invalid search level");
	if (*T == NULL) {
		*T = (KDTree)malloc(sizeof(struct KDNode));
		if (*T == NULL)
			call_error("Malloc failed in insert_helper");
		(*T)->left = NULL;
		(*T)->right = NULL;
		((*T)->value)(0) = point(0);
		((*T)->value)(1) = point(1);
		((*T)->value)(2) = point(2);
	}
	else {
		if (point(level) < (*T)->value(level)) 
			insert_helper(point, &((*T)->left), (level+1) % 3);
		else 
			insert_helper(point, &((*T)->right), (level+1) % 3);
	}
}

/* insert:
 * 		Input: point (to be inserted into the tree), kd-tree (can't be NULL), 
 		Return: None. Modify the tree in place by inserting the point into the tree
 */
void insert(Vector3ld point, KDTree *T) {
	if (T == NULL)
		call_error("Invalid pointer for kd-tree in insert");
	return insert_helper(point, T, 0);
}

/* find_nearest_helper:
 * 		Input: kd-tree, point (whose closest match needs to be searched in kd-tree), 
 			   the level to search, a storage for current best found, a storage for
 			   current distance
		Requires TreeType to be KDTree or KDNormalTree
 		Return: None. Modify the found storages in place
 */
template <class TreeType>
void find_nearest_helper(TreeType T, Vector3ld target, int level, TreeType bestN, 
						 long double *bestDistance) {
	long double distance, diff;
	// If reaches the leaf of the tree, end search
	if (T == NULL)
		return;
	
	distance = find_distance(T->value, target);
	diff = (T->value)(level) - target(level);
	
	if (distance <= *bestDistance) {
		*bestDistance = distance;
		*bestN = *T;
	}
	//If find exact match, exit early
	if (!*bestDistance) 
		return;
	
	// Our kd-tree search is for x, y, z points only
	level = (level+1) % 3;
	find_nearest_helper(diff > 0 ? T->left : T->right, target, level, bestN, bestDistance);
	//If the candidate hypersphere crosses this splitting plane, look on the
    // other side of the plane by examining the other subtree.
    if (fabs(diff) >= *bestDistance) 
    	return;
    find_nearest_helper(diff > 0 ? T->right : T->left, target, level, bestN, bestDistance);
}

/* find_nearest:
 * 		Input: point (whose closest match needs to be searched in kd-tree), kd-tree
 		Return: The sub-tree whose node is the best match
 */
template <class NodeType>
NodeType* find_nearest(Vector3ld target, NodeType *T) {
	NodeType *bestN = (NodeType*)malloc(sizeof(NodeType));
	
	if (!bestN)
		call_error("Malloc failed in find_nearest");
	long double *distanceResult = (long double*)malloc(sizeof(long double));
	*distanceResult = numeric_limits<long double>::max();

	find_nearest_helper(T, target, 0, bestN, distanceResult);

	//free(distanceResult);
	return bestN;
}

/*
 * kd_search:
 		Input: target point cloud, kd-tree, inlier ratio, Xreg from last iteration
 			   to transform points 
		Return: pc = set of all closest points
				pr = set of all target points in corresponding order with pc
 				res = mean of the sum of all the distances calculated
 */

MatrixXld kd_search_targets(PointCloud *targets_p, KDTree T) {
	// #########################################################################################
	// Update this function in a way that it will return the nearest neighbour points and their indices.
	// #########################################################################################
	int numTargets = (*targets_p).cols();
	// int inlierSize = trunc(numTargets * inlierRatio);	// Round down to int
	PointCloud resultTargets = PointCloud(3, numTargets);
	MatrixXld resultMatches = MatrixXld(4, numTargets);	// First 3 rows = point, 
														// 4th row = distance 

	long double totalDistance = 0;
	PointCloud targetsNew = PointCloud(3, numTargets);

	// Transform the target points before searching
	targetsNew = *targets_p;

	if ((*targets_p).cols() != numTargets) {
		ostringstream errorString;
		errorString << "Pointcloud (" << (*targets_p).cols()<< ") doesn't match target size (" 
					<< numTargets << ")\n";
		call_error(errorString.str());
	}

	// Find numTargets closest points together with corresponding targets
	for (int count = 0; count < numTargets; count++) {
		KDTree nearestPoint = find_nearest(targetsNew.col(count), T);

		(resultMatches.col(count))(0) = (nearestPoint->value)(0);
		(resultMatches.col(count))(1) = (nearestPoint->value)(1);
		(resultMatches.col(count))(2) = (nearestPoint->value)(2);
		(resultMatches.col(count))(3) = find_distance(nearestPoint->value, targetsNew.col(count));
		resultTargets.col(count) = (*targets_p).col(count);	// We want to return the original targets
	}

	return resultMatches;
}

// This function frees the tree
void free_tree(KDTree T) {
	if (T) {
		free_tree(T->left);
		free_tree(T->right);
		free(T);
	}
	return;
}