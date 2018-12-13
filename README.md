# Iterative Closest Point Algorithm

![armadillo dataset results](https://github.com/vinits5/icp_registration_cpp/master/doc/results.jpg)

The code is an implementation of Iterative Closest Point Registration algorithm. This algorithm is used to find the rigid transformation between Sensor Point Cloud and Model Point Cloud.

### Dataset
It also contains two datasets. One of armadillo and another of car. Armadillo dataset has no noise. But Car dataset has a guassian noise.

### Codes
**icp_with_options.cpp** contains the main file to handle point clouds and call ICP method.
**ICP/ICP.h** is the class for Iterative Closest Point Algorithm.
**KDTree/kd_tree.h** is the class for KD-Trees. It has **insert** & **kd_search_targets** methods to insert elements in a KD-Tree and search a nearest neighbour for given target respectively.

### Run the code
> mkdir build
> 
> make
