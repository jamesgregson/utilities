#ifndef SIGNED_DIST_FUNC_H
#define SIGNED_DIST_FUNC_H

#include<vector>

class mesh_signed_distance {
private:
    void *tree;
public:
    mesh_signed_distance( const std::vector<float> &coords, const std::vector<int> &tris );
    mesh_signed_distance( const std::vector<double> &coords, const std::vector<int> &tris );
    double operator()( const float *in ) const;
    double operator()( const double *in ) const;
};

#endif