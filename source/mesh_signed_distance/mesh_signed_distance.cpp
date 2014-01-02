#include<vec3.h>
#include<mesh_signed_distance_cgal.h>

#include<mesh_signed_distance.h>


mesh_signed_distance::mesh_signed_distance( const std::vector<float> &coords, const std::vector<int> &tris ){
    typedef geom::vec3<double> real3;
    std::vector< real3 > tcoords;
    for( int i=0; i<coords.size(); i+=3 ){
        tcoords.push_back( real3( coords[i+0], coords[i+1], coords[i+2] ) );
    }
    this->tree = new mesh_signed_distance_cgal< double, geom::vec3<double> >( tcoords, tris );
}
mesh_signed_distance::mesh_signed_distance( const std::vector<double> &coords, const std::vector<int> &tris ){
    typedef geom::vec3<double> real3;
    std::vector< real3 > tcoords;
    for( int i=0; i<coords.size(); i+=3 ){
        tcoords.push_back( real3( coords[i+0], coords[i+1], coords[i+2] ) );
    }
    this->tree = new mesh_signed_distance_cgal< double, geom::vec3<double> >( tcoords, tris );
}
double mesh_signed_distance::operator()( const float *in ) const {
    geom::vec3<double> p( in[0], in[1], in[2] );
    return ((mesh_signed_distance_cgal<double,geom::vec3<double> >*)tree)->operator()( p );
}

double mesh_signed_distance::operator()( const double *in ) const {
    geom::vec3<double> p( in[0], in[1], in[2] );
    return ((mesh_signed_distance_cgal<double,geom::vec3<double> >*)tree)->operator()( p );
}