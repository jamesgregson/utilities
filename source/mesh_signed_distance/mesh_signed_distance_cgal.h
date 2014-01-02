#ifndef MESH_SIGNED_DISTANCE_CGAL_H
#define MESH_SIGNED_DISTANCE_CGAL_H

#include<vector>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_triangle_primitive.h>

template< typename real, typename real3 >
class mesh_signed_distance_cgal {
private:
    typedef CGAL::Simple_cartesian<double> K;
    typedef typename K::FT FT;
    typedef typename K::Ray_3 Ray;
    typedef typename K::Line_3 Line;
    typedef typename K::Point_3 Point;
    typedef typename K::Triangle_3 Triangle;
    typedef typename std::list<Triangle>::iterator Iterator;
    typedef CGAL::AABB_triangle_primitive<K,Iterator> Primitive;
    typedef CGAL::AABB_traits<K,Primitive> AABB_triangle_traits;
    typedef CGAL::AABB_tree<AABB_triangle_traits> Tree;

    Tree *m_tree;
public:
    mesh_signed_distance_cgal( const std::vector<real3> &coords, const std::vector<int> &tris ){
        std::list<Triangle> triangles;
        for( int i=0; i<tris.size(); i+=3 ){
            Point a( coords[tris[i+0]][0], coords[tris[i+0]][1], coords[tris[i+0]][2] );
            Point b( coords[tris[i+1]][0], coords[tris[i+1]][1], coords[tris[i+1]][2] );
            Point c( coords[tris[i+2]][0], coords[tris[i+2]][1], coords[tris[i+2]][2] );
            triangles.push_back( Triangle( a, b, c ) );
        }
        m_tree = new Tree( triangles.begin(), triangles.end() );
    }
    
    inline real operator()( const real3 &p ) const {
        Point ro( p[0], p[1], p[2] );
        Point rd( p[0]+1e6, p[1], p[2] );
        Ray ray( ro, rd );
        int n = m_tree->number_of_intersected_primitives( ray );
        FT sd = sqrt(m_tree->squared_distance( ro ));
        return (n%2 == 0) ? sd : -sd;
    }
};

#endif
