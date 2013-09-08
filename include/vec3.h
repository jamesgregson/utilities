#ifndef FLIP_VEC2_H
#define FLIP_VEC2_H

#include<cmath>
#include<iostream>

namespace geom {
    
    template< typename real >
    class vec3 {
    private:
        union {
            struct {
                real m_x;
                real m_y;
                real m_z;
            };
            real m_v[3];
        };
    public:
        vec3( ){
            set( 0.0, 0.0, 0.0 );
        }
        
        vec3( const real x, const real y, const real z ){
            set( x, y, z );
        }
        
        inline const real *ptr() const {
            return m_v;
        }
        
        inline void set( const real x, const real y, const real z ){
            m_x = x; m_y = y; m_z = z;
        }
        
        inline real dot( const vec3<real> &in ) const {
            return m_x*in.m_x + m_y*in.m_y + m_z*in.m_z;
        }
        
        inline vec3<real> cross( const vec3<real> &in ) const {
            return vec3<real>( m_y*in.m_z-m_z*in.m_y, m_z*in.m_x-m_x*in.m_z, m_x*in.m_y-m_y*in.m_x );
        }
        
        inline real length_squared() const {
            return m_x*m_x + m_y*m_y + m_z*m_z;
        }
        
        inline real length() const {
            return sqrt( m_x*m_x + m_y*m_y + m_z*m_z );
        }
        
        inline real normalize(){
            real len = length();
            m_x /= len; m_y /= len; m_z /= len;
            return len;
        }
        
        inline vec3<real> component_parallel_to( const vec3<real> &b ) const {
            real blen2 = b.length_squared();
            return b*dot(b)/blen2;
        }
        
        inline vec3<real> component_perpendicular_to( const vec3<real> &b ) const {
            return operator-( component_parallel_to( b ) );
        }
        
        inline real &x(){
            return m_x;
        }
        
        inline const real &x() const {
            return m_x;
        }
        
        inline real &y(){
            return m_y;
        }
        
        inline const real &y() const {
            return m_y;
        }
        
        inline real &z(){
            return m_z;
        }
        
        inline const real &z() const {
            return m_z;
        }
        
        inline real &operator[]( const int id ){
            return m_v[id];
        }
        
        inline const real &operator[]( const int id ) const {
            return m_v[id];
        }
        
        inline vec3<real> operator-() const {
            return vec3<real>( -m_x, -m_y, -m_z );
        }
        
        inline vec3<real> operator+( const vec3<real> b ) const {
            return vec3<real>( m_x+b.m_x, m_y+b.m_y, m_z+b.m_z );
        }
        
        inline vec3<real> operator-( const vec3<real> b ) const {
            return vec3<real>( m_x-b.m_x, m_y-b.m_y, m_z-b.m_z );
        }
        
        inline vec3<real> operator*( const real b ) const {
            return vec3<real>( m_x*b, m_y*b, m_z*b );
        }
        
        inline vec3<real> operator/( const real b ) const {
            return vec3<real>( m_x/b, m_y/b, m_z/b );
        }
        
        inline vec3<real> operator+=( const vec3<real> b ){
            vec3<real> tmp( m_x+b.m_x, m_y+b.m_y, m_z+b.m_z );
            *this=tmp;
            return tmp;
        }
        
        inline vec3<real> operator-=( const vec3<real> b ){
            vec3<real> tmp( m_x-b.m_x, m_y-b.m_y, m_z-b.m_z );
            *this=tmp;
            return tmp;
        }
        
        inline vec3<real> operator*=( const real b ){
            vec3<real> tmp( m_x*b, m_y*b, m_z*b );
            *this=tmp;
            return tmp;
        }
        
        inline vec3<real> operator/=( const real b ){
            vec3<real> tmp( m_x/b, m_y/b, m_z/b );
            *this=tmp;
            return tmp;
        }
    };
    
    template<typename real>
    vec3<real> operator*( const float &b, const vec3<real> &a ){
        return vec3<real>( b*a.x(), b*a.y(), b*a.z() );
    }
    
    template<typename real>
    vec3<real> operator*( const double &b, const vec3<real> &a ){
        return vec3<real>( b*a.x(), b*a.y(), b*a.z() );
    }
    
};

template<typename real>
std::ostream& operator<<( std::ostream &in, const geom::vec3<real> &val ){
    in << "[" << val.x() << ", " << val.y() << ", " << val.z() << "]";
    return in;
}

#endif