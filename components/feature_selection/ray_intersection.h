#ifndef RAY_INTERSECTION_H
#define RAY_INTERSECTION_H

#include "triangle_mesh_type.h"

template<typename T>
class Ray
{
public:

    typedef typename vcg::Point3<T> PointType;

    PointType ori;
    PointType dir;
    PointType dir_inv;
    
    Ray(const PointType& o, const PointType& d): ori(o), dir(d.normalized())
    {
        dir_inv =  PointType(1. / dir[0], 1. / dir[1], 1. / dir[2]);
    } 

    PointType operator()(T t) const { return  ori + dir * t; }
};


template<typename T>
struct Intersection
{        
    T dist;
    bool happened;  
    vcg::Point3<T> pos;
    void* ptr;

    Intersection(): happened(false), pos(), dist(std::numeric_limits<T>::max()), ptr(nullptr) {    pos.SetZero();  }
};

#endif