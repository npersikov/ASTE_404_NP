#ifndef _VEC_H
#define _VEC_H

#include <ostream>
#include <math.h>

template<typename T>
class _vec3
{
    public : 
        _vec3<T>(): d{0,0,0} {}
        _vec3<T>(T a, T b, T c) : d{a,b,c} {}
        T& operator[] (int i) {return d[i];}
        T operator[] (int i) const {return d[i];}

        friend _vec3<T> operator+(const _vec3<T>&a, const _vec3<T>&b)
        {
            return _vec3<T>(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
        }

        friend _vec3<T> operator-(const _vec3<T>&a, const _vec3<T>&b)
        {
            return _vec3<T>(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
        }

        friend T dot(const _vec3<T>&a, const _vec3<T>&b)
        {
            return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
        }

        // Does it make sense to make this function static, since it has more to
        // do with the class itself and not an instance of _vec3? I think I would
        // make a funciton that is implemented like this:
        // _vec3<T> a = _vec3<T>(1,2,3);
        // double maginitude = a.mag();
        // instead of
        // double maginitude = mag(a);
        friend double mag(const _vec3<T>&a)
        {
            return sqrt(dot(a,a));
        }

        friend std::ostream& operator<<(std::ostream &out, const _vec3<T>&a)
        {
            out << a[0] << " " << a[1] << " " << a[2]; 
            return out;
        }

    protected:
        T d[3];

};

using double3 = _vec3<double>; // use vec3 and give it a nickname(?)

#endif