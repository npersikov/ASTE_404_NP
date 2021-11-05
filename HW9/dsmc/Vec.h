#ifndef _VEC_H
#define _VEC_H

#include <iostream>

template<typename T>
struct _vec3 {
    _vec3(T a, T b, T c) :
        d{ a,b,c } {}
    _vec3() : d{ 0,0,0 } {}

    T operator[] (int i) const { return d[i]; }
    T& operator[] (int i) { return d[i]; }

    //vec3*scalar 
    friend _vec3<T> operator*(_vec3<T> a, double s) {
        return _vec3<T>{a[0]*s, a[1]*s, a[2]*s };
    }

    //vec3/scalar
    friend _vec3<T> operator/(_vec3<T> a, double s) {
        return _vec3<T>{a[0]/s, a[1]/s, a[2]/s };
    }

    //scalar*vec3
    friend _vec3<T> operator*(double s,_vec3<T> a) {
        return _vec3<T>{a[0]*s, a[1]*s, a[2]*s };
    }

    //vec3+vec3
    friend _vec3<T> operator+(_vec3<T> a, _vec3<T> b) {
        return _vec3<T>{a[0]+b[0], a[1]+b[1], a[2]+b[2] };
    }

    //vec3-vec3
	friend _vec3<T> operator-(_vec3<T> a, _vec3<T> b) {
		return _vec3<T>{a[0]-b[0], a[1]-b[1], a[2]-b[2] };
	}

	//vec3 - scalar
	template<typename S>
	friend _vec3<T> operator-(_vec3<T> a, S s) {
		return _vec3<T>{a[0]-s, a[1]-s, a[2]-s };
	}


    //vec3*vec3
    template<typename S>
    friend _vec3<T> operator*(_vec3<S> a, _vec3<T> b) {
	   return _vec3<T>{a[0]*b[0], a[1]*b[1], a[2]*b[2] };
	}

    //vec3/vec3, division by zero risk, returns type of A
    template<typename S>
    friend _vec3<T> operator/(_vec3<T> a, _vec3<S> b) {
	   return _vec3<T>{a[0]/b[0], a[1]/b[1], a[2]/b[2] };
	}

    friend std::ostream& operator<<(std::ostream &out, _vec3<T> a) {
	   return out<<a[0]<<" "<<a[1]<<" "<<a[2];
	}

    friend T mag(_vec3<T> a) {return sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);}

    //compound addition
    void operator+=(_vec3<T> a) {
        d[0] += a[0]; d[1] += a[1]; d[2] += a[2];
    }

    //compound subtraction
    void operator-=(_vec3<T> a) {
        d[0] -= a[0]; d[1] -= a[1]; d[2] -= a[2];
    }

    //sets all values to s
    void operator=(T s) {d[0]=s;d[1]=s;d[2]=s;}

	friend _vec3<T> unit(const _vec3<T> &a) {
		double m = mag(a);
		return {a[0]/m,a[1]/m,a[2]/m};
	}

	//cross product
	friend _vec3<T> cross(const _vec3<T> &a, const _vec3<T> &b) {
		return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
	}

	//dot product of two vectors
	friend T dot(const _vec3<T> &a, const _vec3<T> &b) {
		return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];}

protected:
    T d[3];
};

using double3 = _vec3<double>;
using int3 = _vec3<int>;



#endif

