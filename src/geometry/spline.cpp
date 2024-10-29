
#include "../geometry/spline.h"

template<typename T> T Spline<T>::at(float time) const {

	// A4T1b: Evaluate a Catumull-Rom spline

	// Given a time, find the nearest positions & tangent values
	// defined by the control point map.

	// Transform them for use with cubic_unit_spline

	// Be wary of edge cases! What if time is before the first knot,
	// before the second knot, etc...

    if(!any()) {
        return T();
    } 
    if(knots.size() == 1u) {
        return knots.begin()->second;
    }
    float start_T = knots.begin()->first;
    if(time < start_T) {
        return knots.begin()->second;
    }
    float end_T = std::prev(knots.end())->first;
    if(time > end_T) {
        return std::prev(knots.end())->second;
    }
    if(has(time)) {
        return knots.at(time);
    }

    std::vector<float>left_times, right_times;
    for(auto [t, val] : knots) {
        if(t < time) {
            left_times.push_back(t);
        }else {
            right_times.push_back(t);
        }
    }

    float t0, t1, t2 , t3;
    T p0, p1, p2, p3;
    T m0, m1;
	t1 = left_times.back();
    t2 = right_times[0];
    p1 = knots.at(t1);
    p2 = knots.at(t2);

    if(left_times.size() < 2) {
        t0 = t1 - (t2 - t1); 
        p0 = p1 - (p2 - p1);
    }else {
        t0 = left_times[left_times.size() - 2];
        p0 = knots.at(t0);
    }
    if (right_times.size() <= 1){
        t3 = t2 + (t2 - t1);
        p3 = p2 + (p2 - p1);
    }
    else {
        t3 = right_times[1];
        p3 = knots.at(t3);
    }
    m0 = (p2 - p0) / (t2 - t0);
    m1 = (p3 - p1) / (t3 - t1);
    
    float t_ratio = (time - t1) / (t2 - t1);

    return cubic_unit_spline(t_ratio, p1, p2, m0, m1);
}

template<typename T>
T Spline<T>::cubic_unit_spline(float time, const T& position0, const T& position1,
                               const T& tangent0, const T& tangent1) {

	// A4T1a: Hermite Curve over the unit interval

	// Given time in [0,1] compute the cubic spline coefficients and use them to compute
	// the interpolated value at time 'time' based on the positions & tangents

	// Note that Spline is parameterized on type T, which allows us to create splines over
	// any type that supports the * and + operators.
    float time_2 = time * time;
    float time_3 = time_2 * time;
    float h_00 = 2 * time_3 - 3 * time_2 + 1;
    float h_10 = time_3 - 2 * time_2 + time;
    float h_01 = -2 * time_3 + 3 * time_2;
    float h_11 = time_3 - time_2;

	return T(h_00 * position0 + h_10 * tangent0 + h_01 * position1 + h_11 * tangent1);
}

template class Spline<float>;
template class Spline<double>;
template class Spline<Vec4>;
template class Spline<Vec3>;
template class Spline<Vec2>;
template class Spline<Mat4>;
template class Spline<Spectrum>;
