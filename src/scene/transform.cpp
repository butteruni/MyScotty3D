
#include "transform.h"

Mat4 Transform::local_to_parent() const {
	return Mat4::translate(translation) * rotation.to_mat() * Mat4::scale(scale);
}

Mat4 Transform::parent_to_local() const {
	return Mat4::scale(1.0f / scale) * rotation.inverse().to_mat() * Mat4::translate(-translation);
}

Mat4 Transform::local_to_world() const {
	// A1T1: local_to_world
	//don't use Mat4::inverse() in your code.

	// return Mat4::I; //<-- wrong, but here so code will compile
    Mat4 T = Mat4::translate(translation);
    Mat4 R = rotation.to_mat();
    Mat4 S = Mat4::scale(scale);
    Mat4 local_matrix = T * R * S;
    if(auto parentPtr = parent.lock()) {
        return parentPtr->local_to_world() * local_matrix;
    }
    return local_matrix;
}

Mat4 Transform::world_to_local() const {
	// A1T1: world_to_local
	//don't use Mat4::inverse() in your code.
    Mat4 T_inv = Mat4::translate(-translation);
    Mat4 R_inv = rotation.to_mat().T();
    Mat4 S_inv = Mat4::scale(1.0f / scale);
	// return Mat4::I; //<-- wrong, but here so code will compile
    Mat4 local_matrix = S_inv * R_inv * T_inv; 
    if(auto parentPtr = parent.lock()) {
        return local_matrix * parentPtr->world_to_local();
    }
    return local_matrix;
}

bool operator!=(const Transform& a, const Transform& b) {
	return a.parent.lock() != b.parent.lock() || a.translation != b.translation ||
	       a.rotation != b.rotation || a.scale != b.scale;
}
