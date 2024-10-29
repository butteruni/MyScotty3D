#include <unordered_set>
#include "skeleton.h"
#include "test.h"
#include <iostream>

void Skeleton::Bone::compute_rotation_axes(Vec3 *x_, Vec3 *y_, Vec3 *z_) const {
	assert(x_ && y_ && z_);
	auto &x = *x_;
	auto &y = *y_;
	auto &z = *z_;

	//y axis points in the direction of extent:
	y = extent.unit();
	//if extent is too short to normalize nicely, point along the skeleton's 'y' axis:
	if (!y.valid()) {
		y = Vec3{0.0f, 1.0f, 0.0f};
	}

	//x gets skeleton's 'x' axis projected to be orthogonal to 'y':
	x = Vec3{1.0f, 0.0f, 0.0f};
	x = (x - dot(x,y) * y).unit();
	if (!x.valid()) {
		//if y perfectly aligns with skeleton's 'x' axis, x, gets skeleton's z axis:
		x = Vec3{0.0f, 0.0f, 1.0f};
		x = (x - dot(x,y) * y).unit(); //(this should do nothing)
	}

	//z computed from x,y:
	z = cross(x,y);

	//x,z rotated by roll:
	float cr = std::cos(roll / 180.0f * PI_F);
	float sr = std::sin(roll / 180.0f * PI_F);
	// x = cr * x + sr * -z;
	// z = cross(x,y);
	std::tie(x, z) = std::make_pair(cr * x + sr * -z, cr * z + sr * x);
}

std::vector< Mat4 > Skeleton::bind_pose() const {
	//A4T2a: bone-to-skeleton transformations in the bind pose
	//(the bind pose does not rotate by Bone::pose)

	std::vector< Mat4 > bind;
	bind.reserve(bones.size());

	//NOTE: bones is guaranteed to be ordered such that parents appear before child bones.
    Vec3 r = base + base_offset;
	for (auto const &bone : bones) {
		//avoid complaints about unused bone
		//placeholder -- your code should actually compute the correct transform:

        if(bone.parent == -1u) {
            bind.emplace_back(Mat4::translate(r));
            continue;
        }
        auto parentid = bone.parent;
        auto e_par = bones[parentid].extent;
        auto parent_transform = bind[parentid];
        bind.emplace_back(parent_transform * Mat4::translate(e_par));
		
	}

	assert(bind.size() == bones.size()); //should have a transform for every bone.
	return bind;
}

std::vector< Mat4 > Skeleton::current_pose() const {
    //A4T2a: bone-to-skeleton transformations in the current pose

	//Similar to bind_pose(), but takes rotation from Bone::pose into account.
	// (and translation from Skeleton::base_offset!)

	//You'll probably want to write a loop similar to bind_pose().

	//Useful functions:
	//Bone::compute_rotation_axes() will tell you what axes (in local bone space) Bone::pose should rotate around.
	//Mat4::angle_axis(angle, axis) will produce a matrix that rotates angle (in degrees) around a given axis.
    
	// return std::vector< Mat4 >(bones.size(), Mat4::I);

    Vec3 r = base + base_offset;
    std::vector<Mat4>pose;
    pose.reserve(bones.size());
    for(auto const & bone : bones) {
        Vec3 x, y ,z;
        bone.compute_rotation_axes(&x, &y, &z);
        Mat4 rotation = Mat4::angle_axis(bone.pose.z, z) *
                        Mat4::angle_axis(bone.pose.y, y) *
                        Mat4::angle_axis(bone.pose.x, x);
        if(bone.parent == -1u) {
            pose.emplace_back(Mat4::translate(r) * rotation);
            continue;
        }
        auto parentid = bone.parent;
        auto e_par = bones[parentid].extent;
        auto parent_transform = pose[parentid];
        pose.emplace_back(parent_transform * Mat4::translate(e_par) * rotation);
    }
    assert(pose.size() == bones.size());
    return pose;
}
std::vector< Vec3 > Skeleton:: gradient_in_current_pose() const {
    //A4T2b: IK gradient

    // Computes the gradient (partial derivative) of IK energy relative to each bone's Bone::pose, in the current pose.

	//The IK energy is the sum over all *enabled* handles of the squared distance from the tip of Handle::bone to Handle::target
	std::vector< Vec3 > gradient(bones.size(), Vec3{0.0f, 0.0f, 0.0f});
    std::vector< Mat4 > curTransform = current_pose();
	//TODO: loop over handles and over bones in the chain leading to the handle, accumulating gradient contributions.
	//remember bone.compute_rotation_axes() -- should be useful here, too!
    auto getEndPosInWorld = [&](size_t index) {
        return (curTransform[index] * Vec4(bones[index].extent, 1.0)).xyz(); 
    };
    for(auto &handle : handles) {
        if(!handle.enabled)
            continue;
        Vec3 h = handle.target;
        // compute the output point
        Vec3 p = getEndPosInWorld(handle.bone);
        Vec3 delta = p - h;
        // std::cout << "h: " << h << "\n";
        // std::cout << "p: " << p << "\n";
        // walk up the kinematic chain
        for(BoneIndex b = handle.bone; b != -1u && b < bones.size(); b = bones[b].parent) {

            Bone const & bone = bones[b];

            // TODO: compute [linear transform’]
            // Mat4 xf = curTransform[b];
            Mat4 xf = Mat4();
            Vec3 r = Vec3();
            if(bone.parent != -1u) {
                xf = curTransform[bone.parent] * Mat4::translate(bones[bone.parent].extent);
            }else {
                xf = Mat4::translate(base + base_offset);
            }
            r = (xf * Vec4(0.f, 0.f, 0.f, 1.f)).xyz();
            
            // std::cout << "r :" << r << "\n";
            // Vec3 q = bone.pose;
            Vec3 x;// TODO: compute bone’s x-axis in local space
            Vec3 y;// TODO: compute bone’s y-axis in local space
            Vec3 z;// TODO: compute bone’s z-axis in local space
            
            bone.compute_rotation_axes(&x, &y, &z);
            // std::cout << "x :" << x << "\n";
            // std::cout << "y :" << y << "\n";
            // std::cout << "z :" << z << "\n";
            Mat4 Rx = Mat4::angle_axis(bone.pose.x, x);
            Mat4 Ry = Mat4::angle_axis(bone.pose.y, y);
            Mat4 Rz = Mat4::angle_axis(bone.pose.z, z);
            Vec3 theta = bone.pose;
            z = (xf * Rz).cols[2].xyz();
            y = (xf * Rz * Ry).cols[1].xyz();
            x = (xf * Rz * Ry * Rx).cols[0].xyz();
            gradient[b].x += dot(cross(x, p - r), delta);
            gradient[b].y += dot(cross(y, p - r), delta);
            gradient[b].z += dot(cross(z, p - r), delta);  
        }
    }
	assert(gradient.size() == bones.size());
	return gradient;
}

bool Skeleton::solve_ik(uint32_t steps) {
	//A4T2b - gradient descent
	//check which handles are enabled
	//run `steps` iterations
	
	//call gradient_in_current_pose() to compute d loss / d pose
	//add ...
    float tau = 0.1f;
    while(steps--) {
        auto gradients = gradient_in_current_pose();
        for(size_t i = 0;i < bones.size(); ++i) {
            bones[i].pose -= tau * gradients[i]; 
        }
    }
	//if at a local minimum (e.g., gradient is near-zero), return 'true'.
	//if run through all steps, return `false`.
    auto gradients = gradient_in_current_pose();
    float res = 0.f;
    for(auto gradient: gradients)
        res += gradient.norm();
    if(res < EPS_F) {
        return true;
    }
    
	return false;
}

Vec3 Skeleton::closest_point_on_line_segment(Vec3 const &a, Vec3 const &b, Vec3 const &p) {
	//A4T3: bone weight computation (closest point helper)

    // Return the closest point to 'p' on the line segment from a to b
    
	//Efficiency note: you can do this without any sqrt's! (no .unit() or .norm() is needed!)
    Vec3 ab = b - a;
    Vec3 ap = p - a;
    
    float t = dot(ap ,ab) / dot(ab, ab);
    
    t = std::clamp(t, 0.0f, 1.0f);
    
    return a + ab * t;
}

void Skeleton::assign_bone_weights(Halfedge_Mesh *mesh_) const {
	assert(mesh_);
	auto &mesh = *mesh_;
	(void)mesh; //avoid complaints about unused mesh

	//A4T3: bone weight computation

	//visit every vertex and **set new values** in Vertex::bone_weights (don't append to old values)

	//be sure to use bone positions in the bind pose (not the current pose!)

	//you should fill in the helper closest_point_on_line_segment() before working on this function
    auto bindPose = bind_pose();
    for(auto &v: mesh.vertices) {
        Vec3 pos = v.position;
        float sum = 0;
        std::vector<Halfedge_Mesh::Vertex::Bone_Weight> w;
        w.reserve(bones.size());
        for(size_t i = 0;i < bones.size();++i) {
            Vec3 endPoint = bindPose[i].cols[3].xyz();
            Vec3 startPoint = (bindPose[i] * Mat4::translate(bones[i].extent)).cols[3].xyz();
            Vec3 nearPoint = closest_point_on_line_segment(startPoint, endPoint, pos);
            float d = (pos - nearPoint).norm();
            float w_ij = std::max(0.f, bones[i].radius - d) / bones[i].radius;
            w.emplace_back(Halfedge_Mesh::Vertex::Bone_Weight{(uint32_t)i, w_ij});
            sum += w_ij;
        }
        for(size_t i = 0;i < w.size(); ++i) {
            w[i].weight /= sum;
        }
        v.bone_weights = w;
    }
}
Indexed_Mesh Skeleton::skin(Halfedge_Mesh const &mesh, std::vector< Mat4 > const &bind, std::vector< Mat4 > const &current) {
	assert(bind.size() == current.size());


	//A4T3: linear blend skinning

	//one approach you might take is to first compute the skinned positions (at every vertex) and normals (at every corner)
	// then generate faces in the style of Indexed_Mesh::from_halfedge_mesh

	//---- step 1: figure out skinned positions ---

	std::unordered_map< Halfedge_Mesh::VertexCRef, Vec3 > skinned_positions;
	std::unordered_map< Halfedge_Mesh::HalfedgeCRef, Vec3 > skinned_normals;
	//reserve hash table space to (one hopes) avoid re-hashing:
	skinned_positions.reserve(mesh.vertices.size());
	skinned_normals.reserve(mesh.halfedges.size());

	//(you will probably want to precompute some bind-to-current transformation matrices here)
    std::vector<Mat4> B_inv;
    B_inv.reserve(bind.size());
    for(size_t i = 0;i < bind.size(); ++i) {
        B_inv.emplace_back(bind[i].inverse());
    }
	for (auto vi = mesh.vertices.begin(); vi != mesh.vertices.end(); ++vi) {
		// skinned_positions.emplace(vi, vi->position); //PLACEHOLDER! Replace with code that computes the position of the vertex according to vi->position and vi->bone_weights.
		//NOTE: vertices with empty bone_weights should remain in place.
        Vec3 new_pos = Vec3();
        Vec3 old_pos = vi->position;
        for(auto Info: vi->bone_weights) {
            new_pos += (Info.weight * current[Info.bone] * B_inv[Info.bone] * Vec4(old_pos, 0)).xyz();
        }
        skinned_positions.emplace(vi, new_pos);
		//circulate corners at this vertex:
		auto h = vi->halfedge;
		do {
			//NOTE: could skip if h->face->boundary, since such corners don't get emitted
            if(h->face->boundary)
                continue;
			// skinned_normals.emplace(h, h->corner_normal); //PLACEHOLDER! Replace with code that properly transforms the normal vector! Make sure that you normalize correctly.
            Halfedge_Mesh::VertexRef v = h->vertex;
            Mat4 M = Mat4::Zero;
            for(auto Info : v->bone_weights) {
                M += Info.weight * current[Info.bone] * B_inv[Info.bone];
            }
            Mat4 M_inv = M.inverse();
            Mat4 M_inv_T = M_inv.T();
            Vec3 old_normal = h->corner_normal;
            Vec3 new_normal = (M_inv_T * Vec4(old_normal, 0.f)).xyz();
            new_normal.normalize();
            skinned_normals.emplace(h, new_normal);
			h = h->twin->next;
		} while (h != vi->halfedge);
	}

	//---- step 2: transform into an indexed mesh ---

	//Hint: you should be able to use the code from Indexed_Mesh::from_halfedge_mesh (SplitEdges version) pretty much verbatim, you'll just need to fill in the positions and normals.

	// Indexed_Mesh result = Indexed_Mesh::from_halfedge_mesh(mesh, Indexed_Mesh::SplitEdges); //PLACEHOLDER! you'll probably want to copy the SplitEdges case from this function o'er here and modify it to use skinned_positions and skinned_normals.
	std::vector<Indexed_Mesh::Vert> verts;
	std::vector<Indexed_Mesh::Index> idxs;
    for (Halfedge_Mesh::FaceCRef f = mesh.faces.begin(); f != mesh.faces.end(); f++) {
			if (f->boundary) continue;

			//every corner gets its own copy of a vertex:
			uint32_t corners_begin = static_cast<uint32_t>(verts.size());
			Halfedge_Mesh::HalfedgeCRef h = f->halfedge;
			do {
                Indexed_Mesh::Vert vert;
				vert.pos = skinned_positions[h->vertex];
				vert.norm = skinned_normals[h];
				vert.uv = h->corner_uv;
				vert.id = f->id;
				verts.emplace_back(vert);
				h = h->next;
			} while (h != f->halfedge);
			uint32_t corners_end = static_cast<uint32_t>(verts.size());

			//divide face into a triangle fan:
			for (size_t i = corners_begin + 1; i + 1 < corners_end; i++) {
				idxs.emplace_back(corners_begin);
				idxs.emplace_back(static_cast<uint32_t>(i));
				idxs.emplace_back(static_cast<uint32_t>(i+1));
			}
		}
    return Indexed_Mesh(std::move(verts), std::move(idxs));
}

void Skeleton::for_bones(const std::function<void(Bone&)>& f) {
	for (auto& bone : bones) {
		f(bone);
	}
}


void Skeleton::erase_bone(BoneIndex bone) {
	assert(bone < bones.size());
	//update indices in bones:
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (bones[b].parent == -1U) continue;
		if (bones[b].parent == bone) {
			assert(b > bone); //topological sort!
			//keep bone tips in the same place when deleting parent bone:
			bones[b].extent += bones[bone].extent;
			bones[b].parent = bones[bone].parent;
		} else if (bones[b].parent > bone) {
			assert(b > bones[b].parent); //topological sort!
			bones[b].parent -= 1;
		}
	}
	// erase the bone
	bones.erase(bones.begin() + bone);
	//update indices in handles (and erase any handles on this bone):
	for (uint32_t h = 0; h < handles.size(); /* later */) {
		if (handles[h].bone == bone) {
			erase_handle(h);
		} else if (handles[h].bone > bone) {
			handles[h].bone -= 1;
			++h;
		} else {
			++h;
		}
	}
}

void Skeleton::erase_handle(HandleIndex handle) {
	assert(handle < handles.size());

	//nothing internally refers to handles by index so can just delete:
	handles.erase(handles.begin() + handle);
}


Skeleton::BoneIndex Skeleton::add_bone(BoneIndex parent, Vec3 extent) {
	assert(parent == -1U || parent < bones.size());
	Bone bone;
	bone.extent = extent;
	bone.parent = parent;
	//all other parameters left as default.

	//slightly unfortunate hack:
	//(to ensure increasing IDs within an editing session, but reset on load)
	std::unordered_set< uint32_t > used;
	for (auto const &b : bones) {
		used.emplace(b.channel_id);
	}
	while (used.count(next_bone_channel_id)) ++next_bone_channel_id;
	bone.channel_id = next_bone_channel_id++;

	//all other parameters left as default.

	BoneIndex index = BoneIndex(bones.size());
	bones.emplace_back(bone);

	return index;
}

Skeleton::HandleIndex Skeleton::add_handle(BoneIndex bone, Vec3 target) {
	assert(bone < bones.size());
	Handle handle;
	handle.bone = bone;
	handle.target = target;
	//all other parameters left as default.

	//slightly unfortunate hack:
	//(to ensure increasing IDs within an editing session, but reset on load)
	std::unordered_set< uint32_t > used;
	for (auto const &h : handles) {
		used.emplace(h.channel_id);
	}
	while (used.count(next_handle_channel_id)) ++next_handle_channel_id;
	handle.channel_id = next_handle_channel_id++;

	HandleIndex index = HandleIndex(handles.size());
	handles.emplace_back(handle);

	return index;
}


Skeleton Skeleton::copy() {
	//turns out that there aren't any fancy pointer data structures to fix up here.
	return *this;
}

void Skeleton::make_valid() {
	for (uint32_t b = 0; b < bones.size(); ++b) {
		if (!(bones[b].parent == -1U || bones[b].parent < b)) {
			warn("bones[%u].parent is %u, which is not < %u; setting to -1.", b, bones[b].parent, b);
			bones[b].parent = -1U;
		}
	}
	if (bones.empty() && !handles.empty()) {
		warn("Have %u handles but no bones. Deleting handles.", uint32_t(handles.size()));
		handles.clear();
	}
	for (uint32_t h = 0; h < handles.size(); ++h) {
		if (handles[h].bone >= HandleIndex(bones.size())) {
			warn("handles[%u].bone is %u, which is not < bones.size(); setting to 0.", h, handles[h].bone);
			handles[h].bone = 0;
		}
	}
}

//-------------------------------------------------

Indexed_Mesh Skinned_Mesh::bind_mesh() const {
	return Indexed_Mesh::from_halfedge_mesh(mesh, Indexed_Mesh::SplitEdges);
}

Indexed_Mesh Skinned_Mesh::posed_mesh() const {
	return Skeleton::skin(mesh, skeleton.bind_pose(), skeleton.current_pose());
}

Skinned_Mesh Skinned_Mesh::copy() {
	return Skinned_Mesh{mesh.copy(), skeleton.copy()};
}
