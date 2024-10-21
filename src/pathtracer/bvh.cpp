
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"
#include <functional>
#include <stack>

namespace PT {

struct BVHBuildData {
	BVHBuildData(size_t start, size_t range, size_t dst) : start(start), range(range), node(dst) {
	}
	size_t start; ///< start index into the primitive array
	size_t range; ///< range of index into the primitive array
	size_t node;  ///< address to update
};

struct SAHBucketData {
	BBox bb;          ///< bbox of all primitives
	size_t num_prims; ///< number of primitives in the bucket
};
template <typename Primitive>
size_t BVH<Primitive>::build_recursive(size_t max_leaf_size, size_t start, size_t size) {
    size_t node_idx = new_node();
    BBox bound = BBox();
    for(size_t i = start; i < start + size; ++i) {
        bound.enclose(primitives.at(i).bbox());
    }

    nodes.at(node_idx).size = size;
    nodes.at(node_idx).start = start;
    nodes.at(node_idx).bbox = bound;
    bound.max += Vec3(1.f, 1.f, 1.f);
    if(size < max_leaf_size) {
        nodes.at(node_idx).l = node_idx;
        nodes.at(node_idx).r = node_idx;
        return node_idx;
    }
    constexpr size_t bucket_size = 8;
    std::vector<Primitive> Partitions[2];
    float minCost = FLT_MAX;
    for(size_t dim = 0; dim < 3; ++dim) {

        float bucket_len = bound.dia()[dim] / bucket_size;
        float left = bound.min[dim];
        std::vector<std::vector<Primitive>> buckets(bucket_size);
        std::vector<BBox> bboxOfBuckets(bucket_size);
        for(size_t i = start; i < start + size; ++i) {
            BBox box = primitives.at(i).bbox();
            float center = box.center()[dim];
            size_t bucket_id = std::floor((center - left) / bucket_len);
            buckets.at(bucket_id).push_back(std::move(primitives[i]));
            bboxOfBuckets.at(bucket_id).enclose(box);
        }
        
        std::vector<BBox> leftBox(bucket_size);
        std::vector<size_t> leftAccNum(bucket_size);
        leftAccNum[0] = 0;
        std::vector<BBox> rightBox(bucket_size);
        std::vector<size_t> rightAccNum(bucket_size);
        rightAccNum[0] = 0;
        for(size_t i = 1;i < bucket_size;++i) {
            leftBox[i] = leftBox[i - 1];
            leftBox[i].enclose(bboxOfBuckets[i - 1]);
            leftAccNum[i] = leftAccNum[i - 1] + buckets[i - 1].size();
            
            rightBox[i] = rightBox[i - 1];
            rightBox[i].enclose(bboxOfBuckets[bucket_size - i]);
            rightAccNum[i] = rightAccNum[i - 1] + buckets[bucket_size - i].size();
        }

        size_t minLeftIdx = 0;
        float dimBestCost = FLT_MAX;
        for(size_t leftIdx = 1;leftIdx < bucket_size; ++leftIdx) {
            size_t rightIdx = bucket_size - leftIdx;
            float leftS = leftBox.at(leftIdx).surface_area();
            float rightS = rightBox.at(rightIdx).surface_area();
            float leftN = leftAccNum.at(leftIdx);
            float rightN = rightAccNum.at(rightIdx);
        
            float curCost = leftS * leftN + rightS * rightN;
            if(curCost < dimBestCost) {
                dimBestCost = curCost;
                minLeftIdx = leftIdx;
            }
        }

        if(dimBestCost < minCost) {
            minCost = dimBestCost;
            Partitions[0].clear();
            Partitions[1].clear();

            for(size_t i = 0;i < minLeftIdx; ++i) {
                for(size_t j = 0;j < buckets[i].size(); j++)
                    Partitions[0].push_back(std::move(buckets[i][j]));
            }
            for(size_t i = minLeftIdx;i < bucket_size;++i) {
                for(size_t j = 0; j < buckets[i].size(); j++)
                    Partitions[1].push_back(std::move(buckets[i][j]));
            }
        }
    }
    if(Partitions[0].size() == size || Partitions[1].size() == size) {
        int left_size = size / 2;
        size_t l = build_recursive(max_leaf_size, start, left_size);
        size_t r = build_recursive(max_leaf_size, start + left_size, size - left_size);
        nodes.at(node_idx).l = l;
        nodes.at(node_idx).r = r;
        return node_idx;
    }

    for(size_t i = 0;i < Partitions[0].size();++i) {
        primitives.at(i + start) = std::move(Partitions[0][i]);
    }
    for(size_t i = 0;i < Partitions[1].size(); ++i) {
        primitives.at(i + start + Partitions[0].size()) =  std::move(Partitions[1][i]);
    }
    size_t l = build_recursive(max_leaf_size, start, Partitions[0].size());
    size_t r = build_recursive(max_leaf_size, start + Partitions[0].size(), Partitions[1].size());
    nodes.at(node_idx).l = l;
    nodes.at(node_idx).r = r;
    return node_idx;
}
template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    nodes.clear();
    primitives = std::move(prims);
    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.
	//TODO
    build_recursive(max_leaf_size, 0, primitives.size());
 
}

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	//TODO: replace this code with a more efficient traversal:
    Trace ret;
    for(const Primitive& prim : primitives) {
        Trace hit = prim.hit(ray);
        ret = Trace::min(ret, hit);
    }
    return ret;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
