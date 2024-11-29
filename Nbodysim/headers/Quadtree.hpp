#pragma once
#include "Node.hpp"
#include "Quad.hpp"
#include <vector>
#include <algorithm>
#include <mutex>
class Quadtree
{
public:
    static inline const size_t ROOT = 0;
    float t_sq;
    float e_sq;
    size_t leaf_capacity;
    std::vector<Node> nodes;
    std::vector<size_t> parents;
    std::mutex tree_mutex;

    Quadtree(float theta, float epsilon, size_t num_bodies, size_t leaf_capacity = 16)
        : t_sq(theta * theta), e_sq(epsilon * epsilon), leaf_capacity(leaf_capacity)
    {
        // nodes.reserve(10000); // pre-allocate for better performance
        // parents.reserve(1000);
        size_t estimated_nodes = num_bodies * 4 / 3; // perfect quad tree node count
        nodes.reserve(estimated_nodes);
        parents.reserve(num_bodies / leaf_capacity);
    }

    void clear(Quad quad)
    {
        nodes.clear();
        parents.clear();
        nodes.push_back(Node(0, quad));
    }

    void insert(Vec2 pos, float mass)
    {
        std::lock_guard<std::mutex> lock(tree_mutex);
        size_t node = ROOT;

        while (nodes[node].is_branch())
        {
            size_t quadrant = nodes[node].data.quad.find_quadrant(pos); 
            node = nodes[node].children + quadrant;
        }

        if (nodes[node].is_empty())
        {
            nodes[node].data.pos = pos;
            nodes[node].data.mass = mass;
            return;
        }

        Vec2 existing_pos = nodes[node].data.pos;
        float existing_mass = nodes[node].data.mass;

        if (pos == existing_pos)
        {
            nodes[node].data.mass += mass;
            return;
        }

        while (true)
        {
            size_t children = nodes.size();
            nodes[node].children = children;
            parents.push_back(node);

            auto quads = nodes[node].data.quad.subdivide(); 
            for (size_t i = 0; i < 4; ++i)
            {
                nodes.emplace_back(Node(
                    (i < 3) ? children + i + 1 : nodes[node].next,
                    quads[i],
                    nodes[node].depth + 1
                    ));
            }
            size_t q1 = nodes[node].data.quad.find_quadrant(existing_pos);
            size_t q2 = nodes[node].data.quad.find_quadrant(pos);

            if (q1 == q2)
            {
                node = children + q1;
            }
            else
            {
                nodes[children + q1].data.pos = existing_pos;
                nodes[children + q1].data.mass = existing_mass;
                nodes[children + q2].data.pos = pos;
                nodes[children + q2].data.mass = mass;
                return;
            }
        }
    }

    // "Evil Fast" inverse square root using the "Quake III" algorithm
    // inline float fast_inv_sqrt(float x) const
    // {
    //     float half_x = 0.5f * x;
    //     int i = *(int *)&x;           // Interpret float bits as int
    //     i = 0x5f3759df - (i >> 1);    // Initial approximation
    //     x = *(float *)&i;             // Interpret int bits back as float
    //     x *= 1.5f - (half_x * x * x); // Newton-Raphson refinement
    //     return x;
    // }

    inline constexpr float fast_inv_sqrt(float number) const noexcept
    {
        const auto y = std::bit_cast<float>(
            0x5f3759df - (std::bit_cast<uint32_t>(number) >> 1));
        return y * (1.5f - (number * 0.5f * y * y));
    }

Vec2 acc(const Vec2& pos, const std::vector<Body>& bodies) const {
    Vec2 acc = Vec2::zero();
    size_t node = ROOT;

    while(true) {
        const Node& n = nodes[node];
        Vec2 d = n.data.pos - pos;
        float d_sq = d.mag_sq();

        if (n.data.quad.size * n.data.quad.size < d_sq * t_sq) {
            // far enough... use approximation
            if (d_sq > 0) {
                float inv_dist = fast_inv_sqrt(d_sq + e_sq);
                float inv_dist_cubed = inv_dist * inv_dist * inv_dist;
                acc += d * (n.data.mass * inv_dist_cubed);
            }

            if (n.next == 0) break;
            node = n.next;
        }
        else if (n.is_leaf()) {
            for (size_t i = n.bodies.start; i < n.bodies.end; ++i) {
                const Body& body = bodies[i];
                Vec2 r = body.pos - pos;
                float r_sq = r.mag_sq();

                if (r_sq > 0) {
                    float inv_dist = fast_inv_sqrt(r_sq + e_sq);
                    float inv_dist_cubed = inv_dist * inv_dist * inv_dist;
                    acc += r * (body.mass * inv_dist_cubed);
                }
            }

            if (n.next == 0) break;
            node = n.next;
        }
        else {
            node = n.children;
        }
    }

    return acc;
}

    void build(const std::vector<Body> &bodies)
    {
        size_t estimated_nodes = bodies.size() * 2;
        nodes.reserve(estimated_nodes);
        parents.reserve(bodies.size());

        clear(Quad::new_containing(bodies));
        for (const auto &body : bodies)
        {
            insert(body.pos, body.mass);
        }

        propagate();
    }

    void subdivide(size_t node_index, std::vector<Body> &bodies)
    {
        const Vec2 &center = nodes[node_index].data.quad.center;
        std::array<size_t, 5> split;
        split[0] = nodes[node_index].bodies.start;
        split[4] = nodes[node_index].bodies.end;

        auto predicate_y = [&center](const Body &body)
        {
            return body.pos[1] < center[1]; 
        };
        split[2] = partition_in_place(bodies, split[0], split[4], predicate_y);

        auto predicate_x = [&center](const Body &body)
        {
            return body.pos[0] < center[0];
        };

        split[1] = partition_in_place(bodies, split[0], split[2], predicate_x);
        split[3] = partition_in_place(bodies, split[2], split[4], predicate_x);

        parents.push_back(node_index);
        size_t children_start = nodes.size();
        nodes[node_index].children = children_start;

        auto quads = nodes[node_index].data.quad.subdivide();
        for (size_t i = 0; i < 4; ++i)
        {
            nodes.emplace_back(Node(
                (i < 3) ? children_start + i + 1 : nodes[node_index].next,
                quads[i]));
            nodes.back().bodies.start = split[i];
            nodes.back().bodies.end = split[i + 1];
        }
    }


private:
    template <typename Pred>
    static size_t
    partition_in_place(std::vector<Body> &bodies, size_t start, size_t end, Pred pred)
    {
        if (start >= end)
            return start;

        size_t left = start;
        size_t right = end - 1;

        while (true)
        {
            while (left <= right && pred(bodies[left]))
                left++;
            while (left < right && !pred(bodies[right]))
                right--;

            if (left >= right)
                return left;

            std::swap(bodies[left], bodies[right]);
            left++;
            right--;
        }
    }

    void propagate()
    {
        for (auto it = parents.rbegin(); it != parents.rend(); ++it)
        {
            size_t node = *it;
            size_t child = nodes[node].children;

            nodes[node].data.pos = Vec2::zero();
            nodes[node].data.mass = 0;

            for (size_t i = 0; i < 4; ++i)
            {
                const Node &child_node = nodes[child + i];
                nodes[node].data.pos += child_node.data.pos * child_node.data.mass;
                nodes[node].data.mass += child_node.data.mass;
            }

            if (nodes[node].data.mass > 0)
            {
                nodes[node].data.pos /= nodes[node].data.mass;
            }
        }
    }

};