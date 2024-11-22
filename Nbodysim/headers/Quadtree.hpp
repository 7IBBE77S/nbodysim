#pragma once
#include "Node.hpp"
#include "Partition.hpp"
#include <vector>

class Quadtree
{
public:
     static const size_t ROOT = 0;
    float t_sq;
    float e_sq;
    size_t leaf_capacity;
    std::vector<Node> nodes;
    std::vector<size_t> parents;

    Quadtree(float theta, float epsilon, size_t leaf_capacity = 16)
        : t_sq(theta * theta), e_sq(epsilon * epsilon), leaf_capacity(leaf_capacity) 
    {
        nodes.reserve(1000);
        parents.reserve(1000);
    }

    void clear(Quad quad)
    {
        nodes.clear();
        parents.clear();
        nodes.push_back(Node(0, quad));
    }

    void insert(Vec2 pos, float mass)
    {
        size_t node = ROOT;
        while (nodes[node].is_branch())
        {
            size_t quadrant = nodes[node].quad.find_quadrant(pos);
            node = nodes[node].children + quadrant;
        }

        if (nodes[node].is_empty())
        {
            nodes[node].pos = pos;
            nodes[node].mass = mass;
            return;
        }

        Vec2 p = nodes[node].pos;
        float m = nodes[node].mass;
        if (pos == p)
        {
            nodes[node].mass += mass;
            return;
        }

        while (true)
        {
            size_t children = subdivide(node);
            size_t q1 = nodes[node].quad.find_quadrant(p);
            size_t q2 = nodes[node].quad.find_quadrant(pos);

            if (q1 == q2)
            {
                node = children + q1;
            }
            else
            {
                size_t n1 = children + q1;
                size_t n2 = children + q2;
                nodes[n1].pos = p;
                nodes[n1].mass = m;
                nodes[n2].pos = pos;
                nodes[n2].mass = mass;
                return;
            }
        }
    }

    // Vec2 acc(Vec2 pos) const {
    //     Vec2 acc = Vec2::zero();
    //     size_t node = ROOT;

    //     while (node < nodes.size()) {
    //         if (node != ROOT && nodes[node].is_empty()) {
    //             node = nodes[node].next;
    //             continue;
    //         }

    //         Vec2 d = nodes[node].pos - pos;
    //         float d_sq = d.mag_sq();

    //         if (node != ROOT && (nodes[node].quad.size * nodes[node].quad.size < d_sq * t_sq)) {
    //             if (d_sq > 0) {
    //                 float force = nodes[node].mass / (d_sq * std::sqrt(d_sq + e_sq));
    //                 acc += d * force;
    //             }
    //             node = nodes[node].next;
    //         } else if (nodes[node].is_leaf()) {
    //             if (d_sq > 0) {
    //                 float force = nodes[node].mass / (d_sq * std::sqrt(d_sq + e_sq));
    //                 acc += d * force;
    //             }
    //             node = nodes[node].next;
    //         } else {
    //             node = nodes[node].children;
    //         }
    //     }
    //     return acc;
    // }


    Vec2 acc(Vec2 pos) const
    {
        Vec2 acc = Vec2::zero();
        size_t node = ROOT;

        while (true)
        {
            const Node &n = nodes[node];
            Vec2 d = n.pos - pos;
            float d_sq = d.mag_sq() + e_sq;

            if (n.is_leaf() || n.quad.size * n.quad.size < d_sq * t_sq)
            {
                if (d_sq > 0)
                {
                    float force = n.mass / (d_sq * std::sqrt(d_sq));
                    acc += d * force;
                }

                if (n.next == 0)
                    break;
                node = n.next;
            }
            else
            {
                node = n.children;
            }
        }

        return acc;
    }

    void build(const std::vector<Body> &bodies)
    {
        nodes.reserve(bodies.size() * 2);
        parents.reserve(bodies.size());
        clear(Quad::new_containing(bodies));

        for (size_t i = 0; i < bodies.size(); i++)
        {
            insert(bodies[i].pos, bodies[i].mass);
        }
        propagate();
    }

    size_t subdivide(size_t node) {  
        parents.push_back(node);
        size_t children = nodes.size();
        nodes[node].children = children;

        std::array<size_t, 4> nexts = {
            children + 1,
            children + 2,
            children + 3,
            nodes[node].next
        };

        std::array<Quad, 4> quads = nodes[node].quad.subdivide();
        for (size_t i = 0; i < 4; ++i) {
            nodes.push_back(Node(nexts[i], quads[i]));
        }

        return children;
    }

private:
    //UNUSED
    std::array<size_t, 4> partition_bodies(std::vector<Body> &bodies, size_t start, size_t end, const Quad &quad)
    {
        std::array<size_t, 4> indices;

        size_t mid_x = Partition<Body>::partition(bodies, start, end,
                                                  [&quad](const Body &body)
                                                  {
                                                      return body.pos.x < quad.center.x;
                                                  });

        indices[0] = start;
        indices[1] = Partition<Body>::partition(bodies, start, mid_x,
                                                [&quad](const Body &body)
                                                {
                                                    return body.pos.y < quad.center.y;
                                                });
        indices[2] = Partition<Body>::partition(bodies, mid_x, end,
                                                [&quad](const Body &body)
                                                {
                                                    return body.pos.y < quad.center.y;
                                                });
        indices[3] = end;

        return indices;
    }

    void propagate()
    {
        for (auto it = parents.rbegin(); it != parents.rend(); ++it)
        {
            size_t node = *it;
            Vec2 com = Vec2::zero();
            float total_mass = 0.0f;

            for (size_t i = 0; i < 4; i++)
            {
                size_t child = nodes[node].children + i;
                if (nodes[child].mass > 0)
                {
                    com += nodes[child].pos * nodes[child].mass;
                    total_mass += nodes[child].mass;
                }
            }

            if (total_mass > 0)
            {
                nodes[node].pos = com / total_mass;
                nodes[node].mass = total_mass;
            }
        }
    }
};