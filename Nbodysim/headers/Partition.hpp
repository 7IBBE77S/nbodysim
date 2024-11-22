//UNUSED (may incorporate in the future)
#pragma once
#include <vector>
#include <functional>


// TODO:
//figure out a way to combine with quadtree for optimization 
template<typename T>
class Partition {
public:
    // generic partition function that returns dividing index
    static size_t partition(std::vector<T>& data, size_t start, size_t end, 
                          const std::function<bool(const T&)>& predicate) {
        if (start >= end) return start;
        
        size_t low = start;
        size_t high = end - 1;
        
        while (low < high) {
            while (low < high && predicate(data[low])) {
                low++;
            }
            
            while (low < high && !predicate(data[high])) {
                high--;
            }
            
            if (low < high) {
                std::swap(data[low], data[high]);
            }
        }
        
        return predicate(data[low]) ? low + 1 : low;
    }

    //UNUSED
    // specialized version for spatial partitioning
    template<typename Bounds>
    static size_t spatial_partition(std::vector<T>& data, size_t start, size_t end,
                                  const Bounds& bounds, size_t axis) {
        return partition(data, start, end,
            [&](const T& item) {
                return axis == 0 ? 
                    item.pos.x < bounds.center.x : 
                    item.pos.y < bounds.center.y;
            });
    }

    // multi-way partition (ie for quadtree)
    static std::array<size_t, 4> quad_partition(std::vector<T>& data, 
                                              size_t start, size_t end,
                                              const Vec2& center) {
        std::array<size_t, 4> indices;
        
        // split on x axis
        size_t mid_x = partition(data, start, end,
            [&](const T& item) { return item.pos.x < center.x; });
            
        // then split each half on y axis
        indices[0] = partition(data, start, mid_x,
            [&](const T& item) { return item.pos.y < center.y; });
            
        indices[1] = mid_x;
        
        indices[2] = partition(data, mid_x, end,
            [&](const T& item) { return item.pos.y < center.y; });
            
        indices[3] = end;
        
        return indices;
    }
};


// #pragma once
// #include <vector>
// #include "Body.hpp"
// #include "Vec2.hpp"

// class Partition {
// public:
//     // Simplified quad partition matching Rust's efficient approach
//     static std::array<size_t, 5> quad_partition(
//         std::vector<Body>& bodies,
//         size_t start, 
//         size_t end,
//         const Vec2& center) {
        
//         std::array<size_t, 5> split = {start, 0, 0, 0, end};

//         // First split by y (matching Rust implementation)
//         auto predicate_y = [&](const Body& body) { 
//             return body.pos.y < center.y; 
//         };
//         split[2] = partition(bodies, split[0], split[4], predicate_y);

//         // Then split both halves by x
//         auto predicate_x = [&](const Body& body) { 
//             return body.pos.x < center.x; 
//         };
//         split[1] = partition(bodies, split[0], split[2], predicate_x);
//         split[3] = partition(bodies, split[2], split[4], predicate_x);

//         return split;
//     }

// private:
//     template<typename Pred>
//     static size_t partition(std::vector<Body>& bodies, size_t start, size_t end, Pred pred) {
//         size_t l = start;
//         size_t r = end - 1;
        
//         while(l <= r) {
//             while(l <= r && pred(bodies[l])) l++;
//             while(l < r && !pred(bodies[r])) r--;
//             if(l >= r) return l;
//             std::swap(bodies[l], bodies[r]);
//             l++; r--;
//         }
//         return l;
//     }
// };