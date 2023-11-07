//
// Created by Giovanni Bollati on 05/11/23.
//

#ifndef GL_TEST_PROJECT_OCTREE_HPP
#define GL_TEST_PROJECT_OCTREE_HPP

#include <vector>
#include <glm/glm.hpp>
#include <functional>
#include <iostream>

// AF: leaf octree -> value = ?, children = nullptr
//     non leaf octree -> value = ?, children = std::array<octree>*
template <typename T>
class octree {
public:
    // default costructor
    octree() : value() {
        children = nullptr;
    }
    // construct an octree without children and given value
    explicit octree(T t) : value(t) {
        children = nullptr;
    }
    // must delete children, which is dynamically allocated by expand() calls
    ~octree() {
        delete children;
    }

    // getter
    const T& get() const;

    // setter
    void set(T t);

    // return const octree& -> octree children can't be modified
    const octree<T>& operator[](int i) const;

    // check if octree is a leaf octree
    [[nodiscard]] bool isLeaf() const;

    //this is the only function that can recursively expand an octree
    // "execute" set this octree value to parameter "t", then execute "f" on parameter "t" to get 8 next values;
    // if "condition"(8  next values) is true then the octree expands and reiterate with the new 8 values if
    // pre condition: this->isLeaf() == true
    void execute(int times, T t, std::function<std::array<T, 8>(T)> f, std::function<bool(std::array<T, 8>)> condition);

    // return octree size in bytes
    [[nodiscard]] int bytesize() const;

private:
    T value;
    std::array<octree<T>, 8>* children;

    // expand the octree: allocate children array and initialize with default octree
    void expand();
};

#endif //GL_TEST_PROJECT_OCTREE_HPP
