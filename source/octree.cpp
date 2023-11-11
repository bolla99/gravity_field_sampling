//
// Created by Giovanni Bollati on 05/11/23.
//

#include <octree.hpp>
#include <gravity.hpp>

template <typename T>
const T& octree<T>::get() const {
    return value;
}
template <typename T>
void octree<T>::set(T t) {
    value = t;
}
template <typename T>
const octree<T>& octree<T>::operator[](int i) const {
    return children->at(i);
}
template <typename T>
bool octree<T>::isLeaf() const {
    return children == nullptr;
}
template <typename T>
void octree<T>::execute(int times, T t, std::function<std::array<T, 8>(T)> f, std::function<bool(std::array<T, 8>)> condition) {
    this->set(t);
    if(times == 0) return;
    std::array<T, 8> next_values = f(t);
    if(condition(next_values)) {
        times--;
        this->expand();
        for(int i = 0; i < 8; i++) {
            children->at(i).execute(times, next_values[i], f, condition);
        }
    }
}
template <typename T>
int octree<T>::bytesize() const {
    int size = sizeof(T);
    if(!isLeaf()) {
        for(int i = 0; i < 8; i++) {
            size += children->at(i).bytesize();
        }
    }
    return size;
}
template <typename T>
void octree<T>::expand() {
    children = new std::array<octree, 8>();
}


template class octree<gravity::cube>;
template class octree<gravity::gravity_cube>;