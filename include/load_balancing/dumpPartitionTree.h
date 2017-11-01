#pragma once
#include <ostream>
#include "PartitionTree.h"

template <typename T, size_t D>
void printPointCloud(std::ostream& stream, const Point_Ids<T, D>& pointCloud) {
    std::cout << " {";
    for(const auto& point : pointCloud){
        stream << point;
    }
    std::cout << "}";
}


template <typename T, size_t D>
void dumpPartitionTree(std::ostream& stream, const PartitionTree<D, T>& tree, size_t indentation = 0) {
    for(size_t i = 0; i < indentation; ++i){
        stream << "|";
        if(i < indentation - 1)
            stream << "  ";
    }
    
    //std::cout << "-> (" << partition.size() << ")";
    
    if(auto children = std::get_if<PartitionTree<D, T>::ChildContainer>(&partition.attachment)) {
        for(const auto& child : *children)
        {
            stream << "\n";
            dumpPartitionTree(stream, child, indentation + 1);
        }
    } else {
        printPointCloud<T, D>(stream, std::get<PartitionTree<D, T>::ChildContainer>(partition.attachment));
    }
}
