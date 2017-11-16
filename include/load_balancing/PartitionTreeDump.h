#pragma once

#include <iostream>
#include "PartitionTree.h"

namespace LoadBalancing
{
    template <uint D, typename T>
    void printPointCloud(const Point_Ids& pointCloud) {
        std::cout << "(" << pointCloud.size() << ")" << " {";
        for(auto id : pointCloud){
            if(dPoint<D, T>::isFinite(id))
                std::cout << id << " ";
        }
        std::cout << "}";
    }

    template <uint D, typename T>
    void dumpPartition(const PartitionTree<D, T>& tree, size_t indentation = 0) {
        for(size_t i = 0; i < indentation; ++i){
            std::cout << "|";
            if(i < indentation - 1)
                std::cout << "  ";
        }
        
        if(std::holds_alternative<typename PartitionTree<D, T>::ChildContainer>(tree.attachment)) {
            const auto& children = std::get<typename PartitionTree<D, T>::ChildContainer>(tree.attachment);
            for(const auto& child : children) {
                std::cout << "\n";
                dumpPartition(child, indentation + 1);
            }
        } else {
            const auto& ids = std::get<Point_Ids>(tree.attachment);
            printPointCloud<D, T>(ids);
        }
    }
}
    
