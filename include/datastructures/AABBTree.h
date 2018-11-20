#pragma once

#include "Geometry.h"

// forward declaration of AABB Tree
template<uint D, typename Precision>
class AABBTree;

template<uint D, typename Precision>
struct TreeNode {

private:

    friend AABBTree<D, Precision>;

    static constexpr int cChildren = 1u << D;
    using pNode = std::unique_ptr<TreeNode>;
    using tChildren = std::array<pNode, cChildren>;

public:

    TreeNode(const dBox<D, Precision> &bounds_, const Precision cell_size_, const bool isLeaf)
            : bounds(bounds_), cell_size(cell_size_), nChildren(isLeaf ? -1 : 0) {

//#ifndef NDEBUG
//        auto checkSize = [&]() -> bool{
//            bool leaf = true;
//            for(uint d = 0; d < D; ++d) {
//                if (bounds.high[d] - bounds.low[d] > cell_size)
//                    leaf = false;
//            }
//            return leaf;
//        };
//#endif
//
//        ASSERT((nChildren == -1) == checkSize());

        }

    TreeNode(const TreeNode &other) {

        bounds = other.bounds;
        cell_size = other.cell_size;
        nChildren = other.nChildren;

        for (int i = 0; i < other.nChildren; ++i) {
            children[i] = std::make_unique<TreeNode>(*other.children[i]);
        }

    }

	TreeNode(TreeNode&& other) {
		std::swap(bounds, other.bounds);
		std::swap(cell_size, other.cell_size);
		std::swap(nChildren, other.nChildren);
		std::swap(children, other.children);
	}

	TreeNode& operator=(const TreeNode& other) {
        bounds = other.bounds;
        cell_size = other.cell_size;
        nChildren = other.nChildren;

        for (int i = 0; i < other.nChildren; ++i) {
            children[i] = std::make_unique<TreeNode>(*other.children[i]);
        }

		return *this;
	}

	TreeNode& operator=(TreeNode&& other) {
		std::swap(bounds, other.bounds);
		std::swap(cell_size, other.cell_size);
		std::swap(nChildren, other.nChildren);
		std::swap(children, other.children);
	}

    bool contains(const dBox<D, Precision> &b) const {
        return bounds.contains(b) || bounds == b;
    }

    bool contains(const TreeNode &t) const {
        return bounds.contains(t.bounds);
    }

    bool leaf() const {
        return nChildren == -1;
    }

    bool insert(pNode &&b) {

        //box does not contain node to be inserted
        if (!bounds.contains(b->bounds)) {
            return false;
        }

        // box is identical
        if (bounds == b->bounds) {
            //we have identical boxes, if the other box is a leaf and we don't have any children yet, be a leaf
            if(nChildren <= 0 && b->nChildren == -1)
                nChildren = -1;

            return true;
        }

        //test if one of children nodes contains the new node
        for (int i = 0; i < nChildren; ++i) {
            if (children[i]->contains(b->bounds)) {

                //try to insert box into children node
                if (children[i]->insert(std::move(b))) {
                    return true;
                }

                //split the node and insert new node
                children[i]->split();
                return children[i]->insert(std::move(b));
            }
        }

        //there is space in our children array
        if (nChildren < cChildren) {
            children[nChildren] = std::move(b);
            nChildren++;
            return true;
        } else {
            //split the node and insert new node
            split();
            return insert(std::move(b));
        }
    }

    bool intersects(const dSphere<D, Precision> &c) const {

        // test if sphere intersects own bounds
        if (!bounds.intersects(c))
            return false;

        if (leaf())
            return true;

        // test children
        for (int i = 0; i < nChildren; ++i) {
            if (children[i]->intersects(c))
                return true;
        }

        // no child intersects sphere
        return false;

    }

private:

    void split() {
        tChildren old;
        for (int i = 0; i < cChildren; ++i) {
            old[i] = std::move(children[i]);
        }

//        printf("splitting (%.0f %.0f %.0f) (%.0f %.0f %.0f)\n", bounds.low[0], bounds.low[1], bounds.low[2], bounds.high[0], bounds.high[1], bounds.high[2]);
        for (int i = 0; i < cChildren; ++i) {
            dBox<D, Precision> child;
            for (uint d = 0; d < D; ++d) {

                // half side length rounded down to nearest multiple of cell size
                Precision mid = std::floor((bounds.low[d] + (bounds.high[d] - bounds.low[d]) / 2) / cell_size) * cell_size;

                child.low[d] = i & (1 << d) ? mid : bounds.low[d];
                child.high[d] = i & (1 << d) ? bounds.high[d] : mid;
            }
//            printf("\t%u: (%.0f %.0f %.0f) (%.0f %.0f %.0f)\n", i, child.low[0], child.low[1], child.low[2], child.high[0], child.high[1], child.high[2]);
            children[i] = std::make_unique<TreeNode<D, Precision>>(child, cell_size, false);
        }

        for (int i = 0; i < cChildren; ++i) {
//            printf("\told %u: (%.0f %.0f %.0f) (%.0f %.0f %.0f)\n", i, old[i]->bounds.low[0], old[i]->bounds.low[1], old[i]->bounds.low[2],
//                   old[i]->bounds.high[0], old[i]->bounds.high[1], old[i]->bounds.high[2]);
            insert(std::move(old[i]));
        }
    }

private:

    dBox<D, Precision> bounds;
    Precision cell_size;
    tChildren children;
    int nChildren;

};

template<uint D, typename Precision>
class AABBTree {

public:

    AABBTree(const dBox<D, Precision> &bounds, const Precision cell_size) : m_root(round(bounds, cell_size), cell_size, false) {}

    AABBTree(const AABBTree &other) : m_root(other.m_root) {}

    bool insert(const dBox<D, Precision> &a) {
        return m_root.insert(std::make_unique<TreeNode<D, Precision>>(a, m_root.cell_size, true));
    }

    bool intersects(const dSphere<D, Precision> &c) const {
        return m_root.intersects(c);
    }

private:
    TreeNode<D, Precision> m_root;

    static dBox<D, Precision> round(dBox<D, Precision> bounds, const Precision cell_size) {

        // round bounds to cell size
        for(uint d = 0; d < D; ++d){
            bounds.low[d] = std::floor(bounds.low[d] / cell_size) * cell_size;
            bounds.high[d] = std::ceil(bounds.high[d] / cell_size) * cell_size;
        }

        return bounds;
    };
};
