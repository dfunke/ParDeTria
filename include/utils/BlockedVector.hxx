#pragma once

#include <vector>
#include <memory>

template<typename... VArgs>
class BlockedVector {

  public:
    template<class Container, class IT>
    struct iterator : public std::iterator<std::bidirectional_iterator_tag, typename Container::value_type> {

    public:
        iterator(Container & container, std::size_t block, IT it)
        : m_blockedContainer(container),
          m_block(block),
          m_it(it) { }

        iterator &operator++() {
          ++m_it;
          if(m_it == m_blockedContainer.m_blocks[m_block]->end()
             && m_block < m_blockedContainer.m_blocks.size() - 1){
              m_it = m_blockedContainer.m_blocks[++m_block]->begin();
          }
          return *this;
        }

        iterator operator++(int) {
            auto tmp = *this;
            ++(*this);
            return tmp;
        }

        iterator &operator--() {
            if(m_it != m_blockedContainer.m_blocks[m_block]->begin())
              --m_it;
            else if(m_block > 0){
              m_it = --m_blockedContainer.m_blocks[--m_block]->end();
            }
            return *this;
        }

        iterator operator--(int) {
            auto tmp = *this;
            --(*this);
            return tmp;
        }

        bool operator==(iterator j) const {
          return m_block == j.m_block && m_it == j.m_it;
        }

        bool operator!=(iterator j) const { return !(*this == j); }

        const auto &operator*() const { return *m_it; }
        const auto *operator->() const { return m_it; }

        auto &operator*() { return *m_it; }
        auto *operator->() { return m_it; }

    protected:
        Container &  m_blockedContainer;
        std::size_t m_block;
        IT  m_it;
    };

  private:
    typedef BlockedVector<VArgs...> self;
    typedef std::vector<VArgs...> block;
    typedef std::unique_ptr<block> block_ptr;

    typedef typename block::value_type value_type;
    typedef typename block::iterator internal_iterator;
    typedef typename block::const_iterator internal_const_iterator;


  public:
    template <typename... CTorArgs>
    BlockedVector(CTorArgs... args) {
      new_block(args...);
    }

    void append(self && other){
      m_blocks.insert(m_blocks.end(),
          std::make_move_iterator(other.m_blocks.begin()),
          std::make_move_iterator(other.m_blocks.end()));
    }

    template<typename... Args>
    void emplace_back(Args... args){
      m_blocks.back()->emplace_back(args...);
    }

    template <typename... CTorArgs>
    void new_block(CTorArgs... args){
      m_blocks.emplace_back(std::make_unique<block>(args...));
    }

    auto begin() { return iterator<self, internal_iterator>(*this, 0, m_blocks.front()->begin()); }
    auto end() { return iterator<self, internal_iterator>(*this, m_blocks.size() - 1, m_blocks.back()->end()); }

    auto begin() const { return iterator<const self, internal_const_iterator>(*this, 0, m_blocks.front()->begin()); }
    auto end() const { return iterator<const self, internal_const_iterator>(*this, m_blocks.size() - 1, m_blocks.back()->end()); }

  private:
    std::vector<block_ptr> m_blocks;

};
