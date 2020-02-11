# k_dimensional_tree

Fast KD-Tree (a recursive Implementation with fallback to linear search in case of a small number of points) in 2 or 3 dimensions.

The implementation uses recursion for both construction and nn-search and is strictly C++20 (and up).

    namespace sax

    struct vector_tag_t { };
    struct array_tag_t { };

    template<typename T, typename P = point2<T>, typename Type = vector_tag_t, std::size_t N = 0>
    struct two_dimensional_tree;

A `sax::point2` and `sax::point3 class` are provided, or just drop in **SFML**'s `sf::Vector2` or `sf::Vector3`. The parameter N has no effect in case of a vector, which is the default container type.

No map is implemented, but one can return an index of the node that's closest. Simply make lookup-tables for any "mapped" data.
