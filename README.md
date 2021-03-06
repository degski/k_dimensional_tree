# k_dimensional_tree

Fast KD-Tree (a recursive Implementation with fallback to linear search in case of a small number of points) in 1, 2 or 3 dimensions.

The implementation uses recursion for both construction and nn-search and is strictly C++20 (and up) (`std::static_vector`).


    namespace sax

    struct dynamic_tag { };
    struct static_tag { };

    using point2f = point2<float>;
    using point2d = point2<double>;

    template<typename Type, std::size_t MaxStaticSize = 0, typename Point = point2<Type>, typename TagType = dynamic_tag,
         typename Allocator = std::allocator<Point>>
    struct two_dimensional_tree;


A `sax::point1`, `sax::point2` and `sax::point3 class` are provided, one can cast from **SFML**'s `sf::Vector2` or `sf::Vector3` to the respective `pointN`. The parameter `MaxStaticSize` has no effect in case of a `dynamic_tag`, which is the default container type.

No map/node-data is implemented, but one can return an index of the node that's closest, enabling flexible lookup-tables for any 'mapped' data, while keeping the data structure as compact as possible.

Other than default construction, there is no other way to populate the tree than from range-iterators, objects in the 'source' will be permutated after the construction.

### IDEAS:

### TODO:
