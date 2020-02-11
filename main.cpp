
// MIT License
//
// Copyright (c) 2020 degski
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <array>
#include <fstream>
#include <iterator>
#include <list>
#include <map>
#include <random>
#include <sax/iostream.hpp>
#include <string>
#include <type_traits>
#include <vector>

#include <plf/plf_nanotimer.h>

#include <sax/splitmix.hpp>

#include "one_dimensional_tree.hpp"
#include "three_dimensional_tree.hpp"
#include "two_dimensional_tree.hpp"

namespace sax {

template<typename base_type, std::size_t S>
using k_dimensional_tree = typename std::conditional<
    2u == S, two_dimensional_tree<base_type>,
    typename std::conditional<3u == S, three_dimensional_tree<base_type>, detail::message<S>>::type>::type;
}

int main ( ) {

    sax::splitmix64 rng{ [] ( ) {
        std::random_device rdev;
        return ( static_cast<std::size_t> ( rdev ( ) ) << 32 ) | static_cast<std::size_t> ( rdev ( ) );
    }( ) };
    std::uniform_real_distribution<float> disy{ 0.0f, 100.0f };
    std::uniform_real_distribution<float> disx{ 0.0f, 40.0f };

    plf::nanotimer timer;

    constexpr int n = 25'000;

    // std::cout << bin_tree_size ( n ) << nl;

    std::vector<sax::point2f> points;

    for ( int i = 0; i < n; ++i )
        points.emplace_back ( disx ( rng ), disy ( rng ) );

    assert ( sax::detail::median_it ( std::begin ( points ), std::end ( points ) ) ==
             sax::median2 ( std::begin ( points ), std::end ( points ) ) );

    exit ( 0 );

    timer.start ( );

    sax::two_dimensional_tree<float> tree ( std::begin ( points ), std::end ( points ) );

    std::cout << "elapsed construction " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    // std::cout << nl << tree << nl << nl;

    /*

    sax::point2d point{ disx ( rng ), disy ( rng ) };

    auto found_impl     = tree.find_nearest_recursive ( point );
    auto found_preorder = tree.find_nearest_preorder ( point );
    auto found_linear   = tree.find_nearest_linear ( point, points );

    constexpr int cnt = 1'000'000;

    timer.start ( );
    for ( int i = 0; i < cnt; ++i ) {
        sax::point2d point{ disx ( rng ), disy ( rng ) };
        found_impl += tree.find_nearest_recursive ( { disx ( rng ), disy ( rng ) } );
    }
    std::cout << "elapsed im " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    timer.start ( );
    for ( int i = 0; i < cnt; ++i ) {
        found_preorder += tree.find_nearest_preorder ( { disx ( rng ), disy ( rng ) } );
    }
    std::cout << "elapsed po " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    timer.start ( );
    for ( int i = 0; i < cnt; ++i ) {
        // found_linear += tree.find_nearest_linear ( { disx ( rng ), disy ( rng ) }, points );
    }
    std::cout << "elapsed li " << ( std::uint64_t ) timer.get_elapsed_us ( ) << " us" << nl;

    std::cout << "nearest im " << found_impl << nl;
    std::cout << "nearest pr " << found_preorder << nl;
    std::cout << "nearest li " << found_linear << nl;

    */

    return EXIT_SUCCESS;
}

// int main ( ) { return EXIT_SUCCESS; }
