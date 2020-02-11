
// MIT License
//
// Copyright (c) 2018, 2019, 2020 degski
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

#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdlib>

#include <type_traits>

#include <sax/stl.hpp>

namespace sax {
namespace detail {

constexpr std::size_t linear_bound = 64u;

template<typename>
struct signed_double_width_integer {};
template<>
struct signed_double_width_integer<std::uint8_t> {
    using type = std::int16_t;
};
template<>
struct signed_double_width_integer<std::uint16_t> {
    using type = std::int32_t;
};
template<>
struct signed_double_width_integer<std::uint32_t> {
    using type = std::int64_t;
};
template<>
struct signed_double_width_integer<std::uint64_t> {
    using type = std::int64_t;
}; // Inshallah.
template<>
struct signed_double_width_integer<std::int8_t> {
    using type = std::int16_t;
};
template<>
struct signed_double_width_integer<std::int16_t> {
    using type = std::int32_t;
};
template<>
struct signed_double_width_integer<std::int32_t> {
    using type = std::int64_t;
};
template<>
struct signed_double_width_integer<std::int64_t> {
    using type = std::int64_t;
}; // Inshallah.

template<typename TagType>
struct same_sized_int {
    using type = std::make_signed_t<TagType>;
};
template<>
struct same_sized_int<float> {
    using type = std::int32_t;
};
template<>
struct same_sized_int<double> {
    using type = std::int64_t;
};

template<std::size_t StdArraySize>
[[nodiscard]] constexpr std::size_t array_size ( ) noexcept {
    return StdArraySize > linear_bound ? sax::next_power_2 ( StdArraySize + 1 ) - 1 : StdArraySize;
}

template<typename ForwardIt>
[[nodiscard]] ForwardIt median ( ForwardIt const first_, ForwardIt const last_ ) noexcept {
    return std::next ( first_, std::distance ( first_, last_ ) / 2 );
}

template<std::size_t S>
struct message { // needs fixing.

    template<typename... Args>
    message ( Args... ) {
        static_assert ( not( 2 == S or 3 == S ), "2 or 3 dimensions only" );
    }
};
} // namespace detail

struct vector_tag {};
struct array_tag {};

} // namespace sax
