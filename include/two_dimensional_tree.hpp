
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

#include <limits>
#include <sax/iostream.hpp>
#include <type_traits>

#include "detail/kdt.hpp"

namespace sax {

template<typename Type>
struct point2 {

    using value_type = Type;

    value_type x, y;

    point2 ( ) noexcept : x{ std::numeric_limits<value_type>::quiet_NaN ( ) } { };
    point2 ( point2 const & ) noexcept = default;
    point2 ( point2 && ) noexcept      = default;
    point2 ( value_type && x_, value_type && y_ ) noexcept : x{ std::move ( x_ ) }, y{ std::move ( y_ ) } {}

    // template<typename SfmlVec>
    // point2 ( SfmlVec && v_ ) noexcept : x{ std::move ( v_.x ) }, y{ std::move ( v_.y ) } {}

    [[maybe_unused]] point2 & operator= ( point2 const & ) noexcept = default;
    [[maybe_unused]] point2 & operator= ( point2 && ) noexcept = default;

    [[nodiscard]] bool operator== ( point2 const & p_ ) const noexcept { return x == p_.x and y == p_.y; }
    [[nodiscard]] bool operator!= ( point2 const & p_ ) const noexcept { return x != p_.x or y != p_.y; }

    [[maybe_unused]] point2 & operator+= ( point2 const & p_ ) noexcept {
        x += p_.x;
        y += p_.y;
        return *this;
    }
    [[maybe_unused]] point2 & operator-= ( point2 const & p_ ) noexcept {
        x -= p_.x;
        y -= p_.y;
        return *this;
    }

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, point2 const & p_ ) noexcept {
        if ( not std::isnan ( p_.x ) )
            out_ << '<' << p_.x << ' ' << p_.y << '>';
        else
            out_ << "<* *>";
        return out_;
    }
};

using point2f = point2<float>;
using point2d = point2<double>;

// Implicit KD full binary tree of dimension 2.
template<typename Type, typename P = point2<Type>, typename TagType = vector_tag, std::size_t StdArraySize = 0>
struct two_dimensional_tree {

    using value_type = P;
    using base_type  = Type;
    using dist_type =
        std::conditional_t<std::is_floating_point_v<base_type>, base_type, detail::signed_double_width_integer<base_type>>;
    using pointer         = value_type *;
    using reference       = value_type &;
    using const_pointer   = value_type const *;
    using const_reference = value_type const &;

    using container_type = TagType;
    using container      = std::conditional_t<std::is_same_v<container_type, array_tag>,
                                         std::array<P, detail::array_size<StdArraySize> ( )>, std::vector<P>>;

    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:
    template<typename ForwardIt>
    [[nodiscard]] std::size_t get_dimensions_order ( ForwardIt const first_, ForwardIt const last_ ) const noexcept {
        auto const [ min_x, max_x ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.x < b.x; } );
        auto const [ min_y, max_y ] =
            std::minmax_element ( first_, last_, [] ( auto const & a, auto const & b ) { return a.y < b.y; } );
        return ( max_x->x - min_x->x ) < ( max_y->y - min_y->y );
    }

    [[nodiscard]] pointer left ( pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] pointer right ( pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer left ( const_pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer right ( const_pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }

    [[nodiscard]] bool is_leaf ( const_pointer const p_ ) const noexcept {
        return m_leaf_start < p_ or std::isnan ( left ( p_ )->x );
    }

    template<typename RandomIt>
    void kd_construct_xy ( pointer const p_, RandomIt const first_, RandomIt const last_ ) noexcept {
        RandomIt median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_yx ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_yx ( right ( p_ ), median, last_ );
        }
    }
    template<typename RandomIt>
    void kd_construct_yx ( pointer const p_, RandomIt const first_, RandomIt const last_ ) noexcept {
        RandomIt median = std::next ( first_, std::distance ( first_, last_ ) / 2 );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.y < b.y; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_xy ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_xy ( right ( p_ ), median, last_ );
        }
    }

    void nn_search_xy ( const_pointer const p_ ) const noexcept {
        dist_type d = two_dimensional_tree::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_to.x ) > dist_type{ 0 } ) {
            nn_search_yx ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yx ( right ( p_ ) );
        }
        else {
            nn_search_yx ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_yx ( left ( p_ ) );
        }
    }
    void nn_search_yx ( const_pointer const p_ ) const noexcept {
        dist_type d = two_dimensional_tree::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->y - m_to.y ) > dist_type{ 0 } ) {
            nn_search_xy ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xy ( right ( p_ ) );
        }
        else {
            nn_search_xy ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_xy ( left ( p_ ) );
        }
    }

    void nn_search_linear ( const_pointer const ) const noexcept {
        for ( value_type const & v : m_data ) {
            dist_type const d = distance_squared ( m_to, v );
            if ( d < m_min_distance_squared ) {
                m_point                = &v;
                m_min_distance_squared = d;
            }
        }
    }

    container m_data;
    const_pointer m_leaf_start = nullptr;
    void ( two_dimensional_tree::*nn_search ) ( const_pointer const ) const noexcept;

    // These mutable types are class global result types.
    mutable const_pointer m_point = nullptr;
    mutable value_type m_to;
    mutable dist_type m_min_distance_squared = std::numeric_limits<dist_type>::max ( );

    public:
    two_dimensional_tree ( ) noexcept {}
    two_dimensional_tree ( two_dimensional_tree const & ) = delete;
    two_dimensional_tree ( two_dimensional_tree && rhs_ ) noexcept :
        m_data{ std::move ( rhs_.m_data ) }, m_leaf_start{ rhs_.m_leaf_start }, nn_search{ rhs_.nn_search } {}

    two_dimensional_tree ( std::initializer_list<value_type> il_ ) noexcept {
        if constexpr ( std::is_same_v<container_type, array_tag> )
            assert ( il_.size ( ) == StdArraySize );
        if ( il_.size ( ) ) {
            if ( il_.size ( ) > detail::linear_bound ) {
                if constexpr ( std::is_same_v<container_type, array_tag> )
                    std::fill ( std::begin ( m_data ) + m_data.size ( ) / 2 - 1, std::end ( m_data ), value_type{ } );
                else
                    m_data.resize ( capacity<std::size_t> ( il_.size ( ) ) );
                m_leaf_start = m_data.data ( ) + m_data.size ( ) / 2 - 1;
                container points;
                points.reserve ( il_.size ( ) );
                std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( points ) );
                switch ( get_dimensions_order ( std::begin ( il_ ), std::end ( il_ ) ) ) {
                    case 0:
                        kd_construct_xy ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &two_dimensional_tree::nn_search_xy;
                        break;
                    case 1:
                        kd_construct_yx ( m_data.data ( ), std::begin ( points ), std::end ( points ) );
                        nn_search = &two_dimensional_tree::nn_search_yx;
                        break;
                }
            }
            else {
                if constexpr ( std::is_same_v<container_type, array_tag> ) {
                    std::copy_n ( std::begin ( il_ ), il_.size ( ), std::begin ( m_data ) );
                }
                else {
                    m_data.reserve ( il_.size ( ) );
                    std::copy ( std::begin ( il_ ), std::end ( il_ ), std::back_inserter ( m_data ) );
                }
                nn_search = &two_dimensional_tree::nn_search_linear;
            }
        }
    }

    template<typename ForwardIt>
    two_dimensional_tree ( ForwardIt first_, ForwardIt last_ ) noexcept {
        initialize ( first_, last_ );
    }

    // Returns (constexpr) the size of the std::array, or the class template parameter StdArraySize ( = 0).
    [[nodiscard]] static constexpr std::size_t size ( ) noexcept {
        if constexpr ( std::is_same_v<container_type, array_tag> )
            return detail::array_size<StdArraySize> ( );
        else
            return StdArraySize;
    }

    [[nodiscard]] iterator begin ( ) noexcept { return m_data.begin ( ); }
    [[nodiscard]] const_iterator begin ( ) const noexcept { return m_data.cbegin ( ); }
    [[nodiscard]] const_iterator cbegin ( ) const noexcept { return m_data.cbegin ( ); }

    [[nodiscard]] iterator end ( ) noexcept { return m_data.end ( ); }
    [[nodiscard]] const_iterator end ( ) const noexcept { return m_data.cend ( ); }
    [[nodiscard]] const_iterator cend ( ) const noexcept { return m_data.cend ( ); }

    [[nodiscard]] const_reference root ( ) const noexcept { return m_data.front ( ); }

    [[nodiscard]] bool is_valid ( iterator it_ ) noexcept { return not std::isnan ( it_->x ); }
    [[nodiscard]] bool is_valid ( const_iterator it_ ) const noexcept { return not std::isnan ( it_->x ); }
    [[nodiscard]] bool is_not_valid ( iterator it_ ) noexcept { return std::isnan ( it_->x ); }
    [[nodiscard]] bool is_not_valid ( const_iterator it_ ) const noexcept { return std::isnan ( it_->x ); }

    [[nodiscard]] static bool is_valid ( const_reference value_type_ ) noexcept { return not std::isnan ( value_type_.x ); }
    [[nodiscard]] static bool is_not_valid ( const_reference value_type_ ) noexcept { return std::isnan ( value_type_.x ); }

    two_dimensional_tree & operator= ( two_dimensional_tree const & ) = delete;
    two_dimensional_tree & operator                                   = ( two_dimensional_tree && rhs_ ) noexcept {
        m_data       = std::move ( rhs_.m_data );
        m_leaf_start = rhs_.m_leaf_start;
        nn_search    = rhs_.nn_search;
        return *this;
    }

    template<typename size_type>
    [[nodiscard]] reference operator[] ( size_type const i_ ) noexcept {
        return m_data[ i_ ];
    }
    template<typename size_type>
    [[nodiscard]] const_reference operator[] ( size_type const i_ ) const noexcept {
        return m_data[ i_ ];
    }

    template<typename ForwardIt>
    void initialize ( ForwardIt const first_, ForwardIt const last_ ) noexcept {
        if constexpr ( std::is_same_v<container_type, array_tag> )
            assert ( std::distance ( first_, last_ ) == StdArraySize );
        if ( first_ < last_ ) {
            auto const n = std::distance ( first_, last_ );
            if ( n > detail::linear_bound ) {
                if constexpr ( std::is_same_v<container_type, array_tag> )
                    std::fill ( std::begin ( m_data ) + m_data.size ( ) / 2 - 1, std::end ( m_data ), value_type{ } );
                else
                    m_data.resize ( capacity<std::size_t> ( static_cast<std::size_t> ( n ) ) );
                m_leaf_start = m_data.data ( ) + m_data.size ( ) / 2 - 1;
                switch ( get_dimensions_order ( first_, last_ ) ) {
                    case 0:
                        kd_construct_xy ( m_data.data ( ), first_, last_ );
                        nn_search = &two_dimensional_tree::nn_search_xy;
                        break;
                    case 1:
                        kd_construct_yx ( m_data.data ( ), first_, last_ );
                        nn_search = &two_dimensional_tree::nn_search_yx;
                        break;
                }
            }
            else {
                if constexpr ( std::is_same_v<container_type, array_tag> ) {
                    std::copy_n ( first_, n, std::begin ( m_data ) );
                }
                else {
                    m_data.reserve ( n );
                    std::copy ( first_, last_, std::back_inserter ( m_data ) );
                }
                nn_search = &two_dimensional_tree::nn_search_linear;
            }
        }
    }

    [[nodiscard]] const_pointer nn_pointer ( value_type const & point_ ) const noexcept {
        m_to                   = point_;
        m_min_distance_squared = std::numeric_limits<dist_type>::max ( );
        ( this->*nn_search ) ( m_data.data ( ) );
        return m_point;
    }
    [[nodiscard]] sax::pair<const_pointer, base_type> nn_pointer_distance ( value_type const & point_ ) const noexcept {
        return { nn_pointer ( point_ ), static_cast<base_type> ( m_min_distance_squared ) };
    }

    [[nodiscard]] std::ptrdiff_t nn_index ( value_type const & point_ ) const noexcept {
        return nn_pointer ( point_ ) - m_data.data ( );
    }
    [[nodiscard]] sax::pair<std::ptrdiff_t, base_type> nn_index_distance ( value_type const & point_ ) const noexcept {
        return { nn_index ( point_ ), static_cast<base_type> ( m_min_distance_squared ) };
    }

    [[nodiscard]] base_type nn_distance ( value_type const & point_ ) const noexcept {
        assert ( m_min_distance_squared != std::numeric_limits<dist_type>::max ( ) );
        return static_cast<base_type> ( m_min_distance_squared );
    }

    [[nodiscard]] static constexpr dist_type distance_squared ( value_type const & p1_, value_type const & p2_ ) noexcept {
        return ( ( static_cast<dist_type> ( p1_.x ) - static_cast<dist_type> ( p2_.x ) ) *
                 ( static_cast<dist_type> ( p1_.x ) - static_cast<dist_type> ( p2_.x ) ) ) +
               ( ( static_cast<dist_type> ( p1_.y ) - static_cast<dist_type> ( p2_.y ) ) *
                 ( static_cast<dist_type> ( p1_.y ) - static_cast<dist_type> ( p2_.y ) ) );
    }

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, two_dimensional_tree const & tree_ ) noexcept {
        for ( auto const & p : tree_.m_data )
            out_ << p;
        return out_;
    }

    private:
    template<typename U>
    [[nodiscard]] static constexpr U capacity ( U const i_ ) noexcept {
        assert ( i_ > 0 );
        return i_ > detail::linear_bound ? sax::next_power_2 ( i_ + 1 ) - 1 : i_;
    }
};

} // namespace sax
