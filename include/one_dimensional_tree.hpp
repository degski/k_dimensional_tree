
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
struct point1 {

    using value_type = Type;

    value_type x, y;

    constexpr point1 ( ) noexcept : x ( std::numeric_limits<Type>::quiet_NaN ( ) ) {}
    point1 ( point1 const & ) noexcept = default;
    point1 ( point1 && ) noexcept      = default;
    point1 ( value_type && x_ ) noexcept : x{ std::move ( x_ ) } {}

    [[maybe_unused]] point1 & operator= ( point1 const & ) noexcept = default;
    [[maybe_unused]] point1 & operator= ( point1 && ) noexcept = default;

    [[nodiscard]] bool operator== ( point1 const & p_ ) const noexcept { return x == p_.x; }
    [[nodiscard]] bool operator!= ( point1 const & p_ ) const noexcept { not operator== ( p_ ); }

    [[maybe_unused]] point1 & operator+= ( point1 const & p_ ) noexcept {
        x += p_.x;
        return *this;
    }
    [[maybe_unused]] point1 & operator-= ( point1 const & p_ ) noexcept {
        x -= p_.x;
        return *this;
    }

    // For debugging.

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, point1 const & p_ ) noexcept {
        if ( not std::isnan ( p_.x ) )
            out_ << '<' << p_.x << '>';
        else
            out_ << "<*>";
        return out_;
    }
};

using point1f = point1<float>;
using point1d = point1<double>;

// Implicit KD full binary tree of dimension 2.
template<typename Type, typename Point = point1<Type>, typename TagType = dynamic_tag, std::size_t MaxStaticSize = 0>
struct one_dimensional_tree {

    using value_type = Point;
    using base_type  = Type;
    using dist_type =
        std::conditional_t<std::is_floating_point_v<base_type>, base_type, detail::signed_double_width_integer<base_type>>;
    using pointer         = value_type *;
    using reference       = value_type &;
    using const_pointer   = value_type const *;
    using const_reference = value_type const &;

    using container_type = TagType;
    using container      = std::conditional_t<std::is_same_v<container_type, static_tag>, std::static_vector<Point, MaxStaticSize>,
                                         std::vector<Point>>;

    using iterator       = typename container::iterator;
    using const_iterator = typename container::const_iterator;

    private:
    [[nodiscard]] pointer left ( pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] pointer right ( pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer left ( const_pointer const p_ ) const noexcept { return ( p_ + 1 ) + ( p_ - m_data.data ( ) ); }
    [[nodiscard]] const_pointer right ( const_pointer const p_ ) const noexcept { return ( p_ + 2 ) + ( p_ - m_data.data ( ) ); }

    [[nodiscard]] bool is_leaf ( const_pointer const p_ ) const noexcept {
        return m_leaf_start < p_ or std::isnan ( left ( p_ )->x );
    }

    template<typename RandomIt>
    void kd_construct_x ( pointer const p_, RandomIt const first_, RandomIt const last_ ) noexcept {
        RandomIt median = detail::median_it ( first_, last_ );
        std::nth_element ( first_, median, last_, [] ( value_type const & a, value_type const & b ) { return a.x < b.x; } );
        *p_ = *median;
        if ( first_ != median ) {
            kd_construct_x ( left ( p_ ), first_, median );
            if ( ++median != last_ )
                kd_construct_x ( right ( p_ ), median, last_ );
        }
    }

    void nn_search_x ( const_pointer const p_ ) const noexcept {
        dist_type d = one_dimensional_tree::distance_squared ( *p_, m_to );
        if ( d < m_min_distance_squared ) {
            m_min_distance_squared = d;
            m_point                = p_;
        }
        if ( is_leaf ( p_ ) )
            return;
        if ( ( d = p_->x - m_to.x ) > dist_type{ 0 } ) {
            nn_search_x ( left ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_x ( right ( p_ ) );
        }
        else {
            nn_search_x ( right ( p_ ) );
            if ( ( ( d * d ) < m_min_distance_squared ) )
                nn_search_x ( left ( p_ ) );
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
    void ( one_dimensional_tree::*nn_search ) ( const_pointer const ) const noexcept;

    // These mutable types are class global result types.
    mutable const_pointer m_point = nullptr;
    mutable value_type m_to;
    mutable dist_type m_min_distance_squared = std::numeric_limits<dist_type>::max ( );

    public:
    one_dimensional_tree ( ) noexcept {}
    one_dimensional_tree ( one_dimensional_tree const & ) = delete;
    one_dimensional_tree ( one_dimensional_tree && rhs_ ) noexcept :
        m_data{ std::move ( rhs_.m_data ) }, m_leaf_start{ rhs_.m_leaf_start }, nn_search{ rhs_.nn_search } {}

    one_dimensional_tree ( std::initializer_list<value_type> ) noexcept = delete;

    template<typename ForwardIt>
    one_dimensional_tree ( ForwardIt first_, ForwardIt last_ ) noexcept {
        initialize ( first_, last_ );
    }

    template<typename ForwardIt>
    void initialize ( ForwardIt const first_, ForwardIt const last_ ) noexcept {
        auto const n = std::distance ( first_, last_ );
        if ( n ) {
            if ( n > detail::linear_bound ) {
                m_data.resize ( capacity ( static_cast<std::size_t> ( n ) ) );
                m_leaf_start = detail::median_ptr ( m_data.data ( ), m_data.size ( ) );
                kd_construct_x ( m_data.data ( ), first_, last_ );
                nn_search = &one_dimensional_tree::nn_search_x;
            }
            else {
                if constexpr ( std::is_same_v<container_type, dynamic_tag> )
                    m_data.reserve ( n );
                std::copy ( first_, last_, std::back_inserter ( m_data ) );
                nn_search = &one_dimensional_tree::nn_search_linear;
            }
        }
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

    one_dimensional_tree & operator= ( one_dimensional_tree const & ) = delete;
    one_dimensional_tree & operator                                   = ( one_dimensional_tree && rhs_ ) noexcept {
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
        return ( static_cast<dist_type> ( p1_.x ) - static_cast<dist_type> ( p2_.x ) ) *
               ( static_cast<dist_type> ( p1_.x ) - static_cast<dist_type> ( p2_.x ) );
    }

    // For debugging.

    template<typename Stream>
    [[maybe_unused]] friend Stream & operator<< ( Stream & out_, one_dimensional_tree const & tree_ ) noexcept {
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
