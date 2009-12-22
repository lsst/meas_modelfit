define(`GENERAL_ASSIGN',
`
    /// \brief $1 assignment of arrays and array expressions.
    template <typename OtherT>
    Array const &
    operator $1(Expression<OtherT> const & expr) const {
        NDARRAY_ASSERT(expr.getShape() == this->getShape());
        indir(`$3',$1)
        return *this;
    }

    /// \brief $1 assignment of scalars.
    template <typename ScalarT>
#ifndef DOXYGEN
    typename boost::enable_if<boost::is_convertible<ScalarT,Element>, Array const &>::type
#else
    Array const &
#endif
    operator $1(ScalarT const & scalar) const {
        indir(`$2',$1)
        return *this;
    }')dnl
define(`BASIC_ASSIGN_SCALAR',`std::fill(this->begin(),this->end(),scalar);')dnl
define(`BASIC_ASSIGN_EXPR',`std::copy(expr.begin(),expr.end(),this->begin());')dnl
define(`AUGMENTED_ASSIGN_SCALAR',
`Iterator const i_end = this->end();
        for (Iterator i = this->begin(); i != i_end; ++i) (*i) $1 scalar;')dnl
define(`AUGMENTED_ASSIGN_EXPR',
`Iterator const i_end = this->end();
        typename OtherT::Iterator j = expr.begin();
        for (Iterator i = this->begin(); i != i_end; ++i, ++j) (*i) $1 (*j);')dnl
define(`BASIC_ASSIGN',`GENERAL_ASSIGN(`=',`BASIC_ASSIGN_SCALAR',`BASIC_ASSIGN_EXPR')')dnl
define(`AUGMENTED_ASSIGN',`GENERAL_ASSIGN($1,`AUGMENTED_ASSIGN_SCALAR',`AUGMENTED_ASSIGN_EXPR')')dnl
