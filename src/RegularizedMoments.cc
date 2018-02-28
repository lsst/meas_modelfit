template <typename iter>
Moments::Moments(iter begin, iter end) {
    // The first entry in the vector will be zeroth moment
    // vec[0]
    zero = *begin;
    // The next two entries will be the first moment, i.e. x, y
    // vec[1], vec[2]
    one << *(begin++), *(begin++);
    // The next three entries correspond to the second moment
    // i.e. xx, xy, yy: vec[3], vec[4], vec[5]
    two << *(begin++), *(begin++), *begin, *(begin++);
    assert(begin++ == end)
}

Moments::Moments(std::vector<double> moments): Moments(moments.begin(), moments.end()){
}

Moments::Moments(double zero, FirstMoment first, SecondMoment second): zero(zero), one(first), two(second){
}
