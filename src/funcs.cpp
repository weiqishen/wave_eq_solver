#include "funcs.h"
#include <algorithm>
// Method that searches a value in a sorted array without repeated entries and returns position
int get_index(int value, ndarray<int> &in_array)
{
    int size = in_array.get_len();
    int ju, jm, jl;

    jl = 0;
    ju = size - 1;

    //if (!is_sorted(in_array.get_ptr(), in_array.get_ptr() + size))
    //    Fatal_Error("array not sorted");

    if (value == in_array(0))
        return 0;
    else if (value == in_array(size - 1))
        return size - 1;

    while (ju - jl > 1)
    {
        jm = (ju + jl)/2;
        if (value >= in_array(jm))
            jl = jm;
        else
            ju = jm;
    }

    if (value == in_array(jl))
        return jl;
    else
        return -1;
}