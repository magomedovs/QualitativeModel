#ifndef LINSPACE
#define LINSPACE

#include <array>
#include <vector>
#include <cmath>

/* Fill the std::array with values from range [interval_begin, interval_end]. */
template <typename T, size_t SIZE>
void linspace(T interval_begin, T interval_end, std::array<T, SIZE> &arr) 
{
    const T step = (interval_end - interval_begin) / static_cast<T>(SIZE - 1);
    for (size_t i = 0; i < SIZE; ++i) 
    {
        arr[i] = interval_begin + i * step;
    }
}
/* Fill the std::vector with values from range [interval_begin, interval_end]. */
template <typename T>
void linspace(T interval_begin, T interval_end, std::vector<T> &arr) 
{
    const size_t SIZE = arr.size();
    const T step = (interval_end - interval_begin) / static_cast<T>(SIZE - 1);
    for (size_t i = 0; i < SIZE; ++i) 
    {
        arr[i] = interval_begin + i * step;
    }
}
/* Fill the array of size SIZE with values from range [interval_begin, interval_end]. */
template <typename T>
void linspace(T interval_begin, T interval_end, T arr[], size_t SIZE) 
{
    const T step = (interval_end - interval_begin) / static_cast<T>(SIZE - 1);
    for (size_t i = 0; i < SIZE; ++i) 
    {
        arr[i] = interval_begin + i * step;
    }
}

#endif
