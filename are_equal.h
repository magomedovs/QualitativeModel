#ifndef ARE_EQUAL
#define ARE_EQUAL

template <typename T>
bool are_equal(const T &a, const T &b, const double tol = 1e-10)
{
  if (std::abs( a - b ) < tol)
  {
    return true;
  }
  else
  {
    return false;
  }
}

#endif