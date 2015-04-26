#include "CoCoA/ULong2Long.H"
using CoCoA::ULong2Long;

#include <limits>
using std::numeric_limits;

// Check whether ULong2Long works at four negative values:
//   MinLong, MinLong+1, (MinLong+1)/2, and -1
// Exits with 0 iff the fn works correctly at all four points.
// A non-zero exit code indicates which check failed, but this
// information is not currently used anywhere.

// This program is compiled and run by cpp-flags-ulong2long.sh
// trying the various possible values for CoCoA_ULONG2LONG.

int main()
{
  volatile unsigned long ul;
  long l = numeric_limits<long>::min();
  ul = l;
  if (ULong2Long(ul) != l) return 1;

  l = l+1;
  ul = l;
  if (ULong2Long(ul) != l) return 2;

  l = l/2;
  ul = l;
  if (ULong2Long(ul) != l) return 3;

  l = -1;
  ul = l;
  if (ULong2Long(ul) != l) return 4;

  return 0;
}
