// CEED Interface Stub for MFEM WebAssembly Build
// This file provides empty declarations when MFEM_USE_CEED is not defined

#ifndef MFEM_CEED_INTERFACE_OPERATOR_HPP
#define MFEM_CEED_INTERFACE_OPERATOR_HPP

#include "../../../config/config.hpp"

#ifndef MFEM_USE_CEED

namespace mfem
{

namespace ceed
{

// Empty stub classes when CEED is not available
class Operator {};

} // namespace ceed

} // namespace mfem

#endif // MFEM_USE_CEED not defined

#endif // MFEM_CEED_INTERFACE_OPERATOR_HPP
