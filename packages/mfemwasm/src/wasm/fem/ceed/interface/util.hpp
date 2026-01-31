// CEED Interface Stub for MFEM WebAssembly Build
#ifndef MFEM_CEED_INTERFACE_UTIL_HPP
#define MFEM_CEED_INTERFACE_UTIL_HPP

#include "../../../config/config.hpp"

#ifndef MFEM_USE_CEED

namespace mfem
{

// Stub function when CEED is not available
inline bool DeviceCanUseCeed() { return false; }

namespace ceed
{
// Placeholder namespace
} // namespace ceed

} // namespace mfem

#endif // MFEM_USE_CEED

#endif
