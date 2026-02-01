// CEED Interface Stub for MFEM WebAssembly Build
#ifndef MFEM_CEED_INTERFACE_UTIL_HPP
#define MFEM_CEED_INTERFACE_UTIL_HPP

#include "../../../config/config.hpp"

#ifndef MFEM_USE_CEED

// Forward declaration
namespace mfem { class FiniteElementSpace; }

namespace mfem
{

// Stub function when CEED is not available
inline bool DeviceCanUseCeed() { return false; }

namespace ceed
{

// Stub function when CEED is not available - does nothing
inline void RemoveBasisAndRestriction(const FiniteElementSpace *fes) { }

} // namespace ceed

} // namespace mfem

#endif // MFEM_USE_CEED

#endif
