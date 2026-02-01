// CEED Interface Stub for MFEM WebAssembly Build
// This file provides empty declarations when MFEM_USE_CEED is not defined

#ifndef MFEM_CEED_INTERFACE_OPERATOR_HPP
#define MFEM_CEED_INTERFACE_OPERATOR_HPP

#include "../../../config/config.hpp"
#include "../../../linalg/vector.hpp"

#ifndef MFEM_USE_CEED

namespace mfem
{

namespace ceed
{

// Empty stub class when CEED is not available
// Provides the interface that integrators expect
class Operator
{
public:
   virtual ~Operator() { }
   virtual void AddMult(const Vector &x, Vector &y) const { }
   virtual void GetDiagonal(Vector &d) const { }
   virtual void Mult(const Vector &x, Vector &y) const { }
};

} // namespace ceed

} // namespace mfem

#endif // MFEM_USE_CEED not defined

#endif // MFEM_CEED_INTERFACE_OPERATOR_HPP
