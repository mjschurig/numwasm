// CEED Convection Integrator Stub for MFEM WebAssembly Build
// These are empty stub classes - CEED is not available in WASM

#ifndef MFEM_CEED_CONVECTION_HPP
#define MFEM_CEED_CONVECTION_HPP

#include "../../../bilininteg.hpp"

namespace mfem
{
namespace ceed
{

// Stub classes - they should never be instantiated since DeviceCanUseCeed() returns false

class PAConvectionIntegrator : public Operator
{
public:
   PAConvectionIntegrator(const FiniteElementSpace &, const IntegrationRule &,
                          VectorCoefficient *, real_t)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MFConvectionIntegrator : public Operator
{
public:
   MFConvectionIntegrator(const FiniteElementSpace &, const IntegrationRule &,
                          VectorCoefficient *, real_t)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MixedPAConvectionIntegrator : public Operator
{
public:
   MixedPAConvectionIntegrator(ConvectionIntegrator &, const FiniteElementSpace &,
                               VectorCoefficient *, real_t)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MixedMFConvectionIntegrator : public Operator
{
public:
   MixedMFConvectionIntegrator(ConvectionIntegrator &, const FiniteElementSpace &,
                               VectorCoefficient *, real_t)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

} // namespace ceed
} // namespace mfem

#endif // MFEM_CEED_CONVECTION_HPP
