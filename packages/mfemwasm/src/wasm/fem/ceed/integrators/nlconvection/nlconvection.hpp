// CEED NL Convection Integrator Stub for MFEM WebAssembly Build
// These are empty stub classes - CEED is not available in WASM

#ifndef MFEM_CEED_NLCONVECTION_HPP
#define MFEM_CEED_NLCONVECTION_HPP

#include "../../../nonlininteg.hpp"

namespace mfem
{
namespace ceed
{

// Stub classes - they should never be instantiated since DeviceCanUseCeed() returns false

class PAVectorConvectionNLIntegrator : public Operator
{
public:
   PAVectorConvectionNLIntegrator(const FiniteElementSpace &, const IntegrationRule &,
                                   Coefficient *)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
};

class MFVectorConvectionNLIntegrator : public Operator
{
public:
   MFVectorConvectionNLIntegrator(const FiniteElementSpace &, const IntegrationRule &,
                                   Coefficient *)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
};

} // namespace ceed
} // namespace mfem

#endif // MFEM_CEED_NLCONVECTION_HPP
