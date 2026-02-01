// CEED Mass Integrator Stub for MFEM WebAssembly Build
// These are empty stub classes - CEED is not available in WASM

#ifndef MFEM_CEED_MASS_HPP
#define MFEM_CEED_MASS_HPP

#include "../../../bilininteg.hpp"

namespace mfem
{
namespace ceed
{

// Stub classes - they should never be instantiated since DeviceCanUseCeed() returns false

class PAMassIntegrator : public Operator
{
public:
   PAMassIntegrator(const FiniteElementSpace &, const IntegrationRule &, Coefficient *)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MFMassIntegrator : public Operator
{
public:
   MFMassIntegrator(const FiniteElementSpace &, const IntegrationRule &, Coefficient *)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MixedPAMassIntegrator : public Operator
{
public:
   MixedPAMassIntegrator(MassIntegrator &, const FiniteElementSpace &, Coefficient *)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MixedMFMassIntegrator : public Operator
{
public:
   MixedMFMassIntegrator(MassIntegrator &, const FiniteElementSpace &, Coefficient *)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

} // namespace ceed
} // namespace mfem

#endif // MFEM_CEED_MASS_HPP
