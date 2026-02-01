// CEED Diffusion Integrator Stub for MFEM WebAssembly Build
// These are empty stub classes - CEED is not available in WASM

#ifndef MFEM_CEED_DIFFUSION_HPP
#define MFEM_CEED_DIFFUSION_HPP

#include "../../../bilininteg.hpp"

namespace mfem
{
namespace ceed
{

// Stub classes - they should never be instantiated since DeviceCanUseCeed() returns false

class PADiffusionIntegrator : public Operator
{
public:
   PADiffusionIntegrator(const FiniteElementSpace &, const IntegrationRule &,
                         Coefficient *, VectorCoefficient *, MatrixCoefficient *)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MFDiffusionIntegrator : public Operator
{
public:
   MFDiffusionIntegrator(const FiniteElementSpace &, const IntegrationRule &,
                         Coefficient *, VectorCoefficient *, MatrixCoefficient *)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MixedPADiffusionIntegrator : public Operator
{
public:
   MixedPADiffusionIntegrator(DiffusionIntegrator &, const FiniteElementSpace &)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

class MixedMFDiffusionIntegrator : public Operator
{
public:
   MixedMFDiffusionIntegrator(DiffusionIntegrator &, const FiniteElementSpace &)
   { MFEM_ABORT("CEED not available"); }
   void Mult(const Vector &, Vector &) const override { }
   void GetDiagonal(Vector &) const { }
};

} // namespace ceed
} // namespace mfem

#endif // MFEM_CEED_DIFFUSION_HPP
