/**
 * NumJS Polynomial Module
 *
 * Provides polynomial classes for multiple basis systems,
 * compatible with NumPy's numpy.polynomial module.
 *
 * Available polynomial types:
 * - Polynomial: Power series (standard form)
 * - Chebyshev: Chebyshev polynomials of the first kind
 * - Legendre: Legendre polynomials
 * - Hermite: Physicist's Hermite polynomials
 * - HermiteE: Probabilist's Hermite polynomials
 * - Laguerre: Laguerre polynomials
 */

// Utilities
export {
  PolyError,
  PolyDomainWarning,
  trimseq,
  trimcoef,
  as_series,
  getdomain,
  mapdomain,
  mapparms,
} from './polyutils.js';

// Base class
export { ABCPolyBase, maxpower, companionEigenvalues } from './_polybase.js';

// Polynomial (power series)
export {
  Polynomial,
  polyval,
  polyval2d,
  polyval3d,
  polyvander,
  polyvander2d,
  polyder,
  polyint,
  polyfit,
  polyroots,
  polycompanion,
  polyfromroots,
  polyadd,
  polysub,
  polymul,
  polydiv,
  polypow,
} from './polynomial.js';

// Chebyshev
export {
  Chebyshev,
  chebval,
  chebval2d,
  chebvander,
  chebder,
  chebint,
  chebfit,
  chebroots,
  chebcompanion,
  chebfromroots,
  chebinterpolate,
  poly2cheb,
  cheb2poly,
  chebadd,
  chebsub,
  chebmul,
  chebdiv,
  chebpow,
} from './chebyshev.js';

// Legendre
export {
  Legendre,
  legval,
  legvander,
  legder,
  legint,
  legfit,
  legroots,
  legcompanion,
  legfromroots,
  poly2leg,
  leg2poly,
  legadd,
  legsub,
  legmul,
  legdiv,
  legpow,
} from './legendre.js';

// Hermite (Physicist's)
export {
  Hermite,
  hermval,
  hermvander,
  hermder,
  hermint,
  hermfit,
  hermroots,
  hermcompanion,
  hermfromroots,
  poly2herm,
  herm2poly,
  hermadd,
  hermsub,
  hermmul,
  hermdiv,
  hermpow,
} from './hermite.js';

// HermiteE (Probabilist's)
export {
  HermiteE,
  hermeval,
  hermevander,
  hermeder,
  hermeint,
  hermefit,
  hermeroots,
  hermecompanion,
  hermefromroots,
  poly2herme,
  herme2poly,
  hermeadd,
  hermesub,
  hermemul,
  hermediv,
  hermepow,
} from './hermite_e.js';

// Laguerre
export {
  Laguerre,
  lagval,
  lagvander,
  lagder,
  lagint,
  lagfit,
  lagroots,
  lagcompanion,
  lagfromroots,
  poly2lag,
  lag2poly,
  lagadd,
  lagsub,
  lagmul,
  lagdiv,
  lagpow,
} from './laguerre.js';
