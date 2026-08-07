# Third-Party Notices

This repository includes third-party code used by the hypergeometric-function
support in `common/hyp2f1.m` and `common/make_hyp2f1/`.

## Gauss Hypergeometric Function MATLAB/MEX Support

- Source: Siyi Deng, "Gauss hypergeometric function," MATLAB Central File
  Exchange, Version 1.0.0.0, published 11 Oct 2013.
- URL: https://www.mathworks.com/matlabcentral/fileexchange/43865-gauss-hypergeometric-function
- Included files include `common/hyp2f1.m`, `common/make_hyp2f1/mexhyp2f1.c`,
  and `common/make_hyp2f1/make_hyp2f1.m`.

## Cephes Math Library Components

Several C source files in `common/make_hyp2f1/` are based on Cephes Math Library
code by Stephen L. Moshier. The original headers in those files include the
applicable copyright notices. Files with Cephes notices include:

- `gamma.c`
- `hyp2f1.c`
- `mconf.h`
- `polevl.c`
- `psi.c`
- `round.c`

Keep the original source headers intact when modifying or redistributing these
files.
