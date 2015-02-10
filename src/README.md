## Sources for PO project

There is 3 applications in the folder :

- [po_app](po_app.cfg) : algorithm to resolve the problem, includes :
  - [main.cpp](main.cpp) : loads the mesh, defines options and launchs appropriate methods
  - [psi0.hpp](psi0.hpp) : objects to get psi0
  - [psi0.cpp](psi0.cpp) : sources for the object
  - [eigenprob.hpp](eigenprob.hpp) : objects to get the eigenmodes
  - [eigenprob.cpp](eigenprob.hpp) : sources for the object
  - [spectralproblem.hpp](spectralproblem.hpp) : objects to resolve the spectral problem
  - [spectralproblem.cpp](spectralproblem.cpp) : sources for the object
- [po_basischange](basischange.cpp) : creates the matrix to change the basis, see #17
- [po_eigenmixed](eigenmixed.cpp) : compute the mixed eigenmodes, see #19
