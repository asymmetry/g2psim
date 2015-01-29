## G2PSim

Simulation package for Jefferson Lab experiment E08-027 (g2p).

### Usage
Run.C is a root script written for debug and test. It can be used directly in ROOT environment.

Main.cc is written with supports of a simple set of command-line parameters as well as a fully functional configuration file. It is compiled to an executable.

### Authors:
Chao Gu (contact), Min Huang, Jie Liu, Jixie Zhang

### Versions:
* Main program: 1.8.0
* HRSTrans package: 1.3.2
* G2PPhys package: 1.5.0

### History:
* v1.8, 01/29/2015, 7f25161
  * Add geometry classes

* v1.7, 11/12/2013, 0235332
  * Support configuration file
  * Add energy loss and multiple scattering

* v1.6, 07/25/2013, 1cc22c5
  * Add BPM simulation

* v1.5, 03/26/2013, 03d82df
  * Rewrite simulation framework, add G2PProcBase
  * Separate cross section models as G2PPhys package

* v1.4, 02/20/2013, a9cdc62
  * Add target field and ray tracing functions

* v1.3, 01/25/2013, 8382c3a
  * Reorganize HRSTrans package

* v1.2, 01/18/2013, 7f95f23
  * Add focus plane coordinate transformations
  * Add ROOT support

* v1.1, 01/15/2013, c836551
  * Rewrite simulation framework, add G2PAppBase

* v1.0, 12/20/2012, 4ab06ff
  * First version developed from Jixie's TestSNAKE package
  * Add reconstruction function to use optics database

### Acknowledge:
* Authors of [ROOT](https://root.cern.ch).
* O. Hansen et al. for [Hall A analyzer](http://hallaweb.jlab.org/podd/index.html).
* J. Zhang for HRSMC and TestSNAKE package.
* J. Huang et al. for HRS Optics Optimize package.
* X. Zheng et al. for SAMC package.
* M. Friedman for electron elastic cross sections fits of 1H, 4He and 14N.
* L. Cardman for electron elastic cross sections fits of 12C.
* P. Bosted for electron inclusive inelastic cross sections [fits](https://userweb.jlab.org/~bosted/fits.html).
* Authors of QFS model, EPC model and WISER model.
* Authors of [Geant4](http://geant4.cern.ch/).
