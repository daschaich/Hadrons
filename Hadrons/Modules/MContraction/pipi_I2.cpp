// by Chris Culver
// Copying and modifying as appropriate Meson.cpp
// this will compute the correlation function for pion-pion scattering
// in isospin 2.

#include <Hadrons/Modules/MContraction/pipi_I2.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MContraction;

template class Grid::Hadrons::MContraction::TPiPi<FIMPL,FIMPL,FIMPL,FIMPL>;
