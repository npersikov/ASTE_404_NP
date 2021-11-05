#ifndef _OUTPUT_H
#define _OUTPUT_H

#include "World.h"
#include "Species.h"

// namespace is used to collect all output functions, similar to Java static class
namespace Output {
	void saveVTI(World &world, Species &species);
}
#endif
