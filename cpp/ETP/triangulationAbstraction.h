#ifndef DEF_ETP_TRIANGULATION_ABSTRACT
#define DEF_ETP_TRIANGULATION_ABSTRACT

#include "./Triangle.h"
#include "./Edge.h"
#include "./Component.h"
#include <stack>
#include <limits>

namespace ETP {

  void collapseUnrootedTree( std::shared_ptr<ETP::Triangle> t, std::shared_ptr<ETP::Component> c );

}

#endif
