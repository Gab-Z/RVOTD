#ifndef DEF_ETP_SEARCH_NODE
#define DEF_ETP_SEARCH_NODE

#include <memory>
#include "./Triangle.h"

namespace ETP {

template <typename sP = int, typename sD = float>
class Triangle;

template <typename P, typename D>

struct SearchNode{

  std::shared_ptr<ETP::Triangle<P,D>> comeFrom = nullptr;
  D hValue = 0;
  D gValue = 0;

  SearchNode(){};
  /*
  SearchNode( const SearchNode<P,D>& _sn ){
    this = _sn;
  };
  */
  SearchNode( D _hValue, D _gValue ): hValue( _hValue ), gValue( _gValue ){};
  SearchNode( D _hValue, D _gValue, std::shared_ptr<ETP::Triangle<P,D>> _comeFrom ): hValue( _hValue ), gValue( _gValue ), comeFrom( _comeFrom ){};
  ~SearchNode(){};


};//end of class

}// end of namespace

#endif
