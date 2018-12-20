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
  D hValue = 0.0;
  D gValue = 0.0;
  int comeFromIdx = -1;

  SearchNode(){};
  SearchNode( D _hValue, D _gValue ): hValue( _hValue ), gValue( _gValue ){};
  SearchNode( D _hValue, D _gValue, std::shared_ptr<ETP::Triangle<P,D>> _comeFrom ): hValue( _hValue ), gValue( _gValue ), comeFrom( _comeFrom ){};
  SearchNode( D _hValue, D _gValue, std::shared_ptr<ETP::Triangle<P,D>> _comeFrom, int _comeFromIdx ): hValue( _hValue ), gValue( _gValue ), comeFrom( _comeFrom ), comeFromIdx( _comeFromIdx ){};
  ~SearchNode(){};

};//end of class

template <typename P, typename D>
struct searchEndPoint{

  std::shared_ptr<ETP::Triangle<P,D>> triangle;
  D gValue = 0.0;
  D gValueToSearchStart = 0.0;
  std::shared_ptr<ETP::searchEndPoint<P,D>> parent;
  int parentIndex;
  std::shared_ptr<ETP::Edge<P,D>> entryEdge = nullptr;

  searchEndPoint(){};
  searchEndPoint( D _gValue, std::shared_ptr<ETP::Triangle<P,D>> _triangle ): gValue( _gValue ), triangle( _triangle ){
    parent = nullptr;
  };
  searchEndPoint( D _gValue, std::shared_ptr<ETP::Triangle<P,D>> _triangle, std::shared_ptr<ETP::searchEndPoint<P,D>> _parent ): gValue( _gValue ), triangle( _triangle ), parent( _parent ){};
  searchEndPoint( D _gValue, std::shared_ptr<ETP::Triangle<P,D>> _triangle, std::shared_ptr<ETP::searchEndPoint<P,D>> _parent, int _parentIndex ): gValue( _gValue ), triangle( _triangle ), parent( _parent ), parentIndex( _parentIndex ){};
  searchEndPoint( D _gValue, std::shared_ptr<ETP::Triangle<P,D>> _triangle, std::shared_ptr<ETP::searchEndPoint<P,D>> _parent, int _parentIndex, D _gValueToSearchStart ): gValue( _gValue ), triangle( _triangle ), parent( _parent ), parentIndex( _parentIndex ), gValueToSearchStart( _gValueToSearchStart ){};
  ~searchEndPoint(){};

};//end of class

template <typename P, typename D>
class UsedList{
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> list;
  public:
    UsedList(){};

    ~UsedList(){
      for( auto t : list ){
        t->clearSearchNode();
      }
    };

    void push_back( std::shared_ptr<ETP::Triangle<P,D>> _tri ){
      list.push_back( _tri );
    }

};//end of class

template <typename P, typename D>
struct searchResult{
  std::shared_ptr<ETP::UsedList<P,D>> usedList;
  std::shared_ptr<ETP::searchEndPoint<P,D>> goalEndPoint;
  searchResult(){};
  searchResult( std::shared_ptr<ETP::searchEndPoint<P,D>> _goalEndPoint, std::shared_ptr<ETP::UsedList<P,D>> _usedList ): goalEndPoint( _goalEndPoint ), usedList( _usedList ){};
  std::vector<std::shared_ptr<ETP::Triangle<P,D>>> getFunnel(){
    std::vector<std::shared_ptr<ETP::Triangle<P,D>>> ret;
    if( goalEndPoint == nullptr ) return ret;
    std::shared_ptr<ETP::Triangle<P,D>> tCurrent = goalEndPoint->triangle;
    //std::shared_ptr<ETP::SearchNode<P,D>> searchNodeCurrent = goalEndPoint->triangle->getSearchNode();
    for(;;){
      ret.insert( ret.begin(), tCurrent );
      std::shared_ptr<ETP::SearchNode<P,D>> snCurrent = tCurrent->getSearchNode();
      if( snCurrent == nullptr ) break;
      std::shared_ptr<ETP::Triangle<P,D>> tComeFrom = snCurrent->comeFrom;
      if( tComeFrom == nullptr ) break;
      std::shared_ptr<ETP::Triangle<P,D>> intermediateNode = tCurrent->getAdjacentVal( snCurrent->comeFromIdx );
      if( intermediateNode != tComeFrom ){
        for(;;){
          ret.insert( ret.begin(), intermediateNode );
          std::shared_ptr<ETP::Triangle<P,D>> nextNode = intermediateNode->getAdjacentVal( intermediateNode->getConnectedNodeIndex( tComeFrom ) );
          if( nextNode == tComeFrom ) break;
          intermediateNode = nextNode;
        }
      }
      tCurrent = tComeFrom;
    }
    return ret;
  }
};// end of struct

}// end of namespace

#endif
