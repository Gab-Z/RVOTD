#include "./Component.h"
namespace ETP {

  Component::Component(){}
  Component::Component ( int _id ):Id( _id){}
  Component::~Component(){}
  void Component::setId( int _id ){
    Id = _id;
  }
  int Component::id(){
    return Id;
  }
}
