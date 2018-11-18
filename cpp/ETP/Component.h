#ifndef ETP_COMPONENT
#define ETP_COMPONENT

namespace ETP {

class Component {
  int Id;

public:
  Component();
  Component ( int _id );
  ~Component();
  void setId( int _id );
  int id();
};

}

#endif
