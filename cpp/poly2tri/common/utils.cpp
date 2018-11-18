#include "./utils.h"

namespace p2t {

Orientation Orient2d(const Point& pa, const Point& pb, const Point& pc)
{
  double detleft = (pa.x - pc.x) * (pb.y - pc.y);
  double detright = (pa.y - pc.y) * (pb.x - pc.x);
  double val = detleft - detright;
  if (val > -EPSILON && val < EPSILON) {
    return COLLINEAR;
  } else if (val > 0) {
    return CCW;
  }
  return CW;
}

/*
bool InScanArea(Point& pa, Point& pb, Point& pc, Point& pd)
{
  double pdx = pd.x;
  double pdy = pd.y;
  double adx = pa.x - pdx;
  double ady = pa.y - pdy;
  double bdx = pb.x - pdx;
  double bdy = pb.y - pdy;

  double adxbdy = adx * bdy;
  double bdxady = bdx * ady;
  double oabd = adxbdy - bdxady;

  if (oabd <= EPSILON) {
    return false;
  }

  double cdx = pc.x - pdx;
  double cdy = pc.y - pdy;

  double cdxady = cdx * ady;
  double adxcdy = adx * cdy;
  double ocad = cdxady - adxcdy;

  if (ocad <= EPSILON) {
    return false;
  }

  return true;
}

*/

bool InScanArea(const Point& pa, const Point& pb, const Point& pc, const Point& pd)
{
  double oadb = (pa.x - pb.x)*(pd.y - pb.y) - (pd.x - pb.x)*(pa.y - pb.y);
  if (oadb >= -EPSILON) {
    return false;
  }

  double oadc = (pa.x - pc.x)*(pd.y - pc.y) - (pd.x - pc.x)*(pa.y - pc.y);
  if (oadc <= EPSILON) {
    return false;
  }
  return true;
}
}
