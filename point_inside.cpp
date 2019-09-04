/***************************************************************************
 *   Copyright (C) 2016 by Саша Миленковић                                 *
 *   sasa.milenkovic.xyz@gmail.com *
 *                                                                         *
 *   The code in is_inside_sm function is free software; you can           *
 *   redistribute it and/or modify it under the terms of the GNU General   *
 *   Public License as published by the Free Software Foundation; either   *
 *   version 2 of the License, or (at your option) any later version.      * *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *   ( http://www.gnu.org/licenses/gpl-3.0.en.html )                       *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#include <stdio.h>
#include <time.h>
#include <stdlib.h>

//===================================================================
struct Point {
  Point(double x0, double y0) : x(x0), y(y0) {}
  double x;
  double y;
};

//===================================================================
const Point test_polygon2[39] = {
  Point(0, 1),   Point(1, 1),   Point(1, 2),   Point(2, 2),   Point(2, 3),
  Point(3, 3),   Point(3, 2),   Point(4, 0),   Point(5, 9),   Point(6, 0),
  Point(7, 2),   Point(8, 0),   Point(8, -2),  Point(7, -3),  Point(6, -2),
  Point(5, -2),  Point(4, -2),  Point(3, -1),  Point(2, -2),  Point(1, -2),
  Point(0, -3),  Point(-2, -3), Point(-3, -4), Point(-4, -3), Point(-5, -3),
  Point(-5, -2), Point(-4, -2), Point(-4, .5), Point(-5, .5), Point(-5, 1),
  Point(-4, 1),  Point(-4, 2),  Point(-4, 4),  Point(-3, 4),  Point(-3, 2),
  Point(-2, 2),  Point(-2, 0),  Point(-1, 1),  Point(0, 1)
};

//===================================================================
// isLeft(): tests if a point is Left|On|Right of an infinite line.
//    Input:  three points P0, P1, and P2
//    Return: >0 for P2 left of the line through P0 and P1
//            =0 for P2  on the line
//            <0 for P2  right of the line
inline double isLeft(Point P0, Point P1, Point P2) {
  return ((P1.x - P0.x) * (P2.y - P0.y) - (P2.x - P0.x) * (P1.y - P0.y));
}

//===================================================================
// cn_PnPoly(): crossing number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  0 = outside, 1 = inside
// This code is patterned after [Franklin, 2000]
// Taken from: http://geomalgorithms.com/a03-_inclusion.html
int cn_PnPoly(Point P, const Point *V, int n) {
  int cn = 0; // the  crossing number counter

  // loop through all edges of the polygon
  for (int i = 0; i < n; i++) {                     // edge from V[i]  to V[i+1]
    if (((V[i].y <= P.y) && (V[i + 1].y > P.y))     // an upward crossing
        || ((V[i].y > P.y) && (V[i + 1].y <= P.y))) // a downward crossing
    {
      // compute  the actual edge-ray intersect x-coordinate
      double vt = (double)(P.y - V[i].y) / (V[i + 1].y - V[i].y);
      if (P.x < V[i].x + vt * (V[i + 1].x - V[i].x)) // P.x < intersect
        ++cn; // a valid crossing of y=P.y right of P.x
    }
  }

  return (cn & 1); // 0 if even (out), and 1 if  odd (in)
}

//===================================================================
// wn_WindingNumber(): winding number test for a point in a polygon
//      Input:   P = a point,
//               V[] = vertex points of a polygon V[n+1] with V[n]=V[0]
//      Return:  wn = the winding number (=0 only when P is outside)
// Taken from: http://geomalgorithms.com/a03-_inclusion.html
int wn_WindingNumber(Point P, const Point *V, int n) {
  int wn = 0; // the  winding number counter

  // loop through all edges of the polygon
  for (int i = 0; i < n; i++) {            // edge from V[i] to  V[i+1]
    if (V[i].y <= P.y) {                   // start y <= P.y
      if (V[i + 1].y > P.y)                // an upward crossing
        if (isLeft(V[i], V[i + 1], P) > 0) // P left of  edge
          ++wn;                            // need to have a valid up intersect
    } else {                               // start y > P.y (no test needed)
      if (V[i + 1].y <= P.y)               // a downward crossing
        if (isLeft(V[i], V[i + 1], P) < 0) // P right of  edge
          --wn; // need to have a valid down intersect
    }
  }

  return wn;
}

//===================================================================
/*! This function gives answer whether the given point is inside or outside the
 predefined polygon
 *  Unlike standard ray-casting algorithm, this one works on edges! (with
 performance benefit)

 * arguments:
 * Polygon - searched polygon
 * Point - an arbitrary point that can be inside or outside the polygon
 * length - the number of point in polygon (Attention! The list itself has an
 additional member - the last point coincides with the first)
 *
 *
 * return value:
 *  1 - the point is inside the polygon (including the case when point lies on
 some polygon's line)
 *  0 - the point is outside the polygon
 * -1 - the point is on edge
 */
int is_inside_sm(const Point *polygon, const Point point, int length) {
  int intersections = 0;
  double F, /*DY,*/ dy, dy2 = point.y - polygon[0].y;

  for (int ii = 0, jj = 1; ii < length; ++ii, ++jj) {
    dy = dy2;
    dy2 = point.y - polygon[jj].y;

    // consider only lines which are not completely above/bellow/right from the
    // point
    if (dy * dy2 <= 0. &&
        (point.x >= polygon[ii].x || point.x >= polygon[jj].x)) {
      // non-horizontal line
      if ((dy < 0.) || (dy2 < 0.)) // => dy*dy2 != 0
      {
        F = dy * (polygon[jj].x - polygon[ii].x) / (dy - dy2) + polygon[ii].x;

        if (point.x > F) // line is left from the point - the ray moving towards
                         // left, will intersect it
          ++intersections;
        else if (point.x == F) // point on line
          return -1;
      }

      // # point on upper peak (dy2=dx2=0) or horizontal line (dy=dy2=0 and
      // dx*dx2<=0)
      else if (dy2 == 0. &&
               (point.x == polygon[jj].x ||
                (dy == 0. &&
                 (point.x - polygon[ii].x) * (point.x - polygon[jj].x) <= 0.)))
        return -1;
    }
  }

  return intersections & 1;
}

//===================================================================
// This is how I made performance tests.
int main(int argc, char *argv[]) {
  printf("Point-in test\n\n");

  Point testpoint(7.f, 1.5f);

  int maxx = 8, minx = -6;
  int maxy = 4, miny = -4;
  int repeats = 7.5e7;
  clock_t t;

  srand(1389);
  t = clock();
  for (int jj = 0; jj < repeats; jj++) {
    testpoint.x = 0.1f * (rand() % (10 * maxx - 10 * minx + 1) + 10 * minx);
    testpoint.y = 0.1f * (rand() % (10 * maxy - 10 * miny + 1) + 10 * miny);
    cn_PnPoly(testpoint, test_polygon2, 38);
  }
  t = clock() - t;
  int iclock = (int)(t);
  double dsec = (double)(t) / CLOCKS_PER_SEC;
  printf("cn_CrossingNumber:  %d clocks,   %f seconds\n\n", iclock, dsec);

  srand(1389);
  t = clock();
  for (int jj = 0; jj < repeats; jj++) {
    testpoint.x = 0.1f * (rand() % (10 * maxx - 10 * minx + 1) + 10 * minx);
    testpoint.y = 0.1f * (rand() % (10 * maxy - 10 * miny + 1) + 10 * miny);
    wn_WindingNumber(testpoint, test_polygon2, 38);
  }
  t = clock() - t;
  iclock = (int)(t);
  dsec = (double)(t) / CLOCKS_PER_SEC;
  printf("wn_WindingNumber:  %d clocks,   %f seconds\n\n", iclock, dsec);

  srand(1389);
  t = clock();
  for (int jj = 0; jj < repeats; jj++) {
    testpoint.x = 0.1f * (rand() % (10 * maxx - 10 * minx + 1) + 10 * minx);
    testpoint.y = 0.1f * (rand() % (10 * maxy - 10 * miny + 1) + 10 * miny);
    is_inside_sm(test_polygon2, testpoint, 38);
  }
  t = clock() - t;
  iclock = (int)(t);
  dsec = (double)(t) / CLOCKS_PER_SEC;
  printf("is_inside_sm:  %d clocks,   %f seconds\n\n", iclock, dsec);

  return 0;
}
