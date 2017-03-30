/*  geom.cpp
*
*  gift wrapping algorithm for finding the convex hull on a set of points in 3D
*
*  created by bridget went and claire mccarthy 2/26/17
*  last modified 3/1/17
*
*/


#include "geom.h"

using namespace std;

//dequeue<edge3d> frontier;
unordered_set<edge3d, edge_hash, edge_equals> Q;

/* returns 6 times the signed volume of the polyhedron formed by triangle abc
and point d. formula given in Computational Geometry textbook.
*/
int signed_volume(point3d a, point3d b, point3d c, point3d d) {

  return (((-(a.z - d.z)*(b.y - d.y)*(c.x - d.x) + (a.y - d.y)*(b.z - d.z)*(c.x - d.x) +
  (a.z - d.z)*(b.x - d.x)*(c.y - d.y) - (a.x - d.x)*(b.z - d.z)*(c.y - d.y) -
  (a.y - d.y)*(b.x - d.x)*(c.z - d.z) + (a.x - d.x)*(b.y - d.y)*(c.z - d.z))));

}

/* return 1 if p,q,r, t on same plane, and 0 otherwise */
int coplanar(point3d p, point3d q, point3d r, point3d t) {
  return (int)(signed_volume(p,q,r,t) == 0);
}


/* return 1 if d is  strictly left of abc; 0 otherwise */
int left(point3d a, point3d b, point3d c, point3d d) {
  return (int)(signed_volume(a,b,c,d) < 0);
}

int isEqual(point3d a, point3d b) {
  return (a.x == b.x && a.y == b.y && a.z == b.z);
}



/* returns 1 if exactly one face adjacent to edge e has been discovered so far */
int edgeProcessed(edge3d &e) {
  // find the edge in the hash map
  unordered_set<edge3d, edge_hash, edge_equals>::const_iterator got;
  got = Q.find(e);
  if (got != Q.end()) {
    return 1;
  }
  return 0;
}

edge3d makeEdge(point3d *p, point3d *q) {
  edge3d e = {.p = p, .q = q};
  return e;
}

/* returns the point p that makes the minimal angle with edge e */
point3d* findPointWithMinAngle(edge3d &e, vector<point3d> &points) {

  bool extreme = true;
  bool isLeft = false;

  assert(!isEqual(*e.p, *e.q));

  for (int i = 0; i < points.size(); i++) {

    if (e.left != NULL) {

      if (isEqual(points[i], *e.p) || isEqual(points[i], *e.q) ||
      isEqual(points[i], *e.left->p1) || isEqual(points[i], *e.left->p2) ||
      isEqual(points[i], *e.left->p3)) {
        continue;
      } else {
        if (isEqual(points[i], *e.p) || isEqual(points[i], *e.q)) {
          continue;
        }
      }

      extreme = true;

      // see if first point is left or right
      if (left(*e.p, *e.q, points[i], points[0])) {
        isLeft = true;
      } else {
        isLeft = false;
      }


      // now check if plane defined by this edge with the current point i is extreme
      for (int p = 0; p < points.size(); ++p) {

        if (isEqual(points[p], points[i]) || isEqual(points[p], *e.p) || isEqual(points[p], *e.q)) {
          continue;
        }


        if (coplanar(*e.p, *e.q, points[i], points[p])) {
          printf("coplanar case\n");
          break;
        }
        if (isLeft) {
          if (!left(*e.p, *e.q, points[i], points[p])) {
            extreme = false;
            break;
          }
        } else {
          if (left(*e.p, *e.q, points[i], points[p])) {
            extreme = false;
            break;
          }
        }
      }
      if (extreme) {
        return &points[i];
      }
    }

  }

  printf("could not find any point with minimum edge\n");
  return &points[0];
}

//  float minAngle;
//  point3d* minPoint = &points[0];
//
//
//  minPoint->x = points[0].x;
//  minPoint->y = points[0].y;
//
//  float currAngle;
//  float prod;
//  float mod1;
//  float mod2;
//  float cos0;
//  prod = (e.p->x - e.q->x)*(e.q->x - points[0].x) + (e.p->y - e.q->y)*(e.q->y - points[0].y)+ (e.p->z - e.q->z)*(e.q->z - points[0].z);
//  mod1 = sqrt((e.p->x - e.q->x)*(e.p->x - e.q->x) + (e.p->y - e.q->y)*(e.p->y - e.q->y) + (e.p->z - e.q->z)*(e.p->z - e.q->z));
//  mod2 = sqrt((e.q->x - points[0].x)*(e.q->x - points[0].x) + (e.q->y - points[0].y)*(e.q->y - points[0].y) + (e.q->z - points[0].z)*(e.q->z - points[0].z));
//  cos0 = (prod/((mod1)*(mod2)));
//  currAngle = acos(cos0)*180/M_PI;
//  minAngle = currAngle;
//  for (int i = 0; i < points.size(); i++) {
//    //calculate angle
//    if (isEqual(points[i], *e.p) || isEqual(points[i], *e.q)) {
//      continue;
//    }
//    prod = (e.p->x - e.q->x)*(e.q->x - points[i].x) + (e.p->y - e.q->y)*(e.q->y - points[i].y)+ (e.p->z - e.q->z)*(e.q->z - points[i].z);
//    mod1 = sqrt((e.p->x - e.q->x)*(e.p->x - e.q->x) + (e.p->y - e.q->y)*(e.p->y - e.q->y) + (e.p->z - e.q->z)*(e.p->z - e.q->z));
//    mod2 = sqrt((e.q->x - points[i].x)*(e.q->x - points[i].x) + (e.q->y - points[i].y)*(e.q->y - points[i].y) + (e.q->z - points[i].z)*(e.q->z - points[i].z));
//    cos0 = (prod/((mod1)*(mod2)));
//    currAngle = acos(cos0)*180/M_PI;
//
//    if (currAngle < minAngle) {
//      minAngle = currAngle;
//      minPoint = &points[i];
//      minPoint->x = points[i].x;
//      minPoint->y = points[i].y;
//
//    }
//  }
// return minPoint;

//}


/* computes the first face on the hull */
face3d findInitialFace(vector<point3d> &points) {
  /* find point with min-y coord */
  //float minY_val = FLOAT_MAX;
  int minY_index = -1;
  //  int minX_index = -1;
  int minY = BIG_INT;

  for (int i = 0; i < points.size(); i++) {
    if (points[i].y < minY) {
      minY_index = i;
      minY = points[i].y;
    }
  }
  point3d *start = &points[minY_index];

  bool extreme = true;

  for (int j = 0; j < points.size(); ++j) {
    for (int k = j + 1; k < points.size(); ++k) {
      if (isEqual(*start, points[j]) || isEqual(*start, points[k])) {
        continue;
      }
      extreme = true;



      for (int p = 0; p < points.size(); ++p) {

        if(isEqual(points[p], *start) || p==j || p==k) {
          continue;
        }
        printf("j: %d k: %d p: %d\n", j, k, p);
        if (coplanar(*start, points[j], points[k], points[p])) {
          extreme = true;
          break;
        }
        if (!left(*start, points[j], points[k], points[p])) {
          extreme = false;
          break;
        }

      }
      if (extreme) {
        // create a face and add this face to the hull
        //  printf("extreme\n");
        face3d face = {.p1 = start, .p2 = &points[j], .p3 = &points[k]};
        printf("FIRST FACE:\n \t %d %d %d\n \t %d %d %d\n \t %d %d %d\n",
        face.p1->x, face.p1->y, face.p1->z,
        face.p2->x, face.p2->y, face.p2->z,
        face.p3->x, face.p3->y, face.p3->z);
        printf("j: %d k: %d\n", j, k);
        return face;
      }
    }
  }
  printf("can't find extreme face\n");
  exit(1);
  // face3d face = {.p1 = start, .p2 = start, .p3 = start};
  // return face;

}



/* compute and return the convex hull of the points */
vector<face3d> gift_wrapping_hull(vector<point3d> &points) {

  vector<face3d> result;

  assert(result.empty());
  assert(Q.empty());

  // find a face guaranteed to be on the hull
  face3d face = findInitialFace(points);
  result.push_back(face);

  edge3d e12 = {.p = face.p1, .q = face.p2, .left = &face};
  edge3d e23 = {.p = face.p2, .q = face.p3, .left = &face};
  edge3d e31 = {.p = face.p3, .q = face.p1, .left = &face};

  Q.insert(e12);
  Q.insert(e23);
  Q.insert(e31);

  unordered_set<edge3d, edge_hash, edge_equals>::iterator it;


  while (!Q.empty()) {
    //printf("frontier size: %d\n", frontier.size());
    it = Q.begin();
    edge3d e = *it;
    Q.erase(it);

    //edges.erase(e);


    //printf("\n edge = (%d %d %d) \t (%d %d %d)\n", e.p->x, e.p->y, e.p->z, e.q->x, e.q->y, e.q->z);

    assert(!isEqual(*e.p, *e.q));

    point3d* next = findPointWithMinAngle(e, points);
    //printf("next = %d %d %d\n", next->x, next->y, next->z);

    assert(!isEqual(*e.p, *next) && !isEqual(*e.q, *next));

    face3d f = {.p1 = e.p, .p2 = e.q, .p3 = next};
    //printf("adding face: (%d %d %d)\t (%d %d %d)\t (%d %d %d)\n", e.p->x, e.p->y, e.p->z, e.q->x, e.q->y, e.q->z, next->x, next->y, next->z);
    result.push_back(f); // push onto the hull


    // add next 2 edges
    edge3d e31 = {.p = next, .q = f.p1, .left = &f};
    edge3d e23 = {.p = f.p2, .q = next, .left = &f};

    // if both are in the hash, last face
    if (edgeProcessed(e31) && edgeProcessed(e23)) {
      it = Q.find(e31);
      Q.erase(it);
      it = Q.find(e23);
      Q.erase(it);
      //break;

    } else {

      //printf("edge: (%d %d %d), (%d %d %d) ", e31.p->x, e31.p->y, e31.p->z, e31.q->x, e31.q->y, e31.q->z);
      if (!edgeProcessed(e31)) {
        printf("push\n");
        Q.insert(e31);
        // frontier.push(e31);
        // edges.insert(e31);
      } else {
        printf("already seen this edge; removing.\n");
        //edges.erase(e31);
        it = Q.find(e31);
        Q.erase(it);
        //frontier.pop();
      }
      //printf("edge: (%d %d %d), (%d %d %d) ", e23.p->x, e23.p->y, e23.p->z, e23.q->x, e23.q->y, e23.q->z);
      if (!edgeProcessed(e23)) {
        printf("push\n");
        //frontier.push(e23);
        //edges.insert(e23);
        Q.insert(e23);
      } else {
        printf("removing edge: (%d %d %d) (%d %d %d)\n", e23.p->x,e23.p->y,e23.p->z, e23.q->x, e23.q->y, e23.q->z);
        unordered_set<edge3d, edge_hash, edge_equals>::iterator f = Q.find(e23);
        Q.erase(f);
      }
    }


    printf("edges size : %lu\n", Q.size());

  }

  for (int p = 0; p < result.size(); ++p) {
    /*printf("hull face %d:\n \t %d %d %d\n \t %d %d %d\n \t %d %d %d\n",
    p, result[p].p1->x, result[p].p1->y, result[p].p1->z,
    result[p].p2->x, result[p].p2->y, result[p].p2->z,
    result[p].p3->x, result[p].p3->y, result[p].p3->z);
*/
  }

  assert(Q.empty());

  return result;
}
