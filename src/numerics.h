/*****************************************************************************/
/*                                                                           */
/* Written by:  Martin Held                                                  */
/*                                                                           */
/* E-Mail:      held@cs.sbg.ac.at                                            */
/* Fax Mail:    (+43 662) 8044-611                                           */
/* Voice Mail:  (+43 662) 8044-6304                                          */
/* Snail Mail:  Martin Held                                                  */
/*              FB Computerwissenschaften                                    */
/*              Universitaet Salzburg                                        */
/*              A-5020 Salzburg, Austria                                     */
/*                                                                           */
/*****************************************************************************/


#ifndef MIN_CONVEX_DECOMP_NUMERICS_H
#define MIN_CONVEX_DECOMP_NUMERICS_H

#define ZERO 10e-16

inline double Orientation(pnt* p1, pnt* p2, pnt* p3)
{
   return ( (p2->x - p1->x)*(p3->y - p1->y) -
            (p2->y - p1->y)*(p3->x - p1->x) );
}

inline int Sign(double x)
{
   if (x > 0.0)       return  1;
   else if (x < 0.0)  return -1;
   else               return  0;
}

inline boolean collinear(pnt* p1, pnt* p2, pnt* p3)
{
   double delta;

   delta = Orientation(p1,p2,p3);

   return ((delta >= -ZERO)  &&  (delta <= ZERO));
}


inline boolean CCW(pnt* p1, pnt* p2, pnt* p3)
{
   return Orientation(p1,p2,p3) >= 0.0;
}


inline boolean CW(pnt* p1, pnt* p2, pnt* p3)
{
   return Orientation(p1,p2,p3) <= 0.0;
}


inline double Area(pnt* p2, pnt* p3)
{
   return ( p2->x * p3->y -  p2->y * p3->x);
}



#endif
