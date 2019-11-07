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

#ifndef MIN_CONVEX_DECOMP_HEADERS_H
#define MIN_CONVEX_DECOMP_HEADERS_H

#include <list>
#include "data.h"

int p_comp(const void *, const void *);
int p_comp_y(const void *, const void *);
bool p_comp_y(pnt a, pnt b);

void ArgEval(int argc, char *argv[], rt_options *rt_opt);

FILE *OpenFile(const char *file_name, const char *access);

void initRand(rt_options *rt_opt);

void ReadInput(FILE *file, pnt *pnts, int *num_pnt);

void WriteOutputCoord(FILE *output, node *pnt, int num_pnt);

void WriteOutputIndex(FILE *output, node *pnt, int num_pnt);

void StartComputation(Data *data, rt_options &rt_opt, FILE *output);

void OnionLayers(pnt *pnts, int num_pnts, loop *layers, int *num_layers, 
                 node *nodes);

void ConvexHull(pnt *vtx, int num_vtx, int *ch_vtx, int *num_ch_vtx);

void WriteLayers(FILE *output, pnt *pnts, loop *layers, int num_layers, 
                 node *nodes, boolean obj);

int DetermineLowerBound(pnt *pnts, int num_pnts, loop *layers, 
                        int num_layers, node *nodes);

void GetCWmostVertexInHalfplane(int i1, int i2, int *i3, node *nodes, 
                                pnt *pnts);

void GetCCWmostVertexInHalfplane(int i1, int i2, int i_start, int *i3, 
                                 node *nodes, pnt *pnts);

void GetVerticesInHalfplane(int i1, int i2, int *o3, int *o4, node *nodes, 
                            pnt *pnts);

void ComputeApproxDecomp(Data *data, FILE *output, pnt *pnts, int num_pnts,
                         loop *layers, node *nodes,
                         boolean randomized, boolean obj,
						 rt_options &rt_opt);

void ComputeApproxDecompOnion(Data *data, FILE *output, pnt *pnts, int num_pnts,
                              loop *layers, int num_layers, node *nodes,
                              int lower_bound, boolean obj, rt_options &rt_opt);

void StoreAsOnionLayer(int *ch_vtx, int num_ch_vtx, 
                       loop *layers, int layer_id, 
                       pnt *vtx, node *nodes, boolean randomized);

boolean ExtractInnerPoints(pnt *vtx, int *num_vtx, 
                           int *ch_vtx, int num_ch_vtx);

boolean UpdateInnerPoints(pnt *vtx, int *num_vtx);

void FreeHulls(void);

void WriteOneLayer(FILE *output, pnt *pnts, loop *layers, int layer_id,
                   node *nodes, boolean obj);

void WriteOneTriLayer(FILE *output, pnt *pnts, int tri[], node *nodes,
                      boolean obj);

void HandleOnePoint(FILE *output, pnt *pnts, loop *layers, node *nodes,
                    int *num_cvx_areas, int L0, int *convex, boolean obj);

void HandleTwoPoints(FILE *output, pnt *pnts, loop *layers, node *nodes,
                     int *num_cvx_areas, int L0, int *convex, boolean obj);

void GeneralCase(FILE *output, pnt *pnts, loop *layers, node *nodes, 
                 pnt *vtx, int *num_cvx_areas, int *convex, boolean obj);

void SetConvexHullFlags(pnt *vtx, node *nodes, int start, int end);

inline void PrintNode(node* nodes, int nde, char const *s);

void HandleOnionAnnulus(FILE *output, pnt *pnts, loop *layers, node *nodes, 
                        int *num_cvx_areas, int L0, int *convex,
                        boolean obj);

void AddToConvexChain(int *convex, int *num_convex, node *nodes, int start,
                      int end, boolean ccw);

boolean WriteConvexChain(FILE *output, pnt *pnts, int *convex, 
                         int *num_convex, boolean obj);

void HandleDegenerateLoop(FILE *output, pnt *pnts, loop *layers, node *nodes, 
                          int *num_cvx_areas, int L0, int *convex,
                          boolean obj, int *k1, int *k2);

void WriteObjVertices(FILE *output, pnt *pnts, int num_pnts);

#endif
