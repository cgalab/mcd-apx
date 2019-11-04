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

void ArgEval(int argc, char *argv[], rt_options *rt_opt);

FILE *OpenFile(const char *file_name, const char *access);

void ReadInput(FILE *file, pnt *pnts, int *num_pnt);

void WriteOutputCoord(FILE *output, node *pnt, int num_pnt);

void WriteOutputIndex(FILE *output, node *pnt, int num_pnt);

void *ReallocateArray(void *old_ptr, int number, size_t size);

void OnionLayers(pnt *pnts, int num_pnts, loop *layers, int *num_layers, 
                 node *nodes, int *num_nodes);

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

void ComputeApproxDecomp(FILE *output, pnt *pnts, int num_pnts,
                         loop *layers, node *nodes,
                         boolean randomized, boolean obj);

void ComputeApproxDecompOnion(FILE *output, pnt *pnts, int num_pnts,
                              loop *layers, int num_layers, node *nodes,
                              int lower_bound, boolean obj);

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

boolean WriteConvexChain(FILE *output, pnt *pnts, node *nodes, int *convex, 
                         int *num_convex, boolean obj);

void HandleDegenerateLoop(FILE *output, pnt *pnts, loop *layers, node *nodes, 
                          int *num_cvx_areas, int L0, int *convex,
                          boolean obj);

void WriteObjVertices(FILE *output, pnt *pnts, int num_pnts);

#endif
