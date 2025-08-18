/* version 3.696.
   Written by Mary Kuhner, Jon Yamato, Joseph Felsenstein, Akiko Fuseki,
   Sean Lamont, and Andrew Keeffe.

   Copyright (c) 1993-2014, Joseph Felsenstein
   All rights reserved.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

#include <float.h>

#include "phylip.h"
#include "dist.h"

/* Static variable for thread count control */
static int neighbor_requested_threads = 0;

#ifdef OPENMP_ENABLED
#include <omp.h>
#endif

#ifndef OLDC
/* function prototypes */
void neighbor_getoptions(void);
void neighbor_allocrest(void);
void neighbor_doinit(void);
void neighbor_inputoptions(void);
void getinput(void);
void neighbor_describe(node *, double);
void neighbor_summarize(void);
void nodelabel(boolean);
void jointree(void);
void neighbor_maketree(void);
void freerest(void);
#ifdef OPENMP_ENABLED
void init_neighbor_parallel(int requested_threads);
#endif
/* function prototypes */
#endif


Char neighbor_infilename[FNMLNGTH], neighbor_outfilename[FNMLNGTH], neighbor_outtreename[FNMLNGTH];
long neighbor_nonodes2, neighbor_outgrno, neighbor_col, neighbor_datasets, neighbor_ith;
long neighbor_inseed;
vector *neighbor_x;
intvector *neighbor_reps;

boolean neighbor_jumble, neighbor_lower, neighbor_upper, neighbor_outgropt, neighbor_replicates, neighbor_trout,
               neighbor_printdata, neighbor_progress, neighbor_treeprint, neighbor_mulsets, neighbor_njoin;
tree neighbor_curtree;
longer neighbor_seed;
long *neighbor_enterorder;
Char neighbor_progname[20];

/* variables for neighbor_maketree, propagated globally for C version: */
node **neighbor_cluster;

#ifdef OPENMP_ENABLED
/* Parallelization variables */
static int neighbor_num_threads = 1;

void init_neighbor_parallel(int requested_threads)
{
  /* Initialize parallel processing with user-specified thread count */
#ifdef AUTO_THREAD_DETECTION
  int max_threads = omp_get_max_threads();
  /* Use requested threads or auto-detect if 0 */
  if (requested_threads > 0) {
    neighbor_num_threads = (requested_threads <= max_threads) ? requested_threads : max_threads;
  } else {
    neighbor_num_threads = max_threads; /* Auto-detect */
  }
  omp_set_num_threads(neighbor_num_threads);
#else
  neighbor_num_threads = omp_get_max_threads();
#endif
  
  if (neighbor_progress) {
    // printf("Neighbor-joining method: OpenMP parallel processing with %d thread(s)", neighbor_num_threads);
    if (requested_threads > 1) printf(" (user-specified)");
  }
}
#endif


void neighbor_getoptions()
{
  /* interactively set options */
  long inseed0 = 0, loopcount;

  fprintf(outfile, "\nNeighbor-Joining/UPGMA method version %s\n\n",VERSION);
  neighbor_jumble = false;
  neighbor_lower = false;
  neighbor_outgrno = 1;
  neighbor_outgropt = false;
  neighbor_replicates = false;
  neighbor_trout = true;
  neighbor_upper = false;
  neighbor_printdata = false;
  neighbor_progress = true;
  neighbor_treeprint = true;
  neighbor_njoin = true;
  loopcount = 0;
}  /* neighbor_getoptions */


void neighbor_allocrest()
{
  long i;

  neighbor_x = (vector *)Malloc(spp*sizeof(vector));
  for (i = 0; i < spp; i++)
    neighbor_x[i] = (vector)Malloc(spp*sizeof(double));
  neighbor_reps = (intvector *)Malloc(spp*sizeof(intvector));
  for (i = 0; i < spp; i++)
    neighbor_reps[i] = (intvector)Malloc(spp*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  neighbor_enterorder = (long *)Malloc(spp*sizeof(long));
  neighbor_cluster = (node **)Malloc(spp*sizeof(node *));
}  /* neighbor_allocrest */


void freerest()
{
  long i;

  for (i = 0; i < spp; i++)
    free(neighbor_x[i]);
  free(neighbor_x);
  for (i = 0; i < spp; i++)
    free(neighbor_reps[i]);
  free(neighbor_reps);
  free(nayme);
  free(neighbor_enterorder);
  free(neighbor_cluster);
}  /* freerest */


void neighbor_doinit()
{
  /* initializes variables */
  node *p;

  inputnumbers2(&spp, &neighbor_nonodes2, 2);
  neighbor_nonodes2 += (neighbor_njoin ? 0 : 1);
  neighbor_getoptions();
  alloctree(&neighbor_curtree.nodep, neighbor_nonodes2+1);
  p = neighbor_curtree.nodep[neighbor_nonodes2]->next;
  neighbor_curtree.nodep[neighbor_nonodes2]->next = neighbor_curtree.nodep[neighbor_nonodes2];
  free(p->next);
  free(p);
  neighbor_allocrest();

}  /* neighbor_doinit */


void neighbor_inputoptions()
{
  /* read options information */

  if (neighbor_ith != 1)
    samenumsp2(neighbor_ith);
  putc('\n', outfile);
  if (neighbor_njoin)
    fprintf(outfile, " Neighbor-joining method\n");
  else
    fprintf(outfile, " UPGMA method\n");
  fprintf(outfile, "\n Negative branch lengths allowed\n\n");
}  /* neighbor_inputoptions */


void neighbor_describe(node *p, double height)
{
  /* print out information for one branch */
  long i;
  node *q;

  q = p->back;
  if (neighbor_njoin)
    fprintf(outfile, "%4ld          ", q->index - spp);
  else
    fprintf(outfile, "%4ld     ", q->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
    putc(' ', outfile);
  } else {
    if (neighbor_njoin)
      fprintf(outfile, "%4ld       ", p->index - spp);
    else {
      fprintf(outfile, "%4ld       ", p->index - spp);
    }
  }
  if (neighbor_njoin)
    fprintf(outfile, "%12.5f\n", q->v);
  else
    fprintf(outfile, "%10.5f      %10.5f\n", q->v, q->v+height);
  if (!p->tip) {
    neighbor_describe(p->next->back, height+q->v);
    neighbor_describe(p->next->next->back, height+q->v);
  }
}  /* neighbor_describe */


void neighbor_summarize()
{
  /* print out branch lengths etc. */
  putc('\n', outfile);
  if (neighbor_njoin) {
    fprintf(outfile, "remember:");
    if (neighbor_outgropt)
      fprintf(outfile, " (although rooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  if (neighbor_njoin) {
    fprintf(outfile, "\nBetween        And            Length\n");
    fprintf(outfile, "-------        ---            ------\n");
  } else {
    fprintf(outfile, "From     To            Length          Height\n");
    fprintf(outfile, "----     --            ------          ------\n");
  }
  neighbor_describe(neighbor_curtree.start->next->back, 0.0);
  neighbor_describe(neighbor_curtree.start->next->next->back, 0.0);
  if (neighbor_njoin)
    neighbor_describe(neighbor_curtree.start->back, 0.0);
  fprintf(outfile, "\n\n");
}  /* neighbor_summarize */


void nodelabel(boolean isnode)
{
 /*
  if (isnode)
    printf("node");
  else
    printf("species");
    */
}  /* nodelabel */


void jointree()
{
  /* calculate the tree */
  long nc, nextnode, mini=0, minj=0, i, j, ia, ja, ii, jj, nude, iter;
  double fotu2, total, tmin, dio, djo, bi, bj, bk, dmin=0, da;
  long el[3];
  vector av;
  intvector oc;

  double *R;   /* added in revisions by Y. Ina */
  R = (double *)Malloc(spp * sizeof(double));

  for (i = 0; i <= spp - 2; i++) {
#ifdef OPENMP_ENABLED
    #pragma omp parallel for private(j, da) schedule(dynamic)
#endif
    for (j = i + 1; j < spp; j++) {
      da = (neighbor_x[i][j] + neighbor_x[j][i]) / 2.0;
      neighbor_x[i][j] = da;
      neighbor_x[j][i] = da;
    }
  }
  /* First initialization */
  fotu2 = spp - 2.0;
  nextnode = spp + 1;
  av = (vector)Malloc(spp*sizeof(double));
  oc = (intvector)Malloc(spp*sizeof(long));
  for (i = 0; i < spp; i++) {
    av[i] = 0.0;
    oc[i] = 1;
  }
  /* Enter the main cycle */
  if (neighbor_njoin)
    iter = spp - 3;
  else
    iter = spp - 1;
  for (nc = 1; nc <= iter; nc++) {
    for (j = 2; j <= spp; j++) {
      for (i = 0; i <= j - 2; i++)
        neighbor_x[j - 1][i] = neighbor_x[i][j - 1];
    }
    tmin = DBL_MAX;
    /* Compute sij and minimize */
    if (neighbor_njoin) {     /* many revisions by Y. Ina from here ... */
#ifdef OPENMP_ENABLED
      #pragma omp parallel for private(i) schedule(static)
#endif
      for (i = 0; i < spp; i++)
        R[i] = 0.0;
      for (ja = 2; ja <= spp; ja++) {
        jj = neighbor_enterorder[ja - 1];
        if (neighbor_cluster[jj - 1] != NULL) {
#ifdef OPENMP_ENABLED
          #pragma omp parallel for private(ia, ii) reduction(+:R[:spp]) schedule(dynamic)
#endif
          for (ia = 0; ia <= ja - 2; ia++) {
            ii = neighbor_enterorder[ia];
            if (neighbor_cluster[ii - 1] != NULL) {
              R[ii - 1] += neighbor_x[ii - 1][jj - 1];
              R[jj - 1] += neighbor_x[ii - 1][jj - 1];
            }
          }
        }
      }
    } /* ... to here */
    for (ja = 2; ja <= spp; ja++) {
      jj = neighbor_enterorder[ja - 1];
      if (neighbor_cluster[jj - 1] != NULL) {
#ifdef OPENMP_ENABLED
        #pragma omp parallel for private(ia, ii, total) schedule(dynamic)
#endif
        for (ia = 0; ia <= ja - 2; ia++) {
          ii = neighbor_enterorder[ia];
          if (neighbor_cluster[ii - 1] != NULL) {
            if (neighbor_njoin) {
              total = fotu2 * neighbor_x[ii - 1][jj - 1] - R[ii - 1] - R[jj - 1];
               /* this statement part of revisions by Y. Ina */
            } else
              total = neighbor_x[ii - 1][jj - 1];
#ifdef OPENMP_ENABLED
            #pragma omp critical
#endif
            {
              if (total < tmin) {
                tmin = total;
                mini = ii;
                minj = jj;
              }
            }
          }
        }
      }
    }
    /* compute lengths and print */
    if (neighbor_njoin) {
      dio = 0.0;
      djo = 0.0;
#ifdef OPENMP_ENABLED
      #pragma omp parallel for private(i) reduction(+:dio,djo) schedule(static)
#endif
      for (i = 0; i < spp; i++) {
        dio += neighbor_x[i][mini - 1];
        djo += neighbor_x[i][minj - 1];
      }
      dmin = neighbor_x[mini - 1][minj - 1];
      dio = (dio - dmin) / fotu2;
      djo = (djo - dmin) / fotu2;
      bi = (dmin + dio - djo) * 0.5;
      bj = dmin - bi;
      bi -= av[mini - 1];
      bj -= av[minj - 1];
    } else {
      bi = neighbor_x[mini - 1][minj - 1] / 2.0 - av[mini - 1];
      bj = neighbor_x[mini - 1][minj - 1] / 2.0 - av[minj - 1];
      av[mini - 1] += bi;
    }

    if (neighbor_progress) {
    /*  printf("Cycle %3ld: ", iter - nc + 1);*/
      if (neighbor_njoin)
        nodelabel((boolean)(av[mini - 1] > 0.0));
      else
        nodelabel((boolean)(oc[mini - 1] > 1.0));
     /* printf(" %ld (%10.5f) joins ", mini, bi);*/
      if (neighbor_njoin)
        nodelabel((boolean)(av[minj - 1] > 0.0));
      else
        nodelabel((boolean)(oc[minj - 1] > 1.0));
      /*printf(" %ld (%10.5f)\n", minj, bj);*/
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    hookup(neighbor_curtree.nodep[nextnode - 1]->next, neighbor_cluster[mini - 1]);
    hookup(neighbor_curtree.nodep[nextnode - 1]->next->next, neighbor_cluster[minj - 1]);
    neighbor_cluster[mini - 1]->v = bi;
    neighbor_cluster[minj - 1]->v = bj;
    neighbor_cluster[mini - 1]->back->v = bi;
    neighbor_cluster[minj - 1]->back->v = bj;
    neighbor_cluster[mini - 1] = neighbor_curtree.nodep[nextnode - 1];
    neighbor_cluster[minj - 1] = NULL;
    nextnode++;
    if (neighbor_njoin)
      av[mini - 1] = dmin * 0.5;
    /* re-initialization */
    fotu2 -= 1.0;
#ifdef OPENMP_ENABLED
    #pragma omp parallel for private(j, da) schedule(dynamic)
#endif
    for (j = 0; j < spp; j++) {
      if (neighbor_cluster[j] != NULL) {
        if (neighbor_njoin) {
          da = (neighbor_x[mini - 1][j] + neighbor_x[minj - 1][j]) * 0.5;
          if (mini - j - 1 < 0)
            neighbor_x[mini - 1][j] = da;
          if (mini - j - 1 > 0)
            neighbor_x[j][mini - 1] = da;
        } else {
          da = neighbor_x[mini - 1][j] * oc[mini - 1] + neighbor_x[minj - 1][j] * oc[minj - 1];
          da /= oc[mini - 1] + oc[minj - 1];
          neighbor_x[mini - 1][j] = da;
          neighbor_x[j][mini - 1] = da;
        }
      }
    }
#ifdef OPENMP_ENABLED
    #pragma omp parallel for private(j) schedule(static)
#endif
    for (j = 0; j < spp; j++) {
      neighbor_x[minj - 1][j] = 0.0;
      neighbor_x[j][minj - 1] = 0.0;
    }
    oc[mini - 1] += oc[minj - 1];
  }
  /* the last cycle */
  nude = 1;
  for (i = 1; i <= spp; i++) {
    if (neighbor_cluster[i - 1] != NULL) {
      el[nude - 1] = i;
      nude++;
    }
  }
  if (!neighbor_njoin) {
    neighbor_curtree.start = neighbor_cluster[el[0] - 1];
    neighbor_curtree.start->back = NULL;
    free(av);
    free(oc);
    return;
  }
  bi = (neighbor_x[el[0] - 1][el[1] - 1] + neighbor_x[el[0] - 1][el[2] - 1] - neighbor_x[el[1] - 1]
        [el[2] - 1]) * 0.5;
  bj = neighbor_x[el[0] - 1][el[1] - 1] - bi;
  bk = neighbor_x[el[0] - 1][el[2] - 1] - bi;
  bi -= av[el[0] - 1];
  bj -= av[el[1] - 1];
  bk -= av[el[2] - 1];
  if (neighbor_progress) {
    /*printf("last cycle:\n");*/
    putchar(' ');
    nodelabel((boolean)(av[el[0] - 1] > 0.0));
    /*printf(" %ld  (%10.5f) joins ", el[0], bi);*/
    nodelabel((boolean)(av[el[1] - 1] > 0.0));
    /*printf(" %ld  (%10.5f) joins ", el[1], bj);*/
    nodelabel((boolean)(av[el[2] - 1] > 0.0));
   /* printf(" %ld  (%10.5f)\n", el[2], bk);*/
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
  hookup(neighbor_curtree.nodep[nextnode - 1], neighbor_cluster[el[0] - 1]);
  hookup(neighbor_curtree.nodep[nextnode - 1]->next, neighbor_cluster[el[1] - 1]);
  hookup(neighbor_curtree.nodep[nextnode - 1]->next->next, neighbor_cluster[el[2] - 1]);
  neighbor_cluster[el[0] - 1]->v = bi;
  neighbor_cluster[el[1] - 1]->v = bj;
  neighbor_cluster[el[2] - 1]->v = bk;
  neighbor_cluster[el[0] - 1]->back->v = bi;
  neighbor_cluster[el[1] - 1]->back->v = bj;
  neighbor_cluster[el[2] - 1]->back->v = bk;
  neighbor_curtree.start = neighbor_cluster[el[0] - 1]->back;
  free(av);
  free(oc);
  free(R);
}  /* jointree */


void neighbor_maketree()
{
  /* construct the tree */
  long i ;

#ifdef OPENMP_ENABLED
  init_neighbor_parallel(neighbor_requested_threads);
#endif

  inputdata(neighbor_replicates, neighbor_printdata, neighbor_lower, neighbor_upper, neighbor_x, neighbor_reps);
  if (neighbor_njoin && (spp < 3)) {
    printf("\nERROR: Neighbor-Joining runs must have at least 3 species\n\n");
    exxit(-1);
  }
 /* if (neighbor_progress)
    putchar('\n');*/
  if (neighbor_ith == 1)
    setuptree(&neighbor_curtree, neighbor_nonodes2 + 1);
  for (i = 1; i <= spp; i++)
    neighbor_enterorder[i - 1] = i;
  if (neighbor_jumble)
    randumize(neighbor_seed, neighbor_enterorder);
  for (i = 0; i < spp; i++)
    neighbor_cluster[i] = neighbor_curtree.nodep[i];
  jointree();
  if (neighbor_njoin)
    neighbor_curtree.start = neighbor_curtree.nodep[neighbor_outgrno - 1]->back;
  printree(neighbor_curtree.start, neighbor_treeprint, neighbor_njoin, (boolean)(!neighbor_njoin));
  if (neighbor_treeprint)
    neighbor_summarize();
  if (neighbor_trout) {
    neighbor_col = 0;
    if (neighbor_njoin)
      treeout(neighbor_curtree.start, &neighbor_col, 0.43429448222, neighbor_njoin, neighbor_curtree.start);
    else
      neighbor_curtree.root = neighbor_curtree.start,
      treeoutr(neighbor_curtree.start,&neighbor_col,&neighbor_curtree);
  }
 /* if (neighbor_progress) {
    printf("\nOutput written on file \"%s\"\n\n", neighbor_outfilename);
    if (neighbor_trout)
      printf("Tree written on file \"%s\"\n\n", neighbor_outtreename);
  }*/
}  /* neighbor_maketree */


int neighbor_build_tree(const char *path_name_infile, const char *path_name_outfile, int num_threads)
{  /* main program */
  /* Store requested thread count in static variable for use by init_neighbor_parallel */
  neighbor_requested_threads = num_threads;
  
  init();
  char *outfile_path_name=malloc(strlen(path_name_outfile)+strlen("_report")+1);
  char *outtree_path_name=malloc(strlen(path_name_outfile)+strlen("_tree.nwk")+1);
  strcpy(outfile_path_name,path_name_outfile);
  strcat(outfile_path_name,"_report");
  strcpy(outtree_path_name,path_name_outfile);
  strcat(outtree_path_name,"_tree.nwk");
  openfile(&infile,path_name_infile,"input file", "r",neighbor_infilename);
  openfile(&outfile,outfile_path_name,"output file", "w",neighbor_outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  neighbor_mulsets = false;
  neighbor_datasets = 1;
  neighbor_doinit();
  if (neighbor_trout)
    openfile(&outtree,outtree_path_name,"output tree file", "w",neighbor_outtreename);
  neighbor_ith = 1;
  while (neighbor_ith <= neighbor_datasets) {
    if (neighbor_datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n",neighbor_ith);
     /* if (neighbor_progress)
        printf("Data set # %ld:\n",neighbor_ith);*/
    }
    neighbor_inputoptions();
    neighbor_maketree();
    if (eoln(infile) && (neighbor_ith < neighbor_datasets))
      scan_eoln(infile);
    neighbor_ith++;
  }
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
  freerest();
  freetree(&neighbor_curtree.nodep, neighbor_nonodes2+1);
  free(outfile_path_name);
  free(outtree_path_name);
#ifdef MAC
  fixmacfile(neighbor_outfilename);
  fixmacfile(neighbor_outtreename);
#endif
 /* printf("Done.\n\n");*/
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}





