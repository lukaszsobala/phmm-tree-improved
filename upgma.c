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

#ifdef OPENMP_ENABLED
#include <omp.h>
#endif

#ifndef OLDC
/* function prototypes */
void upgma_getoptions(void);
void upgma_allocrest(void);
void upgma_doinit(void);
void upgma_inputoptions(void);
void upgma_getinput(void);
void upgma_describe(node *, double);
void upgma_summarize(void);
void upgma_nodelabel(boolean);
void upgma_jointree(void);
void upgma_maketree(void);
void upgma_freerest(void);
#ifdef OPENMP_ENABLED
void init_upgma_parallel(void);
#endif
/* function prototypes */
#endif


static Char infilename[FNMLNGTH], outfilename[FNMLNGTH], outtreename[FNMLNGTH];
static long nonodes2, outgrno, col, datasets, ith;
static long inseed;
static vector *x;
static intvector *reps;
static boolean jumble, lower, upper, outgropt, replicates, trout,
               printdata, progress, treeprint, mulsets, njoin;
static tree curtree;
static longer seed;
static long *enterorder;
static Char progname[20];

/* OpenMP parallel processing variables */
#ifdef OPENMP_ENABLED
static int upgma_num_threads = 1;
#endif

/* variables for upgma_maketree, propagated globally for C version: */
node **cluster;


void upgma_getoptions()
{
  /* interactively set options */
  long inseed0 = 0, loopcount;

  fprintf(outfile, "\nNeighbor-Joining/UPGMA method version %s\n\n",VERSION);
  jumble = false;
  lower = false;
  outgrno = 1;
  outgropt = false;
  replicates = false;
  trout = true;
  upper = false;
  printdata = false;
  progress = true;
  treeprint = true;
  njoin = false;
  loopcount = 0;
}  /* upgma_getoptions */


#ifdef OPENMP_ENABLED
void init_upgma_parallel(void)
{
  /* Initialize parallel processing */
#ifdef AUTO_THREAD_DETECTION
  upgma_num_threads = omp_get_max_threads();
  omp_set_num_threads(upgma_num_threads);
#else
  upgma_num_threads = omp_get_max_threads();
#endif
  
  if (progress) {
    printf("UPGMA method: OpenMP parallel processing with %d thread(s)\n", upgma_num_threads);
  }
}
#endif


void upgma_allocrest()
{
  long i;

  x = (vector *)Malloc(spp*sizeof(vector));
  for (i = 0; i < spp; i++)
    x[i] = (vector)Malloc(spp*sizeof(double));
  reps = (intvector *)Malloc(spp*sizeof(intvector));
  for (i = 0; i < spp; i++)
    reps[i] = (intvector)Malloc(spp*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
  cluster = (node **)Malloc(spp*sizeof(node *));
}  /* upgma_allocrest */


void upgma_freerest()
{
  long i;

  for (i = 0; i < spp; i++)
    free(x[i]);
  free(x);
  for (i = 0; i < spp; i++)
    free(reps[i]);
  free(reps);
  free(nayme);
  free(enterorder);
  free(cluster);
}  /* upgma_freerest */


void upgma_doinit()
{
  /* initializes variables */
  node *p;

  inputnumbers2(&spp, &nonodes2, 2);
  nonodes2 += (njoin ? 0 : 1);
  upgma_getoptions();
  alloctree(&curtree.nodep, nonodes2+1);
  p = curtree.nodep[nonodes2]->next;
  curtree.nodep[nonodes2]->next = curtree.nodep[nonodes2];
  free(p->next);
  free(p);
  upgma_allocrest();

}  /* upgma_doinit */


void upgma_inputoptions()
{
  /* read options information */

  if (ith != 1)
    samenumsp2(ith);
  putc('\n', outfile);
  if (njoin)
    fprintf(outfile, " Neighbor-joining method\n");
  else
    fprintf(outfile, " UPGMA method\n");
  fprintf(outfile, "\n Negative branch lengths allowed\n\n");
}  /* upgma_inputoptions */


void upgma_describe(node *p, double height)
{
  /* print out information for one branch */
  long i;
  node *q;

  q = p->back;
  if (njoin)
    fprintf(outfile, "%4ld          ", q->index - spp);
  else
    fprintf(outfile, "%4ld     ", q->index - spp);
  if (p->tip) {
    for (i = 0; i < nmlngth; i++)
      putc(nayme[p->index - 1][i], outfile);
    putc(' ', outfile);
  } else {
    if (njoin)
      fprintf(outfile, "%4ld       ", p->index - spp);
    else {
      fprintf(outfile, "%4ld       ", p->index - spp);
    }
  }
  if (njoin)
    fprintf(outfile, "%12.5f\n", q->v);
  else
    fprintf(outfile, "%10.5f      %10.5f\n", q->v, q->v+height);
  if (!p->tip) {
    upgma_describe(p->next->back, height+q->v);
    upgma_describe(p->next->next->back, height+q->v);
  }
}  /* upgma_describe */


void upgma_summarize()
{
  /* print out branch lengths etc. */
  putc('\n', outfile);
  if (njoin) {
    fprintf(outfile, "remember:");
    if (outgropt)
      fprintf(outfile, " (although rooted by outgroup)");
    fprintf(outfile, " this is an unrooted tree!\n");
  }
  if (njoin) {
    fprintf(outfile, "\nBetween        And            Length\n");
    fprintf(outfile, "-------        ---            ------\n");
  } else {
    fprintf(outfile, "From     To            Length          Height\n");
    fprintf(outfile, "----     --            ------          ------\n");
  }
  upgma_describe(curtree.start->next->back, 0.0);
  upgma_describe(curtree.start->next->next->back, 0.0);
  if (njoin)
    upgma_describe(curtree.start->back, 0.0);
  fprintf(outfile, "\n\n");
}  /* upgma_summarize */


void upgma_nodelabel(boolean isnode)
{
 /*
  if (isnode)
    printf("node");
  else
    printf("species");
    */
}  /* upgma_nodelabel */


void upgma_jointree()
{
  /* calculate the tree */
  long nc, nextnode, mini=0, minj=0, i, j, ia, ja, ii, jj, nude, iter;
  double fotu2, total, tmin, dio, djo, bi, bj, bk, dmin=0, da;
  long el[3];
  vector av;
  intvector oc;

  double *R;   /* added in revisions by Y. Ina */
  R = (double *)Malloc(spp * sizeof(double));

#ifdef OPENMP_ENABLED
  #pragma omp parallel for private(j, da) schedule(dynamic)
#endif
  for (i = 0; i <= spp - 2; i++) {
    for (j = i + 1; j < spp; j++) {
      da = (x[i][j] + x[j][i]) / 2.0;
      x[i][j] = da;
      x[j][i] = da;
    }
  }
  /* First initialization */
  fotu2 = spp - 2.0;
  nextnode = spp + 1;
  av = (vector)Malloc(spp*sizeof(double));
  oc = (intvector)Malloc(spp*sizeof(long));
#ifdef OPENMP_ENABLED
  #pragma omp parallel for schedule(static)
#endif
  for (i = 0; i < spp; i++) {
    av[i] = 0.0;
    oc[i] = 1;
  }
  /* Enter the main cycle */
  if (njoin)
    iter = spp - 3;
  else
    iter = spp - 1;
  for (nc = 1; nc <= iter; nc++) {
#ifdef OPENMP_ENABLED
    #pragma omp parallel for private(i) schedule(static)
#endif
    for (j = 2; j <= spp; j++) {
      for (i = 0; i <= j - 2; i++)
        x[j - 1][i] = x[i][j - 1];
    }
    tmin = DBL_MAX;
    /* Compute sij and minimize */
    if (njoin) {     /* many revisions by Y. Ina from here ... */
#ifdef OPENMP_ENABLED
      #pragma omp parallel for schedule(static)
#endif
      for (i = 0; i < spp; i++)
        R[i] = 0.0;
      for (ja = 2; ja <= spp; ja++) {
        jj = enterorder[ja - 1];
        if (cluster[jj - 1] != NULL) {
          for (ia = 0; ia <= ja - 2; ia++) {
            ii = enterorder[ia];
            if (cluster[ii - 1] != NULL) {
              R[ii - 1] += x[ii - 1][jj - 1];
              R[jj - 1] += x[ii - 1][jj - 1];
            }
          }
        }
      }
    } /* ... to here */
    for (ja = 2; ja <= spp; ja++) {
      jj = enterorder[ja - 1];
      if (cluster[jj - 1] != NULL) {
        for (ia = 0; ia <= ja - 2; ia++) {
          ii = enterorder[ia];
          if (cluster[ii - 1] != NULL) {
            if (njoin) {
              total = fotu2 * x[ii - 1][jj - 1] - R[ii - 1] - R[jj - 1];
               /* this statement part of revisions by Y. Ina */
            } else
              total = x[ii - 1][jj - 1];
            if (total < tmin) {
              tmin = total;
              mini = ii;
              minj = jj;
            }
          }
        }
      }
    }
    /* compute lengths and print */
    if (njoin) {
      dio = 0.0;
      djo = 0.0;
#ifdef OPENMP_ENABLED
      #pragma omp parallel for reduction(+:dio,djo) schedule(static)
#endif
      for (i = 0; i < spp; i++) {
        dio += x[i][mini - 1];
        djo += x[i][minj - 1];
      }
      dmin = x[mini - 1][minj - 1];
      dio = (dio - dmin) / fotu2;
      djo = (djo - dmin) / fotu2;
      bi = (dmin + dio - djo) * 0.5;
      bj = dmin - bi;
      bi -= av[mini - 1];
      bj -= av[minj - 1];
    } else {
      bi = x[mini - 1][minj - 1] / 2.0 - av[mini - 1];
      bj = x[mini - 1][minj - 1] / 2.0 - av[minj - 1];
      av[mini - 1] += bi;
    }

    if (progress) {
    /*  printf("Cycle %3ld: ", iter - nc + 1);*/
      if (njoin)
        upgma_nodelabel((boolean)(av[mini - 1] > 0.0));
      else
        upgma_nodelabel((boolean)(oc[mini - 1] > 1.0));
     /* printf(" %ld (%10.5f) joins ", mini, bi);*/
      if (njoin)
        upgma_nodelabel((boolean)(av[minj - 1] > 0.0));
      else
        upgma_nodelabel((boolean)(oc[minj - 1] > 1.0));
      /*printf(" %ld (%10.5f)\n", minj, bj);*/
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    hookup(curtree.nodep[nextnode - 1]->next, cluster[mini - 1]);
    hookup(curtree.nodep[nextnode - 1]->next->next, cluster[minj - 1]);
    cluster[mini - 1]->v = bi;
    cluster[minj - 1]->v = bj;
    cluster[mini - 1]->back->v = bi;
    cluster[minj - 1]->back->v = bj;
    cluster[mini - 1] = curtree.nodep[nextnode - 1];
    cluster[minj - 1] = NULL;
    nextnode++;
    if (njoin)
      av[mini - 1] = dmin * 0.5;
    /* re-initialization */
    fotu2 -= 1.0;
#ifdef OPENMP_ENABLED
    #pragma omp parallel for private(da) schedule(dynamic)
#endif
    for (j = 0; j < spp; j++) {
      if (cluster[j] != NULL) {
        if (njoin) {
          da = (x[mini - 1][j] + x[minj - 1][j]) * 0.5;
          if (mini - j - 1 < 0)
            x[mini - 1][j] = da;
          if (mini - j - 1 > 0)
            x[j][mini - 1] = da;
        } else {
          da = x[mini - 1][j] * oc[mini - 1] + x[minj - 1][j] * oc[minj - 1];
          da /= oc[mini - 1] + oc[minj - 1];
          x[mini - 1][j] = da;
          x[j][mini - 1] = da;
        }
      }
    }
#ifdef OPENMP_ENABLED
    #pragma omp parallel for schedule(static)
#endif
    for (j = 0; j < spp; j++) {
      x[minj - 1][j] = 0.0;
      x[j][minj - 1] = 0.0;
    }
    oc[mini - 1] += oc[minj - 1];
  }
  /* the last cycle */
  nude = 1;
  for (i = 1; i <= spp; i++) {
    if (cluster[i - 1] != NULL) {
      el[nude - 1] = i;
      nude++;
    }
  }
  if (!njoin) {
    curtree.start = cluster[el[0] - 1];
    curtree.start->back = NULL;
    free(av);
    free(oc);
    return;
  }
  bi = (x[el[0] - 1][el[1] - 1] + x[el[0] - 1][el[2] - 1] - x[el[1] - 1]
        [el[2] - 1]) * 0.5;
  bj = x[el[0] - 1][el[1] - 1] - bi;
  bk = x[el[0] - 1][el[2] - 1] - bi;
  bi -= av[el[0] - 1];
  bj -= av[el[1] - 1];
  bk -= av[el[2] - 1];
  if (progress) {
    /*printf("last cycle:\n");*/
    putchar(' ');
    upgma_nodelabel((boolean)(av[el[0] - 1] > 0.0));
    /*printf(" %ld  (%10.5f) joins ", el[0], bi);*/
    upgma_nodelabel((boolean)(av[el[1] - 1] > 0.0));
    /*printf(" %ld  (%10.5f) joins ", el[1], bj);*/
    upgma_nodelabel((boolean)(av[el[2] - 1] > 0.0));
   /* printf(" %ld  (%10.5f)\n", el[2], bk);*/
#ifdef WIN32
    phyFillScreenColor();
#endif
  }
  hookup(curtree.nodep[nextnode - 1], cluster[el[0] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next, cluster[el[1] - 1]);
  hookup(curtree.nodep[nextnode - 1]->next->next, cluster[el[2] - 1]);
  cluster[el[0] - 1]->v = bi;
  cluster[el[1] - 1]->v = bj;
  cluster[el[2] - 1]->v = bk;
  cluster[el[0] - 1]->back->v = bi;
  cluster[el[1] - 1]->back->v = bj;
  cluster[el[2] - 1]->back->v = bk;
  curtree.start = cluster[el[0] - 1]->back;
  free(av);
  free(oc);
  free(R);
}  /* upgma_jointree */


void upgma_maketree()
{
  /* construct the tree */
  long i ;

#ifdef OPENMP_ENABLED
  init_upgma_parallel();
#endif

  inputdata(replicates, printdata, lower, upper, x, reps);
  if (njoin && (spp < 3)) {
    printf("\nERROR: Neighbor-Joining runs must have at least 3 species\n\n");
    exxit(-1);
  }
 /* if (progress)
    putchar('\n');*/
  if (ith == 1)
    setuptree(&curtree, nonodes2 + 1);
  for (i = 1; i <= spp; i++)
    enterorder[i - 1] = i;
  if (jumble)
    randumize(seed, enterorder);
  for (i = 0; i < spp; i++)
    cluster[i] = curtree.nodep[i];
  upgma_jointree();
  if (njoin)
    curtree.start = curtree.nodep[outgrno - 1]->back;
  printree(curtree.start, treeprint, njoin, (boolean)(!njoin));
  if (treeprint)
    upgma_summarize();
  if (trout) {
    col = 0;
    if (njoin)
      treeout(curtree.start, &col, 0.43429448222, njoin, curtree.start);
    else
      curtree.root = curtree.start,
      treeoutr(curtree.start,&col,&curtree);
  }
 /* if (progress) {
    printf("\nOutput written on file \"%s\"\n\n", outfilename);
    if (trout)
      printf("Tree written on file \"%s\"\n\n", outtreename);
  }*/
}  /* upgma_maketree */


int upgma_build_tree(const char *path_name_infile, const char *path_name_outfile)
{  /* main program */
  init();
  char *outfile_path_name=malloc(strlen(path_name_outfile)+strlen("_upgma_outfile")+1);
  char *outtree_path_name=malloc(strlen(path_name_outfile)+strlen("_upgma_outtree")+1);
  strcpy(outfile_path_name,path_name_outfile);
  strcat(outfile_path_name,"_upgma_outfile");
  strcpy(outtree_path_name,path_name_outfile);
  strcat(outtree_path_name,"_upgma_outtree");
  openfile(&infile,path_name_infile,"input file", "r",infilename);
  openfile(&outfile,outfile_path_name,"output file", "w",outfilename);
  ibmpc = IBMCRT;
  ansi = ANSICRT;
  mulsets = false;
  datasets = 1;
  upgma_doinit();
  if (trout)
    openfile(&outtree,outtree_path_name,"output tree file", "w",outtreename);
  ith = 1;
  while (ith <= datasets) {
    if (datasets > 1) {
      fprintf(outfile, "Data set # %ld:\n",ith);
     /* if (progress)
        printf("Data set # %ld:\n",ith);*/
    }
    upgma_inputoptions();
    upgma_maketree();
    if (eoln(infile) && (ith < datasets))
      scan_eoln(infile);
    ith++;
  }
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
  upgma_freerest();
  freetree(&curtree.nodep, nonodes2+1);
  free(outfile_path_name);
  free(outtree_path_name);
#ifdef MAC
  fixmacfile(outfilename);
  fixmacfile(outtreename);
#endif
 /* printf("Done.\n\n");*/
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
  return 0;
}





