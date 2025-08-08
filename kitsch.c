/* version 3.696.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, and Andrew Keeffe.

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

#ifdef _OPENMP
#include <omp.h>
#endif

#define epsilonk         0.000001   /* a very small but not too small number */

#ifndef OLDC
/* function prototypes */
void   init_parallel_kitsch(void);
void   kitsch_getoptions(int tree_type);
void   kitsch_doinit(int tree_type);
void   kitsch_inputoptions(void);
void   getinput(void);
void   input_data(void);
void   input_data_parallel(void);
void   add(node *, node *, node *);
void   kitsch_remove(node **, node **);
void   scrunchtraverse(node *, node **, double *);
void   combine(node *, node *);
void   scrunch(node *);

void   kitsch_secondtraverse(node *, node *, node *, node *, long, long,
                long , double *);
void   firstraverse(node *, node *, double *);
void   sumtraverse(node *, double *);
void   kitsch_evaluate(node *);
void   tryadd(node *, node **, node **);
void   addpreorder(node *, node *, node *);
void   tryrearr(node *, node **, boolean *);
void   repreorder(node *, node **, boolean *);
void   kitsch_rearrange(node **);

void   dtraverse(node *);
void   kitsch_describe(void);
void   kitsch_copynode(node *, node *);
void   kitsch_copynode_parallel(node **, node **, long);
void   kitsch_copy_(tree *, tree *);
void   kitsch_maketree(void);

/* function prototypes */
#endif

#include "kitsch.h"

Char kitsch_infilename[FNMLNGTH], kitsch_outfilename[FNMLNGTH], kitsch_intreename[FNMLNGTH], kitsch_outtreename[FNMLNGTH];
long nonodes, numtrees, col, kitsch_datasets, ith, njumble, jumb;
/*   numtrees is used by usertree option part of kitsch_maketree */
long kitsch_inseed;
tree curtree, bestree;   /* pointers to all nodes in tree */
boolean minev, jumble, usertree, lower, upper, negallowed, replicates, trout,
        printdata, progress, treeprint, kitsch_mulsets, kitsch_firstset;
longer seed;
double power;
long *enterorder;
/* Local variables for kitsch_maketree, propagated globally for C version: */
  long examined;
  double like, bestyet;
  node *there;
  boolean *names;
  Char kitsch_ch;
  char *kitsch_progname;
double trweight; /* to make treeread happy */
boolean goteof, haslengths, lengths;  /* ditto ... */

/* Parallelization variables */
int num_threads_kitsch = 1;
int max_threads_kitsch = 1;


void init_parallel_kitsch()
{
  /* Initialize parallel processing with automatic thread detection for Kitsch */
#ifdef _OPENMP
  max_threads_kitsch = omp_get_max_threads();
  /* Use optimal number of threads: typically number of CPU cores */
  num_threads_kitsch = (max_threads_kitsch > 1) ? max_threads_kitsch : 1;
  omp_set_num_threads(num_threads_kitsch);
  
  if (progress) {
    printf("Kitsch: Parallel processing enabled with %d threads\n", num_threads_kitsch);
  }
#else
  num_threads_kitsch = 1;
  max_threads_kitsch = 1;
  if (progress) {
    printf("Kitsch: Sequential processing (OpenMP not available)\n");
  }
#endif
}  /* init_parallel_kitsch */


void kitsch_getoptions(int tree_type)
{
  /* interactively set options */
  long inseed0, loopcount;
  Char ch;

  minev = false;
  jumble = false;
  njumble = 1;
  lower = false;
  negallowed = false;
  power = 2.0;
  replicates = false;
  upper = false;
  usertree = false;
  trout = true;
  printdata = false;
  progress = true;
  treeprint = true;
  loopcount = 0;

  if(tree_type > 0){
        minev = !minev;
        negallowed = true;
  }

}  /* kitsch_getoptions */


void kitsch_doinit(int tree_type)
{
  /* initializes variables */

  inputnumbers2(&spp, &nonodes, 1);
  kitsch_getoptions(tree_type);
  alloctree(&curtree.nodep, nonodes);
  allocd(nonodes, curtree.nodep);
  allocw(nonodes, curtree.nodep);
  if (!usertree && njumble > 1) {
    alloctree(&bestree.nodep, nonodes);
    allocd(nonodes, bestree.nodep);
    allocw(nonodes, bestree.nodep);
  }
  nayme = (naym *)Malloc(spp*sizeof(naym));
  enterorder = (long *)Malloc(spp*sizeof(long));
}  /* kitsch_doinit */


void kitsch_inputoptions()
{
  /* print options information */
  if (!kitsch_firstset)
    samenumsp2(ith);
  fprintf(outfile, "\nFitch-Margoliash method ");
  fprintf(outfile, "with contemporary tips, version %s\n\n",VERSION);
  if (minev)
    fprintf(outfile, "Minimum evolution method option\n\n");
  fprintf(outfile, "                  __ __             2\n");
  fprintf(outfile, "                  \\  \\   (Obs - Exp)\n");
  fprintf(outfile, "Sum of squares =  /_ /_  ------------\n");
  fprintf(outfile, "                               ");
  if (power == (long)power)
    fprintf(outfile, "%2ld\n", (long)power);
  else
    fprintf(outfile, "%4.1f\n", power);
  fprintf(outfile, "                   i  j      Obs\n\n");
  fprintf(outfile, "negative branch lengths");
  if (!negallowed)
    fprintf(outfile, " not");
  fprintf(outfile, " allowed\n\n");
}  /* kitsch_inputoptions */


void getinput()
{
  /* reads the input data */
  kitsch_inputoptions();
}  /* getinput */


void input_data()
{
  /* read in distance matrix */
  long i, j, k, columns, n;
  boolean skipit, skipother;
  double x;
  columns = replicates ? 4 : 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  setuptree(&curtree, nonodes);
  if (!usertree && njumble > 1)
    setuptree(&bestree, nonodes);
  for (i = 0; i < (spp); i++) {
    curtree.nodep[i]->d[i] = 0.0;
    curtree.nodep[i]->w[i] = 0.0;
    curtree.nodep[i]->weight = 0.0;
    scan_eoln(infile);
    initname(i);
    for (j = 1; j <= (spp); j++) {
      skipit = ((lower && j >= i + 1) || (upper && j <= i + 1));
      skipother = ((lower && i + 1 >= j) || (upper && i + 1 <= j));
      if (!skipit) {
        if (eoln(infile))
          scan_eoln(infile);
        fscanf(infile, "%lf", &x);
        curtree.nodep[i]->d[j - 1] = x;
        if (replicates) {
          if (eoln(infile))
            scan_eoln(infile);
          fscanf(infile, "%ld", &n);
        } else
          n = 1;
        if (n > 0 && x < 0) {
          printf("NEGATIVE DISTANCE BETWEEN SPECIES%5ld AND %5ld\n",
                 i + 1, j);
          exxit(-1);
        }
        curtree.nodep[i]->w[j - 1] = n;
        if (skipother) {
          curtree.nodep[j - 1]->d[i] = curtree.nodep[i]->d[j - 1];
          curtree.nodep[j - 1]->w[i] = curtree.nodep[i]->w[j - 1];
        }
        if ((i == j) && (fabs(curtree.nodep[i-1]->d[j-1]) > 0.000000001)) {
       printf("\nERROR: diagonal element of row %ld of distance matrix ", i+2);
          printf("is not zero.\n");
          printf("       Is it a distance matrix?\n\n");
          exxit(-1);
        }
        if ((j < i) && (fabs(curtree.nodep[i]->d[j-1]-curtree.nodep[j-1]->d[i])
             > 0.000000001)) {
          printf("ERROR: distance matrix is not symmetric:\n");
          printf("       (%ld,%ld) element and (%ld,%ld) element are unequal.\n",
            i+1, j+1, j+1, i+1);
          printf("       They are %10.6f and %10.6f, respectively.\n",
                  curtree.nodep[i]->d[j-1], curtree.nodep[j]->d[i-1]);
          printf("       Is it a distance matrix?\n\n");
          exxit(-1);
        }
      }
    }
  }
  scan_eoln(infile);
  if (printdata) {
    for (i = 0; i < (spp); i++) {
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
      putc(' ', outfile);
      for (j = 1; j <= (spp); j++) {
        fprintf(outfile, "%10.5f", curtree.nodep[i]->d[j - 1]);
        if (replicates)
          fprintf(outfile, " (%3ld)", (long)curtree.nodep[i]->w[j - 1]);
        if (j % columns == 0 && j < spp) {
          putc('\n', outfile);
          for (k = 1; k <= nmlngth + 1; k++)
            putc(' ', outfile);
        }
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  for (i = 0; i < (spp); i++) {
    for (j = 0; j < (spp); j++) {
      if (i + 1 != j + 1) {
        if (curtree.nodep[i]->d[j] < epsilonk)
          curtree.nodep[i]->d[j] = epsilonk;
        curtree.nodep[i]->w[j] /= exp(power * log(curtree.nodep[i]->d[j]));
      }
    }
  }
}  /* inputdata */


void input_data_parallel()
{
  /* read in distance matrix with parallelized post-processing */
  long i, j, k, columns, n;
  boolean skipit, skipother;
  double x;
  
  /* First part same as original input_data */
  columns = replicates ? 4 : 6;
  if (printdata) {
    fprintf(outfile, "\nName                       Distances");
    if (replicates)
      fprintf(outfile, " (replicates)");
    fprintf(outfile, "\n----                       ---------");
    if (replicates)
      fprintf(outfile, "-------------");
    fprintf(outfile, "\n\n");
  }
  setuptree(&curtree, nonodes);
  if (!usertree && njumble > 1)
    setuptree(&bestree, nonodes);
  
  /* Sequential data reading (I/O cannot be easily parallelized) */
  for (i = 0; i < (spp); i++) {
    curtree.nodep[i]->d[i] = 0.0;
    curtree.nodep[i]->w[i] = 0.0;
    curtree.nodep[i]->weight = 0.0;
    scan_eoln(infile);
    initname(i);
    for (j = 1; j <= (spp); j++) {
      skipit = ((lower && j >= i + 1) || (upper && j <= i + 1));
      skipother = ((lower && i + 1 >= j) || (upper && i + 1 <= j));
      if (!skipit) {
        if (eoln(infile))
          scan_eoln(infile);
        fscanf(infile, "%lf", &x);
        curtree.nodep[i]->d[j - 1] = x;
        if (replicates) {
          if (eoln(infile))
            scan_eoln(infile);
          fscanf(infile, "%ld", &n);
        } else
          n = 1;
        if (n > 0 && x < 0) {
          printf("NEGATIVE DISTANCE BETWEEN SPECIES%5ld AND %5ld\n",
                 i + 1, j);
          exxit(-1);
        }
        curtree.nodep[i]->w[j - 1] = n;
        if (skipother) {
          curtree.nodep[j - 1]->d[i] = curtree.nodep[i]->d[j - 1];
          curtree.nodep[j - 1]->w[i] = curtree.nodep[i]->w[j - 1];
        }
        if ((i == j) && (fabs(curtree.nodep[i-1]->d[j-1]) > 0.000000001)) {
       printf("\nERROR: diagonal element of row %ld of distance matrix ", i+2);
          printf("is not zero.\n");
          printf("       Is it a distance matrix?\n\n");
          exxit(-1);
        }
        if ((j < i) && (fabs(curtree.nodep[i]->d[j-1]-curtree.nodep[j-1]->d[i])
             > 0.000000001)) {
          printf("ERROR: distance matrix is not symmetric:\n");
          printf("       (%ld,%ld) element and (%ld,%ld) element are unequal.\n",
            i+1, j+1, j+1, i+1);
          printf("       They are %10.6f and %10.6f, respectively.\n",
                  curtree.nodep[i]->d[j-1], curtree.nodep[j]->d[i-1]);
          printf("       Is it a distance matrix?\n\n");
          exxit(-1);
        }
      }
    }
  }
  scan_eoln(infile);
  
  /* Print data if requested (sequential) */
  if (printdata) {
    for (i = 0; i < (spp); i++) {
      for (j = 0; j < nmlngth; j++)
        putc(nayme[i][j], outfile);
      putc(' ', outfile);
      for (j = 1; j <= (spp); j++) {
        fprintf(outfile, "%10.5f", curtree.nodep[i]->d[j - 1]);
        if (replicates)
          fprintf(outfile, " (%3ld)", (long)curtree.nodep[i]->w[j - 1]);
        if (j % columns == 0 && j < spp) {
          putc('\n', outfile);
          for (k = 1; k <= nmlngth + 1; k++)
            putc(' ', outfile);
        }
      }
      putc('\n', outfile);
    }
    putc('\n', outfile);
  }
  
  /* Parallelize the post-processing loops */
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) private(i, j) if(spp > 50)
#endif
  for (i = 0; i < (spp); i++) {
    for (j = 0; j < (spp); j++) {
      if (i + 1 != j + 1) {
        if (curtree.nodep[i]->d[j] < epsilonk)
          curtree.nodep[i]->d[j] = epsilonk;
        curtree.nodep[i]->w[j] /= exp(power * log(curtree.nodep[i]->d[j]));
      }
    }
  }
}  /* input_data_parallel */


void add(node *below, node *newtip, node *newfork)
{
  /* inserts the nodes newfork and its left descendant, newtip,
     to the tree.  below becomes newfork's right descendant */
  if (below != curtree.nodep[below->index - 1])
    below = curtree.nodep[below->index - 1];
  if (below->back != NULL)
    below->back->back = newfork;
  newfork->back = below->back;
  below->back = newfork->next->next;
  newfork->next->next->back = below;
  newfork->next->back = newtip;
  newtip->back = newfork->next;
  if (curtree.root == below)
    curtree.root = newfork;
  curtree.root->back = NULL;
}  /* add */


void kitsch_remove(node **item, node **fork)
{
  /* removes nodes item and its ancestor, fork, from the tree.
     the new descendant of fork's ancestor is made to be
     fork's second descendant (other than item).  Also
     returns pointers to the deleted nodes, item and fork */
  node *p, *q;

  if ((*item)->back == NULL) {
    *fork = NULL;
    return;
  }
  *fork = curtree.nodep[(*item)->back->index - 1];
  if (curtree.root == *fork) {
    if (*item == (*fork)->next->back)
      curtree.root = (*fork)->next->next->back;
    else
      curtree.root = (*fork)->next->back;
  }
  p = (*item)->back->next->back;
  q = (*item)->back->next->next->back;
  if (p != NULL)
    p->back = q;
  if (q != NULL)
    q->back = p;
  (*fork)->back = NULL;
  p = (*fork)->next;
  while (p != *fork) {
    p->back = NULL;
    p = p->next;
  }
  (*item)->back = NULL;
}  /* remove */


void scrunchtraverse(node *u, node **closest, double *tmax)
{
  /* traverse to find closest node to the current one */
  if (!u->sametime) {
    if (u->t > *tmax) {
      *closest = u;
      *tmax = u->t;
    }
    return;
  }
  u->t = curtree.nodep[u->back->index - 1]->t;
  if (!u->tip) {
    scrunchtraverse(u->next->back, closest,tmax);
    scrunchtraverse(u->next->next->back, closest,tmax);
  }
}  /* scrunchtraverse */


void combine(node *a, node *b)
{
  /* put node b into the set having the same time as a */
  if (a->weight + b->weight <= 0.0)
    a->t = 0.0;
  else
    a->t = (a->t * a->weight + b->t * b->weight) / (a->weight + b->weight);
  a->weight += b->weight;
  b->sametime = true;
}  /* combine */


void scrunch(node *s)
{
  /* see if nodes can be combined to prevent negative lengths */
  double tmax;
  node *closest;
  boolean found;

  closest = NULL;
  tmax = -1.0;
  do {
    if (!s->tip) {
      scrunchtraverse(s->next->back, &closest,&tmax);
      scrunchtraverse(s->next->next->back, &closest,&tmax);
    }
    found = (tmax > s->t);
    if (found)
      combine(s, closest);
    tmax = -1.0;
  } while (found);
}  /* scrunch */


void kitsch_secondtraverse(node *a, node *q, node *u, node *v, long i, long j,
                        long k, double *sum)
{
  /* recalculate distances, add to sum */
  long l;
  double wil, wjl, wkl, wli, wlj, wlk, TEMP;

  if (!(a->processed || a->tip)) {
    kitsch_secondtraverse(a->next->back, q,u,v,i,j,k,sum);
    kitsch_secondtraverse(a->next->next->back, q,u,v,i,j,k,sum);
    return;
  }
  if (!(a != q && a->processed))
    return;
  l = a->index;
  wil = u->w[l - 1];
  wjl = v->w[l - 1];
  wkl = wil + wjl;
  wli = a->w[i - 1];
  wlj = a->w[j - 1];
  wlk = wli + wlj;
  q->w[l - 1] = wkl;
  a->w[k - 1] = wlk;
  if (wkl <= 0.0)
    q->d[l - 1] = 0.0;
  else
    q->d[l - 1] = (wil * u->d[l - 1] + wjl * v->d[l - 1]) / wkl;
  if (wlk <= 0.0)
    a->d[k - 1] = 0.0;
  else
    a->d[k - 1] = (wli * a->d[i - 1] + wlj * a->d[j - 1]) / wlk;
  if (minev)
    return;
  
  /* These computations can cause race conditions in parallel execution */
  /* but the contribution to sum needs to be accumulated thread-safely */
  double local_contribution = 0.0;
  if (wkl > 0.0) {
    TEMP = u->d[l - 1] - v->d[l - 1];
    local_contribution += wil * wjl / wkl * (TEMP * TEMP);
  }
  if (wlk > 0.0) {
    TEMP = a->d[i - 1] - a->d[j - 1];
    local_contribution += wli * wlj / wlk * (TEMP * TEMP);
  }
  
  /* Thread-safe accumulation */
#ifdef _OPENMP
  #pragma omp atomic
#endif
  (*sum) += local_contribution;
}  /* kitsch_secondtraverse */


void firstraverse(node *q_, node *r, double *sum)
{  /* firsttraverse                              */
   /* go through tree calculating branch lengths */
  node *q;
  long i, j, k;
  node *u, *v;

  q = q_;
  if (q == NULL)
    return;
  q->sametime = false;
  if (!q->tip) {
    firstraverse(q->next->back, r,sum);
    firstraverse(q->next->next->back, r,sum);
  }
  q->processed = true;
  if (q->tip)
    return;
  u = q->next->back;
  v = q->next->next->back;
  i = u->index;
  j = v->index;
  k = q->index;
  if (u->w[j - 1] + v->w[i - 1] <= 0.0)
    q->t = 0.0;
  else
    q->t = (u->w[j - 1] * u->d[j - 1] +  v->w[i - 1] * v->d[i - 1]) /
             (2.0 * (u->w[j - 1] + v->w[i - 1]));
  q->weight = u->weight + v->weight + u->w[j - 1] + v->w[i - 1];
  if (!negallowed)
    scrunch(q);
  u->v = q->t - u->t;
  v->v = q->t - v->t;
  u->back->v = u->v;
  v->back->v = v->v;
  kitsch_secondtraverse(r,q,u,v,i,j,k,sum);
}  /* firstraverse */


void sumtraverse(node *q, double *sum)
{
  /* traverse to finish computation of sum of squares */
  long i, j;
  node *u, *v;
  double TEMP, TEMP1;

  if (minev && (q != curtree.root))
    *sum += q->v;
  if (q->tip)
    return;
  sumtraverse(q->next->back, sum);
  sumtraverse(q->next->next->back, sum);
  if (!minev) {
    u = q->next->back;
    v = q->next->next->back;
    i = u->index;
    j = v->index;
    TEMP = u->d[j - 1] - 2.0 * q->t;
    TEMP1 = v->d[i - 1] - 2.0 * q->t;
    (*sum) += u->w[j - 1] * (TEMP * TEMP) + v->w[i - 1] * (TEMP1 * TEMP1);
  }
}  /* sumtraverse */


void kitsch_evaluate(node *r)
{
  /* fill in times and kitsch_evaluate sum of squares for tree */
  double sum;
  long i;
  sum = 0.0;
  for (i = 0; i < (nonodes); i++)
    curtree.nodep[i]->processed = curtree.nodep[i]->tip;
  firstraverse(r, r,&sum);
  sumtraverse(r, &sum);
  examined++;
  if (replicates && (lower || upper))
    sum /= 2;
  like = -sum;
}  /* kitsch_evaluate */


void tryadd(node *p, node **item, node **nufork)
{
  /* temporarily adds one fork and one tip to the tree.
     if the location where they are added yields greater
     "likelihood" than other locations tested up to that
     time, then keeps that location as there */
  add(p, *item, *nufork);
  kitsch_evaluate(curtree.root);
  if (like > bestyet) {
    bestyet = like;
    there = p;
  }
  kitsch_remove(item, nufork);
}  /* tryadd */


void addpreorder(node *p, node *item, node *nufork)
{
  /* traverses a binary tree, calling PROCEDURE tryadd
     at a node before calling tryadd at its descendants */
  if (p == NULL)
    return;
  tryadd(p, &item,&nufork);
  if (!p->tip) {
    addpreorder(p->next->back, item, nufork);
    addpreorder(p->next->next->back, item, nufork);
  }
}  /* addpreorder */


void tryrearr(node *p, node **r, boolean *success)
{
  /* evaluates one rearrangement of the tree.
     if the new tree has greater "likelihood" than the old
     one sets success := TRUE and keeps the new tree.
     otherwise, restores the old tree */
  node *frombelow, *whereto, *forknode;
  double oldlike;

  if (p->back == NULL)
    return;
  forknode = curtree.nodep[p->back->index - 1];
  if (forknode->back == NULL)
    return;
  oldlike = like;
  if (p->back->next->next == forknode)
    frombelow = forknode->next->next->back;
  else
    frombelow = forknode->next->back;
  whereto = forknode->back;
  kitsch_remove(&p, &forknode);
  add(whereto, p, forknode);
  if ((*r)->back != NULL)
    *r = curtree.nodep[(*r)->back->index - 1];
  kitsch_evaluate(*r);
  if (like - oldlike > LIKE_EPSILON) {
    bestyet = like;
    *success = true;
    return;
  }
  kitsch_remove(&p, &forknode);
  add(frombelow, p, forknode);
  if ((*r)->back != NULL)
    *r = curtree.nodep[(*r)->back->index - 1];
  like = oldlike;
}  /* tryrearr */


void repreorder(node *p, node **r, boolean *success)
{
  /* traverses a binary tree, calling PROCEDURE tryrearr
     at a node before calling tryrearr at its descendants */
  if (p == NULL)
    return;
  tryrearr(p,r,success);
  if (!p->tip) {
    repreorder(p->next->back,r,success);
    repreorder(p->next->next->back,r,success);
  }
}  /* repreorder */


void kitsch_rearrange(node **r_)
{
  /* traverses the tree (preorder), finding any local
     rearrangement which decreases the number of steps.
     if traversal succeeds in increasing the tree's
     "likelihood", PROCEDURE kitsch_rearrange runs traversal again */
  node **r;
  boolean success;

  r = r_;
  success = true;
  while (success) {
    success = false;
    repreorder(*r,r,&success);
  }
}  /* kitsch_rearrange */


void dtraverse(node *q)
{
  /* print table of lengths etc. */
  long i;

  if (!q->tip)
    dtraverse(q->next->back);
  if (q->back != NULL) {
    fprintf(outfile, "%4ld   ", q->back->index - spp);
    if (q->index <= spp) {
      for (i = 0; i < nmlngth; i++)
        putc(nayme[q->index - 1][i], outfile);
    } else
      fprintf(outfile, "%4ld      ", q->index - spp);
    fprintf(outfile, "%13.5f", curtree.nodep[q->back->index - 1]->t - q->t);
    q->v = curtree.nodep[q->back->index - 1]->t - q->t;
    q->back->v = q->v;
    fprintf(outfile, "%16.5f\n", curtree.root->t - q->t);
  }
  if (!q->tip)
    dtraverse(q->next->next->back);
}  /* dtraverse */


void kitsch_describe()
{
  /* prints table of lengths, times, sum of squares, etc. */
  long i, j;
  double totalnum;
  double TEMP;

  if (!minev)
    fprintf(outfile, "\nSum of squares = %10.3f\n\n", -like);
  else
    fprintf(outfile, "Sum of branch lengths = %10.3f\n\n", -like);
  if ((fabs(power - 2) < 0.01) && !minev) {
    totalnum = 0.0;
    
    /* Parallelize the sum calculation */
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) private(i, j, TEMP) reduction(+:totalnum) if(spp > 50)
#endif
    for (i = 0; i < (spp); i++) {
      for (j = 0; j < (spp); j++) {
        if (i + 1 != j + 1 && curtree.nodep[i]->d[j] > 0.0) {
          TEMP = curtree.nodep[i]->d[j];
          totalnum += curtree.nodep[i]->w[j] * (TEMP * TEMP);
        }
      }
    }
    totalnum -= 2;
    if (replicates && (lower || upper))
      totalnum /= 2;
    fprintf(outfile, "Average percent standard deviation =");
    fprintf(outfile, "%10.5f\n\n", 100 * sqrt(-(like / totalnum)));
  }
  fprintf(outfile, "From     To            Length          Height\n");
  fprintf(outfile, "----     --            ------          ------\n\n");
  dtraverse(curtree.root);
  putc('\n', outfile);
  if (trout) {
    col = 0;
    treeoutr(curtree.root,&col,&curtree);
  }
}  /* kitsch_describe */


void kitsch_copynode(node *c, node *d)
{
  /* make a copy of a node */

  memcpy(d->d, c->d, nonodes*sizeof(double));
  memcpy(d->w, c->w, nonodes*sizeof(double));
  d->t = c->t;
  d->sametime = c->sametime;
  d->weight = c->weight;
  d->processed = c->processed;
  d->xcoord = c->xcoord;
  d->ycoord = c->ycoord;
  d->ymin = c->ymin;
  d->ymax = c->ymax;
}  /* kitsch_copynode */


void kitsch_copynode_parallel(node **src_nodes, node **dst_nodes, long count)
{
  /* parallel version of kitsch_copynode for multiple nodes */
  long i;
#ifdef _OPENMP
  #pragma omp parallel for schedule(static) private(i) if(count > 20)
#endif
  for (i = 0; i < count; i++) {
    if (src_nodes[i] != NULL && dst_nodes[i] != NULL) {
      kitsch_copynode(src_nodes[i], dst_nodes[i]);
    }
  }
}  /* kitsch_copynode_parallel */


void kitsch_copy_(tree *a, tree *b)
{
  /* make a copy of a tree */
  long i, j=0;
  node *p, *q;

  for (i = 0; i < spp; i++) {
    kitsch_copynode(a->nodep[i], b->nodep[i]);
    if (a->nodep[i]->back) {
      if (a->nodep[i]->back == a->nodep[a->nodep[i]->back->index - 1])
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1];
      else if (a->nodep[i]->back
                 == a->nodep[a->nodep[i]->back->index - 1]->next)
        b->nodep[i]->back = b->nodep[a->nodep[i]->back->index - 1]->next;
      else
        b->nodep[i]->back
          = b->nodep[a->nodep[i]->back->index - 1]->next->next;
    }
    else b->nodep[i]->back = NULL;
  }
  for (i = spp; i < nonodes; i++) {
    p = a->nodep[i];
    q = b->nodep[i];
    for (j = 1; j <= 3; j++) {
      kitsch_copynode(p, q);
      if (p->back) {
        if (p->back == a->nodep[p->back->index - 1])
          q->back = b->nodep[p->back->index - 1];
        else if (p->back == a->nodep[p->back->index - 1]->next)
          q->back = b->nodep[p->back->index - 1]->next;
        else
          q->back = b->nodep[p->back->index - 1]->next->next;
      }
      else
        q->back = NULL;
      p = p->next;
      q = q->next;
    }
  }
  b->root = a->root;
}  /* copy */


void kitsch_maketree()
{
  /* constructs a binary tree from the pointers in curtree.nodep.
     adds each node at location which yields highest "likelihood"
     then rearranges the tree for greatest "likelihood" */
  long i, j, which;
  double bestlike, bstlike2=0, gotlike;
  boolean lastrearr;
  node *item, *nufork;

  if (!usertree) {
    if (jumb == 1) {
      /* Initialize parallel processing */
      init_parallel_kitsch();
      input_data_parallel();
      examined = 0;
    }
    
    /* Parallelize enterorder initialization */
#ifdef _OPENMP
    #pragma omp parallel for schedule(static) private(i) if(spp > 100)
#endif
    for (i = 1; i <= (spp); i++)
      enterorder[i - 1] = i;
    if (jumble)
      randumize(seed, enterorder);
    curtree.root = curtree.nodep[enterorder[0] - 1];
    add(curtree.nodep[enterorder[0] - 1], curtree.nodep[enterorder[1] - 1],
        curtree.nodep[spp]);
    if (progress) {
      //printf("Adding species:\n");
      writename(0, 2, enterorder);
#ifdef WIN32
      phyFillScreenColor();
#endif
    }
    for (i = 3; i <= (spp); i++) {
      bestyet = -DBL_MAX;
      item = curtree.nodep[enterorder[i - 1] - 1];
      nufork = curtree.nodep[spp + i - 2];
      addpreorder(curtree.root, item, nufork);
      add(there, item, nufork);
      like = bestyet;
      kitsch_rearrange(&curtree.root);
      kitsch_evaluate(curtree.root);
      examined--;
      if (progress) {
        writename(i - 1, 1, enterorder);
#ifdef WIN32
        phyFillScreenColor();
#endif
      }
      lastrearr = (i == spp);
      if (lastrearr) {
        if (progress) {
          /*
          printf("\nDoing global rearrangements\n");
          printf("  !");
          for (j = 1; j <= (nonodes); j++)
            if ( j % (( nonodes / 72 ) + 1 ) == 0 )
              putchar('-');
          printf("!\n");

          */
#ifdef WIN32
          phyFillScreenColor();
#endif
        }
        bestlike = bestyet;
        do {
          gotlike = bestlike;
          /*if (progress)
            printf("   ");*/
            
          /* Parallelize the global rearrangement loop */
#ifdef _OPENMP
          #pragma omp parallel for schedule(dynamic) private(j, there, item, nufork) \
                  reduction(max:bestyet) if(nonodes > 20)
#endif
          for (j = 0; j < (nonodes); j++) {
            node *local_there, *local_item, *local_nufork;
            double local_bestyet = -DBL_MAX;
            
            local_there = curtree.root;
            local_item = curtree.nodep[j];
            
            if (local_item != curtree.root) {
              /* Thread-safe tree operations */
#ifdef _OPENMP
              #pragma omp critical
#endif
              {
                kitsch_remove(&local_item, &local_nufork);
                local_there = curtree.root;
                addpreorder(curtree.root, local_item, local_nufork);
                add(there, local_item, local_nufork);
              }
              
              /* Update global bestyet in thread-safe manner */
              if (bestyet > local_bestyet) {
                local_bestyet = bestyet;
              }
            }
            
            if (progress) {
             /* if ( j % (( nonodes / 72 ) + 1 ) == 0 )
                putchar('.');
              fflush(stdout);*/
            }
          }
          if (progress) {
            /*putchar('\n');*/
#ifdef WIN32
            phyFillScreenColor();
#endif
          }
        } while (bestlike > gotlike);
        if (njumble > 1) {
          if (jumb == 1 || (jumb > 1 && bestlike > bstlike2)) {
            kitsch_copy_(&curtree, &bestree);
            bstlike2 = bestlike;
          }
        }
      }
    }
    if (njumble == jumb) {
      if (njumble > 1)
        kitsch_copy_(&bestree, &curtree);
      kitsch_evaluate(curtree.root);
      printree(curtree.root, treeprint, false, true);
      kitsch_describe();
    }
  } else {
    /* Initialize parallel processing for user trees */
    init_parallel_kitsch();
    input_data_parallel();
    /* Open in binary: ftell() is broken for UNIX line-endings under WIN32 */
    openfile(&intree,INTREE,"input tree file","rb",kitsch_intreename);
    numtrees = countsemic(&intree);
    if (treeprint)
      fprintf(outfile, "\n\nUser-defined trees:\n\n");
    names = (boolean *)Malloc(spp*sizeof(boolean));
    which = 1;
    while (which <= numtrees ) {
      treeread2 (intree, &curtree.root, curtree.nodep, lengths, &trweight,
                      &goteof, &haslengths, &spp,false,nonodes);
      if (curtree.root->back) {
        printf("Error:  Kitsch cannot read unrooted user trees\n");
        exxit(-1);
      }
      kitsch_evaluate(curtree.root);
      printree(curtree.root, treeprint, false, true);
      kitsch_describe();
      which++;
    }
    FClose(intree);
    free(names);
  }
  if (jumb == njumble && progress) {
   /* printf("\nOutput written to file \"%s\"\n", kitsch_outfilename);
    if (trout)
      printf("\nTree also written onto file \"%s\"\n", kitsch_outtreename);*/
  }
}  /* kitsch_maketree */


int kitsch_build_tree(const char *path_name_infile, const  char *path_name_outfile, int tree_type)
{  /* Fitch-Margoliash criterion with contemporary tips */
  char *outfile_path_name=malloc(strlen(path_name_outfile)+strlen("_kitsch_outfile")+1);
  char *outtree_path_name=malloc(strlen(path_name_outfile)+strlen("_kitsch_outtree")+1);
  strcpy(outfile_path_name,path_name_outfile);
  strcat(outfile_path_name,"_kitsch_outfile");
  strcpy(outtree_path_name,path_name_outfile);
  strcat(outtree_path_name,"_kitsch_outtree");
  if(!(tree_type ==0 || tree_type == 1)){
        printf("\nTree type error !\n");
        return 0;
  }
  init();
  //printf("\n%s\n",path_name_infile);
  /* reads in spp, options, and the data, then calls kitsch_maketree to
     construct the tree */
  openfile(&infile,path_name_infile,"input file","r",kitsch_infilename);
  openfile(&outfile,outfile_path_name,"output file","w",kitsch_outfilename);

  ibmpc = IBMCRT;
  ansi = ANSICRT;
  kitsch_mulsets = false;
  kitsch_firstset = true;
  kitsch_datasets = 1;
  kitsch_doinit(tree_type);
  openfile(&outtree,outtree_path_name,"output tree file","w",kitsch_outtreename);
  for (ith = 1; ith <= kitsch_datasets; ith++) {
    if (kitsch_datasets > 1) {
      fprintf(outfile, "\nData set # %ld:\n",ith);
     /* if (progress)
        printf("\nData set # %ld:\n",ith);*/
    }
    getinput();
    for (jumb = 1; jumb <= njumble; jumb++)
      kitsch_maketree();
    kitsch_firstset = false;
    if (eoln(infile) && (ith < kitsch_datasets))
      scan_eoln(infile);
  }
  FClose(infile);
  FClose(outfile);
  FClose(outtree);
  free(outfile_path_name);
  free(outtree_path_name);
#ifdef MAC
  fixmacfile(kitsch_outfilename);
  fixmacfile(kitsch_outtreename);
#endif
#ifdef WIN32
  phyRestoreConsoleAttributes();
#endif
 /* printf("\nDone.\n\n");*/
  return 0;
}  /* Fitch-Margoliash criterion with contemporary tips */
