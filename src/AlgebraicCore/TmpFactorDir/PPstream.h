#ifndef PPstream_h
#define PPstream_h

struct PPstream_struct
{
  int not_yet_started;
  int nvars;
  int *last_pp;
  int (*cmp)(const void*, const void*);
  int frontier_nelems, frontier_max_sz;
  int **frontier;
  int avoid_last;
  int avoid_nelems, avoid_max_sz;
  int **avoid;
};

typedef struct PPstream_struct *PPstream;

/***************************************************************************/
PPstream PPstream_ctor(int nvars, int (*cmp)(const void*, const void*));
void PPstream_dtor(PPstream self);
int *PPstream_next(PPstream self);
void PPstream_avoid_last(PPstream self);

#endif
