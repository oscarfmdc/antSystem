#include "btsp.h"
#include "tsp_instance.h"
#include "mymalloc.h"
#include <string.h>

t_number *** bTSP_wdist = NULL;
int *** bTSP_wcandlist = NULL;

t_btsp_instance btsp_read_tsplib_instance (const char * filename)
{
    FILE * tsp_file;
    t_tsp_instance instance1, instance2;
    t_btsp_instance btsp_instance;

    tsp_file = fopen (filename, "r");
    if (tsp_file == NULL) {
        perror (filename);
        exit (1);
    }
    read_tsplib_instance (tsp_file, &instance1);
    read_tsplib_instance (tsp_file, &instance2);
    if (instance1.n != instance2.n) {
        fprintf (stderr, "%s: error reading `%s', different instance sizes\n",
                 __FUNCTION__, filename);
        exit (1);
    }
    /* FIXME: there is some memory leaked in instance1 and
       instance2.  */
    btsp_instance.size = instance1.n;
    btsp_instance.distance1 = instance1.distance;
    btsp_instance.distance2 = instance2.distance;
    return btsp_instance;
}

static int
read_problem_size(FILE *input)
{
    char  buf [1000];
    char     *p;
    int size;

    /* FIXME: this could perhaps be simplified by using fscanf (" %d ").  */
    /* Read empty line.  */
    fgets((char*) &buf, 13, input);
    /* Read integer.  */
    fscanf(input, "%d", &size);
    /* Read empty line. */
    fgets((char*) &buf, 1000, input);
    if ((p=strstr((char*)(&buf), "\n"))!=NULL) {    
        p++;
    }
    return size;
}

int btsp_read_instance(const char * filename)
{
    FILE * tsp_file;

    tsp_file = fopen (filename, "r");
    if (tsp_file == NULL) {
        perror (filename);
        exit (1);
    }

    btsp_instance.size = read_problem_size (tsp_file);

    DEBUG2 (fprintf (stderr, "%s:Problem Size: %d\n",
                     __FUNCTION__, btsp_instance.size));

    btsp_instance.distance1 = matrix_number_read (tsp_file, btsp_instance.size, btsp_instance.size);
    btsp_instance.distance2 = matrix_number_read (tsp_file, btsp_instance.size, btsp_instance.size);

    fclose (tsp_file);

    return btsp_instance.size;
}


