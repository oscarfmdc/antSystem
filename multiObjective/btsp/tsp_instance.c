#include "common.h"
#include "tsp_instance.h"
#include "work.h"
#define TRACE DEBUG2

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

struct point {
  double x;
  double y;
};
static struct point *nodeptr;

#define TSPLIB_PI 3.14159265358979323846264

static inline double
dtrunc (double x)
{
    int k;

    k = (int) x;
    x = (double) k;
    return x;
}

/*    
      FUNCTION: the following four functions implement different ways of 
                computing distances for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
*/

static int round_distance (int i, int j) 
/*    
      FUNCTION: compute Euclidean distances between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
    double xd = nodeptr[i].x - nodeptr[j].x;
    double yd = nodeptr[i].y - nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd) + 0.5;

    return (int) r;
}

static int ceil_distance (int i, int j) 
/*    
      FUNCTION: compute ceiling distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
    double xd = nodeptr[i].x - nodeptr[j].x;
    double yd = nodeptr[i].y - nodeptr[j].y;
    double r  = sqrt(xd*xd + yd*yd) + 0.000000001;

    return (int)r;
}

static int geo_distance (int i, int j) 
/*    
      FUNCTION: compute geometric distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: adapted from concorde code
                for the definition of how to compute this distance see TSPLIB
*/
{
    double deg, min;
    double lati, latj, longi, longj;
    double q1, q2, q3;
    int dd;
    double x1 = nodeptr[i].x, x2 = nodeptr[j].x, 
	y1 = nodeptr[i].y, y2 = nodeptr[j].y;

    deg = dtrunc (x1);
    min = x1 - deg;
    lati = TSPLIB_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (x2);
    min = x2 - deg;
    latj = TSPLIB_PI * (deg + 5.0 * min / 3.0) / 180.0;

    deg = dtrunc (y1);
    min = y1 - deg;
    longi = TSPLIB_PI * (deg + 5.0 * min / 3.0) / 180.0;
    deg = dtrunc (y2);
    min = y2 - deg;
    longj = TSPLIB_PI * (deg + 5.0 * min / 3.0) / 180.0;

    q1 = cos (longi - longj);
    q2 = cos (lati - latj);
    q3 = cos (lati + latj);
    dd = (int) (6378.388 * acos (0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    return dd;

}

static int att_distance (int i, int j) 
/*    
      FUNCTION: compute ATT distance between two nodes rounded to next 
                integer for TSPLIB instances
      INPUT:    two node indices
      OUTPUT:   distance between the two nodes
      COMMENTS: for the definition of how to compute this distance see TSPLIB
*/
{
    double xd = nodeptr[i].x - nodeptr[j].x;
    double yd = nodeptr[i].y - nodeptr[j].y;
    double rij = sqrt ((xd * xd + yd * yd) / 10.0);
    double tij = dtrunc (rij);
    int dij;

    if (tij < rij)
        dij = (int) tij + 1;
    else
        dij = (int) tij;
    return dij;
}



static t_number ** 
compute_distances (int n,
                   int  (*distance)(int, int) /* function pointer */)
/*    
      FUNCTION: computes the matrix of all intercity distances
      INPUT:    none
      OUTPUT:   pointer to distance matrix, has to be freed when program stops
*/
{
    int     i, j;
    t_number     **matrix;

    if((matrix = malloc(sizeof(t_number) * n * n +
			sizeof(t_number *) * n	 )) == NULL){
	fprintf(stderr,"Out of memory, exit.");
	exit(1);
    }
    for ( i = 0 ; i < n ; i++ ) {
	matrix[i] = (t_number *)(matrix + n) + i*n;
	for (j = 0  ; j < n ; j++) {
	    matrix[i][j] = distance(i, j);
	}
    }
    return matrix;
}

/* FIXME: this could return an error and fatal could be configured to not
   abort.  */
void
read_tsplib_instance (FILE *tsp_file,
                      t_tsp_instance *tsp_instance) 
/* FIXME: update the description.  */
/*    
      FUNCTION: parse and read instance file
      INPUT:    instance name
      OUTPUT:   list of coordinates for all nodes
      COMMENTS: Instance files have to be in TSPLIB format, otherwise procedure fails
*/
{
    char buf[LINE_BUF_LEN];
    char edge_weight_type[LINE_BUF_LEN];
    int     i, j;
    int n;
    int  (*distance)(int, int) = NULL;

    /* FIXME: this allows writing pass the end of buf. */
    fscanf (tsp_file,"%s", buf);
    while (strcmp("NODE_COORD_SECTION", buf) != 0) {
        /* Parse "NAME:" or "NAME : ".  */
        /* FIXME: this could be done much better by using 
           fscanf (" NAME : "). */
	if (strcmp("NAME :", buf) == 0) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf) );
	    fscanf(tsp_file, "%s", buf);
	    strcpy(tsp_instance->name, buf);
	    TRACE ( printf("%s \n", tsp_instance->name));
	    buf[0]=0;
	}
	else if ( strcmp("NAME:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    strcpy(tsp_instance->name, buf);
	    TRACE ( printf("%s \n", tsp_instance->name));
	    buf[0]=0;
	}
	else if ( strcmp("COMMENT", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf));
	    buf[0]=0;
	}
	else if ( strcmp("COMMENT:", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf));
	    buf[0]=0;
	}
	else if ( strcmp("TYPE", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf));
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf));
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	    buf[0]=0;
	}
	else if ( strcmp("TYPE:", buf) == 0 ) {
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf));
	    if( strcmp("TSP", buf) != 0 ) {
		fprintf(stderr,"\n Not a TSP instance in TSPLIB format !!\n");
		exit(1);
	    }
	    buf[0]=0;
	}
	else if( strcmp("DIMENSION", buf) == 0 ){
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf));;
	    fscanf(tsp_file, "%d", &n);
	    tsp_instance->n = n;
	    TRACE ( printf("%d\n", n));;
	    assert ( n > 2 && n < 6000);
	    buf[0]=0;
	}
	else if ( strcmp("DIMENSION:", buf) == 0 ) {
	    fscanf(tsp_file, "%d", &n);
	    tsp_instance->n = n;
	    TRACE ( printf("%d\n", n));;
	    assert ( n > 2 && n < 6000);
	    buf[0]=0;
	}
	else if( strcmp("DISPLAY_DATA_TYPE", buf) == 0 ){
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf));;
	    buf[0]=0;
	}
	else if ( strcmp("DISPLAY_DATA_TYPE:", buf) == 0 ) {
	    fgets(buf, LINE_BUF_LEN, tsp_file);
	    TRACE ( printf("%s", buf));;
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE", buf) == 0 ){
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s ", buf));;
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf));;
	    if ( strcmp("EUC_2D", buf) == 0 ) {
		distance = round_distance;
	    }
	    else if ( strcmp("CEIL_2D", buf) == 0 ) {
		distance = ceil_distance;
	    }
	    else if ( strcmp("GEO", buf) == 0 ) {
		distance = geo_distance;
	    }
	    else if ( strcmp("ATT", buf) == 0 ) {
		distance = att_distance;
	    }
	    else {
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
                exit (1);
            }
	    strcpy (edge_weight_type, buf);
	    buf[0]=0;
	}
	else if( strcmp("EDGE_WEIGHT_TYPE:", buf) == 0 ){
	    /* set pointer to appropriate distance function; has to be one of 
	       EUC_2D, CEIL_2D, GEO, or ATT. Everything else fails */
	    buf[0]=0;
	    fscanf(tsp_file, "%s", buf);
	    TRACE ( printf("%s\n", buf));
	    if ( strcmp("EUC_2D", buf) == 0 ) {
		distance = round_distance;
	    }
	    else if ( strcmp("CEIL_2D", buf) == 0 ) {
		distance = ceil_distance;
	    }
	    else if ( strcmp("GEO", buf) == 0 ) {
		distance = geo_distance;
	    }
	    else if ( strcmp("ATT", buf) == 0 ) {
		distance = att_distance;
	    }
	    else {
		fprintf(stderr,"EDGE_WEIGHT_TYPE %s not implemented\n",buf);
		exit(1);
	    }
	    strcpy(edge_weight_type, buf);
	    buf[0]=0;
	}
	buf[0]=0;
	fscanf(tsp_file,"%s", buf);
    }


    if( strcmp("NODE_COORD_SECTION", buf) == 0 ){
	TRACE ( printf("found section contaning the node coordinates\n"));
	    }
    else{
	fprintf(stderr,"\n\nSome error ocurred finding start of coordinates from tsp file !!\n");
	exit(1);
    }

    if( (nodeptr = malloc(sizeof(struct point) * n)) == NULL )
	exit(EXIT_FAILURE);
    else {
	for ( i = 0 ; i < n ; i++ ) {
	    fscanf(tsp_file,"%d %lf %lf", &j, &nodeptr[i].x, &nodeptr[i].y );
	}
    }
    TRACE ( printf("number of cities is %d\n",n);
            printf("\n... done\n"));
    /* FIXME: if there is no EDGE_TYPE section, this crashes.  */
    tsp_instance->distance = compute_distances(tsp_instance->n, distance);
    free (nodeptr);
}

t_tsp_instance *
read_tsplib_file (const char * tsp_filename)
/* FIXME: update the description.  */
{
    t_tsp_instance *tsp_instance;
    FILE *tsp_file = fopen (tsp_filename, "r");
    if (tsp_file == NULL) {
	perror (tsp_filename);
	exit (1);
    }
    tsp_instance = malloc (sizeof(t_tsp_instance));
    read_tsplib_instance (tsp_file, tsp_instance);
    return tsp_instance;
}
