#include "solution.h"
/* FIXME: This should not include solution.h but only btsp.h. 
 *#include "btsp.h"
 */
#include "libmisc.h"

/* FIXME: These two should be set in set_wls_params. */
extern int LocalSearch_type;
extern int LS_candlist_size;

void Sol_print_wls_params (FILE *stream)
{
    fprintf (stream, "LS type : %d", LocalSearch_type);
    if (LS_candlist_size < bTSP_size)
        fprintf (stream, " LS candidate list : %d", LS_candlist_size);
    else
        fprintf (stream, " LS candidate list : not used");
}

/*

       AAAA    CCCC   OOOO   TTTTTT   SSSSS  PPPPP
      AA  AA  CC     OO  OO    TT    SS      PP  PP
      AAAAAA  CC     OO  OO    TT     SSSS   PPPPP
      AA  AA  CC     OO  OO    TT        SS  PP
      AA  AA   CCCC   OOOO     TT    SSSSS   PP

######################################################
##########    ACO algorithms for the TSP    ##########
######################################################

      Version: 1.0
      File:    ls.c
      Author:  Thomas Stuetzle
      Purpose: implementation of local search routines
      Check:   README and gpl.txt
      Copyright (C) 1999  Thomas Stuetzle
*/

/***************************************************************************

    Program's name: acotsp

    Ant Colony Optimization algorithms (AS, ACS, EAS, RAS, MMAS, BWAS) for the 
    symmetric TSP 

    Copyright (C) 2004  Thomas Stuetzle

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    email: stuetzle no@spam informatik.tu-darmstadt.de
    mail address: Universitaet Darmstadt
                  Fachbereich Informatik
                  Hochschulstr. 10
                  D-64283 Darmstadt
		  Germany

***************************************************************************/

static int dlb_flag = TRUE;  /* flag indicating whether don't look bits are used. I recommend 
                                to always use it if local search is applied */

static int *
generate_random_permutation(int n)
/*    
      FUNCTION:       generate a random permutation of the integers 0 .. n-1
      INPUT:          length of the array
      OUTPUT:         pointer to the random permutation
      (SIDE)EFFECTS:  the array holding the random permutation is allocated in this 
                      function. Don't forget to free again the memory!
      COMMENTS:       only needed by the local search procedures
*/
{
   int  *r;

   r = create_int_vector (n);

   Rand_int_permutation (r, n);

   return r;
}

static void
two_opt_first(t_number **distance, int **nn_list, int nn_ls, int *tour)
/*    
      FUNCTION:       2-opt a tour 
      INPUT:          pointer to the tour that undergoes local optimization
      OUTPUT:         none
      (SIDE)EFFECTS:  tour is 2-opt
      COMMENTS:       the neighbourhood is scanned in random order (this need 
                      not be the best possible choice). Concerning the speed-ups used 
		      here consult, for example, Chapter 8 of
		      Holger H. Hoos and Thomas Stuetzle, 
		      Stochastic Local Search---Foundations and Applications, 
		      Morgan Kaufmann Publishers, 2004.
		      or some of the papers online available from David S. Johnson.
*/
{
    int c1, c2;             /* cities considered for an exchange */
    int s_c1, s_c2;         /* successor cities of c1 and c2     */
    int p_c1, p_c2;         /* predecessor cities of c1 and c2   */   
    int pos_c1, pos_c2;     /* positions of cities c1, c2        */
    int i, j, h, l;
    int improvement_flag, improve_node, help, n_improves = 0, n_exchanges=0;
    int h1=0, h2=0, h3=0, h4=0;
    t_number radius;             /* radius of nn-search */
    t_number gain = 0;
    int *random_vector;
    int *pos;               /* positions of cities in tour */ 
    int *dlb;               /* vector containing don't look bits */ 
    int n = bTSP_size;

    pos = malloc(n * sizeof(int));
    dlb = malloc(n * sizeof(int));
    for ( i = 0 ; i < n ; i++ ) {
	pos[tour[i]] = i;
	dlb[i] = FALSE;
    }

    improvement_flag = TRUE;
    random_vector = generate_random_permutation( n );

    while ( improvement_flag ) {

	improvement_flag = FALSE;

	for (l = 0 ; l < n; l++) {

	    c1 = random_vector[l]; 
	    assert ( c1 < n && c1 >= 0);
            if ( dlb_flag && dlb[c1] )
                continue;
	    improve_node = FALSE;
	    pos_c1 = pos[c1];
	    s_c1 = tour[pos_c1+1];
	    radius = distance[c1][s_c1];

	    /* First search for c1's nearest neighbours, use successor of c1 */
	    for ( h = 0 ; h < nn_ls ; h++ ) {
		c2 = nn_list[c1][h]; /* exchange partner, determine its position */
		if ( radius > distance[c1][c2] ) {
		    s_c2 = tour[pos[c2]+1];
		    gain =  - radius + distance[c1][c2] + 
			distance[s_c1][s_c2] - distance[c2][s_c2];
		    if ( gain < 0 ) {
			h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; 
			improve_node = TRUE;
			goto exchange2opt;
		    }
		}
		else 
		    break;
	    }      
	    /* Search one for next c1's h-nearest neighbours, use predecessor c1 */
	    if (pos_c1 > 0)
		p_c1 = tour[pos_c1-1];
	    else 
		p_c1 = tour[n-1];
	    radius = distance[p_c1][c1];
	    for ( h = 0 ; h < nn_ls ; h++ ) {
		c2 = nn_list[c1][h];  /* exchange partner, determine its position */
		if ( radius > distance[c1][c2] ) {
		    pos_c2 = pos[c2];
		    if (pos_c2 > 0)
			p_c2 = tour[pos_c2-1];
		    else 
			p_c2 = tour[n-1];
		    if ( p_c2 == c1 )
			continue;
		    if ( p_c1 == c2 )
			continue;
		    gain =  - radius + distance[c1][c2] + 
			distance[p_c1][p_c2] - distance[p_c2][c2];
		    if ( gain < 0 ) {
			h1 = p_c1; h2 = c1; h3 = p_c2; h4 = c2; 
			improve_node = TRUE;
			goto exchange2opt;
		    }
		}
		else 
		    break;
	    }      
	    if (improve_node) {
	    exchange2opt:
		n_exchanges++;
		improvement_flag = TRUE;
		dlb[h1] = FALSE; dlb[h2] = FALSE;
		dlb[h3] = FALSE; dlb[h4] = FALSE;
		/* Now perform move */
		if ( pos[h3] < pos[h1] ) {
		    help = h1; h1 = h3; h3 = help;
		    help = h2; h2 = h4; h4 = help;
		}
		if ( pos[h3] - pos[h2] < n / 2 + 1) {
		    /* reverse inner part from pos[h2] to pos[h3] */
		    i = pos[h2]; j = pos[h3];
		    while (i < j) {
			c1 = tour[i];
			c2 = tour[j];
			tour[i] = c2;
			tour[j] = c1;
			pos[c1] = j;
			pos[c2] = i;
			i++; j--;
		    }
		}
		else {
		    /* reverse outer part from pos[h4] to pos[h1] */
		    i = pos[h1]; j = pos[h4];
		    if ( j > i )
			help = n - (j - i) + 1;
		    else 
			help = (i - j) + 1;
		    help = help / 2;
		    for ( h = 0 ; h < help ; h++ ) {
			c1 = tour[i];
			c2 = tour[j];
			tour[i] = c2;
			tour[j] = c1;
			pos[c1] = j;
			pos[c2] = i;
			i--; j++;
			if ( i < 0 )
			    i = n-1;
			if ( j >= n )
			    j = 0;
		    }
		    tour[n] = tour[0];
		}
	    } else {
		dlb[c1] = TRUE;
	    }
      
	}
	if ( improvement_flag ) {
	    n_improves++;
	}
    }
    free( random_vector );
    free( dlb );
    free( pos );
}

static void
two_h_opt_first(t_number **distance, int **nn_list, int nn_ls, int *tour)
/*    
      FUNCTION:       2-h-opt a tour
      INPUT:          pointer to the tour that undergoes local optimization
      OUTPUT:         none
      (SIDE)EFFECTS:  tour is 2-h-opt
      COMMENTS:       for details on 2-h-opt see J. L. Bentley. Fast algorithms for geometric 
                      traveling salesman problems. ORSA Journal on Computing, 
		      4(4):387--411, 1992. 
		      The neighbourhood is scanned in random order (this need 
                      not be the best possible choice). Concerning the speed-ups used 
		      here consult, for example, Chapter 8 of
		      Holger H. Hoos and Thomas Stuetzle, 
		      Stochastic Local Search---Foundations and Applications, 
		      Morgan Kaufmann Publishers, 2004.
		      or some of the papers online available from David S. Johnson.
*/
{

    int c1, c2;         /* cities considered for an exchange */
    int s_c1, s_c2;     /* successors of c1 and c2           */
    int p_c1, p_c2;     /* predecessors of c1 and c2         */   
    int pos_c1, pos_c2;     /* positions of cities c1, c2        */
    int i, j, h, l;
    int improvement_flag, improve_node;
    int h1=0, h2=0, h3=0, h4=0, h5=0, help;
    t_number radius;             /* radius of nn-search */
    t_number gain = 0;
    int *random_vector;
    int two_move, node_move;

    int *pos;               /* positions of cities in tour */ 
    int *dlb;               /* vector containing don't look bits */ 
    int n = bTSP_size;
  
    pos = malloc(n * sizeof(int));
    dlb = malloc(n * sizeof(int));
    for ( i = 0 ; i < n ; i++ ) {
	pos[tour[i]] = i;
	dlb[i] = FALSE;
    }

    improvement_flag = TRUE;
    random_vector = generate_random_permutation( n );

    while ( improvement_flag ) {

	improvement_flag = FALSE; two_move = FALSE; node_move = FALSE;

	for (l = 0 ; l < n; l++) {

	    c1 = random_vector[l]; 
	    assert ( c1 < n && c1 >= 0);
		if ( dlb_flag && dlb[c1] )
		    continue;
	    improve_node = FALSE;
	    pos_c1 = pos[c1];
	    s_c1 = tour[pos_c1+1];
	    radius = distance[c1][s_c1];

	    /* First search for c1's nearest neighbours, use successor of c1 */
	    for ( h = 0 ; h < nn_ls ; h++ ) {
		c2 = nn_list[c1][h]; /* exchange partner, determine its position */
		if ( radius > distance[c1][c2] ) {
		    pos_c2 = pos[c2];
		    s_c2 = tour[pos_c2+1];
		    gain =  - radius + distance[c1][c2] + 
			distance[s_c1][s_c2] - distance[c2][s_c2];
		    if ( gain < 0 ) {
			h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; 
			improve_node = TRUE; two_move = TRUE; node_move = FALSE;
			goto exchange;
		    }
		    if (pos_c2 > 0)
			p_c2 = tour[pos_c2-1];
		    else 
			p_c2 = tour[n-1];
		    gain = - radius + distance[c1][c2] + distance[c2][s_c1]
			+ distance[p_c2][s_c2] - distance[c2][s_c2]
			- distance[p_c2][c2]; 
		    if ( c2 == s_c1 )
			gain = 0;
		    if ( p_c2 == s_c1 )
			gain = 0;

		    gain = 0;

		    if ( gain < 0 ) {
			h1 = c1; h2 = s_c1; h3 = c2; h4 = p_c2; h5 = s_c2; 
			improve_node = TRUE; node_move = TRUE; two_move = FALSE;
			goto exchange;
		    }
		}
		else 
		    break;
	    }      
	    /* Second search for c1's nearest neighbours, use predecessor c1 */
	    if (pos_c1 > 0)
		p_c1 = tour[pos_c1-1];
	    else 
		p_c1 = tour[n-1];
	    radius = distance[p_c1][c1];
	    for ( h = 0 ; h < nn_ls ; h++ ) {
		c2 = nn_list[c1][h];  /* exchange partner, determine its position */
		if ( radius > distance[c1][c2] ) {
		    pos_c2 = pos[c2];
		    if (pos_c2 > 0)
			p_c2 = tour[pos_c2-1];
		    else 
			p_c2 = tour[n-1];
		    if ( p_c2 == c1 )
			continue;
		    if ( p_c1 == c2 )
			continue;
		    gain =  - radius + distance[c1][c2] + 
			distance[p_c1][p_c2] - distance[p_c2][c2];
		    if ( gain < 0 ) {
			h1 = p_c1; h2 = c1; h3 = p_c2; h4 = c2; 
			improve_node = TRUE; two_move = TRUE; node_move = FALSE;
			goto exchange;
		    }
		    s_c2 = tour[pos[c2]+1];
		    gain = - radius + distance[c2][c1] + distance[p_c1][c2]
			+ distance[p_c2][s_c2] - distance[c2][s_c2]
			- distance[p_c2][c2]; 
		    if ( p_c1 == c2 )
			gain = 0;
		    if ( p_c1 == s_c2 )
			gain = 0;

		    if ( gain < 0 ) {
			h1 = p_c1; h2 = c1; h3 = c2; h4 = p_c2; h5 = s_c2; 
			improve_node = TRUE; node_move = TRUE; two_move = FALSE;
			goto exchange;
		    }
		}
		else 
		    break;
	    }      
	exchange:
	    if (improve_node) {
		if ( two_move ) {
		    improvement_flag = TRUE;
		    dlb[h1] = FALSE; dlb[h2] = FALSE;
		    dlb[h3] = FALSE; dlb[h4] = FALSE;
		    /* Now perform move */
		    if ( pos[h3] < pos[h1] ) {
			help = h1; h1 = h3; h3 = help;
			help = h2; h2 = h4; h4 = help;
		    }
		    if ( pos[h3] - pos[h2] < n / 2 + 1) {
			/* reverse inner part from pos[h2] to pos[h3] */
			i = pos[h2]; j = pos[h3];
			while (i < j) {
			    c1 = tour[i];
			    c2 = tour[j];
			    tour[i] = c2;
			    tour[j] = c1;
			    pos[c1] = j;
			    pos[c2] = i;
			    i++; j--;
			}
		    }
		    else {
			/* reverse outer part from pos[h4] to pos[h1] */
			i = pos[h1]; j = pos[h4];
			if ( j > i )
			    help = n - (j - i) + 1;
			else 
			    help = (i - j) + 1;
			help = help / 2;
			for ( h = 0 ; h < help ; h++ ) {
			    c1 = tour[i];
			    c2 = tour[j];
			    tour[i] = c2;
			    tour[j] = c1;
			    pos[c1] = j;
			    pos[c2] = i;
			    i--; j++;
			    if ( i < 0 )
				i = n-1;
			    if ( j >= n )
				j = 0;
			}
			tour[n] = tour[0];
		    }
		} else if ( node_move ) {
		    improvement_flag = TRUE;
		    dlb[h1] = FALSE; dlb[h2] = FALSE; dlb[h3] = FALSE; 
		    dlb[h4] = FALSE; dlb[h5] = FALSE;
		    /* Now perform move */
		    if ( pos[h3] < pos[h1] ) {
			help = pos[h1] - pos[h3];
			i = pos[h3];
			for ( h = 0 ; h < help ; h++ ) {
			    c1 = tour[i+1];
			    tour[i] = c1;
			    pos[c1] = i;
			    i++;
			}
			tour[i] = h3;
			pos[h3] = i;
			tour[n] = tour[0];
		    } else {
			/* pos[h3] > pos[h1] */ 
			help = pos[h3] - pos[h1];
			/*  	    if ( help < n / 2 + 1) { */
			i = pos[h3];
			for ( h = 0 ; h < help - 1 ; h++ ) {
			    c1 = tour[i-1];
			    tour[i] = c1;
			    pos[c1] = i;
			    i--;
			}
			tour[i] = h3;
			pos[h3] = i;
			tour[n] = tour[0];
			/*  	  } */
		    }  
		} else {
		    fprintf (stderr,"error: %s: unreacheable condition reached!\n",
                             __FUNCTION__);
		    exit(1);
		}
		two_move = FALSE; node_move = FALSE;
	    } else {
		dlb[c1] = TRUE;
	    }
      
	}
    }
    free( random_vector );
    free( dlb );
    free( pos );
}

static void
three_opt_first(t_number **distance, int **nn_list, int nn_ls, int *tour)

/*    
      FUNCTION:       3-opt the tour
      INPUT:          pointer to the tour that is to optimize
      OUTPUT:         none
      (SIDE)EFFECTS:  tour is 3-opt
      COMMENT:        this is certainly not the best possible implementation of a 3-opt 
                      local search algorithm. In addition, it is very lengthy; the main 
		      reason herefore is that awkward way of making an exchange, where 
		      it is tried to copy only the shortest possible part of a tour. 
		      Whoever improves the code regarding speed or solution quality, please 
		      drop me the code at stuetzle no@spam informatik.tu-darmstadt.de
		      The neighbourhood is scanned in random order (this need 
                      not be the best possible choice). Concerning the speed-ups used 
		      here consult, for example, Chapter 8 of
		      Holger H. Hoos and Thomas Stuetzle, 
		      Stochastic Local Search---Foundations and Applications, 
		      Morgan Kaufmann Publishers, 2004.
		      or some of the papers online available from David S. Johnson.
*/
{
    /* In case a 2-opt move should be performed, we only need to store opt2_move = TRUE,
       as h1, .. h4 are used in such a way that they store the indices of the correct move */

    int   c1, c2, c3;           /* cities considered for an exchange */
    int   s_c1, s_c2, s_c3;     /* successors of these cities        */
    int   p_c1, p_c2, p_c3;     /* predecessors of these cities      */   
    int   pos_c1, pos_c2, pos_c3;     /* positions of cities c1, c2, c3    */
    int   i, j, h, g, l;
    int   improvement_flag, help;
    int   h1=0, h2=0, h3=0, h4=0, h5=0, h6=0; /* memorize cities involved in a move */
    t_number   diffs, diffp;
    int   between = FALSE; 
    int   opt2_flag;  /* = TRUE: perform 2-opt move, otherwise none or 3-opt move */
    int   move_flag;  /* 
			      move_flag = 0 --> no 3-opt move 
			      move_flag = 1 --> between_move (c3 between c1 and c2)
			      move_flag = 2 --> not_between with successors of c2 and c3
			      move_flag = 3 --> not_between with predecessors of c2 and c3
			      move_flag = 4 --> cyclic move 
			   */
    t_number gain, move_value, radius, add1, add2;
    t_number decrease_breaks;    /* Stores decrease by breaking two edges (a,b) (c,d) */
    int val[3];
    int n1, n2, n3;
    int *pos;               /* positions of cities in tour */ 
    int *dlb;               /* vector containing don't look bits */ 
    int *h_tour;            /* help vector for performing exchange move */ 
    int *hh_tour;           /* help vector for performing exchange move */ 
    int *random_vector;
    int n = bTSP_size;

    pos = malloc(n * sizeof(int));
    dlb = malloc(n * sizeof(int));
    h_tour = malloc(n * sizeof(int));
    hh_tour = malloc(n * sizeof(int));

    for ( i = 0 ; i < n ; i++ ) {
	pos[tour[i]] = i;
	dlb[i] = FALSE;
    }
    improvement_flag = TRUE;
    random_vector = generate_random_permutation( n );

    while ( improvement_flag ) {
	move_value = 0;
	improvement_flag = FALSE;

	for ( l = 0 ; l < n ; l++ ) {

	    c1 = random_vector[l];
	    if ( dlb_flag && dlb[c1] )
		continue;
	    opt2_flag = FALSE;

	    move_flag = 0;
	    pos_c1 = pos[c1];
	    s_c1 = tour[pos_c1+1];
	    if (pos_c1 > 0)
		p_c1 = tour[pos_c1-1];
	    else 
		p_c1 = tour[n-1];

	    h = 0;    /* Search for one of the h-nearest neighbours */
	    while ( h < nn_ls ) {

		c2   = nn_list[c1][h];  /* second city, determine its position */
		pos_c2 = pos[c2];
		s_c2 = tour[pos_c2+1];
		if (pos_c2 > 0)
		    p_c2 = tour[pos_c2-1];
		else 
		    p_c2 = tour[n-1];
	  
		diffs = 0; diffp = 0;

		radius = distance[c1][s_c1];
		add1   = distance[c1][c2];

		/* Here a fixed radius neighbour search is performed */
		if ( radius > add1 ) {
		    decrease_breaks = - radius - distance[c2][s_c2];
		    diffs =  decrease_breaks + add1 + distance[s_c1][s_c2];
		    diffp =  - radius - distance[c2][p_c2] + 
			distance[c1][p_c2] + distance[s_c1][c2];
		}
		else 
		    break;
		if ( p_c2 == c1 )  /* in case p_c2 == c1 no exchange is possible */
		    diffp = 0;
		if ( (diffs < move_value) || (diffp < move_value) ) {
		    improvement_flag = TRUE; 
		    if (diffs <= diffp) { 
			h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; 
			move_value = diffs; 
			opt2_flag = TRUE; move_flag = 0;
			/*     	    goto exchange; */
		    } else {
			h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; 
			move_value = diffp;  
			opt2_flag = TRUE; move_flag = 0;
			/*     	    goto exchange; */
		    }
		}
		/* Now perform the innermost search */
		g = 0;
		while (g < nn_ls) {
	  
		    c3   = nn_list[s_c1][g];
		    pos_c3 = pos[c3];
		    s_c3 = tour[pos_c3+1];
		    if (pos_c3 > 0)
			p_c3 = tour[pos_c3-1];
		    else 
			p_c3 = tour[n-1];
		  
		    if ( c3 == c1 ) {
			g++;
			continue;
		    }
		    else {
			add2 = distance[s_c1][c3];
			/* Perform fixed radius neighbour search for innermost search */
			if ( decrease_breaks + add1 < add2 ) {
			  
			    if ( pos_c2 > pos_c1 ) {
				if ( pos_c3 <= pos_c2 && pos_c3 > pos_c1 )
				    between = TRUE;
				else 
				    between = FALSE;
			    }
			    else if ( pos_c2 < pos_c1 )
				if ( pos_c3 > pos_c1 || pos_c3 < pos_c2 )
				    between = TRUE;
				else 
				    between = FALSE;
			    else {
				printf(" Strange !!, pos_1 %d == pos_2 %d, \n", pos_c1, pos_c2);
			    }
			  
			    if ( between ) {
				/* We have to add edges (c1,c2), (c3,s_c1), (p_c3,s_c2) to get 
				   valid tour; it's the only possibility */
			      
				gain = decrease_breaks - distance[c3][p_c3] +
				    add1 + add2 +
				    distance[p_c3][s_c2];
			      
				/* check for improvement by move */
				if ( gain < move_value ) {
				    improvement_flag = TRUE; /* g = neigh_ls + 1; */
				    move_value = gain;
				    opt2_flag = FALSE;
				    move_flag = 1;
				    /* store nodes involved in move */
				    h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; h5 = p_c3; h6 = c3;
				    goto exchange;
				} 
			    }
			    else {   /* not between(pos_c1,pos_c2,pos_c3) */
			      
				/* We have to add edges (c1,c2), (s_c1,c3), (s_c2,s_c3) */
			      
				gain = decrease_breaks - distance[c3][s_c3] +
				    add1 + add2 + 
				    distance[s_c2][s_c3];
			      
				if ( pos_c2 == pos_c3 ) {
				    gain = 20000;
				}
			      
				/* check for improvement by move */
				if ( gain < move_value ) {
				    improvement_flag = TRUE; /* g = neigh_ls + 1; */
				    move_value = gain;
				    opt2_flag = FALSE;
				    move_flag = 2;
				    /* store nodes involved in move */
				    h1 = c1; h2 = s_c1; h3 = c2; h4 = s_c2; h5 = c3; h6 = s_c3;
				    goto exchange;
				}
			      
				/* or add edges (c1,c2), (s_c1,c3), (p_c2,p_c3) */
				gain = - radius - distance[p_c2][c2] 
				    - distance[p_c3][c3] +
				    add1 + add2 + 
				    distance[p_c2][p_c3];
			      
				if ( c3 == c2 || c2 == c1 || c1 == c3 || p_c2 == c1 ) {
				    gain = 2000000;
				}
			      
				if ( gain < move_value ) {
				    improvement_flag = TRUE;
				    move_value = gain;
				    opt2_flag = FALSE;
				    move_flag = 3;
				    h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; h5 = p_c3; h6 = c3;
				    goto exchange;
				}
			      
				/* Or perform the 3-opt move where no subtour inversion is necessary 
				   i.e. delete edges (c1,s_c1), (c2,p_c2), (c3,s_c3) and 
				   add edges (c1,c2), (c3,s_c1), (p_c2,s_c3) */
			      
				gain = - radius - distance[p_c2][c2] - 
				    distance[c3][s_c3]
				    + add1 + add2 + distance[p_c2][s_c3];
			      
				/* check for improvement */
				if ( gain < move_value ) {
				    improvement_flag = TRUE;
				    move_value = gain;
				    opt2_flag = FALSE;
				    move_flag = 4;
				    improvement_flag = TRUE;
				    /* store nodes involved in move */
				    h1 = c1; h2 = s_c1; h3 = p_c2; h4 = c2; h5 = c3; h6 = s_c3; 
				    goto exchange;
				}
			    }
			}
			else
			    g = nn_ls + 1;
		    }
		    g++;
		}
		h++;
	    }
	    if ( move_flag || opt2_flag ) {
	    exchange:
		move_value = 0;

		/* Now make the exchange */
		if ( move_flag ) {
		    dlb[h1] = FALSE; dlb[h2] = FALSE; dlb[h3] = FALSE; 
		    dlb[h4] = FALSE; dlb[h5] = FALSE; dlb[h6] = FALSE;
		    pos_c1 = pos[h1]; pos_c2 = pos[h3]; pos_c3 = pos[h5];
		  
		    if ( move_flag == 4 ) {

			if ( pos_c2 > pos_c1 ) 
			    n1 = pos_c2 - pos_c1;
			else
			    n1 = n - (pos_c1 - pos_c2);
			if ( pos_c3 > pos_c2 ) 
			    n2 = pos_c3 - pos_c2;
			else
			    n2 = n - (pos_c2 - pos_c3);
			if ( pos_c1 > pos_c3 ) 
			    n3 = pos_c1 - pos_c3;
			else
			    n3 = n - (pos_c3 - pos_c1);
		      
			/* n1: length h2 - h3, n2: length h4 - h5, n3: length h6 - h1 */
			val[0] = n1; val[1] = n2; val[2] = n3; 
			/* Now order the partial tours */
			h = 0;
			help = INT_MIN;
			for ( g = 0; g <= 2; g++) {
			    if ( help < val[g] ) {
				help = val[g];
				h = g;
			    }
			}
		      
			/* order partial tours according length */
			if ( h == 0 ) {
			    /* copy part from pos[h4] to pos[h5]
			       direkt kopiert: Teil von pos[h6] to pos[h1], it
			       remains the part from pos[h2] to pos[h3] */
			    j = pos[h4];
			    h = pos[h5];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h) {
				i++;
				j++;
				if ( j  >= n )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    /* First copy partial tour 3 in new position */
			    j = pos[h4];
			    i = pos[h6];
			    tour[j] = tour[i];
			    pos[tour[i]] = j; 
			    while ( i != pos_c1) {
				i++;
				if ( i >= n )
				    i = 0;
				j++;
				if ( j >= n )
				    j = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j; 
			    }
			  
			    /* Now copy stored part from h_tour */
			    j++;
			    if ( j >= n )
				j = 0;
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
			else if ( h == 1 ) {
			  
			    /* copy part from pos[h6] to pos[h1]
			       direkt kopiert: Teil von pos[h2] to pos[h3], it
			       remains the part from pos[h4] to pos[h5] */
			    j = pos[h6];
			    h = pos[h1];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h) {
				i++;
				j++;
				if ( j  >= n )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    /* First copy partial tour 3 in new position */
			    j = pos[h6];
			    i = pos[h2];
			    tour[j] = tour[i];
			    pos[tour[i]] = j; 
			    while ( i != pos_c2) {
				i++;
				if ( i >= n )
				    i = 0;
				j++;
				if ( j >= n )
				    j = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j; 
			    }
			  
			    /* Now copy stored part from h_tour */
			    j++;
			    if ( j >= n )
				j = 0;
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
			else if ( h == 2 ) {
			    /* copy part from pos[h2] to pos[h3]
			       direkt kopiert: Teil von pos[h4] to pos[h5], it
			       remains the part from pos[h6] to pos[h1] */
			    j = pos[h2];
			    h = pos[h3];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h) {
				i++;
				j++;
				if ( j  >= n )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
	      
			    /* First copy partial tour 3 in new position */
			    j = pos[h2];
			    i = pos[h4];
			    tour[j] = tour[i];
			    pos[tour[i]] = j; 
			    while ( i != pos_c3) {
				i++;
				if ( i >= n )
				    i = 0;
				j++;
				if ( j >= n )
				    j = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j; 
			    }
			  
			    /* Now copy stored part from h_tour */
			    j++;
			    if ( j >= n )
				j = 0;
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i]; 
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		    }
		    else if ( move_flag == 1 ) {
		      
			if ( pos_c3 < pos_c2 ) 
			    n1 = pos_c2 - pos_c3;
			else
			    n1 = n - (pos_c3 - pos_c2);
			if ( pos_c3 > pos_c1 ) 
			    n2 = pos_c3 - pos_c1 + 1;
			else
			    n2 = n - (pos_c1 - pos_c3 + 1);
			if ( pos_c2 > pos_c1 ) 
			    n3 = n - (pos_c2 - pos_c1 + 1);
			else
			    n3 = pos_c1 - pos_c2 + 1;
		      
			/* n1: length h6 - h3, n2: length h5 - h2, n2: length h1 - h3 */
			val[0] = n1; val[1] = n2; val[2] = n3; 
			/* Now order the partial tours */
			h = 0;
			help = INT_MIN;
			for ( g = 0; g <= 2; g++) {
			    if ( help < val[g] ) {
				help = val[g];
				h = g;
			    }
			}
			/* order partial tours according length */
		      
			if ( h == 0 ) {
			  
			    /* copy part from pos[h5] to pos[h2]
			       (inverted) and from pos[h4] to pos[h1] (inverted)
			       it remains the part from pos[h6] to pos[h3] */
			    j = pos[h5];
			    h = pos[h2];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h1];
			    h = pos[h4];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }
			  
			    j = pos[h4];
			    for ( i = 0; i< n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if (j >= n)
				    j = 0;
			    }
			  
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i< n1 ; i++ ) {
				tour[j] = h_tour[i]; 
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
			else if ( h == 1 ) {
			  
			    /* copy part from h3 to h6 (wird inverted) erstellen : */
			    j = pos[h3];
			    h = pos[h6];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h) {
				i++;
				j--;
				if ( j  < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h6];
			    i = pos[h4];
			  
			    tour[j] = tour[i];
			    pos[tour[i]] = j; 
			    while ( i != pos_c1) {
				i++;
				j++;
				if ( j >= n)
				    j = 0;
				if ( i >= n)
				    i = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j; 
			    }
			  
			    /* Now copy stored part from h_tour */
			    j++;
			    if ( j >= n )
				j = 0;
			    i = 0;
			    tour[j] = h_tour[i];
			    pos[h_tour[i]] = j; 
			    while ( j != pos_c1 ) {
				j++;
				if ( j >= n )
				    j = 0;
				i++;
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
			    }
			    tour[n] = tour[0];
			}
		      
			else if ( h == 2 ) {
			  
			    /* copy part from pos[h2] to pos[h5] and
			       from pos[h3] to pos[h6] (inverted), it
			       remains the part from pos[h4] to pos[h1] */
			    j = pos[h2];
			    h = pos[h5];
			    i = 0;
			    h_tour[i] =  tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j++;
				if ( j >= n )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
			    j = pos_c2;
			    h = pos[h6];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }
			  
			    j = pos[h2];
			    for ( i = 0; i< n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n)
				    j = 0;
			    }
			  
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i< n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		    }
		    else if ( move_flag == 2 ) {
		      
			if ( pos_c3 < pos_c1 ) 
			    n1 = pos_c1 - pos_c3;
			else
			    n1 = n - (pos_c3 - pos_c1);
			if ( pos_c3 > pos_c2 ) 
			    n2 = pos_c3 - pos_c2;
			else
			    n2 = n - (pos_c2 - pos_c3);
			if ( pos_c2 > pos_c1 ) 
			    n3 = pos_c2 - pos_c1;
			else
			    n3 = n - (pos_c1 - pos_c2);
		      
			val[0] = n1; val[1] = n2; val[2] = n3; 
			/* Determine which is the longest part */
			h = 0;
			help = INT_MIN;
			for ( g = 0; g <= 2; g++) {
			    if ( help < val[g] ) {
				help = val[g];
				h = g;
			    }
			}
			/* order partial tours according length */
		      
			if ( h == 0 ) {
			  
			    /* copy part from pos[h3] to pos[h2]
			       (inverted) and from pos[h5] to pos[h4], it
			       remains the part from pos[h6] to pos[h1] */
			    j = pos[h3];
			    h = pos[h2];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h5];
			    h = pos[h4];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }
			  
			    j = pos[h2];
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i]; 
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
	      
			    for ( i = 0; i < n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			    /*  	      getchar(); */
			}
			else if ( h == 1 ) {
			  
			    /* copy part from pos[h2] to pos[h3] and
			       from pos[h1] to pos[h6] (inverted), it
			       remains the part from pos[h4] to pos[h5] */
			    j = pos[h2];
			    h = pos[h3];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j++;
				if ( j >= n  )
				    j = 0;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h1];
			    h = pos[h6];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j =  n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }
			    j = pos[h6];
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i]; 
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    for ( i = 0; i < n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		      
			else if ( h == 2 ) {
			  
			    /* copy part from pos[h1] to pos[h6]
			       (inverted) and from pos[h4] to pos[h5],
			       it remains the part from pos[h2] to
			       pos[h3] */
			    j = pos[h1];
			    h = pos[h6];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }

			    j = pos[h4];
			    h = pos[h5];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j++;
				if ( j >= n  )
				    j = 0;
				hh_tour[i] = tour[j];
				n2++;
			    }

			    j = pos[h4];
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			  
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i < n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		    }
		    else if ( move_flag == 3 ) {
		      
			if ( pos_c3 < pos_c1 ) 
			    n1 = pos_c1 - pos_c3;
			else
			    n1 = n - (pos_c3 - pos_c1);
			if ( pos_c3 > pos_c2 ) 
			    n2 = pos_c3 - pos_c2;
			else
			    n2 = n - (pos_c2 - pos_c3);
			if ( pos_c2 > pos_c1 ) 
			    n3 = pos_c2 - pos_c1;
			else
			    n3 = n - (pos_c1 - pos_c2);
			/* n1: length h6 - h1, n2: length h4 - h5, n2: length h2 - h3 */
		      
			val[0] = n1; val[1] = n2; val[2] = n3; 
			/* Determine which is the longest part */
			h = 0;
			help = INT_MIN;
			for ( g = 0; g <= 2; g++) {
			    if ( help < val[g] ) {
				help = val[g];
				h = g;
			    }
			}
			/* order partial tours according length */
		      
			if ( h == 0 ) {
			  
			    /* copy part from pos[h2] to pos[h3]
			       (inverted) and from pos[h4] to pos[h5]
			       it remains the part from pos[h6] to pos[h1] */
			    j = pos[h3];
			    h = pos[h2];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }
			  
			    j = pos[h2];
			    h = pos[h5];
			    i = pos[h4];
			    tour[j] = h4;
			    pos[h4] = j;
			    while ( i != h ) {
				i++;
				if ( i >= n )
				    i = 0;
				j++;
				if ( j >= n )
				    j = 0;
				tour[j] = tour[i];
				pos[tour[i]] = j;
			    }
			    j++;
			    if ( j >= n )
				j = 0;
			    for ( i = 0; i < n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
			else if ( h == 1 ) {

			    /* copy part from pos[h3] to pos[h2]
			       (inverted) and from  pos[h6] to pos[h1],
			       it remains the part from pos[h4] to pos[h5] */
			    j = pos[h3];
			    h = pos[h2];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0  )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }

			    j = pos[h6];
			    h = pos[h1];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j++;
				if ( j >= n )
				    j = 0;
				hh_tour[i] = tour[j];
				n2++;
			    }
			  
			    j = pos[h6];
			    for ( i = 0; i<n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }

			    for ( i = 0 ; i < n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		      
			else if ( h == 2 ) {
			  
			    /* copy part from pos[h4] to pos[h5]
			       (inverted) and from pos[h6] to pos[h1] (inverted)
			       it remains the part from pos[h2] to pos[h3] */
			    j = pos[h5];
			    h = pos[h4];
			    i = 0;
			    h_tour[i] = tour[j];
			    n1 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				h_tour[i] = tour[j];
				n1++;
			    }

			    j = pos[h1];
			    h = pos[h6];
			    i = 0;
			    hh_tour[i] = tour[j];
			    n2 = 1;
			    while ( j != h ) {
				i++;
				j--;
				if ( j < 0 )
				    j = n-1;
				hh_tour[i] = tour[j];
				n2++;
			    }

			    j = pos[h4];
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i< n1 ; i++ ) {
				tour[j] = h_tour[i];
				pos[h_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    /* Now copy stored part from h_tour */
			    for ( i = 0; i< n2 ; i++ ) {
				tour[j] = hh_tour[i];
				pos[hh_tour[i]] = j; 
				j++;
				if ( j >= n )
				    j = 0;
			    }
			    tour[n] = tour[0];
			}
		    }
		    else {
			printf(" Some very strange error must have occurred !!!\n\n");
			exit(0);
		    }
		}
		if (opt2_flag) {

		    /* Now perform move */
		    dlb[h1] = FALSE; dlb[h2] = FALSE;
		    dlb[h3] = FALSE; dlb[h4] = FALSE;
		    if ( pos[h3] < pos[h1] ) {
			help = h1; h1 = h3; h3 = help;
			help = h2; h2 = h4; h4 = help;
		    }
		    if ( pos[h3]-pos[h2] < n / 2 + 1) {
			/* reverse inner part from pos[h2] to pos[h3] */
			i = pos[h2]; j = pos[h3];
			while (i < j) {
			    c1 = tour[i];
			    c2 = tour[j];
			    tour[i] = c2;
			    tour[j] = c1;
			    pos[c1] = j;
			    pos[c2] = i;
			    i++; j--;
			}
		    }
		    else {
			/* reverse outer part from pos[h4] to pos[h1] */
			i = pos[h1]; j = pos[h4];
			if ( j > i )
			    help = n - (j - i) + 1;
			else 
			    help = (i - j) + 1;
			help = help / 2;
			for ( h = 0 ; h < help ; h++ ) {
			    c1 = tour[i];
			    c2 = tour[j];
			    tour[i] = c2;
			    tour[j] = c1;
			    pos[c1] = j;
			    pos[c2] = i;
			    i--; j++;
			    if ( i < 0 )
				i = n - 1;
			    if ( j >= n )
				j = 0;
			}
			tour[n] = tour[0];
		    }
		}
	    }
	    else {
		dlb[c1] = TRUE;
	    }
	}
    }
    free( random_vector );
    free( h_tour );
    free( hh_tour );
    free( pos );
    free( dlb );
}


t_number
SolWeighted_local_search(t_solution *solution, t_number lambda,
                         t_number max_weight)
/*    
      FUNCTION:       manage the local search phase; apply local search to ALL ants; in 
                      dependence of ls_flag one of 2-opt, 2.5-opt, and 3-opt local search
		      is chosen.
      INPUT:          none
      OUTPUT:         none
      (SIDE)EFFECTS:  all ants of the colony have locally optimal tours
      COMMENTS:       typically, best performance is obtained by applying local search 
                      to all ants. It is known that some improvements (e.g. convergence 
		      speed towards high quality solutions) may be obtained for some 
		      ACO algorithms by applying local search to only some of the ants.
		      Overall best performance is typcially obtained by using 3-opt.
*/
{
    int n = btsp_instance.size;
    t_number previous_value, new_value;
    static t_number **wdist = NULL;
    static t_number prev_lambda;
    assert (bTSP_wcandlist);
    if (wdist == NULL) {
        prev_lambda = lambda + 1;
        wdist = create_number_matrix (n, n);
    }

    previous_value = (max_weight - lambda) * SolGetObjective (solution, 1) 
        + lambda * SolGetObjective (solution, 2);

    
    if (prev_lambda != lambda) {
        calc_weighted_matrix (wdist,
                              btsp_instance.distance1,
                              btsp_instance.distance2,
                              btsp_instance.size,
                              btsp_instance.weights[lambda],
                              btsp_instance.num_weights - 1);
    }

    switch (LocalSearch_type) {
    case 1: /* 2-opt local search */
        two_opt_first (wdist,
                       bTSP_wcandlist[lambda],
                       LS_candlist_size,
                       SolGetVector(solution));    
        break;
    case 2: /* 2.5-opt local search */
        two_h_opt_first (wdist,
                         bTSP_wcandlist[lambda],
                         LS_candlist_size,
                         SolGetVector(solution)); 
        break;
    case 3: /* 3-opt local search */
        three_opt_first (wdist,
                         bTSP_wcandlist[lambda],
                         LS_candlist_size,
                         SolGetVector(solution));  
        break;
    default:
        eprintf ("%s(): "
                 "type of local search procedure (%d) not implemented!\n",
                 __FUNCTION__, LocalSearch_type);
    }

    SolEvaluate (solution);
    DEBUG1 (SolCheck (solution));

    new_value = (max_weight - lambda) * SolGetObjective (solution, 1) 
        + lambda * SolGetObjective (solution, 2);

    DEBUG2_FUNPRINT ("previous_value = %g  new_value = %g"
                     "  improvement = %g\n", 
                     (double) previous_value, (double) new_value,
                     (double) (new_value - previous_value));

    prev_lambda = lambda;

    return new_value - previous_value;
}

