//most updated pthread program
#include <stdio.h>
#include <pthread.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define         X_RESN  200       /* x resolution */
#define         Y_RESN  200       /* y resolution */
#define         no_threads 3
#define         t_interval 0.01
#define         num_n_body 500
#define         g_const 1.0
#define			total_times 1000

void   initRandomSeed(); 
double randomReal(double low, double high);
double *randomRealSets (int n, double low, double high);
void   printInfo(double runningTime);

pthread_mutex_t mutex1;
double up_boundary = Y_RESN*(2.0/3);
double bottom_boundary = Y_RESN*(1.0/3);
double left_boundary = X_RESN*(1.0/3);
double right_boundary = X_RESN*(2.0/3);
double *x_coordinate = randomRealSets(num_n_body, left_boundary, right_boundary);
double *y_coordinate = randomRealSets(num_n_body, bottom_boundary, up_boundary);
double *x_volocity   = randomRealSets(num_n_body, 0, 0);
double *y_volocity   = randomRealSets(num_n_body, 0, 0);
double *n_mass       = randomRealSets(num_n_body, 0.5, 1.5);

int    collision [num_n_body];
double x_factor  [num_n_body];
double y_factor  [num_n_body];

int    average_tasks =  int( floor(num_n_body/(no_threads)) );
int    global_index  = 0;

void *slave(void *ignored){

	int local_index, end, has_collided;
    double mass, mass2, Fx, Fy, fx, fy, r,
           x_pos_1, y_pos_1, x_vol_1, y_vol_1,
           x_pos_2, y_pos_2;

	pthread_mutex_lock(&mutex1);
	local_index = global_index;
	global_index += average_tasks;
	pthread_mutex_unlock(&mutex1);

	end = (local_index+average_tasks < num_n_body) ? 
		  (local_index+average_tasks) : (num_n_body-1);

    for(int i = local_index; i < end; i++) {

        mass = n_mass[i];
        Fx = Fy = 0;   
        has_collided = 0;
        collision[i] = -1;
        x_pos_1 = x_coordinate[i];
        y_pos_1 = y_coordinate[i];
        x_vol_1 = x_volocity[i];
        y_vol_1 = y_volocity[i];

        if (x_pos_1 < left_boundary) {
            x_vol_1 = -x_vol_1;
            x_pos_1 = left_boundary + (left_boundary - x_pos_1);
        } 
        else if (x_pos_1 > right_boundary) {
            x_vol_1 = -x_vol_1;
            x_pos_1 = right_boundary - (x_pos_1 - right_boundary);                  
        } 
        if (y_pos_1 < bottom_boundary) {
            y_vol_1 = -y_vol_1;
            y_pos_1 = bottom_boundary + (bottom_boundary - y_pos_1);  
        }
        else if (y_pos_1 > up_boundary) {
            y_vol_1 = -y_vol_1;
            y_pos_1 = up_boundary - (y_pos_1 - up_boundary); 
        }
            
        for(int k = 0; k < num_n_body; k++) { 
            if (k == i) continue;

            // pthread_mutex_lock(&mutex1);
            x_pos_2 = x_coordinate[k];
            y_pos_2 = y_coordinate[k];                   
            // pthread_mutex_unlock(&mutex1);

            r = sqrt((x_pos_1 - x_pos_2)*
                    ( x_pos_1 - x_pos_2)+
                    ( y_pos_1 - y_pos_2)*
                    ( y_pos_1 - y_pos_2));
            
            if (r < 2) {
                if (i < k && has_collided == 0) {
                    collision[i] = k;
                    has_collided = 1;
                }
                r = 2;
            }        
            mass2 = n_mass[k];
            fx = mass * mass2 * (x_pos_2- x_pos_1)/(r*r*r);
            fy = mass * mass2 * (y_pos_2- y_pos_1)/(r*r*r);                                 
            Fx += fx;
            Fy += fy;
        }

        x_factor[i] = Fx/mass;
        y_factor[i] = Fy/mass;

        x_volocity[i]   =  x_vol_1;
        y_volocity[i]   =  y_vol_1;
        x_coordinate[i] =  x_pos_1;
        y_coordinate[i] =  y_pos_1;

    }
}

int main() {

	Window          win;                            /* initialization for a window */
    unsigned
    int             width, height,                  /* window size */
                    x, y,                           /* window position */
                    border_width,                   /*border width in pixels */
                    display_width, display_height,  /* size of screen */
                    screen;                         /* which screen */
    char            *window_name = "Mandelbrot Set", *display_name = NULL;
    GC              gc;
    unsigned
    long            valuemask = 0;
    XGCValues       values;
    Display         *display;
    XSizeHints      size_hints;
    Pixmap          bitmap;
    XPoint          points[800];
    FILE            *fp, *fopen ();
    char            str[100];
    XSetWindowAttributes attr[1];
    timeval t_start, t_end;
   
    /* connect to Xserver */
    if (  (display = XOpenDisplay (display_name)) == NULL ) {
        fprintf (stderr, "drawon: cannot connect to X server %s\n",
                    XDisplayName (display_name) );
    exit (-1);
    }
        
    /* get screen size */

    screen = DefaultScreen (display);
    display_width = DisplayWidth (display, screen);
    display_height = DisplayHeight (display, screen);

    /* set window size */

    width = X_RESN;
    height = Y_RESN;

    /* set window position */

    x = 0;
    y = 0;

    /* create opaque window */

    border_width = 4;
    win = XCreateSimpleWindow (display, RootWindow (display, screen),
                            x, y, width, height, border_width, 
                            BlackPixel (display, screen), WhitePixel (display, screen));

    size_hints.flags = USPosition|USSize;
    size_hints.x = x;
    size_hints.y = y;
    size_hints.width = width;
    size_hints.height = height;
    size_hints.min_width = 300;
    size_hints.min_height = 300;
    
    XSetNormalHints (display, win, &size_hints);
    XStoreName(display, win, window_name);

    /* create graphics context */

    gc = XCreateGC (display, win, valuemask, &values);

    XSetBackground (display, gc, WhitePixel (display, screen));
    XSetForeground (display, gc, BlackPixel (display, screen));
    XSetLineAttributes (display, gc, 1, LineSolid, CapRound, JoinRound);

    attr[0].backing_store = Always;
    attr[0].backing_planes = 1;
    attr[0].backing_pixel = BlackPixel(display, screen);

    XChangeWindowAttributes(display, win, CWBackingStore | CWBackingPlanes | CWBackingPixel, attr);

    XMapWindow (display, win);
    XSync(display, 0);

    /*create slaves */
	int i, k;
	pthread_t thread[no_threads];
	pthread_mutex_init(&mutex1, NULL);

	int times = 1;
    double mass, mass2, temp_x_vol, temp_y_vol, 
           x_vol_1, y_vol_1, x_vol_2, y_vol_2, 
           x_pos_1, y_pos_1;

	gettimeofday(&t_start, NULL);
	while(times < total_times) {
		XClearWindow(display, win);
		global_index = 0;
		for (i = 0; i < no_threads; i++) {
			if (pthread_create(&thread[i], NULL, slave, NULL) != 0)
				perror("Pthread_create fails");
		}
		for (i = 0; i < no_threads; i++) {
			if (pthread_join(thread[i], NULL) != 0)
				perror("Pthread_join fails");
		}

        for(int i = 0; i < num_n_body; i++) {
            printf("collision %d\n", collision[i]);
        }

		for(int i = 0; i < num_n_body; i++) {

            x_vol_1 = x_volocity  [i];
            y_vol_1 = y_volocity  [i];
            x_pos_1 = x_coordinate[i];
            y_pos_1 = y_coordinate[i];

            k = collision[i];
            
            if (k != -1) {

                mass  = n_mass[i];               
                mass2 = n_mass[k]; 
                x_vol_2 = x_volocity[k];
                y_vol_2 = y_volocity[k];              
                temp_x_vol = x_vol_1;
                temp_y_vol = y_vol_1;

                x_vol_1 = ( (mass-mass2)*x_vol_1+2*mass2*x_vol_2 ) /
                                  (mass+mass2);
                y_vol_1 = ( (mass-mass2)*y_vol_1+2*mass2*y_vol_2 ) /
                                  (mass+mass2);
                x_vol_2 = ( (mass2-mass)*x_vol_2+2*mass*temp_x_vol ) /
                                  (mass+mass2);
                y_vol_2 = ( (mass2-mass)*y_vol_2+2*mass*temp_y_vol ) /
                                 (mass+mass2);

                x_volocity[k] = x_vol_2;
                y_volocity[k] = y_vol_2;

                
            }

            x_vol_1 = x_vol_1 + x_factor[i]*t_interval;
            y_vol_1 = y_vol_1 + y_factor[i]*t_interval;
            x_pos_1 = x_pos_1 + x_vol_1*t_interval;
            y_pos_1 = y_pos_1 + y_vol_1*t_interval;

            x_coordinate[i] = x_pos_1;
            y_coordinate[i] = y_pos_1;
            x_volocity  [i] = x_vol_1;
            y_volocity  [i] = y_vol_1;

            XDrawPoint (display, win, gc, x_pos_1, y_pos_1);
            // usleep(1);
            printf("");
        }

		XFlush (display);
		times++;
	}
	gettimeofday(&t_end, NULL);
    double runningTime = (t_end.tv_sec - t_start.tv_sec) + 
                         (t_end.tv_usec - t_start.tv_usec) / 1000000.0;
    printInfo(runningTime);
    //printf("times: %d\n", times);
 }

void initRandomSeed() {
    static bool initialized = false;
    if (!initialized) {
        srand(int(time(NULL)));
        initialized = true;
    }
}

double randomReal(double low, double high) {
    // initRandomSeed();
    double d = rand() / (double(RAND_MAX)+1);
    double s = d*(high-low);
    return low + s;
}

double *randomRealSets (int n, double low, double high) {
    double *real_sets;
    real_sets = (double*)malloc((n) * sizeof(double));
    for (int i = 0; i < n; i++) {
        real_sets[i] = randomReal(low, high);
    }
    return real_sets;
}

void printInfo(double runningTime) {
    printf("\nName: %s\n", "Liang Jialu");
    printf("Student ID: %d\n", 118010164);
    printf("Assignment 3, N Body, Pthread implementation.\n");
    printf("runTime is %f\n", runningTime);
    printf("number of threads = %d", no_threads);
}