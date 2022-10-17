#include <omp.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <stdio.h>

#define         X_RESN  200       /* x resolution */
#define         Y_RESN  200       /* y resolution */
#define         t_interval 0.01
#define         num_n_body 500
#define         g_const 1.0
#define         total_times 1000
#define         n_thread 2


void   initRandomSeed(); 
double randomReal(double low, double high);
double *randomRealSets (int n, double low, double high);
void   printInfo(double runningTime);
int    collision[num_n_body];
double x_factor [num_n_body];
double y_factor [num_n_body];

int main (int argc, char **argv)
{
    Window          win;                            /* initialization for a window */
    unsigned
    int             width, height,                  /* window size */
                    x, y,                           /* window position */
                    border_width,                   /*border width in pixels */
                    display_width, display_height,  /* size of screen */
                    screen;                         /* which screen */

    char            *window_name = "openMP N Body", *display_name = NULL;
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

    timeval t_start, t_end;
    double up_boundary = Y_RESN*(2.0/3);
    double bottom_boundary = Y_RESN*(1.0/3);
    double left_boundary = X_RESN*(1.0/3);
    double right_boundary = X_RESN*(2.0/3);
    double *x_coordinate = randomRealSets(num_n_body, left_boundary, right_boundary);
    double *y_coordinate = randomRealSets(num_n_body, bottom_boundary, up_boundary);
    double *x_volocity   = randomRealSets(num_n_body, 0, 0);
    double *y_volocity   = randomRealSets(num_n_body, 0, 0);
    double *n_mass       = randomRealSets(num_n_body, 0.5, 1.5);
    double mass, Fx, Fy;

    int times = 1;  

    gettimeofday(&t_start, NULL);  
    while (times < total_times) {

        XClearWindow(display, win);

        int A;

        #pragma omp parallel for num_threads(n_thread) private(mass, Fx, Fy)
        for(int i = 0; i < num_n_body; i++) {
        	// printf("Thread ID: %d, %d\n", omp_get_thread_num(), i);
            Fx = Fy = 0;
            collision[i] = -1;
            mass = n_mass[i];

            if (x_coordinate[i] < left_boundary) {
                x_volocity[i] = -x_volocity[i];
                x_coordinate[i] = left_boundary + (left_boundary - x_coordinate[i]);
            } 
            else if (x_coordinate[i] > right_boundary) {
                x_volocity[i] = -x_volocity[i];
                x_coordinate[i] = right_boundary - (x_coordinate[i] - right_boundary);                  
            } 
            if (y_coordinate[i] < bottom_boundary) {
                y_volocity[i] = -y_volocity[i];
                y_coordinate[i] = bottom_boundary + (bottom_boundary - y_coordinate[i]);  
            }
            else if (y_coordinate[i] > up_boundary) {
                y_volocity[i] = -y_volocity[i];
                y_coordinate[i] = up_boundary - (y_coordinate[i] - up_boundary); 
            }

            double r, fx, fy, mass2;
            int has_collided = 0;
            // #pragma omp parallel num_threads(n_thread) reduction(+:Fx) reduction(+:Fy) private(r, fx, fy, mass2)
            for(int k = 0; k < num_n_body; k++) {
            	// printf("Thread ID: %d, %d\n", omp_get_thread_num(), k);
             //    has_collided = 0;
                if (k == i) continue;
                r = sqrt((x_coordinate[i]- x_coordinate[k])*(x_coordinate[i]- x_coordinate[k])+
                        ( y_coordinate[i]- y_coordinate[k])*(y_coordinate[i]- y_coordinate[k]));
                mass2 = n_mass[k];

                if (r < 2) {
                    if (i < k && has_collided == 0) {
                        collision[i] = k;
                        has_collided = 1;                    
                    }
                    r = 2;
                }
                fx = mass * mass2 * (x_coordinate[k]- x_coordinate[i])/(r*r*r);
                fy = mass * mass2 * (y_coordinate[k]- y_coordinate[i])/(r*r*r);
                Fx += fx;
                Fy += fy;
                
            }
            
            x_factor[i] = Fx/mass;
            y_factor[i] = Fy/mass;

            //see if collision happens
        }

        int k, has_collided;
        double mass, mass2, temp_x_vol, temp_y_vol, r;
        for(int i = 0; i < num_n_body; i++) {
            has_collided = 0;
            k = collision[i];

             XDrawPoint (display, win, gc, x_coordinate[i], y_coordinate[i]);
            // // usleep(1);
             printf("");
        
            if (k != -1) {

                mass = n_mass[i];
                mass2 = n_mass[k];

                temp_x_vol = x_volocity[i];
                temp_y_vol = y_volocity[i];
                              
                x_volocity[i] = ( (mass-mass2)*x_volocity[i]+2*mass2*x_volocity[k] ) /
                                  (mass+mass2);
                y_volocity[i] = ( (mass-mass2)*y_volocity[i]+2*mass2*y_volocity[k] ) /
                                  (mass+mass2);
                x_volocity[k] = ( (mass2-mass)*x_volocity[k]+2*mass*temp_x_vol ) /
                                  (mass+mass2);
                y_volocity[k] = ( (mass2-mass)*y_volocity[k]+2*mass*temp_y_vol ) /
                                  (mass+mass2);
            }

            x_volocity[i] = x_volocity[i] + x_factor[i]*t_interval;
            y_volocity[i] = y_volocity[i] + y_factor[i]*t_interval;
            x_coordinate[i] = x_coordinate[i] + x_volocity[i]*t_interval;
            y_coordinate[i] = y_coordinate[i] + y_volocity[i]*t_interval;
        }

        XFlush (display);
        times++;
 
    }

    gettimeofday(&t_end, NULL);
    double runningTime = (t_end.tv_sec - t_start.tv_sec) + 
                         (t_end.tv_usec - t_start.tv_usec) / 1000000.0;
    printInfo(runningTime);

    return 0;
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
    printf("Assignment 3, N Body, openMP implementation.\n");
    printf("runTime is %f\n", runningTime);
    printf("number of threads = %d", n_thread);
}