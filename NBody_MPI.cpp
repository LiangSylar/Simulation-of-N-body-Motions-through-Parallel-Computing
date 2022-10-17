/* This MPI program is written to solve the bug for collision */
/* The program need to make sure each body collide with another body within the same processor */
#include "mpi.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xos.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>

#define         X_RESN  200       /* x resolution */
#define         Y_RESN  200       /* y resolution */
#define         t_interval 1
#define         num_n_body 500
#define         g_const 1.0
#define         MASTER 0
#define         total_times 1000

void   initRandomSeed(); 
double randomReal(double low, double high);
double *randomRealSets (int n, double low, double high);
void   printInfo(double runningTime);

int main (int argc, char* argv[])
{

    MPI_Init(&argc, & argv);
    int taskid, 
        numtasks, 
        len;
    char hostname[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(hostname, &len);
    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &taskid);

    int average_tasks1 =  int( floor(num_n_body/(numtasks-1)) );
    int average_tasks2 = num_n_body - (numtasks - 2)*average_tasks1; 
    int average_tasks = (taskid == numtasks-1) ? average_tasks2 : average_tasks1;

    double up_boundary = Y_RESN*(2.0/3);
    double bottom_boundary = Y_RESN*(1.0/3);
    double left_boundary = X_RESN*(1.0/3);
    double right_boundary = X_RESN*(2.0/3);

    if (taskid == MASTER) { // Master receives position of points.
        
        Window          win;                            /* initialization for a window */
        unsigned
        int             width, height,                  /* window size */
                        x, y,                           /* window position */
                        border_width,                   /*border width in pixels */
                        display_width, display_height,  /* size of screen */
                        screen;                         /* which screen */

        char            *window_name = "N Body", *display_name = NULL;
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
        gettimeofday(&t_start, NULL);  

        double up_boundary = Y_RESN*(2.0/3);
        double bottom_boundary = Y_RESN*(1.0/3);
        double left_boundary = X_RESN*(1.0/3);
        double right_boundary = X_RESN*(2.0/3);
        double *x_coordinate = randomRealSets(num_n_body, left_boundary, right_boundary);
        double *y_coordinate = randomRealSets(num_n_body, bottom_boundary, up_boundary);
        double *n_mass       = randomRealSets(num_n_body, 0.5, 1.5);
        double *x_volocity   = randomRealSets(num_n_body, 0, 0);
        double *y_volocity   = randomRealSets(num_n_body, 0, 0);

        double x_factor[num_n_body]; 
        double y_factor[num_n_body];
        int collision [num_n_body];

        MPI_Bcast(n_mass, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

        int recv_average_tasks;
        MPI_Status status;
        int times = 1;
        while (times < total_times) { 

            // printf("times: %d\n", times);

            XClearWindow(display, win);
            MPI_Bcast(x_coordinate, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(y_coordinate, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(x_volocity, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(y_volocity, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            
            int start;
            for(int i = 1; i < numtasks; i++) {
                
                recv_average_tasks = (i == numtasks-1) ? average_tasks2 : average_tasks1;
                
                double recv_x_coordinate[recv_average_tasks];
                double recv_y_coordinate[recv_average_tasks];
                double recv_x_factor    [recv_average_tasks];
                double recv_y_factor    [recv_average_tasks];
                double recv_x_volocity  [recv_average_tasks];
                double recv_y_volocity  [recv_average_tasks];
                int    recv_collision   [recv_average_tasks];

                MPI_Recv(recv_x_coordinate, recv_average_tasks, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(recv_y_coordinate, recv_average_tasks, MPI_DOUBLE, i, 2, MPI_COMM_WORLD, &status);
                MPI_Recv(recv_x_volocity, recv_average_tasks, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
                MPI_Recv(recv_y_volocity, recv_average_tasks, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &status);
                MPI_Recv(recv_x_factor, recv_average_tasks, MPI_DOUBLE, i, 6, MPI_COMM_WORLD, &status);
                MPI_Recv(recv_y_factor, recv_average_tasks, MPI_DOUBLE, i, 7, MPI_COMM_WORLD, &status);
                
                MPI_Recv(recv_collision, recv_average_tasks, MPI_INT, i, 5, MPI_COMM_WORLD, &status);
                start = (i-1)*average_tasks1;
                for(int k = 0; k < recv_average_tasks; k++) {                   
                    x_coordinate[k+start] = recv_x_coordinate[k];
                    y_coordinate[k+start] = recv_y_coordinate[k];
                    x_volocity  [k+start] = recv_x_volocity[k];
                    y_volocity  [k+start] = recv_y_volocity[k];
                    // printf("%f, %f \n", x_volocity[k+start], recv_x_volocity[k]);
                    collision   [k+start] = recv_collision[k];
                    x_factor    [k+start] = recv_x_factor [k];
                    y_factor    [k+start] = recv_y_factor [k];
                }
            }

            int k;
            double mass, mass2, temp_x_vol, temp_y_vol;
            for(int i = 0; i < num_n_body; i++) {

                k = collision[i];
                
                if (k != -1) {

                    mass  = n_mass[i];                    
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

                XDrawPoint (display, win, gc, x_coordinate[i], y_coordinate[i]);
                // // usleep(1);
                printf("");
            }
    
            XFlush(display);
            times++;
        }
        gettimeofday(&t_end, NULL);
        double runningTime = (t_end.tv_sec - t_start.tv_sec) + 
                             (t_end.tv_usec - t_start.tv_usec) / 1000000.0;
        printInfo(runningTime);

    }

    else { // slaves send back position of the points.

        double Srecv_x_coordinate [num_n_body]; 
        double Srecv_y_coordinate [num_n_body];
        double Srecv_n_mass       [num_n_body];
        double Srecv_x_volocity   [num_n_body];
        double Srecv_y_volocity   [num_n_body];

        double Ssend_x_coordinate [average_tasks];
        double Ssend_y_coordinate [average_tasks];
        double Ssend_x_volocity   [average_tasks];
        double Ssend_y_volocity   [average_tasks];
        double Ssend_x_factor     [average_tasks];
        double Ssend_y_factor     [average_tasks];
        int    Ssend_collision    [average_tasks];

        MPI_Bcast(Srecv_n_mass, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

        double mass,  mass2, Fx, Fy, fx, fy, r, has_collided;

        int start = (taskid-1)*average_tasks1;
        int Stimes = 1;
        while (Stimes < total_times) {

            MPI_Bcast(Srecv_x_coordinate, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(Srecv_y_coordinate, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(Srecv_x_volocity, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);
            MPI_Bcast(Srecv_y_volocity, num_n_body, MPI_DOUBLE, MASTER, MPI_COMM_WORLD);

            for(int i = start; i < average_tasks+start; i++) {
                Fx = Fy = 0;   
                has_collided = 0;
                Ssend_collision[i-start] = -1;
                Ssend_x_volocity[i-start] = Srecv_x_volocity[i];
                Ssend_y_volocity[i-start] = Srecv_y_volocity[i];
                Ssend_x_coordinate[i-start] = Srecv_x_coordinate[i];
                Ssend_y_coordinate[i-start] = Srecv_y_coordinate[i];
                mass      = Srecv_n_mass[i];

                if (Srecv_x_coordinate[i] < left_boundary) {
                    Ssend_x_volocity[i-start] = -Srecv_x_volocity[i];
                    Ssend_x_coordinate[i-start] = left_boundary + (left_boundary - Srecv_x_coordinate[i]);
                } 
                else if (Srecv_x_coordinate[i] > right_boundary) {
                    Ssend_x_volocity[i-start] = -Srecv_x_volocity[i];
                    Ssend_x_coordinate[i-start] = right_boundary - (Srecv_x_coordinate[i] - right_boundary);                  
                } 
                if (Srecv_y_coordinate[i] < bottom_boundary) {
                    Ssend_y_volocity[i-start] = -Srecv_y_volocity[i];
                    Ssend_y_coordinate[i-start] = bottom_boundary + (bottom_boundary - Srecv_y_coordinate[i]);
                }
                else if (Srecv_y_coordinate[i] > up_boundary) {
                    Ssend_y_volocity[i-start] = -Srecv_y_volocity[i];
                    Ssend_y_coordinate[i-start] = up_boundary - (Srecv_y_coordinate[i] - up_boundary); 
                }

                /* collision between particles */
            
                for(int k = 0; k < num_n_body; k++) { 
                    if (k == i) continue;
                    r = sqrt((Srecv_x_coordinate[i]- Srecv_x_coordinate[k])*
                            ( Srecv_x_coordinate[i]- Srecv_x_coordinate[k])+
                            ( Srecv_y_coordinate[i]- Srecv_y_coordinate[k])*
                            ( Srecv_y_coordinate[i]- Srecv_y_coordinate[k]));
                    mass2 = Srecv_n_mass[k];
                    if (r < 2) {
                        if (i < k && has_collided == 0) {
                            Ssend_collision[i-start] = k;
                            has_collided = 1;
                        }
                        r = 2;
                    }
                    fx = mass * mass2 * (Srecv_x_coordinate[k]- Srecv_x_coordinate[i])/(r*r*r);
                    fy = mass * mass2 * (Srecv_y_coordinate[k]- Srecv_y_coordinate[i])/(r*r*r);                                 
                    Fx += fx;
                    Fy += fy;
                }

                Ssend_x_factor[i-start] = Fx/mass;
                Ssend_y_factor[i-start] = Fy/mass;
            }

            Stimes++;
            MPI_Send(Ssend_x_coordinate, average_tasks, MPI_DOUBLE, MASTER, 1, MPI_COMM_WORLD);
            MPI_Send(Ssend_y_coordinate, average_tasks, MPI_DOUBLE, MASTER, 2, MPI_COMM_WORLD);
            MPI_Send(Ssend_x_volocity, average_tasks, MPI_DOUBLE, MASTER, 3, MPI_COMM_WORLD);
            MPI_Send(Ssend_y_volocity, average_tasks, MPI_DOUBLE, MASTER, 4, MPI_COMM_WORLD);
            MPI_Send(Ssend_collision, average_tasks, MPI_INT, MASTER, 5, MPI_COMM_WORLD);
            MPI_Send(Ssend_x_factor, average_tasks, MPI_DOUBLE, MASTER, 6, MPI_COMM_WORLD);
            MPI_Send(Ssend_y_factor, average_tasks, MPI_DOUBLE, MASTER, 7, MPI_COMM_WORLD);

        }

    }
    printf("number of cores = %d", numtasks);
    MPI_Finalize();
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
    initRandomSeed();
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
    printf("Assignment 3, N Body, MPI implementation.\n");
    printf("runTime is %f\n", runningTime);
}
