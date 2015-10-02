#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <math.h>
#define VECSIZE 64000
#define ITERATIONS 100
#define NUMDIM 5
double When()
{
  struct timeval tp;
  gettimeofday(&tp, NULL);
  return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

struct INFO
{
        double val;
        int rank;
};

void Reduce(double param[VECSIZE], int numdim, int rank, MPI_Status status)
{
        int curdim, i;
        int notparticipating =0;
        int bitmask =1;
        double compare[VECSIZE];
        for(curdim=0;curdim<numdim;curdim++)
        {
                if((rank & notparticipating) == 0){
                        if((rank & bitmask) != 0){
                                int mesg_dest = rank^bitmask;
                                MPI_Send(param,VECSIZE, MPI_INT, mesg_dest,0,MPI_COMM_WORLD);
                        }
                        else{
                                int msg_src = rank^bitmask;
                                MPI_Recv(compare,VECSIZE,MPI_INT, msg_src,0,MPI_COMM_WORLD,&status);
                                for(i=0; i<VECSIZE;i++)
                                {
                                        if(compare[i]>param[i]) param[i] = compare[i];
                                }
                        }
                }
                notparticipating = notparticipating^bitmask;
                bitmask<<=1;
        }
}

void Broadcast(double param[VECSIZE], int numdim, int rank, MPI_Status status)
{
        double dnumdim = numdim;
        int notparticipating = pow(2,dnumdim-1)-1;
        int bitmask = pow(2,dnumdim-1);
        double compare[VECSIZE];
        int curdim,i;
        for(curdim=0;curdim<numdim;curdim++)
        {
                if((rank & notparticipating) ==0){
                        if((rank & bitmask) ==0){
                                int msg_dest = rank^bitmask;
                                MPI_Send(param,VECSIZE, MPI_INT, msg_dest,0,MPI_COMM_WORLD);
                        }else{
                                int msg_src = rank^bitmask;
                                MPI_Recv(compare,VECSIZE,MPI_INT, msg_src,0,MPI_COMM_WORLD,&status);
                                for(i =0; i<VECSIZE;i++)
                                {
                                        param[i] = compare[i];
                                }
                        }
           }
                notparticipating >>=1;
                bitmask >>=1;
        }

}

main(int argc, char *argv[])
{
        int iproc, nproc,i, iter;
        char host[255], message[55];
        MPI_Status status;

        MPI_Init(&argc, &argv);
        MPI_Comm_size(MPI_COMM_WORLD, &nproc);
        MPI_Comm_rank(MPI_COMM_WORLD, &iproc);

        gethostname(host,253);

        double vector[VECSIZE];
        int myrank, root = 0;

        MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
        // Start time here
        srand(myrank+5);
        double start = When();
        for(iter = 0; iter < ITERATIONS; iter++) {
          for(i = 0; i < VECSIZE; i++) {
            vector[i] = rand();
          }
          Reduce(vector,NUMDIM, myrank,status);
          //MPI_Reduce( in, out, VECSIZE, MPI_DOUBLE_INT, MPI_MAXLOC, root, MPI_COMM_WORLD);
          // At this point, the answer resides on process root
          if (myrank == root) {
              /* read ranks out
               */
              //printf("the root reduced to %f\n",vector[0]);
          }
          // Now broadcast this max vector to everyone else.
          Broadcast(vector, NUMDIM, myrank,status);
          //MPI_Bcast(vector, VECSIZE, MPI_DOUBLE_INT, root, MPI_COMM_WORLD);
          for(i = 0; i < VECSIZE; i++) {
          //printf("final proc %d [%d]=%f from %d\n",myrank,i,out[i].val,out[i].rank);
          }
        }
        MPI_Finalize();
        double end = When();
        if(myrank == root) {
          printf("Time %f\n",end-start);
        }
}
                                                        
