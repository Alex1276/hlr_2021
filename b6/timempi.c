#include <stdio.h>
#include <mpi.h>

int
main (int argc, char** argv) {
  int rank;
  int size;
  char name[MPI_MAX_PROCESSOR_NAME];
  int length;

  // INIT MPI environment 
  MPI_Init(&argc,&argv);

  // get tis process' rank
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);

  // get number of ranks 
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  //get name of processor 
  MPI_Get_processor_name(name,&length);

  int buffer_len = 150;
  char buffer[buffer_len];

  // print from rank 0
  // else send to rank 0 and print there
  if (rank == 0) {
    //sprintf(buffer," Rank %d beendet jetzt!\n",rank);
    // print rank 0 message
    printf("%s",buffer);
    // for each rank != 0
    for (int i = 1; i<size;i++) {
      // receive message from rank
      MPI_Recv(buffer,buffer_len,MPI_CHAR,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);  
      //print recieved message
      printf("%s\n",buffer);
    }
  }else{
    double time = MPI_Wtime();
    sprintf(buffer,"%s: %fl",name,time);
    // send message to rank 0
    MPI_Send(buffer,buffer_len,MPI_CHAR,0,rank,MPI_COMM_WORLD);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  printf("RAnk: %d  beendet jetzt! \n",rank);
  MPI_Finalize();
  return 0;
}