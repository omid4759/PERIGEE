#include "GenBC_UserLPM.hpp"

GenBC_UserLPM::GenBC_UserLPM( const char * const &lpn_filename, const int &in_N,
    const double &dt3d )
: N( in_N ), h( dt3d/static_cast<double>(N) ),
  absTol( 1.0e-8 ), relTol( 1.0e-5 )
{
  // Now read the lpn files for num_ebc, Ra, Ca, Ra_micro, cim, Rv and Pd
  std::string temp_name( lpn_filename );
  SYS_T::file_check( temp_name ); // make sure the file is on the disk

  std::ifstream reader;
  reader.open( lpn_filename, std::ifstream::in );

  std::istringstream sstrm;
  std::string sline;
  std::string bc_type;
  //const PetscMPIInt rank = SYS_T::get_MPI_rank();
  tstart=0.0;
  tend=N*h;
  // The first non-commented lien should be
  // RCR num_ebc
  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str(sline);
      sstrm >> bc_type;
      sstrm >> num_ebc;
      sstrm >> num_LPM_unknowns;
      sstrm.clear();
      break;
    }
  }

  // Check the file's bc_type matches RCR
  if( bc_type.compare("UserLPM") ==0
      || bc_type.compare("USERLPM") == 0
      || bc_type.compare("userlpm") == 0 )
  {
    Surface_To_LPM.resize(num_ebc);
    face_type.resize(num_ebc);
    prev_0D_sol.resize(num_LPM_unknowns);
    Xi0.resize(num_LPM_unknowns);
    Q0.resize(num_ebc);
    P0.resize(num_ebc);

  }
  else SYS_T::print_fatal("Error: the outflow model in %s does not match GenBC_UserLPM.\n", lpn_filename);




  // Read files for each ebc to set the values of Ra, Ca, Ra_micro, Cim,Rv, and Pd
  int counter = 0;

  while( std::getline(reader, sline) )
  {
    if( sline[0] != '#' && !sline.empty() )
    {
      sstrm.str( sline );
      int face_id;
      sstrm >> face_id;

      // Make sure the face_id, the first column in the file are listed
      // from 0 to ebc_id - 1
      if(face_id != counter) SYS_T::print_fatal("Error: GenBC_UserLPM the input file %s has wrong format in the face id column (the first column). \n", lpn_filename);

      sstrm >> face_type[ counter ];

      sstrm >> Surface_To_LPM[counter];

      sstrm.clear();

      counter += 1;
    }
  }


  if(counter != num_ebc ) SYS_T::print_fatal("Error: GenBC_UserLPM the input file %s does not contain complete data for GenBC faces. \n", lpn_filename);

  reader.close();

  SYS_T::commPrint( "===> GenBC_UserLPM data are read in from %s.\n", lpn_filename );

  // Set a zero initial value. They should be reset based on the initial
  // 3D solutions.
  for(int ii=0; ii<num_ebc; ++ii)
  {

      Q0[ii] = 0.0;

  }

  for(int ii=0;ii<num_LPM_unknowns;++ii){
     Xi0[ii]=0.0;
     prev_0D_sol[ii]=0.0;
  }

  if(SYS_T::get_MPI_rank() == 0){
   remove( myfifo1) ;

   remove( myfifo2) ;

	int creation = mkfifo(myfifo1, 0666);
    if(creation==-1){
      perror("mkfifo error");
      SYS_T::print_fatal("mkfifo error.\n");
    }
    creation = mkfifo(myfifo2, 0666);
    if(creation==-1){
      perror("mkfifo error");
      SYS_T::print_fatal("mkfifo error.\n");
    }
  }


}

GenBC_UserLPM::~GenBC_UserLPM()
{}

void GenBC_UserLPM::print_info() const
{
  SYS_T::commPrint( "     User defined LPM: N = %d, h = %e, num_ebc = %d \n", N, h, num_ebc );

}


void GenBC_UserLPM::get_m( double * const &in_dot_Q, double * const &in_Q,
 double * const &m ) const
{


  double * in_Q_1=new double[num_ebc];
  double * in_Q_2=new double[num_ebc];
  double * left=new double[num_ebc];
  double * right=new double[num_ebc];
  double * diff=new double[num_ebc];
  for (int ii=0;ii<num_ebc;ii++){
   diff[ii] = std::abs(in_Q[ii]) * relTol;
   if( diff[ii] < absTol ) diff[ii] = absTol;
   in_Q_1[ii]=in_Q[ii]+0.5*diff[ii];
   in_Q_2[ii]=in_Q[ii]-0.5*diff[ii];
  }
   get_P(in_dot_Q,in_Q_1,left);
   get_P(in_dot_Q,in_Q_2,right);

  for(int ii =0; ii<num_ebc;++ii){

   m[ii]=(left[ii] - right[ii]) / diff[ii];
  }

  delete [] in_Q_1;
  in_Q_1=nullptr;
  delete [] in_Q_2;
  in_Q_2=nullptr;
  delete [] left;
  left=nullptr;
  delete [] right;
  right=nullptr;
  delete [] diff;
  diff=nullptr;

}


void GenBC_UserLPM::get_P0(double * const &Pn) const
{
  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  for ( int myrank = 0; myrank < size; ++myrank ){
   if(myrank==rank){

   std::string str = "./get_face_value &";
  // str = str + std::to_string(ii);

    // Convert string to const char * as system requires
    // parameter of type const char *
    const char *command = str.c_str();


   system(command);

    int fd1,fd2;

    fd1 = open(myfifo1, O_WRONLY ); //open(myfifo, O_WRONLY | O_NONBLOCK);
  //  fd2 = open(myfifo2, O_WRONLY );
    if(fd1==-1){
      perror("open error \n");
	  SYS_T::print_fatal("open myfifo1 failed \n");
    }

    double * LPM_sols=new double[num_LPM_unknowns];
    double * Q_faces=new double[num_ebc];

    for(int ii=0;ii<num_LPM_unknowns;++ii){
     LPM_sols[ii]=Xi0[ii];
    }
    for(int ii=0;ii<num_ebc;++ii){
     Q_faces[ii]=Q0[ii];
    }

   // printf("get P0 send LPM_sols=%lf %lf Q=%lf %lf \n",LPM_sols[0],LPM_sols[1],Q_faces[0],Q_faces[1]);

    write(fd1,LPM_sols,sizeof(double)*num_LPM_unknowns);
    write(fd1,Q_faces,sizeof(double)*num_ebc);
    close(fd1);

    fd2 = open(myfifo2, O_RDONLY);

     // Read from FIFO
    read(fd2,Pn,sizeof(double)*num_ebc);

    // Print the read message

    close(fd2);
 //   printf("get P0, Pn= %lf %lf \n",Pn[0],Pn[1]);
   // printf("for comparison QR+x=%lf %lf \n",Q0[0]*1000.0+Xi0[0],Q0[1]*1000.0+Xi0[1]);
   // printf("get_P0 values received going to delete [] rank=%d num_LPM_unknows=%d \n",rank,num_LPM_unknowns);

    delete [] LPM_sols; LPM_sols=nullptr;
    delete [] Q_faces; Q_faces=nullptr;

    }

    MPI_Barrier(PETSC_COMM_WORLD);
   }



   // Pn[0]=Q0[0]*1000.0+Xi0[0];
   // Pn[1]=Q0[1]*1000.0+Xi0[2];
   // printf("usr LPM get_P0 received =%lf %lf  \n",Pn[0],Pn[1]);
}


void GenBC_UserLPM::get_P( double *const &in_dot_Q,
   double * const &in_Q, double * const &P ) const
{


  const PetscMPIInt rank = SYS_T::get_MPI_rank();
  const PetscMPIInt size = SYS_T::get_MPI_size();

  for ( int myrank = 0; myrank < size; ++myrank ){
   if(myrank==rank){

//    clock_t clock_start = clock();
//    clock_t clock_stop;
//    double tsec;
    double *Xi_m=new double[num_LPM_unknowns];
    double * in_Q0=new double[num_ebc];

    for(int ii=0;ii<num_LPM_unknowns;++ii){
      Xi_m[ii]= Xi0[ii];
    }
    for(int ii=0;ii<num_ebc;++ii){
     in_Q0[ii]=Q0[ii];
    }


    std::string str = "./GenBC_User &";

    // Convert string to const char * as system requires
    // parameter of type const char *
    const char *command = str.c_str();

    system(command);


    int fd1,fd2;

    /*
	// FIFO file path
	const char * myfifo1 = "/tmp/myfifo1";
    const char * myfifo2 = "/tmp/myfifo2";
	// Creating the named file(FIFO)
	// mkfifo(<pathname>, <permission>)
	mkfifo(myfifo1, 0666);
	mkfifo(myfifo2, 0666);
*/
    fd1 = open(myfifo1, O_WRONLY ); //open(myfifo, O_WRONLY | O_NONBLOCK);
   // fd2 = open(myfifo2, O_WRONLY );

   // printf("output LPM sols \n");
  //  printf("get P send Xi_m=%lf %lf inQ=%lf %lf \n",Xi_m[0],Xi_m[1],in_Q[0],in_Q[1]);
    write(fd1,Xi_m,sizeof(double)*num_LPM_unknowns);
   //  write(fd,num_ebc,sizeof(int));

    write(fd1,in_Q0,sizeof(double)*num_ebc);
    write(fd1,in_Q,sizeof(double)*num_ebc);
    write(fd1,&tstart,sizeof(double));
    write(fd1,&tend,sizeof(double));



    close(fd1);


    fd2 = open(myfifo2, O_RDONLY);

     // Read from FIFO


    read(fd2, Xi_m, sizeof(double)*num_LPM_unknowns);
    read(fd2, P, sizeof(double)*num_ebc);
    close(fd2);

    for (int ii=0; ii<num_LPM_unknowns;++ii){
     prev_0D_sol[ii]=Xi_m[ii];
    }
    delete [] Xi_m; Xi_m=nullptr;
    delete [] in_Q0; in_Q0=nullptr;
 //   clock_stop = clock();
 //   tsec = ((double) (clock_stop-clock_start)/CLOCKS_PER_SEC );
 //   printf("UsrLPM CPU time: %lf secs\n", tsec);
     }
     MPI_Barrier(PETSC_COMM_WORLD);
    }


}


void GenBC_UserLPM::reset_initial_sol( double *const &in_Q_0,
    double * const  &in_P_0, const double &curr_time )
{


  //clock_t clock_start = clock();
  //clock_t clock_stop;
  //double tsec;
  for(int ii=0;ii<num_ebc;++ii){
   Q0[ii]  = in_Q_0[ii];
   }

  //when reset_initial_sol is called, use the last converged 0D solition as initial solutions.

  for(int ii=0;ii<num_LPM_unknowns;++ii){
    Xi0[ii]= prev_0D_sol[ii];
  }

  tstart=curr_time;
  tend=curr_time+N*h;


 // clock_stop = clock();
 // tsec = ((double) (clock_stop-clock_start)/CLOCKS_PER_SEC );
 // printf("reset initial_sol =%d CPU time: %lf secs\n",ii, tsec);
}

double GenBC_UserLPM:: F(double * const &pi, double * const &q, double * const &K)const
{

   std::string str = "./GenBC_User & ";
  // str = str + std::to_string(ii);

    // Convert string to const char * as system requires
    // parameter of type const char *
    const char *command = str.c_str();

 //  printf("F Calling external program  \n");
    system(command);




    int fd1,fd2;
    /*
	// FIFO file path
	const char * myfifo1 = "/tmp/myfifo1";
    const char * myfifo2 = "/tmp/myfifo2";
	// Creating the named file(FIFO)
	// mkfifo(<pathname>, <permission>)
	mkfifo(myfifo1, 0666);
    mkfifo(myfifo2, 0666);

    */
    fd1 = open(myfifo1, O_WRONLY ); //open(myfifo, O_WRONLY | O_NONBLOCK);
  //  fd2 = open(myfifo2, O_WRONLY );

   // printf("F  pi=%lf %lf q=%lf %lf \n",pi[0],pi[1],q[0],q[1]);
    write(fd1,pi,sizeof(double)*num_LPM_unknowns);

    write(fd1,q,sizeof(double)*num_ebc);
  //   printf("F finish writing q  \n");
    close(fd1);


    mkfifo(myfifo2, 0666);
    fd2 = open(myfifo2, O_RDONLY);

     // Read from FIFO
/*
    for(int ii=0;ii<num_LPM_unknowns;++ii){
      read(fd, &K[ii], sizeof(double));
    }
*/

    read(fd2,K,sizeof(double)*num_LPM_unknowns);
    close(fd2);


   // printf("F complete \n");


      return 0.0;
}


// EOF
