typedef std::vector< double > DVec

// take a matrix and return its sqare root
void corset( const TMatrixD& V, TMatrixD& C );

// Given the sqrt of the covariance matrix, generate correlated random numbers.  This assumes the mean of all variables is 0.  If it is not, you'll 
void corgen( const TMatrixD& C, DVec& x );

void yourRoutine (const TMatrix &covarMatrix)
{
      TMatrixD sqrtCov( NP, NP );
      corset( cov, sqrtCov );

      for (int loop = 0; loop < 10000; ++loop)
      {
          vector< double > values (0, NP);
          corgen (sqrtCov, values);
          // values assumes all means are at 0.  If they are not, you'll have too add the means to each element now.

         // use values
      } // for loop
}


void corset( const TMatrixD& V, TMatrixD& C )
{ 
  // calculate sqrt(V) as lower diagonal matrix
  for( int i = 0; i < NP; ++i ) {
    for( int j = 0; j < NP; ++j ) {
      C[i][j] = 0;
    }
  }

  for( int j = 0; j < NP; ++j ) {
    // diagonal terms first
    double Ck = 0;
    for( int k = 0; k < j; ++k ) {
      Ck += C[j][k] * C[j][k];
    } // k 
    C[j][j] = sqrt( fabs( V[j][j] - Ck ) );

    // off-diagonal terms
    for( int i = j+1; i < NP; ++i ) {
      Ck = 0;
      for( int k = 0; k < j; ++k ) {
	Ck += C[i][k] * C[j][k];
      } //k
      C[i][j] = ( V[i][j] - Ck ) / C[j][j];
    }// i 
  } // j
}

void corgen( const TMatrixD& C, DVec& x )
{
  DVec z (0, NP);

  // np random numbers from unit Gaussian
  for( int i = 0; i < NP; ++i ) {
    z[i] = gRandom->Gaus( 0.0, 1.0 );
  }

  for( int i = 0; i < NP; ++i ) {
    x[i] = 0;
    for( int j = 0; j <= i; ++j ) {
      x[i] += C[i][j] * z[j];
    } // j
  } // i

}