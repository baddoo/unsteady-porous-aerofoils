function x = myJacobiNodes ( n, alpha, beta )

%*****************************************************************************80
%
%% J_POLYNOMIAL_ZEROS: zeros of Jacobi polynomial J(n,a,b,x).
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    19 March 2012
%
%  Author:
%
%    John Burkardt.
%
%  Reference:
%
%    Sylvan Elhay, Jaroslav Kautsky,
%    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
%    Interpolatory Quadrature,
%    ACM Transactions on Mathematical Software,
%    Volume 13, Number 4, December 1987, pages 399-415.
%
%  Parameters:
%
%    Input, integer N, the order.
%
%    Input, real ALPHA, BETA, the parameters.
%    -1 < ALPHA, BETA.
%
%    Output, real X(N), the zeros.
%
  ab = alpha + beta;
  abi = 2.0 + ab;
%
%  Define the zero-th moment.
%
  zemu = 2.0^( ab + 1.0 ) * double ( gamma ( sym( alpha + 1.0 ))) ...
    * double (gamma (sym( beta + 1.0 ))) / double ( gamma ( sym ( abi )));
%
%  Define the Jacobi matrix.
%
  x = zeros ( n, 1 );
  bj = zeros ( n, 1 );

  x(1) = ( beta - alpha ) / abi;
  bj(1) = 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) ...
    / ( ( abi + 1.0 ) * abi * abi );
  a2b2 = beta * beta - alpha * alpha;

  for j = 2 : n
    abi = 2.0 * j + ab;
    x(j) = a2b2 / ( ( abi - 2.0 ) * abi );
    abi = abi^2;
    bj(j) = 4.0 * j * ( j + alpha ) * ( j + beta ) * ( j + ab ) ...
      / ( ( abi - 1.0 ) * abi );
  end
  bj(1:n) =  sqrt ( bj(1:n) );

  w = zeros ( n, 1 );
  w(1) = sqrt ( zemu );
%
%  Diagonalize the Jacobi matrix.
%
  [ x, w ] = imtqlx ( n, x, bj, w );

  return
end