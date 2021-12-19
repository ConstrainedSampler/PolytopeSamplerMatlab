function value = r8_lgit ( a, x, algap1 )

%*****************************************************************************80
%
%% R8_LGIT evaluates the log of Tricomi's incomplete gamma function.
%
%  Discussion:
%
%    Perron's continued fraction is used for large X and X <= A.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    28 September 2011
%
%  Author:
%
%    Original FORTRAN77 version by Wayne Fullerton.
%    MATLAB version by John Burkardt.
%
%  Reference:
%
%    Wayne Fullerton,
%    Portable Special Function Routines,
%    in Portability of Numerical Software,
%    edited by Wayne Cowell,
%    Lecture Notes in Computer Science, Volume 57,
%    Springer 1977,
%    ISBN: 978-3-540-08446-4,
%    LC: QA297.W65.
%
%  Input:
%
%    real A, the parameter.
%
%    real X, the argument.
%
%    real ALGAP1, the logarithm of the gamma function of A+1.
%
%  Output:
%
%    real VALUE, the log of Tricomi's incomplete gamma function.
%
  persistent eps
  persistent sqeps

  if ( isempty ( eps ) )
    eps = 0.5 * r8_mach ( 3 );
    sqeps = sqrt ( r8_mach ( 4 ) );
  end

  if ( x <= 0.0 )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8_LGIT - Fatal error!\n' );
    fprintf ( 1, '  X <= 0.\n' );
    error ( 'R8_LGIT - Fatal error!' )
  end

  if ( a < x )
    fprintf ( 1, '\n' );
    fprintf ( 1, 'R8_LGIT - Fatal error!\n' );
    fprintf ( 1, '  A < X.\n' );
    error ( 'R8_LGIT - Fatal error!' )
  end

  ax = a + x;
  a1x = ax + 1.0;
  r = 0.0;
  p = 1.0;
  s = p;

  for k = 1 : 200
    fk = k;
    t = ( a + fk ) * x * ( 1.0 + r );
    r = t / ( ( ax + fk ) * ( a1x + fk ) - t );
    p = r * p;
    s = s + p;
    if ( abs ( p ) < eps * s )
      hstar = 1.0 - x * s / a1x;
      value = - x - algap1 - log ( hstar );
      return
    end
  end

  fprintf ( 1, '\n' );
  fprintf ( 1, 'R8_LGIT - Fatal error!\n' );
  fprintf ( 1, '  No convergence after 200 iterations.\n' );
  error ( 'R8_LGIT - Fatal error!' )
end
