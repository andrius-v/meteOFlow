<?

define ("BOUNDARY_CONDITION", 0);
//0 Neumann
//1 Periodic
//2 Symmetric

define ("ZOOM_SIGMA_ZERO", 0.6);

define ("BOUNDARY_CONDITION_DIRICHLET", 0);
define ("BOUNDARY_CONDITION_REFLECTING", 1);
define ("BOUNDARY_CONDITION_PERIODIC", 2);

define ("DEFAULT_GAUSSIAN_WINDOW_SIZE", 5);
define ("DEFAULT_BOUNDARY_CONDITION", BOUNDARY_CONDITION_REFLECTING);

define ("MAX_ITERATIONS", 300);
define ("PRESMOOTHING_SIGMA", 0.8);
define ("GRAD_IS_ZERO", 1E-10);

define ("PAR_DEFAULT_NPROC", 0);
define ("PAR_DEFAULT_TAU", 0.25);
define ("PAR_DEFAULT_LAMBDA", 0.15);
define ("PAR_DEFAULT_THETA", 0.3);
define ("PAR_DEFAULT_NSCALES", 100);
define ("PAR_DEFAULT_ZFACTOR", 0.5);
define ("PAR_DEFAULT_NWARPS", 5);
define ("PAR_DEFAULT_EPSILON", 0.01);
define ("PAR_DEFAULT_VERBOSE", 0);

class TVL1{
	
	private $u;
	private $v;
	public $time;
	
	function __construct(&$img0, &$img1, $scl)
	{
		$nproc   = PAR_DEFAULT_NPROC;
		$tau     = PAR_DEFAULT_TAU;
		$lambda  = PAR_DEFAULT_LAMBDA;
		$theta   = PAR_DEFAULT_THETA;
		$nscales = PAR_DEFAULT_NSCALES;
		$zfactor = PAR_DEFAULT_ZFACTOR;
		$nwarps  = PAR_DEFAULT_NWARPS;
		$epsilon = PAR_DEFAULT_EPSILON;
		$verbose = PAR_DEFAULT_VERBOSE;
		
		// read the input images
		//int    nx, ny, nx2, ny2;
		$nx = imagesx($img0);
		$ny = imagesy($img0);
		$nx2 = imagesx($img1);
		$ny2 = imagesy($img1);
		$this->RGB2Y($img0, $I0, $nx, $ny);
		$this->RGB2Y($img1, $I1, $nx2, $ny2);
		
		$N = 1 + log(hypot($nx, $ny) / 16.0) / log(1 / $zfactor);
		if ($N < $nscales) $nscales = intval($N);
		
		//allocate memory for the flow
		//float *u = xmalloc(2 * nx * ny * sizeof*u);
		//float *v = u + nx*ny;;
		
		$this->u = new SplFixedArray($nx * $ny);
		$this->v = new SplFixedArray($nx * $ny);
		$time_start = microtime(true);
		$this->Dual_TVL1_optic_flow_multiscale($I0, $I1, $this->u, $this->v, $nx, $ny, $tau, $lambda, $theta, $nscales, $zfactor, $nwarps, $epsilon, $verbose, $scl);
		$time_end = microtime(true);
		$this->time = $time_end - $time_start;
		//echo ('time = '); echo($time);
	}
	
	private function RGB2Y(&$img, &$I, $x, $y)
	{
		for ($j = 0; $j < $y; $j++)
		{
			for ($i = 0; $i < $x; $i++)
			{
				$c = imagecolorat($img, $i, $j);
				$r = ($c >> 16) & 0xFF;
				$g = ($c >> 8) & 0xFF;
				$b = $c & 0xFF;
				$I[$j * $x + $i] = 0.299 * $r + 0.587 * $g + 0.114 * $b;
			}
		}
	}
	
	public function Y2RGB(&$img, $scl, $th, $clr, $stp)
	{
		
		$w = imagesx($img);
		$h = imagesy($img);
		imageantialias ($img , true);
		for ($j = 0; $j < $h; $j = $j + $stp)
		{
			for ($i = 0; $i < $w; $i = $i + $stp)
			{
				$d = hypot($scl * $this->u[$j * $w + $i], $scl * $this->v[$j * $w + $i]);
				if (($d > 1) && ($d > ($scl / 2)))
				{
					$c = imagecolorallocate($img, $clr, $clr, $clr);
					//imagesetthickness($img, $th);
					$dx = $scl * $this->u[$j * $w + $i];
					$dy = $scl * $this->v[$j * $w + $i];
					imageline($img, $i, $j, $i + $dx, $j + $dy, $c);
					//$c = imagecolorallocate($img, $clr, 0, 0);
					//imagesetpixel($img, $i, $j, $c);
					imagesetpixel($img, $i + $dx - 1, $j + $dy, $c);
					imagesetpixel($img, $i + $dx , $j + + $dy - 1, $c);
					imagesetpixel($img, $i + $dx  + 1, $j + $dy, $c);
					imagesetpixel($img, $i + $dx , $j + $dy + 1, $c);
				}
			}
		}
		
		//return ($img);
	}
	
	public function isDone()
	{
		echo('Woot woot - we are done');
	}
	
	/*public function Y2RGB()
	{
		$img = imagecreatetruecolor($this->w, $this->h);
		for ($j = 0; $j < $this->h; $j++)
		{
			for ($i = 0; $i < $this->w; $i++)
			{
				$l = $this->I1[$j * $this->w + $i];
				$c = imagecolorallocate($img, $l, $l, $l);
				imagesetpixel($img, $i, $j, $c);
				//imageline($img, $i, $j, $i + $this->u[$j * $this->w + $i], $j + $this->v[$j * $this->w + $i]);
			}
		}
		
		for ($j = 0; $j < $this->h; $j = $j + 10)
		{
			for ($i = 0; $i < $this->w; $i = $i + 10)
			{
				//$l = $this->I1[$j * $this->w + $i];
				//$c = imagecolorallocate($img, $l, $l, $l);
				$c = imagecolorallocate($img, 0, 0, 0);
				//imagesetpixel($img, $i+$this->u, $j + $this->v, $c);
				imageline($img, $i, $j, $i + $this->u[$j * $this->w + $i]/10000, $j + $this->v[$j * $this->w + $i]/10000, $c);
			}
		}
		
		return ($img);
	}*/
	
	/**
	  *
	  * Neumann boundary condition test
	  *
	**/	
	public function neumann_bc($x, $nx, &$out)
	{
		if($x < 0)
		{
			$x = 0;
			$out = true;
		}
		else if ($x >= $nx)
		{
			$x = $nx - 1;
			$out = true;
		}

		return ($x);
	}
	
	/**
	  *
	  * Periodic boundary condition test
	  *
	**/
	public function periodic_bc($x, $nx, &$out)
	{
		if($x < 0)
		{
			$n   = 1 - intval($x / ($nx + 1));
			$ixx = $x + $n * $nx;

			$x = $ixx % $nx;
			$out = true;
		}
		else if($x >= $nx)
		{
			$x = $x % $nx;
			$out = true;
		}

		return ($x);
	}
	
	/**
	  *
	  * Symmetric boundary condition test
	  *
	**/
	public function symmetric_bc($x, $nx, &$out)
	{
		if($x < 0)
		{
			$borde = $nx - 1;
			$xx = -$x;
			$n  = intval($xx / $borde) % 2;

			if ($n) $x = $borde - ($xx % $borde );
			else $x = $xx % $borde;
			$out = true;
		}
		else if ($x >= $nx)
		{
			$borde = $nx - 1;
			$n = intval($x / $borde) % 2;

			if ($n) $x = $borde - ($x % $borde);
			else $x = $x % $borde;
			$out = true;
		}

		return ($x);
	}
	
	/**
	  *
	  * Cubic interpolation in one dimension
	  *
	**/
	public function cubic_interpolation_cell(
		&$v,  	//interpolation points
		$x      //point to be interpolated
	)
	{
		return  ($v[1] + 0.5 * $x * ($v[2] - $v[0] + $x * (2.0 *  $v[0] - 5.0 * $v[1] + 4.0 * $v[2] - $v[3] + $x * (3.0 * ($v[1] - $v[2]) + $v[3] - $v[0]))));
	}
	
	/**
	  *
	  * Bicubic interpolation in two dimensions
	  *
	**/
	public function bicubic_interpolation_cell(
		&$p, //array containing the interpolation points
		$x,       //x position to be interpolated
		$y        //y position to be interpolated
	)
	{
		$v = new SplFixedArray(4);
		$v[0] = $this->cubic_interpolation_cell($p[0], $y);
		$v[1] = $this->cubic_interpolation_cell($p[1], $y);
		$v[2] = $this->cubic_interpolation_cell($p[2], $y);
		$v[3] = $this->cubic_interpolation_cell($p[3], $y);
		return ($this->cubic_interpolation_cell($v, $x));
		unset($v);
	}
	
	/**
	  *
	  * Compute the bicubic interpolation of a point in an image.
	  * Detect if the point goes outside the image domain.
	  *
	**/
	public function bicubic_interpolation_at(
		&$input, //image to be interpolated
		$uu,    //x component of the vector field
		$vv,    //y component of the vector field
		$nx,    //image width
		$ny,    //image height
		$border_out //if true, return zero outside the region
	)
	{
		$sx = ($uu < 0)? -1: 1;
		$sy = ($vv < 0)? -1: 1;

		//int x, y, mx, my, dx, dy, ddx, ddy;
		$x; $y; $mx; $my; $dx; $dy; $ddx; $ddy;
		$out = false;

		//apply the corresponding boundary conditions
		switch(BOUNDARY_CONDITION) {

			case 0: $x   = $this->neumann_bc(intval($uu), $nx, $out);
				$y   = $this->neumann_bc(intval($vv), $ny, $out);
				$mx  = $this->neumann_bc(intval($uu) - $sx, $nx, $out);
				$my  = $this->neumann_bc(intval($vv) - $sx, $ny, $out); // CHCK sx??? i think should be sy
				$dx  = $this->neumann_bc(intval($uu) + $sx, $nx, $out);
				$dy  = $this->neumann_bc(intval($vv) + $sy, $ny, $out);
				$ddx = $this->neumann_bc(intval($uu) + 2 * $sx, $nx, $out);
				$ddy = $this->neumann_bc(intval($vv) + 2 * $sy, $ny, $out);
				break;

			case 1: $x   = $this->periodic_bc(intval($uu), $nx, $out);
				$y   = $this->periodic_bc(intval($vv), $ny, $out);
				$mx  = $this->periodic_bc(intval($uu) - $sx, $nx, $out);
				$my  = $this->periodic_bc(intval($vv) - $sx, $ny, $out); // CHCK same as above
				$dx  = $this->periodic_bc(intval($uu) + $sx, $nx, $out);
				$dy  = $this->periodic_bc(intval($vv) + $sy, $ny, $out);
				$ddx = $this->periodic_bc(intval($uu) + 2 * $sx, $nx, $out);
				$ddy = $this->periodic_bc(intval($vv) + 2 * $sy, $ny, $out);
				break;

			case 2: $x   = $this->symmetric_bc(intval($uu), $nx, $out);
				$y   = $this->symmetric_bc(intval($vv), $ny, $out);
				$mx  = $this->symmetric_bc(intval($uu) - $sx, $nx, $out);
				$my  = $this->symmetric_bc(intval($vv) - $sx, $ny, $out); // CHCK same as above
				$dx  = $this->symmetric_bc(intval($uu) + $sx, $nx, $out);
				$dy  = $this->symmetric_bc(intval($vv) + $sy, $ny, $out);
				$ddx = $this->symmetric_bc(intval($uu) + 2 * $sx, $nx, $out);
				$ddy = $this->symmetric_bc(intval($vv) + 2 * $sy, $ny, $out);
				break;

			default:$x   = $this->neumann_bc(intval($uu), $nx, $out);
				$y   = $this->neumann_bc(intval($vv), $ny, $out);
				$mx  = $this->neumann_bc(intval($uu) - $sx, $nx, $out);
				$my  = $this->neumann_bc(intval($vv) - $sx, $ny, $out); // CHCK same as above
				$dx  = $this->neumann_bc(intval($uu) + $sx, $nx, $out);
				$dy  = $this->neumann_bc(intval($vv) + $sy, $ny, $out);
				$ddx = $this->neumann_bc(intval($uu) + 2 * $sx, $nx, $out);
				$ddy = $this->neumann_bc(intval($vv) + 2 * $sy, $ny, $out);
				break;
		}

		if($out && $border_out)
			return (0.0);

		else
		{
			//obtain the interpolation points of the image
			$p11 = $input[$mx  + $nx * $my];
			$p12 = $input[$x   + $nx * $my];
			$p13 = $input[$dx  + $nx * $my];
			$p14 = $input[$ddx + $nx * $my];

			$p21 = $input[$mx  + $nx * $y];
			$p22 = $input[$x   + $nx * $y];
			$p23 = $input[$dx  + $nx * $y];
			$p24 = $input[$ddx + $nx * $y];

			$p31 = $input[$mx  + $nx * $dy];
			$p32 = $input[$x   + $nx * $dy];
			$p33 = $input[$dx  + $nx * $dy];
			$p34 = $input[$ddx + $nx * $dy];

			$p41 = $input[$mx  + $nx * $ddy];
			$p42 = $input[$x   + $nx * $ddy];
			$p43 = $input[$dx  + $nx * $ddy];
			$p44 = $input[$ddx + $nx * $ddy];

			//create array
			$pol = array(
				array($p11, $p21, $p31, $p41),
				array($p12, $p22, $p32, $p42),
				array($p13, $p23, $p33, $p43),
				array($p14, $p24, $p34, $p44)
				);

			//return interpolation
			return ($this->bicubic_interpolation_cell($pol, $uu - $x, $vv - $y));
			unset($pol);
		}
	}
	
	/**
	  *
	  * Compute the bicubic interpolation of an image.
	  *
	**/
	public function bicubic_interpolation_warp(
		&$input,     // image to be warped
		&$u,         // x component of the vector field
		&$v,         // y component of the vector field
		&$output,    // image warped with bicubic interpolation
		$nx,        // image width
		$ny,        // image height
		$border_out // if true, put zeros outside the region
	)
	{
		for($i = 0; $i < $ny; $i++)
			for($j = 0; $j < $nx; $j++)
			{
				$p  = $i * $nx + $j;
				$uu = $j + $u[$p];
				$vv = $i + $v[$p];

				// obtain the bicubic interpolation at position (uu, vv)
				$output[$p] = $this->bicubic_interpolation_at($input, $uu, $vv, $nx, $ny, $border_out);
			}
	}
	
	/**
	  *
	  * Compute the size of a zoomed image from the zoom factor
	  *
	**/
	public function zoom_size(
		$nx,      // width of the orignal image
		$ny,      // height of the orignal image
		&$nxx,    // width of the zoomed image
		&$nyy,    // height of the zoomed image
		$factor // zoom factor between 0 and 1
	)
	{
		//compute the new size corresponding to factor
		//we add 0.5 for rounding off to the closest number
		$nxx = intval($nx * $factor + 0.5);
		$nyy = intval($ny * $factor + 0.5);
	}
	
	/**
	  *
	  * Downsample an image
	  *
	**/
	public function zoom_out(
		&$I,    // input image
		&$Iout,       // output image
		$nx,      // image width
		$ny,      // image height
		$factor // zoom factor between 0 and 1
	)
	{
		// temporary working image
		$Is = new SplFixedArray($nx * $ny);
		for($i = 0; $i < $nx * $ny; $i++)
			$Is[$i] = $I[$i];

		// compute the size of the zoomed image
		$nxx; $nyy;
		$this->zoom_size($nx, $ny, $nxx, $nyy, $factor);

		// compute the Gaussian sigma for smoothing
		$sigma = ZOOM_SIGMA_ZERO * sqrt(1.0 / ($factor * $factor) - 1.0);

		// pre-smooth the image
		$this->gaussian($Is, $nx, $ny, $sigma);

		// re-sample the image using bicubic interpolation
		for ($i1 = 0; $i1 < $nyy; $i1++)
		{
			for ($j1 = 0; $j1 < $nxx; $j1++)
				{
					$i2  = intval($i1 / $factor);
					$j2  = intval($j1 / $factor);

					$g = $this->bicubic_interpolation_at($Is, $j2, $i2, $nx, $ny, false);
					$Iout[$i1 * $nxx + $j1] = $g;
				}
		}

		unset($Is);
	}
	
	/**
	  *
	  * Function to upsample the image
	  *
	**/
	public function zoom_in(
		&$I, // input image
		&$Iout,    // output image
		$nx,         // width of the original image
		$ny,         // height of the original image
		$nxx,        // width of the zoomed image
		$nyy         // height of the zoomed image
	)
	{
		// compute the zoom factor
		$factorx = $nxx / $nx;
		$factory = $nyy / $ny;

		// re-sample the image using bicubic interpolation
		for ($i1 = 0; $i1 < $nyy; $i1++)
		for ($j1 = 0; $j1 < $nxx; $j1++)
		{
			$i2 = $i1 / $factory;
			$j2 = $j1 / $factorx;

			$g = $this->bicubic_interpolation_at($I, $j2, $i2, $nx, $ny, false);
			$Iout[$i1 * $nxx + $j1] = $g;
		}
	}
	
	/**
	 *
	 * Details on how to compute the divergence and the grad(u) can be found in:
	 * [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
	 * Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
	 *
	 **/


	/**
	 *
	 * Function to compute the divergence with backward differences
	 * (see [2] for details)
	 *
	 **/
	public function divergence(
			&$v1, // x component of the vector field
			&$v2, // y component of the vector field
			&$div,      // output divergence
			$nx,    // image width
			$ny     // image height
		)
	{
		// compute the divergence on the central body of the image
		for ($i = 1; $i < $ny - 1; $i++)
		{
			for($j = 1; $j < $nx - 1; $j++)
			{
				$p  = $i * $nx + $j;
				$p1 = $p - 1;
				$p2 = $p - $nx;

				$v1x = $v1[$p] - $v1[$p1];
				$v2y = $v2[$p] - $v2[$p2];

				$div[$p] = $v1x + $v2y;
			}
		}

		// compute the divergence on the first and last rows
		for ($j = 1; $j < $nx - 1; $j++)
		{
			$p = ($ny - 1) * $nx + $j;

			$div[$j] = $v1[$j] - $v1[$j - 1] + $v2[$j];
			$div[$p] = $v1[$p] - $v1[$p - 1] - $v2[$p - $nx];
		}

		// compute the divergence on the first and last columns
		for ($i = 1; $i < $ny - 1; $i++)
		{
			$p1 = $i * $nx;
			$p2 = ($i + 1) * $nx - 1;

			$div[$p1] =  $v1[$p1]   + $v2[$p1] - $v2[$p1 - $nx];
			$div[$p2] = -$v1[$p2 - 1] + $v2[$p2] - $v2[$p2 - $nx];

		}

		$div[0] =  $v1[0] + $v2[0];
		$div[$nx - 1] = -$v1[$nx - 2] + $v2[$nx - 1];
		$div[($ny - 1) * $nx] = $v1[($ny - 1) * $nx] - $v2[($ny - 2) * $nx];
		$div[$ny * $nx - 1] = -$v1[$ny * $nx - 2] - $v2[($ny - 1) * $nx - 1];
	}
	
	/**
	 *
	 * Function to compute the gradient with forward differences
	 * (see [2] for details)
	 *
	 **/
	public function forward_gradient(
			&$f, //input image
			&$fx,      //computed x derivative
			&$fy,      //computed y derivative
			$nx,   //image width
			$ny    //image height
		)
	{
		// compute the gradient on the central body of the image
		for ($i = 0; $i < $ny - 1; $i++)
		{
			for($j = 0; $j < $nx - 1; $j++)
			{
				$p  = $i * $nx + $j;
				$p1 = $p + 1;
				$p2 = $p + $nx;

				$fx[$p] = $f[$p1] - $f[$p];
				$fy[$p] = $f[$p2] - $f[$p];
			}
		}

		// compute the gradient on the last row
		for ($j = 0; $j < $nx - 1; $j++)
		{
			$p = ($ny - 1) * $nx + $j;

			$fx[$p] = $f[$p + 1] - $f[$p];
			$fy[$p] = 0;
		}

		// compute the gradient on the last column
		for ($i = 1; $i < $ny; $i++)
		{
			$p = $i * $nx - 1;

			$fx[$p] = 0;
			$fy[$p] = $f[$p + $nx] - $f[$p];
		}

		$fx[$ny * $nx - 1] = 0;
		$fy[$ny * $nx - 1] = 0;
	}
	
	/**
	 *
	 * Function to compute the gradient with centered differences
	 *
	 **/
	public function centered_gradient(
			&$input,  //input image
			&$dx,           //computed x derivative
			&$dy,           //computed y derivative
			$nx,        //image width
			$ny         //image height
		)
	{
		// compute the gradient on the center body of the image
		for ($i = 1; $i < $ny - 1; $i++)
		{
			for($j = 1; $j < $nx - 1; $j++)
			{
				$k = $i * $nx + $j;
				$dx[$k] = 0.5 * ($input[$k + 1] - $input[$k - 1]);
				$dy[$k] = 0.5 * ($input[$k + $nx] - $input[$k - $nx]);
			}
		}

		// compute the gradient on the first and last rows
		for ($j = 1; $j < $nx - 1; $j++)
		{
			$dx[$j] = 0.5 * ($input[$j + 1] - $input[$j - 1]);
			$dy[$j] = 0.5 * ($input[$j + $nx] - $input[$j]);

			$k = ($ny - 1) * $nx + $j;

			$dx[$k] = 0.5 * ($input[$k + 1] - $input[$k - 1]);
			$dy[$k] = 0.5 * ($input[$k] - $input[$k - $nx]);
		}

		// compute the gradient on the first and last columns
		for($i = 1; $i < $ny - 1; $i++)
		{
			$p = $i * $nx;
			$dx[$p] = 0.5 * ($input[$p + 1] - $input[$p]);
			$dy[$p] = 0.5 * ($input[$p + $nx] - $input[$p - $nx]);

			$k = ($i + 1) * $nx - 1;

			$dx[$k] = 0.5 * ($input[$k] - $input[$k - 1]);
			$dy[$k] = 0.5 * ($input[$k + $nx] - $input[$k - $nx]);
		}

		// compute the gradient at the four corners
		$dx[0] = 0.5 * ($input[1] - $input[0]);
		$dy[0] = 0.5 * ($input[$nx] - $input[0]);

		$dx[$nx - 1] = 0.5 * ($input[$nx - 1] - $input[$nx - 2]);
		$dy[$nx - 1] = 0.5 * ($input[2 * $nx - 1] - $input[$nx - 1]);

		$dx[($ny - 1) * $nx] = 0.5 * ($input[($ny - 1) * $nx + 1] - $input[($ny - 1) * $nx]);
		$dy[($ny - 1) * $nx] = 0.5 * ($input[($ny - 1) * $nx] - $input[($ny - 2) * $nx]);

		$dx[$ny * $nx - 1] = 0.5 * ($input[$ny* $nx - 1] - $input[$ny * $nx - 1 - 1]);
		$dy[$ny * $nx - 1] = 0.5 * ($input[$ny * $nx - 1] - $input[($ny - 1) * $nx - 1]);
	}
	
	/**
	 *
	 * In-place Gaussian smoothing of an image
	 *
	 */
	public function gaussian(
		&$I,             // input/output image
		$xdim,       // image width
		$ydim,       // image height
		$sigma    // Gaussian sigma
	)
	{
		$boundary_condition = DEFAULT_BOUNDARY_CONDITION;
		$window_size = DEFAULT_GAUSSIAN_WINDOW_SIZE;

		$den  = 2 * $sigma * $sigma;
		$size = intval($window_size * $sigma) + 1;
		$bdx  = $xdim + $size;
		$bdy  = $ydim + $size;

		if ($boundary_condition && $size > $xdim) {
			//fprintf(stderr, "GaussianSmooth: sigma too large\n");
			//abort();
			//break; // some error
		}

		// compute the coefficients of the 1D convolution kernel
		$B = new SplFixedArray($size);
		for($i = 0; $i < $size; $i++)
		{
			$B[$i] = 1 / ($sigma * sqrt(2.0 * 3.1415926)) * exp(-$i * $i / $den);
		}

		// normalize the 1D convolution kernel
		$norm = 0;
		for($i = 0; $i < $size; $i++)
		{
			$norm += $B[$i];
		}
		
		$norm *= 2;
		$norm -= $B[0];
		for($i = 0; $i < $size; $i++)
		{
			$B[$i] /= $norm;
		}

		// convolution of each line of the input image
		$R = new SplFixedArray($size + $xdim + $size);

		for ($k = 0; $k < $ydim; $k++)
		{
			for ($i = $size; $i < $bdx; $i++)
			{
				$R[$i] = $I[$k * $xdim + $i - $size];
			}

			switch ($boundary_condition)
			{
			case BOUNDARY_CONDITION_DIRICHLET:
				for($i = 0, $j = $bdx; $i < $size; $i++, $j++)
				{
					$R[$i] = $R[$j] = 0;
				}
				break;

			case BOUNDARY_CONDITION_REFLECTING:
				for($i = 0, $j = $bdx; $i < $size; $i++, $j++)
				{
					$R[$i] = $I[$k * $xdim + $size - $i];
					$R[$j] = $I[$k * $xdim + $xdim - $i - 1];
				}
				break;

			case BOUNDARY_CONDITION_PERIODIC:
				for($i = 0, $j = $bdx; $i < $size; $i++, $j++)
				{
					$R[$i] = $I[$k * $xdim + $xdim - $size + $i];
					$R[$j] = $I[$k * $xdim + $i];
				}
				break;
			}

			for ($i = $size; $i < $bdx; $i++)
			{
				$sum = $B[0] * $R[$i];
				for ($j = 1; $j < $size; $j++ )
				{
					$sum += $B[$j] * ($R[$i - $j] + $R[$i + $j]);
				}
				$I[$k * $xdim + $i - $size] = $sum;
			}
		}

		// convolution of each column of the input image
		$T = new SplFixedArray($size + $ydim + $size);

		for ($k = 0; $k < $xdim; $k++)
		{
			for ($i = $size; $i < $bdy; $i++)
			{
				$T[$i] = $I[($i - $size) * $xdim + $k];
			}

			switch ($boundary_condition)
			{
			case BOUNDARY_CONDITION_DIRICHLET:
				for ($i = 0, $j = $bdy; $i < $size; $i++, $j++)
				{
					$T[$i] = $T[$j] = 0;
				}
				break;

			case BOUNDARY_CONDITION_REFLECTING:
				for ($i = 0, $j = $bdy; $i < $size; $i++, $j++)
				{
					$T[$i] = $I[($size - $i) * $xdim + $k];
					$T[$j] = $I[($ydim - $i - 1) * $xdim + $k];
				}
				break;

			case BOUNDARY_CONDITION_PERIODIC:
				for($i = 0, $j = $bdx; $i < $size; $i++, $j++)
				{
					$T[$i] = $I[($ydim - $size + $i) * $xdim + $k];
					$T[$j] = $I[$i * $xdim + $k];
				}
				break;
			}

			for ($i = $size; $i < $bdy; $i++)
			{
				$sum = $B[0] * $T[$i];
				for ($j = 1; $j < $size; $j++)
				{
					$sum += $B[$j] * ($T[$i - $j] + $T[$i + $j]);
				}
				$I[($i - $size) * $xdim + $k] = $sum;
			}
		}

		unset($R);
		unset($T);
		unset($B);
	}
	
	/**
	 * Implementation of the Zach, Pock and Bischof dual TV-L1 optic flow method
	 *
	 * see reference:
	 *  [1] C. Zach, T. Pock and H. Bischof, "A Duality Based Approach for Realtime
	 *      TV-L1 Optical Flow", In Proceedings of Pattern Recognition (DAGM),
	 *      Heidelberg, Germany, pp. 214-223, 2007
	 *
	 *
	 * Details on the total variation minimization scheme can be found in:
	 *  [2] A. Chambolle, "An Algorithm for Total Variation Minimization and
	 *      Applications", Journal of Mathematical Imaging and Vision, 20: 89-97, 2004
	 **/


	/**
	 *
	 * Function to compute the optical flow in one scale
	 *
	 **/
	public function Dual_TVL1_optic_flow(
			&$I0,           // source image
			&$I1,           // target image
			&$u1,           // x component of the optical flow
			&$u2,           // y component of the optical flow
			$nx,      // image width
			$ny,      // image height
			$tau,     // time step
			$lambda,  // weight parameter for the data term
			$theta,   // weight parameter for (u - v)²
			$warps,   // number of warpings per scale
			$epsilon, // tolerance for numerical convergence
			$verbose  // enable/disable the verbose mode
		)
	{
		$size = $nx * $ny;
		$l_t = $lambda * $theta;

		//size_t sf = sizeof(float);
		$I1x    = new SplFixedArray($size);
		$I1y    = new SplFixedArray($size);
		//$I1w    = new SplFixedArray($size);
		$I1wx   = new SplFixedArray($size);
		$I1wy   = new SplFixedArray($size);
		$rho_c  = new SplFixedArray($size);
		//$v1     = new SplFixedArray($size);
		//$v2     = new SplFixedArray($size);
		$p11    = new SplFixedArray($size);
		$p12    = new SplFixedArray($size);
		$p21    = new SplFixedArray($size);
		$p22    = new SplFixedArray($size);
		//$div    = new SplFixedArray($size);
		$grad   = new SplFixedArray($size);
		//$div_p1 = new SplFixedArray($size);
		//$div_p2 = new SplFixedArray($size);
		//$u1x    = new SplFixedArray($size);
		//$u1y    = new SplFixedArray($size);
		//$u2x    = new SplFixedArray($size);
		//$u2y    = new SplFixedArray($size);

		$this->centered_gradient($I1, $I1x, $I1y, $nx, $ny);

		// initialization of p
		for ($i = 0; $i < $size; $i++)
		{
			$p11[$i] = $p12[$i] = 0.0;
			$p21[$i] = $p22[$i] = 0.0;
		}

		for ($warpings = 0; $warpings < $warps; $warpings++)
		{
			$I1w    = new SplFixedArray($size);
			// compute the warping of the target image and its derivatives
			$this->bicubic_interpolation_warp($I1, $u1, $u2, $I1w,  $nx, $ny, true);
			$this->bicubic_interpolation_warp($I1x, $u1, $u2, $I1wx, $nx, $ny, true);
			$this->bicubic_interpolation_warp($I1y, $u1, $u2, $I1wy, $nx, $ny, true);

			for ($i = 0; $i < $size; $i++)
			{
				$Ix2 = $I1wx[$i] * $I1wx[$i];
				$Iy2 = $I1wy[$i] * $I1wy[$i];

				// store the |Grad(I1)|^2
				$grad[$i] = ($Ix2 + $Iy2);

				// compute the constant part of the rho function
				$rho_c[$i] = ($I1w[$i] - $I1wx[$i] * $u1[$i] - $I1wy[$i] * $u2[$i] - $I0[$i]);
			}
			
			unset($I1w);

			$n = 0;
			$error = INF;
			while ($error > $epsilon * $epsilon && $n < MAX_ITERATIONS)
			{
				$n++;
				// estimate the values of the variable (v1, v2)
				// (thresholding opterator TH)
				$v1     = new SplFixedArray($size);
				$v2     = new SplFixedArray($size);
				for ($i = 0; $i < $size; $i++)
				{
					$rho = $rho_c[$i] + ($I1wx[$i] * $u1[$i] + $I1wy[$i] * $u2[$i]);

					$d1; $d2;

					if ($rho < -$l_t * $grad[$i])
					{
						$d1 = $l_t * $I1wx[$i];
						$d2 = $l_t * $I1wy[$i];
					}
					else
					{
						if ($rho > $l_t * $grad[$i])
						{
							$d1 = -$l_t * $I1wx[$i];
							$d2 = -$l_t * $I1wy[$i];
						}
						else
						{
							if ($grad[$i] < GRAD_IS_ZERO)
								$d1 = $d2 = 0;
							else
							{
								$fi = -$rho / $grad[$i];
								$d1 = $fi * $I1wx[$i];
								$d2 = $fi * $I1wy[$i];
							}
						}
					}

					$v1[$i] = $u1[$i] + $d1;
					$v2[$i] = $u2[$i] + $d2;
				}
				
				$div_p1 = new SplFixedArray($size);
				$div_p2 = new SplFixedArray($size);

				// compute the divergence of the dual variable (p1, p2)
				$this->divergence($p11, $p12, $div_p1, $nx ,$ny);
				$this->divergence($p21, $p22, $div_p2, $nx ,$ny);

				// estimate the values of the optical flow (u1, u2)
				$error = 0.0;
				for ($i = 0; $i < $size; $i++)
				{
					$u1k = $u1[$i];
					$u2k = $u2[$i];

					$u1[$i] = $v1[$i] + $theta * $div_p1[$i];
					$u2[$i] = $v2[$i] + $theta * $div_p2[$i];

					$error += ($u1[$i] - $u1k) * ($u1[$i] - $u1k) + ($u2[$i] - $u2k) * ($u2[$i] - $u2k);
				}
				$error /= $size;
				unset($v1);
				unset($v2);
				unset($div_p1);
				unset($div_p2);
				
				$u1x    = new SplFixedArray($size);
				$u1y    = new SplFixedArray($size);
				$u2x    = new SplFixedArray($size);
				$u2y    = new SplFixedArray($size);

				// compute the gradient of the optical flow (Du1, Du2)
				$this->forward_gradient($u1, $u1x, $u1y, $nx , $ny);
				$this->forward_gradient($u2, $u2x, $u2y, $nx , $ny);

				// estimate the values of the dual variable (p1, p2)
				for ($i = 0; $i < $size; $i++)
				{
					$taut = $tau / $theta;
					$g1   = hypot($u1x[$i], $u1y[$i]);
					$g2   = hypot($u2x[$i], $u2y[$i]);
					$ng1  = 1.0 + $taut * $g1;
					$ng2  = 1.0 + $taut * $g2;

					$p11[$i] = ($p11[$i] + $taut * $u1x[$i]) / $ng1;
					$p12[$i] = ($p12[$i] + $taut * $u1y[$i]) / $ng1;
					$p21[$i] = ($p21[$i] + $taut * $u2x[$i]) / $ng2;
					$p22[$i] = ($p22[$i] + $taut * $u2y[$i]) / $ng2;
				}
				
				unset($u1x);
				unset($u1y);
				unset($u2x);
				unset($u2y);
			}
		}

		// delete allocated memory
		unset($I1x);
		unset($I1y);
		//unset($I1w);
		unset($I1wx);
		unset($I1wy);
		unset($rho_c);
		//unset($v1);
		//unset($v2);
		unset($p11);
		unset($p12);
		unset($p21);
		unset($p22);
		//unset($div);
		unset($grad);
		//unset($div_p1);
		//unset($div_p2);
		//unset($u1x);
		//unset($u1y);
		//unset($u2x);
		//unset($u2y);
	}
	
	/**
	 *
	 * Compute the max and min of an array
	 *
	 **/
	public function getminmax(
		&$min,     // output min
		&$max,     // output max
		&$x, // input array
		$n           // array size
	)
	{
		$min = $max = $x[0];
		for ($i = 1; $i < $n; $i++) {
			if ($x[$i] < $min)
				$min = $x[$i];
			if ($x[$i] > $max)
				$max = $x[$i];
		}
	}
	
	/**
	 *
	 * Function to normalize the images between 0 and 255
	 *
	 **/
	public function image_normalization(
			&$I0,  // input image0
			&$I1,  // input image1
			&$I0n,       // normalized output image0
			&$I1n,       // normalized output image1
			$size          // size of the image
		)
	{
		$max0; $max1; $min0; $min1;

		// obtain the max and min of each image
		$this->getminmax($min0, $max0, $I0, $size);
		$this->getminmax($min1, $max1, $I1, $size);

		// obtain the max and min of both images
		$max = ($max0 > $max1)? $max0 : $max1;
		$min = ($min0 < $min1)? $min0 : $min1;
		$den = $max - $min;

		if ($den > 0)
			// normalize both images
			for ($i = 0; $i < $size; $i++)
			{
				$I0n[$i] = 255.0 * ($I0[$i] - $min) / $den;
				$I1n[$i] = 255.0 * ($I1[$i] - $min) / $den;
			}

		else
			// copy the original images
			for ($i = 0; $i < $size; $i++)
			{
				$I0n[$i] = $I0[$i];
				$I1n[$i] = $I1[$i];
			}
	}
	
	/**
	 *
	 * Function to compute the optical flow using multiple scales
	 *
	 **/
	public function Dual_TVL1_optic_flow_multiscale(
			&$I0,           // source image
			&$I1,           // target image
			&$u1,           // x component of the optical flow
			&$u2,           // y component of the optical flow
			$nxx,     // image width
			$nyy,     // image height
			$tau,     // time step
			$lambda,  // weight parameter for the data term
			$theta,   // weight parameter for (u - v)²
			$nscales, // number of scales
			$zfactor, // factor for building the image piramid
			$warps,   // number of warpings per scale
			$epsilon, // tolerance for numerical convergence
			$verbose,  // enable/disable the verbose mode
			$scl
		)
	{
		$size = $nxx * $nyy;

		// allocate memory for the pyramid structure
		$I0s = new SplFixedArray($nscales);
		$I1s = new SplFixedArray($nscales);
		$u1s = new SplFixedArray($nscales);
		$u2s = new SplFixedArray($nscales);
		$nx  = new SplFixedArray($nscales);
		$ny  = new SplFixedArray($nscales);

		$I0s[0] = new SplFixedArray($size);
		$I1s[0] = new SplFixedArray($size);

		$u1s[0] = $u1;
		$u2s[0] = $u2;
		$nx [0] = $nxx;
		$ny[0] = $nyy;

		// normalize the images between 0 and 255
		$this->image_normalization($I0, $I1, $I0s[0], $I1s[0], $size);

		// pre-smooth the original images
		$this->gaussian($I0s[0], $nx[0], $ny[0], PRESMOOTHING_SIGMA);
		$this->gaussian($I1s[0], $nx[0], $ny[0], PRESMOOTHING_SIGMA);

		// create the scales
		for ($s = 1; $s < $nscales; $s++)
		{
			$temp1; $temp2;
			$this->zoom_size($nx[$s - 1], $ny[$s - 1], /*$nx[$s]*/$temp1, /*$ny[$s]*/ $temp2, $zfactor);
			$nx[$s] = $temp1; $ny[$s] = $temp2;
			$sizes = $nx[$s] * $ny[$s];

			// allocate memory
			$I0s[$s] = new SplFixedArray($sizes);
			$I1s[$s] = new SplFixedArray($sizes);
			$u1s[$s] = new SplFixedArray($sizes);
			$u2s[$s] = new SplFixedArray($sizes);

			// zoom in the images to create the pyramidal structure
			$this->zoom_out($I0s[$s - 1], $I0s[$s], $nx[$s - 1], $ny[$s - 1], $zfactor);
			$this->zoom_out($I1s[$s - 1], $I1s[$s], $nx[$s - 1], $ny[$s - 1], $zfactor);
		}

		// initialize the flow at the coarsest scale
		for ($i = 0; $i < $nx[$nscales - 1] * $ny[$nscales - 1]; $i++)
		{
			$u1s[$nscales - 1][$i] = $u2s[$nscales - 1][$i] = 0.0;
		}

		// pyramidal structure for computing the optical flow
		for ($s = $nscales - 1; $s >= 1; $s--)
		{
			//echo(' iteracyja: ');
			//var_dump($s);
			//echo(' NX: ');
			//var_dump($nx[$s]);
			//echo(' NY: ');
			//var_dump($ny[$s]);
			// compute the optical flow at the current scale
			$this->Dual_TVL1_optic_flow($I0s[$s], $I1s[$s], $u1s[$s], $u2s[$s], $nx[$s], $ny[$s], $tau, $lambda, $theta, $warps, $epsilon, $verbose);

			// if this was the last scale, finish now
			if ($s==$scl){
				//echo ('NXX = '); var_dump($nxx);
				//echo ('NYY = '); var_dump($nyy);
				//echo ('NX = '); var_dump($nx[$s]);
				//echo ('NY = '); var_dump($ny[$s]);
				for ($j = 0; $j < $nyy; $j++)
				{
					for ($i = 0; $i < $nxx; $i++)
					{
						//$ii = intval($i/2);
						//$jj = intval($j/2);
						//echo (" X= "); var_dump($ii);
						//echo (" Y= "); var_dump($jj);
						//$u1[$j * $nxx + $i] = 0;
						//$u2[$j * $nxx + $i] = 0;
						$u1[$j * $nxx + $i] = $u1s[$s][intval($j / pow(2, $scl)) * $nx[$s] + intval($i / pow(2, $scl))] / $zfactor;
						$u2[$j * $nxx + $i] = $u2s[$s][intval($j / pow(2, $scl)) * $nx[$s] + intval($i / pow(2, $scl))] / $zfactor;
					}
				}
				break;
			}

			// otherwise, upsample the optical flow

			// zoom the optical flow for the next finer scale
			$this->zoom_in($u1s[$s], $u1s[$s - 1], $nx[$s], $ny[$s], $nx[$s - 1], $ny[$s - 1]);
			$this->zoom_in($u2s[$s], $u2s[$s - 1], $nx[$s], $ny[$s], $nx[$s - 1], $ny[$s - 1]);
			
			unset($I0s[$s]);
			unset($I1s[$s]);
			unset($u1s[$s]);
			unset($u2s[$s]);

			// scale the optical flow with the appropriate zoom factor
			for ($i = 0; $i < $nx[$s - 1] * $ny[$s - 1]; $i++)
			{
				$u1s[$s - 1][$i] *= 1.0 / $zfactor;
				$u2s[$s - 1][$i] *= 1.0 / $zfactor;
			}
		}

		// delete allocated memory
		/*for ($i = 1; $i < $nscales; $i++)
		{
			unset($I0s[$i]);
			unset($I1s[$i]);
			unset($u1s[$i]);
			unset($u2s[$i]);
		}*/
		unset($I0s[0]);
		unset($I1s[0]);

		unset($I0s);
		unset($I1s);
		unset($u1s);
		unset($u2s);
		unset($nx);
		unset($ny);
	}

}

?>