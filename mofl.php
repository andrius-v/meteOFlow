<?
	include_once('TVL1.php');
	gc_enable();
	
	function LoadPNG($imgname)
	{
		/* Attempt to open */
		$im = @imagecreatefrompng($imgname);

		/* See if it failed */
		if(!$im)
		{
			/* Create a blank image */
			$im  = imagecreatetruecolor(150, 30);
			$bgc = imagecolorallocate($im, 255, 255, 255);
			$tc  = imagecolorallocate($im, 0, 0, 0);

			imagefilledrectangle($im, 0, 0, 150, 30, $bgc);

			/* Output an error message */
			imagestring($im, 1, 5, 5, 'Error loading ' . $imgname, $tc);
		}

		return $im;
	}
	
	// Lets check if data is updated
	$img0 = LoadPNG('check.png');
	$img1 = LoadPNG('http://meteo.lt/dokumentai/operatyvi_inf/radaras/radarlarge+36.png');
	$imgT = imagecreatetruecolor(imagesx($img0), imagesy($img0));
	
	imagecopy($imgT, $img1, 0, 0, 0, imagesy($img1) - imagesy($img0), imagesx($img0),  imagesy($img0));
	
	$diff = false;
	for ($y = 0; $y < imagesy($imgT); $y++)
	{
		for ($x = 0; $x < imagesx($imgT); $x++)
		{
			$c0 = imagecolorat($img0, $x, $y);
			$cT = imagecolorat($imgT, $x, $y);
			if (($c0 - $cT) != 0)
			{
				$diff = true;
				break;
			}
		}
	}
	
	//imagedestroy($img0);
	//imagedestroy($imgT);
	
	if ($diff)
	{
		$img0 = LoadPNG('http://meteo.lt/dokumentai/operatyvi_inf/radaras/radarlarge+35.png');
		$img2 = imagecreatetruecolor(imagesx($imgT), imagesy($imgT));
		imagecopy($img2, $img1, 0, 0, 0, imagesy($img1) - imagesy($img2), imagesx($img2),  imagesy($img2));
		$hs = new TVL1($img0, $img1, 1);
		$hs->Y2RGB($img1, 2, 1, 255, 8);
		imagecopy($img1, $img2, 0, imagesy($img1) - imagesy($img2), 0, 0, imagesx($img2),  imagesy($img2));
		header('Content-Disposition: Attachment;filename=mOFl.png');
		header('Content-Type: image/png');
		imagepng($img1);
		imagepng($img1, "mOFl.png");
		imagepng($imgT, 'check.png');
	}
	else
	{
		$imgT = loadPNG('mOFl.png');
		header('Content-Disposition: Attachment;filename=mOFl.png');
		header('Content-Type: image/png');
		imagepng($imgT);
	}
	
	imagedestroy($img2);
	imagedestroy($imgT);
	imagedestroy($img0);
	imagedestroy($img1);
	gc_disable();
?>