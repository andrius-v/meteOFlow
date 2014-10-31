meteOFlow
=========

Optical flow in PHP for meteo.lt radar pictures. Requires GD.

Usage
=====

To initialize class you need two images. Optical flow calculation starts on class construction.

`$opticalFlow = new TVL1($firstImage, $secondImage, $fineTune);`

Dimensions of images must be the same. As script uses tons of memory `$fineTune` parameter tells the script when to stop.
My sugesstion start at `$fineTune = 0` (finest scale), and increase it until memory/time limits are satisfying.

To draw optical flow as lines onto image use method Y2RGB:

`$opticalFlow->Y2RGB($image, $scale, $thickness, $color, $step);`

where:
- `$image` - GD image on which we draw optical flow
- `$scale` - (float) optical flow scale factor, scales flow vectors and is used as drawing threshold
- `$thickness` - thickness of lines (not used)
- `$step` - defines flow drawing grid.

=====
Original C code and article can be found here: http://www.ipol.im/pub/art/2013/26/
