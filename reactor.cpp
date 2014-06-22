#include <Magick++.h>
#include <Magick++/STL.h>
#include "common.h"

using Magick::InitializeMagick;
using Magick::Image;
using Magick::writeImages;
using Magick::DrawableCircle;
using Magick::DrawableText;

int main(int argc,char **argv) 
{ 
  InitializeMagick(*argv);
  vector<Image> v;
  for (int i = 0; i < 100; ++i)
  {
      v.push_back(Image("1000x1000", "white"));
      v[i].strokeColor("red"); // Outline color
      v[i].fillColor("green"); // Fill color
      v[i].strokeWidth(2);
      v[i].draw( DrawableCircle(50 + i ,50 + i, 150 + i,100 + i));
      v[i].draw( DrawableText(300, 300, "I love Russ пупкин"));
      v[i].pixelColor(49 + i/3, 49, "red");
  }
  writeImages(v.begin(), v.end(), "ttt.gif");
  //Image image( "100x100", "white" ); 
  //image.pixelColor( 49, 49, "red" ); 
  //image.write( "red_pixel.gif" ); 
  return 0; 
}
