// declare some constants
#declare data = "data.csv"

// declare variables
#declare var1 = 0.0;
#declare var2 = 0.0;
#declare psi4 = 0.0;
#declare psi6 = 0.0;

// File with predefined colors
#include "colors.inc"

// Place the camera
camera {
  location <50, 50, 100>
  look_at  <50, 50, 0>
}

// create a lightsource
light_source { 
  <50, 50, 100> 
  color White
}

// Set a background color
background { color White }

// open a data file
#fopen file data read

// read as long as there is data
#while (defined(file))
  
  #read (file,var1,var2,psi4,psi6)

  sphere {
    <var1*100/sqrt(500), var2*100/sqrt(500), 0>, 0.5
    texture {
      pigment { color Gray30 }
    }
  }

  sphere {
    <var1*100/sqrt(500), var2*100/sqrt(500), 0>, 0.5
    texture {
      pigment { color rgbf <0, 0, 1, psi6> }
    }
  }

  sphere {
    <var1*100/sqrt(500), var2*100/sqrt(500), 0>, 0.5
    texture {
      pigment { color rgbf <1, 0, 0, psi4> }
    }
  }

#end