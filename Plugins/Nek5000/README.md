This is a plugin source code for the NEK5000 data format, for use with version 5.7 of ParaView

How to compile:

  Use ParaView version 5.7. Depending on the optional package you have built into your source (e.g. OSPRay, OIDN, VisRTX), you will need some extra env flags to help cmake. In my case, I explicitly set ospray_DIR, OpenImageDenoise_DIR, and VisRTX_DIR before running cmake)
   
  cd Nek5000

  mkdir build

  cd build

  cmake ../

  make
