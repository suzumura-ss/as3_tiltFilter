#include <AS3/AS3.h>
#include <sys/types.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>


class Vector3 {
public:
  double x, y, z;
  Vector3(double x0=0, double y0=0, double z0=0)
  {
    x = x0;  y = y0;  z = z0;
  }
  ~Vector3() {}
};

class Matrix3x3 {
public:
  double m11, m12, m13;
  double m21, m22, m23;
  double m31, m32, m33;
  inline Matrix3x3()
  {
    identity();
  }
  inline ~Matrix3x3() {}
  inline void identity()
  {
    m11 = 1;  m12 = 0;  m13 = 0;
    m21 = 0;  m22 = 1;  m23 = 0;
    m31 = 0;  m32 = 0;  m33 = 1;
  }
  inline Vector3 operator *(const Vector3& vec) const
  {
    Vector3 ret;
    ret.x = m11*vec.x + m12*vec.y + m13*vec.z;
    ret.y = m21*vec.x + m22*vec.y + m23*vec.z;
    ret.z = m31*vec.x + m32*vec.y + m33*vec.z;
    return ret;
  }
};

class TiltFilter {
protected:
  Matrix3x3 _tilt;
  double _yaw;
  int _width;
  int _rowInBytes;
  int _height;
  const uint8_t* _src;

public:
  TiltFilter()
  {
    _tilt.identity();
    _src = NULL;
  }
  ~TiltFilter()
  {
  }

  void setTilt(double y, double p, double r)
  {
    _yaw = y;
    double cosP = cos(p);
    double sinP = sin(p);
    double cosR = cos(r);
    double sinR = sin(r);
    _tilt.m11 =         cosR;
    _tilt.m12 =        -sinR;
    _tilt.m13 =         0.0;
    _tilt.m21 =  cosP * sinR;
    _tilt.m22 =  cosP * cosR;
    _tilt.m23 =  sinP;
    _tilt.m31 = -sinP * sinR;
    _tilt.m32 = -sinP * cosR;
    _tilt.m33 =  cosP;
    
  }

  void setSource(const uint8_t* src, int width, int height)
  {
    _src = src;
    _width = width;
    _height = height;
    _rowInBytes = width*4;
  }

protected:
  static inline float crop(float v)
  {
    return (v<0)? 0: ((v>255)? 255: v);
  }

  inline void colorAdd(float* rgb, float weight, int x, int y) const
  {
    size_t ofs = (x%_width)*4+(y%_height)*_rowInBytes;
    rgb[1] += weight*_src[ofs+1];
    rgb[2] += weight*_src[ofs+2];
    rgb[3] += weight*_src[ofs+3];
  }

  inline void sampler(uint8_t* ret, double x, double y) const
  {
    int ix = x;
    int iy = y;
    double ax = x-ix;
    double ay = y-iy;
    float rgb[4] = {0,0,0,0};
    colorAdd(rgb, (1-ax)*(1-ay), ix  , iy  );
    colorAdd(rgb, (  ax)*(1-ay), ix+1, iy  );
    colorAdd(rgb, (1-ax)*(  ay), ix  , iy+1);
    colorAdd(rgb, (  ax)*(  ay), ix+1, iy+1);

    ret[0] = 255;
    ret[1] = crop(rgb[1]);
    ret[2] = crop(rgb[2]);
    ret[3] = crop(rgb[3]);
  }

public:
  void apply(uint8_t* dst) const
  {
#ifdef _OPENMP
    int t = omp_get_num_procs();
    omp_set_num_threads(t/t);
#endif

    double dx = _width * _yaw / (2.0*M_PI);

#pragma omp parallel for
    for (int i = 0; i<_height; ++i) {
      double theta0 = M_PI/2.0 - M_PI*i/_height;
      double y0 = sin(theta0);
      double cosTheta = cos(theta0);
      for (int j = 0; j<_width; ++j) {
        double phi0 = M_PI*2.0*j/_width;
        Vector3 p(cosTheta * cos(phi0), y0, cosTheta * sin(phi0));
        Vector3 q = _tilt*p;
        double theta = asin(q.y);
        double phi = atan2(q.z, q.x);
        double u = phi * _width / (2.0 * M_PI);
        double v = _height / 2.0 - theta * _height / M_PI;
        int o = (((int)(j+dx)+_width)%_width)*4 + i*_rowInBytes;
        sampler(&dst[o], u, v);
      }
    }
  }

public:
  void as3_return() const
  {
    AS3_Return((int)this);
  }

  static TiltFilter* as3_getSelf()
  {
    int self;
    AS3_GetScalarFromVar(self, as3_self);
    return (TiltFilter*)self;
  }
};



void TiltFilter_init() __attribute__((used,
  annotate("as3sig:public function TiltFilter_init(as3_yaw:Number, as3_pitch:Number, as3_roll:Number):int"),
  annotate("as3package:info.smoche.TiltFilter")
));
void TiltFilter_init()
{
  double yaw, pitch, roll;
  AS3_GetScalarFromVar(yaw, as3_yaw);
  AS3_GetScalarFromVar(pitch, as3_pitch);
  AS3_GetScalarFromVar(roll, as3_roll);

  TiltFilter* filter = new TiltFilter();
  filter->setTilt(yaw, pitch, roll);
  filter->as3_return();
}


void TiltFilter_free() __attribute__((used,
  annotate("as3sig:public function TiltFilter_free(as3_self:int):void"),
  annotate("as3package:info.smoche.TiltFilter")
));
void TiltFilter_free()
{
  TiltFilter* self = TiltFilter::as3_getSelf();
  delete self;
}


void TiltFilter_apply() __attribute__((used,
  annotate("as3sig:public function TiltFilter_apply(as3_self:int, as3_target:BitmapData, root:Sprite):void"),
  annotate("as3package:info.smoche.TiltFilter"),
  annotate("as3import:flash.display.BitmapData"),
  annotate("as3import:flash.display.Sprite"),
  annotate("as3import:flash.geom.Rectangle"),
  annotate("as3import:flash.utils.ByteArray")
));
void TiltFilter_apply()
{
  TiltFilter* self = TiltFilter::as3_getSelf();

  inline_as3("\
    var as3_width:uint = as3_target.width;\
    var as3_height:uint = as3_target.height;\
    var as3_rect:Rectangle = new Rectangle(0, 0, as3_width, as3_height);\
    \
    var ba:ByteArray = as3_target.getPixels(as3_rect);\
    var as3_src:int = CModule.malloc(ba.length);\
    ba.position = 0;\
    CModule.writeBytes(as3_src, ba.length, ba);\
    var as3_dst:int = CModule.malloc(ba.length);\
    \
    CModule.rootSprite = root;\
  ");
  uint8_t* src;
  uint8_t* dst;
  size_t width, height;
  AS3_GetScalarFromVar(src, as3_src);
  AS3_GetScalarFromVar(dst, as3_dst);
  AS3_GetScalarFromVar(width, as3_width);
  AS3_GetScalarFromVar(height, as3_height);

  self->setSource(src, width, height);
  self->apply(dst);

  inline_as3("\
    ba = new ByteArray();\
    CModule.readBytes(as3_dst, as3_width*as3_height*4, ba);\
    ba.position = 0;\
    as3_target.setPixels(as3_rect, ba);\
    CModule.free(as3_src);\
    CModule.free(as3_dst);\
  ");
}

