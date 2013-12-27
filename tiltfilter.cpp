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

  Matrix3x3(double a11, double a12, double a13,
            double a21, double a22, double a23,
            double a31, double a32, double a33)
  {
    m11 = a11;  m12 = a12;  m13 = a13;
    m21 = a21;  m22 = a22;  m23 = a23;
    m31 = a31;  m32 = a32;  m33 = a33;
  }

  ~Matrix3x3() {}

  void identity()
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

  Matrix3x3 operator *(const Matrix3x3& mat) const
  {
    Matrix3x3 ret;
    ret.m11 = m11*mat.m11 + m12*mat.m21 + m13*mat.m31;
    ret.m12 = m11*mat.m12 + m12*mat.m22 + m13*mat.m32;
    ret.m13 = m11*mat.m13 + m12*mat.m23 + m13*mat.m33;

    ret.m21 = m21*mat.m11 + m22*mat.m21 + m23*mat.m31;
    ret.m22 = m21*mat.m12 + m22*mat.m22 + m23*mat.m32;
    ret.m23 = m21*mat.m13 + m22*mat.m23 + m23*mat.m33;

    ret.m31 = m31*mat.m11 + m32*mat.m21 + m33*mat.m31;
    ret.m32 = m31*mat.m12 + m32*mat.m22 + m33*mat.m32;
    ret.m33 = m31*mat.m13 + m32*mat.m23 + m33*mat.m33;
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
    Matrix3x3 rotYp(0, 0, -1, 0, 1, 0,  1, 0, 0);
    Matrix3x3 rotP (1, 0, 0, 0, cosP, -sinP, 0, sinP, cosP);
    Matrix3x3 rotR (cosR, sinR, 0, -sinR, cosR, 0, 0, 0, 1);
    Matrix3x3 rotYn(0, 0,  1, 0, 1, 0, -1, 0, 0);
    _tilt = rotYp * rotR * rotP * rotYn;
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
    y = (y<0)? 0: ((y>=_height)? _height-1: y);
    size_t ofs = (x%_width)*4+y*_rowInBytes;
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
  void apply(uint8_t* dst, int ps, int ph) const
  {
    int dx = _width * _yaw / (2.0*M_PI);
    double w = _width;
    double h = _height;

#ifdef _OPENMP
    int procs = omp_get_num_procs();
    omp_set_num_threads(procs);
#pragma pragma omp parallel for
#endif
    for (int i = ps; i<ps+ph; ++i) {
      double theta0 = M_PI / 2.0 - M_PI * i / h;
      double y0 = sin(theta0);
      double cosTheta = cos(theta0);
      int o = i*_rowInBytes;
      for (int j = 0; j<_width; ++j) {
        double phi0 = M_PI * 2.0 * j / w;
        Vector3 p(cosTheta * cos(phi0), y0, cosTheta * sin(phi0));
        Vector3 q = _tilt*p;
        double theta = asin(q.y);
        double phi = atan2(q.z, q.x);
        double u = phi * w / (2.0 * M_PI);
        double v = h / 2.0 - theta * h / M_PI;
        if (u<0) u+=_width;
        if (v<0) v+=_height;
        sampler(&dst[o + ((j + dx + _width)%_width)*4], u, v);
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
  annotate("as3package:info.smoche.TiltFilterBase")
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
  annotate("as3package:info.smoche.TiltFilterBase")
));
void TiltFilter_free()
{
  TiltFilter* self = TiltFilter::as3_getSelf();
  delete self;
}


void TiltFilter_apply() __attribute__((used,
  annotate("as3sig:public function TiltFilter_apply(as3_self:int, as3_target:BitmapData, as3_procStart:uint, as3_procHeight:uint, root:Sprite):void"),
  annotate("as3package:info.smoche.TiltFilterBase"),
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
  size_t width, height, procStart, procHeight;
  AS3_GetScalarFromVar(src, as3_src);
  AS3_GetScalarFromVar(dst, as3_dst);
  AS3_GetScalarFromVar(width, as3_width);
  AS3_GetScalarFromVar(height, as3_height);
  AS3_GetScalarFromVar(procStart, as3_procStart);
  AS3_GetScalarFromVar(procHeight, as3_procHeight);
  if (procStart>=height) return;
  if (procStart+procHeight>height) {
    procHeight = height - procStart;
  }

  self->setSource(src, width, height);
  self->apply(dst, procStart, procHeight);

  inline_as3("\
    ba = new ByteArray();\
    CModule.readBytes(as3_dst, as3_width*as3_height*4, ba);\
    ba.position = 0;\
    as3_target.setPixels(as3_rect, ba);\
    CModule.free(as3_src);\
    CModule.free(as3_dst);\
  ");
}

