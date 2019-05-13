/** @file sift.h
 ** @brief SIFT (@ref sift)
 ** @author Andrea Vedaldi
 **/

/*
Copyright (C) 2007-12 Andrea Vedaldi and Brian Fulkerson.
All rights reserved.

This file is part of the VLFeat library and is made available under
the terms of the BSD license (see the COPYING file).
*/

#ifndef VL_SIFT_H
#define VL_SIFT_H

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
// add this to support sse2
#include <immintrin.h>

/* All additional defines are here*/
// define inline and export, see vlfeat/vl/host.h
#  define VL_INLINE static __inline__
//#  ifdef VL_BUILD_DLL
//#    ifdef __cplusplus
//#      define VL_EXPORT __attribute__((visibility ("default"))) extern "C"
//#    else
#      define VL_EXPORT __attribute__((visibility ("default"))) extern
//#    endif
//#  else
//#    ifdef __cplusplus
//#      define VL_EXPORT extern "C"
//#    else
//#      define VL_EXPORT extern
//#    endif
//#  endif

// define min, max, VL_SHIFT_LEFT, see vlfeat/vl/generic.h
#define VL_MIN(x,y) (((x)<(y))?(x):(y))
#define VL_MAX(x,y) (((x)>(y))?(x):(y))
#define VL_SHIFT_LEFT(x,n) (((n)>=0)?((x)<<(n)):((x)>>-(n)))
// define VL_PAD_BY_CONTINUITY and VL_TRANSPOSE, see vlfeat/vl/imopv.h
#define VL_PAD_BY_ZERO       (0x0 << 0) /**< @brief Pad with zeroes. */
#define VL_PAD_MASK          (0x3)      /**< @brief Padding field selector. */
#define VL_PAD_BY_CONTINUITY (0x1 << 0) /**< @brief Pad by continuity. */
#define VL_TRANSPOSE         (0x1 << 2) /**< @brief Transpose result. */
// define constant
#define VL_TRUE 1   /**< @brief @c true (1) constant */
#define VL_FALSE 0  /**< @brief @c false (0) constant */
#define VL_PI 3.141592653589793
#define VL_EPSILON_F 1.19209290E-07F
#define VL_EPSILON_D 2.220446049250313e-16

// Define error constant
#define VL_ERR_OK       0  /**< No error */
#define VL_ERR_OVERFLOW 1  /**< Buffer overflow error */
#define VL_ERR_ALLOC    2  /**< Resource allocation error */
#define VL_ERR_BAD_ARG  3  /**< Bad argument or illegal data error */
#define VL_ERR_IO       4  /**< Input/output error */
#define VL_ERR_EOF      5  /**< End-of-file or end-of-sequence error */
#define VL_LOG_OF_2 0.693147180559945

// ref: vl/float.th
// vl_imconvcol_vf_sse2 take vl_sift_pix (float) as input, which is T
// So FLT should be VL_TYPE_FLOAT, SFX should be f, then VSIZE and VTYPE is ensured
#  define VSIZE  4
#  define VSFX   s
#  define VTYPE  __m128

// define vl_uindex, see vlfeat/vl/host.h (suppose LP64, LLP64)
typedef long long unsigned  vl_uint64 ;
typedef vl_uint64           vl_uindex ;
typedef long long           vl_int64 ;
typedef vl_uint64           vl_size ;
typedef vl_int64            vl_index ;
typedef vl_uint64           vl_uindex ;
typedef vl_uint64           vl_uintptr ;
typedef int                 vl_int32 ;
typedef int                 vl_bool ;
// if ILP32, use two lines below
//typedef long long           vl_int64 ;
//typedef int       unsigned  vl_uint32 ;
//typedef vl_uint32           vl_uindex ;
//typedef vl_uint32           vl_size ;
//typedef vl_int32            vl_index ;
//typedef vl_uint32           vl_uindex ;
//typedef vl_uint32           vl_uintptr ;
//typedef int                 vl_int32 ;
//typedef int                 vl_bool ;


/** @brief SIFT filter pixel type */
typedef float vl_sift_pix ;

/** ------------------------------------------------------------------
 ** @brief SIFT filter keypoint
 **
 ** This structure represent a keypoint as extracted by the SIFT
 ** filter ::VlSiftFilt.
 **/

typedef struct _VlSiftKeypoint
{
  int o ;           /**< o coordinate (octave). */

  int ix ;          /**< Integer unnormalized x coordinate. */
  int iy ;          /**< Integer unnormalized y coordinate. */
  int is ;          /**< Integer s coordinate. */

  float x ;     /**< x coordinate. */
  float y ;     /**< y coordinate. */
  float s ;     /**< s coordinate. */
  float sigma ; /**< scale. */
} VlSiftKeypoint ;

/** ------------------------------------------------------------------
 ** @brief SIFT filter
 **
 ** This filter implements the SIFT detector and descriptor.
 **/

typedef struct _VlSiftFilt
{
  double sigman ;       /**< nominal image smoothing. */
  double sigma0 ;       /**< smoothing of pyramid base. */
  double sigmak ;       /**< k-smoothing */
  double dsigma0 ;      /**< delta-smoothing. */

  int width ;           /**< image width. */
  int height ;          /**< image height. */
  int O ;               /**< number of octaves. */
  int S ;               /**< number of levels per octave. */
  int o_min ;           /**< minimum octave index. */
  int s_min ;           /**< minimum level index. */
  int s_max ;           /**< maximum level index. */
  int o_cur ;           /**< current octave. */

  vl_sift_pix *temp ;   /**< temporary pixel buffer. */
  vl_sift_pix *octave ; /**< current GSS data. */
  vl_sift_pix *dog ;    /**< current DoG data. */
  int octave_width ;    /**< current octave width. */
  int octave_height ;   /**< current octave height. */

  vl_sift_pix *gaussFilter ;  /**< current Gaussian filter */
  double gaussFilterSigma ;   /**< current Gaussian filter std */
  vl_size gaussFilterWidth ;  /**< current Gaussian filter width */

  VlSiftKeypoint* keys ;/**< detected keypoints. */
  int nkeys ;           /**< number of detected keypoints. */
  int keys_res ;        /**< size of the keys buffer. */

  double peak_thresh ;  /**< peak threshold. */
  double edge_thresh ;  /**< edge threshold. */
  double norm_thresh ;  /**< norm threshold. */
  double magnif ;       /**< magnification factor. */
  double windowSize ;   /**< size of Gaussian window (in spatial bins) */

  vl_sift_pix *grad ;   /**< GSS gradient data. */
  int grad_o ;          /**< GSS gradient data octave. */

} VlSiftFilt ;

/** @name Create and destroy
 ** @{
 **/
VL_EXPORT
VlSiftFilt*  vl_sift_new    (int width, int height,
                             int noctaves, int nlevels,
                             int o_min) ;
VL_EXPORT
void         vl_sift_delete (VlSiftFilt *f) ;
/** @} */

/** @name Process data
 ** @{
 **/

VL_EXPORT
int   vl_sift_process_first_octave       (VlSiftFilt *f,
                                          vl_sift_pix const *im) ;

VL_EXPORT
int   vl_sift_process_next_octave        (VlSiftFilt *f) ;

VL_EXPORT
void  vl_sift_detect                     (VlSiftFilt *f) ;

VL_EXPORT
int   vl_sift_calc_keypoint_orientations (VlSiftFilt *f,
                                          double *angle,
//                                          double angles [4],
                                          VlSiftKeypoint const*k);
VL_EXPORT
void  vl_sift_calc_keypoint_descriptor   (VlSiftFilt *f,
                                          vl_sift_pix *descr,
                                          VlSiftKeypoint const* k,
                                          double angle) ;

VL_EXPORT
void  vl_sift_calc_raw_descriptor        (VlSiftFilt const *f,
                                          vl_sift_pix const* image,
                                          vl_sift_pix *descr,
                                          int widht, int height,
                                          double x, double y,
                                          double s, double angle0) ;

VL_EXPORT
void  vl_sift_keypoint_init              (VlSiftFilt const *f,
                                          VlSiftKeypoint *k,
                                          double x,
                                          double y,
                                          double sigma) ;

VL_EXPORT
void  vl_imconvcol_vf (vl_sift_pix* dst, vl_size dst_stride,
                       vl_sift_pix const* src,
                       vl_size src_width, vl_size src_height, vl_size src_stride,
                       vl_sift_pix const* filt, vl_index filt_begin, vl_index filt_end,
                       int step, unsigned int flags);
/** @} */

/** @name Retrieve data and parameters
 ** @{
 **/
VL_INLINE int    vl_sift_get_octave_index   (VlSiftFilt const *f) ;
VL_INLINE int    vl_sift_get_noctaves       (VlSiftFilt const *f) ;
VL_INLINE int    vl_sift_get_octave_first   (VlSiftFilt const *f) ;
VL_INLINE int    vl_sift_get_octave_width   (VlSiftFilt const *f) ;
VL_INLINE int    vl_sift_get_octave_height  (VlSiftFilt const *f) ;
VL_INLINE int    vl_sift_get_nlevels        (VlSiftFilt const *f) ;
VL_INLINE int    vl_sift_get_nkeypoints     (VlSiftFilt const *f) ;
VL_INLINE double vl_sift_get_peak_thresh    (VlSiftFilt const *f) ;
VL_INLINE double vl_sift_get_edge_thresh    (VlSiftFilt const *f) ;
VL_INLINE double vl_sift_get_norm_thresh    (VlSiftFilt const *f) ;
VL_INLINE double vl_sift_get_magnif         (VlSiftFilt const *f) ;
VL_INLINE double vl_sift_get_window_size    (VlSiftFilt const *f) ;

VL_INLINE vl_sift_pix *vl_sift_get_octave  (VlSiftFilt const *f, int s) ;
VL_INLINE VlSiftKeypoint const *vl_sift_get_keypoints (VlSiftFilt const *f) ;
/** @} */

/** @name Set parameters
 ** @{
 **/
VL_INLINE void vl_sift_set_peak_thresh (VlSiftFilt *f, double t) ;
VL_INLINE void vl_sift_set_edge_thresh (VlSiftFilt *f, double t) ;
VL_INLINE void vl_sift_set_norm_thresh (VlSiftFilt *f, double t) ;
VL_INLINE void vl_sift_set_magnif      (VlSiftFilt *f, double m) ;
VL_INLINE void vl_sift_set_window_size (VlSiftFilt *f, double m) ;
/** @} */

/* -------------------------------------------------------------------
 *                                     Inline functions implementation
 * ---------------------------------------------------------------- */

/** ------------------------------------------------------------------
 ** @brief Get current octave index.
 ** @param f SIFT filter.
 ** @return index of the current octave.
 **/

VL_INLINE int
vl_sift_get_octave_index (VlSiftFilt const *f)
{
  return f-> o_cur ;
}

/** ------------------------------------------------------------------
 ** @brief Get number of octaves.
 ** @param f SIFT filter.
 ** @return number of octaves.
 **/

VL_INLINE int
vl_sift_get_noctaves (VlSiftFilt const *f)
{
  return f-> O ;
}

/**-------------------------------------------------------------------
 ** @brief Get first octave.
 ** @param f SIFT filter.
 ** @return index of the first octave.
 **/

VL_INLINE int
vl_sift_get_octave_first (VlSiftFilt const *f)
{
  return f-> o_min ;
}

/** ------------------------------------------------------------------
 ** @brief Get current octave width
 ** @param f SIFT filter.
 ** @return current octave width.
 **/

VL_INLINE int
vl_sift_get_octave_width (VlSiftFilt const *f)
{
  return f-> octave_width ;
}

/** ------------------------------------------------------------------
 ** @brief Get current octave height
 ** @param f SIFT filter.
 ** @return current octave height.
 **/

VL_INLINE int
vl_sift_get_octave_height (VlSiftFilt const *f)
{
  return f-> octave_height ;
}

/** ------------------------------------------------------------------
 ** @brief Get current octave data
 ** @param f SIFT filter.
 ** @param s level index.
 **
 ** The level index @a s ranges in the interval <tt>s_min = -1</tt>
 ** and <tt> s_max = S + 2</tt>, where @c S is the number of levels
 ** per octave.
 **
 ** @return pointer to the octave data for level @a s.
 **/

VL_INLINE vl_sift_pix *
vl_sift_get_octave (VlSiftFilt const *f, int s)
{
  int w = vl_sift_get_octave_width  (f) ;
  int h = vl_sift_get_octave_height (f) ;
  return f->octave + w * h * (s - f->s_min) ;
}

/** ------------------------------------------------------------------
 ** @brief Get number of levels per octave
 ** @param f SIFT filter.
 ** @return number of leves per octave.
 **/

VL_INLINE int
vl_sift_get_nlevels (VlSiftFilt const *f)
{
  return f-> S ;
}

/** ------------------------------------------------------------------
 ** @brief Get number of keypoints.
 ** @param f SIFT filter.
 ** @return number of keypoints.
 **/

VL_INLINE int
vl_sift_get_nkeypoints (VlSiftFilt const *f)
{
  return f-> nkeys ;
}

/** ------------------------------------------------------------------
 ** @brief Get keypoints.
 ** @param f SIFT filter.
 ** @return pointer to the keypoints list.
 **/

VL_INLINE VlSiftKeypoint const *
vl_sift_get_keypoints (VlSiftFilt const *f)
{
  return f-> keys ;
}

/** ------------------------------------------------------------------
 ** @brief Get peaks treashold
 ** @param f SIFT filter.
 ** @return threshold ;
 **/

VL_INLINE double
vl_sift_get_peak_thresh (VlSiftFilt const *f)
{
  return f -> peak_thresh ;
}

/** ------------------------------------------------------------------
 ** @brief Get edges threshold
 ** @param f SIFT filter.
 ** @return threshold.
 **/

VL_INLINE double
vl_sift_get_edge_thresh (VlSiftFilt const *f)
{
  return f -> edge_thresh ;
}

/** ------------------------------------------------------------------
 ** @brief Get norm threshold
 ** @param f SIFT filter.
 ** @return threshold.
 **/

VL_INLINE double
vl_sift_get_norm_thresh (VlSiftFilt const *f)
{
  return f -> norm_thresh ;
}

/** ------------------------------------------------------------------
 ** @brief Get the magnification factor
 ** @param f SIFT filter.
 ** @return magnification factor.
 **/

VL_INLINE double
vl_sift_get_magnif (VlSiftFilt const *f)
{
  return f -> magnif ;
}

/** ------------------------------------------------------------------
 ** @brief Get the Gaussian window size.
 ** @param f SIFT filter.
 ** @return standard deviation of the Gaussian window (in spatial bin units).
 **/

VL_INLINE double
vl_sift_get_window_size (VlSiftFilt const *f)
{
  return f -> windowSize ;
}



/** ------------------------------------------------------------------
 ** @brief Set peaks threshold
 ** @param f SIFT filter.
 ** @param t threshold.
 **/

VL_INLINE void
vl_sift_set_peak_thresh (VlSiftFilt *f, double t)
{
  f -> peak_thresh = t ;
}

/** ------------------------------------------------------------------
 ** @brief Set edges threshold
 ** @param f SIFT filter.
 ** @param t threshold.
 **/

VL_INLINE void
vl_sift_set_edge_thresh (VlSiftFilt *f, double t)
{
  f -> edge_thresh = t ;
}

/** ------------------------------------------------------------------
 ** @brief Set norm threshold
 ** @param f SIFT filter.
 ** @param t threshold.
 **/

VL_INLINE void
vl_sift_set_norm_thresh (VlSiftFilt *f, double t)
{
  f -> norm_thresh = t ;
}

/** ------------------------------------------------------------------
 ** @brief Set the magnification factor
 ** @param f SIFT filter.
 ** @param m magnification factor.
 **/

VL_INLINE void
vl_sift_set_magnif (VlSiftFilt *f, double m)
{
  f -> magnif = m ;
}

/** ------------------------------------------------------------------
 ** @brief Set the Gaussian window size
 ** @param f SIFT filter.
 ** @param x Gaussian window size (in units of spatial bin).
 **
 ** This is the parameter @f$ \hat \sigma_{\text{win}} @f$ of
 ** the standard SIFT descriptor @ref sift-tech-descriptor-std.
 **/

VL_INLINE void
vl_sift_set_window_size (VlSiftFilt *f, double x)
{
  f -> windowSize = x ;
}

/* Additional inline functions are here */
// vl_floor_d vl_floor_f, see vlfeat/vl/mathop.h
VL_INLINE long int
vl_floor_d (double x)
{
  long int xi = (long int) x ;
  if (x >= 0 || (double) xi == x) return xi ;
  else return xi - 1 ;
}
VL_INLINE long int
vl_floor_f (float x)
{
  long int xi = (long int) x ;
  if (x >= 0 || (float) xi == x) return xi ;
  else return xi - 1 ;
}

// vl_abs_d, see mathop.h
VL_INLINE double
vl_abs_d (double x)
{
  return __builtin_fabs (x) ;
  // if not using gnuc, use below
  //return fabs(x) ;
}

// vl_fast_resqrt_f, see mathop.h
VL_INLINE float
vl_fast_resqrt_f (float x)
{
  /* 32-bit version */
  union {
    float x ;
    vl_int32  i ;
  } u ;

  float xhalf = (float) 0.5 * x ;

  /* convert floating point value in RAW integer */
  u.x = x ;

  /* gives initial guess y0 */
  u.i = 0x5f3759df - (u.i >> 1);
  /*u.i = 0xdf59375f - (u.i>>1);*/

  /* two Newton steps */
  u.x = u.x * ( (float) 1.5  - xhalf*u.x*u.x) ;
  u.x = u.x * ( (float) 1.5  - xhalf*u.x*u.x) ;
  return u.x ;
}

// vl_fast_sqrt_f, see mathop.h
VL_INLINE float
vl_fast_sqrt_f (float x)
{
  return (x < 1e-8) ? 0 : x * vl_fast_resqrt_f (x) ;
}

// vl_mod_2pi_f, see mathop.h
VL_INLINE float
vl_mod_2pi_f (float x)
{
  while (x > (float)(2 * VL_PI)) x -= (float) (2 * VL_PI) ;
  while (x < 0.0F) x += (float) (2 * VL_PI);
  return x ;
}

// vl_abs_f, see mathop.h
VL_INLINE float
vl_abs_f (float x)
{
  return __builtin_fabsf (x) ;
  // if not using gnuc, use below
  //return fabsf(x) ;
}

// vl_fast_atan2_f, see mathop.h
VL_INLINE float
vl_fast_atan2_f (float y, float x)
{
  float angle, r ;
  float const c3 = 0.1821F ;
  float const c1 = 0.9675F ;
  float abs_y    = vl_abs_f (y) + VL_EPSILON_F ;

  if (x >= 0) {
    r = (x - abs_y) / (x + abs_y) ;
    angle = (float) (VL_PI / 4) ;
  } else {
    r = (x + abs_y) / (abs_y - x) ;
    angle = (float) (3 * VL_PI / 4) ;
  }
  angle += (c3*r*r - c1) * r ;
  return (y < 0) ? - angle : angle ;
}


/* VL_SIFT_H */
#endif
