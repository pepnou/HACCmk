//
// Argonne Leadership Computing Facility
// Vitali Morozov morozov@anl.gov
//


#include <math.h>
#include <stdint.h>
#include <immintrin.h>
#include <stdio.h>

float Q_rsqrt( float number )
{	
	const float x2 = number * 0.5F;
	const float threehalfs = 1.5F;

	union {
		float f;
		uint32_t i;
	} conv = {number}; // member 'f' set to value of 'number'.
	conv.i  = 0x5f3759df - ( conv.i >> 1 );
	conv.f  *= ( threehalfs - ( x2 * conv.f * conv.f ) );
	return conv.f;
}

/*void Step10_orig( int count1, float xxi, float yyi, float zzi, float fsrrmax2, float mp_rsm2, float *xx1, float *yy1, float *zz1, float *mass1, float *dxi, float *dyi, float *dzi )
{

    const float ma0 = 0.269327, ma1 = -0.0750978, ma2 = 0.0114808, ma3 = -0.00109313, ma4 = 0.0000605491, ma5 = -0.00000147177;
    
    float dxc, dyc, dzc, m, r2, f, xi, yi, zi;
    int j;

    xi = 0.; yi = 0.; zi = 0.;

    for ( j = 0; j < count1; j++ ) 
    {
        dxc = xx1[j] - xxi;
        dyc = yy1[j] - yyi;
        dzc = zz1[j] - zzi;
  
        r2 = dxc * dxc + dyc * dyc + dzc * dzc;
       
        m = ( r2 < fsrrmax2 ) ? mass1[j] : 0.0f;

        f =  pow( r2 + mp_rsm2, -1.5 ) - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));
        
        f = ( r2 > 0.0f ) ? m * f : 0.0f;

        xi = xi + f * dxc;
        yi = yi + f * dyc;
        zi = zi + f * dzc;
    }

    *dxi = xi;
    *dyi = yi;
    *dzi = zi;
}*/

/*void Step10_orig( int count1, int i, float fsrrmax2, float mp_rsm2, float *xx1, float *yy1, float *zz1, float *mass1, float *dxi, float *dyi, float *dzi )
{

  const float ma0 = 0.269327, ma1 = -0.0750978, ma2 = 0.0114808, ma3 = -0.00109313, ma4 = 0.0000605491, ma5 = -0.00000147177;

  float dxc, dyc, dzc, m, r2, f, xi, yi, zi;
  int j;

  xi = 0.; yi = 0.; zi = 0.;

  __m256 r_xxi = _mm256_broadcast_ss(&(xx1[i])),
         r_yyi = _mm256_broadcast_ss(&(yy1[i])),
         r_zzi = _mm256_broadcast_ss(&(zz1[i])),
         r_ma0 = _mm256_broadcast_ss(&ma0),
         r_ma1 = _mm256_broadcast_ss(&ma1),
         r_ma2 = _mm256_broadcast_ss(&ma2),
         r_ma3 = _mm256_broadcast_ss(&ma3),
         r_ma4 = _mm256_broadcast_ss(&ma4),
         r_ma5 = _mm256_broadcast_ss(&ma5),
         r_mp_rsm2 = _mm256_broadcast_ss(&mp_rsm2),
         r_fssrmax2 = _mm256_broadcast_ss(&fsrrmax2);

  float massi = mass1[i];
  mass1[i] = 0.f;

  float vectx[8], vecty[8], vectz[8];

  for ( j = 0; j < count1; j += 8 ) 
  {
    //dxc = xx1[j] - xxi;
    //dyc = yy1[j] - yyi;
    //dzc = zz1[j] - zzi;
    __m256 r_dxc = _mm256_sub_ps(
        _mm256_load_ps(
          &(xx1[j])) , 
        r_xxi);
    __m256 r_dyc = _mm256_sub_ps(
        _mm256_load_ps(
          &(yy1[j])) , 
        r_yyi);
    __m256 r_dzc = _mm256_sub_ps(
        _mm256_load_ps(
          &(zz1[j])) , 
        r_zzi);

    //r2 = dxc * dxc + dyc * dyc + dzc * dzc;
    __m256 r_r2 = _mm256_add_ps(
        _mm256_mul_ps(
          r_dxc,
          r_dxc),
        _mm256_add_ps(
          _mm256_mul_ps(
            r_dyc,
            r_dyc),
          _mm256_mul_ps(
            r_dzc,
            r_dzc)));

    //m = ( r2 < fsrrmax2 ) ? mass1[j] : 0.0f;
    __m256i mask = _mm256_cvtps_epi32(_mm256_cmp_ps(r_r2, r_fssrmax2, _CMP_LT_OQ));
    __m256 r_m = _mm256_maskload_ps(&(mass1[j]), mask);


    //f =  pow( r2 + mp_rsm2, -1.5 ) - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));
    __m256 r_r2m = _mm256_add_ps(r_r2, r_mp_rsm2);
    __m256 r_f = _mm256_mul_ps(
        r_m ,
        _mm256_sub_ps( 
          _mm256_div_ps(
            _mm256_rsqrt_ps(r_r2m), 
            r_r2m),
          _mm256_add_ps( 
            r_ma0, 
            _mm256_mul_ps(
              r_r2,
              _mm256_add_ps(
                r_ma1,
                _mm256_mul_ps(
                  r_r2,
                  _mm256_add_ps(
                    r_ma2,
                    _mm256_mul_ps(
                      r_r2,
                      _mm256_add_ps(
                        r_ma3,
                        _mm256_mul_ps(
                          r_r2,
                          _mm256_add_ps(
                            r_ma4,
                            _mm256_mul_ps(
                              r_r2,
                              r_ma5))))))))))));

    //f = ( r2 > 0.0f ) ? m * f : 0.0f;
    
    _mm256_store_ps(
        vectx,
        _mm256_mul_ps(
          r_dxc,
          r_f));

    _mm256_store_ps(
        vecty,
        _mm256_mul_ps(
          r_dyc,
          r_f));

    _mm256_store_ps(
        vectz,
        _mm256_mul_ps(
          r_dzc,
          r_f));
    
    for(int l = 0; l < 8; l++) {
      xi += vectx[l];
      yi += vecty[l];
      zi += vectz[l];
    }
    

    //xi = xi + f * dxc;
    //yi = yi + f * dyc;
    //zi = zi + f * dzc;
  }

  mass1[i] = massi;

  *dxi = xi;
  *dyi = yi;
  *dzi = zi;
}*/


/*void Step10_orig( int count1, int i, float fsrrmax2, float mp_rsm2, float *xx1, float *mass1, float *dxi )
{

    const float ma0 = 0.269327, ma1 = -0.0750978, ma2 = 0.0114808, ma3 = -0.00109313, ma4 = 0.0000605491, ma5 = -0.00000147177;
    
    float dxc, m, r2, f, xi, invsqrt, r2m;
    int j;

    xi = 0.;

    for ( j = 0; j < count1; j++ ) 
    {
        dxc = xx1[j] - xx1[i];
  
        r2 = (dxc * dxc);
       
        m = ( r2 < fsrrmax2 ) ? mass1[j] : 0.0f;

        //f =  pow( r2 + mp_rsm2, -1.5 ) - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));
        //f =  1 / ((r2 + mp_rsm2) * sqrt(r2 + mp_rsm2)) - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));
        //f =  pow( r2 + mp_rsm2, -1.5 ) - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));
        //f =  Q_rsqrt(r2 + mp_rsm2) / (r2 + mp_rsm2) - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));
        
        //r2m = r2 + mp_rsm2;
        //__m128 in = _mm_load_ss( &r2m );
        //_mm_store_ss( &invsqrt, _mm_rsqrt_ss( in ) );
        
        //f =  invsqrt / r2m - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));

        f = ( r2 > 0.0f ) ? m * f : 0.0f;

        xi = xi + f * dxc;
    }

    *dxi = xi;
}*/

void Step10_orig( int count1, int i, float fsrrmax2, float mp_rsm2, float *xx1, float *mass1, float *dxi)
{

  const float ma0 = 0.269327, ma1 = -0.0750978, ma2 = 0.0114808, ma3 = -0.00109313, ma4 = 0.0000605491, ma5 = -0.00000147177;

  float dxc, m, r2, f, xi, invsqrt, r2m;
  int j;

  xi = 0.;

  __m256 r_xxi = _mm256_broadcast_ss(&(xx1[i])),
         r_ma0 = _mm256_broadcast_ss(&ma0),
         r_ma1 = _mm256_broadcast_ss(&ma1),
         r_ma2 = _mm256_broadcast_ss(&ma2),
         r_ma3 = _mm256_broadcast_ss(&ma3),
         r_ma4 = _mm256_broadcast_ss(&ma4),
         r_ma5 = _mm256_broadcast_ss(&ma5),
         r_mp_rsm2 = _mm256_broadcast_ss(&mp_rsm2),
         r_fsrrmax2 = _mm256_broadcast_ss(&fsrrmax2);

  float massi = mass1[i];
  mass1[i] = 0.;

  float vect[8];
   

  for ( j = 0; j < count1; j += 8 ) 
  {
    //dxc = xx1[j] - xxi;
    __m256 r_dxc = _mm256_sub_ps(
        _mm256_load_ps(
          &(xx1[j])) , 
        r_xxi);

    /*if(iter == 400 && j == 0 && i == 0) {
      _mm256_store_ps(vect, r_dxc);
    
      for(int l = 0; l < 8; l++) {
        printf("%f\n",vect[l]);
      }

      printf("\n");
    }*/


    //r2 = (dxc * dxc);
    __m256 r_r2 = _mm256_mul_ps(r_dxc, r_dxc);

    /*if(iter == 400 && j == 0 && i == 0) {
      _mm256_store_ps(vect, r_r2);
    
      for(int l = 0; l < 8; l++) {
        printf("%f\n",vect[l]);
      }

      printf("\n");
    }*/

    //m = ( r2 < fsrrmax2 ) ? mass1[j] : 0.0f;
    __m256i mask = _mm256_cvtps_epi32(_mm256_cmp_ps(r_r2, r_fsrrmax2, _CMP_LT_OQ));
    __m256 r_m = _mm256_maskload_ps(&(mass1[j]), mask);

    //f =  Q_rsqrt(r2 + mp_rsm2) / (r2 + mp_rsm2) - ( ma0 + r2*(ma1 + r2*(ma2 + r2*(ma3 + r2*(ma4 + r2*ma5)))));
    
    __m256 r_r2m = _mm256_add_ps(r_r2, r_mp_rsm2);

    /*if(iter == 400 && j == 0 && i == 0) {
      _mm256_store_ps(vect, r_r2m);
    
      for(int l = 0; l < 8; l++) {
        printf("%f\n",vect[l]);
      }

      printf("\n");
    }*/

    __m256 r_f = _mm256_sub_ps( 
      _mm256_div_ps(
        _mm256_rsqrt_ps(r_r2m), 
        r_r2m),
      _mm256_add_ps( 
        r_ma0, 
        _mm256_mul_ps(
          r_r2,
          _mm256_add_ps(
            r_ma1,
            _mm256_mul_ps(
              r_r2,
              _mm256_add_ps(
                r_ma2,
                _mm256_mul_ps(
                  r_r2,
                  _mm256_add_ps(
                    r_ma3,
                    _mm256_mul_ps(
                      r_r2,
                      _mm256_add_ps(
                        r_ma4,
                        _mm256_mul_ps(
                          r_r2,
                          r_ma5)))))))))));

    /*if(iter == 400 && j == 0 && i == 0) {
      _mm256_store_ps(vect, r_f);
    
      for(int l = 0; l < 8; l++) {
        printf("%f\n",vect[l]);
      }

      printf("\n");
    }*/

    

    //f = ( r2 > 0.0f ) ? m * f : 0.0f;
    
    _mm256_store_ps(
        vect, 
        _mm256_mul_ps(
          r_dxc, 
          _mm256_mul_ps(
            r_f, 
            r_m)));
    
    /*if(iter == 400 && j == 0 && i == 0) {
      for(int l = 0; l < 8; l++) {
        printf("%f\n",vect[l]);
      }

      printf("\n");
    }*/

    for(int l = 0; l < 8; l++) {
      xi += vect[l];
    }


    //xi = xi + f * dxc;
  }

  mass1[i] = massi;

  *dxi = xi;


}


