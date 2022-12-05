/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: FFTImplementationCallback.c
 *
 * MATLAB Coder version            : 5.4
 * C/C++ source code generated on  : 05-Dec-2022 10:00:34
 */

/* Include Files */
#include "FFTImplementationCallback.h"
#include "NEDwaves_emxutil.h"
#include "NEDwaves_types.h"
#include "rt_nonfinite.h"
#include <math.h>

/* Function Declarations */
static void c_FFTImplementationCallback_r2b(const emxArray_creal_T *x,
                                            int unsigned_nRows,
                                            const emxArray_real_T *costab,
                                            const emxArray_real_T *sintab,
                                            emxArray_creal_T *y);

static void d_FFTImplementationCallback_doH(
    const emxArray_real_T *x, int xoffInit, emxArray_creal_T *y, int nrowsx,
    int nRows, int nfft, const emxArray_creal_T *wwc,
    const emxArray_real_T *costab, const emxArray_real_T *sintab,
    const emxArray_real_T *costabinv, const emxArray_real_T *sintabinv);

/* Function Definitions */
/*
 * Arguments    : const emxArray_creal_T *x
 *                int unsigned_nRows
 *                const emxArray_real_T *costab
 *                const emxArray_real_T *sintab
 *                emxArray_creal_T *y
 * Return Type  : void
 */
static void c_FFTImplementationCallback_r2b(const emxArray_creal_T *x,
                                            int unsigned_nRows,
                                            const emxArray_real_T *costab,
                                            const emxArray_real_T *sintab,
                                            emxArray_creal_T *y)
{
  const creal_T *x_data;
  creal_T *y_data;
  const double *costab_data;
  const double *sintab_data;
  double temp_im;
  double temp_re;
  double temp_re_tmp;
  double twid_re;
  int i;
  int iDelta2;
  int iheight;
  int iy;
  int j;
  int ju;
  int k;
  int nRowsD2;
  sintab_data = sintab->data;
  costab_data = costab->data;
  x_data = x->data;
  iy = y->size[0];
  y->size[0] = unsigned_nRows;
  emxEnsureCapacity_creal_T(y, iy);
  y_data = y->data;
  if (unsigned_nRows > x->size[0]) {
    iy = y->size[0];
    y->size[0] = unsigned_nRows;
    emxEnsureCapacity_creal_T(y, iy);
    y_data = y->data;
    for (iy = 0; iy < unsigned_nRows; iy++) {
      y_data[iy].re = 0.0;
      y_data[iy].im = 0.0;
    }
  }
  iDelta2 = x->size[0];
  if (iDelta2 > unsigned_nRows) {
    iDelta2 = unsigned_nRows;
  }
  iheight = unsigned_nRows - 2;
  nRowsD2 = unsigned_nRows / 2;
  k = nRowsD2 / 2;
  iy = 0;
  ju = 0;
  for (i = 0; i <= iDelta2 - 2; i++) {
    bool tst;
    y_data[iy] = x_data[i];
    iy = unsigned_nRows;
    tst = true;
    while (tst) {
      iy >>= 1;
      ju ^= iy;
      tst = ((ju & iy) == 0);
    }
    iy = ju;
  }
  y_data[iy] = x_data[iDelta2 - 1];
  if (unsigned_nRows > 1) {
    for (i = 0; i <= iheight; i += 2) {
      temp_re_tmp = y_data[i + 1].re;
      temp_im = y_data[i + 1].im;
      temp_re = y_data[i].re;
      twid_re = y_data[i].im;
      y_data[i + 1].re = temp_re - temp_re_tmp;
      y_data[i + 1].im = twid_re - temp_im;
      y_data[i].re = temp_re + temp_re_tmp;
      y_data[i].im = twid_re + temp_im;
    }
  }
  iy = 2;
  iDelta2 = 4;
  iheight = ((k - 1) << 2) + 1;
  while (k > 0) {
    int b_temp_re_tmp;
    for (i = 0; i < iheight; i += iDelta2) {
      b_temp_re_tmp = i + iy;
      temp_re = y_data[b_temp_re_tmp].re;
      temp_im = y_data[b_temp_re_tmp].im;
      y_data[b_temp_re_tmp].re = y_data[i].re - temp_re;
      y_data[b_temp_re_tmp].im = y_data[i].im - temp_im;
      y_data[i].re += temp_re;
      y_data[i].im += temp_im;
    }
    ju = 1;
    for (j = k; j < nRowsD2; j += k) {
      double twid_im;
      int ihi;
      twid_re = costab_data[j];
      twid_im = sintab_data[j];
      i = ju;
      ihi = ju + iheight;
      while (i < ihi) {
        b_temp_re_tmp = i + iy;
        temp_re_tmp = y_data[b_temp_re_tmp].im;
        temp_im = y_data[b_temp_re_tmp].re;
        temp_re = twid_re * temp_im - twid_im * temp_re_tmp;
        temp_im = twid_re * temp_re_tmp + twid_im * temp_im;
        y_data[b_temp_re_tmp].re = y_data[i].re - temp_re;
        y_data[b_temp_re_tmp].im = y_data[i].im - temp_im;
        y_data[i].re += temp_re;
        y_data[i].im += temp_im;
        i += iDelta2;
      }
      ju++;
    }
    k /= 2;
    iy = iDelta2;
    iDelta2 += iDelta2;
    iheight -= iy;
  }
}

/*
 * Arguments    : const emxArray_real_T *x
 *                int xoffInit
 *                emxArray_creal_T *y
 *                int nrowsx
 *                int nRows
 *                int nfft
 *                const emxArray_creal_T *wwc
 *                const emxArray_real_T *costab
 *                const emxArray_real_T *sintab
 *                const emxArray_real_T *costabinv
 *                const emxArray_real_T *sintabinv
 * Return Type  : void
 */
static void d_FFTImplementationCallback_doH(
    const emxArray_real_T *x, int xoffInit, emxArray_creal_T *y, int nrowsx,
    int nRows, int nfft, const emxArray_creal_T *wwc,
    const emxArray_real_T *costab, const emxArray_real_T *sintab,
    const emxArray_real_T *costabinv, const emxArray_real_T *sintabinv)
{
  emxArray_creal_T *fv;
  emxArray_creal_T *fy;
  emxArray_creal_T *reconVar1;
  emxArray_creal_T *reconVar2;
  emxArray_creal_T *ytmp;
  emxArray_int32_T *wrapIndex;
  emxArray_real_T *b_costab;
  emxArray_real_T *b_sintab;
  emxArray_real_T *costab1q;
  emxArray_real_T *hcostabinv;
  emxArray_real_T *hsintab;
  emxArray_real_T *hsintabinv;
  const creal_T *wwc_data;
  creal_T *fv_data;
  creal_T *fy_data;
  creal_T *reconVar1_data;
  creal_T *reconVar2_data;
  creal_T *y_data;
  creal_T *ytmp_data;
  const double *costab_data;
  const double *costabinv_data;
  const double *sintab_data;
  const double *sintabinv_data;
  const double *x_data;
  double b_temp_re_tmp;
  double re_tmp;
  double temp_im;
  double temp_re;
  double twid_im;
  double twid_re;
  double z_tmp;
  double *b_costab_data;
  double *b_sintab_data;
  double *costab1q_data;
  double *hcostabinv_data;
  double *hsintab_data;
  double *hsintabinv_data;
  int hnRows;
  int hszCostab;
  int i;
  int istart;
  int j;
  int ju;
  int k;
  int nRowsD2;
  int nd2;
  int temp_re_tmp;
  int *wrapIndex_data;
  bool tst;
  sintabinv_data = sintabinv->data;
  costabinv_data = costabinv->data;
  sintab_data = sintab->data;
  costab_data = costab->data;
  wwc_data = wwc->data;
  y_data = y->data;
  x_data = x->data;
  emxInit_creal_T(&ytmp, 1);
  hnRows = nRows / 2;
  istart = ytmp->size[0];
  ytmp->size[0] = hnRows;
  emxEnsureCapacity_creal_T(ytmp, istart);
  ytmp_data = ytmp->data;
  if (hnRows > nrowsx) {
    istart = ytmp->size[0];
    ytmp->size[0] = hnRows;
    emxEnsureCapacity_creal_T(ytmp, istart);
    ytmp_data = ytmp->data;
    for (istart = 0; istart < hnRows; istart++) {
      ytmp_data[istart].re = 0.0;
      ytmp_data[istart].im = 0.0;
    }
  }
  if ((x->size[0] & 1) == 0) {
    tst = true;
    ju = x->size[0];
  } else if (x->size[0] >= nRows) {
    tst = true;
    ju = nRows;
  } else {
    tst = false;
    ju = x->size[0] - 1;
  }
  emxInit_real_T(&costab1q, 2);
  if (ju > nRows) {
    ju = nRows;
  }
  nd2 = nRows << 1;
  temp_im = 6.2831853071795862 / (double)nd2;
  j = nd2 / 2 / 2;
  istart = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = j + 1;
  emxEnsureCapacity_real_T(costab1q, istart);
  costab1q_data = costab1q->data;
  costab1q_data[0] = 1.0;
  nd2 = j / 2 - 1;
  for (k = 0; k <= nd2; k++) {
    costab1q_data[k + 1] = cos(temp_im * ((double)k + 1.0));
  }
  istart = nd2 + 2;
  nd2 = j - 1;
  for (k = istart; k <= nd2; k++) {
    costab1q_data[k] = sin(temp_im * (double)(j - k));
  }
  emxInit_real_T(&b_costab, 2);
  emxInit_real_T(&b_sintab, 2);
  costab1q_data[j] = 0.0;
  j = costab1q->size[1] - 1;
  nd2 = (costab1q->size[1] - 1) << 1;
  istart = b_costab->size[0] * b_costab->size[1];
  b_costab->size[0] = 1;
  b_costab->size[1] = nd2 + 1;
  emxEnsureCapacity_real_T(b_costab, istart);
  b_costab_data = b_costab->data;
  istart = b_sintab->size[0] * b_sintab->size[1];
  b_sintab->size[0] = 1;
  b_sintab->size[1] = nd2 + 1;
  emxEnsureCapacity_real_T(b_sintab, istart);
  b_sintab_data = b_sintab->data;
  b_costab_data[0] = 1.0;
  b_sintab_data[0] = 0.0;
  for (k = 0; k < j; k++) {
    b_costab_data[k + 1] = costab1q_data[k + 1];
    b_sintab_data[k + 1] = -costab1q_data[(j - k) - 1];
  }
  istart = costab1q->size[1];
  for (k = istart; k <= nd2; k++) {
    b_costab_data[k] = -costab1q_data[nd2 - k];
    b_sintab_data[k] = -costab1q_data[k - j];
  }
  emxInit_real_T(&hsintab, 2);
  emxInit_real_T(&hcostabinv, 2);
  emxInit_real_T(&hsintabinv, 2);
  hszCostab = costab->size[1] / 2;
  istart = costab1q->size[0] * costab1q->size[1];
  costab1q->size[0] = 1;
  costab1q->size[1] = hszCostab;
  emxEnsureCapacity_real_T(costab1q, istart);
  costab1q_data = costab1q->data;
  istart = hsintab->size[0] * hsintab->size[1];
  hsintab->size[0] = 1;
  hsintab->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hsintab, istart);
  hsintab_data = hsintab->data;
  istart = hcostabinv->size[0] * hcostabinv->size[1];
  hcostabinv->size[0] = 1;
  hcostabinv->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hcostabinv, istart);
  hcostabinv_data = hcostabinv->data;
  istart = hsintabinv->size[0] * hsintabinv->size[1];
  hsintabinv->size[0] = 1;
  hsintabinv->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hsintabinv, istart);
  hsintabinv_data = hsintabinv->data;
  for (i = 0; i < hszCostab; i++) {
    nd2 = ((i + 1) << 1) - 2;
    costab1q_data[i] = costab_data[nd2];
    hsintab_data[i] = sintab_data[nd2];
    hcostabinv_data[i] = costabinv_data[nd2];
    hsintabinv_data[i] = sintabinv_data[nd2];
  }
  emxInit_int32_T(&wrapIndex, 2);
  emxInit_creal_T(&reconVar1, 1);
  emxInit_creal_T(&reconVar2, 1);
  istart = reconVar1->size[0];
  reconVar1->size[0] = hnRows;
  emxEnsureCapacity_creal_T(reconVar1, istart);
  reconVar1_data = reconVar1->data;
  istart = reconVar2->size[0];
  reconVar2->size[0] = hnRows;
  emxEnsureCapacity_creal_T(reconVar2, istart);
  reconVar2_data = reconVar2->data;
  istart = wrapIndex->size[0] * wrapIndex->size[1];
  wrapIndex->size[0] = 1;
  wrapIndex->size[1] = hnRows;
  emxEnsureCapacity_int32_T(wrapIndex, istart);
  wrapIndex_data = wrapIndex->data;
  for (i = 0; i < hnRows; i++) {
    istart = i << 1;
    temp_im = b_sintab_data[istart];
    temp_re = b_costab_data[istart];
    reconVar1_data[i].re = temp_im + 1.0;
    reconVar1_data[i].im = -temp_re;
    reconVar2_data[i].re = 1.0 - temp_im;
    reconVar2_data[i].im = temp_re;
    if (i + 1 != 1) {
      wrapIndex_data[i] = (hnRows - i) + 1;
    } else {
      wrapIndex_data[0] = 1;
    }
  }
  emxFree_real_T(&b_sintab);
  emxFree_real_T(&b_costab);
  z_tmp = (double)ju / 2.0;
  istart = (int)((double)ju / 2.0);
  for (hszCostab = 0; hszCostab < istart; hszCostab++) {
    temp_re_tmp = (hnRows + hszCostab) - 1;
    temp_re = wwc_data[temp_re_tmp].re;
    temp_im = wwc_data[temp_re_tmp].im;
    nd2 = xoffInit + (hszCostab << 1);
    twid_re = x_data[nd2];
    twid_im = x_data[nd2 + 1];
    ytmp_data[hszCostab].re = temp_re * twid_re + temp_im * twid_im;
    ytmp_data[hszCostab].im = temp_re * twid_im - temp_im * twid_re;
  }
  if (!tst) {
    temp_re_tmp = (hnRows + (int)z_tmp) - 1;
    temp_re = wwc_data[temp_re_tmp].re;
    temp_im = wwc_data[temp_re_tmp].im;
    twid_re = x_data[xoffInit + ((int)z_tmp << 1)];
    ytmp_data[(int)z_tmp].re = temp_re * twid_re + temp_im * 0.0;
    ytmp_data[(int)z_tmp].im = temp_re * 0.0 - temp_im * twid_re;
    if ((int)z_tmp + 2 <= hnRows) {
      istart = (int)((double)ju / 2.0) + 2;
      for (i = istart; i <= hnRows; i++) {
        ytmp_data[i - 1].re = 0.0;
        ytmp_data[i - 1].im = 0.0;
      }
    }
  } else if ((int)z_tmp + 1 <= hnRows) {
    istart = (int)((double)ju / 2.0) + 1;
    for (i = istart; i <= hnRows; i++) {
      ytmp_data[i - 1].re = 0.0;
      ytmp_data[i - 1].im = 0.0;
    }
  }
  emxInit_creal_T(&fy, 1);
  z_tmp = (double)nfft / 2.0;
  nd2 = (int)z_tmp;
  istart = fy->size[0];
  fy->size[0] = (int)z_tmp;
  emxEnsureCapacity_creal_T(fy, istart);
  fy_data = fy->data;
  if ((int)z_tmp > ytmp->size[0]) {
    istart = fy->size[0];
    fy->size[0] = (int)z_tmp;
    emxEnsureCapacity_creal_T(fy, istart);
    fy_data = fy->data;
    for (istart = 0; istart < nd2; istart++) {
      fy_data[istart].re = 0.0;
      fy_data[istart].im = 0.0;
    }
  }
  ju = ytmp->size[0];
  istart = (int)z_tmp;
  if (ju <= istart) {
    istart = ju;
  }
  hszCostab = (int)z_tmp - 2;
  nRowsD2 = (int)z_tmp / 2;
  k = nRowsD2 / 2;
  nd2 = 0;
  ju = 0;
  for (i = 0; i <= istart - 2; i++) {
    fy_data[nd2] = ytmp_data[i];
    j = (int)z_tmp;
    tst = true;
    while (tst) {
      j >>= 1;
      ju ^= j;
      tst = ((ju & j) == 0);
    }
    nd2 = ju;
  }
  fy_data[nd2] = ytmp_data[istart - 1];
  if ((int)z_tmp > 1) {
    for (i = 0; i <= hszCostab; i += 2) {
      b_temp_re_tmp = fy_data[i + 1].re;
      temp_re = fy_data[i + 1].im;
      re_tmp = fy_data[i].re;
      twid_re = fy_data[i].im;
      fy_data[i + 1].re = re_tmp - b_temp_re_tmp;
      fy_data[i + 1].im = twid_re - temp_re;
      fy_data[i].re = re_tmp + b_temp_re_tmp;
      fy_data[i].im = twid_re + temp_re;
    }
  }
  nd2 = 2;
  hszCostab = 4;
  ju = ((k - 1) << 2) + 1;
  while (k > 0) {
    for (i = 0; i < ju; i += hszCostab) {
      temp_re_tmp = i + nd2;
      temp_re = fy_data[temp_re_tmp].re;
      temp_im = fy_data[temp_re_tmp].im;
      fy_data[temp_re_tmp].re = fy_data[i].re - temp_re;
      fy_data[temp_re_tmp].im = fy_data[i].im - temp_im;
      fy_data[i].re += temp_re;
      fy_data[i].im += temp_im;
    }
    istart = 1;
    for (j = k; j < nRowsD2; j += k) {
      int ihi;
      twid_re = costab1q_data[j];
      twid_im = hsintab_data[j];
      i = istart;
      ihi = istart + ju;
      while (i < ihi) {
        temp_re_tmp = i + nd2;
        b_temp_re_tmp = fy_data[temp_re_tmp].im;
        temp_im = fy_data[temp_re_tmp].re;
        temp_re = twid_re * temp_im - twid_im * b_temp_re_tmp;
        temp_im = twid_re * b_temp_re_tmp + twid_im * temp_im;
        fy_data[temp_re_tmp].re = fy_data[i].re - temp_re;
        fy_data[temp_re_tmp].im = fy_data[i].im - temp_im;
        fy_data[i].re += temp_re;
        fy_data[i].im += temp_im;
        i += hszCostab;
      }
      istart++;
    }
    k /= 2;
    nd2 = hszCostab;
    hszCostab += hszCostab;
    ju -= nd2;
  }
  emxInit_creal_T(&fv, 1);
  c_FFTImplementationCallback_r2b(wwc, (int)z_tmp, costab1q, hsintab, fv);
  fv_data = fv->data;
  nd2 = fy->size[0];
  emxFree_real_T(&costab1q);
  emxFree_real_T(&hsintab);
  for (istart = 0; istart < nd2; istart++) {
    re_tmp = fy_data[istart].re;
    twid_im = fv_data[istart].im;
    temp_im = fy_data[istart].im;
    temp_re = fv_data[istart].re;
    fy_data[istart].re = re_tmp * temp_re - temp_im * twid_im;
    fy_data[istart].im = re_tmp * twid_im + temp_im * temp_re;
  }
  c_FFTImplementationCallback_r2b(fy, (int)z_tmp, hcostabinv, hsintabinv, fv);
  fv_data = fv->data;
  emxFree_creal_T(&fy);
  emxFree_real_T(&hsintabinv);
  emxFree_real_T(&hcostabinv);
  if (fv->size[0] > 1) {
    temp_im = 1.0 / (double)fv->size[0];
    nd2 = fv->size[0];
    for (istart = 0; istart < nd2; istart++) {
      fv_data[istart].re *= temp_im;
      fv_data[istart].im *= temp_im;
    }
  }
  istart = wwc->size[0];
  for (k = hnRows; k <= istart; k++) {
    temp_im = wwc_data[k - 1].re;
    temp_re = fv_data[k - 1].im;
    twid_re = wwc_data[k - 1].im;
    twid_im = fv_data[k - 1].re;
    nd2 = k - hnRows;
    ytmp_data[nd2].re = temp_im * twid_im + twid_re * temp_re;
    ytmp_data[nd2].im = temp_im * temp_re - twid_re * twid_im;
  }
  emxFree_creal_T(&fv);
  for (i = 0; i < hnRows; i++) {
    double ytmp_re_tmp;
    istart = wrapIndex_data[i];
    temp_im = ytmp_data[i].re;
    temp_re = reconVar1_data[i].im;
    twid_re = ytmp_data[i].im;
    twid_im = reconVar1_data[i].re;
    re_tmp = ytmp_data[istart - 1].re;
    b_temp_re_tmp = -ytmp_data[istart - 1].im;
    z_tmp = reconVar2_data[i].im;
    ytmp_re_tmp = reconVar2_data[i].re;
    y_data[i].re = 0.5 * ((temp_im * twid_im - twid_re * temp_re) +
                          (re_tmp * ytmp_re_tmp - b_temp_re_tmp * z_tmp));
    y_data[i].im = 0.5 * ((temp_im * temp_re + twid_re * twid_im) +
                          (re_tmp * z_tmp + b_temp_re_tmp * ytmp_re_tmp));
    istart = hnRows + i;
    y_data[istart].re = 0.5 * ((temp_im * ytmp_re_tmp - twid_re * z_tmp) +
                               (re_tmp * twid_im - b_temp_re_tmp * temp_re));
    y_data[istart].im = 0.5 * ((temp_im * z_tmp + twid_re * ytmp_re_tmp) +
                               (re_tmp * temp_re + b_temp_re_tmp * twid_im));
  }
  emxFree_creal_T(&reconVar2);
  emxFree_creal_T(&reconVar1);
  emxFree_int32_T(&wrapIndex);
  emxFree_creal_T(&ytmp);
}

/*
 * Arguments    : const emxArray_real_T *x
 *                int xoffInit
 *                emxArray_creal_T *y
 *                int unsigned_nRows
 *                const emxArray_real_T *costab
 *                const emxArray_real_T *sintab
 * Return Type  : void
 */
void c_FFTImplementationCallback_doH(const emxArray_real_T *x, int xoffInit,
                                     emxArray_creal_T *y, int unsigned_nRows,
                                     const emxArray_real_T *costab,
                                     const emxArray_real_T *sintab)
{
  emxArray_creal_T *reconVar1;
  emxArray_creal_T *reconVar2;
  emxArray_int32_T *bitrevIndex;
  emxArray_int32_T *wrapIndex;
  emxArray_real_T *hcostab;
  emxArray_real_T *hsintab;
  creal_T *reconVar1_data;
  creal_T *reconVar2_data;
  creal_T *y_data;
  const double *costab_data;
  const double *sintab_data;
  const double *x_data;
  double b_y_re_tmp;
  double c_y_re_tmp;
  double d_y_re_tmp;
  double temp2_im;
  double temp2_re;
  double temp_im;
  double temp_im_tmp;
  double temp_re;
  double temp_re_tmp;
  double y_re_tmp;
  double z_tmp;
  double *hcostab_data;
  double *hsintab_data;
  int b_j1;
  int hszCostab;
  int i;
  int iDelta;
  int ihi;
  int istart;
  int ju;
  int k;
  int nRows;
  int nRowsD2;
  int *bitrevIndex_data;
  int *wrapIndex_data;
  bool tst;
  sintab_data = sintab->data;
  costab_data = costab->data;
  y_data = y->data;
  x_data = x->data;
  emxInit_real_T(&hcostab, 2);
  emxInit_real_T(&hsintab, 2);
  nRows = unsigned_nRows / 2;
  istart = y->size[0];
  if (istart > nRows) {
    istart = nRows;
  }
  ihi = nRows - 2;
  nRowsD2 = nRows / 2;
  k = nRowsD2 / 2;
  hszCostab = costab->size[1] / 2;
  b_j1 = hcostab->size[0] * hcostab->size[1];
  hcostab->size[0] = 1;
  hcostab->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hcostab, b_j1);
  hcostab_data = hcostab->data;
  b_j1 = hsintab->size[0] * hsintab->size[1];
  hsintab->size[0] = 1;
  hsintab->size[1] = hszCostab;
  emxEnsureCapacity_real_T(hsintab, b_j1);
  hsintab_data = hsintab->data;
  for (i = 0; i < hszCostab; i++) {
    iDelta = ((i + 1) << 1) - 2;
    hcostab_data[i] = costab_data[iDelta];
    hsintab_data[i] = sintab_data[iDelta];
  }
  emxInit_int32_T(&wrapIndex, 2);
  emxInit_creal_T(&reconVar1, 1);
  emxInit_creal_T(&reconVar2, 1);
  b_j1 = reconVar1->size[0];
  reconVar1->size[0] = nRows;
  emxEnsureCapacity_creal_T(reconVar1, b_j1);
  reconVar1_data = reconVar1->data;
  b_j1 = reconVar2->size[0];
  reconVar2->size[0] = nRows;
  emxEnsureCapacity_creal_T(reconVar2, b_j1);
  reconVar2_data = reconVar2->data;
  b_j1 = wrapIndex->size[0] * wrapIndex->size[1];
  wrapIndex->size[0] = 1;
  wrapIndex->size[1] = nRows;
  emxEnsureCapacity_int32_T(wrapIndex, b_j1);
  wrapIndex_data = wrapIndex->data;
  for (i = 0; i < nRows; i++) {
    temp_re = sintab_data[i];
    temp2_re = costab_data[i];
    reconVar1_data[i].re = temp_re + 1.0;
    reconVar1_data[i].im = -temp2_re;
    reconVar2_data[i].re = 1.0 - temp_re;
    reconVar2_data[i].im = temp2_re;
    if (i + 1 != 1) {
      wrapIndex_data[i] = (nRows - i) + 1;
    } else {
      wrapIndex_data[0] = 1;
    }
  }
  emxInit_int32_T(&bitrevIndex, 1);
  z_tmp = (double)unsigned_nRows / 2.0;
  ju = 0;
  hszCostab = 1;
  iDelta = (int)z_tmp;
  b_j1 = bitrevIndex->size[0];
  bitrevIndex->size[0] = (int)z_tmp;
  emxEnsureCapacity_int32_T(bitrevIndex, b_j1);
  bitrevIndex_data = bitrevIndex->data;
  for (b_j1 = 0; b_j1 < iDelta; b_j1++) {
    bitrevIndex_data[b_j1] = 0;
  }
  for (b_j1 = 0; b_j1 <= istart - 2; b_j1++) {
    bitrevIndex_data[b_j1] = hszCostab;
    iDelta = (int)z_tmp;
    tst = true;
    while (tst) {
      iDelta >>= 1;
      ju ^= iDelta;
      tst = ((ju & iDelta) == 0);
    }
    hszCostab = ju + 1;
  }
  bitrevIndex_data[istart - 1] = hszCostab;
  if ((x->size[0] & 1) == 0) {
    tst = true;
    istart = x->size[0];
  } else if (x->size[0] >= unsigned_nRows) {
    tst = true;
    istart = unsigned_nRows;
  } else {
    tst = false;
    istart = x->size[0] - 1;
  }
  if (istart <= unsigned_nRows) {
    iDelta = istart;
  } else {
    iDelta = unsigned_nRows;
  }
  temp_re = (double)iDelta / 2.0;
  if (istart > unsigned_nRows) {
    istart = unsigned_nRows;
  }
  b_j1 = (int)((double)istart / 2.0);
  for (i = 0; i < b_j1; i++) {
    hszCostab = xoffInit + (i << 1);
    y_data[bitrevIndex_data[i] - 1].re = x_data[hszCostab];
    y_data[bitrevIndex_data[i] - 1].im = x_data[hszCostab + 1];
  }
  if (!tst) {
    b_j1 = bitrevIndex_data[(int)temp_re] - 1;
    y_data[b_j1].re = x_data[xoffInit + ((int)temp_re << 1)];
    y_data[b_j1].im = 0.0;
  }
  emxFree_int32_T(&bitrevIndex);
  if (nRows > 1) {
    for (i = 0; i <= ihi; i += 2) {
      temp_re_tmp = y_data[i + 1].re;
      temp_im_tmp = y_data[i + 1].im;
      y_data[i + 1].re = y_data[i].re - temp_re_tmp;
      y_data[i + 1].im = y_data[i].im - y_data[i + 1].im;
      y_data[i].re += temp_re_tmp;
      y_data[i].im += temp_im_tmp;
    }
  }
  iDelta = 2;
  hszCostab = 4;
  ju = ((k - 1) << 2) + 1;
  while (k > 0) {
    for (i = 0; i < ju; i += hszCostab) {
      b_j1 = i + iDelta;
      temp_re = y_data[b_j1].re;
      temp_im = y_data[b_j1].im;
      y_data[b_j1].re = y_data[i].re - temp_re;
      y_data[b_j1].im = y_data[i].im - temp_im;
      y_data[i].re += temp_re;
      y_data[i].im += temp_im;
    }
    istart = 1;
    for (nRows = k; nRows < nRowsD2; nRows += k) {
      temp2_re = hcostab_data[nRows];
      temp2_im = hsintab_data[nRows];
      i = istart;
      ihi = istart + ju;
      while (i < ihi) {
        b_j1 = i + iDelta;
        temp_re_tmp = y_data[b_j1].im;
        y_re_tmp = y_data[b_j1].re;
        temp_re = temp2_re * y_re_tmp - temp2_im * temp_re_tmp;
        temp_im = temp2_re * temp_re_tmp + temp2_im * y_re_tmp;
        y_data[b_j1].re = y_data[i].re - temp_re;
        y_data[b_j1].im = y_data[i].im - temp_im;
        y_data[i].re += temp_re;
        y_data[i].im += temp_im;
        i += hszCostab;
      }
      istart++;
    }
    k /= 2;
    iDelta = hszCostab;
    hszCostab += hszCostab;
    ju -= iDelta;
  }
  emxFree_real_T(&hsintab);
  emxFree_real_T(&hcostab);
  iDelta = (int)z_tmp / 2;
  temp_re_tmp = y_data[0].re;
  temp_im_tmp = y_data[0].im;
  y_data[0].re = 0.5 * ((temp_re_tmp * reconVar1_data[0].re -
                         temp_im_tmp * reconVar1_data[0].im) +
                        (temp_re_tmp * reconVar2_data[0].re -
                         -temp_im_tmp * reconVar2_data[0].im));
  y_data[0].im = 0.5 * ((temp_re_tmp * reconVar1_data[0].im +
                         temp_im_tmp * reconVar1_data[0].re) +
                        (temp_re_tmp * reconVar2_data[0].im +
                         -temp_im_tmp * reconVar2_data[0].re));
  y_data[(int)z_tmp].re = 0.5 * ((temp_re_tmp * reconVar2_data[0].re -
                                  temp_im_tmp * reconVar2_data[0].im) +
                                 (temp_re_tmp * reconVar1_data[0].re -
                                  -temp_im_tmp * reconVar1_data[0].im));
  y_data[(int)z_tmp].im = 0.5 * ((temp_re_tmp * reconVar2_data[0].im +
                                  temp_im_tmp * reconVar2_data[0].re) +
                                 (temp_re_tmp * reconVar1_data[0].im +
                                  -temp_im_tmp * reconVar1_data[0].re));
  for (i = 2; i <= iDelta; i++) {
    double temp2_im_tmp;
    temp_re_tmp = y_data[i - 1].re;
    temp_im_tmp = y_data[i - 1].im;
    b_j1 = wrapIndex_data[i - 1];
    temp2_im = y_data[b_j1 - 1].re;
    temp2_im_tmp = y_data[b_j1 - 1].im;
    y_re_tmp = reconVar1_data[i - 1].im;
    b_y_re_tmp = reconVar1_data[i - 1].re;
    c_y_re_tmp = reconVar2_data[i - 1].im;
    d_y_re_tmp = reconVar2_data[i - 1].re;
    y_data[i - 1].re =
        0.5 * ((temp_re_tmp * b_y_re_tmp - temp_im_tmp * y_re_tmp) +
               (temp2_im * d_y_re_tmp - -temp2_im_tmp * c_y_re_tmp));
    y_data[i - 1].im =
        0.5 * ((temp_re_tmp * y_re_tmp + temp_im_tmp * b_y_re_tmp) +
               (temp2_im * c_y_re_tmp + -temp2_im_tmp * d_y_re_tmp));
    hszCostab = ((int)z_tmp + i) - 1;
    y_data[hszCostab].re =
        0.5 * ((temp_re_tmp * d_y_re_tmp - temp_im_tmp * c_y_re_tmp) +
               (temp2_im * b_y_re_tmp - -temp2_im_tmp * y_re_tmp));
    y_data[hszCostab].im =
        0.5 * ((temp_re_tmp * c_y_re_tmp + temp_im_tmp * d_y_re_tmp) +
               (temp2_im * y_re_tmp + -temp2_im_tmp * b_y_re_tmp));
    temp_im = reconVar1_data[b_j1 - 1].im;
    temp_re = reconVar1_data[b_j1 - 1].re;
    y_re_tmp = reconVar2_data[b_j1 - 1].im;
    temp2_re = reconVar2_data[b_j1 - 1].re;
    y_data[b_j1 - 1].re =
        0.5 * ((temp2_im * temp_re - temp2_im_tmp * temp_im) +
               (temp_re_tmp * temp2_re - -temp_im_tmp * y_re_tmp));
    y_data[b_j1 - 1].im =
        0.5 * ((temp2_im * temp_im + temp2_im_tmp * temp_re) +
               (temp_re_tmp * y_re_tmp + -temp_im_tmp * temp2_re));
    b_j1 = (b_j1 + (int)z_tmp) - 1;
    y_data[b_j1].re = 0.5 * ((temp2_im * temp2_re - temp2_im_tmp * y_re_tmp) +
                             (temp_re_tmp * temp_re - -temp_im_tmp * temp_im));
    y_data[b_j1].im = 0.5 * ((temp2_im * y_re_tmp + temp2_im_tmp * temp2_re) +
                             (temp_re_tmp * temp_im + -temp_im_tmp * temp_re));
  }
  emxFree_int32_T(&wrapIndex);
  if (iDelta != 0) {
    temp_re_tmp = y_data[iDelta].re;
    temp_im_tmp = y_data[iDelta].im;
    y_re_tmp = reconVar1_data[iDelta].im;
    b_y_re_tmp = reconVar1_data[iDelta].re;
    c_y_re_tmp = reconVar2_data[iDelta].im;
    d_y_re_tmp = reconVar2_data[iDelta].re;
    temp_re = temp_re_tmp * d_y_re_tmp;
    temp2_re = temp_re_tmp * b_y_re_tmp;
    y_data[iDelta].re = 0.5 * ((temp2_re - temp_im_tmp * y_re_tmp) +
                               (temp_re - -temp_im_tmp * c_y_re_tmp));
    temp2_im = temp_re_tmp * c_y_re_tmp;
    temp_im = temp_re_tmp * y_re_tmp;
    y_data[iDelta].im = 0.5 * ((temp_im + temp_im_tmp * b_y_re_tmp) +
                               (temp2_im + -temp_im_tmp * d_y_re_tmp));
    b_j1 = (int)z_tmp + iDelta;
    y_data[b_j1].re = 0.5 * ((temp_re - temp_im_tmp * c_y_re_tmp) +
                             (temp2_re - -temp_im_tmp * y_re_tmp));
    y_data[b_j1].im = 0.5 * ((temp2_im + temp_im_tmp * d_y_re_tmp) +
                             (temp_im + -temp_im_tmp * b_y_re_tmp));
  }
  emxFree_creal_T(&reconVar2);
  emxFree_creal_T(&reconVar1);
}

/*
 * Arguments    : const emxArray_real_T *x
 *                int n2blue
 *                int nfft
 *                const emxArray_real_T *costab
 *                const emxArray_real_T *sintab
 *                const emxArray_real_T *sintabinv
 *                emxArray_creal_T *y
 * Return Type  : void
 */
void c_FFTImplementationCallback_dob(const emxArray_real_T *x, int n2blue,
                                     int nfft, const emxArray_real_T *costab,
                                     const emxArray_real_T *sintab,
                                     const emxArray_real_T *sintabinv,
                                     emxArray_creal_T *y)
{
  emxArray_creal_T *fv;
  emxArray_creal_T *fy;
  emxArray_creal_T *r;
  emxArray_creal_T *wwc;
  creal_T *fv_data;
  creal_T *fy_data;
  creal_T *r1;
  creal_T *wwc_data;
  creal_T *y_data;
  const double *costab_data;
  const double *sintab_data;
  const double *x_data;
  double nt_im;
  double nt_re;
  int chan;
  int i;
  int ihi;
  int istart;
  int iy;
  int k;
  int nChan;
  int nInt2;
  int nInt2m1;
  int rt;
  bool ishalflength;
  sintab_data = sintab->data;
  costab_data = costab->data;
  x_data = x->data;
  nChan = x->size[1];
  emxInit_creal_T(&wwc, 1);
  if ((nfft != 1) && ((nfft & 1) == 0)) {
    istart = nfft / 2;
    nInt2m1 = (istart + istart) - 1;
    ihi = wwc->size[0];
    wwc->size[0] = nInt2m1;
    emxEnsureCapacity_creal_T(wwc, ihi);
    wwc_data = wwc->data;
    rt = 0;
    wwc_data[istart - 1].re = 1.0;
    wwc_data[istart - 1].im = 0.0;
    nInt2 = istart << 1;
    for (k = 0; k <= istart - 2; k++) {
      iy = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= iy) {
        rt += iy - nInt2;
      } else {
        rt += iy;
      }
      nt_im = -3.1415926535897931 * (double)rt / (double)istart;
      if (nt_im == 0.0) {
        nt_re = 1.0;
        nt_im = 0.0;
      } else {
        nt_re = cos(nt_im);
        nt_im = sin(nt_im);
      }
      ihi = (istart - k) - 2;
      wwc_data[ihi].re = nt_re;
      wwc_data[ihi].im = -nt_im;
    }
    ihi = nInt2m1 - 1;
    for (k = ihi; k >= istart; k--) {
      wwc_data[k] = wwc_data[(nInt2m1 - k) - 1];
    }
  } else {
    nInt2m1 = (nfft + nfft) - 1;
    ihi = wwc->size[0];
    wwc->size[0] = nInt2m1;
    emxEnsureCapacity_creal_T(wwc, ihi);
    wwc_data = wwc->data;
    rt = 0;
    wwc_data[nfft - 1].re = 1.0;
    wwc_data[nfft - 1].im = 0.0;
    nInt2 = nfft << 1;
    for (k = 0; k <= nfft - 2; k++) {
      iy = ((k + 1) << 1) - 1;
      if (nInt2 - rt <= iy) {
        rt += iy - nInt2;
      } else {
        rt += iy;
      }
      nt_im = -3.1415926535897931 * (double)rt / (double)nfft;
      if (nt_im == 0.0) {
        nt_re = 1.0;
        nt_im = 0.0;
      } else {
        nt_re = cos(nt_im);
        nt_im = sin(nt_im);
      }
      ihi = (nfft - k) - 2;
      wwc_data[ihi].re = nt_re;
      wwc_data[ihi].im = -nt_im;
    }
    ihi = nInt2m1 - 1;
    for (k = ihi; k >= nfft; k--) {
      wwc_data[k] = wwc_data[(nInt2m1 - k) - 1];
    }
  }
  ihi = y->size[0] * y->size[1];
  y->size[0] = nfft;
  y->size[1] = x->size[1];
  emxEnsureCapacity_creal_T(y, ihi);
  y_data = y->data;
  if (nfft > x->size[0]) {
    ihi = y->size[0] * y->size[1];
    y->size[0] = nfft;
    y->size[1] = x->size[1];
    emxEnsureCapacity_creal_T(y, ihi);
    y_data = y->data;
    iy = nfft * x->size[1];
    for (ihi = 0; ihi < iy; ihi++) {
      y_data[ihi].re = 0.0;
      y_data[ihi].im = 0.0;
    }
  }
  if (x->size[1] - 1 >= 0) {
    if ((n2blue != 1) && ((nfft & 1) == 0)) {
      ishalflength = true;
    } else {
      ishalflength = false;
    }
  }
  emxInit_creal_T(&r, 1);
  emxInit_creal_T(&fy, 1);
  emxInit_creal_T(&fv, 1);
  for (chan = 0; chan < nChan; chan++) {
    iy = chan * x->size[0];
    ihi = r->size[0];
    r->size[0] = nfft;
    emxEnsureCapacity_creal_T(r, ihi);
    r1 = r->data;
    if (nfft > x->size[0]) {
      ihi = r->size[0];
      r->size[0] = nfft;
      emxEnsureCapacity_creal_T(r, ihi);
      r1 = r->data;
      for (ihi = 0; ihi < nfft; ihi++) {
        r1[ihi].re = 0.0;
        r1[ihi].im = 0.0;
      }
    }
    if (ishalflength) {
      d_FFTImplementationCallback_doH(x, iy, r, x->size[0], nfft, n2blue, wwc,
                                      costab, sintab, costab, sintabinv);
      r1 = r->data;
    } else {
      double b_nt_re_tmp;
      double twid_im;
      double twid_re;
      int nRowsD2;
      int nt_re_tmp;
      nInt2m1 = x->size[0];
      if (nfft <= nInt2m1) {
        nInt2m1 = nfft;
      }
      for (k = 0; k < nInt2m1; k++) {
        nt_re_tmp = (nfft + k) - 1;
        ihi = iy + k;
        r1[k].re = wwc_data[nt_re_tmp].re * x_data[ihi];
        r1[k].im = wwc_data[nt_re_tmp].im * -x_data[ihi];
      }
      ihi = nInt2m1 + 1;
      for (k = ihi; k <= nfft; k++) {
        r1[k - 1].re = 0.0;
        r1[k - 1].im = 0.0;
      }
      ihi = fy->size[0];
      fy->size[0] = n2blue;
      emxEnsureCapacity_creal_T(fy, ihi);
      fy_data = fy->data;
      if (n2blue > r->size[0]) {
        ihi = fy->size[0];
        fy->size[0] = n2blue;
        emxEnsureCapacity_creal_T(fy, ihi);
        fy_data = fy->data;
        for (ihi = 0; ihi < n2blue; ihi++) {
          fy_data[ihi].re = 0.0;
          fy_data[ihi].im = 0.0;
        }
      }
      istart = r->size[0];
      if (istart > n2blue) {
        istart = n2blue;
      }
      rt = n2blue - 2;
      nRowsD2 = n2blue / 2;
      k = nRowsD2 / 2;
      iy = 0;
      nInt2 = 0;
      for (i = 0; i <= istart - 2; i++) {
        bool tst;
        fy_data[iy] = r1[i];
        nInt2m1 = n2blue;
        tst = true;
        while (tst) {
          nInt2m1 >>= 1;
          nInt2 ^= nInt2m1;
          tst = ((nInt2 & nInt2m1) == 0);
        }
        iy = nInt2;
      }
      fy_data[iy] = r1[istart - 1];
      if (n2blue > 1) {
        for (i = 0; i <= rt; i += 2) {
          b_nt_re_tmp = fy_data[i + 1].re;
          nt_im = fy_data[i + 1].im;
          twid_im = fy_data[i].re;
          nt_re = fy_data[i].im;
          fy_data[i + 1].re = twid_im - b_nt_re_tmp;
          fy_data[i + 1].im = nt_re - nt_im;
          fy_data[i].re = twid_im + b_nt_re_tmp;
          fy_data[i].im = nt_re + nt_im;
        }
      }
      nInt2m1 = 2;
      rt = 4;
      nInt2 = ((k - 1) << 2) + 1;
      while (k > 0) {
        for (i = 0; i < nInt2; i += rt) {
          nt_re_tmp = i + nInt2m1;
          nt_re = fy_data[nt_re_tmp].re;
          nt_im = fy_data[nt_re_tmp].im;
          fy_data[nt_re_tmp].re = fy_data[i].re - nt_re;
          fy_data[nt_re_tmp].im = fy_data[i].im - nt_im;
          fy_data[i].re += nt_re;
          fy_data[i].im += nt_im;
        }
        istart = 1;
        for (iy = k; iy < nRowsD2; iy += k) {
          twid_re = costab_data[iy];
          twid_im = sintab_data[iy];
          i = istart;
          ihi = istart + nInt2;
          while (i < ihi) {
            nt_re_tmp = i + nInt2m1;
            b_nt_re_tmp = fy_data[nt_re_tmp].im;
            nt_im = fy_data[nt_re_tmp].re;
            nt_re = twid_re * nt_im - twid_im * b_nt_re_tmp;
            nt_im = twid_re * b_nt_re_tmp + twid_im * nt_im;
            fy_data[nt_re_tmp].re = fy_data[i].re - nt_re;
            fy_data[nt_re_tmp].im = fy_data[i].im - nt_im;
            fy_data[i].re += nt_re;
            fy_data[i].im += nt_im;
            i += rt;
          }
          istart++;
        }
        k /= 2;
        nInt2m1 = rt;
        rt += rt;
        nInt2 -= nInt2m1;
      }
      c_FFTImplementationCallback_r2b(wwc, n2blue, costab, sintab, fv);
      fv_data = fv->data;
      iy = fy->size[0];
      for (ihi = 0; ihi < iy; ihi++) {
        twid_im = fy_data[ihi].re;
        nt_im = fv_data[ihi].im;
        nt_re = fy_data[ihi].im;
        twid_re = fv_data[ihi].re;
        fy_data[ihi].re = twid_im * twid_re - nt_re * nt_im;
        fy_data[ihi].im = twid_im * nt_im + nt_re * twid_re;
      }
      c_FFTImplementationCallback_r2b(fy, n2blue, costab, sintabinv, fv);
      fv_data = fv->data;
      if (fv->size[0] > 1) {
        nt_im = 1.0 / (double)fv->size[0];
        iy = fv->size[0];
        for (ihi = 0; ihi < iy; ihi++) {
          fv_data[ihi].re *= nt_im;
          fv_data[ihi].im *= nt_im;
        }
      }
      ihi = wwc->size[0];
      for (k = nfft; k <= ihi; k++) {
        nt_im = wwc_data[k - 1].re;
        nt_re = fv_data[k - 1].im;
        twid_re = wwc_data[k - 1].im;
        twid_im = fv_data[k - 1].re;
        iy = k - nfft;
        r1[iy].re = nt_im * twid_im + twid_re * nt_re;
        r1[iy].im = nt_im * nt_re - twid_re * twid_im;
      }
    }
    iy = r->size[0];
    for (ihi = 0; ihi < iy; ihi++) {
      y_data[ihi + y->size[0] * chan] = r1[ihi];
    }
  }
  emxFree_creal_T(&fv);
  emxFree_creal_T(&fy);
  emxFree_creal_T(&r);
  emxFree_creal_T(&wwc);
}

/*
 * File trailer for FFTImplementationCallback.c
 *
 * [EOF]
 */
