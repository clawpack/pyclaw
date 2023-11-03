
#define PY_ARRAY_UNIQUE_SYMBOL reconstruct_ARRAY_API

#include <stdio.h>

#include <Python.h>
#include <numpy/ndarrayobject.h>

void
smoothness_k3 (const double *restrict f, int n, int fsi,
	       double *restrict sigma, int ssi, int ssr)
{
  int i;
  double sigma0, sigma1, sigma2;
  for (i = 2; i < n - 2; i++)
    {
      sigma0 =
	+3.33333333333333 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	10.3333333333333 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	3.66666666666667 * f[(i + 0) * fsi] * f[(i + 2) * fsi] +
	8.33333333333333 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	6.33333333333333 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	1.33333333333333 * f[(i + 2) * fsi] * f[(i + 2) * fsi];
      sigma1 =
	+1.33333333333333 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	4.33333333333333 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	1.66666666666667 * f[(i - 1) * fsi] * f[(i + 1) * fsi] +
	4.33333333333333 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	4.33333333333333 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	1.33333333333333 * f[(i + 1) * fsi] * f[(i + 1) * fsi];
      sigma2 =
	+1.33333333333333 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	6.33333333333333 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	3.66666666666667 * f[(i - 2) * fsi] * f[(i + 0) * fsi] +
	8.33333333333333 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	10.3333333333333 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	3.33333333333333 * f[(i + 0) * fsi] * f[(i + 0) * fsi];
      sigma[i * ssi + 0 * ssr] = sigma0;
      sigma[i * ssi + 1 * ssr] = sigma1;
      sigma[i * ssi + 2 * ssr] = sigma2;
    }
}


PyObject *
py_smoothness_k3 (PyObject * self, PyObject * args)
{
  double *sigma, *f;
  PyArrayObject *f_py, *sigma_py;

  long int n;
  int ssi, ssr, fsi;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &f_py, &sigma_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  f = (double *) PyArray_DATA (f_py);

  n = PyArray_DIM (f_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);

  smoothness_k3 (f, n, fsi, sigma, ssi, ssr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_left_k3 (const double *restrict sigma, int n, int ssi, int ssr,
		 double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2;
  double omega0, omega1, omega2;
  for (i = 2; i < n - 2; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      accumulator = 0.0;
      omega0 = +0.1 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.6 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.3 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
    }
}



PyObject *
py_weights_left_k3 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_left_k3 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_left_k3 (const double *restrict f, int n, int fsi,
		     const double *restrict omega, int wsi, int wsl, int wsr,
		     double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2;
  double omega0, omega1, omega2;
  for (i = 2; i < n - 2; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      fr0 =
	+1.83333333333333 * f[(i + 0) * fsi] +
	-1.16666666666667 * f[(i + 1) * fsi] +
	+0.333333333333333 * f[(i + 2) * fsi];
      fr1 =
	+0.333333333333333 * f[(i - 1) * fsi] +
	+0.833333333333333 * f[(i + 0) * fsi] +
	-0.166666666666667 * f[(i + 1) * fsi];
      fr2 =
	-0.166666666666667 * f[(i - 2) * fsi] +
	+0.833333333333333 * f[(i - 1) * fsi] +
	+0.333333333333333 * f[(i + 0) * fsi];

      fr[i * frsi + 0 * frsl] = fr0 * omega0 + fr1 * omega1 + fr2 * omega2;

      /* debugging */
      /* printf("i %d f %lf %lf %lf %lf %lf w %lf %lf %lf fr %lf\n", i, */
      /*        f[(i - 2) * fsi], */
      /*        f[(i - 1) * fsi], */
      /*        f[(i - 0) * fsi], */
      /*        f[(i + 1) * fsi], */
      /*        f[(i + 2) * fsi], */
      /*        omega0, omega1, omega2, */
      /*        fr[i * frsi + 0 * frsl]); */
    }
}

PyObject *
py_reconstruct_left_k3 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_left_k3 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_right_k3 (const double *restrict sigma, int n, int ssi, int ssr,
		  double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2;
  double omega0, omega1, omega2;
  for (i = 2; i < n - 2; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      accumulator = 0.0;
      omega0 = +0.3 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.6 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.1 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
    }
}



PyObject *
py_weights_right_k3 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_right_k3 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_right_k3 (const double *restrict f, int n, int fsi,
		      const double *restrict omega, int wsi, int wsl, int wsr,
		      double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2;
  double omega0, omega1, omega2;
  for (i = 2; i < n - 2; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      fr0 =
	+0.333333333333333 * f[(i + 0) * fsi] +
	+0.833333333333333 * f[(i + 1) * fsi] +
	-0.166666666666667 * f[(i + 2) * fsi];
      fr1 =
	-0.166666666666667 * f[(i - 1) * fsi] +
	+0.833333333333333 * f[(i + 0) * fsi] +
	+0.333333333333333 * f[(i + 1) * fsi];
      fr2 =
	+0.333333333333333 * f[(i - 2) * fsi] +
	-1.16666666666667 * f[(i - 1) * fsi] +
	+1.83333333333333 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] = fr0 * omega0 + fr1 * omega1 + fr2 * omega2;
    }
}

PyObject *
py_reconstruct_right_k3 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_right_k3 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
smoothness_k4 (const double *restrict f, int n, int fsi,
	       double *restrict sigma, int ssi, int ssr)
{
  int i;
  double sigma0, sigma1, sigma2, sigma3;
  for (i = 3; i < n - 3; i++)
    {
      sigma0 =
	+8.77916666666667 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	39.175 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	29.3416666666667 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	7.725 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	45.8458333333333 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	71.8583333333333 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	19.3416666666667 * f[(i + 1) * fsi] * f[(i + 3) * fsi] +
	29.3458333333333 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	16.175 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	2.27916666666667 * f[(i + 3) * fsi] * f[(i + 3) * fsi];
      sigma1 =
	+2.27916666666667 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	10.5083333333333 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	8.00833333333333 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	2.05833333333333 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	14.3458333333333 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	24.8583333333333 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	6.675 * f[(i + 0) * fsi] * f[(i + 2) * fsi] +
	11.8458333333333 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	6.84166666666667 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	1.1125 * f[(i + 2) * fsi] * f[(i + 2) * fsi];
      sigma2 =
	+1.1125 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	6.84166666666667 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	6.675 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	2.05833333333333 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	11.8458333333333 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	24.8583333333333 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	8.00833333333333 * f[(i - 1) * fsi] * f[(i + 1) * fsi] +
	14.3458333333333 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	10.5083333333333 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	2.27916666666667 * f[(i + 1) * fsi] * f[(i + 1) * fsi];
      sigma3 =
	+2.27916666666667 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	16.175 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	19.3416666666667 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	7.725 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	29.3458333333333 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	71.8583333333333 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	29.3416666666667 * f[(i - 2) * fsi] * f[(i + 0) * fsi] +
	45.8458333333333 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	39.175 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	8.77916666666667 * f[(i + 0) * fsi] * f[(i + 0) * fsi];
      sigma[i * ssi + 0 * ssr] = sigma0;
      sigma[i * ssi + 1 * ssr] = sigma1;
      sigma[i * ssi + 2 * ssr] = sigma2;
      sigma[i * ssi + 3 * ssr] = sigma3;
    }
}


PyObject *
py_smoothness_k4 (PyObject * self, PyObject * args)
{
  double *sigma, *f;
  PyArrayObject *f_py, *sigma_py;

  long int n;
  int ssi, ssr, fsi;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &f_py, &sigma_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  f = (double *) PyArray_DATA (f_py);

  n = PyArray_DIM (f_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);

  smoothness_k4 (f, n, fsi, sigma, ssi, ssr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_left_k4 (const double *restrict sigma, int n, int ssi, int ssr,
		 double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3;
  double omega0, omega1, omega2, omega3;
  for (i = 3; i < n - 3; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      accumulator = 0.0;
      omega0 = +0.0285714285714286 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.342857142857143 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.514285714285714 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.114285714285714 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
    }
}



PyObject *
py_weights_left_k4 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_left_k4 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_left_k4 (const double *restrict f, int n, int fsi,
		     const double *restrict omega, int wsi, int wsl, int wsr,
		     double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3;
  double omega0, omega1, omega2, omega3;
  for (i = 3; i < n - 3; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      fr0 =
	+2.08333333333333 * f[(i + 0) * fsi] +
	-1.91666666666667 * f[(i + 1) * fsi] +
	+1.08333333333333 * f[(i + 2) * fsi] + -0.25 * f[(i + 3) * fsi];
      fr1 =
	+0.25 * f[(i - 1) * fsi] + +1.08333333333333 * f[(i + 0) * fsi] +
	-0.416666666666667 * f[(i + 1) * fsi] +
	+0.0833333333333333 * f[(i + 2) * fsi];
      fr2 =
	-0.0833333333333333 * f[(i - 2) * fsi] +
	+0.583333333333333 * f[(i - 1) * fsi] +
	+0.583333333333333 * f[(i + 0) * fsi] +
	-0.0833333333333333 * f[(i + 1) * fsi];
      fr3 =
	+0.0833333333333333 * f[(i - 3) * fsi] +
	-0.416666666666667 * f[(i - 2) * fsi] +
	+1.08333333333333 * f[(i - 1) * fsi] + +0.25 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3;
    }
}

PyObject *
py_reconstruct_left_k4 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_left_k4 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_right_k4 (const double *restrict sigma, int n, int ssi, int ssr,
		  double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3;
  double omega0, omega1, omega2, omega3;
  for (i = 3; i < n - 3; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      accumulator = 0.0;
      omega0 = +0.114285714285714 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.514285714285714 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.342857142857143 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.0285714285714286 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
    }
}



PyObject *
py_weights_right_k4 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_right_k4 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_right_k4 (const double *restrict f, int n, int fsi,
		      const double *restrict omega, int wsi, int wsl, int wsr,
		      double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3;
  double omega0, omega1, omega2, omega3;
  for (i = 3; i < n - 3; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      fr0 =
	+0.25 * f[(i + 0) * fsi] + +1.08333333333333 * f[(i + 1) * fsi] +
	-0.416666666666667 * f[(i + 2) * fsi] +
	+0.0833333333333333 * f[(i + 3) * fsi];
      fr1 =
	-0.0833333333333333 * f[(i - 1) * fsi] +
	+0.583333333333333 * f[(i + 0) * fsi] +
	+0.583333333333333 * f[(i + 1) * fsi] +
	-0.0833333333333333 * f[(i + 2) * fsi];
      fr2 =
	+0.0833333333333333 * f[(i - 2) * fsi] +
	-0.416666666666667 * f[(i - 1) * fsi] +
	+1.08333333333333 * f[(i + 0) * fsi] + +0.25 * f[(i + 1) * fsi];
      fr3 =
	-0.25 * f[(i - 3) * fsi] + +1.08333333333333 * f[(i - 2) * fsi] +
	-1.91666666666667 * f[(i - 1) * fsi] +
	+2.08333333333333 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3;
    }
}

PyObject *
py_reconstruct_right_k4 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_right_k4 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
smoothness_k5 (const double *restrict f, int n, int fsi,
	       double *restrict sigma, int ssi, int ssr)
{
  int i;
  double sigma0, sigma1, sigma2, sigma3, sigma4;
  for (i = 4; i < n - 4; i++)
    {
      sigma0 =
	+21.4123015873016 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	128.869246031746 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	150.560119047619 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	81.644246031746 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	17.1287698412698 * f[(i + 0) * fsi] * f[(i + 4) * fsi] +
	202.492658730159 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	488.507142857143 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	269.535317460317 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	57.144246031746 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	301.86369047619 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	338.173809523809 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	72.3934523809524 * f[(i + 2) * fsi] * f[(i + 4) * fsi] +
	95.8259920634921 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	41.369246031746 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	4.49563492063492 * f[(i + 4) * fsi] * f[(i + 4) * fsi];
      sigma1 =
	+4.49563492063492 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	27.8275793650794 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	32.7684523809524 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	17.519246031746 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	3.58710317460317 * f[(i - 1) * fsi] * f[(i + 3) * fsi] +
	48.1593253968254 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	121.42380952381 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	66.8686507936508 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	13.9359126984127 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	80.6136904761905 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	92.2571428571429 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	19.685119047619 * f[(i + 1) * fsi] * f[(i + 3) * fsi] +
	27.4926587301587 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	12.0775793650794 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	1.37063492063492 * f[(i + 3) * fsi] * f[(i + 3) * fsi];
      sigma2 =
	+1.37063492063492 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	10.119246031746 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	13.4767857142857 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	7.72757936507937 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	1.62876984126984 * f[(i - 2) * fsi] * f[(i + 2) * fsi] +
	20.8259920634921 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	59.3404761904762 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	35.5353174603175 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	7.72757936507937 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	45.8636904761905 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	59.3404761904762 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	13.4767857142857 * f[(i + 0) * fsi] * f[(i + 2) * fsi] +
	20.8259920634921 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	10.119246031746 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	1.37063492063492 * f[(i + 2) * fsi] * f[(i + 2) * fsi];
      sigma3 =
	+1.37063492063492 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	12.0775793650794 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	19.685119047619 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	13.9359126984127 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	3.58710317460317 * f[(i - 3) * fsi] * f[(i + 1) * fsi] +
	27.4926587301587 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	92.2571428571429 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	66.8686507936508 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	17.519246031746 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	80.6136904761905 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	121.42380952381 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	32.7684523809524 * f[(i - 1) * fsi] * f[(i + 1) * fsi] +
	48.1593253968254 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	27.8275793650794 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	4.49563492063492 * f[(i + 1) * fsi] * f[(i + 1) * fsi];
      sigma4 =
	+4.49563492063492 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	41.369246031746 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	72.3934523809524 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	57.144246031746 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	17.1287698412698 * f[(i - 4) * fsi] * f[(i + 0) * fsi] +
	95.8259920634921 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	338.173809523809 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	269.535317460317 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	81.644246031746 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	301.86369047619 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	488.507142857143 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	150.560119047619 * f[(i - 2) * fsi] * f[(i + 0) * fsi] +
	202.492658730159 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	128.869246031746 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	21.4123015873016 * f[(i + 0) * fsi] * f[(i + 0) * fsi];
      sigma[i * ssi + 0 * ssr] = sigma0;
      sigma[i * ssi + 1 * ssr] = sigma1;
      sigma[i * ssi + 2 * ssr] = sigma2;
      sigma[i * ssi + 3 * ssr] = sigma3;
      sigma[i * ssi + 4 * ssr] = sigma4;
    }
}


PyObject *
py_smoothness_k5 (PyObject * self, PyObject * args)
{
  double *sigma, *f;
  PyArrayObject *f_py, *sigma_py;

  long int n;
  int ssi, ssr, fsi;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &f_py, &sigma_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  f = (double *) PyArray_DATA (f_py);

  n = PyArray_DIM (f_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);

  smoothness_k5 (f, n, fsi, sigma, ssi, ssr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_left_k5 (const double *restrict sigma, int n, int ssi, int ssr,
		 double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4;
  double omega0, omega1, omega2, omega3, omega4;
  for (i = 4; i < n - 4; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      accumulator = 0.0;
      omega0 = +0.00793650793650794 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.158730158730159 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.476190476190476 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.317460317460317 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.0396825396825397 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
    }
}



PyObject *
py_weights_left_k5 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_left_k5 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_left_k5 (const double *restrict f, int n, int fsi,
		     const double *restrict omega, int wsi, int wsl, int wsr,
		     double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4;
  double omega0, omega1, omega2, omega3, omega4;
  for (i = 4; i < n - 4; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      fr0 =
	+2.28333333333333 * f[(i + 0) * fsi] +
	-2.71666666666667 * f[(i + 1) * fsi] +
	+2.28333333333333 * f[(i + 2) * fsi] + -1.05 * f[(i + 3) * fsi] +
	+0.2 * f[(i + 4) * fsi];
      fr1 =
	+0.2 * f[(i - 1) * fsi] + +1.28333333333333 * f[(i + 0) * fsi] +
	-0.716666666666667 * f[(i + 1) * fsi] +
	+0.283333333333333 * f[(i + 2) * fsi] + -0.05 * f[(i + 3) * fsi];
      fr2 =
	-0.05 * f[(i - 2) * fsi] + +0.45 * f[(i - 1) * fsi] +
	+0.783333333333333 * f[(i + 0) * fsi] +
	-0.216666666666667 * f[(i + 1) * fsi] +
	+0.0333333333333333 * f[(i + 2) * fsi];
      fr3 =
	+0.0333333333333333 * f[(i - 3) * fsi] +
	-0.216666666666667 * f[(i - 2) * fsi] +
	+0.783333333333333 * f[(i - 1) * fsi] + +0.45 * f[(i + 0) * fsi] +
	-0.05 * f[(i + 1) * fsi];
      fr4 =
	-0.05 * f[(i - 4) * fsi] + +0.283333333333333 * f[(i - 3) * fsi] +
	-0.716666666666667 * f[(i - 2) * fsi] +
	+1.28333333333333 * f[(i - 1) * fsi] + +0.2 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4;
    }
}

PyObject *
py_reconstruct_left_k5 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_left_k5 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_right_k5 (const double *restrict sigma, int n, int ssi, int ssr,
		  double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4;
  double omega0, omega1, omega2, omega3, omega4;
  for (i = 4; i < n - 4; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      accumulator = 0.0;
      omega0 = +0.0396825396825397 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.317460317460317 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.476190476190476 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.158730158730159 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.00793650793650794 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
    }
}



PyObject *
py_weights_right_k5 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_right_k5 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_right_k5 (const double *restrict f, int n, int fsi,
		      const double *restrict omega, int wsi, int wsl, int wsr,
		      double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4;
  double omega0, omega1, omega2, omega3, omega4;
  for (i = 4; i < n - 4; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      fr0 =
	+0.2 * f[(i + 0) * fsi] + +1.28333333333333 * f[(i + 1) * fsi] +
	-0.716666666666667 * f[(i + 2) * fsi] +
	+0.283333333333333 * f[(i + 3) * fsi] + -0.05 * f[(i + 4) * fsi];
      fr1 =
	-0.05 * f[(i - 1) * fsi] + +0.45 * f[(i + 0) * fsi] +
	+0.783333333333333 * f[(i + 1) * fsi] +
	-0.216666666666667 * f[(i + 2) * fsi] +
	+0.0333333333333333 * f[(i + 3) * fsi];
      fr2 =
	+0.0333333333333333 * f[(i - 2) * fsi] +
	-0.216666666666667 * f[(i - 1) * fsi] +
	+0.783333333333333 * f[(i + 0) * fsi] + +0.45 * f[(i + 1) * fsi] +
	-0.05 * f[(i + 2) * fsi];
      fr3 =
	-0.05 * f[(i - 3) * fsi] + +0.283333333333333 * f[(i - 2) * fsi] +
	-0.716666666666667 * f[(i - 1) * fsi] +
	+1.28333333333333 * f[(i + 0) * fsi] + +0.2 * f[(i + 1) * fsi];
      fr4 =
	+0.2 * f[(i - 4) * fsi] + -1.05 * f[(i - 3) * fsi] +
	+2.28333333333333 * f[(i - 2) * fsi] +
	-2.71666666666667 * f[(i - 1) * fsi] +
	+2.28333333333333 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4;
    }
}

PyObject *
py_reconstruct_right_k5 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_right_k5 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
smoothness_k6 (const double *restrict f, int n, int fsi,
	       double *restrict sigma, int ssi, int ssr)
{
  int i;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5;
  for (i = 5; i < n - 5; i++)
    {
      sigma0 =
	+50.8449983465608 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	392.364947089947 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	630.016005291005 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	524.091633597884 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	223.711722883598 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	38.9611441798942 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	784.153745039683 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	2577.47390873016 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	2173.45958994709 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	935.902678571429 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	163.974454365079 * f[(i + 1) * fsi] * f[(i + 5) * fsi] +
	2153.15287698413 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	3670.6671957672 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	1592.23273809524 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	280.413392857143 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	1577.03019179894 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	1376.16603835979 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	243.404894179894 * f[(i + 3) * fsi] * f[(i + 5) * fsi] +
	301.592981150794 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	107.061706349206 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	9.52844742063492 * f[(i + 5) * fsi] * f[(i + 5) * fsi];
      sigma1 =
	+9.52844742063492 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	75.3802248677249 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	121.878968253968 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	100.724503968254 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	42.4485284391534 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	7.2796626984127 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	160.102240410053 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	539.221593915344 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	455.140145502646 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	194.365641534392 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	33.622833994709 * f[(i + 0) * fsi] * f[(i + 4) * fsi] +
	468.437599206349 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	808.852380952381 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	350.570701058201 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	61.2508928571429 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	356.263988095238 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	313.436871693122 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	55.3456349206349 * f[(i + 2) * fsi] * f[(i + 4) * fsi] +
	69.8574487433862 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	24.9316137566138 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	2.2468501984127 * f[(i + 4) * fsi] * f[(i + 4) * fsi];
      sigma2 =
	+2.2468501984127 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	19.6825396825397 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	33.782671957672 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	28.6231150793651 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	12.059871031746 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	2.03058862433862 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	46.7370783730159 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	168.881316137566 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	148.024404761905 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	63.8887896825397 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	10.954083994709 * f[(i - 1) * fsi] * f[(i + 3) * fsi] +
	161.301025132275 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	296.11164021164 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	131.695701058201 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	23.0874669312169 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	142.159821428571 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	131.286408730159 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	23.6771164021164 * f[(i + 1) * fsi] * f[(i + 3) * fsi] +
	31.6207589285714 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	11.8218915343915 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	1.15437334656085 * f[(i + 3) * fsi] * f[(i + 3) * fsi];
      sigma3 =
	+1.15437334656085 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	11.8218915343915 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	23.6771164021164 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	23.0874669312169 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	10.954083994709 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	2.03058862433862 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	31.6207589285714 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	131.286408730159 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	131.695701058201 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	63.8887896825397 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	12.059871031746 * f[(i - 2) * fsi] * f[(i + 2) * fsi] +
	142.159821428571 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	296.11164021164 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	148.024404761905 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	28.6231150793651 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	161.301025132275 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	168.881316137566 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	33.782671957672 * f[(i + 0) * fsi] * f[(i + 2) * fsi] +
	46.7370783730159 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	19.6825396825397 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	2.2468501984127 * f[(i + 2) * fsi] * f[(i + 2) * fsi];
      sigma4 =
	+2.2468501984127 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	24.9316137566138 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	55.3456349206349 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	61.2508928571429 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	33.622833994709 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	7.2796626984127 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	69.8574487433862 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	313.436871693122 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	350.570701058201 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	194.365641534392 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	42.4485284391534 * f[(i - 3) * fsi] * f[(i + 1) * fsi] +
	356.263988095238 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	808.852380952381 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	455.140145502646 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	100.724503968254 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	468.437599206349 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	539.221593915344 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	121.878968253968 * f[(i - 1) * fsi] * f[(i + 1) * fsi] +
	160.102240410053 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	75.3802248677249 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	9.52844742063492 * f[(i + 1) * fsi] * f[(i + 1) * fsi];
      sigma5 =
	+9.52844742063492 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	107.061706349206 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	243.404894179894 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	280.413392857143 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	163.974454365079 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	38.9611441798942 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	301.592981150794 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	1376.16603835979 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	1592.23273809524 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	935.902678571429 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	223.711722883598 * f[(i - 4) * fsi] * f[(i + 0) * fsi] +
	1577.03019179894 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	3670.6671957672 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	2173.45958994709 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	524.091633597884 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	2153.15287698413 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	2577.47390873016 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	630.016005291005 * f[(i - 2) * fsi] * f[(i + 0) * fsi] +
	784.153745039683 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	392.364947089947 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	50.8449983465608 * f[(i + 0) * fsi] * f[(i + 0) * fsi];
      sigma[i * ssi + 0 * ssr] = sigma0;
      sigma[i * ssi + 1 * ssr] = sigma1;
      sigma[i * ssi + 2 * ssr] = sigma2;
      sigma[i * ssi + 3 * ssr] = sigma3;
      sigma[i * ssi + 4 * ssr] = sigma4;
      sigma[i * ssi + 5 * ssr] = sigma5;
    }
}


PyObject *
py_smoothness_k6 (PyObject * self, PyObject * args)
{
  double *sigma, *f;
  PyArrayObject *f_py, *sigma_py;

  long int n;
  int ssi, ssr, fsi;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &f_py, &sigma_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  f = (double *) PyArray_DATA (f_py);

  n = PyArray_DIM (f_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);

  smoothness_k6 (f, n, fsi, sigma, ssi, ssr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_left_k6 (const double *restrict sigma, int n, int ssi, int ssr,
		 double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5;
  double omega0, omega1, omega2, omega3, omega4, omega5;
  for (i = 5; i < n - 5; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      accumulator = 0.0;
      omega0 = +0.00216450216450216 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.0649350649350649 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.324675324675325 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.432900432900433 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.162337662337662 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega5 = +0.012987012987013 / (1e-36 + sigma5) / (1e-36 + sigma5);
      accumulator += omega5;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega5 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr] = omega5;
    }
}



PyObject *
py_weights_left_k6 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_left_k6 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_left_k6 (const double *restrict f, int n, int fsi,
		     const double *restrict omega, int wsi, int wsl, int wsr,
		     double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4, fr5;
  double omega0, omega1, omega2, omega3, omega4, omega5;
  for (i = 5; i < n - 5; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      omega5 = omega[i * wsi + 0 * wsl + 5 * wsr];
      fr0 =
	+2.45 * f[(i + 0) * fsi] + -3.55 * f[(i + 1) * fsi] +
	+3.95 * f[(i + 2) * fsi] + -2.71666666666667 * f[(i + 3) * fsi] +
	+1.03333333333333 * f[(i + 4) * fsi] +
	-0.166666666666667 * f[(i + 5) * fsi];
      fr1 =
	+0.166666666666667 * f[(i - 1) * fsi] + +1.45 * f[(i + 0) * fsi] +
	-1.05 * f[(i + 1) * fsi] + +0.616666666666667 * f[(i + 2) * fsi] +
	-0.216666666666667 * f[(i + 3) * fsi] +
	+0.0333333333333333 * f[(i + 4) * fsi];
      fr2 =
	-0.0333333333333333 * f[(i - 2) * fsi] +
	+0.366666666666667 * f[(i - 1) * fsi] + +0.95 * f[(i + 0) * fsi] +
	-0.383333333333333 * f[(i + 1) * fsi] +
	+0.116666666666667 * f[(i + 2) * fsi] +
	-0.0166666666666667 * f[(i + 3) * fsi];
      fr3 =
	+0.0166666666666667 * f[(i - 3) * fsi] +
	-0.133333333333333 * f[(i - 2) * fsi] +
	+0.616666666666667 * f[(i - 1) * fsi] +
	+0.616666666666667 * f[(i + 0) * fsi] +
	-0.133333333333333 * f[(i + 1) * fsi] +
	+0.0166666666666667 * f[(i + 2) * fsi];
      fr4 =
	-0.0166666666666667 * f[(i - 4) * fsi] +
	+0.116666666666667 * f[(i - 3) * fsi] +
	-0.383333333333333 * f[(i - 2) * fsi] + +0.95 * f[(i - 1) * fsi] +
	+0.366666666666667 * f[(i + 0) * fsi] +
	-0.0333333333333333 * f[(i + 1) * fsi];
      fr5 =
	+0.0333333333333333 * f[(i - 5) * fsi] +
	-0.216666666666667 * f[(i - 4) * fsi] +
	+0.616666666666667 * f[(i - 3) * fsi] + -1.05 * f[(i - 2) * fsi] +
	+1.45 * f[(i - 1) * fsi] + +0.166666666666667 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4 + fr5 * omega5;
    }
}

PyObject *
py_reconstruct_left_k6 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_left_k6 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_right_k6 (const double *restrict sigma, int n, int ssi, int ssr,
		  double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5;
  double omega0, omega1, omega2, omega3, omega4, omega5;
  for (i = 5; i < n - 5; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      accumulator = 0.0;
      omega0 = +0.012987012987013 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.162337662337662 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.432900432900433 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.324675324675325 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.0649350649350649 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega5 = +0.00216450216450216 / (1e-36 + sigma5) / (1e-36 + sigma5);
      accumulator += omega5;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega5 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr] = omega5;
    }
}



PyObject *
py_weights_right_k6 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_right_k6 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_right_k6 (const double *restrict f, int n, int fsi,
		      const double *restrict omega, int wsi, int wsl, int wsr,
		      double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4, fr5;
  double omega0, omega1, omega2, omega3, omega4, omega5;
  for (i = 5; i < n - 5; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      omega5 = omega[i * wsi + 0 * wsl + 5 * wsr];
      fr0 =
	+0.166666666666667 * f[(i + 0) * fsi] + +1.45 * f[(i + 1) * fsi] +
	-1.05 * f[(i + 2) * fsi] + +0.616666666666667 * f[(i + 3) * fsi] +
	-0.216666666666667 * f[(i + 4) * fsi] +
	+0.0333333333333333 * f[(i + 5) * fsi];
      fr1 =
	-0.0333333333333333 * f[(i - 1) * fsi] +
	+0.366666666666667 * f[(i + 0) * fsi] + +0.95 * f[(i + 1) * fsi] +
	-0.383333333333333 * f[(i + 2) * fsi] +
	+0.116666666666667 * f[(i + 3) * fsi] +
	-0.0166666666666667 * f[(i + 4) * fsi];
      fr2 =
	+0.0166666666666667 * f[(i - 2) * fsi] +
	-0.133333333333333 * f[(i - 1) * fsi] +
	+0.616666666666667 * f[(i + 0) * fsi] +
	+0.616666666666667 * f[(i + 1) * fsi] +
	-0.133333333333333 * f[(i + 2) * fsi] +
	+0.0166666666666667 * f[(i + 3) * fsi];
      fr3 =
	-0.0166666666666667 * f[(i - 3) * fsi] +
	+0.116666666666667 * f[(i - 2) * fsi] +
	-0.383333333333333 * f[(i - 1) * fsi] + +0.95 * f[(i + 0) * fsi] +
	+0.366666666666667 * f[(i + 1) * fsi] +
	-0.0333333333333333 * f[(i + 2) * fsi];
      fr4 =
	+0.0333333333333333 * f[(i - 4) * fsi] +
	-0.216666666666667 * f[(i - 3) * fsi] +
	+0.616666666666667 * f[(i - 2) * fsi] + -1.05 * f[(i - 1) * fsi] +
	+1.45 * f[(i + 0) * fsi] + +0.166666666666667 * f[(i + 1) * fsi];
      fr5 =
	-0.166666666666667 * f[(i - 5) * fsi] +
	+1.03333333333333 * f[(i - 4) * fsi] +
	-2.71666666666667 * f[(i - 3) * fsi] + +3.95 * f[(i - 2) * fsi] +
	-3.55 * f[(i - 1) * fsi] + +2.45 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4 + fr5 * omega5;
    }
}

PyObject *
py_reconstruct_right_k6 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_right_k6 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
smoothness_k7 (const double *restrict f, int n, int fsi,
	       double *restrict sigma, int ssi, int ssr)
{
  int i;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6;
  for (i = 6; i < n - 6; i++)
    {
      sigma0 =
	+119.876965822244 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	1140.52691383077 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	2345.30742098565 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	2648.53826913313 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	1719.4383816338 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	605.48066052389 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	90.0461092238523 * f[(i + 0) * fsi] * f[(i + 6) * fsi] +
	2787.97471636003 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	11665.8977583874 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	13315.7065438111 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	8706.93798656205 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	3081.76169462482 * f[(i + 1) * fsi] * f[(i + 5) * fsi] -
	460.055012375742 * f[(i + 1) * fsi] * f[(i + 6) * fsi] +
	12350.3314303752 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	28424.0145572791 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	18693.1184907107 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	6644.20048656205 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	995.024029781946 * f[(i + 2) * fsi] * f[(i + 6) * fsi] +
	16453.1759122308 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	21738.2182609828 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	7752.80284010742 * f[(i + 3) * fsi] * f[(i + 5) * fsi] -
	1164.09012098498 * f[(i + 3) * fsi] * f[(i + 6) * fsi] +
	7205.30018037518 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	5153.46025838745 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	775.459272837502 * f[(i + 4) * fsi] * f[(i + 6) * fsi] +
	923.494716360029 * f[(i + 5) * fsi] * f[(i + 5) * fsi] -
	278.412561978916 * f[(i + 5) * fsi] * f[(i + 6) * fsi] +
	21.0141417481695 * f[(i + 6) * fsi] * f[(i + 6) * fsi];
      sigma1 =
	+21.0141417481695 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	204.151875250521 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	422.538941047379 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	475.96589258992 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	306.899801386885 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	107.134680585618 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	15.7854224954572 * f[(i - 1) * fsi] * f[(i + 5) * fsi] +
	519.247146915584 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	2207.33120746152 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	2525.45484628026 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	1645.22305600649 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	578.412852032227 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	85.6558534251243 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	2394.05596741222 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	5559.25606962482 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	3658.67693978475 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	1295.61101896946 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	192.87048039923 * f[(i + 1) * fsi] * f[(i + 5) * fsi] +
	3266.81402951472 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	4339.66656345198 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	1547.32768578644 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	231.522065429427 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	1452.34531926407 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	1042.03954079485 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	156.661780553551 * f[(i + 3) * fsi] * f[(i + 5) * fsi] +
	187.891961730399 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	56.7392209295334 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	4.29972816792261 * f[(i + 5) * fsi] * f[(i + 5) * fsi];
      sigma2 =
	+4.29972816792261 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	44.4107718554594 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	94.9327296276255 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	108.110491355352 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	69.4589063251563 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	23.9268024991983 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	3.45697342138314 * f[(i - 2) * fsi] * f[(i + 4) * fsi] +
	121.202864508177 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	537.187110239298 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	626.822593193843 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	409.688449525012 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	142.893546476671 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	20.8355370670996 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	616.654347041847 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	1479.69665604457 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	986.137009229197 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	348.912986562049 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	51.4183199054032 * f[(i + 0) * fsi] * f[(i + 4) * fsi] +
	910.756159144354 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	1239.85097703223 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	445.834938872856 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	66.5117259232537 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	430.708745189995 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	315.141276905964 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	47.66729752886 * f[(i + 2) * fsi] * f[(i + 4) * fsi] +
	58.6280496933622 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	18.0035187690396 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	1.40409545187323 * f[(i + 4) * fsi] * f[(i + 4) * fsi];
      sigma3 =
	+1.40409545187323 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	16.2003629048421 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	38.1364719115761 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	46.8683617257228 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	31.7749557078724 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	11.3047114498156 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	1.65381755718561 * f[(i - 3) * fsi] * f[(i + 3) * fsi] +
	48.9015913600289 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	238.769633387446 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	302.017191959275 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	209.541111562049 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	75.9954446248196 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	11.3047114498156 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	302.86268037518 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	792.178909130992 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	564.852865710678 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	209.541111562049 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	31.7749557078724 * f[(i - 1) * fsi] * f[(i + 3) * fsi] +
	537.03007889744 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	792.178909130992 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	302.017191959275 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	46.8683617257228 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	302.86268037518 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	238.769633387446 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	38.1364719115761 * f[(i + 1) * fsi] * f[(i + 3) * fsi] +
	48.9015913600289 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	16.2003629048421 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	1.40409545187323 * f[(i + 3) * fsi] * f[(i + 3) * fsi];
      sigma4 =
	+1.40409545187323 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	18.0035187690396 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	47.66729752886 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	66.5117259232537 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	51.4183199054032 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	20.8355370670996 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	3.45697342138314 * f[(i - 4) * fsi] * f[(i + 2) * fsi] +
	58.6280496933622 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	315.141276905964 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	445.834938872856 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	348.912986562049 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	142.893546476671 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	23.9268024991983 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	430.708745189995 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	1239.85097703223 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	986.137009229197 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	409.688449525012 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	69.4589063251563 * f[(i - 2) * fsi] * f[(i + 2) * fsi] +
	910.756159144354 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	1479.69665604457 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	626.822593193843 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	108.110491355352 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	616.654347041847 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	537.187110239298 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	94.9327296276255 * f[(i + 0) * fsi] * f[(i + 2) * fsi] +
	121.202864508177 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	44.4107718554594 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	4.29972816792261 * f[(i + 2) * fsi] * f[(i + 2) * fsi];
      sigma5 =
	+4.29972816792261 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	56.7392209295334 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	156.661780553551 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	231.522065429427 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	192.87048039923 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	85.6558534251243 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	15.7854224954572 * f[(i - 5) * fsi] * f[(i + 1) * fsi] +
	187.891961730399 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	1042.03954079485 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	1547.32768578644 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	1295.61101896946 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	578.412852032227 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	107.134680585618 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	1452.34531926407 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	4339.66656345198 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	3658.67693978475 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	1645.22305600649 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	306.899801386885 * f[(i - 3) * fsi] * f[(i + 1) * fsi] +
	3266.81402951472 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	5559.25606962482 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	2525.45484628026 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	475.96589258992 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	2394.05596741222 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	2207.33120746152 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	422.538941047379 * f[(i - 1) * fsi] * f[(i + 1) * fsi] +
	519.247146915584 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	204.151875250521 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	21.0141417481695 * f[(i + 1) * fsi] * f[(i + 1) * fsi];
      sigma6 =
	+21.0141417481695 * f[(i - 6) * fsi] * f[(i - 6) * fsi] -
	278.412561978916 * f[(i - 6) * fsi] * f[(i - 5) * fsi] +
	775.459272837502 * f[(i - 6) * fsi] * f[(i - 4) * fsi] -
	1164.09012098498 * f[(i - 6) * fsi] * f[(i - 3) * fsi] +
	995.024029781946 * f[(i - 6) * fsi] * f[(i - 2) * fsi] -
	460.055012375742 * f[(i - 6) * fsi] * f[(i - 1) * fsi] +
	90.0461092238523 * f[(i - 6) * fsi] * f[(i + 0) * fsi] +
	923.494716360029 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	5153.46025838745 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	7752.80284010742 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	6644.20048656205 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	3081.76169462482 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	605.48066052389 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	7205.30018037518 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	21738.2182609828 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	18693.1184907107 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	8706.93798656205 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	1719.4383816338 * f[(i - 4) * fsi] * f[(i + 0) * fsi] +
	16453.1759122308 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	28424.0145572791 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	13315.7065438111 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	2648.53826913313 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	12350.3314303752 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	11665.8977583874 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	2345.30742098565 * f[(i - 2) * fsi] * f[(i + 0) * fsi] +
	2787.97471636003 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	1140.52691383077 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	119.876965822244 * f[(i + 0) * fsi] * f[(i + 0) * fsi];
      sigma[i * ssi + 0 * ssr] = sigma0;
      sigma[i * ssi + 1 * ssr] = sigma1;
      sigma[i * ssi + 2 * ssr] = sigma2;
      sigma[i * ssi + 3 * ssr] = sigma3;
      sigma[i * ssi + 4 * ssr] = sigma4;
      sigma[i * ssi + 5 * ssr] = sigma5;
      sigma[i * ssi + 6 * ssr] = sigma6;
    }
}


PyObject *
py_smoothness_k7 (PyObject * self, PyObject * args)
{
  double *sigma, *f;
  PyArrayObject *f_py, *sigma_py;

  long int n;
  int ssi, ssr, fsi;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &f_py, &sigma_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  f = (double *) PyArray_DATA (f_py);

  n = PyArray_DIM (f_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);

  smoothness_k7 (f, n, fsi, sigma, ssi, ssr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_left_k7 (const double *restrict sigma, int n, int ssi, int ssr,
		 double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6;
  for (i = 6; i < n - 6; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      sigma6 = sigma[i * ssi + 6 * ssr];
      accumulator = 0.0;
      omega0 = +0.000582750582750583 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.0244755244755245 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.183566433566434 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.407925407925408 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.305944055944056 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega5 = +0.0734265734265734 / (1e-36 + sigma5) / (1e-36 + sigma5);
      accumulator += omega5;
      omega6 = +0.00407925407925408 / (1e-36 + sigma6) / (1e-36 + sigma6);
      accumulator += omega6;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega5 /= accumulator;
      omega6 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr] = omega5;
      omega[i * wsi + 0 * wsl + 6 * wsr] = omega6;
    }
}



PyObject *
py_weights_left_k7 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_left_k7 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_left_k7 (const double *restrict f, int n, int fsi,
		     const double *restrict omega, int wsi, int wsl, int wsr,
		     double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4, fr5, fr6;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6;
  for (i = 6; i < n - 6; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      omega5 = omega[i * wsi + 0 * wsl + 5 * wsr];
      omega6 = omega[i * wsi + 0 * wsl + 6 * wsr];
      fr0 =
	+2.59285714285714 * f[(i + 0) * fsi] +
	-4.40714285714286 * f[(i + 1) * fsi] +
	+6.09285714285714 * f[(i + 2) * fsi] +
	-5.57380952380952 * f[(i + 3) * fsi] +
	+3.17619047619048 * f[(i + 4) * fsi] +
	-1.02380952380952 * f[(i + 5) * fsi] +
	+0.142857142857143 * f[(i + 6) * fsi];
      fr1 =
	+0.142857142857143 * f[(i - 1) * fsi] +
	+1.59285714285714 * f[(i + 0) * fsi] +
	-1.40714285714286 * f[(i + 1) * fsi] +
	+1.09285714285714 * f[(i + 2) * fsi] +
	-0.573809523809524 * f[(i + 3) * fsi] +
	+0.176190476190476 * f[(i + 4) * fsi] +
	-0.0238095238095238 * f[(i + 5) * fsi];
      fr2 =
	-0.0238095238095238 * f[(i - 2) * fsi] +
	+0.30952380952381 * f[(i - 1) * fsi] +
	+1.09285714285714 * f[(i + 0) * fsi] +
	-0.573809523809524 * f[(i + 1) * fsi] +
	+0.25952380952381 * f[(i + 2) * fsi] +
	-0.0738095238095238 * f[(i + 3) * fsi] +
	+0.00952380952380952 * f[(i + 4) * fsi];
      fr3 =
	+0.00952380952380952 * f[(i - 3) * fsi] +
	-0.0904761904761905 * f[(i - 2) * fsi] +
	+0.509523809523809 * f[(i - 1) * fsi] +
	+0.759523809523809 * f[(i + 0) * fsi] +
	-0.24047619047619 * f[(i + 1) * fsi] +
	+0.0595238095238095 * f[(i + 2) * fsi] +
	-0.00714285714285714 * f[(i + 3) * fsi];
      fr4 =
	-0.00714285714285714 * f[(i - 4) * fsi] +
	+0.0595238095238095 * f[(i - 3) * fsi] +
	-0.24047619047619 * f[(i - 2) * fsi] +
	+0.759523809523809 * f[(i - 1) * fsi] +
	+0.509523809523809 * f[(i + 0) * fsi] +
	-0.0904761904761905 * f[(i + 1) * fsi] +
	+0.00952380952380952 * f[(i + 2) * fsi];
      fr5 =
	+0.00952380952380952 * f[(i - 5) * fsi] +
	-0.0738095238095238 * f[(i - 4) * fsi] +
	+0.25952380952381 * f[(i - 3) * fsi] +
	-0.573809523809524 * f[(i - 2) * fsi] +
	+1.09285714285714 * f[(i - 1) * fsi] +
	+0.30952380952381 * f[(i + 0) * fsi] +
	-0.0238095238095238 * f[(i + 1) * fsi];
      fr6 =
	-0.0238095238095238 * f[(i - 6) * fsi] +
	+0.176190476190476 * f[(i - 5) * fsi] +
	-0.573809523809524 * f[(i - 4) * fsi] +
	+1.09285714285714 * f[(i - 3) * fsi] +
	-1.40714285714286 * f[(i - 2) * fsi] +
	+1.59285714285714 * f[(i - 1) * fsi] +
	+0.142857142857143 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4 + fr5 * omega5 + fr6 * omega6;
    }
}

PyObject *
py_reconstruct_left_k7 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_left_k7 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_right_k7 (const double *restrict sigma, int n, int ssi, int ssr,
		  double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6;
  for (i = 6; i < n - 6; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      sigma6 = sigma[i * ssi + 6 * ssr];
      accumulator = 0.0;
      omega0 = +0.00407925407925408 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.0734265734265734 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.305944055944056 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.407925407925408 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.183566433566434 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega5 = +0.0244755244755245 / (1e-36 + sigma5) / (1e-36 + sigma5);
      accumulator += omega5;
      omega6 = +0.000582750582750583 / (1e-36 + sigma6) / (1e-36 + sigma6);
      accumulator += omega6;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega5 /= accumulator;
      omega6 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr] = omega5;
      omega[i * wsi + 0 * wsl + 6 * wsr] = omega6;
    }
}



PyObject *
py_weights_right_k7 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_right_k7 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_right_k7 (const double *restrict f, int n, int fsi,
		      const double *restrict omega, int wsi, int wsl, int wsr,
		      double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4, fr5, fr6;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6;
  for (i = 6; i < n - 6; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      omega5 = omega[i * wsi + 0 * wsl + 5 * wsr];
      omega6 = omega[i * wsi + 0 * wsl + 6 * wsr];
      fr0 =
	+0.142857142857143 * f[(i + 0) * fsi] +
	+1.59285714285714 * f[(i + 1) * fsi] +
	-1.40714285714286 * f[(i + 2) * fsi] +
	+1.09285714285714 * f[(i + 3) * fsi] +
	-0.573809523809524 * f[(i + 4) * fsi] +
	+0.176190476190476 * f[(i + 5) * fsi] +
	-0.0238095238095238 * f[(i + 6) * fsi];
      fr1 =
	-0.0238095238095238 * f[(i - 1) * fsi] +
	+0.30952380952381 * f[(i + 0) * fsi] +
	+1.09285714285714 * f[(i + 1) * fsi] +
	-0.573809523809524 * f[(i + 2) * fsi] +
	+0.25952380952381 * f[(i + 3) * fsi] +
	-0.0738095238095238 * f[(i + 4) * fsi] +
	+0.00952380952380952 * f[(i + 5) * fsi];
      fr2 =
	+0.00952380952380952 * f[(i - 2) * fsi] +
	-0.0904761904761905 * f[(i - 1) * fsi] +
	+0.509523809523809 * f[(i + 0) * fsi] +
	+0.759523809523809 * f[(i + 1) * fsi] +
	-0.24047619047619 * f[(i + 2) * fsi] +
	+0.0595238095238095 * f[(i + 3) * fsi] +
	-0.00714285714285714 * f[(i + 4) * fsi];
      fr3 =
	-0.00714285714285714 * f[(i - 3) * fsi] +
	+0.0595238095238095 * f[(i - 2) * fsi] +
	-0.24047619047619 * f[(i - 1) * fsi] +
	+0.759523809523809 * f[(i + 0) * fsi] +
	+0.509523809523809 * f[(i + 1) * fsi] +
	-0.0904761904761905 * f[(i + 2) * fsi] +
	+0.00952380952380952 * f[(i + 3) * fsi];
      fr4 =
	+0.00952380952380952 * f[(i - 4) * fsi] +
	-0.0738095238095238 * f[(i - 3) * fsi] +
	+0.25952380952381 * f[(i - 2) * fsi] +
	-0.573809523809524 * f[(i - 1) * fsi] +
	+1.09285714285714 * f[(i + 0) * fsi] +
	+0.30952380952381 * f[(i + 1) * fsi] +
	-0.0238095238095238 * f[(i + 2) * fsi];
      fr5 =
	-0.0238095238095238 * f[(i - 5) * fsi] +
	+0.176190476190476 * f[(i - 4) * fsi] +
	-0.573809523809524 * f[(i - 3) * fsi] +
	+1.09285714285714 * f[(i - 2) * fsi] +
	-1.40714285714286 * f[(i - 1) * fsi] +
	+1.59285714285714 * f[(i + 0) * fsi] +
	+0.142857142857143 * f[(i + 1) * fsi];
      fr6 =
	+0.142857142857143 * f[(i - 6) * fsi] +
	-1.02380952380952 * f[(i - 5) * fsi] +
	+3.17619047619048 * f[(i - 4) * fsi] +
	-5.57380952380952 * f[(i - 3) * fsi] +
	+6.09285714285714 * f[(i - 2) * fsi] +
	-4.40714285714286 * f[(i - 1) * fsi] +
	+2.59285714285714 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4 + fr5 * omega5 + fr6 * omega6;
    }
}

PyObject *
py_reconstruct_right_k7 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_right_k7 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
smoothness_k8 (const double *restrict f, int n, int fsi,
	       double *restrict sigma, int ssi, int ssr)
{
  int i;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7;
  for (i = 7; i < n - 7; i++)
    {
      sigma0 =
	+282.837600612977 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	3217.68658751845 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	8098.17149379591 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	11602.0447459035 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	10159.0931716562 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	5415.63416825253 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	1622.88916698582 * f[(i + 0) * fsi] * f[(i + 6) * fsi] -
	210.463531989358 * f[(i + 0) * fsi] * f[(i + 7) * fsi] +
	9343.02132742786 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	47645.872787025 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	68840.1294128134 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	60634.3990483282 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	32462.762767691 * f[(i + 1) * fsi] * f[(i + 5) * fsi] -
	9759.93192303144 * f[(i + 1) * fsi] * f[(i + 6) * fsi] +
	1268.95551054292 * f[(i + 1) * fsi] * f[(i + 7) * fsi] +
	61294.8370166772 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	178245.759975439 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	157723.978487162 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	84736.2897924522 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	25544.3501239796 * f[(i + 2) * fsi] * f[(i + 6) * fsi] -
	3328.25158337599 * f[(i + 2) * fsi] * f[(i + 7) * fsi] +
	130199.124980546 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	231245.307361434 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	124579.67848041 * f[(i + 3) * fsi] * f[(i + 5) * fsi] -
	37637.4314325875 * f[(i + 3) * fsi] * f[(i + 6) * fsi] +
	4912.48566104661 * f[(i + 3) * fsi] * f[(i + 7) * fsi] +
	102966.44021251 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	111189.450476982 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	33651.8387772039 * f[(i + 4) * fsi] * f[(i + 6) * fsi] -
	4398.63397429859 * f[(i + 4) * fsi] * f[(i + 7) * fsi] +
	30071.0784359481 * f[(i + 5) * fsi] * f[(i + 5) * fsi] -
	18228.7647006052 * f[(i + 5) * fsi] * f[(i + 6) * fsi] +
	2385.54101829437 * f[(i + 5) * fsi] * f[(i + 7) * fsi] +
	2765.84444133604 * f[(i + 6) * fsi] * f[(i + 6) * fsi] -
	724.638894617214 * f[(i + 6) * fsi] * f[(i + 7) * fsi] +
	47.5028971986251 * f[(i + 7) * fsi] * f[(i + 7) * fsi];
      sigma1 =
	+47.5028971986251 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	549.582823188643 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	1391.20673258008 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	1992.07290287002 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	1737.9199467609 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	921.690511947415 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	274.621224828637 * f[(i - 1) * fsi] * f[(i + 5) * fsi] -
	35.407460560787 * f[(i - 1) * fsi] * f[(i + 6) * fsi] +
	1639.31476541011 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	8454.36155245706 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	12248.796925352 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	10772.9570807356 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	5746.65947583143 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	1719.62507117959 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	222.440595557253 * f[(i + 0) * fsi] * f[(i + 6) * fsi] +
	11054.5384359481 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	32362.4054769819 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	28675.0021841139 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	15380.2247924521 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	4621.402767691 * f[(i + 1) * fsi] * f[(i + 5) * fsi] -
	599.696734390095 * f[(i + 1) * fsi] * f[(i + 6) * fsi] +
	23881.8339625102 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	42591.6661577299 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	22956.5584871618 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	6924.03404832822 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	901.155248375756 * f[(i + 2) * fsi] * f[(i + 6) * fsi] +
	19089.3249805465 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	20664.4461791424 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	6253.56570910972 * f[(i + 3) * fsi] * f[(i + 5) * fsi] -
	816.068383469668 * f[(i + 3) * fsi] * f[(i + 6) * fsi] +
	5612.02326667724 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	3406.48778702496 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	445.58477421919 * f[(i + 4) * fsi] * f[(i + 6) * fsi] +
	518.201327427861 * f[(i + 5) * fsi] * f[(i + 5) * fsi] -
	135.845449952311 * f[(i + 5) * fsi] * f[(i + 6) * fsi] +
	8.91870511033141 * f[(i + 6) * fsi] * f[(i + 6) * fsi];
      sigma2 =
	+8.91870511033141 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	107.291821204516 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	277.006890621306 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	399.198237967023 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	347.463467070642 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	182.82658888745 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	53.8627119593691 * f[(i - 2) * fsi] * f[(i + 4) * fsi] -
	6.85383181299153 * f[(i - 2) * fsi] * f[(i + 5) * fsi] +
	335.04033977354 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	1774.22905245706 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	2601.97484491219 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	2293.25840018007 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	1217.71486645643 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	361.183311920333 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	46.1921948462738 * f[(i - 1) * fsi] * f[(i + 5) * fsi] +
	2403.24289630687 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	7175.23886432754 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	6406.93231432228 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	3435.42821837807 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	1026.47873509069 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	132.007597485334 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	5440.58053610203 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	9841.58822563111 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	5330.2740359658 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	1605.02809925414 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	207.643474097758 * f[(i + 1) * fsi] * f[(i + 5) * fsi] +
	4502.6216168312 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	4924.83347080903 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	1494.60136979645 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	194.560288231573 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	1358.55473224437 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	830.843311716319 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	108.833222879904 * f[(i + 3) * fsi] * f[(i + 5) * fsi] +
	127.914395039744 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	33.7168840352035 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	2.23485487058274 * f[(i + 5) * fsi] * f[(i + 5) * fsi];
      sigma3 =
	+2.23485487058274 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	28.9038461163322 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	78.9596779063593 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	118.296148019933 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	105.236207783825 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	55.7434572736934 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	16.3186498727289 * f[(i - 3) * fsi] * f[(i + 3) * fsi] -
	2.04079389412028 * f[(i - 3) * fsi] * f[(i + 4) * fsi] +
	97.1187623236942 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	547.061953691627 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	839.561493253242 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	761.319673328215 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	409.596543732663 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	121.468497105518 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	15.3584086083991 * f[(i - 2) * fsi] * f[(i + 4) * fsi] +
	793.785102614737 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	2499.75828562384 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	2315.13502362012 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	1267.31229245215 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	381.255607197169 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	48.78798218551 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	2019.32231127564 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	3827.93467624839 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	2136.14046247043 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	653.059881661548 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	84.7024132787544 * f[(i + 0) * fsi] * f[(i + 4) * fsi] +
	1856.32621511438 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	2115.5956853152 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	658.5622523196 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	86.7358790604971 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	615.75035001057 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	390.9897931978 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	52.4035220146045 * f[(i + 2) * fsi] * f[(i + 4) * fsi] +
	63.3507101439102 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	17.3197577124522 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	1.21003447541078 * f[(i + 4) * fsi] * f[(i + 4) * fsi];
      sigma4 =
	+1.21003447541078 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	17.3197577124522 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	52.4035220146045 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	86.7358790604971 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	84.7024132787544 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	48.78798218551 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	15.3584086083991 * f[(i - 4) * fsi] * f[(i + 2) * fsi] -
	2.04079389412028 * f[(i - 4) * fsi] * f[(i + 3) * fsi] +
	63.3507101439102 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	390.9897931978 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	658.5622523196 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	653.059881661548 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	381.255607197169 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	121.468497105518 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	16.3186498727289 * f[(i - 3) * fsi] * f[(i + 3) * fsi] +
	615.75035001057 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	2115.5956853152 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	2136.14046247043 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	1267.31229245215 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	409.596543732663 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	55.7434572736934 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	1856.32621511438 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	3827.93467624839 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	2315.13502362012 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	761.319673328215 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	105.236207783825 * f[(i - 1) * fsi] * f[(i + 3) * fsi] +
	2019.32231127564 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	2499.75828562384 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	839.561493253242 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	118.296148019933 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	793.785102614737 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	547.061953691627 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	78.9596779063593 * f[(i + 1) * fsi] * f[(i + 3) * fsi] +
	97.1187623236942 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	28.9038461163322 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	2.23485487058274 * f[(i + 3) * fsi] * f[(i + 3) * fsi];
      sigma5 =
	+2.23485487058274 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	33.7168840352035 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	108.833222879904 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	194.560288231573 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	207.643474097758 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	132.007597485334 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	46.1921948462738 * f[(i - 5) * fsi] * f[(i + 1) * fsi] -
	6.85383181299153 * f[(i - 5) * fsi] * f[(i + 2) * fsi] +
	127.914395039744 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	830.843311716319 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	1494.60136979645 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	1605.02809925414 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	1026.47873509069 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	361.183311920333 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	53.8627119593691 * f[(i - 4) * fsi] * f[(i + 2) * fsi] +
	1358.55473224437 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	4924.83347080903 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	5330.2740359658 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	3435.42821837807 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	1217.71486645643 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	182.82658888745 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	4502.6216168312 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	9841.58822563111 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	6406.93231432228 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	2293.25840018007 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	347.463467070642 * f[(i - 2) * fsi] * f[(i + 2) * fsi] +
	5440.58053610203 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	7175.23886432754 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	2601.97484491219 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	399.198237967023 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	2403.24289630687 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	1774.22905245706 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	277.006890621306 * f[(i + 0) * fsi] * f[(i + 2) * fsi] +
	335.04033977354 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	107.291821204516 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	8.91870511033141 * f[(i + 2) * fsi] * f[(i + 2) * fsi];
      sigma6 =
	+8.91870511033141 * f[(i - 6) * fsi] * f[(i - 6) * fsi] -
	135.845449952311 * f[(i - 6) * fsi] * f[(i - 5) * fsi] +
	445.58477421919 * f[(i - 6) * fsi] * f[(i - 4) * fsi] -
	816.068383469668 * f[(i - 6) * fsi] * f[(i - 3) * fsi] +
	901.155248375756 * f[(i - 6) * fsi] * f[(i - 2) * fsi] -
	599.696734390095 * f[(i - 6) * fsi] * f[(i - 1) * fsi] +
	222.440595557253 * f[(i - 6) * fsi] * f[(i + 0) * fsi] -
	35.407460560787 * f[(i - 6) * fsi] * f[(i + 1) * fsi] +
	518.201327427861 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	3406.48778702496 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	6253.56570910972 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	6924.03404832822 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	4621.402767691 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	1719.62507117959 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	274.621224828637 * f[(i - 5) * fsi] * f[(i + 1) * fsi] +
	5612.02326667724 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	20664.4461791424 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	22956.5584871618 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	15380.2247924521 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	5746.65947583143 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	921.690511947415 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	19089.3249805465 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	42591.6661577299 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	28675.0021841139 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	10772.9570807356 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	1737.9199467609 * f[(i - 3) * fsi] * f[(i + 1) * fsi] +
	23881.8339625102 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	32362.4054769819 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	12248.796925352 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	1992.07290287002 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	11054.5384359481 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	8454.36155245706 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	1391.20673258008 * f[(i - 1) * fsi] * f[(i + 1) * fsi] +
	1639.31476541011 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	549.582823188643 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	47.5028971986251 * f[(i + 1) * fsi] * f[(i + 1) * fsi];
      sigma7 =
	+47.5028971986251 * f[(i - 7) * fsi] * f[(i - 7) * fsi] -
	724.638894617214 * f[(i - 7) * fsi] * f[(i - 6) * fsi] +
	2385.54101829437 * f[(i - 7) * fsi] * f[(i - 5) * fsi] -
	4398.63397429859 * f[(i - 7) * fsi] * f[(i - 4) * fsi] +
	4912.48566104661 * f[(i - 7) * fsi] * f[(i - 3) * fsi] -
	3328.25158337599 * f[(i - 7) * fsi] * f[(i - 2) * fsi] +
	1268.95551054292 * f[(i - 7) * fsi] * f[(i - 1) * fsi] -
	210.463531989358 * f[(i - 7) * fsi] * f[(i + 0) * fsi] +
	2765.84444133604 * f[(i - 6) * fsi] * f[(i - 6) * fsi] -
	18228.7647006052 * f[(i - 6) * fsi] * f[(i - 5) * fsi] +
	33651.8387772039 * f[(i - 6) * fsi] * f[(i - 4) * fsi] -
	37637.4314325875 * f[(i - 6) * fsi] * f[(i - 3) * fsi] +
	25544.3501239796 * f[(i - 6) * fsi] * f[(i - 2) * fsi] -
	9759.93192303144 * f[(i - 6) * fsi] * f[(i - 1) * fsi] +
	1622.88916698582 * f[(i - 6) * fsi] * f[(i + 0) * fsi] +
	30071.0784359481 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	111189.450476982 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	124579.67848041 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	84736.2897924522 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	32462.762767691 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	5415.63416825253 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	102966.44021251 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	231245.307361434 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	157723.978487162 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	60634.3990483282 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	10159.0931716562 * f[(i - 4) * fsi] * f[(i + 0) * fsi] +
	130199.124980546 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	178245.759975439 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	68840.1294128134 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	11602.0447459035 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	61294.8370166772 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	47645.872787025 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	8098.17149379591 * f[(i - 2) * fsi] * f[(i + 0) * fsi] +
	9343.02132742786 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	3217.68658751845 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	282.837600612977 * f[(i + 0) * fsi] * f[(i + 0) * fsi];
      sigma[i * ssi + 0 * ssr] = sigma0;
      sigma[i * ssi + 1 * ssr] = sigma1;
      sigma[i * ssi + 2 * ssr] = sigma2;
      sigma[i * ssi + 3 * ssr] = sigma3;
      sigma[i * ssi + 4 * ssr] = sigma4;
      sigma[i * ssi + 5 * ssr] = sigma5;
      sigma[i * ssi + 6 * ssr] = sigma6;
      sigma[i * ssi + 7 * ssr] = sigma7;
    }
}


PyObject *
py_smoothness_k8 (PyObject * self, PyObject * args)
{
  double *sigma, *f;
  PyArrayObject *f_py, *sigma_py;

  long int n;
  int ssi, ssr, fsi;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &f_py, &sigma_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  f = (double *) PyArray_DATA (f_py);

  n = PyArray_DIM (f_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);

  smoothness_k8 (f, n, fsi, sigma, ssi, ssr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_left_k8 (const double *restrict sigma, int n, int ssi, int ssr,
		 double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7;
  for (i = 7; i < n - 7; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      sigma6 = sigma[i * ssi + 6 * ssr];
      sigma7 = sigma[i * ssi + 7 * ssr];
      accumulator = 0.0;
      omega0 = +0.000155400155400155 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.0087024087024087 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.0913752913752914 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.304584304584305 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.380730380730381 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega5 = +0.182750582750583 / (1e-36 + sigma5) / (1e-36 + sigma5);
      accumulator += omega5;
      omega6 = +0.0304584304584305 / (1e-36 + sigma6) / (1e-36 + sigma6);
      accumulator += omega6;
      omega7 = +0.00124320124320124 / (1e-36 + sigma7) / (1e-36 + sigma7);
      accumulator += omega7;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega5 /= accumulator;
      omega6 /= accumulator;
      omega7 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr] = omega5;
      omega[i * wsi + 0 * wsl + 6 * wsr] = omega6;
      omega[i * wsi + 0 * wsl + 7 * wsr] = omega7;
    }
}



PyObject *
py_weights_left_k8 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_left_k8 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_left_k8 (const double *restrict f, int n, int fsi,
		     const double *restrict omega, int wsi, int wsl, int wsr,
		     double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4, fr5, fr6, fr7;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7;
  for (i = 7; i < n - 7; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      omega5 = omega[i * wsi + 0 * wsl + 5 * wsr];
      omega6 = omega[i * wsi + 0 * wsl + 6 * wsr];
      omega7 = omega[i * wsi + 0 * wsl + 7 * wsr];
      fr0 =
	+2.71785714285714 * f[(i + 0) * fsi] +
	-5.28214285714286 * f[(i + 1) * fsi] +
	+8.71785714285714 * f[(i + 2) * fsi] +
	-9.94880952380952 * f[(i + 3) * fsi] +
	+7.55119047619048 * f[(i + 4) * fsi] +
	-3.64880952380952 * f[(i + 5) * fsi] +
	+1.01785714285714 * f[(i + 6) * fsi] + -0.125 * f[(i + 7) * fsi];
      fr1 =
	+0.125 * f[(i - 1) * fsi] + +1.71785714285714 * f[(i + 0) * fsi] +
	-1.78214285714286 * f[(i + 1) * fsi] +
	+1.71785714285714 * f[(i + 2) * fsi] +
	-1.19880952380952 * f[(i + 3) * fsi] +
	+0.551190476190476 * f[(i + 4) * fsi] +
	-0.148809523809524 * f[(i + 5) * fsi] +
	+0.0178571428571429 * f[(i + 6) * fsi];
      fr2 =
	-0.0178571428571429 * f[(i - 2) * fsi] +
	+0.267857142857143 * f[(i - 1) * fsi] +
	+1.21785714285714 * f[(i + 0) * fsi] +
	-0.782142857142857 * f[(i + 1) * fsi] +
	+0.467857142857143 * f[(i + 2) * fsi] +
	-0.198809523809524 * f[(i + 3) * fsi] +
	+0.0511904761904762 * f[(i + 4) * fsi] +
	-0.00595238095238095 * f[(i + 5) * fsi];
      fr3 =
	+0.00595238095238095 * f[(i - 3) * fsi] +
	-0.0654761904761905 * f[(i - 2) * fsi] +
	+0.43452380952381 * f[(i - 1) * fsi] +
	+0.884523809523809 * f[(i + 0) * fsi] +
	-0.36547619047619 * f[(i + 1) * fsi] +
	+0.13452380952381 * f[(i + 2) * fsi] +
	-0.0321428571428571 * f[(i + 3) * fsi] +
	+0.00357142857142857 * f[(i + 4) * fsi];
      fr4 =
	-0.00357142857142857 * f[(i - 4) * fsi] +
	+0.0345238095238095 * f[(i - 3) * fsi] +
	-0.16547619047619 * f[(i - 2) * fsi] +
	+0.634523809523809 * f[(i - 1) * fsi] +
	+0.634523809523809 * f[(i + 0) * fsi] +
	-0.16547619047619 * f[(i + 1) * fsi] +
	+0.0345238095238095 * f[(i + 2) * fsi] +
	-0.00357142857142857 * f[(i + 3) * fsi];
      fr5 =
	+0.00357142857142857 * f[(i - 5) * fsi] +
	-0.0321428571428571 * f[(i - 4) * fsi] +
	+0.13452380952381 * f[(i - 3) * fsi] +
	-0.36547619047619 * f[(i - 2) * fsi] +
	+0.884523809523809 * f[(i - 1) * fsi] +
	+0.43452380952381 * f[(i + 0) * fsi] +
	-0.0654761904761905 * f[(i + 1) * fsi] +
	+0.00595238095238095 * f[(i + 2) * fsi];
      fr6 =
	-0.00595238095238095 * f[(i - 6) * fsi] +
	+0.0511904761904762 * f[(i - 5) * fsi] +
	-0.198809523809524 * f[(i - 4) * fsi] +
	+0.467857142857143 * f[(i - 3) * fsi] +
	-0.782142857142857 * f[(i - 2) * fsi] +
	+1.21785714285714 * f[(i - 1) * fsi] +
	+0.267857142857143 * f[(i + 0) * fsi] +
	-0.0178571428571429 * f[(i + 1) * fsi];
      fr7 =
	+0.0178571428571429 * f[(i - 7) * fsi] +
	-0.148809523809524 * f[(i - 6) * fsi] +
	+0.551190476190476 * f[(i - 5) * fsi] +
	-1.19880952380952 * f[(i - 4) * fsi] +
	+1.71785714285714 * f[(i - 3) * fsi] +
	-1.78214285714286 * f[(i - 2) * fsi] +
	+1.71785714285714 * f[(i - 1) * fsi] + +0.125 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4 + fr5 * omega5 + fr6 * omega6 + fr7 * omega7;
    }
}

PyObject *
py_reconstruct_left_k8 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_left_k8 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_right_k8 (const double *restrict sigma, int n, int ssi, int ssr,
		  double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7;
  for (i = 7; i < n - 7; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      sigma6 = sigma[i * ssi + 6 * ssr];
      sigma7 = sigma[i * ssi + 7 * ssr];
      accumulator = 0.0;
      omega0 = +0.00124320124320124 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.0304584304584305 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.182750582750583 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.380730380730381 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.304584304584305 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega5 = +0.0913752913752914 / (1e-36 + sigma5) / (1e-36 + sigma5);
      accumulator += omega5;
      omega6 = +0.0087024087024087 / (1e-36 + sigma6) / (1e-36 + sigma6);
      accumulator += omega6;
      omega7 = +0.000155400155400155 / (1e-36 + sigma7) / (1e-36 + sigma7);
      accumulator += omega7;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega5 /= accumulator;
      omega6 /= accumulator;
      omega7 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr] = omega5;
      omega[i * wsi + 0 * wsl + 6 * wsr] = omega6;
      omega[i * wsi + 0 * wsl + 7 * wsr] = omega7;
    }
}



PyObject *
py_weights_right_k8 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_right_k8 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_right_k8 (const double *restrict f, int n, int fsi,
		      const double *restrict omega, int wsi, int wsl, int wsr,
		      double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4, fr5, fr6, fr7;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7;
  for (i = 7; i < n - 7; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      omega5 = omega[i * wsi + 0 * wsl + 5 * wsr];
      omega6 = omega[i * wsi + 0 * wsl + 6 * wsr];
      omega7 = omega[i * wsi + 0 * wsl + 7 * wsr];
      fr0 =
	+0.125 * f[(i + 0) * fsi] + +1.71785714285714 * f[(i + 1) * fsi] +
	-1.78214285714286 * f[(i + 2) * fsi] +
	+1.71785714285714 * f[(i + 3) * fsi] +
	-1.19880952380952 * f[(i + 4) * fsi] +
	+0.551190476190476 * f[(i + 5) * fsi] +
	-0.148809523809524 * f[(i + 6) * fsi] +
	+0.0178571428571429 * f[(i + 7) * fsi];
      fr1 =
	-0.0178571428571429 * f[(i - 1) * fsi] +
	+0.267857142857143 * f[(i + 0) * fsi] +
	+1.21785714285714 * f[(i + 1) * fsi] +
	-0.782142857142857 * f[(i + 2) * fsi] +
	+0.467857142857143 * f[(i + 3) * fsi] +
	-0.198809523809524 * f[(i + 4) * fsi] +
	+0.0511904761904762 * f[(i + 5) * fsi] +
	-0.00595238095238095 * f[(i + 6) * fsi];
      fr2 =
	+0.00595238095238095 * f[(i - 2) * fsi] +
	-0.0654761904761905 * f[(i - 1) * fsi] +
	+0.43452380952381 * f[(i + 0) * fsi] +
	+0.884523809523809 * f[(i + 1) * fsi] +
	-0.36547619047619 * f[(i + 2) * fsi] +
	+0.13452380952381 * f[(i + 3) * fsi] +
	-0.0321428571428571 * f[(i + 4) * fsi] +
	+0.00357142857142857 * f[(i + 5) * fsi];
      fr3 =
	-0.00357142857142857 * f[(i - 3) * fsi] +
	+0.0345238095238095 * f[(i - 2) * fsi] +
	-0.16547619047619 * f[(i - 1) * fsi] +
	+0.634523809523809 * f[(i + 0) * fsi] +
	+0.634523809523809 * f[(i + 1) * fsi] +
	-0.16547619047619 * f[(i + 2) * fsi] +
	+0.0345238095238095 * f[(i + 3) * fsi] +
	-0.00357142857142857 * f[(i + 4) * fsi];
      fr4 =
	+0.00357142857142857 * f[(i - 4) * fsi] +
	-0.0321428571428571 * f[(i - 3) * fsi] +
	+0.13452380952381 * f[(i - 2) * fsi] +
	-0.36547619047619 * f[(i - 1) * fsi] +
	+0.884523809523809 * f[(i + 0) * fsi] +
	+0.43452380952381 * f[(i + 1) * fsi] +
	-0.0654761904761905 * f[(i + 2) * fsi] +
	+0.00595238095238095 * f[(i + 3) * fsi];
      fr5 =
	-0.00595238095238095 * f[(i - 5) * fsi] +
	+0.0511904761904762 * f[(i - 4) * fsi] +
	-0.198809523809524 * f[(i - 3) * fsi] +
	+0.467857142857143 * f[(i - 2) * fsi] +
	-0.782142857142857 * f[(i - 1) * fsi] +
	+1.21785714285714 * f[(i + 0) * fsi] +
	+0.267857142857143 * f[(i + 1) * fsi] +
	-0.0178571428571429 * f[(i + 2) * fsi];
      fr6 =
	+0.0178571428571429 * f[(i - 6) * fsi] +
	-0.148809523809524 * f[(i - 5) * fsi] +
	+0.551190476190476 * f[(i - 4) * fsi] +
	-1.19880952380952 * f[(i - 3) * fsi] +
	+1.71785714285714 * f[(i - 2) * fsi] +
	-1.78214285714286 * f[(i - 1) * fsi] +
	+1.71785714285714 * f[(i + 0) * fsi] + +0.125 * f[(i + 1) * fsi];
      fr7 =
	-0.125 * f[(i - 7) * fsi] + +1.01785714285714 * f[(i - 6) * fsi] +
	-3.64880952380952 * f[(i - 5) * fsi] +
	+7.55119047619048 * f[(i - 4) * fsi] +
	-9.94880952380952 * f[(i - 3) * fsi] +
	+8.71785714285714 * f[(i - 2) * fsi] +
	-5.28214285714286 * f[(i - 1) * fsi] +
	+2.71785714285714 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4 + fr5 * omega5 + fr6 * omega6 + fr7 * omega7;
    }
}

PyObject *
py_reconstruct_right_k8 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_right_k8 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
smoothness_k9 (const double *restrict f, int n, int fsi,
	       double *restrict sigma, int ssi, int ssr)
{
  int i;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7,
    sigma8;
  for (i = 8; i < n - 8; i++)
    {
      sigma0 =
	+669.714981108808 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	8893.78045641284 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	26542.9748247279 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	46202.7564018809 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	51067.48172473 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	36585.0576509742 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	16552.0014823105 * f[(i + 0) * fsi] * f[(i + 6) * fsi] -
	4316.36463814691 * f[(i + 0) * fsi] * f[(i + 7) * fsi] +
	496.071153428775 * f[(i + 0) * fsi] * f[(i + 8) * fsi] +
	29991.6199268498 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	180813.861301392 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	316864.946394454 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	351925.412140196 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	253037.274363286 * f[(i + 1) * fsi] * f[(i + 5) * fsi] -
	114802.39231254 * f[(i + 1) * fsi] * f[(i + 6) * fsi] +
	30004.6106074921 * f[(i + 1) * fsi] * f[(i + 7) * fsi] -
	3454.6250083913 * f[(i + 1) * fsi] * f[(i + 8) * fsi] +
	274437.463974042 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	966727.207640071 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	1077735.75959481 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	777141.668605927 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	353390.255601705 * f[(i + 2) * fsi] * f[(i + 6) * fsi] -
	92530.8427411113 * f[(i + 2) * fsi] * f[(i + 7) * fsi] +
	10669.6623191712 * f[(i + 2) * fsi] * f[(i + 8) * fsi] +
	854591.512480352 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	1911046.21916467 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	1381212.31577771 * f[(i + 3) * fsi] * f[(i + 5) * fsi] -
	629240.182223211 * f[(i + 3) * fsi] * f[(i + 6) * fsi] +
	165006.507929424 * f[(i + 3) * fsi] * f[(i + 7) * fsi] -
	19050.4296324557 * f[(i + 3) * fsi] * f[(i + 8) * fsi] +
	1070854.47449811 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	1550800.19015232 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	707565.248792344 * f[(i + 4) * fsi] * f[(i + 6) * fsi] -
	185776.153674587 * f[(i + 4) * fsi] * f[(i + 7) * fsi] +
	21470.5360236718 * f[(i + 4) * fsi] * f[(i + 8) * fsi] +
	562311.328233438 * f[(i + 5) * fsi] * f[(i + 5) * fsi] -
	513755.446640071 * f[(i + 5) * fsi] * f[(i + 6) * fsi] +
	135029.257900627 * f[(i + 5) * fsi] * f[(i + 7) * fsi] -
	15619.1414592001 * f[(i + 5) * fsi] * f[(i + 8) * fsi] +
	117469.122961696 * f[(i + 6) * fsi] * f[(i + 6) * fsi] -
	61801.7019274941 * f[(i + 6) * fsi] * f[(i + 7) * fsi] +
	7153.9713035639 * f[(i + 6) * fsi] * f[(i + 8) * fsi] +
	8134.55939472587 * f[(i + 7) * fsi] * f[(i + 7) * fsi] -
	1884.43224565446 * f[(i + 7) * fsi] * f[(i + 8) * fsi] +
	109.193772932945 * f[(i + 8) * fsi] * f[(i + 8) * fsi];
      sigma1 =
	+109.193772932945 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	1469.41675936423 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	4407.32664278072 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	7674.89153356352 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	8466.40114664635 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	6046.29475543025 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	2725.41239353462 * f[(i - 1) * fsi] * f[(i + 5) * fsi] -
	707.980347608125 * f[(i - 1) * fsi] * f[(i + 6) * fsi] +
	81.0556671385418 * f[(i - 1) * fsi] * f[(i + 7) * fsi] +
	5049.77020781835 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	30701.1587180034 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	53947.0217387824 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	59895.4013896723 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	42979.1691915766 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	19443.7923047687 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	5065.26308734773 * f[(i + 0) * fsi] * f[(i + 6) * fsi] -
	581.225261534807 * f[(i + 0) * fsi] * f[(i + 7) * fsi] +
	47140.2493458593 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	166921.455804809 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	186372.636616418 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	134309.272278381 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	60963.6198256632 * f[(i + 1) * fsi] * f[(i + 5) * fsi] -
	15923.1845243454 * f[(i + 1) * fsi] * f[(i + 6) * fsi] +
	1830.9895489579 * f[(i + 1) * fsi] * f[(i + 7) * fsi] +
	148657.090978519 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	333527.451742794 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	241247.066835384 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	109824.674852455 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	28752.7613106242 * f[(i + 2) * fsi] * f[(i + 6) * fsi] -
	3312.55790820803 * f[(i + 2) * fsi] * f[(i + 7) * fsi] +
	187797.717874362 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	272525.224659461 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	124378.188236788 * f[(i + 3) * fsi] * f[(i + 5) * fsi] -
	32630.2392534289 * f[(i + 3) * fsi] * f[(i + 6) * fsi] +
	3765.65529677863 * f[(i + 3) * fsi] * f[(i + 7) * fsi] +
	99127.2745988898 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	90677.1257492533 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	23831.4757387824 * f[(i + 4) * fsi] * f[(i + 6) * fsi] -
	2754.34352099738 * f[(i + 4) * fsi] * f[(i + 7) * fsi] +
	20774.7074754889 * f[(i + 5) * fsi] * f[(i + 5) * fsi] -
	10936.7033079505 * f[(i + 5) * fsi] * f[(i + 6) * fsi] +
	1265.66080746326 * f[(i + 5) * fsi] * f[(i + 7) * fsi] +
	1441.28575449258 * f[(i + 6) * fsi] * f[(i + 6) * fsi] -
	333.964212406558 * f[(i + 6) * fsi] * f[(i + 7) * fsi] +
	19.3647914042221 * f[(i + 7) * fsi] * f[(i + 7) * fsi];
      sigma2 =
	+19.3647914042221 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	267.510578137456 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	813.039719569186 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	1422.29540695141 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	1567.36952565594 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	1114.27213708535 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	498.941434911939 * f[(i - 2) * fsi] * f[(i + 4) * fsi] -
	128.604173640737 * f[(i - 2) * fsi] * f[(i + 5) * fsi] +
	14.6020328694406 * f[(i - 2) * fsi] * f[(i + 6) * fsi] +
	948.240872428061 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	5868.7702184994 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	10399.3092657059 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	11568.2032050107 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	8281.83632095821 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	3728.0916300002 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	964.845939313753 * f[(i - 1) * fsi] * f[(i + 5) * fsi] -
	109.897639186214 * f[(i - 1) * fsi] * f[(i + 6) * fsi] +
	9222.43045243719 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	33080.8713993305 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	37137.9417090107 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	26774.8153713592 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	12118.1388794826 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	3149.43314058221 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	359.908916834338 * f[(i + 0) * fsi] * f[(i + 6) * fsi] +
	29975.0953815866 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	67875.8127912121 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	49266.8129628913 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	22417.2306985197 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	5850.74664900054 * f[(i + 1) * fsi] * f[(i + 5) * fsi] -
	670.849344757305 * f[(i + 1) * fsi] * f[(i + 6) * fsi] +
	38710.2228777378 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	56543.4445813355 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	25858.8233448132 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	6776.71603569857 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	779.596278301446 * f[(i + 2) * fsi] * f[(i + 6) * fsi] +
	20760.5788136853 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	19076.8005289601 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	5018.72343149075 * f[(i + 3) * fsi] * f[(i + 5) * fsi] -
	579.197723970813 * f[(i + 3) * fsi] * f[(i + 6) * fsi] +
	4400.38698330138 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	2323.5095791696 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	268.954810839028 * f[(i + 4) * fsi] * f[(i + 6) * fsi] +
	307.688066683541 * f[(i + 5) * fsi] * f[(i + 5) * fsi] -
	71.4292240810191 * f[(i + 5) * fsi] * f[(i + 6) * fsi] +
	4.15594657554992 * f[(i + 6) * fsi] * f[(i + 6) * fsi];
      sigma3 =
	+4.15594657554992 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	60.2050054904579 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	189.330514253379 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	338.290107858048 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	376.449192281274 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	267.702258737133 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	119.001300721573 * f[(i - 3) * fsi] * f[(i + 3) * fsi] -
	30.2733426005663 * f[(i - 3) * fsi] * f[(i + 4) * fsi] +
	3.3778142788794 * f[(i - 3) * fsi] * f[(i + 5) * fsi] +
	224.5781681988 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	1445.81202311801 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	2631.07992925861 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	2970.48199593336 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	2136.83371274062 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	958.713082546492 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	245.728335017173 * f[(i - 2) * fsi] * f[(i + 4) * fsi] -
	27.586206325686 * f[(i - 2) * fsi] * f[(i + 5) * fsi] +
	2378.03262363703 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	8815.81240974712 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	10104.3776503687 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	7358.38198208451 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	3334.48145529283 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	861.616952916863 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	97.3685006774814 * f[(i - 1) * fsi] * f[(i + 5) * fsi] +
	8314.44047543303 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	19354.0800298309 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	14276.4079218034 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	6538.57561788698 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	1704.37865189406 * f[(i + 0) * fsi] * f[(i + 4) * fsi] -
	193.989288499038 * f[(i + 0) * fsi] * f[(i + 5) * fsi] +
	11427.8857755966 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	17079.2799526704 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	7909.63189419591 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	2081.09545492807 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	238.707145323602 * f[(i + 1) * fsi] * f[(i + 5) * fsi] +
	6460.89964518612 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	6051.52109493231 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	1607.7984936325 * f[(i + 2) * fsi] * f[(i + 4) * fsi] -
	185.954130124362 * f[(i + 2) * fsi] * f[(i + 5) * fsi] +
	1432.32903721728 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	768.643244458397 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	89.6803151793054 * f[(i + 3) * fsi] * f[(i + 5) * fsi] +
	104.120555009079 * f[(i + 4) * fsi] * f[(i + 4) * fsi] -
	24.5175956580063 * f[(i + 4) * fsi] * f[(i + 5) * fsi] +
	1.45672257391222 * f[(i + 5) * fsi] * f[(i + 5) * fsi];
      sigma4 =
	+1.45672257391222 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	22.8431920515406 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	77.2978189959941 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	147.360891739772 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	173.104800126842 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	128.386943302279 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	58.7752622928914 * f[(i - 4) * fsi] * f[(i + 2) * fsi] -
	15.2037101423747 * f[(i - 4) * fsi] * f[(i + 3) * fsi] +
	1.70341067241367 * f[(i - 4) * fsi] * f[(i + 4) * fsi] +
	91.7501465525254 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	634.284062414746 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	1231.84214048546 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	1470.62870986083 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	1106.32708286298 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	512.943219947286 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	134.233377963287 * f[(i - 3) * fsi] * f[(i + 3) * fsi] -
	15.2037101423747 * f[(i - 3) * fsi] * f[(i + 4) * fsi] +
	1119.38719626435 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	4433.56279439219 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	5386.9907367885 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	4116.54995777895 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	1935.50182392709 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	512.943219947286 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	58.7752622928914 * f[(i - 2) * fsi] * f[(i + 4) * fsi] +
	4477.71304825325 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	11088.1845350392 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	8620.4498023975 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	4116.54995777895 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	1106.32708286298 * f[(i - 1) * fsi] * f[(i + 3) * fsi] -
	128.386943302279 * f[(i - 1) * fsi] * f[(i + 4) * fsi] +
	6998.71770798472 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	11088.1845350392 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	5386.9907367885 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	1470.62870986083 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	173.104800126842 * f[(i + 0) * fsi] * f[(i + 4) * fsi] +
	4477.71304825325 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	4433.56279439219 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	1231.84214048546 * f[(i + 1) * fsi] * f[(i + 3) * fsi] -
	147.360891739772 * f[(i + 1) * fsi] * f[(i + 4) * fsi] +
	1119.38719626435 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	634.284062414746 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	77.2978189959941 * f[(i + 2) * fsi] * f[(i + 4) * fsi] +
	91.7501465525254 * f[(i + 3) * fsi] * f[(i + 3) * fsi] -
	22.8431920515406 * f[(i + 3) * fsi] * f[(i + 4) * fsi] +
	1.45672257391222 * f[(i + 4) * fsi] * f[(i + 4) * fsi];
      sigma5 =
	+1.45672257391222 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	24.5175956580063 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	89.6803151793054 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	185.954130124362 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	238.707145323602 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	193.989288499038 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	97.3685006774814 * f[(i - 5) * fsi] * f[(i + 1) * fsi] -
	27.586206325686 * f[(i - 5) * fsi] * f[(i + 2) * fsi] +
	3.3778142788794 * f[(i - 5) * fsi] * f[(i + 3) * fsi] +
	104.120555009079 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	768.643244458397 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	1607.7984936325 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	2081.09545492807 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	1704.37865189406 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	861.616952916863 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	245.728335017173 * f[(i - 4) * fsi] * f[(i + 2) * fsi] -
	30.2733426005663 * f[(i - 4) * fsi] * f[(i + 3) * fsi] +
	1432.32903721728 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	6051.52109493231 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	7909.63189419591 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	6538.57561788698 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	3334.48145529283 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	958.713082546492 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	119.001300721573 * f[(i - 3) * fsi] * f[(i + 3) * fsi] +
	6460.89964518612 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	17079.2799526704 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	14276.4079218034 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	7358.38198208451 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	2136.83371274062 * f[(i - 2) * fsi] * f[(i + 2) * fsi] -
	267.702258737133 * f[(i - 2) * fsi] * f[(i + 3) * fsi] +
	11427.8857755966 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	19354.0800298309 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	10104.3776503687 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	2970.48199593336 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	376.449192281274 * f[(i - 1) * fsi] * f[(i + 3) * fsi] +
	8314.44047543303 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	8815.81240974712 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	2631.07992925861 * f[(i + 0) * fsi] * f[(i + 2) * fsi] -
	338.290107858048 * f[(i + 0) * fsi] * f[(i + 3) * fsi] +
	2378.03262363703 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	1445.81202311801 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	189.330514253379 * f[(i + 1) * fsi] * f[(i + 3) * fsi] +
	224.5781681988 * f[(i + 2) * fsi] * f[(i + 2) * fsi] -
	60.2050054904579 * f[(i + 2) * fsi] * f[(i + 3) * fsi] +
	4.15594657554992 * f[(i + 3) * fsi] * f[(i + 3) * fsi];
      sigma6 =
	+4.15594657554992 * f[(i - 6) * fsi] * f[(i - 6) * fsi] -
	71.4292240810191 * f[(i - 6) * fsi] * f[(i - 5) * fsi] +
	268.954810839028 * f[(i - 6) * fsi] * f[(i - 4) * fsi] -
	579.197723970813 * f[(i - 6) * fsi] * f[(i - 3) * fsi] +
	779.596278301446 * f[(i - 6) * fsi] * f[(i - 2) * fsi] -
	670.849344757305 * f[(i - 6) * fsi] * f[(i - 1) * fsi] +
	359.908916834338 * f[(i - 6) * fsi] * f[(i + 0) * fsi] -
	109.897639186214 * f[(i - 6) * fsi] * f[(i + 1) * fsi] +
	14.6020328694406 * f[(i - 6) * fsi] * f[(i + 2) * fsi] +
	307.688066683541 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	2323.5095791696 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	5018.72343149075 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	6776.71603569857 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	5850.74664900054 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	3149.43314058221 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	964.845939313753 * f[(i - 5) * fsi] * f[(i + 1) * fsi] -
	128.604173640737 * f[(i - 5) * fsi] * f[(i + 2) * fsi] +
	4400.38698330138 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	19076.8005289601 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	25858.8233448132 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	22417.2306985197 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	12118.1388794826 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	3728.0916300002 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	498.941434911939 * f[(i - 4) * fsi] * f[(i + 2) * fsi] +
	20760.5788136853 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	56543.4445813355 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	49266.8129628913 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	26774.8153713592 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	8281.83632095821 * f[(i - 3) * fsi] * f[(i + 1) * fsi] -
	1114.27213708535 * f[(i - 3) * fsi] * f[(i + 2) * fsi] +
	38710.2228777378 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	67875.8127912121 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	37137.9417090107 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	11568.2032050107 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	1567.36952565594 * f[(i - 2) * fsi] * f[(i + 2) * fsi] +
	29975.0953815866 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	33080.8713993305 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	10399.3092657059 * f[(i - 1) * fsi] * f[(i + 1) * fsi] -
	1422.29540695141 * f[(i - 1) * fsi] * f[(i + 2) * fsi] +
	9222.43045243719 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	5868.7702184994 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	813.039719569186 * f[(i + 0) * fsi] * f[(i + 2) * fsi] +
	948.240872428061 * f[(i + 1) * fsi] * f[(i + 1) * fsi] -
	267.510578137456 * f[(i + 1) * fsi] * f[(i + 2) * fsi] +
	19.3647914042221 * f[(i + 2) * fsi] * f[(i + 2) * fsi];
      sigma7 =
	+19.3647914042221 * f[(i - 7) * fsi] * f[(i - 7) * fsi] -
	333.964212406558 * f[(i - 7) * fsi] * f[(i - 6) * fsi] +
	1265.66080746326 * f[(i - 7) * fsi] * f[(i - 5) * fsi] -
	2754.34352099738 * f[(i - 7) * fsi] * f[(i - 4) * fsi] +
	3765.65529677863 * f[(i - 7) * fsi] * f[(i - 3) * fsi] -
	3312.55790820803 * f[(i - 7) * fsi] * f[(i - 2) * fsi] +
	1830.9895489579 * f[(i - 7) * fsi] * f[(i - 1) * fsi] -
	581.225261534807 * f[(i - 7) * fsi] * f[(i + 0) * fsi] +
	81.0556671385418 * f[(i - 7) * fsi] * f[(i + 1) * fsi] +
	1441.28575449258 * f[(i - 6) * fsi] * f[(i - 6) * fsi] -
	10936.7033079505 * f[(i - 6) * fsi] * f[(i - 5) * fsi] +
	23831.4757387824 * f[(i - 6) * fsi] * f[(i - 4) * fsi] -
	32630.2392534289 * f[(i - 6) * fsi] * f[(i - 3) * fsi] +
	28752.7613106242 * f[(i - 6) * fsi] * f[(i - 2) * fsi] -
	15923.1845243454 * f[(i - 6) * fsi] * f[(i - 1) * fsi] +
	5065.26308734773 * f[(i - 6) * fsi] * f[(i + 0) * fsi] -
	707.980347608125 * f[(i - 6) * fsi] * f[(i + 1) * fsi] +
	20774.7074754889 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	90677.1257492533 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	124378.188236788 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	109824.674852455 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	60963.6198256632 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	19443.7923047687 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	2725.41239353462 * f[(i - 5) * fsi] * f[(i + 1) * fsi] +
	99127.2745988898 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	272525.224659461 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	241247.066835384 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	134309.272278381 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	42979.1691915766 * f[(i - 4) * fsi] * f[(i + 0) * fsi] -
	6046.29475543025 * f[(i - 4) * fsi] * f[(i + 1) * fsi] +
	187797.717874362 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	333527.451742794 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	186372.636616418 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	59895.4013896723 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	8466.40114664635 * f[(i - 3) * fsi] * f[(i + 1) * fsi] +
	148657.090978519 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	166921.455804809 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	53947.0217387824 * f[(i - 2) * fsi] * f[(i + 0) * fsi] -
	7674.89153356352 * f[(i - 2) * fsi] * f[(i + 1) * fsi] +
	47140.2493458593 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	30701.1587180034 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	4407.32664278072 * f[(i - 1) * fsi] * f[(i + 1) * fsi] +
	5049.77020781835 * f[(i + 0) * fsi] * f[(i + 0) * fsi] -
	1469.41675936423 * f[(i + 0) * fsi] * f[(i + 1) * fsi] +
	109.193772932945 * f[(i + 1) * fsi] * f[(i + 1) * fsi];
      sigma8 =
	+109.193772932945 * f[(i - 8) * fsi] * f[(i - 8) * fsi] -
	1884.43224565446 * f[(i - 8) * fsi] * f[(i - 7) * fsi] +
	7153.9713035639 * f[(i - 8) * fsi] * f[(i - 6) * fsi] -
	15619.1414592001 * f[(i - 8) * fsi] * f[(i - 5) * fsi] +
	21470.5360236718 * f[(i - 8) * fsi] * f[(i - 4) * fsi] -
	19050.4296324557 * f[(i - 8) * fsi] * f[(i - 3) * fsi] +
	10669.6623191712 * f[(i - 8) * fsi] * f[(i - 2) * fsi] -
	3454.6250083913 * f[(i - 8) * fsi] * f[(i - 1) * fsi] +
	496.071153428775 * f[(i - 8) * fsi] * f[(i + 0) * fsi] +
	8134.55939472587 * f[(i - 7) * fsi] * f[(i - 7) * fsi] -
	61801.7019274941 * f[(i - 7) * fsi] * f[(i - 6) * fsi] +
	135029.257900627 * f[(i - 7) * fsi] * f[(i - 5) * fsi] -
	185776.153674587 * f[(i - 7) * fsi] * f[(i - 4) * fsi] +
	165006.507929424 * f[(i - 7) * fsi] * f[(i - 3) * fsi] -
	92530.8427411113 * f[(i - 7) * fsi] * f[(i - 2) * fsi] +
	30004.6106074921 * f[(i - 7) * fsi] * f[(i - 1) * fsi] -
	4316.36463814691 * f[(i - 7) * fsi] * f[(i + 0) * fsi] +
	117469.122961696 * f[(i - 6) * fsi] * f[(i - 6) * fsi] -
	513755.446640071 * f[(i - 6) * fsi] * f[(i - 5) * fsi] +
	707565.248792344 * f[(i - 6) * fsi] * f[(i - 4) * fsi] -
	629240.182223211 * f[(i - 6) * fsi] * f[(i - 3) * fsi] +
	353390.255601705 * f[(i - 6) * fsi] * f[(i - 2) * fsi] -
	114802.39231254 * f[(i - 6) * fsi] * f[(i - 1) * fsi] +
	16552.0014823105 * f[(i - 6) * fsi] * f[(i + 0) * fsi] +
	562311.328233438 * f[(i - 5) * fsi] * f[(i - 5) * fsi] -
	1550800.19015232 * f[(i - 5) * fsi] * f[(i - 4) * fsi] +
	1381212.31577771 * f[(i - 5) * fsi] * f[(i - 3) * fsi] -
	777141.668605927 * f[(i - 5) * fsi] * f[(i - 2) * fsi] +
	253037.274363286 * f[(i - 5) * fsi] * f[(i - 1) * fsi] -
	36585.0576509742 * f[(i - 5) * fsi] * f[(i + 0) * fsi] +
	1070854.47449811 * f[(i - 4) * fsi] * f[(i - 4) * fsi] -
	1911046.21916467 * f[(i - 4) * fsi] * f[(i - 3) * fsi] +
	1077735.75959481 * f[(i - 4) * fsi] * f[(i - 2) * fsi] -
	351925.412140196 * f[(i - 4) * fsi] * f[(i - 1) * fsi] +
	51067.48172473 * f[(i - 4) * fsi] * f[(i + 0) * fsi] +
	854591.512480352 * f[(i - 3) * fsi] * f[(i - 3) * fsi] -
	966727.207640071 * f[(i - 3) * fsi] * f[(i - 2) * fsi] +
	316864.946394454 * f[(i - 3) * fsi] * f[(i - 1) * fsi] -
	46202.7564018809 * f[(i - 3) * fsi] * f[(i + 0) * fsi] +
	274437.463974042 * f[(i - 2) * fsi] * f[(i - 2) * fsi] -
	180813.861301392 * f[(i - 2) * fsi] * f[(i - 1) * fsi] +
	26542.9748247279 * f[(i - 2) * fsi] * f[(i + 0) * fsi] +
	29991.6199268498 * f[(i - 1) * fsi] * f[(i - 1) * fsi] -
	8893.78045641284 * f[(i - 1) * fsi] * f[(i + 0) * fsi] +
	669.714981108808 * f[(i + 0) * fsi] * f[(i + 0) * fsi];
      sigma[i * ssi + 0 * ssr] = sigma0;
      sigma[i * ssi + 1 * ssr] = sigma1;
      sigma[i * ssi + 2 * ssr] = sigma2;
      sigma[i * ssi + 3 * ssr] = sigma3;
      sigma[i * ssi + 4 * ssr] = sigma4;
      sigma[i * ssi + 5 * ssr] = sigma5;
      sigma[i * ssi + 6 * ssr] = sigma6;
      sigma[i * ssi + 7 * ssr] = sigma7;
      sigma[i * ssi + 8 * ssr] = sigma8;
    }
}


PyObject *
py_smoothness_k9 (PyObject * self, PyObject * args)
{
  double *sigma, *f;
  PyArrayObject *f_py, *sigma_py;

  long int n;
  int ssi, ssr, fsi;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &f_py, &sigma_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  f = (double *) PyArray_DATA (f_py);

  n = PyArray_DIM (f_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);

  smoothness_k9 (f, n, fsi, sigma, ssi, ssr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_left_k9 (const double *restrict sigma, int n, int ssi, int ssr,
		 double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7,
    sigma8;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7,
    omega8;
  for (i = 8; i < n - 8; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      sigma6 = sigma[i * ssi + 6 * ssr];
      sigma7 = sigma[i * ssi + 7 * ssr];
      sigma8 = sigma[i * ssi + 8 * ssr];
      accumulator = 0.0;
      omega0 = +4.11353352529823e-05 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.00296174413821473 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.0414644179350062 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.193500617030029 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.362813656931304 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega5 = +0.290250925545043 / (1e-36 + sigma5) / (1e-36 + sigma5);
      accumulator += omega5;
      omega6 = +0.0967503085150144 / (1e-36 + sigma6) / (1e-36 + sigma6);
      accumulator += omega6;
      omega7 = +0.0118469765528589 / (1e-36 + sigma7) / (1e-36 + sigma7);
      accumulator += omega7;
      omega8 = +0.000370218017276841 / (1e-36 + sigma8) / (1e-36 + sigma8);
      accumulator += omega8;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega5 /= accumulator;
      omega6 /= accumulator;
      omega7 /= accumulator;
      omega8 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr] = omega5;
      omega[i * wsi + 0 * wsl + 6 * wsr] = omega6;
      omega[i * wsi + 0 * wsl + 7 * wsr] = omega7;
      omega[i * wsi + 0 * wsl + 8 * wsr] = omega8;
    }
}



PyObject *
py_weights_left_k9 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_left_k9 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_left_k9 (const double *restrict f, int n, int fsi,
		     const double *restrict omega, int wsi, int wsl, int wsr,
		     double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7,
    omega8;
  for (i = 8; i < n - 8; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      omega5 = omega[i * wsi + 0 * wsl + 5 * wsr];
      omega6 = omega[i * wsi + 0 * wsl + 6 * wsr];
      omega7 = omega[i * wsi + 0 * wsl + 7 * wsr];
      omega8 = omega[i * wsi + 0 * wsl + 8 * wsr];
      fr0 =
	+2.82896825396825 * f[(i + 0) * fsi] +
	-6.17103174603175 * f[(i + 1) * fsi] +
	+11.8289682539683 * f[(i + 2) * fsi] +
	-16.1710317460317 * f[(i + 3) * fsi] +
	+15.3289682539683 * f[(i + 4) * fsi] +
	-9.87103174603175 * f[(i + 5) * fsi] +
	+4.12896825396825 * f[(i + 6) * fsi] +
	-1.01388888888889 * f[(i + 7) * fsi] +
	+0.111111111111111 * f[(i + 8) * fsi];
      fr1 =
	+0.111111111111111 * f[(i - 1) * fsi] +
	+1.82896825396825 * f[(i + 0) * fsi] +
	-2.17103174603175 * f[(i + 1) * fsi] +
	+2.49563492063492 * f[(i + 2) * fsi] +
	-2.17103174603175 * f[(i + 3) * fsi] +
	+1.32896825396825 * f[(i + 4) * fsi] +
	-0.537698412698413 * f[(i + 5) * fsi] +
	+0.128968253968254 * f[(i + 6) * fsi] +
	-0.0138888888888889 * f[(i + 7) * fsi];
      fr2 =
	-0.0138888888888889 * f[(i - 2) * fsi] +
	+0.236111111111111 * f[(i - 1) * fsi] +
	+1.32896825396825 * f[(i + 0) * fsi] +
	-1.00436507936508 * f[(i + 1) * fsi] +
	+0.745634920634921 * f[(i + 2) * fsi] +
	-0.421031746031746 * f[(i + 3) * fsi] +
	+0.162301587301587 * f[(i + 4) * fsi] +
	-0.0376984126984127 * f[(i + 5) * fsi] +
	+0.00396825396825397 * f[(i + 6) * fsi];
      fr3 =
	+0.00396825396825397 * f[(i - 3) * fsi] +
	-0.0496031746031746 * f[(i - 2) * fsi] +
	+0.378968253968254 * f[(i - 1) * fsi] +
	+0.995634920634921 * f[(i + 0) * fsi] +
	-0.504365079365079 * f[(i + 1) * fsi] +
	+0.245634920634921 * f[(i + 2) * fsi] +
	-0.0876984126984127 * f[(i + 3) * fsi] +
	+0.0194444444444444 * f[(i + 4) * fsi] +
	-0.00198412698412698 * f[(i + 5) * fsi];
      fr4 =
	-0.00198412698412698 * f[(i - 4) * fsi] +
	+0.0218253968253968 * f[(i - 3) * fsi] +
	-0.121031746031746 * f[(i - 2) * fsi] +
	+0.545634920634921 * f[(i - 1) * fsi] +
	+0.745634920634921 * f[(i + 0) * fsi] +
	-0.254365079365079 * f[(i + 1) * fsi] +
	+0.078968253968254 * f[(i + 2) * fsi] +
	-0.0162698412698413 * f[(i + 3) * fsi] +
	+0.00158730158730159 * f[(i + 4) * fsi];
      fr5 =
	+0.00158730158730159 * f[(i - 5) * fsi] +
	-0.0162698412698413 * f[(i - 4) * fsi] +
	+0.078968253968254 * f[(i - 3) * fsi] +
	-0.254365079365079 * f[(i - 2) * fsi] +
	+0.745634920634921 * f[(i - 1) * fsi] +
	+0.545634920634921 * f[(i + 0) * fsi] +
	-0.121031746031746 * f[(i + 1) * fsi] +
	+0.0218253968253968 * f[(i + 2) * fsi] +
	-0.00198412698412698 * f[(i + 3) * fsi];
      fr6 =
	-0.00198412698412698 * f[(i - 6) * fsi] +
	+0.0194444444444444 * f[(i - 5) * fsi] +
	-0.0876984126984127 * f[(i - 4) * fsi] +
	+0.245634920634921 * f[(i - 3) * fsi] +
	-0.504365079365079 * f[(i - 2) * fsi] +
	+0.995634920634921 * f[(i - 1) * fsi] +
	+0.378968253968254 * f[(i + 0) * fsi] +
	-0.0496031746031746 * f[(i + 1) * fsi] +
	+0.00396825396825397 * f[(i + 2) * fsi];
      fr7 =
	+0.00396825396825397 * f[(i - 7) * fsi] +
	-0.0376984126984127 * f[(i - 6) * fsi] +
	+0.162301587301587 * f[(i - 5) * fsi] +
	-0.421031746031746 * f[(i - 4) * fsi] +
	+0.745634920634921 * f[(i - 3) * fsi] +
	-1.00436507936508 * f[(i - 2) * fsi] +
	+1.32896825396825 * f[(i - 1) * fsi] +
	+0.236111111111111 * f[(i + 0) * fsi] +
	-0.0138888888888889 * f[(i + 1) * fsi];
      fr8 =
	-0.0138888888888889 * f[(i - 8) * fsi] +
	+0.128968253968254 * f[(i - 7) * fsi] +
	-0.537698412698413 * f[(i - 6) * fsi] +
	+1.32896825396825 * f[(i - 5) * fsi] +
	-2.17103174603175 * f[(i - 4) * fsi] +
	+2.49563492063492 * f[(i - 3) * fsi] +
	-2.17103174603175 * f[(i - 2) * fsi] +
	+1.82896825396825 * f[(i - 1) * fsi] +
	+0.111111111111111 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4 + fr5 * omega5 + fr6 * omega6 + fr7 * omega7 +
	fr8 * omega8;
    }
}

PyObject *
py_reconstruct_left_k9 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_left_k9 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

void
weights_right_k9 (const double *restrict sigma, int n, int ssi, int ssr,
		  double *restrict omega, int wsi, int wsl, int wsr)
{
  int i;
  double accumulator;
  double sigma0, sigma1, sigma2, sigma3, sigma4, sigma5, sigma6, sigma7,
    sigma8;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7,
    omega8;
  for (i = 8; i < n - 8; i++)
    {
      sigma0 = sigma[i * ssi + 0 * ssr];
      sigma1 = sigma[i * ssi + 1 * ssr];
      sigma2 = sigma[i * ssi + 2 * ssr];
      sigma3 = sigma[i * ssi + 3 * ssr];
      sigma4 = sigma[i * ssi + 4 * ssr];
      sigma5 = sigma[i * ssi + 5 * ssr];
      sigma6 = sigma[i * ssi + 6 * ssr];
      sigma7 = sigma[i * ssi + 7 * ssr];
      sigma8 = sigma[i * ssi + 8 * ssr];
      accumulator = 0.0;
      omega0 = +0.000370218017276841 / (1e-36 + sigma0) / (1e-36 + sigma0);
      accumulator += omega0;
      omega1 = +0.0118469765528589 / (1e-36 + sigma1) / (1e-36 + sigma1);
      accumulator += omega1;
      omega2 = +0.0967503085150144 / (1e-36 + sigma2) / (1e-36 + sigma2);
      accumulator += omega2;
      omega3 = +0.290250925545043 / (1e-36 + sigma3) / (1e-36 + sigma3);
      accumulator += omega3;
      omega4 = +0.362813656931304 / (1e-36 + sigma4) / (1e-36 + sigma4);
      accumulator += omega4;
      omega5 = +0.193500617030029 / (1e-36 + sigma5) / (1e-36 + sigma5);
      accumulator += omega5;
      omega6 = +0.0414644179350062 / (1e-36 + sigma6) / (1e-36 + sigma6);
      accumulator += omega6;
      omega7 = +0.00296174413821473 / (1e-36 + sigma7) / (1e-36 + sigma7);
      accumulator += omega7;
      omega8 = +4.11353352529823e-05 / (1e-36 + sigma8) / (1e-36 + sigma8);
      accumulator += omega8;
      omega0 /= accumulator;
      omega1 /= accumulator;
      omega2 /= accumulator;
      omega3 /= accumulator;
      omega4 /= accumulator;
      omega5 /= accumulator;
      omega6 /= accumulator;
      omega7 /= accumulator;
      omega8 /= accumulator;
      omega[i * wsi + 0 * wsl + 0 * wsr] = omega0;
      omega[i * wsi + 0 * wsl + 1 * wsr] = omega1;
      omega[i * wsi + 0 * wsl + 2 * wsr] = omega2;
      omega[i * wsi + 0 * wsl + 3 * wsr] = omega3;
      omega[i * wsi + 0 * wsl + 4 * wsr] = omega4;
      omega[i * wsi + 0 * wsl + 5 * wsr] = omega5;
      omega[i * wsi + 0 * wsl + 6 * wsr] = omega6;
      omega[i * wsi + 0 * wsl + 7 * wsr] = omega7;
      omega[i * wsi + 0 * wsl + 8 * wsr] = omega8;
    }
}



PyObject *
py_weights_right_k9 (PyObject * self, PyObject * args)
{
  double *sigma, *omega;
  PyArrayObject *sigma_py, *omega_py;

  long int n;
  int ssi, ssr, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OO", &sigma_py, &omega_py))
    return NULL;

  if (PyArray_NDIM(sigma_py) != 2 || PyArray_TYPE(sigma_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "sigma must be two-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  sigma = (double *) PyArray_DATA (sigma_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  ssi = PyArray_STRIDES(sigma_py)[0] / sizeof (double);
  ssr = PyArray_STRIDES(sigma_py)[1] / sizeof (double);

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  weights_right_k9 (sigma, n, ssi, ssr, omega, wsi, wsl, wsr);

  Py_INCREF (Py_None);
  return Py_None;
}

void
reconstruct_right_k9 (const double *restrict f, int n, int fsi,
		      const double *restrict omega, int wsi, int wsl, int wsr,
		      double *restrict fr, int frsi, int frsl)
{
  int i;
  double fr0, fr1, fr2, fr3, fr4, fr5, fr6, fr7, fr8;
  double omega0, omega1, omega2, omega3, omega4, omega5, omega6, omega7,
    omega8;
  for (i = 8; i < n - 8; i++)
    {
      omega0 = omega[i * wsi + 0 * wsl + 0 * wsr];
      omega1 = omega[i * wsi + 0 * wsl + 1 * wsr];
      omega2 = omega[i * wsi + 0 * wsl + 2 * wsr];
      omega3 = omega[i * wsi + 0 * wsl + 3 * wsr];
      omega4 = omega[i * wsi + 0 * wsl + 4 * wsr];
      omega5 = omega[i * wsi + 0 * wsl + 5 * wsr];
      omega6 = omega[i * wsi + 0 * wsl + 6 * wsr];
      omega7 = omega[i * wsi + 0 * wsl + 7 * wsr];
      omega8 = omega[i * wsi + 0 * wsl + 8 * wsr];
      fr0 =
	+0.111111111111111 * f[(i + 0) * fsi] +
	+1.82896825396825 * f[(i + 1) * fsi] +
	-2.17103174603175 * f[(i + 2) * fsi] +
	+2.49563492063492 * f[(i + 3) * fsi] +
	-2.17103174603175 * f[(i + 4) * fsi] +
	+1.32896825396825 * f[(i + 5) * fsi] +
	-0.537698412698413 * f[(i + 6) * fsi] +
	+0.128968253968254 * f[(i + 7) * fsi] +
	-0.0138888888888889 * f[(i + 8) * fsi];
      fr1 =
	-0.0138888888888889 * f[(i - 1) * fsi] +
	+0.236111111111111 * f[(i + 0) * fsi] +
	+1.32896825396825 * f[(i + 1) * fsi] +
	-1.00436507936508 * f[(i + 2) * fsi] +
	+0.745634920634921 * f[(i + 3) * fsi] +
	-0.421031746031746 * f[(i + 4) * fsi] +
	+0.162301587301587 * f[(i + 5) * fsi] +
	-0.0376984126984127 * f[(i + 6) * fsi] +
	+0.00396825396825397 * f[(i + 7) * fsi];
      fr2 =
	+0.00396825396825397 * f[(i - 2) * fsi] +
	-0.0496031746031746 * f[(i - 1) * fsi] +
	+0.378968253968254 * f[(i + 0) * fsi] +
	+0.995634920634921 * f[(i + 1) * fsi] +
	-0.504365079365079 * f[(i + 2) * fsi] +
	+0.245634920634921 * f[(i + 3) * fsi] +
	-0.0876984126984127 * f[(i + 4) * fsi] +
	+0.0194444444444444 * f[(i + 5) * fsi] +
	-0.00198412698412698 * f[(i + 6) * fsi];
      fr3 =
	-0.00198412698412698 * f[(i - 3) * fsi] +
	+0.0218253968253968 * f[(i - 2) * fsi] +
	-0.121031746031746 * f[(i - 1) * fsi] +
	+0.545634920634921 * f[(i + 0) * fsi] +
	+0.745634920634921 * f[(i + 1) * fsi] +
	-0.254365079365079 * f[(i + 2) * fsi] +
	+0.078968253968254 * f[(i + 3) * fsi] +
	-0.0162698412698413 * f[(i + 4) * fsi] +
	+0.00158730158730159 * f[(i + 5) * fsi];
      fr4 =
	+0.00158730158730159 * f[(i - 4) * fsi] +
	-0.0162698412698413 * f[(i - 3) * fsi] +
	+0.078968253968254 * f[(i - 2) * fsi] +
	-0.254365079365079 * f[(i - 1) * fsi] +
	+0.745634920634921 * f[(i + 0) * fsi] +
	+0.545634920634921 * f[(i + 1) * fsi] +
	-0.121031746031746 * f[(i + 2) * fsi] +
	+0.0218253968253968 * f[(i + 3) * fsi] +
	-0.00198412698412698 * f[(i + 4) * fsi];
      fr5 =
	-0.00198412698412698 * f[(i - 5) * fsi] +
	+0.0194444444444444 * f[(i - 4) * fsi] +
	-0.0876984126984127 * f[(i - 3) * fsi] +
	+0.245634920634921 * f[(i - 2) * fsi] +
	-0.504365079365079 * f[(i - 1) * fsi] +
	+0.995634920634921 * f[(i + 0) * fsi] +
	+0.378968253968254 * f[(i + 1) * fsi] +
	-0.0496031746031746 * f[(i + 2) * fsi] +
	+0.00396825396825397 * f[(i + 3) * fsi];
      fr6 =
	+0.00396825396825397 * f[(i - 6) * fsi] +
	-0.0376984126984127 * f[(i - 5) * fsi] +
	+0.162301587301587 * f[(i - 4) * fsi] +
	-0.421031746031746 * f[(i - 3) * fsi] +
	+0.745634920634921 * f[(i - 2) * fsi] +
	-1.00436507936508 * f[(i - 1) * fsi] +
	+1.32896825396825 * f[(i + 0) * fsi] +
	+0.236111111111111 * f[(i + 1) * fsi] +
	-0.0138888888888889 * f[(i + 2) * fsi];
      fr7 =
	-0.0138888888888889 * f[(i - 7) * fsi] +
	+0.128968253968254 * f[(i - 6) * fsi] +
	-0.537698412698413 * f[(i - 5) * fsi] +
	+1.32896825396825 * f[(i - 4) * fsi] +
	-2.17103174603175 * f[(i - 3) * fsi] +
	+2.49563492063492 * f[(i - 2) * fsi] +
	-2.17103174603175 * f[(i - 1) * fsi] +
	+1.82896825396825 * f[(i + 0) * fsi] +
	+0.111111111111111 * f[(i + 1) * fsi];
      fr8 =
	+0.111111111111111 * f[(i - 8) * fsi] +
	-1.01388888888889 * f[(i - 7) * fsi] +
	+4.12896825396825 * f[(i - 6) * fsi] +
	-9.87103174603175 * f[(i - 5) * fsi] +
	+15.3289682539683 * f[(i - 4) * fsi] +
	-16.1710317460317 * f[(i - 3) * fsi] +
	+11.8289682539683 * f[(i - 2) * fsi] +
	-6.17103174603175 * f[(i - 1) * fsi] +
	+2.82896825396825 * f[(i + 0) * fsi];
      fr[i * frsi + 0 * frsl] =
	fr0 * omega0 + fr1 * omega1 + fr2 * omega2 + fr3 * omega3 +
	fr4 * omega4 + fr5 * omega5 + fr6 * omega6 + fr7 * omega7 +
	fr8 * omega8;
    }
}

PyObject *
py_reconstruct_right_k9 (PyObject * self, PyObject * args)
{
  double *f, *omega, *fr;
  PyArrayObject *f_py, *omega_py, *fr_py;

  long int n;
  int fsi, frsi, frsl, wsi, wsl, wsr;

  /* parse options */

  if (!PyArg_ParseTuple (args, "OOO", &f_py, &omega_py, &fr_py))
    return NULL;

  if (PyArray_NDIM(f_py) != 1 || PyArray_TYPE(f_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError,
		       "f must be one-dimensional and of type float");
      return NULL;
    }

  if (PyArray_TYPE(fr_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "fr must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(fr_py) == 1 || PyArray_NDIM(fr_py) == 2))
    {
      PyErr_SetString (PyExc_ValueError, "fr must be one or two dimensional");
      return NULL;
    }

  if (PyArray_TYPE(omega_py) != NPY_DOUBLE)
    {
      PyErr_SetString (PyExc_ValueError, "omega must be of type float");
      return NULL;
    }

  if (!(PyArray_NDIM(omega_py) >= 2 && PyArray_NDIM(omega_py) <= 4))
    {
      PyErr_SetString (PyExc_ValueError,
		       "omega must be two, three, or four dimensional");
      return NULL;
    }

  /* get data, n, strides */
  f = (double *) PyArray_DATA (f_py);
  fr = (double *) PyArray_DATA (fr_py);
  omega = (double *) PyArray_DATA (omega_py);

  n = PyArray_DIM (omega_py, 0);

  fsi = PyArray_STRIDES(f_py)[0] / sizeof (double);
  frsi = PyArray_STRIDES(fr_py)[0] / sizeof (double);

  if (n == 1)
    {
      frsl = 0;
    }
  else
    {
      frsl = PyArray_STRIDES(fr_py)[1] / sizeof (double);
    }

  wsi = PyArray_STRIDES(omega_py)[0] / sizeof (double);
  if (PyArray_NDIM(omega_py) == 3)
    {
      wsl = PyArray_STRIDES(omega_py)[1] / sizeof (double);
      wsr = PyArray_STRIDES(omega_py)[2] / sizeof (double);
    }
  else
    {
      wsl = 0;
      wsr = PyArray_STRIDES(omega_py)[1] / sizeof (double);
    }

  reconstruct_right_k9 (f, n, fsi, omega, wsi, wsl, wsr, fr, frsi, frsl);

  Py_INCREF (Py_None);
  return Py_None;
}

static PyMethodDef reconstructmethods[] = {
  {"smoothness_k3", py_smoothness_k3, METH_VARARGS, ""},
  {"weights_left_k3", py_weights_left_k3, METH_VARARGS, ""},
  {"reconstruct_left_k3", py_reconstruct_left_k3, METH_VARARGS, ""},
  {"weights_right_k3", py_weights_right_k3, METH_VARARGS, ""},
  {"reconstruct_right_k3", py_reconstruct_right_k3, METH_VARARGS, ""},
  {"smoothness_k4", py_smoothness_k4, METH_VARARGS, ""},
  {"weights_left_k4", py_weights_left_k4, METH_VARARGS, ""},
  {"reconstruct_left_k4", py_reconstruct_left_k4, METH_VARARGS, ""},
  {"weights_right_k4", py_weights_right_k4, METH_VARARGS, ""},
  {"reconstruct_right_k4", py_reconstruct_right_k4, METH_VARARGS, ""},
  {"smoothness_k5", py_smoothness_k5, METH_VARARGS, ""},
  {"weights_left_k5", py_weights_left_k5, METH_VARARGS, ""},
  {"reconstruct_left_k5", py_reconstruct_left_k5, METH_VARARGS, ""},
  {"weights_right_k5", py_weights_right_k5, METH_VARARGS, ""},
  {"reconstruct_right_k5", py_reconstruct_right_k5, METH_VARARGS, ""},
  {"smoothness_k6", py_smoothness_k6, METH_VARARGS, ""},
  {"weights_left_k6", py_weights_left_k6, METH_VARARGS, ""},
  {"reconstruct_left_k6", py_reconstruct_left_k6, METH_VARARGS, ""},
  {"weights_right_k6", py_weights_right_k6, METH_VARARGS, ""},
  {"reconstruct_right_k6", py_reconstruct_right_k6, METH_VARARGS, ""},
  {"smoothness_k7", py_smoothness_k7, METH_VARARGS, ""},
  {"weights_left_k7", py_weights_left_k7, METH_VARARGS, ""},
  {"reconstruct_left_k7", py_reconstruct_left_k7, METH_VARARGS, ""},
  {"weights_right_k7", py_weights_right_k7, METH_VARARGS, ""},
  {"reconstruct_right_k7", py_reconstruct_right_k7, METH_VARARGS, ""},
  {"smoothness_k8", py_smoothness_k8, METH_VARARGS, ""},
  {"weights_left_k8", py_weights_left_k8, METH_VARARGS, ""},
  {"reconstruct_left_k8", py_reconstruct_left_k8, METH_VARARGS, ""},
  {"weights_right_k8", py_weights_right_k8, METH_VARARGS, ""},
  {"reconstruct_right_k8", py_reconstruct_right_k8, METH_VARARGS, ""},
  {"smoothness_k9", py_smoothness_k9, METH_VARARGS, ""},
  {"weights_left_k9", py_weights_left_k9, METH_VARARGS, ""},
  {"reconstruct_left_k9", py_reconstruct_left_k9, METH_VARARGS, ""},
  {"weights_right_k9", py_weights_right_k9, METH_VARARGS, ""},
  {"reconstruct_right_k9", py_reconstruct_right_k9, METH_VARARGS, ""},
  {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC
PyInit_reconstruct (void)
{
  static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "reconstruct",       /* m_name */
    NULL,                /* m_doc */
    -1,                  /* m_size */
    reconstructmethods,  /* m_methods */
    NULL,                /* m_reload */
    NULL,                /* m_traverse */
    NULL,                /* m_clear */
    NULL                 /* m_free */
  };
  import_array();
  return PyModule_Create(&moduledef);
}
