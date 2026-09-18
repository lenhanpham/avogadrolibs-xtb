/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_eg.f90 (egbond, egbend, egtors, gfnffdampa/dampt),
  src/gfnff/math.f (vsub, vlen, valijklffPBC), src/basic_geo.f90
  (crossprod, vecnorm, impsc) and src/constr.f90 (dphidrPBC, omegaPBC,
  domegadrPBC), Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffegbond.h"

#include "constants.h"
#include "environment.h"
#include "gfnffdata.h"
#include "gfnffparams.h"
#include "gfnffrab.h"

#include <algorithm>
#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

void vecSub(const double a[3], const double b[3], double c[3])
{
  for (int i = 0; i < 3; ++i)
    c[i] = a[i] - b[i];
}

void vecCross(const double a[3], const double b[3], double c[3])
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

double vecLen(const double a[3])
{
  double tot = a[0] * a[0] + a[1] * a[1] + a[2] * a[2];
  double vlen = 0.0;
  if (tot > 0.0)
    vlen = std::sqrt(tot);
  return vlen;
}

// gfortran NORM2 replication (see gfnffoop): scale starts at 1.0, the
// ratio is squared before folding in. Bitwise on the validation battery.
double norm2(const double a[3])
{
  double scale = 1.0, ssq = 0.0;
  for (int k = 0; k < 3; ++k) {
    if (a[k] == 0.0)
      continue;
    double ax = std::fabs(a[k]);
    if (ax > scale) {
      double r = scale / ax;
      ssq = 1.0 + ssq * (r * r);
      scale = ax;
    } else {
      double q = ax / scale;
      ssq = ssq + q * q;
    }
  }
  return scale * std::sqrt(ssq);
}

double vecCos(const double a[3], const double b[3])
{
  double rimp = 0.0;
  for (int i = 0; i < 3; ++i)
    rimp = rimp + a[i] * b[i];
  double al = norm2(a);
  double bl = norm2(b);
  if (al > 0.0 && bl > 0.0)
    return rimp / (al * bl);
  return 0.0;
}

double vecNorm(double r[3], int inorm)
{
  double sp = 0.0;
  for (int i = 0; i < 3; ++i)
    sp = sp + r[i] * r[i];
  double rn = std::sqrt(sp);
  if (inorm > 0 && std::fabs(rn) > 1.0e-14) {
    double or_ = 1.0 / rn;
    for (int i = 0; i < 3; ++i)
      r[i] = or_ * r[i];
  }
  return rn;
}

void bondDamp(double r2, double atcut, const std::vector<double>& rcov,
              int ia, int ja, double& damp, double& ddamp)
{
  double sum = rcov[ia - 1] + rcov[ja - 1];
  double rcut = atcut * (sum * sum);
  double q = r2 / rcut;
  double rr = q * q;
  damp = 1.0 / (1.0 + rr);
  double d = 1.0 + rr;
  ddamp = (-4.0 * rr) / (r2 * (d * d));
}

double dihedralAngle(const std::vector<double>& xyz, int i, int j, int k,
                     int l)
{
  double ra[3], rb[3], rc[3], na[3], nb[3];
  for (int c = 0; c < 3; ++c) {
    ra[c] = xyz[3 * j + c] - xyz[3 * i + c];
    rb[c] = xyz[3 * k + c] - xyz[3 * j + c];
    rc[c] = xyz[3 * l + c] - xyz[3 * k + c];
  }
  // Determinant and bend calls are pure with unused outputs; skipped.
  vecCross(ra, rb, na);
  vecCross(rb, rc, nb);
  vecNorm(na, 1);
  vecNorm(nb, 1);
  double snanb = 0.0;
  for (int c = 0; c < 3; ++c)
    snanb = snanb + na[c] * nb[c];
  if (std::fabs(std::fabs(snanb) - 1.0) < 1.0e-14)
    snanb = snanb >= 0.0 ? 1.0 : -1.0;
  return std::acos(snanb);
}

void dihedralDerivatives(const std::vector<double>& xyz, int i, int j, int k,
                         int l, double phi, double da[3], double db[3],
                         double dc[3], double dd[3])
{
  double ra[3], rb[3], rc[3];
  for (int c = 0; c < 3; ++c) {
    // Zero-shift reduction of dphidrPBC mode 2 (bitwise-equal).
    ra[c] = -xyz[3 * i + c] + xyz[3 * j + c];
    rb[c] = -xyz[3 * j + c] + xyz[3 * k + c];
    rc[c] = -xyz[3 * k + c] + xyz[3 * l + c];
  }
  double rapb[3], rbpc[3], na[3], nb[3];
  for (int c = 0; c < 3; ++c) {
    rapb[c] = ra[c] + rb[c];
    rbpc[c] = rb[c] + rc[c];
  }
  vecCross(ra, rb, na);
  vecCross(rb, rc, nb);
  double nan = vecNorm(na, 0);
  double nbn = vecNorm(nb, 0);
  double cosphi = std::cos(phi), sinphi = std::sin(phi);
  double nenner = nan * nbn * sinphi;
  double onenner = 0.0;
  if (std::fabs(nenner) < 1.0e-14) {
    da[0] = da[1] = da[2] = 0.0;
    db[0] = db[1] = db[2] = 0.0;
    dc[0] = dc[1] = dc[2] = 0.0;
    dd[0] = dd[1] = dd[2] = 0.0;
    if (std::fabs(nan * nbn) > 1.0e-14)
      onenner = 1.0 / (nan * nbn);
    else
      onenner = 0.0;
  } else {
    onenner = 1.0 / nenner;
  }
  double rab[3], rba[3], rac[3], rbb[3], rbc[3], raa[3];
  double rapba[3], rapbb[3], rbpca[3], rbpcb[3];
  vecCross(na, rb, rab);
  vecCross(nb, ra, rba);
  vecCross(na, rc, rac);
  vecCross(nb, rb, rbb);
  vecCross(nb, rc, rbc);
  vecCross(na, ra, raa);
  vecCross(rapb, na, rapba);
  vecCross(rapb, nb, rapbb);
  vecCross(rbpc, na, rbpca);
  vecCross(rbpc, nb, rbpcb);
  if (std::fabs(onenner) > 1.0e-14) {
    for (int c = 0; c < 3; ++c) {
      da[c] = onenner * (cosphi * nbn / nan * rab[c] - rbb[c]);
      db[c] = onenner * (cosphi * (nbn / nan * rapba[c] + nan / nbn * rbc[c]) -
                         (rac[c] + rapbb[c]));
      dc[c] = onenner * (cosphi * (nbn / nan * raa[c] + nan / nbn * rbpcb[c]) -
                         (rba[c] + rbpca[c]));
      dd[c] = onenner * (cosphi * nan / nbn * rbb[c] - rab[c]);
    }
  } else {
    da[0] = da[1] = da[2] = 0.0;
    db[0] = db[1] = db[2] = 0.0;
    dc[0] = dc[1] = dc[2] = 0.0;
    dd[0] = dd[1] = dd[2] = 0.0;
  }
}

double inversionAngle(const std::vector<double>& xyz, int i, int j, int k,
                      int l)
{
  double re[3], rd[3], rv[3], rn[3];
  for (int c = 0; c < 3; ++c) {
    re[c] = xyz[3 * i + c] - xyz[3 * j + c];
    rd[c] = xyz[3 * k + c] - xyz[3 * j + c];
    rv[c] = xyz[3 * l + c] - xyz[3 * i + c];
  }
  vecCross(re, rd, rn);
  vecNorm(rn, 1);
  vecNorm(rv, 1);
  double rnv = rn[0] * rv[0] + rn[1] * rv[1] + rn[2] * rv[2];
  return std::asin(rnv);
}

void inversionDerivatives(const std::vector<double>& xyz, int i, int j, int k,
                          int l, double omega, double da[3], double db[3],
                          double dc[3], double dd[3])
{
  double re[3], rd[3], rv[3], rdme[3];
  for (int c = 0; c < 3; ++c) {
    re[c] = xyz[3 * i + c] - xyz[3 * j + c];
    rd[c] = xyz[3 * k + c] - xyz[3 * j + c];
    rv[c] = xyz[3 * l + c] - xyz[3 * i + c];
    rdme[c] = rd[c] - re[c];
  }
  double rn[3], rve[3], rne[3], rdv[3], rdn[3], rvdme[3], rndme[3];
  vecCross(re, rd, rn);
  double rvn = vecNorm(rv, 0);
  double rnn = vecNorm(rn, 0);
  vecCross(rv, re, rve);
  vecCross(rn, re, rne);
  vecCross(rd, rv, rdv);
  vecCross(rd, rn, rdn);
  vecCross(rv, rdme, rvdme);
  vecCross(rn, rdme, rndme);
  double sinomega = std::sin(omega);
  double nenner = rnn * rvn * std::cos(omega);
  if (std::fabs(nenner) > 1.0e-14) {
    double onenner = 1.0 / nenner;
    for (int c = 0; c < 3; ++c) {
      da[c] = onenner * (rdv[c] - rn[c] -
                         sinomega * (rvn / rnn * rdn[c] - rnn / rvn * rv[c]));
      db[c] = onenner * (rvdme[c] - sinomega * rvn / rnn * rndme[c]);
      dc[c] = onenner * (rve[c] - sinomega * rvn / rnn * rne[c]);
      dd[c] = onenner * (rn[c] - sinomega * rnn / rvn * rv[c]);
    }
  } else {
    for (int c = 0; c < 3; ++c) {
      da[c] = 0.0;
      db[c] = 0.0;
      dc[c] = 0.0;
      dd[c] = 0.0;
    }
  }
}

void bondEnergyGradient(int bond, int iat, int jat, double rab, double rij,
                        const std::vector<double>& drij,
                        const double drijdcn[2], int n,
                        const std::vector<double>& xyz,
                        const std::vector<double>& vbondShift,
                        const std::vector<double>& vbondSteep,
                        const std::vector<double>& vbondPref, double& energy,
                        std::vector<double>& gradient,
                        std::vector<double>& dEdcn, double sigma[3][3],
                        Environment& env)
{
  if (n <= 0) {
    env.error("empty bond energy setup", "bondEnergyGradient");
    return;
  }
  double t8 = vbondSteep[bond];
  double dr = rab - rij;
  double dum = vbondPref[bond] * std::exp(-(t8 * (dr * dr)));
  energy = energy + dum;
  double yy = 2.0 * t8 * dr * dum;
  double dx = -xyz[3 * jat] + xyz[3 * iat];
  double dy = -xyz[3 * jat + 1] + xyz[3 * iat + 1];
  double dz = -xyz[3 * jat + 2] + xyz[3 * iat + 2];
  double vrab[3] = { dx, dy, dz };
  double t4 = -yy * dx / rab;
  double t5 = -yy * dy / rab;
  double t6 = -yy * dz / rab;
  gradient[3 * iat] += t4;
  gradient[3 * iat + 1] += t5;
  gradient[3 * iat + 2] += t6;
  dEdcn[iat] = dEdcn[iat] + yy * drijdcn[0];
  sigma[0][0] += t4 * vrab[0];
  sigma[1][0] += t4 * vrab[1];
  sigma[2][0] += t4 * vrab[2];
  sigma[0][1] += t5 * vrab[0];
  sigma[1][1] += t5 * vrab[1];
  sigma[2][1] += t5 * vrab[2];
  sigma[0][2] += t6 * vrab[0];
  sigma[1][2] += t6 * vrab[1];
  sigma[2][2] += t6 * vrab[2];
  t4 = yy * (dx / rab);
  t5 = yy * (dy / rab);
  t6 = yy * (dz / rab);
  gradient[3 * jat] += t4;
  gradient[3 * jat + 1] += t5;
  gradient[3 * jat + 2] += t6;
  dEdcn[jat] = dEdcn[jat] + yy * drijdcn[1];
  for (int k = 0; k < n; ++k) {
    double dg0 = drij[3 * k] * yy;
    double dg1 = drij[3 * k + 1] * yy;
    double dg2 = drij[3 * k + 2] * yy;
    gradient[3 * k] += dg0;
    gradient[3 * k + 1] += dg1;
    gradient[3 * k + 2] += dg2;
    sigma[0][0] += dg0 * vrab[0];
    sigma[1][0] += dg0 * vrab[1];
    sigma[2][0] += dg0 * vrab[2];
    sigma[0][1] += dg1 * vrab[0];
    sigma[1][1] += dg1 * vrab[1];
    sigma[2][1] += dg1 * vrab[2];
    sigma[0][2] += dg2 * vrab[0];
    sigma[1][2] += dg2 * vrab[1];
    sigma[2][2] += dg2 * vrab[2];
  }
  (void)vbondShift;
}

void bendEnergyGradient(int center, int first, int second, double equilibrium,
                        double forceConstant,
                        const std::vector<int>& numbers,
                        const std::vector<double>& xyz,
                        const std::vector<double>& rcov, double angleCutA,
                        double& energy, double atomGrad[3][3],
                        double stress[3][3])
{
  double va[3] = { xyz[3 * first], xyz[3 * first + 1], xyz[3 * first + 2] };
  double vb[3] = { xyz[3 * center], xyz[3 * center + 1],
                   xyz[3 * center + 2] };
  double vc[3] = { xyz[3 * second], xyz[3 * second + 1],
                   xyz[3 * second + 2] };
  double vab[3], vcb[3], vp[3];
  vecSub(va, vb, vab);
  vecSub(vc, vb, vcb);
  double rab2 = vab[0] * vab[0] + vab[1] * vab[1] + vab[2] * vab[2];
  double rcb2 = vcb[0] * vcb[0] + vcb[1] * vcb[1] + vcb[2] * vcb[2];
  vecCross(vcb, vab, vp);
  double rp = vecLen(vp) + 1.0e-14;
  double cosa = vecCos(vab, vcb);
  cosa = std::min(1.0, std::max(-1.0, cosa));
  double theta = std::acos(cosa);
  double dampij = 0.0, damp2ij = 0.0, dampjk = 0.0, damp2jk = 0.0;
  bondDamp(rab2, angleCutA, rcov, numbers[first], numbers[center], dampij,
           damp2ij);
  bondDamp(rcb2, angleCutA, rcov, numbers[second], numbers[center], dampjk,
           damp2jk);
  double damp = dampij * dampjk;
  double c0 = equilibrium, kijk = forceConstant;
  double dt = 0.0, ea = 0.0, deddt = 0.0;
  if (pi - c0 < 1.0e-6) {
    dt = theta - c0;
    ea = kijk * (dt * dt);
    deddt = 2.0 * kijk * dt;
  } else {
    double diff = cosa - std::cos(c0);
    ea = kijk * (diff * diff);
    deddt = 2.0 * kijk * std::sin(theta) * (std::cos(c0) - cosa);
  }
  double e = ea * damp;
  energy = e;
  double deda[3], dedc[3], dedb[3];
  vecCross(vab, vp, deda);
  double rmul1 = -deddt / (rab2 * rp);
  for (int c = 0; c < 3; ++c)
    deda[c] = deda[c] * rmul1;
  vecCross(vcb, vp, dedc);
  double rmul2 = deddt / (rcb2 * rp);
  for (int c = 0; c < 3; ++c)
    dedc[c] = dedc[c] * rmul2;
  for (int c = 0; c < 3; ++c)
    dedb[c] = deda[c] + dedc[c];
  double term1[3], term2[3];
  for (int c = 0; c < 3; ++c) {
    term1[c] = ea * damp2ij * dampjk * vab[c];
    term2[c] = ea * damp2jk * dampij * vcb[c];
  }
  // Local columns map to (center, first, second) like the caller.
  for (int c = 0; c < 3; ++c) {
    atomGrad[0][c] = -dedb[c] * damp - term1[c] - term2[c];
    atomGrad[1][c] = deda[c] * damp + term1[c];
    atomGrad[2][c] = dedc[c] * damp + term2[c];
  }
  for (int dim1 = 0; dim1 < 3; ++dim1) {
    for (int dim2 = dim1; dim2 < 3; ++dim2) {
      stress[dim1][dim2] = atomGrad[0][dim2] * vb[dim1] +
                           atomGrad[1][dim2] * va[dim1] +
                           atomGrad[2][dim2] * vc[dim1];
    }
  }
  for (int dim1 = 0; dim1 < 3; ++dim1) {
    for (int dim2 = 0; dim2 < dim1; ++dim2)
      stress[dim1][dim2] = stress[dim2][dim1];
  }
}

void torsionEnergyGradient(int outer1, int center1, int center2, int outer2,
                           int kind, double phase, double forceConstant,
                           const std::vector<int>& numbers,
                           const std::vector<double>& xyz,
                           const std::vector<double>& rcov, double angleCutT,
                           bool periodic, double& energy,
                           double atomGrad[4][3], double stress[3][3])
{
  // egtors(m, i=outer1, j=center1, k=center2, l=outer2).
  int i = outer1, j = center1, k = center2, l = outer2;
  double rn = static_cast<double>(kind);
  double phi0 = phase;
  double e = 0.0;
  if (kind > 0) {
    double va[3] = { xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2] };
    double vb[3] = { xyz[3 * j], xyz[3 * j + 1], xyz[3 * j + 2] };
    double vc[3] = { xyz[3 * k], xyz[3 * k + 1], xyz[3 * k + 2] };
    double vd[3] = { xyz[3 * l], xyz[3 * l + 1], xyz[3 * l + 2] };
    double vab[3], vcb[3], vdc[3];
    vecSub(va, vb, vab);
    vecSub(vb, vc, vcb);
    vecSub(vc, vd, vdc);
    double rij = vab[0] * vab[0] + vab[1] * vab[1] + vab[2] * vab[2];
    double rjk = vcb[0] * vcb[0] + vcb[1] * vcb[1] + vcb[2] * vcb[2];
    double rkl = vdc[0] * vdc[0] + vdc[1] * vdc[1] + vdc[2] * vdc[2];
    double dampij = 0.0, damp2ij = 0.0, dampjk = 0.0, damp2jk = 0.0,
           dampkl = 0.0, damp2kl = 0.0;
    bondDamp(rij, angleCutT, rcov, numbers[i], numbers[j], dampij, damp2ij);
    bondDamp(rjk, angleCutT, rcov, numbers[k], numbers[j], dampjk, damp2jk);
    bondDamp(rkl, angleCutT, rcov, numbers[k], numbers[l], dampkl, damp2kl);
    double damp = dampjk * dampij * dampkl;
    double phi = dihedralAngle(xyz, i, j, k, l);
    double dda[3], ddb[3], ddc[3], ddd[3];
    dihedralDerivatives(xyz, i, j, k, l, phi, dda, ddb, ddc, ddd);
    double dphi1 = phi - phi0;
    double c1 = rn * dphi1 + pi;
    double x1cos = std::cos(c1), x1sin = std::sin(c1);
    double et = (1.0 + x1cos) * forceConstant;
    double dij = -rn * x1sin * forceConstant * damp;
    double term1[3], term2[3], term3[3];
    for (int c = 0; c < 3; ++c) {
      term1[c] = et * damp2ij * dampjk * dampkl * vab[c];
      term2[c] = et * damp2jk * dampij * dampkl * vcb[c];
      term3[c] = et * damp2kl * dampij * dampjk * vdc[c];
    }
    for (int c = 0; c < 3; ++c) {
      atomGrad[0][c] = dij * dda[c] + term1[c];
      atomGrad[1][c] = dij * ddb[c] - term1[c] + term2[c];
      atomGrad[2][c] = dij * ddc[c] + term3[c] - term2[c];
      atomGrad[3][c] = dij * ddd[c] - term3[c];
    }
    if (periodic) {
      for (int dim1 = 0; dim1 < 3; ++dim1) {
        for (int dim2 = dim1; dim2 < 3; ++dim2) {
          stress[dim1][dim2] = atomGrad[0][dim2] * va[dim1] +
                               atomGrad[1][dim2] * vb[dim1] +
                               atomGrad[2][dim2] * vc[dim1] +
                               atomGrad[3][dim2] * vd[dim1];
          stress[dim2][dim1] = stress[dim1][dim2];
        }
      }
      for (int dim1 = 0; dim1 < 3; ++dim1) {
        for (int dim2 = 0; dim2 < dim1; ++dim2)
          stress[dim1][dim2] = stress[dim2][dim1];
      }
    }
    e = et * damp;
  } else {
    double va[3] = { xyz[3 * i], xyz[3 * i + 1], xyz[3 * i + 2] };
    double vb[3] = { xyz[3 * j], xyz[3 * j + 1], xyz[3 * j + 2] };
    double vc[3] = { xyz[3 * k], xyz[3 * k + 1], xyz[3 * k + 2] };
    double vd[3] = { xyz[3 * l], xyz[3 * l + 1], xyz[3 * l + 2] };
    double vab[3], vcb[3], vdc[3];
    vecSub(vb, va, vab);
    vecSub(vb, vc, vcb);
    vecSub(vb, vd, vdc);
    double rij = vab[0] * vab[0] + vab[1] * vab[1] + vab[2] * vab[2];
    double rjk = vcb[0] * vcb[0] + vcb[1] * vcb[1] + vcb[2] * vcb[2];
    double rjl = vdc[0] * vdc[0] + vdc[1] * vdc[1] + vdc[2] * vdc[2];
    double dampij = 0.0, damp2ij = 0.0, dampjk = 0.0, damp2jk = 0.0,
           dampjl = 0.0, damp2jl = 0.0;
    bondDamp(rij, angleCutT, rcov, numbers[i], numbers[j], dampij, damp2ij);
    bondDamp(rjk, angleCutT, rcov, numbers[k], numbers[j], dampjk, damp2jk);
    bondDamp(rjl, angleCutT, rcov, numbers[j], numbers[l], dampjl, damp2jl);
    double damp = dampjk * dampij * dampjl;
    double phi = inversionAngle(xyz, i, j, k, l);
    double dda[3], ddb[3], ddc[3], ddd[3];
    inversionDerivatives(xyz, i, j, k, l, phi, dda, ddb, ddc, ddd);
    double et = 0.0, dij = 0.0;
    if (kind == 0) {
      double dphi1 = phi - phi0;
      double c1 = dphi1 + pi;
      double x1cos = std::cos(c1), x1sin = std::sin(c1);
      et = (1.0 + x1cos) * forceConstant;
      dij = -x1sin * forceConstant * damp;
    } else {
      double diff = std::cos(phi) - std::cos(phi0);
      et = forceConstant * (diff * diff);
      dij = 2.0 * forceConstant * std::sin(phi) * (std::cos(phi0) -
                                                  std::cos(phi)) *
            damp;
    }
    double term1[3], term2[3], term3[3];
    for (int c = 0; c < 3; ++c) {
      term1[c] = et * damp2ij * dampjk * dampjl * vab[c];
      term2[c] = et * damp2jk * dampij * dampjl * vcb[c];
      term3[c] = et * damp2jl * dampij * dampjk * vdc[c];
    }
    for (int c = 0; c < 3; ++c) {
      atomGrad[0][c] = dij * dda[c] - term1[c];
      atomGrad[1][c] = dij * ddb[c] + term1[c] + term2[c] + term3[c];
      atomGrad[2][c] = dij * ddc[c] - term2[c];
      atomGrad[3][c] = dij * ddd[c] - term3[c];
    }
    if (periodic) {
      for (int dim1 = 0; dim1 < 3; ++dim1) {
        for (int dim2 = dim1; dim2 < 3; ++dim2) {
          stress[dim1][dim2] = atomGrad[0][dim2] * va[dim1] +
                               atomGrad[1][dim2] * vb[dim1] +
                               atomGrad[2][dim2] * vc[dim1] +
                               atomGrad[3][dim2] * vd[dim1];
          stress[dim2][dim1] = stress[dim1][dim2];
        }
      }
      for (int dim1 = 0; dim1 < 3; ++dim1) {
        for (int dim2 = 0; dim2 < dim1; ++dim2)
          stress[dim1][dim2] = stress[dim2][dim1];
      }
    }
    e = et * damp;
  }
  energy = e;
}

bool bondLengthEstimates(int n, const std::vector<int>& numbers,
                         const std::vector<double>& cn,
                         const std::vector<double>& dlogCn,
                         const std::vector<std::pair<int, int>>& blist,
                         const std::vector<double>& rabIn,
                         std::vector<double>& rabOut,
                         std::vector<double>& grabOut,
                         std::vector<double>& rabdcnOut, Environment& env)
{
  int nbond = static_cast<int>(blist.size());
  if (n <= 0 || nbond <= 0) {
    env.error("empty bond length setup", "bondLengthEstimates");
    return false;
  }
  // p table and scale factors: unsuffixed literals, single widened on use
  // (same values as gfnffBondGuesses).
  const float p[6][2] = { { 29.84522887f, -8.87843763f },
                          { -1.70549806f, 2.10878369f },
                          { 6.54013762f, 0.08009374f },
                          { 6.39169003f, -0.85808076f },
                          { 6.00000000f, -1.15000000f },
                          { 5.60000000f, -1.30000000f } };
  float scaleF[103][103];
  for (int i = 0; i < 103; ++i) {
    for (int j = 0; j < 103; ++j)
      scaleF[i][j] = 1.0f;
  }
  for (int z = 57; z <= 71; ++z) {
    scaleF[z - 1][9 - 1] = 0.93533568f;
    scaleF[9 - 1][z - 1] = 0.93533568f;
    scaleF[z - 1][17 - 1] = 1.0190114f;
    scaleF[17 - 1][z - 1] = 1.0190114f;
    scaleF[z - 1][35 - 1] = 1.0425532f;
    scaleF[35 - 1][z - 1] = 1.0425532f;
    scaleF[z - 1][53 - 1] = 1.0551948f;
    scaleF[53 - 1][z - 1] = 1.0551948f;
  }
  for (int z = 89; z <= 103; ++z) {
    scaleF[z - 1][9 - 1] = 0.93533568f;
    scaleF[9 - 1][z - 1] = 0.93533568f;
    scaleF[z - 1][17 - 1] = 1.0190114f;
    scaleF[17 - 1][z - 1] = 1.0190114f;
    scaleF[z - 1][35 - 1] = 1.0425532f;
    scaleF[35 - 1][z - 1] = 1.0425532f;
    scaleF[z - 1][53 - 1] = 1.0551948f;
    scaleF[53 - 1][z - 1] = 1.0551948f;
  }
  rabOut.assign(nbond, 0.0);
  grabOut.assign(3 * n * nbond, 0.0);
  rabdcnOut.assign(2 * nbond, 0.0);
  for (int k = 0; k < nbond; ++k) {
    int j = blist[k].first;
    int i = blist[k].second;
    int ati = numbers[i], atj = numbers[j];
    int ir = elementRow6(ati) - 1, jr = elementRow6(atj) - 1;
    double ra = static_cast<double>(rabR0[ati - 1]) +
                static_cast<double>(rabCnfak[ati - 1]) * cn[i];
    double rb = static_cast<double>(rabR0[atj - 1]) +
                static_cast<double>(rabCnfak[atj - 1]) * cn[j];
    double den = std::abs(static_cast<double>(rabEn[ati - 1]) -
                          static_cast<double>(rabEn[atj - 1]));
    double k1 =
      0.005 * (static_cast<double>(p[ir][0]) + static_cast<double>(p[jr][0]));
    double k2 =
      0.005 * (static_cast<double>(p[ir][1]) + static_cast<double>(p[jr][1]));
    double den2 = den * den;
    double ff = (1.0 - (k1 * den)) - ((k2 * den2));
    double sf = static_cast<double>(scaleF[ati - 1][atj - 1]);
    rabOut[k] = (((ra + rb) + rabIn[k]) * ff) * sf;
    rabdcnOut[2 * k] =
      (static_cast<double>(rabCnfak[ati - 1]) * ff) * sf;
    rabdcnOut[2 * k + 1] =
      (static_cast<double>(rabCnfak[atj - 1]) * ff) * sf;
    double pre = (sf * ff);
    double ci = static_cast<double>(rabCnfak[ati - 1]);
    double cj = static_cast<double>(rabCnfak[atj - 1]);
    // dcn(c,m,i) = d(cn[i])/d(x[m]) lives at [(m*n+i)*3+c] (the Fortran
    // dlogCN(c,m,i) element order).
    for (int m = 0; m < n; ++m) {
      for (int c = 0; c < 3; ++c)
        grabOut[(k * n + m) * 3 + c] =
          (pre * (((ci * dlogCn[(m * n + i) * 3 + c]) +
                   ((cj * dlogCn[(m * n + j) * 3 + c])))));
    }
  }
  return true;
}

bool bondedRepulsion(int n, const std::vector<int>& numbers,
                     const std::vector<double>& xyz,
                     const std::vector<std::pair<int, int>>& blist,
                     const std::vector<double>& repa,
                     const std::vector<double>& repz, double repScaleB,
                     double& energy, std::vector<double>& gradient,
                     double sigma[3][3], Environment& env)
{
  int nbond = static_cast<int>(blist.size());
  if (n <= 0 || nbond <= 0) {
    env.error("empty bonded repulsion setup", "bondedRepulsion");
    return false;
  }
  for (const auto& pr : blist) {
    int jat = pr.first;
    int iat = pr.second;
    double dx = xyz[3 * iat] - xyz[3 * jat];
    double dy = xyz[3 * iat + 1] - xyz[3 * jat + 1];
    double dz = xyz[3 * iat + 2] - xyz[3 * jat + 2];
    double d[3] = { dx, dy, dz };
    double rab = norm2(d);
    double r2 = rab * rab;
    int ati = numbers[iat], atj = numbers[jat];
    double alpha = std::sqrt(repa[ati - 1] * repa[atj - 1]);
    double repab = (repz[ati - 1] * repz[atj - 1]) * repScaleB;
    double t16 = std::pow(r2, 0.75);
    double t19 = t16 * t16;
    double t26 = (std::exp(((-alpha) * t16))) * repab;
    energy = energy + (t26 / rab);
    double t27 = ((t26 * ((((1.5 * alpha) * t16)) + 1.0))) / t19;
    gradient[3 * iat] = gradient[3 * iat] - (dx * t27);
    gradient[3 * iat + 1] = gradient[3 * iat + 1] - (dy * t27);
    gradient[3 * iat + 2] = gradient[3 * iat + 2] - (dz * t27);
    gradient[3 * jat] = gradient[3 * jat] + (dx * t27);
    gradient[3 * jat + 1] = gradient[3 * jat + 1] + (dy * t27);
    gradient[3 * jat + 2] = gradient[3 * jat + 2] + (dz * t27);
    double sxx = ((dx * t27) * dx), sxy = ((dx * t27) * dy),
           sxz = ((dx * t27) * dz);
    double syx = ((dy * t27) * dx), syy = ((dy * t27) * dy),
           syz = ((dy * t27) * dz);
    double szx = ((dz * t27) * dx), szy = ((dz * t27) * dy),
           szz = ((dz * t27) * dz);
    sigma[0][0] = sigma[0][0] - sxx;
    sigma[0][1] = sigma[0][1] - sxy;
    sigma[0][2] = sigma[0][2] - sxz;
    sigma[1][0] = sigma[1][0] - syx;
    sigma[1][1] = sigma[1][1] - syy;
    sigma[1][2] = sigma[1][2] - syz;
    sigma[2][0] = sigma[2][0] - szx;
    sigma[2][1] = sigma[2][1] - szy;
    sigma[2][2] = sigma[2][2] - szz;
  }
  return true;
}

void batmEnergyGradient(int iat, int jat, int kat,
                        const std::vector<int>& numbers,
                        const std::vector<double>& xyz,
                        const std::vector<double>& charges,
                        const std::vector<double>& zb3atm, double& energy,
                        double atomGrad[3][3], Environment& env)
{
  int n = static_cast<int>(numbers.size());
  if (n <= 0) {
    env.error("empty ATM setup", "batmEnergyGradient");
    return;
  }
  for (int row = 0; row < 3; ++row) {
    for (int c = 0; c < 3; ++c)
      atomGrad[row][c] = 0.0;
  }
  energy = 0.0;
  double fi = 1.0 - ((3.0 * charges[iat]));
  fi = std::min(std::max(fi, -4.0), 4.0);
  double fj = 1.0 - ((3.0 * charges[jat]));
  fj = std::min(std::max(fj, -4.0), 4.0);
  double fk = 1.0 - ((3.0 * charges[kat]));
  fk = std::min(std::max(fk, -4.0), 4.0);
  double ff = (fi * fj) * fk;
  double c9 = (((ff * zb3atm[numbers[iat] - 1])) * zb3atm[numbers[jat] - 1]) *
              zb3atm[numbers[kat] - 1];
  double dxij = xyz[3 * iat] - xyz[3 * jat];
  double dyij = xyz[3 * iat + 1] - xyz[3 * jat + 1];
  double dzij = xyz[3 * iat + 2] - xyz[3 * jat + 2];
  double dxik = xyz[3 * iat] - xyz[3 * kat];
  double dyik = xyz[3 * iat + 1] - xyz[3 * kat + 1];
  double dzik = xyz[3 * iat + 2] - xyz[3 * kat + 2];
  double dxjk = xyz[3 * jat] - xyz[3 * kat];
  double dyjk = xyz[3 * jat + 1] - xyz[3 * kat + 1];
  double dzjk = xyz[3 * jat + 2] - xyz[3 * kat + 2];
  double r2ij = ((dxij * dxij) + (dyij * dyij)) + (dzij * dzij);
  double r2ik = ((dxik * dxik) + (dyik * dyik)) + (dzik * dzik);
  double r2jk = ((dxjk * dxjk) + (dyjk * dyjk)) + (dzjk * dzjk);
  double sr2ij = std::sqrt(r2ij);
  double sr2ik = std::sqrt(r2ik);
  double sr2jk = std::sqrt(r2jk);
  double invsr2ij = 1.0 / sr2ij;
  double invsr2ik = 1.0 / sr2ik;
  double invsr2jk = 1.0 / sr2jk;
  double mijk = ((-r2ij) + r2jk) + r2ik;
  double imjk = ((r2ij - r2jk)) + r2ik;
  double ijmk = ((r2ij + r2jk)) - r2ik;
  double rijk3 = (r2ij * r2jk) * r2ik;
  double rav3 = ((rijk3 * sr2ij) * sr2jk) * sr2ik;
  double ang = ((((0.375 * ijmk)) * imjk) * mijk) / rijk3;
  double angr9 = (ang + 1.0) / rav3;
  energy = c9 * angr9;
  // dang chains: -0.375*(A+B+C-D)/((sr2*rijk3)*rav3) per pair.
  double r2ij2 = r2ij * r2ij, r2jk2 = r2jk * r2jk, r2ik2 = r2ik * r2ik;
  double aij = ((r2ij2 * r2ij)) + ((r2ij2 * (r2jk + r2ik))) +
               ((r2ij * ((((3.0 * r2jk2)) + (((2.0 * r2jk) * r2ik))) +
                         ((3.0 * r2ik2))))) -
               ((((5.0 * (((r2jk - r2ik)) * ((r2jk - r2ik)))))) * ((r2jk + r2ik)));
  double dangij = (((-0.375) * aij)) / (((sr2ij * rijk3)) * rav3);
  double drij = ((-dangij)) * c9;
  double ajk = ((r2jk2 * r2jk)) + ((r2jk2 * (r2ik + r2ij))) +
               ((r2jk * ((((3.0 * r2ik2)) + (((2.0 * r2ik) * r2ij))) +
                         ((3.0 * r2ij2))))) -
               ((((5.0 * (((r2ik - r2ij)) * ((r2ik - r2ij)))))) * ((r2ik + r2ij)));
  double dangjk = (((-0.375) * ajk)) / (((sr2jk * rijk3)) * rav3);
  double drjk = ((-dangjk)) * c9;
  double aik = ((r2ik2 * r2ik)) + ((r2ik2 * (r2jk + r2ij))) +
               ((r2ik * ((((3.0 * r2jk2)) + (((2.0 * r2jk) * r2ij))) +
                         ((3.0 * r2ij2))))) -
               ((((5.0 * (((r2jk - r2ij)) * ((r2jk - r2ij)))))) * ((r2jk + r2ij)));
  double dangik = (((-0.375) * aik)) / (((sr2ik * rijk3)) * rav3);
  double drik = ((-dangik)) * c9;
  double rij[3], rik[3], rjk[3];
  // Reference vectors (0d image branch): jat-iat, kat-iat, kat-jat.
  for (int c = 0; c < 3; ++c) {
    rij[c] = xyz[3 * jat + c] - xyz[3 * iat + c];
    rik[c] = xyz[3 * kat + c] - xyz[3 * iat + c];
    rjk[c] = xyz[3 * kat + c] - xyz[3 * jat + c];
  }
  for (int c = 0; c < 3; ++c) {
    atomGrad[0][c] = (((drij * rij[c]) * invsr2ij)) +
                     (((drik * rik[c]) * invsr2ik));
    atomGrad[1][c] = (((drjk * rjk[c]) * invsr2jk)) -
                     (((drij * rij[c]) * invsr2ij));
    atomGrad[2][c] = ((((-drik) * rik[c]) * invsr2ik)) -
                     (((drjk * rjk[c]) * invsr2jk));
  }
}

} // namespace Xtb
} // namespace Avogadro
