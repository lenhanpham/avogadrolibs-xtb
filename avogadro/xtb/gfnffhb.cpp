/******************************************************************************
  This source file is part of the Avogadro project.
  This source code is released under the GNU Lesser General Public License,
  see "COPYING.LESSER" in this directory.

  Ported from xtb (https://github.com/grimme-lab/xtb), originally
  src/gfnff/gfnff_eg.f90, src/constr.f90 and src/basic_geo.f90,
  Copyright (C) 2019-2020 Sebastian Ehlert.
******************************************************************************/

#include "gfnffhb.h"

#include "constants.h"
#include "eeq.h"
#include "environment.h"
#include "gfnffegbond.h"

#include <cmath>
#include <vector>

namespace Avogadro {
namespace Xtb {

void hbStrengths(int a, int b, const std::vector<double>& hbBas,
                 const std::vector<double>& hbAci, double* ca, double* cb)
{
  ca[0] = hbBas[a];
  cb[0] = hbBas[b];
  ca[1] = hbAci[a];
  cb[1] = hbAci[b];
}

void crossProd(const double* a, const double* b, double* c)
{
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

double cosAngle(const double* a, const double* b)
{
  double rimp = a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  double da[3] = { a[0], a[1], a[2] };
  double db[3] = { b[0], b[1], b[2] };
  double al = norm2(da);
  double bl = norm2(db);
  if (al > 0.0 && bl > 0.0)
    return rimp / (al * bl);
  return 0.0;
}

double vecNorm10(double* v, int normalise)
{
  double sp = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
  double rn = std::sqrt(sp);
  if (normalise > 0 && std::abs(rn) > 1.0e-14) {
    double or_ = 1.0 / rn;
    v[0] = or_ * v[0];
    v[1] = or_ * v[1];
    v[2] = or_ * v[2];
  }
  return rn;
}

namespace {

inline void sub3(const std::vector<double>& xyz, int a, int b, double* out)
{
  out[0] = xyz[3 * a] - xyz[3 * b];
  out[1] = xyz[3 * a + 1] - xyz[3 * b + 1];
  out[2] = xyz[3 * a + 2] - xyz[3 * b + 2];
}

inline double norm3(const double* v)
{
  double d[3] = { v[0], v[1], v[2] };
  return norm2(d);
}

inline double sumSq(const double* v)
{
  return v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
}

inline double chargeFactor(double q, double st, double sf, bool donor)
{
  double ex1 = std::exp(donor ? -st * q : st * q);
  return ex1 / (ex1 + sf);
}

} // namespace

double dihedralValue0d(int r, int b, int c, int h,
                       const std::vector<double>& xyz)
{
  // valijklffPBC mode 2 with zero shifts: ra = B-R, rb = C-B, rc = H-C.
  // The deter/thab/thbc/valijkPBC evaluations are dead (their results
  // never enter acos(snanb)) and are skipped; vecnorm(.,3,1) normalises
  // na/nb in place when their norm exceeds 1e-14 (verbatim side effect).
  double ra[3], rb[3], rc[3], na[3], nb[3];
  sub3(xyz, b, r, ra);
  sub3(xyz, c, b, rb);
  sub3(xyz, h, c, rc);
  crossProd(ra, rb, na);
  crossProd(rb, rc, nb);
  vecNorm10(na, 1);
  vecNorm10(nb, 1);
  double snanb = na[0] * nb[0] + na[1] * nb[1] + na[2] * nb[2];
  if (std::abs(std::abs(snanb) - 1.0) < 1.0e-14)
    snanb = std::copysign(1.0, snanb);
  return std::acos(snanb);
}

void dihedralDerivs0d(int r, int b, int c, int h, double phi,
                      const std::vector<double>& xyz, double* dda,
                      double* ddb, double* ddc, double* ddd)
{
  // dphidrPBC mode 1 with zero shifts.
  const double eps = 1.0e-14;
  double cosphi = std::cos(phi);
  double ra[3], rb[3], rc[3], rapb[3], rbpc[3];
  sub3(xyz, b, r, ra);
  sub3(xyz, c, b, rb);
  sub3(xyz, h, c, rc);
  for (int k = 0; k < 3; ++k) {
    rapb[k] = ra[k] + rb[k];
    rbpc[k] = rb[k] + rc[k];
  }
  double na[3], nb[3];
  crossProd(ra, rb, na);
  crossProd(rb, rc, nb);
  double nan = vecNorm10(na, 0);
  double nbn = vecNorm10(nb, 0);
  double nenner = nan * nbn * std::sin(phi);
  double onenner;
  if (std::abs(nenner) < eps) {
    for (int k = 0; k < 3; ++k)
      dda[k] = ddb[k] = ddc[k] = ddd[k] = 0.0;
    if (std::abs(nan * nbn) > eps)
      onenner = 1.0 / (nan * nbn);
    else
      onenner = 0.0;
  } else {
    onenner = 1.0 / nenner;
  }
  double rab[3], rba[3], rac[3], rbb[3], rbc[3], raa[3];
  double rapba[3], rapbb[3], rbpca[3], rbpcb[3];
  crossProd(na, rb, rab);
  crossProd(nb, ra, rba);
  crossProd(na, rc, rac);
  crossProd(nb, rb, rbb);
  crossProd(nb, rc, rbc);
  crossProd(na, ra, raa);
  crossProd(rapb, na, rapba);
  crossProd(rapb, nb, rapbb);
  crossProd(rbpc, na, rbpca);
  crossProd(rbpc, nb, rbpcb);
  if (std::abs(onenner) > eps) {
    for (int k = 0; k < 3; ++k) {
      dda[k] = onenner * ((((cosphi * (nbn / nan)) * rab[k])) - rbb[k]);
      ddb[k] = onenner * ((cosphi * ((((nbn / nan)) * rapba[k]) +
                                     (((nan / nbn)) * rbc[k]))) -
                          (rac[k] + rapbb[k]));
      ddc[k] = onenner * ((cosphi * ((((nbn / nan)) * raa[k]) +
                                     (((nan / nbn)) * rbpcb[k]))) -
                          (rba[k] + rbpca[k]));
      ddd[k] = onenner * ((((cosphi * (nan / nbn)) * rbb[k])) - rab[k]);
    }
  } else {
    for (int k = 0; k < 3; ++k)
      dda[k] = ddb[k] = ddc[k] = ddd[k] = 0.0;
  }
}

void bendNciMul(int b, int c, int h, double c0, double fc,
                const std::vector<double>& xyz, double& e, double g[3][3])
{
  // egbend_nci_mul for (j,i,k) = (B,C,H) with zero shifts.
  double c0c = std::cos(c0);
  double dc = 1.0 - c0c;
  double kijk = fc / (dc * dc);
  double va[3], vb[3], vc[3], vab[3], vcb[3];
  va[0] = xyz[3 * c];
  va[1] = xyz[3 * c + 1];
  va[2] = xyz[3 * c + 2];
  vb[0] = xyz[3 * b];
  vb[1] = xyz[3 * b + 1];
  vb[2] = xyz[3 * b + 2];
  vc[0] = xyz[3 * h];
  vc[1] = xyz[3 * h + 1];
  vc[2] = xyz[3 * h + 2];
  for (int k = 0; k < 3; ++k) {
    vab[k] = va[k] - vb[k];
    vcb[k] = vc[k] - vb[k];
  }
  double rab2 = sumSq(vab);
  double rcb2 = sumSq(vcb);
  double vp[3];
  crossProd(vcb, vab, vp);
  double tot = sumSq(vp);
  double rp = (tot > 0.0 ? std::sqrt(tot) : 0.0) + 1.0e-14;
  double cosa = cosAngle(vab, vcb);
  if (cosa > 1.0)
    cosa = 1.0;
  if (cosa < -1.0)
    cosa = -1.0;
  double theta = std::acos(cosa);
  double dt, ea, deddt;
  if (pi - c0 < 1.0e-6) {
    dt = theta - c0;
    ea = kijk * (dt * dt);
    deddt = (2.0 * kijk) * dt;
  } else {
    double dca = cosa - c0c;
    ea = kijk * (dca * dca);
    deddt = (((2.0 * kijk) * std::sin(theta))) * (c0c - cosa);
  }
  e = 1.0 - ea;
  double deda[3], dedc[3], dedb[3];
  crossProd(vab, vp, deda);
  double rmul1 = ((-deddt) / (rab2 * rp));
  for (int k = 0; k < 3; ++k)
    deda[k] = deda[k] * rmul1;
  crossProd(vcb, vp, dedc);
  double rmul2 = ((deddt) / (rcb2 * rp));
  for (int k = 0; k < 3; ++k)
    dedc[k] = dedc[k] * rmul2;
  for (int k = 0; k < 3; ++k)
    dedb[k] = deda[k] + dedc[k];
  for (int k = 0; k < 3; ++k) {
    g[k][0] = dedb[k];
    g[k][1] = -deda[k];
    g[k][2] = -dedc[k];
  }
}

void torsNciMul(int r, int b, int c, int h, int rn, double phi, double phi0,
                double tshift, const std::vector<double>& xyz, double& e,
                double g[3][4])
{
  // egtors_nci_mul with zero shifts (rij/rjk/rkl are dead and skipped).
  double fc = (1.0 - tshift) / 2.0;
  double dda[3], ddb[3], ddc[3], ddd[3];
  dihedralDerivs0d(r, b, c, h, phi, xyz, dda, ddb, ddc, ddd);
  double dphi1 = phi - phi0;
  double c1 = ((rn * dphi1) + pi);
  double x1cos = std::cos(c1);
  double x1sin = std::sin(c1);
  double et = (((1.0 + x1cos) * fc) + tshift);
  double dij = ((((-rn) * x1sin)) * fc);
  for (int k = 0; k < 3; ++k) {
    g[k][0] = dij * dda[k];
    g[k][1] = dij * ddb[k];
    g[k][2] = dij * ddc[k];
    g[k][3] = dij * ddd[k];
  }
  e = et;
}

bool abhEg1(int a, int b, int h, const std::vector<int>& numbers,
            const std::vector<double>& xyz, const std::vector<double>& qa,
            const std::vector<double>& hbBas, const std::vector<double>& hbAci,
            const std::vector<double>& rad, const HbParams& par, double mcf,
            double& energy, double gdr[3][3], double sigma[3][3],
            Environment& env)
{
  int n = static_cast<int>(numbers.size());
  if (n <= 0 || a < 0 || b < 0 || h < 0 || a >= n || b >= n || h >= n) {
    env.error("empty HB setup", "abhEg1");
    return false;
  }
  for (int c = 0; c < 3; ++c) {
    gdr[c][0] = 0.0;
    gdr[c][1] = 0.0;
    gdr[c][2] = 0.0;
  }
  energy = 0.0;
  double ca[2], cb[2];
  hbStrengths(a, b, hbBas, hbAci, ca, cb);

  double drah[3], drbh[3], drab[3];
  sub3(xyz, a, h, drah);
  sub3(xyz, b, h, drbh);
  sub3(xyz, a, b, drab);
  double rab = norm3(drab);
  double rab2 = rab * rab;
  double rah = norm3(drah);
  double rah2 = rah * rah;
  double rbh = norm3(drbh);
  double rbh2 = rbh * rbh;

  double rahprbh = rah + rbh + 1.0e-12;
  double radab = rad[numbers[a] - 1] + rad[numbers[b] - 1];

  double expo = (par.hbacut / radab) * (rahprbh / rab - 1.0);
  if (expo > 15.0)
    return true;
  double ratio2 = std::exp(expo);
  double outl = 2.0 / (1.0 + ratio2);

  double ratio1 = std::pow(rab2 / par.hblongcut, par.hbalp);
  double dampl = 1.0 / (1.0 + ratio1);
  double shortcut = par.hbscut * radab;
  double ratio3 = std::pow(shortcut / rab2, par.hbalp);
  double damps = 1.0 / (1.0 + ratio3);
  double damp = damps * dampl;
  double rdamp = (damp / rab2) / rab;

  double qh = chargeFactor(qa[h], par.hbst, par.hbsf, false);
  double qaa = chargeFactor(qa[a], par.hbst, par.hbsf, true);
  double qb = chargeFactor(qa[b], par.hbst, par.hbsf, true);

  double rah4 = rah2 * rah2;
  double rbh4 = rbh2 * rbh2;
  double denom = 1.0 / (rah4 + rbh4);
  double caa = qaa * ca[0];
  double cbb = qb * cb[0];
  double qhoutl = qh * outl;
  double bas = ((caa * rah4) + (cbb * rbh4)) * denom;
  double aci = ((cb[1] * rah4) + (ca[1] * rbh4)) * denom;

  double rterm = (((-aci) * rdamp) * qhoutl);
  energy = bas * rterm;

  double aterm = (((( -aci) * bas) * rdamp) * qh);
  double sterm = (((-rdamp) * bas) * qhoutl);
  double dterm = (((-aci) * bas) * qhoutl);
  double tmp = ((denom * denom) * 4.0);
  double dd24a = ((rah2 * rbh4) * tmp);
  double dd24b = ((rbh2 * rah4) * tmp);

  double ga[3], gb[3], gh[3], dg[3], dga[3], dgb[3], dgh[3], gi;
  gi = (((caa - cbb) * dd24a) * rterm);
  for (int c = 0; c < 3; ++c)
    ga[c] = gi * drah[c];
  gi = (((cbb - caa) * dd24b) * rterm);
  for (int c = 0; c < 3; ++c)
    gb[c] = gi * drbh[c];
  for (int c = 0; c < 3; ++c)
    gh[c] = -ga[c] - gb[c];

  gi = ((cb[1] - ca[1]) * dd24a);
  for (int c = 0; c < 3; ++c)
    dga[c] = (gi * drah[c]) * sterm;
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dga[c];
  gi = ((ca[1] - cb[1]) * dd24b);
  for (int c = 0; c < 3; ++c)
    dgb[c] = (gi * drbh[c]) * sterm;
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] + dgb[c];
  for (int c = 0; c < 3; ++c)
    dgh[c] = -dga[c] - dgb[c];
  for (int c = 0; c < 3; ++c)
    gh[c] = gh[c] + dgh[c];

  double dA = (((2.0 * par.hbalp) * ratio1) / (1.0 + ratio1));
  double dB = (((2.0 * par.hbalp) * ratio3) / (1.0 + ratio3));
  double o2 = 1.0 + ratio2;
  gi = (((rdamp * (((-dA) + dB) - 3.0))) / rab2);
  gi = gi * dterm;
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drab[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dg[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] - dg[c];

  gi = ((((((((aterm * 2.0) * ratio2) * expo) * rahprbh) / (o2 * o2)) /
         (rahprbh - rab))) /
        rab2);
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drab[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dg[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] - dg[c];

  tmp = ((((((-2.0 * aterm) * ratio2) * expo) / (o2 * o2)) / (rahprbh - rab)));
  for (int c = 0; c < 3; ++c)
    dga[c] = (drah[c] * tmp) / rah;
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dga[c];
  for (int c = 0; c < 3; ++c)
    dgb[c] = (drbh[c] * tmp) / rbh;
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] + dgb[c];
  for (int c = 0; c < 3; ++c)
    dgh[c] = -dga[c] - dgb[c];
  for (int c = 0; c < 3; ++c)
    gh[c] = gh[c] + dgh[c];

  // sigma columns for A, B, H (shifts zero in 0d): each column c gets
  // mcf*gX[c]*xyz[X](:) outer products, verbatim.
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      sigma[r][c] = sigma[r][c] + ((mcf * ga[c]) * xyz[3 * a + r]);
      sigma[r][c] = sigma[r][c] + ((mcf * gb[c]) * xyz[3 * b + r]);
      sigma[r][c] = sigma[r][c] + ((mcf * gh[c]) * xyz[3 * h + r]);
    }
  }
  for (int c = 0; c < 3; ++c) {
    gdr[c][0] = ga[c];
    gdr[c][1] = gb[c];
    gdr[c][2] = gh[c];
  }
  return true;
}

bool abhEg2new(int a, int b, int h, const std::vector<int>& nbrsB,
               const std::vector<int>& numbers,
               const std::vector<double>& xyz,
               const std::vector<double>& qa,
               const std::vector<double>& hbBas,
               const std::vector<double>& hbAci,
               const std::vector<double>& rad, const HbParams& par,
               double mcf, double& energy, std::vector<double>& gdr,
               double sigma[3][3], Environment& env)
{
  int n = static_cast<int>(numbers.size());
  int nbb = static_cast<int>(nbrsB.size());
  if (n <= 0 || a < 0 || b < 0 || h < 0 || a >= n || b >= n || h >= n) {
    env.error("empty HB setup", "abhEg2new");
    return false;
  }
  double pBh = 1.0 + par.hbabmix;
  double pAb = -par.hbabmix;
  energy = 0.0;
  double ca[2], cb[2];
  hbStrengths(a, b, hbBas, hbAci, ca, cb);

  std::vector<double> dranb(3 * nbb, 0.0), drbnb(3 * nbb, 0.0);
  std::vector<double> ranb(nbb, 0.0), ranb2(nbb, 0.0), rbnb(nbb, 0.0),
    rbnb2(nbb, 0.0);
  for (int i = 0; i < nbb; ++i) {
    int inb = nbrsB[i];
    for (int c = 0; c < 3; ++c) {
      dranb[i * 3 + c] = xyz[3 * a + c] - xyz[3 * inb + c];
      drbnb[i * 3 + c] = xyz[3 * b + c] - xyz[3 * inb + c];
    }
    ranb2[i] = sumSq(&dranb[i * 3]);
    ranb[i] = std::sqrt(ranb2[i]);
    rbnb2[i] = sumSq(&drbnb[i * 3]);
    rbnb[i] = std::sqrt(rbnb2[i]);
  }

  double drab[3], drah[3], drbh[3];
  sub3(xyz, a, b, drab);
  sub3(xyz, a, h, drah);
  sub3(xyz, b, h, drbh);
  double rab = norm3(drab);
  double rab2 = rab * rab;
  double rah = norm3(drah);
  double rah2 = rah * rah;
  (void)rah2; // computed like the reference; unused below
  double rbh = norm3(drbh);
  double rbh2 = rbh * rbh;
  double rahprbh = rah + rbh + 1.0e-12;
  double radab = rad[numbers[a] - 1] + rad[numbers[b] - 1];

  double expo = (par.hbacut / radab) * (rahprbh / rab - 1.0);
  if (expo > 15.0)
    return true;
  double ratio2 = std::exp(expo);
  double outl = 2.0 / (1.0 + ratio2);

  double hbnbcutSave = (numbers[b] == 7 && nbb == 1) ? 2.0 : par.hbnbcut;
  std::vector<double> ranbprbnb(nbb), expoNb(nbb), ratio2Nb(nbb), outlNb(nbb);
  for (int i = 0; i < nbb; ++i) {
    ranbprbnb[i] = ranb[i] + rbnb[i] + 1.0e-12;
    expoNb[i] = (hbnbcutSave / radab) * (ranbprbnb[i] / rab - 1.0);
    ratio2Nb[i] = std::pow(std::exp(-expoNb[i]), 1.0);
    outlNb[i] = (2.0 / (1.0 + ratio2Nb[i])) - 1.0;
  }
  double outlNbTot = 1.0;
  for (int i = 0; i < nbb; ++i)
    outlNbTot = outlNbTot * outlNb[i];

  double ratio1 = std::pow(rab2 / par.hblongcut, par.hbalp);
  double dampl = 1.0 / (1.0 + ratio1);
  double shortcut = par.hbscut * radab;
  double ratio3 = std::pow(shortcut / rab2, par.hbalp);
  double damps = 1.0 / (1.0 + ratio3);
  double damp = damps * dampl;
  double ddamp = ((((-2.0 * par.hbalp) * ratio1) / (1.0 + ratio1))) +
                 ((((2.0 * par.hbalp) * ratio3) / (1.0 + ratio3)));
  double rbhdamp = damp * (((pBh / rbh2) / rbh));
  double rabdamp = damp * (((pAb / rab2) / rab));
  double rdamp = rbhdamp + rabdamp;

  double qh = chargeFactor(qa[h], par.hbst, par.hbsf, false);
  double qaa = chargeFactor(qa[a], par.hbst, par.hbsf, true);
  double qb = chargeFactor(qa[b], par.hbst, par.hbsf, true);
  double qhoutl = ((qh * outl) * outlNbTot);
  double const_ = ((((ca[1] * qaa) * cb[0]) * qb) * par.xhaciGlobAbh);
  energy = (((-rdamp) * qhoutl) * const_);

  double aterm = (((( -rdamp) * qh) * outlNbTot) * const_);
  double nbterm = (((( -rdamp) * qh) * outl) * const_);
  double dterm = ((-qhoutl) * const_);

  double ga[3] = { 0.0, 0.0, 0.0 }, gb[3] = { 0.0, 0.0, 0.0 },
         gh[3] = { 0.0, 0.0, 0.0 };
  double dg[3], dga[3], dgb[3], dgh[3], gi;
  gi = ((((rabdamp + rbhdamp) * ddamp) - ((3.0) * rabdamp)) / rab2);
  gi = gi * dterm;
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drab[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = dg[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = -dg[c];

  gi = ((((-3.0) * rbhdamp)) / rbh2);
  gi = gi * dterm;
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drbh[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] + dg[c];
  for (int c = 0; c < 3; ++c)
    gh[c] = -dg[c];

  double o2 = 1.0 + ratio2;
  double tmp1 = (((((((-2.0 * aterm) * ratio2) * expo) / (o2 * o2))) /
                 (rahprbh - rab)));
  gi = (((-tmp1) * rahprbh) / rab2);
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drab[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dg[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] - dg[c];
  gi = (tmp1 / rah);
  for (int c = 0; c < 3; ++c)
    dga[c] = gi * drah[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dga[c];
  gi = (tmp1 / rbh);
  for (int c = 0; c < 3; ++c)
    dgb[c] = gi * drbh[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] + dgb[c];
  for (int c = 0; c < 3; ++c)
    dgh[c] = -dga[c] - dgb[c];
  for (int c = 0; c < 3; ++c)
    gh[c] = gh[c] + dgh[c];

  std::vector<double> tmp2(nbb, 0.0), giNb(nbb, 0.0);
  std::vector<double> gnb(3 * nbb, 0.0), dgnb(3, 0.0);
  dgnb.assign(3, 0.0);
  std::vector<char> mask(nbb, 1);
  for (int i = 0; i < nbb; ++i) {
    mask[i] = 0;
    double prod = 1.0;
    for (int k = 0; k < nbb; ++k) {
      if (mask[k])
        prod = prod * outlNb[k];
    }
    double o2nb = 1.0 + ratio2Nb[i];
    tmp2[i] = ((((((((2.0 * nbterm) * prod) * ratio2Nb[i]) * expoNb[i]) /
                  (o2nb * o2nb))) /
                (ranbprbnb[i] - rab)));
    giNb[i] = (((-tmp2[i]) * ranbprbnb[i]) / rab2);
    for (int c = 0; c < 3; ++c)
      dg[c] = giNb[i] * drab[c];
    for (int c = 0; c < 3; ++c)
      ga[c] = ga[c] + dg[c];
    for (int c = 0; c < 3; ++c)
      gb[c] = gb[c] - dg[c];
    mask[i] = 1;
  }
  for (int i = 0; i < nbb; ++i) {
    giNb[i] = (tmp2[i] / ranb[i]);
    for (int c = 0; c < 3; ++c)
      dga[c] = giNb[i] * dranb[i * 3 + c];
    for (int c = 0; c < 3; ++c)
      ga[c] = ga[c] + dga[c];
    giNb[i] = (tmp2[i] / rbnb[i]);
    for (int c = 0; c < 3; ++c)
      dgb[c] = giNb[i] * drbnb[i * 3 + c];
    for (int c = 0; c < 3; ++c)
      gb[c] = gb[c] + dgb[c];
    for (int c = 0; c < 3; ++c)
      dgnb[c] = -dga[c] - dgb[c];
    for (int c = 0; c < 3; ++c)
      gnb[i * 3 + c] = dgnb[c];
  }

  if (nbb < 1) {
    for (int c = 0; c < 3; ++c) {
      gdr[3 * a + c] = gdr[3 * a + c] + ga[c];
      gdr[3 * b + c] = gdr[3 * b + c] + gb[c];
      gdr[3 * h + c] = gdr[3 * h + c] + gh[c];
    }
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c) {
        sigma[r][c] = sigma[r][c] + ((mcf * ga[c]) * xyz[3 * a + r]);
        sigma[r][c] = sigma[r][c] + ((mcf * gb[c]) * xyz[3 * b + r]);
        sigma[r][c] = sigma[r][c] + ((mcf * gh[c]) * xyz[3 * h + r]);
      }
    }
    return true;
  }

  for (int c = 0; c < 3; ++c) {
    gdr[3 * a + c] = gdr[3 * a + c] + ga[c];
    gdr[3 * b + c] = gdr[3 * b + c] + gb[c];
    gdr[3 * h + c] = gdr[3 * h + c] + gh[c];
  }
  for (int i = 0; i < nbb; ++i) {
    int inb = nbrsB[i];
    for (int c = 0; c < 3; ++c)
      gdr[3 * inb + c] = gdr[3 * inb + c] + gnb[i * 3 + c];
  }
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      sigma[r][c] = sigma[r][c] + ((mcf * ga[c]) * xyz[3 * a + r]);
      sigma[r][c] = sigma[r][c] + ((mcf * gb[c]) * xyz[3 * b + r]);
      sigma[r][c] = sigma[r][c] + ((mcf * gh[c]) * xyz[3 * h + r]);
    }
  }
  for (int i = 0; i < nbb; ++i) {
    int inb = nbrsB[i];
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c)
        sigma[r][c] =
          sigma[r][c] + ((mcf * gnb[i * 3 + c]) * xyz[3 * inb + r]);
    }
  }
  return true;
}

bool abhEg2rnr(int a, int b, int h, const std::vector<int>& nbrsB,
               const std::vector<int>& numbers,
               const std::vector<double>& xyz,
               const std::vector<double>& qa,
               const std::vector<double>& hbBas,
               const std::vector<double>& hbAci,
               const std::vector<double>& rad,
               const std::vector<double>& repz, const HbParams& par,
               double mcf, double& energy, std::vector<double>& gdr,
               double sigma[3][3], Environment& env)
{
  int n = static_cast<int>(numbers.size());
  int nbb = static_cast<int>(nbrsB.size());
  if (n <= 0 || a < 0 || b < 0 || h < 0 || a >= n || b >= n || h >= n ||
      nbb != 2) {
    env.error("rNR HB needs exactly two B neighbours", "abhEg2rnr");
    return false;
  }
  double pBh = 1.0 + par.hbabmix;
  double pAb = -par.hbabmix;
  energy = 0.0;
  double ca[2], cb[2];
  hbStrengths(a, b, hbBas, hbAci, ca, cb);
  double lpDist = 0.50 - static_cast<double>(0.018f) * repz[numbers[b] - 1];
  const double hblpcut = 56.0;

  double dranb[2][3], drbnb[2][3], drnb[2][3];
  double ranb[2], rbnb[2];
  double vector[3] = { 0.0, 0.0, 0.0 };
  for (int i = 0; i < 2; ++i) {
    int inb = nbrsB[i];
    for (int c = 0; c < 3; ++c) {
      dranb[i][c] = xyz[3 * a + c] - xyz[3 * inb + c];
      drbnb[i][c] = xyz[3 * b + c] - xyz[3 * inb + c];
      drnb[i][c] = xyz[3 * inb + c] - xyz[3 * b + c];
      vector[c] = vector[c] + drnb[i][c];
    }
    ranb[i] = std::sqrt(sumSq(dranb[i]));
    rbnb[i] = std::sqrt(sumSq(drbnb[i]));
  }
  double vnorm = norm3(vector);
  double lp[3];
  if (vnorm > 1.0e-10) {
    for (int c = 0; c < 3; ++c)
      lp[c] = xyz[3 * b + c] - lpDist * (vector[c] / vnorm);
  } else {
    for (int c = 0; c < 3; ++c)
      lp[c] = xyz[3 * b + c];
    nbb = 0;
  }
  double drab[3], drah[3], drbh[3], dralp[3], drblp[3];
  sub3(xyz, a, b, drab);
  sub3(xyz, a, h, drah);
  sub3(xyz, b, h, drbh);
  for (int c = 0; c < 3; ++c) {
    dralp[c] = xyz[3 * a + c] - lp[c];
    drblp[c] = xyz[3 * b + c] - lp[c];
  }
  double rab = norm3(drab);
  double rab2 = rab * rab;
  double rah = norm3(drah);
  double rbh = norm3(drbh);
  double rbh2 = rbh * rbh;
  double rahprbh = rah + rbh + 1.0e-12;
  double radab = rad[numbers[a] - 1] + rad[numbers[b] - 1];

  double expo = (par.hbacut / radab) * (rahprbh / rab - 1.0);
  if (expo > 15.0)
    return true;
  double ratio2 = std::exp(expo);
  double outl = 2.0 / (1.0 + ratio2);

  double rblp2 = sumSq(drblp);
  double rblp = std::sqrt(rblp2);
  double ralp2 = sumSq(dralp);
  double ralp = std::sqrt(ralp2);
  double ralpprblp = ralp + rblp + 1.0e-12;
  double expoLp = (hblpcut / radab) * (ralpprblp / rab - 1.0);
  double ratio2Lp = std::exp(expoLp);
  double outlLp = 2.0 / (1.0 + ratio2Lp);

  double ranbprbnb[2], expoNb[2], ratio2Nb[2], outlNb[2];
  for (int i = 0; i < nbb; ++i) {
    ranbprbnb[i] = ranb[i] + rbnb[i] + 1.0e-12;
    expoNb[i] = (par.hbnbcut / radab) * (ranbprbnb[i] / rab - 1.0);
    ratio2Nb[i] = std::pow(std::exp(-expoNb[i]), 1.0);
    outlNb[i] = (2.0 / (1.0 + ratio2Nb[i])) - 1.0;
  }
  double outlNbTot = 1.0;
  for (int i = 0; i < nbb; ++i)
    outlNbTot = outlNbTot * outlNb[i];

  double ratio1 = std::pow(rab2 / par.hblongcut, par.hbalp);
  double dampl = 1.0 / (1.0 + ratio1);
  double shortcut = par.hbscut * radab;
  double ratio3 = std::pow(shortcut / rab2, par.hbalp);
  double damps = 1.0 / (1.0 + ratio3);
  double damp = damps * dampl;
  double ddamp = ((((-2.0 * par.hbalp) * ratio1) / (1.0 + ratio1))) +
                 ((((2.0 * par.hbalp) * ratio3) / (1.0 + ratio3)));
  double rbhdamp = damp * (((pBh / rbh2) / rbh));
  double rabdamp = damp * (((pAb / rab2) / rab));
  double rdamp = rbhdamp + rabdamp;

  double qh = chargeFactor(qa[h], par.hbst, par.hbsf, false);
  double qaa = chargeFactor(qa[a], par.hbst, par.hbsf, true);
  double qb = chargeFactor(qa[b], par.hbst, par.hbsf, true);
  double qhoutl = ((((qh * outl) * outlNbTot)) * outlLp);
  double const_ = ((((ca[1] * qaa) * cb[0]) * qb) * par.xhaciGlobAbh);
  energy = (((-rdamp) * qhoutl) * const_);

  double aterm = ((((((-rdamp) * qh) * outlNbTot)) * outlLp) * const_);
  double nbterm = ((((((-rdamp) * qh) * outl)) * outlLp) * const_);
  double lpterm = ((((((-rdamp) * qh) * outl)) * outlNbTot) * const_);
  double dterm = ((-qhoutl) * const_);

  double ga[3] = { 0.0, 0.0, 0.0 }, gb[3] = { 0.0, 0.0, 0.0 },
         gh[3] = { 0.0, 0.0, 0.0 };
  double dg[3], dga[3], dgb[3], dgh[3], gi;
  gi = ((((rabdamp + rbhdamp) * ddamp) - ((3.0) * rabdamp)) / rab2);
  gi = gi * dterm;
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drab[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = dg[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = -dg[c];
  gi = ((((-3.0) * rbhdamp)) / rbh2);
  gi = gi * dterm;
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drbh[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] + dg[c];
  for (int c = 0; c < 3; ++c)
    gh[c] = -dg[c];

  double o2 = 1.0 + ratio2;
  double tmp1 = (((((((-2.0 * aterm) * ratio2) * expo) / (o2 * o2))) /
                 (rahprbh - rab)));
  gi = (((-tmp1) * rahprbh) / rab2);
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drab[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dg[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] - dg[c];
  gi = (tmp1 / rah);
  for (int c = 0; c < 3; ++c)
    dga[c] = gi * drah[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dga[c];
  gi = (tmp1 / rbh);
  for (int c = 0; c < 3; ++c)
    dgb[c] = gi * drbh[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] + dgb[c];
  for (int c = 0; c < 3; ++c)
    dgh[c] = -dga[c] - dgb[c];
  for (int c = 0; c < 3; ++c)
    gh[c] = gh[c] + dgh[c];

  double o2lp = 1.0 + ratio2Lp;
  double tmp3 = (((((((-2.0 * lpterm) * ratio2Lp) * expoLp) / (o2lp * o2lp))) /
                 (ralpprblp - rab)));
  gi = (((-tmp3) * ralpprblp) / rab2);
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drab[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dg[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] - dg[c];
  gi = (tmp3 / ralp);
  for (int c = 0; c < 3; ++c)
    dga[c] = gi * dralp[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dga[c];
  gi = (tmp3 / (rblp + 1.0e-12));
  for (int c = 0; c < 3; ++c)
    dgb[c] = gi * drblp[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] - dga[c]; // verbatim (uses dga, not dgb)
  double glp[3];
  for (int c = 0; c < 3; ++c)
    glp[c] = -dga[c];

  double vsq = sumSq(vector);
  double vpw = std::pow(vsq, 1.5);
  double gii[3][3];
  for (int col = 0; col < 3; ++col) {
    for (int row = 0; row < 3; ++row) {
      double unit = (row == col) ? -1.0 : 0.0;
      gii[row][col] =
        (((-lpDist) * static_cast<double>(2))) *
        ((unit / vnorm) + ((vector[row] * vector[col]) / vpw));
    }
  }
  double gnbLp[3] = { 0.0, 0.0, 0.0 };
  for (int row = 0; row < 3; ++row) {
    for (int col = 0; col < 3; ++col)
      gnbLp[row] = gnbLp[row] + gii[row][col] * glp[col];
  }

  std::vector<double> tmp2(2, 0.0), giNb(2, 0.0);
  double gnb[2][3] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
  std::vector<char> mask(2, 1);
  for (int i = 0; i < nbb; ++i) {
    mask[i] = 0;
    double prod = 1.0;
    for (int k = 0; k < nbb; ++k) {
      if (mask[k])
        prod = prod * outlNb[k];
    }
    double o2nb = 1.0 + ratio2Nb[i];
    tmp2[i] = ((((((((2.0 * nbterm) * prod) * ratio2Nb[i]) * expoNb[i]) /
                  (o2nb * o2nb))) /
                (ranbprbnb[i] - rab)));
    giNb[i] = (((-tmp2[i]) * ranbprbnb[i]) / rab2);
    for (int c = 0; c < 3; ++c)
      dg[c] = giNb[i] * drab[c];
    for (int c = 0; c < 3; ++c)
      ga[c] = ga[c] + dg[c];
    for (int c = 0; c < 3; ++c)
      gb[c] = gb[c] - dg[c];
    mask[i] = 1;
  }
  for (int i = 0; i < nbb; ++i) {
    giNb[i] = (tmp2[i] / ranb[i]);
    for (int c = 0; c < 3; ++c)
      dga[c] = giNb[i] * dranb[i][c];
    for (int c = 0; c < 3; ++c)
      ga[c] = ga[c] + dga[c];
    giNb[i] = (tmp2[i] / rbnb[i]);
    for (int c = 0; c < 3; ++c)
      dgb[c] = giNb[i] * drbnb[i][c];
    for (int c = 0; c < 3; ++c)
      gb[c] = gb[c] + dgb[c];
    for (int c = 0; c < 3; ++c)
      gnb[i][c] = -dga[c] - dgb[c];
  }

  for (int c = 0; c < 3; ++c) {
    gdr[3 * a + c] = gdr[3 * a + c] + ga[c];
    gdr[3 * b + c] = gdr[3 * b + c] + gb[c] + gnbLp[c];
    gdr[3 * h + c] = gdr[3 * h + c] + gh[c];
  }
  for (int i = 0; i < nbb; ++i) {
    int inb = nbrsB[i];
    for (int c = 0; c < 3; ++c)
      gdr[3 * inb + c] =
        gdr[3 * inb + c] + (gnb[i][c] - (gnbLp[c] / 2.0));
  }
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      sigma[r][c] = sigma[r][c] + ((mcf * ga[c]) * xyz[3 * a + r]);
      sigma[r][c] = sigma[r][c] + ((mcf * gb[c]) * xyz[3 * b + r]);
      sigma[r][c] = sigma[r][c] + ((mcf * gh[c]) * xyz[3 * h + r]);
      sigma[r][c] = sigma[r][c] + ((mcf * gnbLp[c]) * xyz[3 * b + r]);
    }
  }
  for (int c = 0; c < 3; ++c)
    gnbLp[c] = gnbLp[c] / 2.0;
  for (int i = 0; i < nbb; ++i) {
    int inb = nbrsB[i];
    for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c)
        sigma[r][c] = sigma[r][c] +
                      ((mcf * (gnb[i][c] - gnbLp[c])) * xyz[3 * inb + r]);
    }
  }
  return true;
}

bool abhEg3(int a, int b, int h, int c,
            const std::vector<int>& cNbrs, const std::vector<int>& numbers,
            const std::vector<double>& xyz, const std::vector<double>& qa,
            const std::vector<double>& hbBas,
            const std::vector<double>& hbAci,
            const std::vector<double>& sqrab, const std::vector<double>& srab,
            const std::vector<double>& rad, const HbParams& par, double mcf,
            double& energy, std::vector<double>& gdr, double sigma[3][3],
            Environment& env)
{
  int n = static_cast<int>(numbers.size());
  (void)sqrab;
  (void)srab; // unused in the 0d body (kept for signature parity)
  if (n <= 0 || a < 0 || b < 0 || h < 0 || c < 0 || a >= n || b >= n ||
      h >= n || c >= n) {
    env.error("empty HB setup", "abhEg3");
    return false;
  }
  double pBh = 1.0 + par.hbabmix;
  double pAb = -par.hbabmix;
  energy = 0.0;
  double ca[2], cb[2];
  hbStrengths(a, b, hbBas, hbAci, ca, cb);

  double dranb[3], drbnb[3];
  sub3(xyz, a, c, dranb);
  sub3(xyz, b, c, drbnb);
  double ranb2 = sumSq(dranb);
  double ranb = std::sqrt(ranb2);
  double rbnb2 = sumSq(drbnb);
  double rbnb = std::sqrt(rbnb2);

  double drab[3], drah[3], drbh[3];
  sub3(xyz, a, b, drab);
  sub3(xyz, a, h, drah);
  sub3(xyz, b, h, drbh);
  double rab = norm3(drab);
  double rab2 = rab * rab;
  double rah = norm3(drah);
  double rbh = norm3(drbh);
  double rbh2 = rbh * rbh;
  double rahprbh = rah + rbh + 1.0e-12;
  double radab = rad[numbers[a] - 1] + rad[numbers[b] - 1];

  double expo = (par.hbacut / radab) * (rahprbh / rab - 1.0);
  if (expo > 15.0)
    return true;
  double ratio2 = std::exp(expo);
  double outl = 2.0 / (1.0 + ratio2);

  double ranbprbnb = ranb + rbnb + 1.0e-12;
  double expoNb = (par.hbnbcut / radab) * (ranbprbnb / rab - 1.0);
  double ratio2Nb = std::exp(-expoNb);
  double outlNbTot = (2.0 / (1.0 + ratio2Nb)) - 1.0;

  double ratio1 = std::pow(rab2 / par.hblongcut, par.hbalp);
  double dampl = 1.0 / (1.0 + ratio1);
  double shortcut = par.hbscut * radab;
  double rb6 = shortcut / rab2;
  double rb62 = rb6 * rb6;
  double ratio3 = rb62 * (rb62 * rb62); // **6 as x2*x4
  double damps = 1.0 / (1.0 + ratio3);
  double damp = damps * dampl;
  double ddamp = ((((-2.0 * par.hbalp) * ratio1) / (1.0 + ratio1))) +
                 ((((2.0 * par.hbalp) * ratio3) / (1.0 + ratio3)));
  double rbhdamp = damp * (((pBh / rbh2) / rbh));
  double rabdamp = damp * (((pAb / rab2) / rab));
  double rdamp = rbhdamp + rabdamp;

  // Torsion list over C's neighbours except B (0d: no image cycles).
  struct TorsEntry
  {
    int r;
  };
  std::vector<TorsEntry> tlist;
  for (int nb : cNbrs) {
    if (nb == b)
      continue;
    tlist.push_back({ nb });
  }
  int ntors = static_cast<int>(tlist.size());
  std::vector<double> etmp(ntors, 0.0);
  std::vector<double> g4tmp(12 * (ntors > 0 ? ntors : 1), 0.0);
  for (int i = 0; i < ntors; ++i) {
    int ii = tlist[i].r;
    double phi = dihedralValue0d(ii, b, c, h, xyz);
    double e = 0.0;
    double g4[3][4];
    torsNciMul(ii, b, c, h, 2, phi, pi / 2.0, par.torsHb, xyz, e, g4);
    etmp[i] = e;
    for (int col = 0; col < 4; ++col) {
      for (int k = 0; k < 3; ++k)
        g4tmp[(i * 4 + col) * 3 + k] = g4[k][col];
    }
  }
  double etors = 1.0;
  for (int i = 0; i < ntors; ++i)
    etors = etors * etmp[i];
  std::vector<double> gtors(3 * n, 0.0);
  std::vector<char> tmask(ntors > 0 ? ntors : 1, 1);
  for (int i = 0; i < ntors; ++i) {
    tmask[i] = 0;
    double prod = 1.0;
    for (int k = 0; k < ntors; ++k) {
      if (tmask[k])
        prod = prod * etmp[k];
    }
    int ii = tlist[i].r;
    for (int k = 0; k < 3; ++k) {
      gtors[3 * ii + k] = gtors[3 * ii + k] + g4tmp[(i * 4 + 0) * 3 + k] * prod;
      gtors[3 * b + k] = gtors[3 * b + k] + g4tmp[(i * 4 + 1) * 3 + k] * prod;
      gtors[3 * c + k] = gtors[3 * c + k] + g4tmp[(i * 4 + 2) * 3 + k] * prod;
      gtors[3 * h + k] = gtors[3 * h + k] + g4tmp[(i * 4 + 3) * 3 + k] * prod;
    }
    tmask[i] = 1;
  }

  double r0 = 120.0;
  double phi0 = ((r0 * pi) / 180.0);
  double bshift = par.bendHb;
  double fc = 1.0 - bshift;
  double eangl = 0.0;
  double g3[3][3];
  // jj/kk/ll are B/C/H for every torsion entry by construction.
  bendNciMul(b, c, h, phi0, fc, xyz, eangl, g3);
  std::vector<double> gangl(3 * n, 0.0);
  for (int k = 0; k < 3; ++k) {
    gangl[3 * b + k] = gangl[3 * b + k] + g3[k][0];
    gangl[3 * c + k] = gangl[3 * c + k] + g3[k][1];
    gangl[3 * h + k] = gangl[3 * h + k] + g3[k][2];
  }

  double qh = chargeFactor(qa[h], par.hbst, par.hbsf, false);
  double qaa = chargeFactor(qa[a], par.hbst, par.hbsf, true);
  double qb = chargeFactor(qa[b], par.hbst, par.hbsf, true);
  double qhoutl = ((qh * outl) * outlNbTot);
  double const_ = ((((ca[1] * qaa) * cb[0]) * qb) * par.xhaciCoh);
  energy = ((((((-rdamp) * qhoutl)) * eangl) * etors) * const_);

  double aterm = ((((((( -rdamp) * qh) * outlNbTot)) * eangl) * etors) * const_);
  double nbterm =
    ((((((( -rdamp) * qh) * outl)) * eangl) * etors) * const_);
  double dterm = (((((-qhoutl)) * eangl) * etors) * const_);
  double tterm = (((((-rdamp) * qhoutl)) * eangl) * const_);
  double bterm = (((((-rdamp) * qhoutl)) * etors) * const_);

  double ga[3] = { 0.0, 0.0, 0.0 }, gb[3] = { 0.0, 0.0, 0.0 },
         gh[3] = { 0.0, 0.0, 0.0 };
  double dg[3], dga[3], dgb[3], dgh[3], gi;
  gi = ((((rabdamp + rbhdamp) * ddamp) - ((3.0) * rabdamp)) / rab2);
  gi = gi * dterm;
  for (int k = 0; k < 3; ++k)
    dg[k] = gi * drab[k];
  for (int k = 0; k < 3; ++k)
    ga[k] = dg[k];
  for (int k = 0; k < 3; ++k)
    gb[k] = -dg[k];
  gi = ((((-3.0) * rbhdamp)) / rbh2);
  gi = gi * dterm;
  for (int k = 0; k < 3; ++k)
    dg[k] = gi * drbh[k];
  for (int k = 0; k < 3; ++k)
    gb[k] = gb[k] + dg[k];
  for (int k = 0; k < 3; ++k)
    gh[k] = -dg[k];

  double o2 = 1.0 + ratio2;
  double tmp1 = (((((((-2.0 * aterm) * ratio2) * expo) / (o2 * o2))) /
                 (rahprbh - rab)));
  gi = (((-tmp1) * rahprbh) / rab2);
  for (int k = 0; k < 3; ++k)
    dg[k] = gi * drab[k];
  for (int k = 0; k < 3; ++k)
    ga[k] = ga[k] + dg[k];
  for (int k = 0; k < 3; ++k)
    gb[k] = gb[k] - dg[k];
  gi = (tmp1 / rah);
  for (int k = 0; k < 3; ++k)
    dga[k] = gi * drah[k];
  for (int k = 0; k < 3; ++k)
    ga[k] = ga[k] + dga[k];
  gi = (tmp1 / rbh);
  for (int k = 0; k < 3; ++k)
    dgb[k] = gi * drbh[k];
  for (int k = 0; k < 3; ++k)
    gb[k] = gb[k] + dgb[k];
  for (int k = 0; k < 3; ++k)
    dgh[k] = -dga[k] - dgb[k];
  for (int k = 0; k < 3; ++k)
    gh[k] = gh[k] + dgh[k];

  double o2nb = 1.0 + ratio2Nb;
  double tmp2 = (((((((2.0 * nbterm) * ratio2Nb) * expoNb) / (o2nb * o2nb))) /
                 (ranbprbnb - rab)));
  double giNb = (((-tmp2) * ranbprbnb) / rab2);
  for (int k = 0; k < 3; ++k)
    dg[k] = giNb * drab[k];
  for (int k = 0; k < 3; ++k)
    ga[k] = ga[k] + dg[k];
  for (int k = 0; k < 3; ++k)
    gb[k] = gb[k] - dg[k];
  giNb = (tmp2 / ranb);
  double gnb[3];
  for (int k = 0; k < 3; ++k)
    dga[k] = giNb * dranb[k];
  for (int k = 0; k < 3; ++k)
    ga[k] = ga[k] + dga[k];
  giNb = (tmp2 / rbnb);
  for (int k = 0; k < 3; ++k)
    dgb[k] = giNb * drbnb[k];
  for (int k = 0; k < 3; ++k)
    gb[k] = gb[k] + dgb[k];
  for (int k = 0; k < 3; ++k)
    gnb[k] = -dga[k] - dgb[k];

  for (int i = 0; i < ntors; ++i) {
    int ii = tlist[i].r;
    for (int k = 0; k < 3; ++k)
      gdr[3 * ii + k] = gdr[3 * ii + k] + gtors[3 * ii + k] * tterm;
  }
  for (int k = 0; k < 3; ++k) {
    gdr[3 * b + k] = gdr[3 * b + k] + gtors[3 * b + k] * tterm;
    gdr[3 * c + k] = gdr[3 * c + k] + gtors[3 * c + k] * tterm;
    gdr[3 * h + k] = gdr[3 * h + k] + gtors[3 * h + k] * tterm;
  }
  for (int k = 0; k < 3; ++k) {
    gdr[3 * b + k] = gdr[3 * b + k] + gangl[3 * b + k] * bterm;
    gdr[3 * c + k] = gdr[3 * c + k] + gangl[3 * c + k] * bterm;
    gdr[3 * h + k] = gdr[3 * h + k] + gangl[3 * h + k] * bterm;
  }
  for (int r = 0; r < 3; ++r) {
    for (int k = 0; k < 3; ++k) {
      sigma[r][k] = sigma[r][k] + ((mcf * ga[k]) * xyz[3 * a + r]);
      sigma[r][k] = sigma[r][k] + ((mcf * gb[k]) * xyz[3 * b + r]);
      sigma[r][k] = sigma[r][k] + ((mcf * gh[k]) * xyz[3 * h + r]);
      sigma[r][k] = sigma[r][k] + ((mcf * gnb[k]) * xyz[3 * c + r]);
    }
  }
  for (int i = 0; i < ntors; ++i) {
    int ii = tlist[i].r;
    for (int r = 0; r < 3; ++r) {
      for (int k = 0; k < 3; ++k)
        sigma[r][k] = sigma[r][k] +
                      ((((mcf * tterm) * gtors[3 * ii + k])) * xyz[3 * ii + r]);
    }
  }
  for (int r = 0; r < 3; ++r) {
    for (int k = 0; k < 3; ++k) {
      sigma[r][k] = sigma[r][k] +
                    ((((mcf * tterm) * gtors[3 * b + k])) * xyz[3 * b + r]);
      sigma[r][k] = sigma[r][k] +
                    ((((mcf * tterm) * gtors[3 * c + k])) * xyz[3 * c + r]);
      sigma[r][k] = sigma[r][k] +
                    ((((mcf * tterm) * gtors[3 * h + k])) * xyz[3 * h + r]);
      sigma[r][k] = sigma[r][k] +
                    ((((mcf * bterm) * gangl[3 * b + k])) * xyz[3 * b + r]);
      sigma[r][k] = sigma[r][k] +
                    ((((mcf * bterm) * gangl[3 * c + k])) * xyz[3 * c + r]);
      sigma[r][k] = sigma[r][k] +
                    ((((mcf * bterm) * gangl[3 * h + k])) * xyz[3 * h + r]);
    }
  }

  for (int k = 0; k < 3; ++k) {
    gdr[3 * a + k] = gdr[3 * a + k] + ga[k];
    gdr[3 * b + k] = gdr[3 * b + k] + gb[k];
    gdr[3 * h + k] = gdr[3 * h + k] + gh[k];
    gdr[3 * c + k] = gdr[3 * c + k] + gnb[k];
  }
  return true;
}

bool rbxEg(int a, int b, int x, double cb, double cx,
           const std::vector<int>& numbers, const std::vector<double>& xyz,
           const std::vector<double>& qa, const std::vector<double>& rad,
           const HbParams& par, double& energy, double gdr[3][3],
           double sigma[3][3], Environment& env)
{
  int n = static_cast<int>(numbers.size());
  if (n <= 0 || a < 0 || b < 0 || x < 0 || a >= n || b >= n || x >= n) {
    env.error("empty XB setup", "rbxEg");
    return false;
  }
  for (int c = 0; c < 3; ++c) {
    gdr[c][0] = 0.0;
    gdr[c][1] = 0.0;
    gdr[c][2] = 0.0;
  }
  energy = 0.0;

  double drax[3], drbx[3], drab[3];
  sub3(xyz, a, x, drax);
  sub3(xyz, b, x, drbx);
  sub3(xyz, a, b, drab);
  double rab2 = sumSq(drab);
  double rab = std::sqrt(rab2);
  double rax2 = sumSq(drax);
  double rax = std::sqrt(rax2) + 1.0e-12;
  double rbx2 = sumSq(drbx);
  double rbx = std::sqrt(rbx2) + 1.0e-12;

  double expo = par.xbacut * ((rax + rbx) / rab - 1.0);
  if (expo > 15.0)
    return true;
  double ratio2 = std::exp(expo);
  double outl = 2.0 / (1.0 + ratio2);

  double ratio1 = std::pow(rbx2 / par.hblongcutXb, par.hbalp);
  double dampl = 1.0 / (1.0 + ratio1);
  double shortcut = par.xbscut * (rad[numbers[a] - 1] + rad[numbers[b] - 1]);
  double ratio3 = std::pow(shortcut / rbx2, par.hbalp);
  double damps = 1.0 / (1.0 + ratio3);
  double damp = damps * dampl;
  double rdamp = (damp / rbx2) / rbx;

  double qx = chargeFactor(qa[x], par.xbst, par.xbsf, false);
  double qb = chargeFactor(qa[b], par.xbst, par.xbsf, true);
  double const_ = (((cb * qb) * cx) * qx);
  double aterm = ((-rdamp) * const_);
  double dterm = ((-outl) * const_);
  energy = ((((-rdamp) * outl) * const_));

  double ga[3] = { 0.0, 0.0, 0.0 }, gb[3] = { 0.0, 0.0, 0.0 },
         gx[3] = { 0.0, 0.0, 0.0 };
  double dg[3], dga[3], dgb[3], dgx[3], gi;
  double dA = (((2.0 * par.hbalp) * ratio1) / (1.0 + ratio1));
  double dB = (((2.0 * par.hbalp) * ratio3) / (1.0 + ratio3));
  double o2 = 1.0 + ratio2;
  gi = (((rdamp * (((-dA) + dB) - 3.0))) / rbx2);
  gi = gi * dterm;
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drbx[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = dg[c];
  for (int c = 0; c < 3; ++c)
    gx[c] = -dg[c];

  double raxprbx = rax + rbx;
  gi = ((((((((2.0 * ratio2) * expo) * raxprbx) / (o2 * o2)) /
           (raxprbx - rab))) /
         rab2));
  gi = gi * aterm;
  for (int c = 0; c < 3; ++c)
    dg[c] = gi * drab[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dg[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] - dg[c];

  gi = (((((((-2.0 * ratio2) * expo) / (o2 * o2)) / (raxprbx - rab))) / rax));
  gi = gi * aterm;
  for (int c = 0; c < 3; ++c)
    dga[c] = gi * drax[c];
  for (int c = 0; c < 3; ++c)
    ga[c] = ga[c] + dga[c];
  gi = (((((((-2.0 * ratio2) * expo) / (o2 * o2)) / (raxprbx - rab))) / rbx));
  gi = gi * aterm;
  for (int c = 0; c < 3; ++c)
    dgb[c] = gi * drbx[c];
  for (int c = 0; c < 3; ++c)
    gb[c] = gb[c] + dgb[c];
  for (int c = 0; c < 3; ++c)
    dgx[c] = -dga[c] - dgb[c];
  for (int c = 0; c < 3; ++c)
    gx[c] = gx[c] + dgx[c];

  for (int c = 0; c < 3; ++c) {
    gdr[c][0] = ga[c];
    gdr[c][1] = gb[c];
    gdr[c][2] = gx[c];
  }
  // sigma is NOT mcf-scaled here (verbatim).
  for (int r = 0; r < 3; ++r) {
    for (int c = 0; c < 3; ++c) {
      sigma[r][c] = sigma[r][c] + (ga[c] * xyz[3 * a + r]);
      sigma[r][c] = sigma[r][c] + (gb[c] * xyz[3 * b + r]);
      sigma[r][c] = sigma[r][c] + (gx[c] * xyz[3 * x + r]);
    }
  }
  return true;
}

bool isXAtom(int atomicNumber)
{
  return atomicNumber == 17 || atomicNumber == 35 || atomicNumber == 53 ||
         atomicNumber == 16 || atomicNumber == 34 || atomicNumber == 52 ||
         atomicNumber == 15 || atomicNumber == 33 || atomicNumber == 51;
}

bool hbDonorLists(int n, const std::vector<int>& numbers,
                  const std::vector<int>& hyb, const std::vector<int>& piFlags,
                  const std::vector<double>& qa,
                  const std::vector<std::vector<int>>& neighbours,
                  const std::vector<int>& firstNb, const int* group,
                  const std::vector<double>& xhbas, double hQaThr,
                  double qaBThr, const std::vector<int>& bpair,
                  const std::vector<double>& hbBas,
                  const std::vector<double>& hbAci, HbDonorLists& lists,
                  Environment& env)
{
  if (n <= 0 || static_cast<int>(numbers.size()) != n ||
      static_cast<int>(qa.size()) != n) {
    env.error("empty HB perception setup", "hbDonorLists");
    return false;
  }
  lists.hatH.clear();
  lists.hatAB.clear();
  lists.xbatA.clear();
  lists.xbatB.clear();
  lists.xbatX.clear();
  // H atoms with positive topology charges.
  for (int i = 0; i < n; ++i) {
    if (numbers[i] != 1)
      continue;
    if (hyb[i] == 1)
      continue;
    double ff = hQaThr;
    int j = firstNb[i];
    if (j < 0)
      continue;
    if (numbers[j] > 10)
      ff -= static_cast<double>(0.20f);
    if (numbers[j] == 6 && hyb[j] == 3)
      ff += static_cast<double>(0.05f);
    if (qa[i] > ff)
      lists.hatH.push_back(i);
  }
  // AB pairs with negative charges and HB strength.
  for (int i = 0; i < n; ++i) {
    if (numbers[i] == 6 && piFlags[i] == 0)
      continue;
    double ff = qaBThr;
    if (numbers[i] > 10)
      ff += static_cast<double>(0.2f);
    if (qa[i] > ff)
      continue;
    for (int j = 0; j < i; ++j) {
      ff = qaBThr;
      if (numbers[j] > 10)
        ff += static_cast<double>(0.2f);
      if (qa[j] > ff)
        continue;
      double ci[2], cj[2];
      hbStrengths(i, j, hbBas, hbAci, ci, cj);
      if (ci[0] * cj[1] < 1.0e-6 && ci[1] * cj[0] < 1.0e-6)
        continue;
      if (numbers[j] == 6 && piFlags[j] == 0)
        continue;
      lists.hatAB.emplace_back(i, j);
    }
  }
  // A-X...B halogen triples (0d: all shifts central, iTrDum = 1).
  for (int i = 0; i < n; ++i) {
    for (int x : neighbours[i]) {
      if (!isXAtom(numbers[x]))
        continue;
      if (numbers[x] == 16 &&
          static_cast<int>(neighbours[x].size()) > 2)
        continue;
      for (int j = 0; j < n; ++j) {
        if (i == j || j == x)
          continue;
        if (bpair[j * n + x] <= 3)
          continue;
        if (xhbas[numbers[j] - 1] < 1.0e-6)
          continue;
        if (group[numbers[j] - 1] == 4) {
          if (piFlags[j] == 0 || qa[j] > static_cast<double>(0.05f))
            continue;
        }
        lists.xbatA.push_back(i);
        lists.xbatB.push_back(j);
        lists.xbatX.push_back(x);
      }
    }
  }
  return true;
}

bool hbBondTriplets(int n, const std::vector<std::pair<int, int>>& hatAB,
                     const std::vector<int>& hatH,
                     const std::vector<double>& xyz, double hbThr1,
                     const std::vector<int>& bpair,
                     std::vector<HbBondTriplet>& triplets, Environment& env)
{
  if (n <= 0 || static_cast<int>(xyz.size()) != 3 * n) {
    env.error("empty HB triplet setup", "hbBondTriplets");
    return false;
  }
  triplets.clear();
  for (const auto& ab : hatAB) {
    int i = ab.first, j = ab.second;
    double dx = xyz[3 * i] - xyz[3 * j];
    double dy = xyz[3 * i + 1] - xyz[3 * j + 1];
    double dz = xyz[3 * i + 2] - xyz[3 * j + 2];
    double rab = dx * dx + dy * dy + dz * dz;
    if (rab > hbThr1)
      continue;
    // 0d: iTrDum = 1 sits in the central cell, so pairs are non-bonded
    // unless directly bonded.
    bool ijnonbond = bpair[i * n + j] != 1;
    for (int nh : hatH) {
      if (bpair[nh * n + i] == 1 && ijnonbond) {
        HbBondTriplet t;
        t.a = i;
        t.b = j;
        t.h = nh;
        triplets.push_back(t);
      }
      if (bpair[nh * n + j] == 1 && ijnonbond) {
        HbBondTriplet t;
        t.a = j;
        t.b = i;
        t.h = nh;
        triplets.push_back(t);
      }
    }
  }
  return true;
}

bool hbAhbMaps(int n, const std::vector<int>& numbers,
               const std::vector<std::pair<int, int>>& bonds,
               const std::vector<HbBondTriplet>& triplets, HbBondMaps& maps,
               Environment& env)
{
  if (n <= 0 || static_cast<int>(numbers.size()) != n) {
    env.error("empty HB map setup", "hbAhbMaps");
    return false;
  }
  int nbond = static_cast<int>(bonds.size());
  int nbTri = static_cast<int>(triplets.size());
  maps.ahA.clear();
  maps.ahH.clear();
  maps.bAtoms.clear();
  maps.isAbh.assign(n, 0);
  maps.nrHb.assign(nbond, 0);
  maps.mapAbh.assign(n, 0);
  maps.mapNab = 0;
  maps.mapNh = 0;
  maps.bMax = 1;
  // Sizing pass mirroring bond_hb_AHB_set1: lin table (1-based, slot 0
  // zeroed), AH count, B max, participation flags and per-bond counts.
  // B_count starts at 1 and nr_hb takes it unconditionally, so bonds
  // without a match keep the stale running count (verbatim).
  std::vector<int> linA, linH, linTrA, linTrH;
  linA.assign(nbTri + 1, 0);
  linH.assign(nbTri + 1, 0);
  linTrA.assign(nbTri + 1, 0);
  linTrH.assign(nbTri + 1, 0);
  int totCount = 0, ahCount = 0, bCount = 1, bMax = 1, linDiff = 0;
  for (int i = 0; i < nbond; ++i) {
    int jj = bonds[i].first, ii = bonds[i].second;
    int hbH = -1, hbA = -1;
    if (numbers[ii] == 1) {
      hbH = ii;
      hbA = jj;
    } else if (numbers[jj] == 1) {
      hbH = jj;
      hbA = ii;
    } else {
      continue;
    }
    if (!(numbers[hbA] == 7 || numbers[hbA] == 8))
      continue;
    for (int j = 0; j < nbTri; ++j) {
      int bat = triplets[j].b;
      int atB = numbers[bat];
      int aat = triplets[j].a;
      int hat = triplets[j].h;
      if (hbA == aat && hbH == hat) {
        if (atB == 7 || atB == 8) {
          ++totCount;
          linA[totCount] = hbA;
          linH[totCount] = hbH;
          linTrA[totCount] = 1;
          linTrH[totCount] = 1;
          maps.isAbh[bat] = 1;
          if (linA[totCount] - linA[totCount - 1] == 0 &&
              linH[totCount] - linH[totCount - 1] == 0)
            linDiff = 0;
          else
            linDiff = 1;
          bool itrSame = true; // 0d: all shifts central
          if (linDiff == 0 && itrSame)
            ++bCount;
          if (linDiff != 0 || !itrSame) {
            ++ahCount;
            bCount = 1;
          }
          if (bCount > bMax)
            bMax = bCount;
        }
      }
    }
    maps.isAbh[hbA] = 1;
    maps.isAbh[hbH] = 1;
    maps.nrHb[i] = bCount;
  }
  // Fill pass mirroring bond_hb_AHB_set.
  std::vector<int> ahA, ahH;
  std::vector<std::vector<int>> bAtoms;
  std::vector<int> bCounts;
  totCount = 0;
  ahCount = 0;
  bCount = 0;
  linDiff = 0;
  linA.assign(nbTri + 1, 0);
  linH.assign(nbTri + 1, 0);
  linTrA.assign(nbTri + 1, 0);
  linTrH.assign(nbTri + 1, 0);
  for (int i = 0; i < nbond; ++i) {
    int jj = bonds[i].first, ii = bonds[i].second;
    int hbH = -1, hbA = -1;
    if (numbers[ii] == 1) {
      hbH = ii;
      hbA = jj;
    } else if (numbers[jj] == 1) {
      hbH = jj;
      hbA = ii;
    } else {
      continue;
    }
    if (!(numbers[hbA] == 7 || numbers[hbA] == 8))
      continue;
    for (int j = 0; j < nbTri; ++j) {
      int bat = triplets[j].b;
      int atB = numbers[bat];
      int aat = triplets[j].a;
      int hat = triplets[j].h;
      if (hbA == aat && hbH == hat) {
        if (atB == 7 || atB == 8) {
          ++totCount;
          linA[totCount] = hbA;
          linH[totCount] = hbH;
          linTrA[totCount] = 1;
          linTrH[totCount] = 1;
          if (linA[totCount] - linA[totCount - 1] == 0 &&
              linH[totCount] - linH[totCount - 1] == 0)
            linDiff = 0;
          else
            linDiff = 1;
          bool itrSame = true; // 0d: all shifts central
          if (linDiff == 0 && itrSame)
            ++bCount;
          if (linDiff != 0 || !itrSame) {
            ++ahCount;
            ahA.push_back(hbA);
            ahH.push_back(hbH);
            bAtoms.emplace_back();
            bCounts.push_back(0);
            bCount = 1;
          }
          bCounts[ahCount - 1] = bCount;
          bAtoms[ahCount - 1].push_back(bat);
        }
      }
    }
  }
  maps.ahA = ahA;
  maps.ahH = ahH;
  maps.bAtoms = bAtoms;
  maps.bMax = bMax;
  // Compact atom->HB-index maps (0-based within the AB and H classes).
  int jc = 0, kc = 0;
  for (int i = 0; i < n; ++i) {
    if (!maps.isAbh[i])
      continue;
    if (numbers[i] == 1)
      maps.mapAbh[i] = jc++;
    else
      maps.mapAbh[i] = kc++;
  }
  maps.mapNh = jc;
  maps.mapNab = kc;
  return true;
}

bool hbBasicity(int n, const std::vector<int>& numbers,
                const std::vector<int>& itag,
                const std::vector<std::vector<int>>& neighbours,
                const std::vector<int>& counts,
                const std::vector<double>& xhbas, std::vector<double>& hbBas,
                Environment& env)
{
  if (n <= 0 || static_cast<int>(numbers.size()) != n) {
    env.error("empty HB basicity setup", "hbBasicity");
    return false;
  }
  hbBas.assign(n, 1.0);
  for (int i = 0; i < n; ++i) {
    int nn = counts[i];
    int ati = numbers[i];
    hbBas[i] = xhbas[ati - 1];
    if (ati == 6 && nn == 2 && itag[i] == 1)
      hbBas[i] = static_cast<double>(1.46f);
    if (!(ati == 8 && nn == 1))
      continue;
    // 0d: single cell, so the neighbour-cell lookup always hits it and
    // nb(nn, i) is the last neighbour in nb order.
    int last = neighbours[i][nn - 1];
    if (numbers[last] == 6)
      hbBas[i] = static_cast<double>(0.68f);
    if (numbers[last] == 7)
      hbBas[i] = static_cast<double>(0.47f);
  }
  return true;
}

bool hbAcidity(int n, const std::vector<int>& numbers,
               const std::vector<std::vector<int>>& neighbours,
               const std::vector<int>& hyb, const std::vector<int>& piFlags,
               const std::vector<double>& xhaci, std::vector<double>& hbAci,
               Environment& env)
{
  if (n <= 0 || static_cast<int>(numbers.size()) != n) {
    env.error("empty HB acidity setup", "hbAcidity");
    return false;
  }
  hbAci.assign(n, 1.0);
  for (int i = 0; i < n; ++i)
    hbAci[i] = xhaci[numbers[i] - 1];
  for (int i = 0; i < n; ++i) {
    if (neighbours[i].empty())
      continue;
    int nn = neighbours[i][0];
    hbAci[i] = xhaci[numbers[i] - 1];
    if (isAmideHydrogen(n, numbers, hyb, neighbours, piFlags, i))
      hbAci[nn] = hbAci[nn] * static_cast<double>(0.80f);
  }
  return true;
}

bool hbErfCoordination(int n, const std::vector<int>& numbers,
                       const std::vector<double>& xyz,
                       const std::vector<double>& rcov, const HbBondMaps& maps,
                       std::vector<double>& hbCn, std::vector<double>& hbDcn,
                       Environment& env)
{
  if (n <= 0 || static_cast<int>(numbers.size()) != n ||
      static_cast<int>(xyz.size()) != 3 * n) {
    env.error("empty HB erf CN setup", "hbErfCoordination");
    return false;
  }
  hbCn.assign(n, 0.0);
  hbDcn.assign(3 * n * n, 0.0);
  const double kn = 27.5;
  const double hlfosqrtpi = 1.0 / 1.77245385091;
  const double rcovScal = static_cast<double>(1.78f);
  int nr = static_cast<int>(maps.ahA.size());
  for (int i = 0; i < nr; ++i) {
    int iat = maps.ahH[i];
    int ati = numbers[iat];
    for (int b : maps.bAtoms[i]) {
      int jat = b;
      int atj = numbers[jat];
      double rx = xyz[3 * jat] - xyz[3 * iat];
      double ry = xyz[3 * jat + 1] - xyz[3 * iat + 1];
      double rz = xyz[3 * jat + 2] - xyz[3 * iat + 2];
      double r2 = rx * rx + ry * ry + rz * rz;
      if (r2 > 900.0)
        continue;
      double r = std::sqrt(r2);
      double rcovij = rcovScal * (rcov[ati - 1] + rcov[atj - 1]);
      double dr = r - rcovij;
      double tmp = 0.5 * (1.0 + std::erf((-kn * dr) / rcovij));
      double k2 = kn * kn;
      double d2 = dr * dr;
      double rc2 = rcovij * rcovij;
      double dtmp =
        -hlfosqrtpi * kn * std::exp((-k2 * d2) / rc2) / rcovij;
      hbCn[iat] += tmp;
      hbCn[jat] += tmp;
      double dx = dtmp * rx / r;
      double dy = dtmp * ry / r;
      double dz = dtmp * rz / r;
      hbDcn[(jat * n + jat) * 3] += dx;
      hbDcn[(jat * n + jat) * 3 + 1] += dy;
      hbDcn[(jat * n + jat) * 3 + 2] += dz;
      hbDcn[(iat * n + jat) * 3] += dx;
      hbDcn[(iat * n + jat) * 3 + 1] += dy;
      hbDcn[(iat * n + jat) * 3 + 2] += dz;
      hbDcn[(jat * n + iat) * 3] -= dx;
      hbDcn[(jat * n + iat) * 3 + 1] -= dy;
      hbDcn[(jat * n + iat) * 3 + 2] -= dz;
      hbDcn[(iat * n + iat) * 3] -= dx;
      hbDcn[(iat * n + iat) * 3 + 1] -= dy;
      hbDcn[(iat * n + iat) * 3 + 2] -= dz;
    }
  }
  return true;
}

bool hbTripletLists(int n, const std::vector<double>& xyz,
                    const std::vector<std::pair<int, int>>& hatAB,
                    const std::vector<int>& hatH,
                    const HbDonorLists& donors,
                    const std::vector<int>& bpair, double hbThr1,
                    double hbThr2, std::vector<HbTriple>& hb1,
                    std::vector<HbTriple>& hb2, std::vector<XbTriple>& xb,
                    Environment& env)
{
  if (n <= 0 || static_cast<int>(xyz.size()) != 3 * n) {
    env.error("empty HB list setup", "hbTripletLists");
    return false;
  }
  hb1.clear();
  hb2.clear();
  xb.clear();
  for (const auto& ab : hatAB) {
    int i = ab.first, j = ab.second;
    double dx = xyz[3 * i] - xyz[3 * j];
    double dy = xyz[3 * i + 1] - xyz[3 * j + 1];
    double dz = xyz[3 * i + 2] - xyz[3 * j + 2];
    double rab = dx * dx + dy * dy + dz * dz;
    if (rab > hbThr1)
      continue;
    // 0d: iTrDum = 1 sits in the central cell.
    bool ijnonbond = bpair[i * n + j] != 1;
    for (int nh : hatH) {
      double ax = xyz[3 * nh] - xyz[3 * i];
      double ay = xyz[3 * nh] - xyz[3 * i + 1];
      double az = xyz[3 * nh] - xyz[3 * i + 2];
      double rih = ax * ax + ay * ay + az * az;
      double bx = xyz[3 * nh] - xyz[3 * j];
      double by = xyz[3 * nh] - xyz[3 * j + 1];
      double bz = xyz[3 * nh] - xyz[3 * j + 2];
      double rjh = bx * bx + by * by + bz * bz;
      if (bpair[nh * n + i] == 1 && ijnonbond) {
        HbTriple t;
        t.a = i;
        t.b = j;
        t.h = nh;
        hb2.push_back(t);
        continue;
      }
      if (bpair[nh * n + j] == 1 && ijnonbond) {
        HbTriple t;
        t.a = j;
        t.b = i;
        t.h = nh;
        hb2.push_back(t);
        continue;
      }
      if (rab + rih + rjh < hbThr2) {
        HbTriple t;
        t.a = i;
        t.b = j;
        t.h = nh;
        hb1.push_back(t);
      }
    }
  }
  for (size_t ix = 0; ix < donors.xbatA.size(); ++ix) {
    int i = donors.xbatA[ix], j = donors.xbatB[ix];
    double dx = xyz[3 * j] - xyz[3 * i];
    double dy = xyz[3 * j + 1] - xyz[3 * i + 1];
    double dz = xyz[3 * j + 2] - xyz[3 * i + 2];
    double rab = dx * dx + dy * dy + dz * dz;
    if (rab > hbThr2)
      continue;
    XbTriple t;
    t.a = i;
    t.b = j;
    t.x = donors.xbatX[ix];
    xb.push_back(t);
  }
  return true;
}

bool egbondHb(int bondIdx, int iat, int jat, double rab, double rij,
              const std::vector<double>& drij, const double drijdcn[2],
              const std::vector<int>& numbers,
              const std::vector<double>& hbCn, const std::vector<double>& hbDcn,
              const std::vector<double>& xyz, double vbondScale,
              double steepness, double prefactor, const HbBondMaps& maps,
              std::vector<char>& consideredAbh, std::vector<double>& dEdcn,
              double& e, std::vector<double>& g, Environment& env)
{
  int n = static_cast<int>(hbCn.size());
  if (n <= 0 || static_cast<int>(g.size()) != 3 * n ||
      static_cast<int>(numbers.size()) != n) {
    env.error("empty HB bond setup", "egbondHb");
    return false;
  }
  int hbH = -1, hbA = -1;
  if (numbers[iat] == 1) {
    hbH = iat;
    hbA = jat;
  } else if (numbers[jat] == 1) {
    hbH = jat;
    hbA = iat;
  } else {
    return true; // no H endpoint (reference prints and returns)
  }
  double t1 = 1.0 - vbondScale;
  double t8 = (-t1 * hbCn[hbH] + 1.0) * steepness;
  double dr = rab - rij;
  double dum = prefactor * std::exp(-t8 * dr * dr);
  e += dum;
  double yy = 2.0 * t8 * dr * dum;
  double dx = -xyz[3 * jat] + xyz[3 * iat];
  double dy = -xyz[3 * jat + 1] + xyz[3 * iat + 1];
  double dz = -xyz[3 * jat + 2] + xyz[3 * iat + 2];
  double t4 = -yy * dx / rab;
  double t5 = -yy * dy / rab;
  double t6 = -yy * dz / rab;
  g[3 * iat] += t4;
  g[3 * iat + 1] += t5;
  g[3 * iat + 2] += t6;
  dEdcn[iat] += yy * drijdcn[0];
  t4 = yy * (dx / rab);
  t5 = yy * (dy / rab);
  t6 = yy * (dz / rab);
  g[3 * jat] += t4;
  g[3 * jat + 1] += t5;
  g[3 * jat + 2] += t6;
  dEdcn[jat] += yy * drijdcn[1];
  for (int k = 0; k < n; ++k) {
    g[3 * k] += drij[3 * k] * yy;
    g[3 * k + 1] += drij[3 * k + 1] * yy;
    g[3 * k + 2] += drij[3 * k + 2] * yy;
  }
  double zz = dum * steepness * dr * dr * t1;
  int nr = static_cast<int>(maps.ahA.size());
  for (int j = 0; j < nr; ++j) {
    if (!(maps.ahH[j] == hbH && maps.ahA[j] == hbA))
      continue;
    // 0d: all recorded shifts are central, so the shift test passes.
    for (int c = 0; c < 3; ++c)
      g[3 * hbH + c] += hbDcn[(hbH * n + hbH) * 3 + c] * zz;
    for (int hbB : maps.bAtoms[j]) {
      int mapA = maps.mapAbh[hbA];
      int mapB = maps.mapAbh[hbB];
      int mapH = maps.mapAbh[hbH];
      int idx = (mapA * maps.mapNab + mapB) * maps.mapNh + mapH;
      if (!consideredAbh[idx]) {
        consideredAbh[idx] = 1;
        for (int c = 0; c < 3; ++c)
          g[3 * hbB + c] -= hbDcn[(hbB * n + hbH) * 3 + c] * zz;
      }
    }
  }
  (void)bondIdx;
  return true;
}

} // namespace Xtb
} // namespace Avogadro
