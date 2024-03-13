
/***************************************************************

 create_TP06_lookup_OpSplit_2D.c

 Version       2.1

 Date          23-oct-2023

 Author        R.H.Clayton (r.h.clayton@sheffield.ac.uk)
 
 This file is part of VentricularFibrosis.

 Copyright (c) Richard Clayton, 
 Department of Computer Science, 
 University of Sheffield, 2006, 2009, 2016, 2023

 VentricularFibrosis is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.

*********************************************************************/
#include <TP06_OpSplit_2D.h>

int create_TP06_lookup_OpSplit_2D( double **lookup )
{

/* This function precalculates variables for lookup tables  */
/* we assume that the ODEs can be approximated by the equation */
/* described by Rush and Larsen                                */
/*                                                             */
/* n = n_inf - (n_inf - n) * exp( -dt / tau_n )                */
/*                                                             */
/* and precalculate n_inf and n_exp = exp( -dt / tau_n )       */
/* for Vm = -100 mV to + 100 mV with a step of 0.1 mV          */

  /* indices for lookup table */
  int na_h_exp       = 1;
  int na_h_inf       = 2;
  int na_j_exp       = 3;
  int na_j_inf       = 4;
  int na_m_exp       = 5;
  int na_m_inf       = 6;
  int ca_d_exp       = 7;
  int ca_d_inf       = 8;
  int ca_f_exp       = 9;
  int ca_f_inf       = 10;
  int ca_f2_exp      = 11;
  int ca_f2_inf      = 12;
  int k_xr1_exp      = 13;
  int k_xr1_inf      = 14;
  int k_xr2_exp      = 15;
  int k_xr2_inf      = 16;
  int k_xs_exp       = 17;
  int k_xs_inf       = 18;
  int to_r_epi_inf   = 19;
  int to_r_epi_exp   = 20;
  int to_s_epi_inf   = 21;
  int to_s_epi_exp   = 22;
  int to_r_endo_inf  = 23;
  int to_r_endo_exp  = 24;
  int to_s_endo_inf  = 25;
  int to_s_endo_exp  = 26;
  int to_r_M_inf     = 27;
  int to_r_M_exp     = 28;
  int to_s_M_inf     = 29;
  int to_s_M_exp     = 30;

  float Vm;
  int Vindex;
  int gain = GAIN;
  int offset = VMOFFSET;

  double alpha_h;  /* Na alpha-h rate constant (ms^-1) */
  double beta_h;   /* Na beta-h rate constant (ms^-1) */
  double alpha_j;  /* Na alpha-j rate constant (ms^-1) */
  double beta_j;   /* Na beta-j rate constant (ms^-1) */
  double alpha_m;  /* Na alpha-m rate constant (ms^-1) */
  double beta_m;   /* Na beta-m rate constant (ms^-1) */
  double tau_h, h_inf;
  double tau_j, j_inf;
  double tau_m, m_inf;

  double d_inf;   /* Steady-state value of activation gate d  */
  double tau_d;   /* Time constant of gate d (ms^-1) */
  double ad, bd, cd;
  double f_inf;   /* Steady-state value of inactivation gate f */
  double tau_f;   /* Time constant of gate f (ms^-1) */
  double af, bf, cf;
  double f2_inf;   /* Steady-state value of inactivation gate f */
  double tau_f2;   /* Time constant of gate f (ms^-1) */
  double af2, bf2, cf2;

  double xr1_inf;  /* Steady-state value of inactivation gate xr1 */
  double tau_xr1;  /* Time constant of gate xr1 (ms^-1) */
  double axr1, bxr1;
  double xr2_inf;  /* Steady-state value of inactivation gate xr2 */
  double tau_xr2;  /* Time constant of gate xr2 (ms^-1) */
  double axr2, bxr2;

  double xs_inf; /* Steady-state value of inactivation gate xs */
  double tau_xs; /* Time constant of gate xs (ms^-1) */
  double axs, bxs;

  double r_inf_epi; /* Steady-state value of ito gate r for epi cells */
  double tau_r_epi; /* Time constant of gate r (ms^-1) */
  double s_inf_epi; /* Steady-state value of ito gate s for epi cells */
  double tau_s_epi; /* Time constant of gate s (ms^-1) */
  double r_inf_endo; /* Steady-state value of ito gate r for endo cells */
  double tau_r_endo; /* Time constant of gate r (ms^-1) */
  double s_inf_endo; /* Steady-state value of ito gate s for endo cells */
  double tau_s_endo; /* Time constant of gate s (ms^-1) */
  double r_inf_M; /* Steady-state value of ito gate r for M cells */
  double tau_r_M; /* Time constant of gate r (ms^-1) */
  double s_inf_M; /* Steady-state value of ito gate s for M cells */
  double tau_s_M; /* Time constant of gate s (ms^-1) */


  for (Vm = -100.0; Vm <= 100.0; Vm += 0.1)
    {
    /* 2001 voltage steps: gain should be 10, offset 1001 */

    Vindex = (int) Vm * gain + offset;

    /* Na current */
    alpha_m = 1.0/(1.0+exp((-60.0-Vm)/5.0));
	beta_m = 0.1/(1.0+exp((Vm+35.0)/5.0))+0.10/(1.0+exp((Vm-50.0)/200.0));
	tau_m = alpha_m * beta_m;
	m_inf = 1.0/((1.0+exp((-56.86-Vm)/9.03))*(1.0+exp((-56.86-Vm)/9.03)));

	if (Vm >= -40.0)
	    {
	    alpha_h = 0.0;
	    beta_h = (0.77/(0.13*(1.0 + exp(-(Vm + 10.66)/11.1))));
	    }
	else
	    {
	    alpha_h = (0.057 * exp(-(Vm + 80.0)/6.8));
	    beta_h = (2.7 * exp(0.079 * Vm)+(3.1e5) * exp(0.3485 * Vm));
	    }
	tau_h = 1.0/(alpha_h + beta_h);
	h_inf = 1.0/((1.0 + exp((Vm + 71.55)/7.43))*(1.0 + exp((Vm + 71.55)/7.43)));

	if(Vm >= -40.0)
	    {
	    alpha_j = 0.;
	    beta_j = (0.6 * exp((0.057) * Vm)/(1.0 + exp(-0.1 * (Vm + 32.0))));
	    }
	else
	    {
		alpha_j = (((-2.5428e4)*exp(0.2444*Vm)-(6.948e-6)*exp(-0.04391*Vm))*(Vm+37.78)/(1.0+exp(0.311*(Vm+79.23))));
		beta_j = (0.02424*exp(-0.01052*Vm)/(1.0+exp(-0.1378*(Vm+40.14))));
	    }
	tau_j = 1.0/(alpha_j + beta_j);
	j_inf = h_inf;

	lookup[na_h_exp][Vindex] = tau_h;
    lookup[na_h_inf][Vindex] = h_inf;

    lookup[na_j_exp][Vindex] = tau_j;
    lookup[na_j_inf][Vindex] = j_inf;

    lookup[na_m_exp][Vindex] = tau_m;
    lookup[na_m_inf][Vindex] = m_inf;

    /* Ca currents */

    d_inf = 1.0/(1.0+exp((-8.0-Vm)/7.5));
    ad = 1.4/(1.0+exp((-35.0-Vm)/13.0))+0.25;
    bd=1.4/(1.0+exp((Vm+5.0)/5.0));
	cd=1.0/(1.0+exp((50.0-Vm)/20.0));
	// bug fixed
    tau_d = ad * bd + cd;

    f_inf = 1.0/(1.0+exp((Vm+20.0)/7.0));
    af=1102.5*exp(-(Vm+27.0)*(Vm+27.0)/225.0);
	bf=200.0/(1.0+exp((13.0-Vm)/10.0));
	cf=(180.0/(1.0+exp((Vm+30.0)/10.0)))+20.0;
    tau_f = af + bf + cf;

    f2_inf = 0.67/(1.0+exp((Vm+35.0)/7.0))+0.33;
    af2=562.0*exp(-(Vm+27.0)*(Vm+27.0)/240.0);
	bf2=31.0/(1.0+exp((25.0-Vm)/10.0));
	cf2=80.0/(1.0+exp((Vm+30.0)/10.0));
	// bug fixed
	tau_f2 = af2 + bf2 + cf2;

    lookup[ca_d_exp][Vindex] = tau_d;
    lookup[ca_d_inf][Vindex] = d_inf;

    lookup[ca_f_exp][Vindex] = tau_f;
    lookup[ca_f_inf][Vindex] = f_inf;

	lookup[ca_f2_exp][Vindex] = tau_f2;
    lookup[ca_f2_inf][Vindex] = f2_inf;

    /* Rapidly inactivating K current IKr */

    xr1_inf = 1.0/(1.0+exp((-26.0-Vm)/7.0));
    axr1 = 450.0/(1.0+exp((-45.0-Vm)/10.0));
    bxr1 = 6.0/(1.0+exp((Vm+30.0)/11.5));
    tau_xr1 = axr1 * bxr1;

    xr2_inf = 1.0/(1.0+exp((Vm+88.0)/24.0));
    axr2 = 3.0/(1.0+exp((-60.0-Vm)/20.0));
    bxr2 = 1.12/(1.0+exp((Vm-60.0)/20.0));
    tau_xr2 = axr2 * bxr2;

    lookup[k_xr1_exp][Vindex] = tau_xr1;
    lookup[k_xr1_inf][Vindex] = xr1_inf;
    lookup[k_xr2_exp][Vindex] = tau_xr2;
    lookup[k_xr2_inf][Vindex] = xr2_inf;

	/* Slowly inactivating K current IKs */

	xs_inf = 1.0/(1.0+exp((-5.0-Vm)/14.0));
	axs = (1400.0/(sqrt(1.0+exp((5.0-Vm)/6.0))));
	bxs = (1.0/(1.0+exp((Vm-35.0)/15.0)));
	tau_xs = axs * bxs + 80.0;

    lookup[k_xs_exp][Vindex] = tau_xs;
    lookup[k_xs_inf][Vindex] = xs_inf;

    /* Ito */

    // Code for EPI cells
	r_inf_epi = 1.0/(1.0+exp((20.0-Vm)/6.0));
	s_inf_epi = 1.0/(1.0+exp((Vm+20.0)/5.0));
	tau_r_epi = 9.5*exp(-(Vm+40.0)*(Vm+40.0)/1800.0)+0.8;
	tau_s_epi = 85.0*exp(-(Vm+45.0)*(Vm+45.0)/320.0)+5.0/(1.0+exp((Vm-20.0)/5.0))+3.0;

	// Code for ENDO cells
	r_inf_endo = 1.0/(1.0+exp((20.0-Vm)/6.0));
	s_inf_endo = 1.0/(1.0+exp((Vm+28.0)/5.0));
	tau_r_endo = 9.5*exp(-(Vm+40.0)*(Vm+40.0)/1800.0)+0.8;
	tau_s_endo = 1000.0*exp(-(Vm+67.0)*(Vm+67.0)/1000.0)+8.0;

	//code for M cells
	r_inf_M = 1.0/(1.0+exp((20.0-Vm)/6.0));
	s_inf_M = 1.0/(1.0+exp((Vm+20.0)/5.0));
	tau_r_M = 9.5*exp(-(Vm+40.0)*(Vm+40.0)/1800.0)+0.8;
	tau_s_M = 85.0*exp(-(Vm+45.0)*(Vm+45.0)/320.0)+5.0/(1.0+exp((Vm-20.0)/5.0))+3.0;

    lookup[to_r_epi_exp][Vindex] = tau_r_epi;
    lookup[to_r_epi_inf][Vindex] = r_inf_epi;
    lookup[to_s_epi_exp][Vindex] = tau_s_epi;
    lookup[to_s_epi_inf][Vindex] = s_inf_epi;

    lookup[to_r_endo_exp][Vindex] = tau_r_endo;
    lookup[to_r_endo_inf][Vindex] = r_inf_endo;
    lookup[to_s_endo_exp][Vindex] = tau_s_endo;
    lookup[to_s_endo_inf][Vindex] = s_inf_endo;

    lookup[to_r_M_exp][Vindex] = tau_r_M;
    lookup[to_r_M_inf][Vindex] = r_inf_M;
    lookup[to_s_M_exp][Vindex] = tau_s_M;
    lookup[to_s_M_inf][Vindex] = s_inf_M;

    //prints lookup table on stdout
    //printf("%5.1f %d ",Vm,Vindex);
    //for(xx=1;xx<=37;xx++)
    //  printf("%7.6f ",lookup[xx][Vindex]);
    //printf("\n");

  } // end of for Vm loop

  return (1); // success
}
