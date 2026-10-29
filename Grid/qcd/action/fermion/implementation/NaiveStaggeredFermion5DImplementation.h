/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: ./lib/qcd/action/fermion/implementation/NaiveStaggeredFermion5DImplementation.h

    Copyright (C) 2015

Author: Peter Boyle <paboyle@ph.ed.ac.uk>
Author: AzusaYamaguchi <ayamaguc@staffmail.ed.ac.uk>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program; if not, write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
#include <Grid/qcd/action/fermion/FermionCore.h>
#include <Grid/qcd/action/fermion/NaiveStaggeredFermion5D.h>
#include <Grid/perfmon/PerfCount.h>

#pragma once

NAMESPACE_BEGIN(Grid);

// 5d lattice for multiple right hand sides; the s direction is not checkerboarded
template<class Impl>
NaiveStaggeredFermion5D<Impl>::NaiveStaggeredFermion5D(GridCartesian         &FiveDimGrid,
                                                       GridRedBlackCartesian &FiveDimRedBlackGrid,
                                                       GridCartesian         &FourDimGrid,
                                                       GridRedBlackCartesian &FourDimRedBlackGrid,
                                                       RealD _mass,
                                                       RealD _c1, RealD _u0,
                                                       const ImplParams &p) :
  Kernels(p),
  _FiveDimGrid        (&FiveDimGrid),
  _FiveDimRedBlackGrid(&FiveDimRedBlackGrid),
  _FourDimGrid        (&FourDimGrid),
  _FourDimRedBlackGrid(&FourDimRedBlackGrid),
  Stencil    (&FiveDimGrid,npoint,Even,directions,displacements,p),
  StencilEven(&FiveDimRedBlackGrid,npoint,Even,directions,displacements,p), // source is Even
  StencilOdd (&FiveDimRedBlackGrid,npoint,Odd ,directions,displacements,p), // source is Odd
  mass(_mass),
  c1(_c1),
  u0(_u0),
  Umu(&FourDimGrid),
  UmuEven(&FourDimRedBlackGrid),
  UmuOdd (&FourDimRedBlackGrid),
  _tmp(&FiveDimRedBlackGrid)
{
  GRID_ASSERT(FiveDimGrid._ndimension==5);
  GRID_ASSERT(FourDimGrid._ndimension==4);
  GRID_ASSERT(FourDimRedBlackGrid._ndimension==4);
  GRID_ASSERT(FiveDimRedBlackGrid._ndimension==5);
  GRID_ASSERT(FiveDimRedBlackGrid._checker_dim==1); // Don't checker the s direction

  // extent of fifth dim and not spread out
  Ls=FiveDimGrid._fdimensions[0];
  GRID_ASSERT(FiveDimRedBlackGrid._fdimensions[0]==Ls);
  GRID_ASSERT(FiveDimGrid._processors[0]         ==1);
  GRID_ASSERT(FiveDimRedBlackGrid._processors[0] ==1);

  // Other dimensions must match the decomposition of the four-D fields
  for(int d=0;d<4;d++){
    GRID_ASSERT(FiveDimGrid._processors[d+1]         ==FourDimGrid._processors[d]);
    GRID_ASSERT(FiveDimRedBlackGrid._processors[d+1] ==FourDimGrid._processors[d]);
    GRID_ASSERT(FourDimRedBlackGrid._processors[d]   ==FourDimGrid._processors[d]);

    GRID_ASSERT(FiveDimGrid._fdimensions[d+1]        ==FourDimGrid._fdimensions[d]);
    GRID_ASSERT(FiveDimRedBlackGrid._fdimensions[d+1]==FourDimGrid._fdimensions[d]);
    GRID_ASSERT(FourDimRedBlackGrid._fdimensions[d]  ==FourDimGrid._fdimensions[d]);

    GRID_ASSERT(FiveDimGrid._simd_layout[d+1]        ==FourDimGrid._simd_layout[d]);
    GRID_ASSERT(FiveDimRedBlackGrid._simd_layout[d+1]==FourDimGrid._simd_layout[d]);
    GRID_ASSERT(FourDimRedBlackGrid._simd_layout[d]  ==FourDimGrid._simd_layout[d]);
  }

  if (Impl::LsVectorised) {

    int nsimd = Simd::Nsimd();

    // Dimension zero of the five-d is the Ls direction
    GRID_ASSERT(FiveDimGrid._simd_layout[0]        ==nsimd);
    GRID_ASSERT(FiveDimRedBlackGrid._simd_layout[0]==nsimd);

    for(int d=0;d<4;d++){
      GRID_ASSERT(FourDimGrid._simd_layout[d]==1);
      GRID_ASSERT(FourDimRedBlackGrid._simd_layout[d]==1);
      GRID_ASSERT(FiveDimRedBlackGrid._simd_layout[d+1]==1);
    }

  } else {

    // Dimension zero of the five-d is the Ls direction
    GRID_ASSERT(FiveDimRedBlackGrid._simd_layout[0]==1);
    GRID_ASSERT(FiveDimGrid._simd_layout[0]        ==1);

  }
  int LLs = FiveDimGrid._rdimensions[0];
  int vol4= FourDimGrid.oSites();
  Stencil.BuildSurfaceList(LLs,vol4);

  vol4=FourDimRedBlackGrid.oSites();
  StencilEven.BuildSurfaceList(LLs,vol4);
  StencilOdd.BuildSurfaceList(LLs,vol4);

  ZeroCounters();
}
template<class Impl>
NaiveStaggeredFermion5D<Impl>::NaiveStaggeredFermion5D(GaugeField &_U,
                                                       GridCartesian         &FiveDimGrid,
                                                       GridRedBlackCartesian &FiveDimRedBlackGrid,
                                                       GridCartesian         &FourDimGrid,
                                                       GridRedBlackCartesian &FourDimRedBlackGrid,
                                                       RealD _mass,
                                                       RealD _c1, RealD _u0,
                                                       const ImplParams &p) :
  NaiveStaggeredFermion5D(FiveDimGrid,FiveDimRedBlackGrid,
                          FourDimGrid,FourDimRedBlackGrid,
                          _mass,_c1,_u0,p)
{
  ImportGauge(_U);
}

template <class Impl>
void NaiveStaggeredFermion5D<Impl>::CopyGaugeCheckerboards(void)
{
  pickCheckerboard(Even, UmuEven, Umu);
  pickCheckerboard(Odd,  UmuOdd , Umu);
}
template<class Impl>
void NaiveStaggeredFermion5D<Impl>::ImportGauge(const GaugeField &_U)
{
  ////////////////////////////////////////////////////////
  // Double store as NaiveStaggeredFermion: the Naik slot is discarded
  ////////////////////////////////////////////////////////
  DoubledGaugeField _UUU(GaugeGrid());
  Impl::DoubleStore(GaugeGrid(), _UUU, Umu, _U, _U );

  ////////////////////////////////////////////////////////
  // Apply scale factors to get the right fermion Kinetic term
  // 0.5 ( U p(x+mu) - Udag(x-mu) p(x-mu) )
  ////////////////////////////////////////////////////////
  for (int mu = 0; mu < Nd; mu++) {

    auto U = PeekIndex<LorentzIndex>(Umu, mu);
    PokeIndex<LorentzIndex>(Umu, U*( 0.5*c1/u0), mu );

    U = PeekIndex<LorentzIndex>(Umu, mu+4);
    PokeIndex<LorentzIndex>(Umu, U*(-0.5*c1/u0), mu+4);
  }

  CopyGaugeCheckerboards();
}
template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopDir(const FermionField &in, FermionField &out,int dir5,int disp)
{
  // As NaiveStaggeredFermion: no single-direction kernel for the naive action
  GRID_ASSERT(0);
};

template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DerivInternal(StencilImpl & st,
                                                  DoubledGaugeField & U,
                                                  GaugeField &mat,
                                                  const FermionField &A,
                                                  const FermionField &B,
                                                  int dag)
{
  // No force terms in multi-rhs solver staggered
  GRID_ASSERT(0);
}

template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopDeriv(GaugeField &mat, const FermionField &A, const FermionField &B, int dag)
{
  GRID_ASSERT(0);
}

template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopDerivEO(GaugeField &mat, const FermionField &A, const FermionField &B, int dag)
{
  GRID_ASSERT(0);
}

template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopDerivOE(GaugeField &mat, const FermionField &A, const FermionField &B, int dag)
{
  GRID_ASSERT(0);
}

template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopInternal(StencilImpl & st,
                                                 DoubledGaugeField & U,
                                                 const FermionField &in, FermionField &out,int dag)
{
  DhopCalls++;
  DhopTotalTime -= usecond();
  if ( StaggeredKernelsStatic::Comms == StaggeredKernelsStatic::CommsAndCompute )
    DhopInternalOverlappedComms(st,U,in,out,dag);
  else
    DhopInternalSerialComms(st,U,in,out,dag);
  DhopTotalTime += usecond();
}

template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopInternalOverlappedComms(StencilImpl & st,
                                                                DoubledGaugeField & U,
                                                                const FermionField &in, FermionField &out,int dag)
{
  Compressor compressor;

  DhopFaceTime -= usecond();
  st.Prepare();
  st.HaloGather(in,compressor);

  std::vector<std::vector<CommsRequest_t> > requests;
  st.CommunicateBegin(requests);

  st.CommsMergeSHM(compressor);
  DhopFaceTime += usecond();

  DhopComputeTime -= usecond();
  {
    int interior=1;
    int exterior=0;
    Kernels::DhopNaive(st,U,in,out,dag,interior,exterior);
  }
  DhopComputeTime += usecond();

  DhopCommTime -= usecond();
  st.CommsMerge(compressor);
  st.CommunicateComplete(requests);
  DhopCommTime += usecond();

  DhopComputeTime2 -= usecond();
  {
    int interior=0;
    int exterior=1;
    Kernels::DhopNaive(st,U,in,out,dag,interior,exterior);
  }
  DhopComputeTime2 += usecond();
}

template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopInternalSerialComms(StencilImpl & st,
                                                            DoubledGaugeField & U,
                                                            const FermionField &in, FermionField &out,int dag)
{
  Compressor compressor;

  DhopCommTime -= usecond();
  st.HaloExchange(in,compressor);
  DhopCommTime += usecond();

  // Dhop takes the 4d grid from U, and makes a 5d index for fermion
  DhopComputeTime -= usecond();
  {
    int interior=1;
    int exterior=1;
    Kernels::DhopNaive(st,U,in,out,dag,interior,exterior);
  }
  DhopComputeTime += usecond();
}

template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopOE(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionRedBlackGrid());    // verifies half grid
  conformable(in.Grid(),out.Grid()); // drops the cb check

  GRID_ASSERT(in.Checkerboard()==Even);
  out.Checkerboard() = Odd;

  DhopInternal(StencilEven,UmuOdd,in,out,dag);
}
template<class Impl>
void NaiveStaggeredFermion5D<Impl>::DhopEO(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionRedBlackGrid());    // verifies half grid
  conformable(in.Grid(),out.Grid()); // drops the cb check

  GRID_ASSERT(in.Checkerboard()==Odd);
  out.Checkerboard() = Even;

  DhopInternal(StencilOdd,UmuEven,in,out,dag);
}
template<class Impl>
void NaiveStaggeredFermion5D<Impl>::Dhop(const FermionField &in, FermionField &out,int dag)
{
  conformable(in.Grid(),FermionGrid()); // verifies full grid
  conformable(in.Grid(),out.Grid());

  out.Checkerboard() = in.Checkerboard();

  DhopInternal(Stencil,Umu,in,out,dag);
}

/////////////////////////////////////////////////////////////////////////
// Performance monitoring
/////////////////////////////////////////////////////////////////////////
template<class Impl>
void NaiveStaggeredFermion5D<Impl>::Report(void)
{
  if ( DhopCalls > 0 ) {
    std::cout << GridLogMessage << "#### NaiveStaggeredFermion5D Dhop: " << DhopCalls << " calls, Ls " << Ls << std::endl;
    std::cout << GridLogMessage << "     total   " << DhopTotalTime / DhopCalls << " us/call" << std::endl;
    std::cout << GridLogMessage << "     face    " << DhopFaceTime / DhopCalls << " us/call" << std::endl;
    std::cout << GridLogMessage << "     comms   " << DhopCommTime / DhopCalls << " us/call" << std::endl;
    std::cout << GridLogMessage << "     compute " << (DhopComputeTime + DhopComputeTime2) / DhopCalls << " us/call" << std::endl;
  }
}
template<class Impl>
void NaiveStaggeredFermion5D<Impl>::ZeroCounters(void)
{
  DhopCalls        = 0;
  DhopTotalTime    = 0;
  DhopCommTime     = 0;
  DhopComputeTime  = 0;
  DhopComputeTime2 = 0;
  DhopFaceTime     = 0;
}

/////////////////////////////////////////////////////////////////////////
// Implement the general interface. Here we use SAME mass on all slices
/////////////////////////////////////////////////////////////////////////
template <class Impl>
void NaiveStaggeredFermion5D<Impl>::Mdir(const FermionField &in, FermionField &out, int dir, int disp)
{
  DhopDir(in, out, dir, disp);
}
template <class Impl>
void NaiveStaggeredFermion5D<Impl>::MdirAll(const FermionField &in, std::vector<FermionField> &out)
{
  GRID_ASSERT(0);
}
template <class Impl>
void NaiveStaggeredFermion5D<Impl>::M(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerNo);
  axpy(out, mass, in, out);
}

template <class Impl>
void NaiveStaggeredFermion5D<Impl>::Mdag(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  Dhop(in, out, DaggerYes);
  axpy(out, mass, in, out);
}

template <class Impl>
void NaiveStaggeredFermion5D<Impl>::Meooe(const FermionField &in, FermionField &out)
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerNo);
  } else {
    DhopOE(in, out, DaggerNo);
  }
}
template <class Impl>
void NaiveStaggeredFermion5D<Impl>::MeooeDag(const FermionField &in, FermionField &out)
{
  if (in.Checkerboard() == Odd) {
    DhopEO(in, out, DaggerYes);
  } else {
    DhopOE(in, out, DaggerYes);
  }
}

template <class Impl>
void NaiveStaggeredFermion5D<Impl>::Mooee(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  typename FermionField::scalar_type scal(mass);
  out = scal * in;
}

template <class Impl>
void NaiveStaggeredFermion5D<Impl>::MooeeDag(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  Mooee(in, out);
}

template <class Impl>
void NaiveStaggeredFermion5D<Impl>::MooeeInv(const FermionField &in, FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  out = (1.0 / (mass)) * in;
}

template <class Impl>
void NaiveStaggeredFermion5D<Impl>::MooeeInvDag(const FermionField &in,FermionField &out)
{
  out.Checkerboard() = in.Checkerboard();
  MooeeInv(in, out);
}

////////////////////////////////////////////////////////
// Conserved current - not yet implemented.
////////////////////////////////////////////////////////
template <class Impl>
void NaiveStaggeredFermion5D<Impl>::ContractConservedCurrent(PropagatorField &q_in_1,
                                                             PropagatorField &q_in_2,
                                                             PropagatorField &q_out,
                                                             PropagatorField &src,
                                                             Current curr_type,
                                                             unsigned int mu)
{
  GRID_ASSERT(0);
}

template <class Impl>
void NaiveStaggeredFermion5D<Impl>::SeqConservedCurrent(PropagatorField &q_in,
                                                        PropagatorField &q_out,
                                                        PropagatorField &src,
                                                        Current curr_type,
                                                        unsigned int mu,
                                                        unsigned int tmin,
                                                        unsigned int tmax,
                                                        ComplexField &lattice_cmplx)
{
  GRID_ASSERT(0);
}

NAMESPACE_END(Grid);
