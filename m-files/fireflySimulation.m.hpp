// Automatically translated using m2cpp 2.0 on 2019-04-12 16:39:20

#ifndef FIREFLYSIMULATION_M_HPP
#define FIREFLYSIMULATION_M_HPP

#include "mconvert.h"
#include <iostream>
#include "SPlot.h"
#include <armadillo>
using namespace arma ;

struct _P
{
  TYPE KeepUnmatched, Results ;
} ;

void fireflySimulation(TYPE varargin, TYPE& phaseOut, TYPE& omegaOut)
{
  SPlot _plot ;
  TYPE Alpha, Beta, DoublingThreshold, E, FilterLength, PlotLength, SimulationLength, animSpeed, animate, defaultAlpha, defaultAnimSpeed, defaultBeta, defaultDoublingThreshold, defaultFilterLength, defaultFilterType, defaultNodes, defaultPlotLength, defaultSampleRate, defaultSeed, defaultSimulationLength, defaulttref, domega, f, filterType, fire, i, in, interactions, j, k, k2, k3, k4, l1, l3, l4, n, omega, omegas, p, p1, p2, pf, phi, phiref, plotSpectrogram, q, r, randomSeed, sortedOmegas, sr, tmp, tref ;
  rowvec _aux_rowvec_1, _aux_rowvec_10, _aux_rowvec_11, _aux_rowvec_12, _aux_rowvec_13, _aux_rowvec_14, _aux_rowvec_2, _aux_rowvec_3, _aux_rowvec_4, _aux_rowvec_5, _aux_rowvec_6, _aux_rowvec_7, _aux_rowvec_8, _aux_rowvec_9 ;
  wall_clock _timer ;
  p = inputParser ;
  p.KeepUnmatched = 1 ;
  defaultSimulationLength = 60000 ;
  defaultFilterLength = 8 ;
  defaultNodes = 5 ;
  defaultSampleRate = 1000 ;
  defaultPlotLength = 10000 ;
  defaulttref = 5 ;
  defaultDoublingThreshold = 20 ;
  defaultSeed = 0 ;
  defaultAnimSpeed = 0.2 ;
  defaultAlpha = 0.4 ;
  defaultBeta = 0.4 ;
  defaultFilterType = "median" ;
  addOptional(p, "SimulationLength", defaultSimulationLength, @isscalar) ;
  addOptional(p, "FilterLength", defaultFilterLength, @isscalar) ;
  addOptional(p, "Nodes", defaultNodes, @isscalar) ;
  addOptional(p, "SampleRate", defaultSampleRate, @isscalar) ;
  addOptional(p, "PlotLength", defaultPlotLength, @isscalar) ;
  addOptional(p, "tref", defaulttref, @isscalar) ;
  addOptional(p, "DoublingThreshold", defaultDoublingThreshold, @isscalar) ;
  addOptional(p, "alpha", defaultAlpha, @isscalar) ;
  addOptional(p, "beta", defaultBeta, @isscalar) ;
  addOptional(p, "PlotSpectrogram", 0, @isscalar) ;
  addOptional(p, "RandomSeed", defaultSeed, @isscalar) ;
  addOptional(p, "AnimSpeed", defaultAnimSpeed, @isscalar) ;
  addOptional(p, "FilterType", defaultFilterType) ;
  addOptional(p, "Animate", 0) ;
  parse(p, varargin{m2cpp::span<uvec>(0, varargin.n_rows-1)}) ;
  SimulationLength = p.Results ;
  FilterLength = p.Results ;
  n = p.Results ;
  sr = p.Results ;
  PlotLength = p.Results ;
  tref = p.Results ;
  DoublingThreshold = p.Results ;
  Alpha = p.Results ;
  Beta = p.Results ;
  plotSpectrogram = p.Results ;
  randomSeed = p.Results ;
  filterType = str2func(p.Results) ;
  animate = p.Results ;
  animSpeed = p.Results*200 ;
  if (randomSeed!=0)
  {
    rng(randomSeed) ;
  }
  omegas = arma::zeros<mat>(SimulationLength, n) ;
  phi = arma::zeros<mat>(SimulationLength, n) ;
  phi.col(0).rows(m2cpp::fspan(1, 1, n)) = arma::randu<mat>(n, 1) ;
  omega = 2*arma::pow(2., ((arma::randu<mat>(1, n)-0.5)*2)) ;
  E = arma::ones<mat>(SimulationLength, n) ;
  f = arma::zeros<rowvec>(n) ;
  domega = arma::zeros<mat>(n, 2) ;
  in = 0 ;
  for (i=2; i<=SimulationLength; i++)
  {
    fire = 0 ;
    for (j=1; j<=n; j++)
    {
      if (phi(i-1, j)<1)
      {
        phi(i, j) = phi(i-1, j)+omega(j)/sr ;
      }
      else
      {
        in = in+1 ;
        interactions(in, arma::strans(arma::join_rows(m2cpp::fspan(1, 1, 3), m2cpp::srow<double>(11)))-1) = {arma::join_rows(arma::join_rows(arma::join_rows(m2cpp::srow<sword>(2), i), j), f(j))} ;
        f(j) = mod(f(j)+1, 2) ;
      }
    }
    phiref = tref*omega*0.0005 ;
    for (int _j=0; _j<length(randperm(n)); _j++)
    {
      j = randperm(n)[_j] ;
      if (phi(i, j)>=1)
      {
        if (f(j)==0)
        {
          in = in+1 ;
          interactions(in, m2cpp::fspan(1, 1, 11)) = {arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(m2cpp::srow<sword>(-1), i), j), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0))} ;
        }
        if (f(j)==1 && fire==1)
        {
          in = in+1 ;
          interactions(in, m2cpp::fspan(1, 1, 3)) = {arma::join_rows(arma::join_rows(m2cpp::srow<sword>(0), i), j)} ;
        }
        if (f(j)==1 && fire==0)
        {
          fire = 1 ;
          phi(i, j) = 1 ;
          for (int _k=0; _k<length(setdiff(randperm(n), j, "stable")); _k++)
          {
            k = setdiff(randperm(n), j, "stable")[_k] ;
            if (phi(i, k)>phiref(k))
            {
              pf = 5 ;
              if (1 == pf)
              {
                phi(i, k) = phi(i, k)+Alpha*phi(i, k) ;
              }
              else if (2 == pf)
              {
                phi(i, k) = phi(i, k)+Alpha*phi(i, k)*(phi(i, k)-0.5)/abs(phi(i, k)-0.5) ;
              }
              else if (3 == pf)
              {
                phi(i, k) = phi(i, k)+Alpha*phi(i, k)-arma::sin(phi(i, k)*2*datum::pi) ;
              }
              else if (4 == pf)
              {
                phi(i, k) = phi(i, k)+Alpha-arma::sin(phi(i, k)*2*datum::pi) ;
              }
              else if (5 == pf)
              {
                phi(i, k) = phi(i, k)+Alpha-arma::sin(phi(i, k)*2*datum::pi)*abs(-arma::sin(phi(i, k)*2*datum::pi)) ;
              }
              else if (6 == pf)
              {
                phi(i, k) = phi(i, k)+Alpha*phi(i, k)-arma::sin(phi(i, k)*2*datum::pi)*abs(-arma::sin(phi(i, k)*2*datum::pi)) ;
              }
              else if (7 == pf)
              {
                phi(i, k) = phi(i, k)+Alpha-arma::cos(phi(i, k)*datum::pi)*abs(-arma::cos(phi(i, k)*datum::pi)) ;
              }
              else if (8 == pf)
              {
                phi(i, k) = phi(i, k)+Alpha*pow((phi(i, k)-0.5), 3) ;
              }
              if (phi(i, k)<0)
              {
                phi(i, k) = 0 ;
              }
              else if (phi(i, k)>=1)
              {
                phi(i, k) = 1 ;
              }
              if (phi(i, k)<phiref(k) || phi(i, k)>1-phiref(k))
              {
                E(m2cpp::fspan(i, 1, SimulationLength), k) = 0 ;
              }
              else
              {
                E(m2cpp::fspan(i, 1, SimulationLength), k) = (1-arma::cos(2*datum::pi*phi(i, k)))/2.0 ;
              }
              domega(k, 0) = ((omega(k)*pow(2, (Beta-arma::sin(2*datum::pi*phi(i, k))*(filterType(E(m2cpp::fspan(m2cpp::length(E.cols(k))-FilterLength+1, 1, m2cpp::length(E.cols(k))), k))))))+prod(domega(k, m2cpp::fspan(1, 1, 2))))/(domega(k, 2)+1) ;
              domega(k, 1) = domega(k, 2)+1 ;
              if (domega(k, 2)>DoublingThreshold)
              {
                omega(k) = 2*omega(k) ;
                domega.rows(k) = {0, 0} ;
                std::cout << {arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(m2cpp::srow(std::string("Someone''s stuck. Doubling frequency of fly ")), num2str(k)), std::string("")), Time), std::string("")), num2str(i/sr)), m2cpp::srow(std::string("s")))} << std::endl ;
              }
              in = in+1 ;
              interactions(in, m2cpp::fspan(1, 1, 9)) = {arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(m2cpp::srow<sword>(1), i), j), k), m2cpp::srow<sword>(1)), filterType(E(m2cpp::fspan(m2cpp::length(E.cols(k))-FilterLength+1, 1, m2cpp::length(E.cols(k))), k))), phi(i, k)-phi(i-1, k)), Beta*arma::sin(2*datum::pi*phi(i, k))), domega(k, 1))} ;
            }
            else
            {
              in = in+1 ;
              interactions(in, m2cpp::fspan(1, 1, 6)) = {arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(m2cpp::srow<sword>(1), i), j), k), m2cpp::srow<sword>(0)), filterType(E(m2cpp::fspan(m2cpp::length(E.cols(k))-FilterLength+1, 1, m2cpp::length(E.cols(k))), k)))} ;
            }
          }
        }
        if (domega(j, 2)>0)
        {
          omega(j) = domega(j, 1) ;
          domega.rows(j) = {0, 0} ;
          if (omega(j)==0)
          {
            std::cout << "a node has died... omega = zero" << std::endl ;
          }
          in = in+1 ;
          interactions(in, m2cpp::fspan(1, 1, 10)) = {arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(m2cpp::srow<sword>(-2), i), j), j), m2cpp::srow<sword>(1)), filterType(E(m2cpp::fspan(m2cpp::length(E.cols(j))-FilterLength+1, 1, m2cpp::length(E.cols(j))), j))), m2cpp::srow<sword>(-1)), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0)), omega(j))} ;
        }
        if (arma::sum(E(m2cpp::fspan(m2cpp::length(E.cols(j))-FilterLength+1, 1, m2cpp::length(E.cols(j))), j))==0 && m2cpp::length(interactions(interactions.col(2)==j && interactions.col(10)==1, 1))>8)
        {
          sortedOmegas = arma::trans(sort(interactions(interactions.col(2)==j && interactions.col(10)==1, 2), "descend")) ;
          omega(j) = -2000*1.0/mean(diff(sortedOmegas(m2cpp::fspan(1, 1, FilterLength)))) ;
          in = in+1 ;
          interactions(in, m2cpp::fspan(1, 1, 10)) = {arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(arma::join_rows(m2cpp::srow<sword>(-3), i), j), j), m2cpp::srow<sword>(1)), filterType(E(m2cpp::fspan(m2cpp::length(E.cols(j))-FilterLength+1, 1, m2cpp::length(E.cols(j))), j))), m2cpp::srow<sword>(-1)), m2cpp::srow<sword>(0)), m2cpp::srow<sword>(0)), omega(j))} ;
        }
      }
    }
    for (j=1; j<=n; j++)
    {
      if (phi(i, j)>=1)
      {
        phi(i, j) = 1 ;
      }
      omegas(i, j) = omega(j) ;
    }
  }
  if (plotSpectrogram==1)
  {
    _plot.figure(1) ;
    _plot.clf(1) ;
    _plot.subplot(5, 4, m2cpp::fspan(1, 1, 2)) ;
    _plot.plot((m2cpp::fspan(1, 1, PlotLength))/sr, phi.rows(m2cpp::fspan(1, 1, PlotLength))) ;
    _plot.title("phase: beginning") ;
    _plot.subplot(5, 4, m2cpp::fspan(3, 1, 4)) ;
    _plot.plot((m2cpp::fspan(SimulationLength-PlotLength, 1, SimulationLength))/sr, phi.rows(m2cpp::fspan(SimulationLength-PlotLength, 1, SimulationLength))) ;
    _plot.title("phase end") ;
    _plot.subplot(5, 4, {5, 10}) ;
    _plot.plot((m2cpp::fspan(1, 1, SimulationLength))/sr, phi.rows(m2cpp::fspan(1, 1, SimulationLength))+1.05*cumsum(arma::ones<umat>(n, SimulationLength))) ;
    _plot.title("phase") ;
    _plot.xlim({arma::join_rows(m2cpp::srow<sword>(0), PlotLength)}/sr) ;
    _plot.hold(1) ;
    scatter(interactions(interactions.col(10)==1, 2)/1000.0, (interactions(interactions.col(10)==1, 3)-1)*1.05, "k.") ;
    _plot.subplot(5, 4, {7, 12}) ;
    double __aux_rowvec_1 [] = {1, -1} ;
    _aux_rowvec_1 = rowvec(__aux_rowvec_1, 2, false) ;
    _plot.plot((m2cpp::fspan(1, 1, SimulationLength))/sr, filter(_aux_rowvec_1, 1, phi.rows(m2cpp::fspan(1, 1, SimulationLength)))+1.15*cumsum(arma::ones<umat>(n, SimulationLength))) ;
    _plot.title("delta phase") ;
    _plot.xlim({arma::join_rows(m2cpp::srow<sword>(0), PlotLength)}/sr) ;
    _plot.subplot(5, 4, {13, 18}) ;
    _plot.plot((m2cpp::fspan(1, 1, SimulationLength))/1000.0, omegas) ;
    _plot.title("freq") ;
    _plot.subplot(5, 4, {15, 20}) ;
    double __aux_rowvec_2 [] = {1, -1} ;
    _aux_rowvec_2 = rowvec(__aux_rowvec_2, 2, false) ;
    tmp = filter(_aux_rowvec_2, 1, omegas) ;
    tmp.rows(m2cpp::fspan(1, 1, n)).reset() ;
    _plot.plot((m2cpp::fspan(n+1, 1, SimulationLength))/1000.0, tmp) ;
    _plot.title("delta freq") ;
    _plot.figure(2) ;
    for (i=1; i<=n; i++)
    {
      _plot.subplot(n, 1, i) ;
      spectrogram(arma::sin(2*datum::pi*phi.cols(i)), 4096, 2048, 10000, sr, "yaxis") ;
      _plot.ylim(0, 8) ;
    }
  }
  if (nargout>0)
  {
    phaseOut = phi ;
    if (nargout>1 || animate==1)
    {
      omegaOut = omegas ;
    }
  }
  if (animate==1)
  {
    _plot.subplot(6, 2, m2cpp::fspan(1, 1, 4)) ;
    _plot.plot((m2cpp::fspan(1, 1, m2cpp::length(phaseOut.col(0))))/sr, phaseOut) ;
    _plot.xlabel("time (s)") ;
    _plot.ylabel("phase") ;
    double __aux_rowvec_3 [] = {0, 0} ;
    _aux_rowvec_3 = rowvec(__aux_rowvec_3, 2, false) ;
    double __aux_rowvec_4 [] = {0, 1} ;
    _aux_rowvec_4 = rowvec(__aux_rowvec_4, 2, false) ;
    double __aux_rowvec_5 [] = {0, 0, 0} ;
    _aux_rowvec_5 = rowvec(__aux_rowvec_5, 3, false) ;
    l1 = line(_aux_rowvec_3, _aux_rowvec_4, "LineWidth", 2, "color", _aux_rowvec_5) ;
    k3 = _plot.subplot(6, 2, {6, 8}) ;
    _plot.plot((m2cpp::fspan(1, 1, m2cpp::length(phaseOut.col(0))))/sr, phaseOut) ;
    _plot.xlabel("time (s)") ;
    _plot.ylabel("phase") ;
    _plot.hold(1) ;
    double __aux_rowvec_6 [] = {0, 0} ;
    _aux_rowvec_6 = rowvec(__aux_rowvec_6, 2, false) ;
    double __aux_rowvec_7 [] = {0, 1} ;
    _aux_rowvec_7 = rowvec(__aux_rowvec_7, 2, false) ;
    double __aux_rowvec_8 [] = {0, 0, 0} ;
    _aux_rowvec_8 = rowvec(__aux_rowvec_8, 3, false) ;
    l3 = line(_aux_rowvec_6, _aux_rowvec_7, "LineWidth", 2, "color", _aux_rowvec_8) ;
    k4 = _plot.subplot(6, 2, {10, 12}) ;
    _plot.plot((m2cpp::fspan(1, 1, m2cpp::length(omegaOut.col(0))))/sr, omegaOut) ;
    _plot.xlabel("time (s)") ;
    _plot.ylabel("frequency") ;
    _plot.hold(1) ;
    double __aux_rowvec_9 [] = {0, 0} ;
    _aux_rowvec_9 = rowvec(__aux_rowvec_9, 2, false) ;
    double __aux_rowvec_10 [] = {0, 4} ;
    _aux_rowvec_10 = rowvec(__aux_rowvec_10, 2, false) ;
    double __aux_rowvec_11 [] = {0, 0, 0} ;
    _aux_rowvec_11 = rowvec(__aux_rowvec_11, 3, false) ;
    l4 = line(_aux_rowvec_9, _aux_rowvec_10, "LineWidth", 2, "color", _aux_rowvec_11) ;
    k2 = _plot.subplot(6, 2, m2cpp::fspan(5, 2, 11)) ;
    for (i=1+5*animSpeed; i<=m2cpp::length(phaseOut.col(0)); i+=animSpeed)
    {
      m2cpp::tic() ;
      _plot.xlim(k3, {arma::join_rows(i-sr, i+sr)}/sr) ;
      _plot.xlim(k4, {arma::join_rows(i-sr, i+sr)}/sr) ;
      polar(datum::pi, 1) ;
      _plot.hold(1) ;
      p1 = polar(k2, phaseOut.rows(m2cpp::fspan(i-5*animSpeed, 1, i))*2*datum::pi, 1/omegaOut.rows(m2cpp::fspan(i-5*animSpeed, 1, i)), "-") ;
      double __aux_rowvec_12 [] = {0, 3} ;
      _aux_rowvec_12 = rowvec(__aux_rowvec_12, 2, false) ;
      double __aux_rowvec_13 [] = {0, 0} ;
      _aux_rowvec_13 = rowvec(__aux_rowvec_13, 2, false) ;
      double __aux_rowvec_14 [] = {0, 0, 0} ;
      _aux_rowvec_14 = rowvec(__aux_rowvec_14, 3, false) ;
      line(_aux_rowvec_12, _aux_rowvec_13, "LineWidth", 3, "color", _aux_rowvec_14) ;
      p2 = polar(k2, phaseOut.rows(m2cpp::fspan(i, 1, i+1))*2*datum::pi, 1/omegaOut.rows(m2cpp::fspan(i, 1, i+1)), ".") ;
      set(l1, "xdata", {arma::join_rows(i/sr, i/sr)}) ;
      set(l3, "xdata", {arma::join_rows(i/sr, i/sr)}) ;
      set(l4, "xdata", {arma::join_rows(i/sr, i/sr)}) ;
      set(p1, "LineWidth", 1) ;
      for (q=1; q<=n; q++)
      {
        r = interactions(interactions.col(0)==2 && interactions.col(1)>i && interactions.col(2)==q, 11) ;
        if (!m2cpp::isempty(r))
        {
          if (r(1)==1)
          {
            set(p2(q), "MarkerSize", 5, "Marker", "o") ;
          }
          else
          {
            set(p2(q), "MarkerSize", 18) ;
          }
        }
      }
      _plot.hold(0) ;
      drawnow ;
      pause(0.2-m2cpp::toc()) ;
    }
  }
  _plot.show() ;
}
#endif