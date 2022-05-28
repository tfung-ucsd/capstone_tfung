function [init_pos, exit_success] = init_ml(true_pos, kpos, emiss, kemiss, angles, kangles,...
                                         A, S, toa_meas, aoa_meas, rss_meas, O_toa, O_aoa, O_rss, ...
                                         sigma_t, sigma_p, sigma_r, c, P0, np, d0);

% Description:  This function provides the initial cordinate estimates to
%   the ML algorithm which will then iterate to find the optimal solution.
%   This GUI supports a large number of options, including time-of-arrival
%   estimates with and without known emission time, angle-of-arrival
%   estimates with and without known orientations, received signal strength
%   measurements, and combinations thereof plus support for sparce networks.
%
%   A single initialization routine supporting all these scenarios is
%   unknown, and, for demonstration purposes, we initialize with true values.
%   The user can replace this initialization routine with their own more
%   realistic version utilizing the input parameters as described below.
%
%
% Inputs:
%   true_pos,       [x1 y2; x2 y2; ... ; xn yn]  true node positions.
%                   elements of 'true_pos' indexed by 'kpos' are
%                   are legitimately known to the init routine.
%   kpos,           indicies of known positions in true_pos vector.
%                   max length is 2n.
%   emiss,          [t1 t2 t3 ...]  vector of known emission times for
%                   sources.  (seconds)
%   kemiss,         indicies of the nodes whose emission times are given in
%                   'emiss'.  that is, emission time emiss(i) corresponds
%                   to kemiss(i).  note, kemiss is an index to the
%                   total node array, NOT just the sources, S.
%   angles,         [a1 a2 a3 ...]  vector of known angles for arrays.
%                   (radians)
%   kangles,        indicies of the nodes whose known angles are given in
%                   'angles'.
%   A,              indicies of the arrays. max of n.  the length of A is
%                   denoted |A|.
%   S,              indicies of the sources. max of n.  if a node is a
%                   source and an array (it is colocated), then its index
%                   will be a member of both A and S.
%   toa_meas,       these are the measured arrival times
%                   (emission time + flight time + noise) from each
%                   array due to each source.  size is |S|x|A|. (seconds)
%   aoa_meas,       measured angles of arrival, |S|x|A|. (radians)
%   rss_meas,       measured rss values, |S|x|A|. (dBm)
%   O_toa,          |A|x|S| matrix indicating which toa's were observed.
%                   1 if observed, 0 if not.
%   O_aoa,          |A|x|S| matrix indicating which aoa's were observed.
%                   1 if observed, 0 if not.
%   O_rss,          |A|x|S| matrix indicating which rss's were observed.
%                   1 if observed, 0 if not.
%   sigma_t,        |A|x|S| matrix of the standard deviation of each toa
%                   measurement. (sec)
%   sigma_p,        |A|x|S| matrix of the standard deviation of each aoa
%                   measurement. (radians)
%   sigma_r,        |A|x|S| matrix of the standard deviation of each rss
%                   measurement. (dB)
%   c,              propagation velocity used in timing measurements (meters/sec)
%   P0,             power at reference distance, d0.  used in rss
%                   measurements.  (dBm)
%   np,             pathloss exponent used in rss measurements.  (unitless)
%   d0,             reference distance used in rss measurements.  (dBm)
%
%
% Outputs:
%   init_pos,       initial estimate of the nodes' positions.  [x1 y2; x2 y2; ... ; xn yn]
%   exit_success,   exit_success=1 if the initialization routine was
%                   successful.  exit_cussess=0 if init failed.  GUI will
%                   give warning and abort localization in this case.
%
%
% Author :  Josh Ash
%           The Ohio State Universtiy
%
% $Id: init_ml.m,v 1.1 2006/04/13 23:00:08 ashj Exp $

init_pos = true_pos;
exit_success = 1;

