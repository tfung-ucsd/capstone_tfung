function [out_pos, out_orien, exitflag] = user_localization(pos, kpos, emiss, kemiss, angles, kangles,...
                                                  A, S, toa_meas, aoa_meas, rss_meas, O_toa, O_aoa, O_rss, ...
                                                  sigma_t, sigma_p, sigma_r, place1, place2, c, P0, np, d0);
% Description:
%   this user-defined m-file can be used to implemnt a localization
%   algorithm of the user's choice.  by default, SeNeLEx uses a maximum
%   likelihood estimator which uses all available information to compute
%   the ML position estimates.  alternatively, SeNeLEx can call this
%   user-defined function.
%
%                                                                   
% inputs:
%   pos,            [x1 x2 ... xn y1 y2 ... yn]'  initial estimates of node
%                   positions.  if true values are known, these are used
%                   instead and indicies added to kpos.  (meters) 
%                   [note: the user-defined function init_ml() is 
%                   responsible for initialization.  by default, this
%                   function passes (unrealistically) the true coordinates 
%                   as initial values.]
%   kpos,           indicies of known positions in pos vector.  max length
%                   is 2n.
%   emiss,          [t1 t2 t3 ...]  vector of known emission times for
%                   sources.  (seconds)
%   kemiss,         indicies of the nodes whose emission times are given in
%                   'emiss'.  that is, emission time emiss(i) corresponds
%                   to node( kemiss(i) ).  note, kemiss is an index to the
%                   total node array, NOT just the sources, S.
%   angles,         [a1 a2 a3 ...]  vector of known angles (orientations) 
%                   for arrays.  (radians)
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
%                   O_toa(i,j) = 1 if array i makes a time-of-arrival
%                   measurement of an emission from source j.  if "use time
%                   information" is not checked, this matrix will be all 
%                   zeros.
%   O_aoa,          |A|x|S| matrix indicating which aoa's were observed.
%                   1 if observed, 0 if not.  if "use angle information" is
%                   not checked, this matrix will be all zeros.
%   O_rss,          |A|x|S| matrix indicating which rss's were observed.
%                   1 if observed, 0 if not.  if "use RSS information" is
%                   not checked, this matrix will be all zeros.
%   sigma_t,        |A|x|S| matrix of the standard deviation of each toa
%                   measurement. (sec)
%   sigma_p,        |A|x|S| matrix of the standard deviation of each aoa
%                   measurement. (radians)
%   sigma_r,        |A|x|S| matrix of the standard deviation of each rss
%                   measurement. (dB)
%   place1          place holder, unused.
%   place2          place holder, unused.
%   c,              propagation velocity used in timing measurements (meters/sec)
%   P0,             power at reference distance, d0.  used in rss
%                   measurements.  (dBm)
%   np,             pathloss exponent used in rss measurements.  (unitless)
%   d0,             reference distance used in rss measurements.  (m)
%
% outputs:
%   out_pos,        estimate of the nodes' positions.  [x1 y1; x2 y2; ... xn yn]
%   out_orien,      estimate of the nodes' orientations [a1 a2 ... an]'
%                   (radians).  return [] if orientations cannot be
%                   determined.
%   exit_flag,      1=successful estimate returned, 0=estimation failed
%
%
% Author :  Josh Ash
%           The Ohio State Universtiy
%
% $Id: user_localization.m,v 1.1 2006/06/04 19:44:50 ashj Exp $

A_len = length(A);  % number of arrays (receivers)
S_len = length(S);  % number of sources (transmitters)
n = length(pos)/2;    % number of total nodes



% return the initial values slightly perturbed
out_pos = reshape(pos,n,2) + 2*randn(n,2);
out_orien = [];  % if not returning orientation, de-select angle info on main gui
exitflag = 1;



