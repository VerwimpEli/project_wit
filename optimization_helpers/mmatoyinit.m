%-------------------------------------------------------------
%
%    Copyright (C) 2009 Krister Svanberg
%
%    This file, mmatoyinit.m, is part of GCMMA-MMA-code.
%    
%    GCMMA-MMA-code is free software; you can redistribute it and/or
%    modify it under the terms of the GNU General Public License as 
%    published by the Free Software Foundation; either version 3 of 
%    the License, or (at your option) any later version.
%    
%    This code is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%    
%    You should have received a copy of the GNU General Public License
%    (file COPYING) along with this file.  If not, see 
%    <http://www.gnu.org/licenses/>.
%    
%    You should have received a file README along with this file,
%    containing contact information.  If not, see
%    <http://www.smoptit.se/> or e-mail mmainfo@smoptit.se or krille@math.kth.se.
%
%------
%
%  Version September 2009.

m = 1;
n = 25*25;
epsimin = 0.0000001;
xval    = rand(n, 1) * 0.4;
xold1   = xval;
xold2   = xval;
xmin    = zeros(n, 1);
xmax    = ones(n, 1);
low     = xmin;
upp     = xmax;
c       = 1e8 * [1]';
d       = [1]';
a0      = 1;
a       = [0]';
outeriter = 0;
maxoutit  = 15;
kkttol  = 1e-8;
%
%---------------------------------------------------------------------
