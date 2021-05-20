function [x, w, n, xR] = getCubeCubaturePointsAndWeights(k)
    % Get cube cubature points and weights for a cubature of precision k,
    % taken from the book "P. Solin, K. Segeth and I. Dolezel: Higher-Order
    % Finite Element Methods", Chapman & Hall/CRC Press, 2003.
    %
    % SYNOPSIS:
    %
    %   [x, w, n, xR] = getCubeCubaturePointsAndWeights(k)
    %
    % PARAMETERS:
    %   k - cubature prescision
    %
    % RETURNS:
    %   x  - Cubature points
    %   w  - Cubature weights
    %   n  - Number of cubature points
    %   xR - Coordinates of reference cube

%{
Copyright 2009-2021 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}

    xR = [-1, -1, -1;
           1, -1, -1;
           1,  1, -1;
          -1,  1, -1;
          -1, -1,  1;
           1, -1,  1;
           1,  1,  1;
          -1,  1,  1];
    
    if k <= 1

        xw = [0.000000000000000  0.000000000000000  0.000000000000000  8.000000000000000];

    elseif k <= 3

        xw = [ 1.000000000000000   0.000000000000000   0.000000000000000  1.333333333333333
              1   0.000000000000000   0.000000000000000  1.333333333333333
               0.000000000000000   1.000000000000000   0.000000000000000  1.333333333333333
               0.000000000000000  -1.000000000000000   0.000000000000000  1.333333333333333
               0.000000000000000   0.000000000000000   1.000000000000000  1.333333333333333
               0.000000000000000   0.000000000000000  -1.000000000000000  1.333333333333333];

    elseif k <= 5

        xw = [ 0.795822425754221   0.000000000000000   0.000000000000000  0.886426592797784
              -0.795822425754221   0.000000000000000   0.000000000000000  0.886426592797784
               0.000000000000000   0.795822425754221   0.000000000000000  0.886426592797784
               0.000000000000000  -0.795822425754221   0.000000000000000  0.886426592797784
               0.000000000000000   0.000000000000000   0.795822425754221  0.886426592797784
               0.000000000000000   0.000000000000000  -0.795822425754221  0.886426592797784
               0.758786910639328   0.758786910639328   0.758786910639328  0.335180055401662
               0.758786910639328  -0.758786910639328   0.758786910639328  0.335180055401662
               0.758786910639328   0.758786910639328  -0.758786910639328  0.335180055401662
               0.758786910639328  -0.758786910639328  -0.758786910639328  0.335180055401662
              -0.758786910639328   0.758786910639328   0.758786910639328  0.335180055401662
              -0.758786910639328  -0.758786910639328   0.758786910639328  0.335180055401662
              -0.758786910639328   0.758786910639328  -0.758786910639328  0.335180055401662
              -0.758786910639328  -0.758786910639328  -0.758786910639328  0.335180055401662];

    elseif k <= 7

        xw = [ 0.000000000000000   0.000000000000000   0.000000000000000  0.788073482744211
               0.848418011472252   0.000000000000000   0.000000000000000  0.499369002307720
              -0.848418011472252   0.000000000000000   0.000000000000000  0.499369002307720
               0.000000000000000   0.848418011472252   0.000000000000000  0.499369002307720
               0.000000000000000  -0.848418011472252   0.000000000000000  0.499369002307720
               0.000000000000000   0.000000000000000   0.848418011472252  0.499369002307720
               0.000000000000000   0.000000000000000  -0.848418011472252  0.499369002307720
               0.652816472101691   0.652816472101691   0.652816472101691  0.478508449425127
               0.652816472101691  -0.652816472101691   0.652816472101691  0.478508449425127
               0.652816472101691   0.652816472101691  -0.652816472101691  0.478508449425127
               0.652816472101691  -0.652816472101691  -0.652816472101691  0.478508449425127
              -0.652816472101691   0.652816472101691   0.652816472101691  0.478508449425127
              -0.652816472101691  -0.652816472101691   0.652816472101691  0.478508449425127
              -0.652816472101691   0.652816472101691  -0.652816472101691  0.478508449425127
              -0.652816472101691  -0.652816472101691  -0.652816472101691  0.478508449425127
               0.000000000000000   1.106412898626718   1.106412898626718  0.032303742334037
               0.000000000000000  -1.106412898626718   1.106412898626718  0.032303742334037
               0.000000000000000   1.106412898626718  -1.106412898626718  0.032303742334037
               0.000000000000000  -1.106412898626718  -1.106412898626718  0.032303742334037
               1.106412898626718   0.000000000000000   1.106412898626718  0.032303742334037
              -1.106412898626718   0.000000000000000   1.106412898626718  0.032303742334037
               1.106412898626718   0.000000000000000  -1.106412898626718  0.032303742334037
              -1.106412898626718   0.000000000000000  -1.106412898626718  0.032303742334037
               1.106412898626718   1.106412898626718   0.000000000000000  0.032303742334037
              -1.106412898626718   1.106412898626718   0.000000000000000  0.032303742334037
               1.106412898626718  -1.106412898626718   0.000000000000000  0.032303742334037
              -1.106412898626718  -1.106412898626718   0.000000000000000  0.032303742334037];

    elseif k <= 9

        xw = [ 0.000000000000000   0.000000000000000   0.000000000000000   0.588405321380412
               1.064082230328777   0.000000000000000   0.000000000000000  -0.152097068487023
              -1.064082230328777   0.000000000000000   0.000000000000000  -0.152097068487023
               0.000000000000000   1.064082230328777   0.000000000000000  -0.152097068487023
               0.000000000000000  -1.064082230328777   0.000000000000000  -0.152097068487023
               0.000000000000000   0.000000000000000   1.064082230328777  -0.152097068487023
               0.000000000000000   0.000000000000000  -1.064082230328777  -0.152097068487023
               0.905830033000216   0.000000000000000   0.000000000000000   0.369012523996709
              -0.905830033000216   0.000000000000000   0.000000000000000   0.369012523996709
               0.000000000000000   0.905830033000216   0.000000000000000   0.369012523996709
               0.000000000000000  -0.905830033000216   0.000000000000000   0.369012523996709
               0.000000000000000   0.000000000000000   0.905830033000216   0.369012523996709
               0.000000000000000   0.000000000000000  -0.905830033000216   0.369012523996709
               0.817286490798906   0.817286490798906   0.817286490798906   0.104007450974435
               0.817286490798906  -0.817286490798906   0.817286490798906   0.104007450974435
               0.817286490798906   0.817286490798906  -0.817286490798906   0.104007450974435
               0.817286490798906  -0.817286490798906  -0.817286490798906   0.104007450974435
              -0.817286490798906   0.817286490798906   0.817286490798906   0.104007450974435
              -0.817286490798906  -0.817286490798906   0.817286490798906   0.104007450974435
              -0.817286490798906   0.817286490798906  -0.817286490798906   0.104007450974435
              -0.817286490798906  -0.817286490798906  -0.817286490798906   0.104007450974435
               0.501292956337400   0.501292956337400   0.501292956337400   0.380660357224238
               0.501292956337400  -0.501292956337400   0.501292956337400   0.380660357224238
               0.501292956337400   0.501292956337400  -0.501292956337400   0.380660357224238
               0.501292956337400  -0.501292956337400  -0.501292956337400   0.380660357224238
              -0.501292956337400   0.501292956337400   0.501292956337400   0.380660357224238
              -0.501292956337400  -0.501292956337400   0.501292956337400   0.380660357224238
              -0.501292956337400   0.501292956337400  -0.501292956337400   0.380660357224238
              -0.501292956337400  -0.501292956337400  -0.501292956337400   0.380660357224238
               1.017168937265364   0.650007853956632   0.000000000000000   0.0930316449988371
               1.017168937265364  -0.650007853956632   0.000000000000000   0.0930316449988371
               0.650007853956632   1.017168937265364   0.000000000000000   0.0930316449988371
              -0.650007853956632   1.017168937265364   0.000000000000000   0.0930316449988371
              -1.017168937265364   0.650007853956632   0.000000000000000   0.0930316449988371
              -1.017168937265364  -0.650007853956632   0.000000000000000   0.0930316449988371
               0.650007853956632  -1.017168937265364   0.000000000000000   0.0930316449988371
              -0.650007853956632  -1.017168937265364   0.000000000000000   0.0930316449988371
               0.000000000000000   0.650007853956632   1.017168937265364   0.0930316449988371
               0.000000000000000  -0.650007853956632   1.017168937265364   0.0930316449988371
               0.000000000000000   1.017168937265364   0.650007853956632   0.0930316449988371
               0.000000000000000   1.017168937265364  -0.650007853956632   0.0930316449988371
               0.000000000000000   0.650007853956632  -1.017168937265364   0.0930316449988371
               0.000000000000000  -0.650007853956632  -1.017168937265364   0.0930316449988371
               0.000000000000000  -1.017168937265364   0.650007853956632   0.0930316449988371
               0.000000000000000  -1.017168937265364  -0.650007853956632   0.0930316449988371
               1.017168937265364   0.000000000000000   0.650007853956632   0.0930316449988371
               1.017168937265364   0.000000000000000  -0.650007853956632   0.0930316449988371
               0.650007853956632   0.000000000000000   1.017168937265364   0.0930316449988371
              -0.650007853956632   0.000000000000000   1.017168937265364   0.0930316449988371
              -1.017168937265364   0.000000000000000   0.650007853956632   0.0930316449988371
              -1.017168937265364   0.000000000000000  -0.650007853956632   0.0930316449988371
               0.650007853956632   0.000000000000000  -1.017168937265364   0.0930316449988371
              -0.650007853956632   0.000000000000000  -1.017168937265364   0.0930316449988371];

    else

        error('Precision not supported!')

    end

    x = xw(:,1:3);
    w = xw(:,4)/8;
    n = numel(w); 

end
