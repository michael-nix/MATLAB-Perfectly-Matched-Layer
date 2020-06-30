% setupPML   Helper function to create a Perfectly Matched Layer for 1D or
%   2D domains that pretty much just works when you need absorbing
%   boundary conditions for a standard wave equation, (du/dt).^2 = del2(u).
%
%   sigmax = SETUPPML(x, dx) if x is one-dimensional returns a row or
%   column vector with a symmetric perfectly matched layer.
%
%   sigmax = SETUPPML(x, dx) if x is two-dimensional assumes x was created
%   using meshgrid such that the dimension traverses the columns.
%
%   [sigmax, sigmay] = SETUPPML(x, dx) if x is two-dimensional also returns
%   an identical perfectly matched layer, sigmay, that is the transpose of
%   sigmax.  Assumes you have an identical spatial resolution for both
%   dimensions, rows and columns.
%
%   If you need two different PMLs, one for x and one for y, just call this
%   function twice with a single output, transposing as needed.
%
%   NOTE: this assumes your x-dimension goes along the columns, i.e. 
%   consists of a single row vector, or for 2D, stacked identical row 
%   vectors.
function [sigmax, sigmay] = setupPML(x, dx)

    xiscol = iscolumn(x);
    if xiscol
        x = x.';
    end

    PMLwidth = 15;   m = 4;
    sigmax_max = log(1e4)*(m+1) / (PMLwidth*dx)^(m+1);

    sigmax = sigmax_max*((1:PMLwidth)*dx).^m;
    sigmax = [sigmax(end:-1:1), ...
        zeros(1, size(x, 2) - 2 * PMLwidth), sigmax];

    sigmax = sparse(sigmax);
    sigmax = repmat(sigmax, size(x, 1), 1);

    if xiscol
        sigmax = sigmax.';
    end

    if nargout > 1
        sigmay = sigmax.';
    end
end