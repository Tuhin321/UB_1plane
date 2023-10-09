function [MX,freq]=fft_hanning(X,NFFT,Fs)
% Function calculates the FFT from X data. Number of points is NFFT and
% sampling frequency is Fs
% Function returns
% MX = vector of amplitudes
% freq = frequency vector


Fn=Fs/2; % Nyquist frequency

%Hanning windowing, removes cutting errors
W=hanning(length(X));
X=W.*X;

% Take fft, padding with zeros, length(FFTX)==NFFT
FFTX=fft(X,NFFT);
NumUniquePts = ceil((NFFT+1)/2);

% fft is symmetric, throw away second half
FFTX=FFTX(1:NumUniquePts);
MX=abs(FFTX);  % Take magnitude of X
% Multiply by 2 to take into account the fact that we threw out second half
% of FFTX above
% Multiplying by 2 for taking into account the impact of Hanning windowing
MX=MX*2*2; 
MX(1)=MX(1)/2;   % Account for endpoint uniqueness
if ~rem(NFFT,2),
    MX(length(MX))=MX(length(MX))/2;
end
MX(length(MX))=MX(length(MX))/2;  % We know NFFT is even

% Scale the FFT so that it is not a function of the length of X.
MX=MX/length(X); % Amplitudit
freq=(0:NumUniquePts-1)*2*Fn/NFFT; %taajuusvektori
% FFT loppuu---------------------------------------------------------------


% Local functions (in case if not included in Matlab installation
% seems to require Signal Processing ToolBoxi

function w = hanning(varargin)
%HANNING   Hanning window.
%   HANNING(N) returns the N-point symmetric Hanning window in a column
%   vector.  Note that the first and last zero-weighted window samples
%   are not included.
%
%   HANNING(N,'symmetric') returns the same result as HANNING(N).
%
%   HANNING(N,'periodic') returns the N-point periodic Hanning window,
%   and includes the first zero-weighted window sample.
%
%   NOTE: Use the HANN function to get a Hanning window which has the 
%          first and last zero-weighted samples. 
%
%   See also BARTLETT, BLACKMAN, BOXCAR, CHEBWIN, HAMMING, HANN, KAISER
%   and TRIANG.

%   Copyright 1988-2000 The MathWorks, Inc.
%   $Revision: 1.9 $  $Date: 2000/06/09 22:04:47 $

% Check number of inputs
error(nargchk(1,2,nargin));

% Check for trivial order
[n,w,trivialwin] = check_order(varargin{1});
if trivialwin, return, end

% Select the sampling option
if nargin == 1,
   sflag = 'symmetric';
else
   sflag = lower(varargin{2});
end

% Allow partial strings for sampling options
allsflags = {'symmetric','periodic'};
sflagindex = strmatch(sflag, allsflags);
if length(sflagindex)~=1         % catch 0 or 2 matches
   error('Sampling flag must be either ''symmetric'' or ''periodic''.');
end
sflag = allsflags{sflagindex};

% Evaluate the window
switch sflag,
case 'periodic'
   % Includes the first zero sample
   w = [0; sym_hanning(n-1)];
case 'symmetric'
   % Does not include the first and last zero sample
   w = sym_hanning(n);
end

%---------------------------------------------------------------------
function w = sym_hanning(n)
%SYM_HANNING   Symmetric Hanning window. 
%   SYM_HANNING Returns an exactly symmetric N point window by evaluating
%   the first half and then flipping the same samples over the other half.

if ~rem(n,2)
   % Even length window
   half = n/2;
   w = calc_hanning(half,n);
   w = [w; w(end:-1:1)];
else
   % Odd length window
   half = (n+1)/2;
   w = calc_hanning(half,n);
   w = [w; w(end-1:-1:1)];
end

%---------------------------------------------------------------------
function w = calc_hanning(m,n)
%CALC_HANNING   Calculates Hanning window samples.
%   CALC_HANNING Calculates and returns the first M points of an N point
%   Hanning window.

w = .5*(1 - cos(2*pi*(1:m)'/(n+1))); 

% [EOF] hanning.m


function [n_out, w, trivalwin] = check_order(n_in)
%CHECK_ORDER Checks the order passed to the window functions.
% [N,W,TRIVALWIN] = CHECK_ORDER(N_ESTIMATE) will round N_ESTIMATE to the
% nearest integer if it is not alreay an integer. In special cases (N is [],
% 0, or 1), TRIVALWIN will be set to flag that W has been modified.

%   Copyright 1988-2000 The MathWorks, Inc.
%   $Revision: 1.4 $  $Date: 2000/06/09 20:50:37 $

w = [];
trivalwin = 0;

% Special case of negative orders:
if n_in < 0,
   error('Order cannot be less than zero.');
end

% Check if order is already an integer or empty
% If not, round to nearest integer.
if isempty(n_in) | n_in == floor(n_in),
   n_out = n_in;
else
   n_out = round(n_in);
   warning('Rounding order to nearest integer.');
end

% Special cases:
if isempty(n_out) | n_out == 0,
   w = zeros(0,1);               % Empty matrix: 0-by-1
   trivalwin = 1; 
elseif n_out == 1,
   w = 1;
   trivalwin = 1;   
end

