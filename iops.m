classdef iops < hgsetget
   %% Create an IOPS class
   %
   % ip = iops(telescope, wavelength); returns a default iops object
   % ip = iops(telescope, wavelength,'useJacket',true,'mCPU',true); returns
   % a default iops object for computation over GPU and multiple CPU
   
   properties
       % refractive index contrast
       deltaN;
       % refractive index
       n;
       % resolution in pixels
       resolution;
       % absorbing boundary, logical
       absorbingBoundary;
       % x/y scale
       xs;
       % z scale
       zs;
       % device length in um
       sLength;
       % segmented telescope
       tel;
       % tag
       tag = 'IOPS';
       % use GPU?
       useGPU = false;
       % multiple cpu useage, logical
       mCPU = false;
       % use the jacket toolbox?
       useJacket = false;
       
       %output coeficient paris normalised to the first
       coefs; % coefs - calibrated and normalised coefficients of output
       inputCoefs; % inputCoefs - how much light is coupled into each fibre
       outputCoefs; % outputCoefs - how much light is at the output of each fibre
       referenceCoefs; % referenceCoefs - normalised coefficients from a flat (0 segment piston) input
       
       % input and outputs
       inputE;
       outputE;
       unmaskedInput; % input before the mask
       
       % segment pair vector
       segPair;
       
       % masks for each waveguide input and output
       mask1;
       mask2;
       mask; % total mask
       maskD; %mask diameter
   end
   
   properties (Dependent, SetObservable=true)
       % wavelength in um
       %wavelength;
       
   end
   
   properties (Access=private)
       % wavelength in um
       wavelength;
       
       p_log;
       % lens for imaging
       p_lens;
       
       
       
       % waveguide core size, um
       core;
       % waveguide spacing, um
       spacing;
       
       % waveguide output
       waveguideOutput;
       
       % refractive index profile
       index;
       
       % bpm constants
       F;
       D;
       MD;
       kx;
       ky;
       
       % initialised already?
       init = false;
       
       % telescope segments
       segments;
       
       % x scale
       p_xs;
   end
   
   methods
   
       % Constructor
       function obj = iops(tel,wavelength,varargin)
           p = inputParser;
           p.addParamValue('n', 1.46, @isnumeric);
           p.addParamValue('deltaN',  0.02678, @isnumeric); %for 1550 nm  0.02678, for 1310 nm 0.025875
           p.addParamValue('sLength', 1e4, @isnumeric);
           p.addRequired('tel', @(x) isa(x,'telescopeAbstract'));
           p.addRequired('wavelength', @(x) isnumeric(x) || isa(x,'photometry') );
           p.addParamValue('resolution', 128, @isnumeric);
           p.addParamValue('absorbingBoundary', false, @islogical);
           p.addParamValue('xs', 0.2, @isnumeric);
           p.addParamValue('zs', 0.2, @isnumeric);
           p.addParamValue('useGPU', false, @islogical);
           p.addParamValue('mCPU', false, @islogical);
           p.addParamValue('useJacket', false, @islogical);
           p.addParamValue('segPair', [1:6;2:7], @isnumeric); % a vector of inputs to compare

           p.addParamValue('logging', true, @islogical);

           % parse input args
           p.parse(tel, wavelength, varargin{:});
           if isa(p.Results.wavelength,'photometry')
               obj.wavelength = p.Results.wavelength.wavelength;
           else
               obj.wavelength = p.Results.wavelength;
           end

           obj.deltaN = p.Results.deltaN;
           obj.n = p.Results.n;
           
           
           if isempty(p.Results.resolution)
               obj.resolution = tel.resolution;
           else
               obj.resolution = p.Results.resolution;
           end
           
           obj.absorbingBoundary = p.Results.absorbingBoundary;
           obj.xs = p.Results.xs;
           obj.zs = p.Results.zs;
           
           obj.sLength = p.Results.sLength;
           obj.tel = p.Results.tel;
           obj.segments = obj.tel.segment;

           if p.Results.logging
               obj.p_log = logBook.checkIn(obj);
               display(obj);
           end
           
           obj.useGPU = p.Results.useGPU;
           obj.mCPU = p.Results.mCPU;
           obj.useJacket = p.Results.useJacket;
           
           obj.p_lens = lens(1);
           obj.p_lens.nyquistSampling = 2; %must match pixel scale of BPM!!
           obj.p_lens.fieldStopSize = 64;
           
           obj.segPair = p.Results.segPair;
           
           % I like even numbers!
           if rem(obj.tel.resolution,2) % odd number of pixels in telescope
               fprintf('IOPS WARNING: Odd number of pixels in Telescope object, IOPS may not run');
           end
           
           obj.core = 4;
           obj.spacing = 8;
           obj.segments = tel.segment();
           
           obj.initIOPS();
           obj.init = true;

       end

       % Deconstructor
       function delete(obj)
           if ~isempty(obj.p_log)
               checkOut(obj.p_log,obj)
           end
       end

       function initIOPS(obj)
           %% initIOPS
           % initialise everything for iops
           obj.calibrateXS();
           obj.createMask();
           obj.initBPM();
           obj.resetCoefs();
           obj.resetReferenceCoefs();
           obj.resetImages();
           obj.calibrateReferenceCoefs();
       end
       
       function display(obj)
           %% DISPLAY Display object information
           %
           % disp(obj) prints information about the iops object
           
           fprintf('___ %s ___\n',obj.tag);
           fprintf(' %d resolution with x and z scale %d and %d over %d um\n',obj.resolution,obj.p_xs,obj.zs,obj.sLength);
           fprintf('----------------------------------------------------\n')
       end
       
       function val = get.MD(obj)
           val = obj.MD;
       end
     
       
%        function varargout = uplus(obj)
%            if nargout>0
%                varargout{1} = obj;
%            end
%        end
       
       function relay(obj,src)
           %% RELAY iops to source and sets coeficients based on output
           %
           % relay(obj,src) propagates an image of the source from
           % each segment through iops and sets coefs to the output value 
           % pairs.
           
           % check if matlabpool is already open, and if we want to use
           % multiple cpu
           poolWasAlreadyOpen = true;
           if matlabpool('size') == 0 && obj.mCPU
               try
                   matlabpool('open')
                   poolWasAlreadyOpen = false;
               catch MErr
                   fprintf('Error starting matlabpool in iops %s',MErr);
               end
           else
               if ~obj.mCPU
                   poolWasAlreadyOpen = false;
               end
           end
           
           % create local version of everything - useful in parfor
           localCoefs = zeros(2,length(obj.coefs));
           localInputCoefs = zeros(2,length(obj.inputCoefs));
           
           res = obj.resolution;
           
           % input and output complex amplitudes
           gE = zeros(res,res,length(obj.segPair));
           iE = zeros(res,res,length(obj.segPair));
           uE = zeros(res,res,length(obj.segPair));
           
           % _try_ to supress some oomao output!
           obj.p_log.verbose = false;
           
           % loop through each segment and compare in pairs
           parfor ii = 1:length(obj.segPair)
               
               % create the input
               input = generateInput(obj,obj.segPair(1,ii),obj.segPair(2,ii),src);
               uE(:,:,ii) = input; %unmasked input
               input = input.*obj.mask; %masked input
               iE(:,:,ii) = input; %store for later

                
               gE(:,:,ii) = bpm(obj,input); %propagate input and keep output
               

               localCoefs(:,ii) = obj.measureOutput(abs(gE(:,:,ii)).^2); %measure output integral and keep
               localInputCoefs(:,ii) = obj.measureOutput(abs(iE(:,:,ii)).^2);% measure input integral and keep
           end
           
           if obj.mCPU && ~poolWasAlreadyOpen
               matlabpool('close');
           end
           
           obj.inputE = iE;
           obj.outputE = gE;
           obj.unmaskedInput = uE;
           obj.outputCoefs = localCoefs;
           obj.inputCoefs = localInputCoefs;
           
           obj.coefs = (localCoefs./localInputCoefs)./obj.referenceCoefs; %normalise coefficients to input and reference
           
           obj.p_log.verbose = true;
           
           
       end
       
       function input = generateInput(obj,seg1,seg2,src)
           %% Generates the input image with phase
           
           % keep a few pixels outside mask for rounding error
           padding = 2; % take an extra few pixels
           maskR = round(obj.maskD/2); % input diameter in px
           offset = round(obj.spacing/obj.p_xs/2); % input spacing in px
           
           
           res = obj.resolution;
           %tel pixel scale
           gmtPx = obj.tel.resolution/obj.tel.D; %GMT pixel scale, pix/m
           
           telRes = obj.tel.resolution;
           % padding for fft - image size matches physical object size
           padRes = round(telRes*(4*0.2043/obj.xs)+1); % nyquist sampling *2
           % create circular pupils
           [X Y] = meshgrid(-padRes/2+0.5:padRes/2-0.5,-padRes/2+0.5:padRes/2-0.5);
           % segment 1 is different to the rest
           
           % subtract one pixel to try and match tel pupil better
           segPupil1 = (hypot(X,Y) < (obj.tel.segmentD*gmtPx/2 - 1)) - (hypot(X,Y) < (obj.tel.centralObscurationD*gmtPx/2 + 1));
           
           segPupil = hypot(X,Y) < (obj.tel.segmentD*gmtPx/2 - 1);
           % noise in the complex part of the image - an fft artifact only
           seg1FftNoise = angle(fftshift(fft2(segPupil1)));
           segFftNoise = angle(fftshift(fft2(segPupil)));
           %radius of segment pupil in px
           segPupilR = round(gmtPx*obj.tel.segmentD/2); 
           
           segCoords = [round(gmtPx*real(obj.tel.segmentCoordinate))+telRes/2; round(gmtPx*imag(obj.tel.segmentCoordinate))+telRes/2];
           % create an image from the segment

           % IMAGE CREATION WITH CENTRED SEGMENT PUPIL

           % some variables to make different pupil geometry easier later
           % on
           pup1 = segPupil;
           pup2 = segPupil;
           noise1 = segFftNoise;
           noise2 = segFftNoise;

           % use central segment pupil
           if seg1 == 1
               pup1 = segPupil1;
               noise1 = seg1FftNoise;
           end

           % use central segment pupil if comparing seg 7 to 1
           if seg2 == 1
               pup2 = segPupil1;
               noise2 = seg1FftNoise;
           end


           % create two images, one for each segment to be compared
           % image for seg1
           phase = zeros(padRes);
           pr2 = round(padRes/2);
           % take the phase from src.phase
           phase(pr2-segPupilR+1:pr2+segPupilR,...
                pr2-segPupilR+1:pr2+segPupilR) =...
                    src.phase(segCoords(2,seg1)-(segPupilR-1):segCoords(2,seg1)+segPupilR,...
                    segCoords(1,seg1)-(segPupilR-1):segCoords(1,seg1)+segPupilR);
                
           phase = phase.*segPupil; %remove any clipping of other segments - done use segPupil1 here!

           % create a complex amplitude from this phase and pupil
           wave = pup1.*exp(1i.*phase);
           
           % create an image from this complex amplitude
           input1 = fftshift(fft2(wave));
           % recreate the image, this time removing noise generated by fft
           input1 = abs(input1).*exp(1i.*(angle(input1)-noise1));

           % same thing for the other segment
           phase2 = zeros(padRes);
           phase2(pr2-segPupilR+1:pr2+segPupilR,...
                pr2-segPupilR+1:pr2+segPupilR) =...
                    src.phase(segCoords(2,seg2)-(segPupilR-1):segCoords(2,seg2)+segPupilR,...
                    segCoords(1,seg2)-(segPupilR-1):segCoords(1,seg2)+segPupilR);
                
           phase2 = phase2.*segPupil;

           wave2 = pup2.*exp(1i.*phase2);

           input2 = fftshift(fft2(wave2));

           input2 = abs(input2).*exp(1i.*(angle(input2)-noise2));
           
           %normalise inputs: take reduced real part (due to exp(i phase shift) into account
           
           %normilisation constants
%            norm1 = max(abs(input1(:)))/real(max(wave2(:)));
%            norm2 = max(abs(input1(:)))/real(max(wave(:)));
%            
%            
%            % sometimes the segment will be slightly off (check
%            % imagesc(real(wave)) - you will see some peaks around the edge
%            % with a value of 1 - you will only see this when real(wave) <1.
%            % to fix this check if the real part over the pupil is different
%            % to the mean over the pupil, of it is there is some error and
%            % the normalisation needs to be different for the other segment
%            if abs(min(real(wave(logical(pup1)))) - mean(real(wave(logical(pup1))))) > 1e-12
%                fprintf('IOPS WARNING: One input (1) unexpectedly large, trying to compensate\n');
%                norm2 = max(abs(input1(:)))/min(real(wave(logical(pup1))));
%            end
%            
%            if seg2 ~=1 && abs(min(real(wave2(logical(pup2)))) - mean(real(wave2(logical(pup2))))) > 1e-12
%                fprintf('IOPS WARNING: One input (2) unexpectedly large, trying to compensate\n');
%                norm1 = max(abs(input1(:)))/min(real(wave2(logical(pup2))));
%            end
%            
%            % normalise inputs
%            input1 = input1./abs(norm1);
%            input2 = input2./abs(norm2);
           norm1 = max(abs(input1(:)));
           input1 = input1./norm1;
           input2 = input2./norm1;
           
           % create a single input by taking the centre of input1 and
           % input2 and placing them in the correct position for coupling
           % into the fibres
           [lx ly] = size(input1);

           % x and y center of image
           cx = round(lx/2);
           cy = round(ly/2);

           % the mask for each input is obj.core/obj.xs pixels wide -
           % need this many pixels from each input + padding


           % clip each input
           input1 = input1(cx - maskR - padding:cx + maskR + padding, ...
                            cy - maskR - padding:cy + maskR + padding);
           input2 = input2(cx - maskR - padding:cx + maskR + padding, ...
                            cy - maskR - padding:cy + maskR + padding);

           input = zeros(res,res);
           % place each input in the correct place
           input(res/2 - maskR - padding:res/2 + maskR + padding, ...
                res/2 - offset - maskR - padding: res/2 - offset + maskR + padding) = input1;


           input(res/2 - maskR - padding:res/2 + maskR + padding, ...
                res/2 + offset - maskR - padding-1: res/2 + offset + maskR + padding-1) = input2;
            
           % both inputs are now normalised to eachother, now normalise the
           % whole thing
           input = input./max(abs(input(:)));


%            input = input.*obj.mask;
       end
       
       function testBPM(obj,offset)
           %% Test the BPM with current parameters
           %


           % test bpm with gaussian inputs offset
           fprintf('Testing BPM with gaussian inputs offset by %d rad\n',offset);

           inputOffset = offset;
           
           [X Y] = meshgrid(-obj.resolution/2+1:obj.resolution/2);
           xw = 2;% gaussian width
           
           input = exp(-((X+obj.spacing/2/obj.p_xs).*obj.p_xs/xw).^2 - ((Y).*obj.p_xs/xw).^2) + exp(-((X-obj.spacing/2/obj.p_xs).*obj.p_xs/xw).^2 - ((Y).*obj.p_xs/xw).^2+1i*inputOffset);
           input = input.*obj.mask;
           
           out = bpm(obj,input);
           
           figure();
           imagesc(abs(out).^2);
           
           coef = measureOutput(obj,abs(out).^2);
           
           fprintf('Gaussian coef: %d %d\n Gaussian ratio: %d\n',coef(1),coef(2),coef(1)/coef(2));
           
           % test with images created from obj.tel segments
           
           src = source();
           src = src.*obj.segments{1}*obj.p_lens;
           width = iops.gaussWidth(abs(src.wave).^2);
           
           input = exp(-((X+obj.spacing/2/obj.p_xs)./width).^2 - ((Y)./width).^2);
           input = input.*obj.mask;
           
           out = bpm(obj,abs(input).^2);
           
           figure();
           imagesc(abs(out).^2);
           
           coef = measureOutput(obj,out);
           
           fprintf('Gaussian coef: %d %d\n Gaussian ratio: %d\n',coef(1),coef(2),coef(1)/coef(2));
           
       end
       
       %% Set/Get xs
       function val = get.xs(obj)
           val = obj.p_xs;
       end
       
       function set.xs(obj,val)
           obj.xs=val;
           obj.p_xs=val;
           if obj.init
               obj.initIOPS();
           end
       end
       
       %% Set/Get zs
       function val = get.zs(obj)
           val = obj.zs;
       end
       
       function set.zs(obj,val)
           obj.zs = val;
           if obj.init
               obj.initBPM();
               obj.calibrateReferenceCoefs();
           end
       end
       
       %% Set/Get n
       function val = get.n(obj)
           val = obj.n;
       end
       
       function set.n(obj,val)
           obj.n = val;
           if obj.init
               obj.initBPM();
               obj.calibrateReferenceCoefs();
           end
       end
       
       %% Set/Get deltaN
       function val = get.deltaN(obj)
           val = obj.deltaN;
       end
       
       function set.deltaN(obj,val)
           obj.deltaN = val;
           if obj.init
               obj.initBPM();
               obj.calibrateReferenceCoefs();
           end
       end
       
       %% Set/Get sLength
       function val = get.sLength(obj)
           val = obj.sLength;
       end
       
       function set.sLength(obj,val)
           obj.sLength = val;
           if obj.init
               obj.calibrateReferenceCoefs();
           end
       end
       
       %% Set/Get tel
       function val = get.tel(obj)
           val = obj.tel;
       end
       
       function set.tel(obj,val)
           obj.tel = val;
           
           obj.segments = obj.tel.segment();
       end
       
       %% Set/Get resolution
       function val = get.resolution(obj)
           val = obj.resolution;
       end
       
       function set.resolution(obj,val)
           obj.resolution = val;
           if obj.init
               obj.initIOPS();
           end
       end
       
       %% Set/Get segPair
       function val = get.segPair(obj)
           val = obj.setPair;
       end
       
       function set.segPair(obj,val)
           obj.setPair = val;
           
           obj.resetCoefs();
           obj.resetReferenceCoefs();
           obj.calibrateReferenceCoefs();
       end
       
       function output = measureOutput(obj,buf)
           %% measureOutput
           % sum each output image around its mask and return the pair of results
           output = [sum(sum(buf.*obj.mask1)); sum(sum(buf.*obj.mask2))];
       end
       
       
       function calibrateReferenceCoefs(obj)
           %% calibrateReferenceCoefs
           % generate reference coefs from flat input
           
           obj.resetCoefs();
           
           src = source('wavelength',photometry.H);
           src=src.*obj.tel.pupil;
           
           relay(obj,src);%ceofs set in here
           
           obj.referenceCoefs = obj.outputCoefs./obj.inputCoefs; %normalise reference coefficients to input
           % reset other coefs
           obj.resetCoefs();
       end
       
       function resetCoefs(obj)
           %% resetCoefs
           % resets obj coefs - not including reference coefs
           obj.coefs = zeros(2,length(obj.segPair));
           obj.inputCoefs = zeros(2,length(obj.segPair));
           obj.outputCoefs = zeros(2,length(obj.segPair));
       end
       
       function resetReferenceCoefs(obj)
           %% resetReferenceCoefs
           % resets reference coefs
           obj.referenceCoefs = zeros(2,length(obj.segPair));
       end
       
       function calibrateXS(obj)
           %% calibrateXS
           % set the x scale (obj.xs) using an image from a telescope
           % segment and a gaussian of known width
           
           % padding for fft - image size matches physical object size
           res = obj.tel.resolution*(4*0.2/obj.xs)+1; % nyquist sampling *2
           
           % pixel scale of tel
           gmtPx = obj.tel.resolution/obj.tel.D; %GMT pixel scale, pix/m
           
           %create a circular pupil in the centre
           [X Y] = meshgrid(-res/2:res/2-1,-res/2:res/2-1);
           
           segPupil1 = (hypot(X,Y) < obj.tel.segmentD*gmtPx/2) - (hypot(X,Y) < obj.tel.centralObscurationD*gmtPx/2);
           
           % create an image from this pupil
           img = fftshift(fft2(segPupil1));
           
           % measure the width of this image
           width = iops.gaussWidth(abs(img).^2);
           
           
           % create a Gaussian of known width
           [X Y] = meshgrid(-obj.resolution/2:obj.resolution/2-1);
           inputWidth = 2;
           inE = exp(-((X).*obj.xs/inputWidth).^2 - ((Y).*obj.xs/inputWidth).^2);
           % measure the width of this Gaussian - in pixels
           inWidth = iops.gaussWidth(abs(inE).^2);
           
           % match the pixel scale to the ratio of these two widths: this
           % matches the image size of correct scale to the fibre size
           obj.p_xs = obj.xs/(width/inWidth); % change the input scale depending on fwhm of image
           
           obj.xs = obj.p_xs;
           
           if obj.p_xs > 0.4
               fprintf('IOPS WARNING: Simulation grid may be too coarse, consider decreasing xs for more accurate results\n');
           end
       end
       
       function createMask(obj)
           %% createMask
           % create and set the fibre profile and mask for input

           
           obj.maskD = round(obj.core/obj.xs);
  
%            disk = fspecial('disk',obj.maskD/2);
%            disk = disk./max(max(disk));
           
           
           
           % want the mask to have an odd number of pixels!!
%            if ~rem(obj.maskD,2)
%                obj.maskD = obj.maskD + 1; 
%            end
           
           % fibre profile made of two discs of diameter core separated by spacing
           obj.maskD = 2*round(obj.maskD/2) + 1; %use an odd number of pixels for the mask
           disk=tools.piston(obj.maskD);
           
           obj.mask = zeros(obj.resolution);
           
           [sx sy] = size(disk);
           
           
           % place the disk fibre profile on each side of centre.
           % not using tools.piston with an offset here in case another
           % profile is needed (eg fspecial('disk',obj.maskD/2); )
           obj.mask(round(obj.resolution/2-sx/2:obj.resolution/2+sx/2-1),...
                    round(obj.resolution/2-obj.spacing/2/obj.xs-sy/2:obj.resolution/2-obj.spacing/2/obj.xs+sy/2-1)) = disk;
           % mask1 and mask2 are the mask for a single fibre
           obj.mask1 = obj.mask;
                
           obj.mask(round(obj.resolution/2-sx/2:obj.resolution/2+sx/2-1),...
                    round(obj.resolution/2+obj.spacing/2/obj.xs-sy/2:obj.resolution/2+obj.spacing/2/obj.xs+sy/2-1)) = disk;
                
           obj.mask2 = obj.mask - obj.mask1;
           
%            obj.mask1 = tools.piston(obj.core/obj.xs,obj.resolution,-obj.spacing/2/obj.xs,0,'type','double','shape','disc');
%            obj.mask2 = tools.piston(obj.core/obj.xs,obj.resolution,obj.spacing/2/obj.xs,0,'type','double','shape','disc');
%            obj.mask = obj.mask1+obj.mask2;
       end
       
       function initBPM(obj)
           %% initBPM
           % parameters for bpm
           obj.index = zeros(obj.resolution,obj.resolution)+obj.mask.*obj.deltaN;
           obj.index = obj.index+obj.n;
           
           
           obj.F = (2*pi/obj.wavelength).*obj.index;% index profile
           obj.D = obj.wavelength/(4*pi*obj.n);% diffraction coeff
           obj.kx = [0:obj.resolution/2-1,(-obj.resolution/2):-1].*(2*pi/(obj.xs*obj.resolution));% Fourier space index - no fftshift required
           obj.ky = [0:obj.resolution/2-1,(-obj.resolution/2):-1].*(2*pi/(obj.xs*obj.resolution));
           % Fourier space diffraction - for BPM
           obj.MD = zeros(obj.resolution,obj.resolution);
           for kyi=1:obj.resolution
                for kxi=1:obj.resolution
                   obj.MD(kxi,kyi) = exp(-1i*obj.D*((obj.kx(kxi)^2)+(obj.ky(kyi)^2))*obj.zs);
                 end
           end
       end
       
       function resetImages(obj)
           % complex amplitude of input and output
           obj.inputE = zeros(obj.resolution,obj.resolution,length(obj.segPair));
           obj.outputE = zeros(obj.resolution,obj.resolution,length(obj.segPair));
           obj.unmaskedInput = zeros(obj.resolution,obj.resolution,length(obj.segPair)); % not multiplied by mask
       end
       
   end
   
   methods (Static)
        function rc = gaussWidth(gauss)
           %% Find the 1/e width of an input 2D gaussian
           % adapted from tools.fitFwhm for 1/e width
           
           C = contourc(gauss./max(gauss(:)),[1/exp(1) 1/exp(1)]);
           rr=hypot(C(1,2:end),C(2,2:end));
           xc = sum(rr.*C(1,2:end))./sum(rr);
           yc = sum(rr.*C(2,2:end))./sum(rr);
           rc = mean(sqrt((C(1,2:end)-xc).^2+ (C(2,2:end)-yc).^2));
        end
       
        function out = generateImages(imgs)
            %% Generate a single image of all input/output images imgs
            [sx, sy, sz] = size(imgs);
            width = round(sx/3/2);
            out = imgs(sx/2-width:sx/2+width-1,:,1);
            
            for ii = 2:sz
                out = [out; imgs(sx/2-width:sx/2+width-1,:,ii)];
            end
            
        end
        
        
   end
   
   methods (Access = protected)
       function out = bpm(obj,input)
           %%  BEAM PROPAGATION MODEL using fourier solution to Schrodinger equation in 2D
           %
           %
           

            if obj.useGPU
                gMD = gpuArray(obj.MD);
                gE = gpuArray(input);
                gZs = gpuArray(obj.zs);
                gF = gpuArray(obj.F);
            else
                if obj.useJacket
                    gMD = gdouble(obj.MD);
                    gE = gdouble(input);
                    gZs = gdouble(obj.zs);
                    gF = gdouble(obj.F);
                else
                    gMD = obj.MD;
                    gE = input;
                    gZs = obj.zs;
                    gF = obj.F;
                end
            end
            
            
            fprintf('IOPS: Propagating segment pair');
            
            
            for j=0:obj.zs:obj.sLength
                
                if rem(j/obj.zs,obj.sLength/(10*obj.zs)) == 0
                    fprintf('.');
                end
                
                e=fft2(gE);
                newe = e.*gMD;
                newE = ifft2(newe);


                gE = newE.*exp(1i*gF*gZs);

                %absorbing boundary
                if (obj.absorbingBoundary)
                    gE(:,1:3) = 0;
                    gE(:,obj.resolution-2:obj.resolution) = 0;
                    gE(1:3,:) = 0;
                    gE(obj.resolution-2:obj.resolution,:) = 0;
                end

            end
            
            fprintf('DONE\n');
           
            if obj.useGPU
                out = gather(gE);
            else
                if obj.useJacket
                    out = double(gE);
                else
                    out = gE;
                end
            end
           
       end
       
       
       
   end
   
end