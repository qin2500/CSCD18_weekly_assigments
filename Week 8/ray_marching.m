0%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% CSC D18 - Computer Graphics - UTSC
%
% This script has the goal of letting you explore ray marching, a general
% technique to process light rays that trades flexibility for computational 
% expense.
%
% The case we are interested in is the intersection test for *impicit* 
% surfaces for which it is not straightforward to solve for the
% intersection point algebraically. 
%
% Such surfaces are therefore tricky to render. Ray marching will allow
% us to render these at a level of accuracy that depends on how much
% computation we want to invest into the ray marching process.
%
% function [img]=ray_marching(imSize,winSize,f,e,l,t)
%
%    imSize - size of the output image in pixels (e.g. [512 512])
%    winSize - size of the viewing window (e.g. [-5 5] gives a window 10 units wide in each direction)
%    f - the camera focal length (a positive number)
%    e - camera center of projection in world coordinates (a vector with [x y z])
%    l - a point in world coordinates the camera is 'looking at' (so we can compute g) 
%    t - camera's up vector
%
%    - Returns the image produced
%
% Example:   ray_marching([512 512],[-5 5],2,[3 3 -7],[0 0 0],[0 1 0]);
%   This produces a 512x512 image, with a viewing window between -5 and 5 in each
%   direction (in camera coordinates). The camera has a focal length of 2, is
%   located at [3 3 -7], looking at the world origin ([0 0 0]) and with an up
%   vector of [0 1 0].
%
% (C) F. Estrada, Oct. 2020
%
% You can find high quality renders of implicit surfaces here:
% https://cindyjs.org/gallery/cindygl/Raytracer/index.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [img]=ray_marching(imSize,winSize,f,e,l,t)

dxy=(winSize(2)-winSize(1))/(imSize(1)-1);   % change in coordinates for each pixel

img=zeros([imSize 3]);              % Initially empty image

% Surface properties - change these to change the surface's appearance
col=[0 .5 1];                    % Surface colour as an RGB triplet in [0,1]
phongCoeff=[.1 .9 .75 3];        % Phong model coefficients Ia, Id, Is, shinyness - Assume Id = Is = [1 1 1] (it's just a totally white light);

% Camera parameters - assumed to be looking toward the origin along the Z axis
%e=[0 0 -7];                     % Camera center
%f=2;                            % Focal length

% Set a point light source somewhere
ls=[-2 2 -9];

% Compute camera coordinate frame
g=l-e;
wc=-g/norm(g);         % w vector
uc=cross(t,wc);
uc=uc/norm(uc);        % u vector
vc=cross(wc,uc);       % v vector


% Like all ray tracing, we will have a loop over the image pixels
fprintf(2,'Processing column: ');
for i=1:imSize(1);
 fprintf(2,'%d, ',i);
 for j=1:imSize(2);
  
    % Like your ray tracer, the first step is to figure out the ray equation for a ray through pixel
    % (i,j). Here we go:
    
    p0=e-(f*wc);                                             % start at e, move -f units along w
    p0=p0+(winSize(1)+(dxy*i))*uc;                           % move to the right x location along the u direction
    p0=p0+(winSize(2)-(dxy*j))*vc;                           % and the correct y location along the v direction        
    d=p0-e;                                                 % Ray direction
    d=d/norm(d);                                            % normalized to unit length
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now that we have the ray, the ray marching process consists of the following:
    %
    % TO DO:
    %
    % - Get the value of the implicit equation for the surface we're rendering at p0 - call this v0
    %   LOOP:
    %    - Take small steps (of size dl, delta lambda) in the direction d by incrementing lambda
    %    - Get the position along the ray given the current lambda
    %    - Evaluate the implicit equation at this point, call the resulting value v
    %    - If we're very lucky and it's zero, that's the intersection! evaluate colour and set pixel colour
    %    - If the *SIGN* of v is different from the sign of v0, we must have passed the intersection
    %         * the intersection is somewhere between the previous point and the current one,
    %           find a location, and evaluate colour at this point to set the pixel to the right colour
    %      ELSE (the sign of v0 is the same as the sign of v) - ccontinue the ray marching loop
    %
    %  Evaluating colour means:
    %    - Get the normal at the intersection (how do we do this from an implicit equation?)
    %    - Evaluate the Phong model to get pixel colour
    %    - Set the pixel colour img(i,j,:)=[r g b];  % where [r,g,b] is the colour you computed in [0,1]
    %
    % If you do this correctly, you should be able to render the surface and it should be properly
    % shaded!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
    % First try out rendering a sphere - you will easily know if your code is working
    % Once you have the sphere working try out the other shapes and have fun rendering those...
 
    p=p0+(lambda*d);
    x=p(1);
    y=p(2);
    z=p(3);
    
    % SELECT which function you want to render by using the corresponding implicit form below - make sure to
    % adjust the window size above to the corresponding range!
    v0=(x^2)+(y^2)+(z^2)-4;                                                         % Sphere (window size -1,1)
%    v0=(x^3)+(y*y*z*z)-(z*z);                                                       % Hummingbird (window size -5,5)
%    v0=(y*y)+(z*z)+(x*x*x)-(x*x);                                                   % Ding Dong  (window size -.5,.5)  
%    v0=(((x*x)+((9/4)*z*z)+(y*y)-1)^3)-((x*x)*(y*y*y))-((9/80)*(z*z)*(y*y*y));      % Sweet (window size -.5,.5)
%    v0=x^4-5*x^2+y^4-5*y^2+z^4-5*z^2+11.8;                                          % Tangle cube (window -1.5,1.5)
%    v0=x^4+2*x^2*y^2+2*x^2*z^2+y^4+2*y^2*z^2+z^4+8*x*y*z-10*x^2-10*y^2-10*z^2+25;    % Tetrahedral (window -1.5, -1.5)
%    v0=x^4+y^4+z^4-x^2-y^2-z^2+0.5;                                                 % Chubbs (window -.5, .5)

    dl=.1;                                                  % Increment in lambda - if your code is missing the surface, make this smaller!
    lambda=0;                                               % Initial lambda - start at p0
    done=0;         
    maxDist = 50;                                        % Set this flag when you have completed the ray marching for this ray to exit the loop    
    while (~done)
     % Implement this loop to do the ray marching, 
     % you'll need to get the value of the implicit equation at the current lambda
     % check if you've reached the intersection (or crossed over it!)
     % and when you've found the intersection, compute the normal
     % so you can do shading using Phong.
        
     % Conditions to end the loop:
     % - You reached the intersection
     % - Lambda is greater than some maximum value (in which case we assume there's no intersection). 
     % You get to decide how large that value should be
     
     % Set the corresponding image pixel (i,j) at the colour that results from evaluating the Phong
     % illumination model at the intersection point returned by the ray marching process. Include
     % an ambient term, a diffuse term, and a specular term, the object colour, phong coefficients,
     % and light location and colour are all defined above, so once you get an intersection point
     % and normal you have all the information needed to compute Phong shading at this point.

        lambda = lambda + dl;

        p = p0 + lambda * d;
        x = p(1);
        y = p(2);
        z = p(3);

        v = x^2 + y^2 + z^2 - 4;

        % Check if v is zero
        if abs(v) < 0.001
            done = 1;

            % calculate normals
            N = 2 * [x, y, z];
            
            
            N = normalize(N); % Normal vector
            L = normalize(ls - p); % Light direction
            V = normalize(e - p); % Viewer direction
            R = reflect(-L, N); % Reflection vector

            % Calculate ambient component
            ambient = 0.1; % Ambient light intensity (can be adjusted)
            ambient_color = ambient * col; % Ambient color contribution

            % Calculate diffuse component
            diff = max(dot(N, L), 0); % Lambert's cosine law
            diffuse_color = diff * col; % Diffuse color contribution

            % Calculate specular component
            shininess = phongCoeff.shininess; % Shininess coefficient
            spec = max(dot(R, V), 0) ^ shininess; % Specular intensity
            specular_color = phongCoeff.specular * spec; % Specular color contribution

            % Combine all components
            color = ambient_color + diffuse_color + specular_color;

            % Clamp the color values to the range [0, 1]
            color = min(max(color, 0), 1);

            img(i,j,:) = color;
            

        elseif lambda > maxDist
            done = 1;  % No intersection, leave the pixel black
        end



    end;
    
 end;
end;



% Let's see your image!
fprintf(2,'... Done!\n');
img=img/max(img(:));
imwrite(uint8(255*img),'output.jpg');       % Save this! in case you want it later.
figure(1);clf;image(img);axis image;
end

function [v] = normalize(v)
    % Normalizes a vector
    v = v / norm(v);
end

function [r] = reflect(I, N)
    % Reflects vector I around normal N
    r = I - 2 * dot(I, N) * N;
end