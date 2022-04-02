moment = 'Zernike';
metric = 'cityblock';
ord = 5;
rep = 5;
          
content = '*';
extension = 'tiff';            
word = strcat(content,extension);

%path = uigetdir;
path = 'D:\ImmaginiLavoro\Medical\MRI\Emphysema-CT';
images = getFileNames(path);
images = sortImages(images);

folderFeatures = [];

for k = 1:length(images)
    image = imread(images{k});
    switch moment
        case 'HU'                  
            folderFeatures(k, :) = Hu_Moments(image); %salvo in un vettore tutti i valori delle feautures                    
        case 'Zernike'
            folderFeatures(k, :) = Zernike_Moments(image, ord, rep); %salvo in un vettore tutti i valori delle feautures                    
        case 'Legendre'
            folderFeatures(k, :) = legmoms_vec(image, ord); %salvo in un vettore tutti i valori delle feautures                    
        case 'Chebyshev'
            folderFeatures(k, :) = dchebmoms_vec(image, ord); %salvo in un vettore tutti i valori delle feautures                    
    end              
end


similarity = 0;

rDB = 0;
gDB = 0;
ARP = 0;
ARR = 0;
eta = 100;

labels = zeros(length(images),1);
for k = 1:length(images)
    classeCartella = extractBetween(images{k}, "_", ".");
    classeCartella = str2num(classeCartella{1:1});
    labels(k) = classeCartella;           
end

     
for i = 1:length(images)
    classeImmagine = labels(i);
    imageFeatures = folderFeatures(i,:);
    classeSearch = labels;
    classeSearch(i,:)=[];
    searchFeatures = folderFeatures;
    searchFeatures(i,:)=[];

    d = pdist2(searchFeatures, imageFeatures, metric);
    [~, index]=sortrows(d);

    classeSearch = classeSearch(index);
    rDB = sum(classeSearch(1:eta)==classeImmagine);   
    gDB = sum(classeSearch==classeImmagine); 
    
    %salvo il valore della sommatoria
    ARP = ARP + rDB/eta;
    ARR = ARR + rDB/gDB;
end         
%calcolo il valore finale
ARP = (100/length(immagini)) * ARP;  
fprintf('\n\nIl valore di ARP è: %4.1f %', ARP);
% ---------------------------------------------------------        

ARR = (100/length(immagini)) * ARR;            
fprintf('\n\nIl valore di ARR è: %4.1f %', ARR);

% ------------ Formula per il calcolo di FScore -----------         
FScore = (2 * ARP * ARR)/(ARP + ARR);
fprintf('\n\nIl valore di FScore è: %4.1f %', FScore);
fprintf('\n');
% ---------------------------------------------------------  

nomiMisure(1,:) = "ARP";
nomiMisure(2,:) = "ARR";
nomiMisure(3,:) = "FScore";

performance(1,:) = ARP;
performance(2,:) = ARR;
performance(3,:) = FScore;         


function fileNames = getFileNames(path, type)
    % Returns all the names of the files contained in a directory
    %
    % INPUT:
    %       path           : the directory to explore
    %       type           : class of file IMG, MAT, SCRITP, ALL(default)
    %       
    % OUTPUT:
    %       fileNames      : names (full paths) of the files in the directory

    if nargin==0
        %Select the directory browsing for the directory.
        path = uigetdir(pwd, 'Select a folder');
        type = 'ALL';
    end

    if nargin==1
        type = 'ALL';
    end

    %Read the list of files in the directory
    dirContents = dir(path);
    nFiles = size(dirContents,1);
    fileNames = cell(nFiles,1);
    idx = false(nFiles,1);
    for i = 1:nFiles
        record = fullfile(path, dirContents(i).name);
        if isdir(record) ~= 1
           switch type
               case 'ALL'
                   idx(i) = 1;
                   fileNames{i} = record;
               case {'IMG','img','IMAGE','image'}
                   [~,~,ext] = fileparts(record);
                   if strcmp(ext, '.png')==1 || strcmp(ext, '.tif')==1 || strcmp(ext, '.bmp')==1 || strcmp(ext, '.jpg')==1 || strcmp(ext, '.jpeg')==1
                       idx(i) = 1;
                       fileNames{i} = record;
                   end
               case {'MAT','mat'}
                   [~,~,ext] = fileparts(record);
                   if strcmp(ext, '.mat')==1
                       idx(i) = 1;
                       fileNames{i} = record;
                   end
              case {'XML','xml'}
                   [~,~,ext] = fileparts(record);
                   if strcmp(ext, '.xml')==1
                       idx(i) = 1;
                       fileNames{i} = record;
                   end
               case {'SCRIPT','M,','m'}
                   [~,~,ext] = fileparts(record);
                   if strcmp(ext, '.m')==1
                       idx(i) = 1;
                       fileNames{i} = record;
                   end
           end
        end
    end

    fileNames = fileNames(idx);
end


function sortedImages = sortImages(images)

    numbers = zeros(1,length(images));
    for ii = 1:length(images)
        C = strsplit(images{ii}, 'patch');    
        C = strsplit(C{2}, '_');
        numbers(ii) = str2num(C{1});
    end
    [~, a_order] = sort(numbers);
    sortedImages = images(a_order);
end


% -------------------------- funzione Hu_Moments --------------------------
function phi = Hu_Moments(F)
    %INVMOMENTS Compute invariant moments of image.
    %   PHI = INVMOMENTS(F) computes the moment invariants of the image
    %   F. PHI is a seven-element row vector containing the moment
    %   invariants as defined in equations (11.3-17) through (11.3-23) of
    %   Gonzalez and Woods, Digital Image Processing, 2nd Ed.
    %
    %   F must be a 2-D, real, nonsparse, numeric or logical matrix.

    %   Copyright 2002-2004 R. C. Gonzalez, R. E. Woods, & S. L. Eddins
    %   Digital Image Processing Using MATLAB, Prentice-Hall, 2004
    %   $Revision: 1.5 $  $Date: 2003/11/21 14:39:19 $

    if (~ismatrix(F)) || issparse(F) || ~isreal(F) || ~(isnumeric(F) || islogical(F))
       error(['F must be a 2-D, real, nonsparse, numeric or logical ' 'matrix.']);                       
    end

    F = double(F);
    phi = compute_phi(compute_eta(compute_m(F)));            


    % -------------------------- funzione compute_m --------------------------
    function m = compute_m(F)
        [M, N] = size(F);
        [x, y] = meshgrid(1:N, 1:M);

        % Turn x, y, and F into column vectors to make the summations a bit
        % easier to compute in the following.
        x = x(:);
        y = y(:);
        F = F(:);

        % DIP equation (11.3-12)
        m.m00 = sum(F);
        % Protect against divide-by-zero warnings.
        if (m.m00 == 0)
           m.m00 = eps;
        end
        % The other central moments:  
        m.m10 = sum(x .* F);
        m.m01 = sum(y .* F);
        m.m11 = sum(x .* y .* F);
        m.m20 = sum(x.^2 .* F);
        m.m02 = sum(y.^2 .* F);
        m.m30 = sum(x.^3 .* F);
        m.m03 = sum(y.^3 .* F);
        m.m12 = sum(x .* y.^2 .* F);
        m.m21 = sum(x.^2 .* y .* F);            
    end

    % -------------------------- funzione compute_eta --------------------------
    function e = compute_eta(m)

        % DIP equations (11.3-14) through (11.3-16).

        xbar = m.m10 / m.m00;
        ybar = m.m01 / m.m00;

        e.eta11 = (m.m11 - ybar*m.m10) / m.m00^2;
        e.eta20 = (m.m20 - xbar*m.m10) / m.m00^2;
        e.eta02 = (m.m02 - ybar*m.m01) / m.m00^2;
        e.eta30 = (m.m30 - 3 * xbar * m.m20 + 2 * xbar^2 * m.m10) / m.m00^2.5;
        e.eta03 = (m.m03 - 3 * ybar * m.m02 + 2 * ybar^2 * m.m01) / m.m00^2.5;
        e.eta21 = (m.m21 - 2 * xbar * m.m11 - ybar * m.m20 + 2 * xbar^2 * m.m01) / m.m00^2.5;                       
        e.eta12 = (m.m12 - 2 * ybar * m.m11 - xbar * m.m02 + 2 * ybar^2 * m.m10) / m.m00^2.5;                               
    end

    % -------------------------- funzione compute_phi --------------------------
    function phi = compute_phi(e)
        % DIP equations (11.3-17) through (11.3-23).

        phi(1) = e.eta20 + e.eta02;
        phi(2) = (e.eta20 - e.eta02)^2 + 4*e.eta11^2;
        phi(3) = (e.eta30 - 3*e.eta12)^2 + (3*e.eta21 - e.eta03)^2;
        phi(4) = (e.eta30 + e.eta12)^2 + (e.eta21 + e.eta03)^2;
        phi(5) = (e.eta30 - 3*e.eta12) * (e.eta30 + e.eta12) * ...
                 ( (e.eta30 + e.eta12)^2 - 3*(e.eta21 + e.eta03)^2 ) + ...
                 (3*e.eta21 - e.eta03) * (e.eta21 + e.eta03) * ...
                 ( 3*(e.eta30 + e.eta12)^2 - (e.eta21 + e.eta03)^2 );
        phi(6) = (e.eta20 - e.eta02) * ( (e.eta30 + e.eta12)^2 - ...
                                         (e.eta21 + e.eta03)^2 ) + ...
                 4 * e.eta11 * (e.eta30 + e.eta12) * (e.eta21 + e.eta03);
        phi(7) = (3*e.eta21 - e.eta03) * (e.eta30 + e.eta12) * ...
                 ( (e.eta30 + e.eta12)^2 - 3*(e.eta21 + e.eta03)^2 ) + ...
                 (3*e.eta12 - e.eta30) * (e.eta21 + e.eta03) * ...
                 ( 3*(e.eta30 + e.eta12)^2 - (e.eta21 + e.eta03)^2 );           
    end
end    



% -------------------------- Funzione Zernike_Moments --------------------------
function [zernike_mom] = Zernike_Moments(f, ord, repet)
    % Function to compute Zernike moments up to order ord.
    %       f = Input image M x N matrix (can be grayscale or color)
    %       ord = The maximum order (scalar)
    %       rep = repetition  number
    %       zernike_mom = Vector containing all the moments up to order ord
    %   
    % EXAMPLE
    %       
    %   RGBfile = imgetfile();
    %   RGB = imread(RGBfile);
    %   if size(RGB, 3)==1
    %       gray = RGB;
    %   else
    %       gray = rgb2gray(RGB);
    %   end
    %   ord = 4;
    %   repet = 2;
    %   mom = Zernike_Moments(gray, ord, repet);
    %   mom2 = Zernike_Moments(imrotate(gray, 90), ord, repet);
    %   mom3 = Zernike_Moments(imresize(gray, 0.5), ord, repet);
    %   MOM = [mom;mom2;mom3];
    if size(f,3)==3  
        f = rgb2gray(f); 
    end

    f = double(f);
    zernike_mom = [];

    for i=0:ord
        for j=0:repet
            if i>=j && mod((i+j),2)==0
                % Computation of Zernike moment with order i and repetition j 
                [~, AOH, ~] = Zernikmoment_MN(f,i,j); 
                zernike_mom = [zernike_mom AOH];
            end
        end
    end  

    % -------------------------- Funzione Zernikemoment_MN --------------------------
    function [Z, A, Phi] = Zernikmoment_MN(f,n,m)
        % Function to find the Zernike moments for an M x N grayscale image
        %       f = Input image M x N matrix 
        %       n = The order of Zernike moment (scalar)
        %       m = The repetition number of Zernike moment (scalar)
        %       Z = Complex Zernike moment 
        %       A = Amplitude of the moment
        %       Phi = phase (angle) of the moment (in degrees)
        %
        [M, N]  = size(f);
        y = 1:M; x = 1:N;  % x row, y col

        [X,Y] = meshgrid(x,y);
        R = sqrt( (2.*X-N-1).^2/N^2 + (2.*Y-M-1).^2/M^2 );
        Theta = atan2((M-1-2.*Y+2)/M,(2.*X-N+1-2)/N);

        R = (R<=1).*R;
        Rad = radialpoly(R,n,m);    % get the radial polynomial
        Product = f.*Rad.*exp(-1i*m*Theta);
        Z = sum(Product(:));        % calculate the moments

        cnt = nnz(R)+1;             % count the number of pixels inside the unit circle
        Z = (n+1)*Z/cnt;            % normalize the amplitude of moments
        A = abs(Z);                 % calculate the amplitude of the moment
        Phi = angle(Z)*180/pi;      % calculate the phase of the moment (in degrees)
    end

    % -------------------------- Funzione radialpoly --------------------------
    function rad = radialpoly(r,n,m)
        % Function to compute Zernike Polynomials:
        %       r = radius
        %       n = the order of Zernike polynomial
        %       m = the repetition of Zernike moment
        rad = zeros(size(r));                     % Initilization
        for s = 0:(n-abs(m))/2
            c = (-1)^s*factorial(n-s)/(factorial(s)*factorial((n+abs(m))/2-s)*factorial((n-abs(m))/2-s));
        rad = rad + c*r.^(n-2*s);

        end
    end
end        

function v = legmoms_vec(F,ord)
    %LEGMOMS_VEC vector of Legendre moments of an image
    %   v=legmoms_vec(F,ord) computes the vector v of the continuous Legendre
    %   moments of the image F, up to order ord.

    M = legmoms(F,ord);
    idx = fliplr(triu(ones(size(M))));
    v = M(logical(idx))';


    function [M, Md, featM, featMd] = legmoms(F, ord, usesimpson)
        %LEGMOMS Legendre moments of an image
        %   [M,Md]=legmoms(F,ord) computes the continuous Legendre moments M and the
        %   discrete moments Md of the image F, up to order ord.

        n = ord+1;
        type = 'Legendre';
        if nargin ==2
         usesimpson = 1;	% use Simpson to compute integrals (otherwise, trapez. rule)
        end

        %%%%%%%%Aggiunto da Lorenzo
        [m1,m2] = size(F);
        if usesimpson && ~(2*round(m1/2)-m1), F = F(1:end-1,:); end % for Simpson
        if usesimpson && ~(2*round(m2/2)-m2), F = F(:,1:end-1); end % for Simpson
        %%%%%%%%

        [m1,m2] = size(F);
        F = im2double(F);
        x = linspace(-1,1,m1)';
        y = linspace(-1,1,m2)';

        [alfa1,beta1] = opcoef(type,n,m1);	% recursion coefficients
        [alfa2,beta2] = opcoef(type,n,m2);	% recursion coefficients
        P1 = opevmat(alfa1,beta1,x);		% values of polynomials on x
        P2 = opevmat(alfa2,beta2,y);		% values of polynomials on y

        % continuous moments
        M = opcmoms(F,P1,P2,usesimpson);
        Md = 4/m1/m2*(P1'*F*P2);	% discrete Legendre polynomials

        idx = zeros(size(Md));
        idx(1:end-1,1:end-1)= fliplr(triu(ones(size(Md)-1)));

        featM = M(logical(idx))';
        featMd = Md(logical(idx))';
    end

    function M = opcmoms(F,P1,P2,usesimpson,fast)
        %OPCMOMS numerical computation of the moments of an image
        %   M=opcmoms(F,P1,P2) computes the moments of the image F with respect to the
        %   basis functions whose values are contained in the columns of P1 and P2.
        %
        %   M=opcmoms(F,P1,P2,method) chooses between Simpson rule (method=1, default)
        %   and trapezoidal rule (method=0).

        if nargin<4 || isempty(usesimpson), usesimpson = 1; end
        if nargin<5 || isempty(fast), fast = 1; end

        [m1,m2] = size(F);
        n1 = size(P1,2)-1;
        n2 = size(P2,2)-1;

        if usesimpson && ~(2*round(m1/2)-m1), error('m1 must be odd.'), end
        if usesimpson && ~(2*round(m2/2)-m2), error('m2 must be odd.'), end

        if fast
            M = zeros(n1+1,n2+1);
            if usesimpson
                w1 = ones(m1,1); w1(2:2:m1-1)=4; w1(3:2:m1-1)=2;
                h1 = 2/3/(m1-1);
                w2 = ones(m2,1); w2(2:2:m2-1)=4; w2(3:2:m2-1)=2;
                h2 = 2/3/(m2-1);
            else
                w1 = [1 2*ones(1,m1-2) 1]';
                h1 = 1/(m1-1);
                w2 = [1 2*ones(1,m2-2) 1]';
                h2 = 1/(m2-1);
            end
            for k = 1:n1+1
                g = h1*(w1'*(F.*repmat(P1(:,k),1,m2)))';
                M(k,:) = h2*(w2'*(repmat(g,1,n2+1).*P2));
            end
        else	% slow
            M = zeros(n1+1,n2+1);
            g = zeros(m1,1);
            for k = 1:n1+1
                for l = 1:n2+1
                    if usesimpson
                        for i = 1:m2
                            g(i) = simpson(F(:,i).*P1(:,k),-1,1);
                        end
                        M(k,l) = simpson(g.*P2(:,l),-1,1);
                    else
                        for i = 1:m2
                            g(i) = trapezi(F(:,i).*P1(:,k),-1,1);
                        end
                        M(k,l) = trapezi(g.*P2(:,l),-1,1);
                    end
                end
            end
        end             
   end

    function P = opevmat(alfa,beta,x,ortho)
        %OPEVMAT evaluation of orthogonal polynomials
        %   P=opevmat(alpha,beta,x) computes the values of the orthogonal polynomials
        %   defined by the recursion coefficients alpha and beta on the vector of
        %   points x. The j-th column of P contains the values of the polynomial of
        %   degree j-1.
        %
        %   P=opevmat(alpha,beta,x,ortho) specifies the normalizazion of the
        %   polynomials: ortho=1 (default) evaluates orthonormal polynomials, ortho=0
        %   is for monic polynomials.

        if nargin<4, ortho = 1; end

        n = length(alfa)-1;
        m = length(x);
        x = x(:);
        P = ones(m,n+1);

        if ortho	% orthonormal
            P(:,1) = P(:,1)/sqrt(beta(1));
            P(:,2) = (x-alfa(1)).*P(:,1);
            P(:,2) = P(:,2) ./ sqrt(beta(2));
            for k = 2:n
                P(:,k+1) = (x-alfa(k)).*P(:,k) - sqrt(beta(k)).*P(:,k-1);
                P(:,k+1) = P(:,k+1) / sqrt(beta(k+1));
            end
        else		% monic
            P(:,2) = x-alfa(1);
            for k = 2:n
                P(:,k+1) = (x-alfa(k)).*P(:,k) - beta(k).*P(:,k-1);
            end
        end             


    end

    function [alfa,beta] = opcoef(type,n,N)
        %OPCOEF	recursion coefficients for orthogonal polynomials
        %   [alpha,beta]=opcoef(type,n) returns the three-term recursion coefficients
        %   for the first n+1 orthogonal polynomials of some families, as specified by
        %   type: 'Legendre', 'DChebychev' (discrete Chebychev), 'Chebychev'.
        %
        %   If type='DChebychev', a third parameter N, specifying the number of points,
        %   is required. This family is orthogonal in [0,N-1], the other ones in
        %   [-1,1].

        if ischar(type)
            switch lower(type(1:5))
            case 'legen'	% Legendre
                type = 1;
            case 'dcheb'	% discrete Chebychev
                type = 2;
            case 'cheby'	% Chebychev
                type = 3;
            otherwise
                error('unknown type.')
            end
        end

        if type == 2
            if nargin<3, error('the number of points is needed.'), end
            if n+1>N, error('n must be lower than N-1.'), end
        end

        vk = [1:n+1]';

        switch type

        case 1		% Legendre
            beta0 = 2;
            alfa = zeros(n+1,1);
            beta = [beta0; vk.^2./(4*vk.^2-1)];

        case 2		% discrete Chebychev
            beta0 = N;
            alfa = N/2*(1-1/N)*ones(n+1,1);
            beta = [beta0; N^2/4*(1-(vk./N).^2)./(4-1./(vk.^2))];

        case 3		% Chebychev
            beta0 = pi;
            alfa = zeros(n+1,1);
            beta = [beta0; 1/2; ones(n,1)/4];

        otherwise
            error('unknown type.')
        end

    end
end

function [v,i,j] = dchebmoms_vec(F,ord)
    %DCHEBMOMSS_VEC vector of discrete Chebyshev moments of an image
    %   v=dchebmoms_vec(F,ord) computes the vector v of the discrete Chebyshev
    %   moments of the image F, up to order ord.

    M = dchebmoms(F,ord);
    idx = fliplr(triu(ones(size(M))));
    v = M(logical(idx))';


    if nargout > 1 
        idx2 = find(idx);
        [i,j] = ind2sub(size(M),idx2);
        %v = M(sub2ind(size(M),j,i));
        i = i-1;
        j = j-1;
    end

    function [M,Mc,P1,P2] = dchebmoms(F,ord)
        %DCHEBMOMS discrete Chebyshev moments of an image
        %   M=dchebmoms(F,ord) computes the matrix M of the discrete Chebyshev moments
        %   of the image F, up to order ord.
        %   [M,Mc]=dchebmoms(F,ord) computes also the continuous Chebyshev moments Mc.
        %   [M,Md,P1,P2]=dchebmoms(F,ord) returns the evaluation of orthogonal
        %   polynomials on both axes.

        n = ord;
        type = 'DChebyshev';
        usesimpson = 1 && (nargin>1);	% use Simpson to compute integrals

        [m1 m2] = size(F);
        if usesimpson && ~(2*round(m1/2)-m1), m1 = m1-1; end	% for Simpson
        if usesimpson && ~(2*round(m2/2)-m2), m2 = m2-1; end	% for Simpson
        F = double(F(1:m1,1:m2));
        F = mat2gray(F);
        x = [0:m1-1]';
        y = [0:m2-1]';

        [alfa1 beta1] = opcoef(type,n,m1);	% recursion coefficients
        [alfa2 beta2] = opcoef(type,n,m2);	% recursion coefficients
        P1 = opevmat(alfa1,beta1,x);		% values of polynomials on x
        P2 = opevmat(alfa2,beta2,y);		% values of polynomials on y

        % moments
        M = P1'*F*P2;
        if nargout>1, Mc = opcmoms(F,P1,P2,usesimpson); end

    end

    function M = opcmoms(F,P1,P2,usesimpson,fast)
        %OPCMOMS numerical computation of the moments of an image
        %   M=opcmoms(F,P1,P2) computes the moments of the image F with respect to the
        %   basis functions whose values are contained in the columns of P1 and P2.
        %
        %   M=opcmoms(F,P1,P2,method) chooses between Simpson rule (method=1, default)
        %   and trapezoidal rule (method=0).

        if nargin<4 || isempty(usesimpson), usesimpson = 1; end
        if nargin<5 || isempty(fast), fast = 1; end

        [m1,m2] = size(F);
        n1 = size(P1,2)-1;
        n2 = size(P2,2)-1;

        if usesimpson && ~(2*round(m1/2)-m1), error('m1 must be odd.'), end
        if usesimpson && ~(2*round(m2/2)-m2), error('m2 must be odd.'), end

        if fast
            M = zeros(n1+1,n2+1);
            if usesimpson
                w1 = ones(m1,1); w1(2:2:m1-1)=4; w1(3:2:m1-1)=2;
                h1 = 2/3/(m1-1);
                w2 = ones(m2,1); w2(2:2:m2-1)=4; w2(3:2:m2-1)=2;
                h2 = 2/3/(m2-1);
            else
                w1 = [1 2*ones(1,m1-2) 1]';
                h1 = 1/(m1-1);
                w2 = [1 2*ones(1,m2-2) 1]';
                h2 = 1/(m2-1);
            end
            for k = 1:n1+1
                g = h1*(w1'*(F.*repmat(P1(:,k),1,m2)))';
                M(k,:) = h2*(w2'*(repmat(g,1,n2+1).*P2));
            end
        else	% slow
            M = zeros(n1+1,n2+1);
            g = zeros(m1,1);
            for k = 1:n1+1
                for l = 1:n2+1
                    if usesimpson
                        for i = 1:m2
                            g(i) = simpson(F(:,i).*P1(:,k),-1,1);
                        end
                        M(k,l) = simpson(g.*P2(:,l),-1,1);
                    else
                        for i = 1:m2
                            g(i) = trapezi(F(:,i).*P1(:,k),-1,1);
                        end
                        M(k,l) = trapezi(g.*P2(:,l),-1,1);
                    end
                end
            end
        end
    end

    function P = opevmat(alfa,beta,x,ortho)
        %OPEVMAT evaluation of orthogonal polynomials
        %   P=opevmat(alpha,beta,x) computes the values of the orthogonal polynomials
        %   defined by the recursion coefficients alpha and beta on the vector of
        %   points x. The j-th column of P contains the values of the polynomial of
        %   degree j-1.
        %
        %   P=opevmat(alpha,beta,x,ortho) specifies the normalizazion of the
        %   polynomials: ortho=1 (default) evaluates orthonormal polynomials, ortho=0
        %   is for monic polynomials.

        if nargin<4, ortho = 1; end

        n = length(alfa)-1;
        m = length(x);
        x = x(:);
        P = ones(m,n+1);

        if ortho	% orthonormal
            P(:,1) = P(:,1)/sqrt(beta(1));
            P(:,2) = (x-alfa(1)).*P(:,1);
            P(:,2) = P(:,2) ./ sqrt(beta(2));
            for k = 2:n
                P(:,k+1) = (x-alfa(k)).*P(:,k) - sqrt(beta(k)).*P(:,k-1);
                P(:,k+1) = P(:,k+1) / sqrt(beta(k+1));
            end
        else		% monic
            P(:,2) = x-alfa(1);
            for k = 2:n
                P(:,k+1) = (x-alfa(k)).*P(:,k) - beta(k).*P(:,k-1);
            end
        end  
    end

    function [alfa,beta] = opcoef(type,n,N)
        %OPCOEF	recursion coefficients for orthogonal polynomials
        %   [alpha,beta]=opcoef(type,n) returns the three-term recursion coefficients
        %   for the first n+1 orthogonal polynomials of some families, as specified by
        %   type: 'Legendre', 'DChebychev' (discrete Chebychev), 'Chebychev'.
        %
        %   If type='DChebychev', a third parameter N, specifying the number of points,
        %   is required. This family is orthogonal in [0,N-1], the other ones in
        %   [-1,1].

        if ischar(type)
            switch lower(type(1:5))
            case 'legen'	% Legendre
                type = 1;
            case 'dcheb'	% discrete Chebychev
                type = 2;
            case 'cheby'	% Chebychev
                type = 3;
            otherwise
                error('unknown type.')
            end
        end

        if type == 2
            if nargin<3, error('the number of points is needed.'), end
            if n+1>N, error('n must be lower than N-1.'), end
        end

        vk = [1:n+1]';

        switch type

        case 1		% Legendre
            beta0 = 2;
            alfa = zeros(n+1,1);
            beta = [beta0; vk.^2./(4*vk.^2-1)];

        case 2		% discrete Chebychev
            beta0 = N;
            alfa = N/2*(1-1/N)*ones(n+1,1);
            beta = [beta0; N^2/4*(1-(vk./N).^2)./(4-1./(vk.^2))];

        case 3		% Chebychev
            beta0 = pi;
            alfa = zeros(n+1,1);
            beta = [beta0; 1/2; ones(n,1)/4];

        otherwise
            error('unknown type.')
        end             
    end

end



            
