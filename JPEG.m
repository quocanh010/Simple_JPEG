k = load('ForemanCIF_Y.mat');
I = k.yseq;

%Encoding
%Perform level shift down
img = I - 128;

%dct2 

imModified = blockproc(img,[8 8],@(blkStruct) dct2(blkStruct.data));
figure()
imshow(imModified, [])

%Quality factors
qF =[75, 50, 25, 12.5];
alpha = zeros(4,1);
for i = 1:4
    if (qF(i) <= 50)
        alpha(i, 1) = 50 / qF(i);
    else
        alpha(i, 1) = 2 - 50 / qF(i);
    end
end
%Import luminance quantization table\
Quv = [16 11 10 16 24 40 51 61;
        12 12 14 19 26 58 60 55;
        14 13 16 24 40 57 69 56;
        14 17 22 29 51 87 80 62;
        18 22 37 56 68 109 103 77;
        24 35 55 64 81 104 113 92;
        49 64 78 87 103 121 120 101;
        72 92 95 98 112 100 103 99;
     ];
%Uniform Quantisatigon 
Squv = zeros(size(img));

Psnr = [];

for a = 1:4
    
    
    Q_a = alpha(a,1) * Quv
    for i = 1: size(img, 1)/ 8 
        for j = 1:size(img, 2) /8 
            Squv(8*(i-1)+1: 8*(i-1)+8, 8*(j-1)+1 :8*(j-1)+8)=  round(imModified(8*(i-1)+1: 8*(i-1)+8, 8*(j-1)+1 :8*(j-1)+8) ./ Q_a);
        end
    end


    %Path 1 Reconstruct image

    %Dequantization
    Ruv = zeros(size(img)); 
    for i = 1: size(img, 1)/ 8 
        for j = 1:size(img, 2) /8 
            Ruv(8*(i-1)+1: 8*(i-1)+8, 8*(j-1)+1 :8*(j-1)+8)=  (Squv(8*(i-1)+1: 8*(i-1)+8, 8*(j-1)+1 :8*(j-1)+8) .* Q_a);
        end
    end

    %Inverse Dct
    img_reconst = blockproc(Ruv,[8 8],@(blkStruct) idct2(blkStruct.data));
    %Level shift up
    img_reconst_test = img_reconst + 128;
    figure;
    imshow(uint8(img_reconst+128))



    MSE = sum((img_reconst_test - (I)) .^ 2, 'all') / (size(I, 1) * size(I, 2));
    %Calculating PSNR
    PSNR = 10 * log10( 255^2 / MSE);
    Psnr = [PSNR, Psnr];
    
end


    %Path 2 generate Bitstream

    %Zig-zag scan 
count = 0
PRED = 0;
for i = 1: size(img, 1)/ 8 
    for j = 1:size(img, 2) /8 
        z = zig_zag_v(Squv(8*(i-1)+1: 8*(i-1)+8, 8*(j-1)+1 :8*(j-1)+8));
        DC = z(1,1);
        AC = z(1,2:end);
    %Processing
    %Differential endcoding
    %DC
    DIFF = DC - PRED;
    PRED = DC;
    
    %AC run-level-coding
    AC_non_zero = find(AC~=0); % none zero index in AC
    
    rlc_pairs = [];
    
    %Encoding the first AC component
    if  (AC_non_zero(1,1) ~= [])
        
        if(AC_non_zero(1,1) == 1)
            rlc_pairs = [rlc_pairs; 0 AC(AC_non_zero(1,1))]
        else
            rlc_pairs = [rlc_pairs; AC_non_zero(1,1) AC(AC_non_zero(1,1))]
        end
        for i = 2:length(AC_non_zero) 
            n_non_zero = AC_non_zero(1,i) - AC_non_zero(1, i-1) - 1;  % # of non-zero entries in a block
            rlc_pairs = [rlc_pairs; n_non_zero AC(AC_non_zero(1,i))]
        end
    else
        rlc_pairs = [rlc_pairs; 0 0]
    end
    if( (size(rlc_pairs, 1) + sum(rlc_pairs(:, 1), 'all')) ~= 63)
        %EOB coding
        rlc_pairs = [rlc_pairs; 0 0]
    end  
    
    %VLC Huffman code
    
    %DC
    if(DIFF  == 0)
        SSSS = 0;
    else
        SSSS = floor(log2(abs(DIFF) + 1));
    end
    
    %DCprefix
    DCprx = lumaDC{SSSS+1, 1}
   
    %Calculate offset v
    if (DIFF > 0)
        v = (dec2bin(DIFF));
    else
        v = dec2bin(abs(DIFF));
        for i = 1:lenght(v)
            if v(i) == '0'
                v(i) = '1';
            else
                v(i) = '0';
            end
        end
    end
    offset = [];
    for i = 1:length(v)
        if(v(i) == '1')
            offset = [offset 1];
        else
            offset = [offset 0];
        end
    end
    %Codewaord for DC
    if (DIFF == 0)
        cw = [DCprx offset];
    else
        cw = [DCprx offset];
    end
    
    
    count = count + 1;
    
    end
end
