function c = bitcount(v)
%BITCOUNT Count the number of bits.
%   BITCOUNT(V) The number of bits set in V is returned
 c= int64(v); 
 % actually, support is only for int32, as 
 % hex2dec fits only into flintmax

c = bitand(c,              hex2dec('55555555')) + ...
    bitand(bitshift(c,-1), hex2dec('55555555'));
c = bitand(c,              hex2dec('33333333')) + ...
    bitand(bitshift(c,-2), hex2dec('33333333'));
c = bitand(c,              hex2dec('0F0F0F0F')) + ...
    bitand(bitshift(c,-4), hex2dec('0F0F0F0F'));
c = bitand(c,              hex2dec('00FF00FF')) + ...
    bitand(bitshift(c,-8), hex2dec('00FF00FF'));
c = bitand(c,              hex2dec('0000FFFF')) + ...
    bitand(bitshift(c,-16),hex2dec('0000FFFF'));
%c = bitand(c,              hex2dec('00000000FFFFFFFF')) + ...
%    bitand(bitshift(c,-32),hex2dec('00000000FFFFFFFF'));
end

