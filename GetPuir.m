% function [ Puir ] = GetPuir(pu_t, pt_i, pt_r)
% 
% Get the probability of user-interest-rating 
% [ M, T ] = size( pu_t );
% N = size( pt_i, 2 );
% R = size( pt_r, 2 );
% 
% for u = 1: M
%     for i = 1: N
%         for r = 1: R
%             p = 0;
%             for t = 1: T
%                 p = p + pu_t( u, t )* pt_i( t, i )* pt_r( t, r );
%             end
%             Puir( u, i, r ) = p;
%         end
%     end
% end
% 
% end