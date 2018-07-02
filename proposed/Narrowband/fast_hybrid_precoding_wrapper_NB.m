function [ FRF,FBB ] = fast_hybrid_precoding_wrapper_NB( Fopt, NRF, FRF_in )

    [Nt, Ns] = size(Fopt);
    y = [];
    
    % FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    if FRF_in == 0
        FRF = exp( 1i*unifrnd(0,2*pi,Nt,NRF) );
    else
        FRF = FRF_in;
    end    
    
    FBB = pinv(FRF) * Fopt;
    
    [FRF,FBB, ~] = fast_hybrid_precoding_NB(Fopt, FRF, FBB);

end