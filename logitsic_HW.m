%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is the programme to compute statistical complexity measures by
% means of ordinal pattern transition network approaches for nonlinear time
% series analysis, which has been shown in
%
% Characterizing dynamical transitions by statistical complexity measures
% based on ordinal pattern transition networks
% Chaos 31, 033127 (2021); doi: 10.1063/5.0038876
% 
% We present the algorithms by showing different control parameters of the
% Logistic map, which has been summarized in Table I in the above paper.
% The embedding parameters are embedD = 6 and embedT = 6, which are only
% for the illustration purposes. For full discussions, please consult with
% the publised paper. 

% Copyright (C) 2021 East China Normal University
% Author: Yong Zou <yzou@phy.ecnu.edu.cn>
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 

clear all;
close all;
clc;

rSpace = [3.83, 3.679, 3.791, 4];

for iR = 1:4
    r = rSpace(iR); 
    disp(['control parameter r = ', num2str(r, '%1.3f'), ' starts...']); 
    
    xdata = [];
    x0 = rand(1, 1);
    for i=1:1:201000
        xnew = r * x0 * (1.0 - x0);
        x0 = xnew;
        if i > 1000
            xdata(i-1000) = xnew;
        end
    end
    size(xdata);
    Nraw = length(xdata);
    embedD = 6; 
    embedT = embedD;
    
    xEm = [];
    shift_mat_ind = reshape(0:embedT:(embedD-1)*embedT,[],1) * ones(1,Nraw-(embedD-1)*embedT) +...
        ones(embedD, 1) * reshape(1:(Nraw-(embedD-1)*embedT),1,[]);
    shift_mat = xdata(shift_mat_ind);
    xEm = shift_mat';
    
    [xa, idx] = sort(xEm');
    oPd = idx';
    
    % unique order patterns
    pU = unique(oPd, 'rows');
    opTme = zeros(length(oPd), 1);
    op=zeros(1,size(pU, 1));
    for i=1:1:size(pU, 1)
        B = pU(i, :);
        indx = find(all(bsxfun(@eq, oPd, B), 2));
        op(i) = length(indx) / length(oPd);
        opTme(indx) = i;
    end
    
    % node entropy
    L = factorial(embedD);
    M = length(op);
    for i=1:1:L
        if(i<=M)
            op1(1,i)=op(1,i);
        else
            op1(1,i)=0;
        end
    end
    
    Pe1(1, 1:1:L) = 1/L;
    Sp1 = -sum(op .* log(op));
    Smax1 = log(L);
    Ho = Sp1 / Smax1  ;
    Spe1 = -sum(Pe1 .* log(Pe1)) ;
    PPe1 = (op1 + Pe1)/2;
    Sppe1 = -sum(PPe1 .* log(PPe1));
    JS1 = Sppe1 - Sp1/2- Spe1/2;
    Q01 = -2/(((L+1)/L)*log(L+1)-2*log(2*L)+log(L));
    QJ1 = Q01 * JS1;
    Co = QJ1 * Ho;
    
    % transition frequencies
    sNode = opTme(1:end-1);
    eNode = opTme(2:end);
    linkNodes = cat(2, sNode, eNode);
    %remove loop
    dNode=eNode-sNode;
    id=find(dNode==0);
    linkNodes(id,:)=[];
    
    Adj = zeros(size(pU,1),size(pU,1));
    for iNode = 1:1:size(pU,1)
        indx = find(linkNodes(:, 1) == iNode);
        toWhom = linkNodes(indx, 2);
        for jNode = 1:1:size(pU,1)
            if isempty(toWhom)
                Adj(iNode, jNode)=0;
            else
                idx = find(toWhom == jNode);
                Adj(iNode, jNode) = length(idx) / length(linkNodes);
            end
        end
    end
    
    w = Adj(:);
    w=w';
    M=length(w);
    He = 0.0;
    for i = 1:M
        if (w(i) > 0.0)
            He = He + w(i) * log(w(i));
        end
    end
    He = -1.0 * He ;
    
    dd = factorial(embedD);
    % remove loop
    L = dd * (dd - 1);
    
    Pe = 1/L;
    w1 = 0;
    PPe1 = (Pe + w1) / 2;
    Sppe = 0;
    for i=1:M
        PPe(i) = (w(i) + Pe) / 2;
        Sppe = Sppe + PPe(i) * log(PPe(i));
    end
    
    Sppe=Sppe+(L-M)*(PPe1)*log(PPe1);
    Sppe=-Sppe;
    Spe=-L*Pe*log(Pe);
    Sp = He;
    Smax = log(L);
    HW = Sp/Smax ;
    JS = Sppe - Sp/2-Spe/2;
    Q0 = -2/(((L+1)/L)*log(L+1)-2*log(2*L)+log(L));
    QJ = Q0 * JS;
    CW = QJ * HW;
    
    disp(['Ho=', num2str(Ho, '%1.3f'), ', Co=', num2str(Co, '%1.3f'), ', Hw=', num2str(HW, '%1.3f'), ', Cw=', num2str(CW, '%1.3f')]);
    
end